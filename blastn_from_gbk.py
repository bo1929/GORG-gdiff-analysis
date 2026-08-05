#!/usr/bin/env python3
"""
blastn_from_gbk.py — BLASTN all annotated genes from a GenBank file against a subject.

Extracts CDS/gene features from a GBK file (reverse-complementing minus-strand
genes), runs ONE high-sensitivity blastn invocation against the subject, and
writes per-gene rows in the standard 15-column block schema used by
blastn_blocks / minimap_blocks / mummer_blocks (align_util.BLOCK_HEADER):

  q_contig q_start q_end q_strand s_contig ssag s_start s_end s_strand
  nident aln_span identity distance tool raw_score

q_* columns are the gene's coordinates from the GBK.  Every gene gets a row:
genes with no hit get a sentinel row (nident=0, aln_span=0, distance=1, '.'
subject columns).  No mapped/unmapped judgement is made here — distances are
reported as-is so any thresholding can be done downstream.

Optional: provide a separate genome FASTA (--genome) to source gene sequences
from instead of the GBK's embedded sequence (contig ids are matched with and
without version suffixes).

Default (chain) mode: per gene, HSPs are grouped by target contig and the best
colinear chain is kept, so every query base pair is assigned to at most one
HSP; nident/aln_span in the output summarize the chained alignment (coverage =
aln_span / gene length).  This avoids the fragmentary best-single-HSP problem
at low ANI.  With --all-hits, all blastn HSPs are dumped instead (no chaining).

Requires: BioPython, NCBI BLAST+ (makeblastdb, blastn)
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent / "pairwise-mapping"))
from align_util import (  # noqa: E402
    BLOCK_HEADER, check_fastas, identity_distance, managed_workdir,
    require_tools, run_cmd, sag_of)

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════════

TOOL = "blastn-gene"
BLAST_TASK = "blastn"
BLAST_DUST = "no"
BLAST_WORD_SIZE = 7
BLAST_EVALUE = "1000"
BLAST_MAX_TARGET_SEQS = 100  # chain mode: candidate targets per gene
BLAST_MAX_HSPS = 50          # per query-target pair
BLAST_ALL_MAX_TARGET_SEQS = 100   # --all-hits raw dump
BLAST_ALL_MAX_HSPS = 10
CHAIN_MAX_OV = 0.2           # max query-overlap fraction allowed between chained HSPs
BLAST_OUTFMT = "6 qseqid sseqid nident length qstart qend sstart send bitscore pident"

# ═══════════════════════════════════════════════════════════════════════════════
# GBK → gene FASTA
# ═══════════════════════════════════════════════════════════════════════════════

def _unversion(seqid: str) -> str:
    """Strip a trailing version suffix: NZ_CP012345.1 -> NZ_CP012345."""
    head, dot, tail = seqid.rpartition(".")
    return head if dot and tail.isdigit() else seqid


def _load_genome_fasta(path: Optional[Path]) -> Dict[str, str]:
    """Load contig sequences, keyed by first header token AND its unversioned form."""
    if path is None:
        return {}
    genomes: Dict[str, str] = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        key = rec.id.split()[0]
        genomes.setdefault(key, str(rec.seq).upper())
        genomes.setdefault(_unversion(key), str(rec.seq).upper())
    print(f"  loaded {len(genomes)} contig key(s) from genome FASTA", file=sys.stderr)
    return genomes


def extract_genes(gbk_path: Path, genome_fasta: Optional[Path] = None
                  ) -> List[dict]:
    """Extract CDS/gene features from a GenBank file.

    Minus-strand genes are reverse-complemented.  Sequences come from
    *genome_fasta* (version-tolerant contig match) when provided, otherwise
    from the GBK's embedded sequence.
    """
    records = list(SeqIO.parse(str(gbk_path), "genbank"))
    if not records:
        raise SystemExit(f"no records found in {gbk_path}")

    genome_seqs = _load_genome_fasta(genome_fasta)

    genes: List[dict] = []
    for rec in records:
        contig = rec.id
        gseq = genome_seqs.get(contig) or genome_seqs.get(_unversion(contig))
        if gseq is None:
            gseq = str(rec.seq).upper()

        for feat in rec.features:
            if feat.type not in ("CDS", "gene"):
                continue
            loc = feat.location
            try:
                start = int(loc.start) + 1  # 1-based
                end = int(loc.end)
            except (TypeError, ValueError):
                continue  # fuzzy location — skip

            strand = "-" if loc.strand == -1 else "+"  # None/compound -> "+"
            qual = feat.qualifiers

            seqid = None
            for tag in ("locus_tag", "gene", "protein_id", "label"):
                vals = qual.get(tag, [])
                if vals and vals[0].strip():
                    seqid = vals[0].strip()
                    break
            if seqid is None:
                seqid = f"{contig}_{start}_{end}"

            seq = gseq[start - 1:end]
            if strand == "-":
                seq = str(Seq(seq).reverse_complement())
            if len(seq) < 1:
                continue

            genes.append({"seqid": seqid, "contig": contig, "start": start,
                          "end": end, "strand": strand, "sequence": seq})

    # Deduplicate by seqid (keep longest)
    seen: Dict[str, dict] = {}
    for g in genes:
        sid = g["seqid"]
        if sid not in seen or len(g["sequence"]) > len(seen[sid]["sequence"]):
            seen[sid] = g

    print(f"  {len(records)} contig(s), {len(genes)} gene features "
          f"({len(seen)} unique seqids)", file=sys.stderr)
    return list(seen.values())


def write_gene_fasta(genes: List[dict], out_fa: Path) -> None:
    with out_fa.open("w") as fh:
        for g in genes:
            fh.write(f">{g['seqid']} {g['contig']}:{g['start']}-{g['end']}({g['strand']})\n")
            fh.write(f"{g['sequence']}\n")


# ═══════════════════════════════════════════════════════════════════════════════
# Block-schema report
# ═══════════════════════════════════════════════════════════════════════════════

def _block_row(g: dict, q_start: int, q_end: int, sseqid: str,
               s_start, s_end, s_strand: str,
               nident: int, aln_span: int, bitscore: float) -> str:
    identity, distance = identity_distance(nident, aln_span)
    return (f"{g['contig']}\t{q_start}\t{q_end}\t{g['strand']}\t"
            f"{sseqid}\t{sag_of(sseqid)}\t{s_start}\t{s_end}\t{s_strand}\t"
            f"{nident}\t{aln_span}\t{identity:.6f}\t{distance:.6f}\t"
            f"{TOOL}\t{bitscore:.1f}")


def _sentinel_row(g: dict) -> str:
    return (f"{g['contig']}\t{g['start']}\t{g['end']}\t{g['strand']}\t"
            f".\t.\t.\t.\t.\t0\t0\t0.000000\t1.000000\t{TOOL}\t.")


def _chain_target(hsps: List[tuple]) -> Optional[tuple]:
    """Best colinear HSP chain for one (gene, target) group.

    hsps: (nident, aln_len, qs, qe, ss, se) with qs <= qe (gene-relative,
    1-based).  The chain maximizes total nident via DP (weighted colinear
    scheduling, O(n^2) on <= BLAST_MAX_HSPS HSPs).  Returns (s_strand,
    s_start, s_end, matches, covered, q_lo, q_hi, score) or None.  Each query
    bp is then assigned to at most one HSP (highest identity wins), giving
    per-bp match assignment over the gene.
    """
    fw = [h for h in hsps if h[4] <= h[5]]
    rc = [h for h in hsps if h[4] > h[5]]
    strand = "+" if sum(h[1] for h in fw) >= sum(h[1] for h in rc) else "-"
    group = sorted(fw if strand == "+" else rc, key=lambda h: h[2])
    if not group:
        return None

    # dp[i] = (score, prev_ix) of the best chain ending at HSP i
    dp = [(h[0], -1) for h in group]
    for i, (nident, aln_len, qs, qe, ss, se) in enumerate(group):
        best = (nident, -1)
        for j in range(i):
            jn, jal, jqs, jqe, jss, jse = group[j]
            if qs <= jqe - CHAIN_MAX_OV * aln_len:
                continue  # too much query overlap
            if strand == "+" and min(ss, se) < min(jss, jse):
                continue  # not colinear
            if strand == "-" and max(ss, se) > max(jss, jse):
                continue
            cand = dp[j][0] + nident
            if cand > best[0]:
                best = (cand, j)
        dp[i] = best
    i = max(range(len(group)), key=lambda k: dp[k][0])
    score = dp[i][0]
    chain = []
    while i >= 0:
        chain.append(group[i])
        i = dp[i][1]
    chain.reverse()

    # per-bp assignment: highest-identity HSP wins each query position
    q_hi = max(h[3] for h in chain)
    iden = [0.0] * (q_hi + 1)
    for nident, aln_len, qs, qe, ss, se in chain:
        per_bp = nident / aln_len
        for pos in range(qs, qe + 1):
            if per_bp > iden[pos]:
                iden[pos] = per_bp
    matches = sum(iden)
    covered = sum(1 for v in iden if v > 0)
    q_lo = min(h[2] for h in chain)
    s_los = [min(h[4], h[5]) for h in chain]
    s_his = [max(h[4], h[5]) for h in chain]
    return (strand, min(s_los), max(s_his), matches, covered, q_lo, q_hi, score)


def chain_gene(hsps: List[tuple]) -> Optional[tuple]:
    """Best target chain for one gene.  hsps: (sseqid, nident, aln_len, qs,
    qe, ss, se).  Returns (sseqid, s_start, s_end, s_strand, matches, covered,
    q_lo, q_hi, score) or None."""
    by_target: Dict[str, List[tuple]] = {}
    for sseqid, nident, aln_len, qs, qe, ss, se in hsps:
        by_target.setdefault(sseqid, []).append((nident, aln_len, qs, qe, ss, se))
    best, best_key = None, None
    for sseqid, group in by_target.items():
        chained = _chain_target(group)
        if chained is None:
            continue
        strand, s_start, s_end, matches, covered, q_lo, q_hi, score = chained
        key = (matches, covered)
        if best_key is None or key > best_key:
            best_key = key
            s_strand = strand
            best = (sseqid, s_start, s_end, s_strand, matches, covered,
                    q_lo, q_hi, score)
    return best


def write_blocks(blast_out: Path, genes: List[dict], out_path: Path,
                 all_hits: bool) -> None:
    """Parse outfmt 6 and write the block TSV.

    Default (chain) mode: HSPs are chained per gene (best colinear chain,
    per-bp assignment); nident/aln_span summarize the chained alignment.
    --all-hits: every HSP gets a row (no chaining)."""
    per_gene: Dict[str, List[tuple]] = {}
    hits: List[tuple] = []
    with blast_out.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            try:
                qseqid, sseqid = f[0], f[1]
                nident, aln_len = int(f[2]), int(f[3])
                qs, qe = int(f[4]), int(f[5])
                ss, se = int(f[6]), int(f[7])
                bitscore = float(f[8])
            except (ValueError, IndexError):
                continue
            if all_hits:
                s_start, s_end = (ss, se) if ss <= se else (se, ss)
                hits.append((qseqid, sseqid, s_start, s_end,
                             "+" if ss <= se else "-", nident, aln_len, bitscore))
            else:
                per_gene.setdefault(qseqid, []).append(
                    (sseqid, nident, aln_len, qs, qe, ss, se))

    by_id = {g["seqid"]: g for g in genes}
    n_mapped, n_nohit = 0, 0
    with out_path.open("w") as oh:
        oh.write(BLOCK_HEADER + "\n")
        if all_hits:
            seen = set()
            for qseqid, sseqid, s_start, s_end, s_strand, nident, aln_len, bitscore in hits:
                g = by_id.get(qseqid)
                if g is None:
                    continue
                seen.add(qseqid)
                oh.write(_block_row(g, g["start"], g["end"], sseqid, s_start,
                                    s_end, s_strand, nident, aln_len, bitscore) + "\n")
                n_mapped += 1
            for g in genes:
                if g["seqid"] not in seen:
                    oh.write(_sentinel_row(g) + "\n")
                    n_nohit += 1
        else:
            for g in genes:
                hsps = per_gene.get(g["seqid"])
                chained = chain_gene(hsps) if hsps else None
                if chained is None:
                    oh.write(_sentinel_row(g) + "\n")
                    n_nohit += 1
                    continue
                (sseqid, s_start, s_end, s_strand, matches, covered,
                 q_lo, q_hi, score) = chained
                if g["strand"] == "+":  # chained coords are gene-relative
                    q_start, q_end = g["start"] + q_lo - 1, g["start"] + q_hi - 1
                else:  # minus-strand genes were reverse-complemented
                    q_start, q_end = g["end"] - q_hi + 1, g["end"] - q_lo + 1
                oh.write(_block_row(g, q_start, q_end, sseqid, s_start, s_end,
                                    s_strand, round(matches), covered, score) + "\n")
                n_mapped += 1
    print(f"wrote {n_mapped} hit rows + {n_nohit} no-hit -> {out_path}",
          file=sys.stderr)


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-g", "--gbk", required=True, type=Path, help="GenBank input file (query genes)")
    p.add_argument("-s", "--subject", required=True, type=Path, help="subject FASTA")
    p.add_argument("-o", "--out", required=True, type=Path, help="output block TSV")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--genome", type=Path, default=None,
                   help="genome FASTA to source gene sequences from (instead of GBK-embedded)")
    p.add_argument("--all-hits", action="store_true",
                   help="dump ALL blastn HSPs (no chaining) instead of one chained row per gene")
    p.add_argument("--workdir", type=Path, default=None)
    p.add_argument("--keep-tmp", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    require_tools("blastn", "makeblastdb", hint="Install NCBI BLAST+.")
    if not args.gbk.is_file():
        raise SystemExit(f"GBK file not found: {args.gbk}")
    check_fastas(args.subject, Path(args.gbk))

    print(f"reading {args.gbk}", file=sys.stderr)
    genes = extract_genes(args.gbk, genome_fasta=args.genome)
    if not genes:
        raise SystemExit("no gene features found in GBK")

    with managed_workdir(args.workdir, prefix="blastn_genes_", keep_tmp=args.keep_tmp) as work:
        gene_fa = work / "genes.fa"
        print(f"writing {len(genes)} genes -> {gene_fa}", file=sys.stderr)
        write_gene_fasta(genes, gene_fa)

        db = work / "subject_db"
        print(f"building BLAST database from {args.subject}", file=sys.stderr)
        run_cmd(["makeblastdb", "-in", str(args.subject.resolve()),
                 "-dbtype", "nucl", "-out", str(db), "-logfile", os.devnull])

        blast_out = work / "genes.blastn6"
        max_target = BLAST_ALL_MAX_TARGET_SEQS if args.all_hits else BLAST_MAX_TARGET_SEQS
        max_hsps = BLAST_ALL_MAX_HSPS if args.all_hits else BLAST_MAX_HSPS
        print(f"running blastn ({'all-hits' if args.all_hits else 'chain'}, "
              f"high sensitivity)", file=sys.stderr)
        run_cmd(["blastn", "-task", BLAST_TASK,
                 "-query", str(gene_fa), "-db", str(db),
                 "-word_size", str(BLAST_WORD_SIZE),
                 "-dust", BLAST_DUST, "-evalue", BLAST_EVALUE,
                 "-max_target_seqs", str(max_target),
                 "-max_hsps", str(max_hsps),
                 "-num_threads", str(args.threads),
                 "-outfmt", BLAST_OUTFMT], stdout_path=blast_out)

        write_blocks(blast_out, genes, args.out, args.all_hits)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
