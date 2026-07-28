#!/usr/bin/env python3
"""
blastn_from_gbk.py — BLASTN all annotated genes from a GenBank file against a reference.

Extracts CDS/gene features from a GBK file, writes them as a single multi-FASTA,
runs ONE high-sensitivity blastn invocation against the reference, and writes a
concise TSV report.  Every gene gets a row (unmapped → identity=0).

  identity = 100 * nident / aln_length
  distance = 1 - identity / 100

Optional: provide a separate genome FASTA (--genome) to source gene sequences
from instead of the GBK's embedded sequence.  With --all-hits, all blastn HSPs
are reported (multiple rows per gene) instead of only the single best hit.

Requires: BioPython, NCBI BLAST+ (makeblastdb, blastn)
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════════════════════════════════════

BLAST_TASK = "blastn"
BLAST_DUST = "no"
BLAST_WORD_SIZE = 7
BLAST_EVALUE = "1000"
BLAST_MAX_TARGET_SEQS = 1    # best-only mode
BLAST_MAX_HSPS = 1
BLAST_ALL_MAX_TARGET_SEQS = 100   # --all-hits mode
BLAST_ALL_MAX_HSPS = 10
BLAST_OUTFMT = "6 qseqid sseqid nident length qstart qend sstart send bitscore pident"

# ═══════════════════════════════════════════════════════════════════════════════
# GBK → gene FASTA
# ═══════════════════════════════════════════════════════════════════════════════

def _load_genome_fasta(path: Optional[Path]) -> Dict[str, str]:
    """Load contig sequences from a FASTA file, keyed by the first token of each header."""
    if path is None:
        return {}
    genomes: Dict[str, str] = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        key = rec.id.split()[0]
        genomes[key] = str(rec.seq).upper()
    print(f"  loaded {len(genomes)} contig(s) from genome FASTA", file=sys.stderr)
    return genomes


def extract_genes(gbk_path: Path, genome_fasta: Optional[Path] = None
                  ) -> List[dict]:
    """Extract CDS/gene features from a GenBank file.

    Sequences come from *genome_fasta* (keyed by contig id) when provided,
    otherwise from the GBK's embedded sequence.

    Returns list of gene dicts with:
      seqid, contig, start, end, strand, sequence
    """
    records = list(SeqIO.parse(str(gbk_path), "genbank"))
    if not records:
        raise SystemExit(f"no records found in {gbk_path}")

    genome_seqs = _load_genome_fasta(genome_fasta)

    genes: List[dict] = []
    for rec in records:
        contig = rec.id
        # Prefer external genome FASTA; fall back to GBK-embedded sequence
        gseq = genome_seqs.get(contig)
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

            strand = "+" if loc.strand >= 0 else "-"
            qual = feat.qualifiers

            # Build a clean seqid: locus_tag > gene > protein_id > label
            seqid = None
            for tag in ("locus_tag", "gene", "protein_id", "label"):
                vals = qual.get(tag, [])
                if vals and vals[0].strip():
                    seqid = vals[0].strip()
                    break
            if seqid is None:
                seqid = f"{contig}_{start}_{end}"

            # Pull sequence from the (possibly external) genome
            seq = gseq[start - 1:end]

            if len(seq) < 1:
                continue

            genes.append({
                "seqid": seqid,
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
                "sequence": seq,
            })

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
    """Write all gene sequences as a single multi-FASTA."""
    with out_fa.open("w") as fh:
        for g in genes:
            fh.write(f">{g['seqid']} {g['contig']}:{g['start']}-{g['end']}({g['strand']})\n")
            fh.write(f"{g['sequence']}\n")


# ═══════════════════════════════════════════════════════════════════════════════
# BLASTN
# ═══════════════════════════════════════════════════════════════════════════════

def makeblastdb(fasta: Path, db: Path) -> None:
    _run(["makeblastdb", "-in", str(fasta.resolve()), "-dbtype", "nucl",
          "-out", str(db), "-logfile", os.devnull])


def run_blastn(query_fa: Path, db: Path, out: Path, threads: int,
               all_hits: bool = False) -> None:
    max_target = BLAST_ALL_MAX_TARGET_SEQS if all_hits else BLAST_MAX_TARGET_SEQS
    max_hsps = BLAST_ALL_MAX_HSPS if all_hits else BLAST_MAX_HSPS
    _run([
        "blastn", "-task", BLAST_TASK,
        "-query", str(query_fa), "-db", str(db),
        "-word_size", str(BLAST_WORD_SIZE),
        "-dust", BLAST_DUST, "-evalue", BLAST_EVALUE,
        "-max_target_seqs", str(max_target),
        "-max_hsps", str(max_hsps),
        "-num_threads", str(threads),
        "-outfmt", BLAST_OUTFMT,
    ], stdout_path=out)


def _run(cmd: List[str], *, stdout_path: Optional[Path] = None) -> None:
    if stdout_path is None:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    else:
        with stdout_path.open("w") as oh:
            proc = subprocess.run(cmd, check=False, stdout=oh,
                                  stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"command failed ({proc.returncode}): {' '.join(cmd)}\n"
                         f"{proc.stderr or ''}")


# ═══════════════════════════════════════════════════════════════════════════════
# Report
# ═══════════════════════════════════════════════════════════════════════════════

REPORT_HEADER = (
    "reference\tquery_name\tquery_seqid\t"
    "identity\tdistance\tnident\taln_length\t"
    "q_start\tq_end\ts_start\ts_end\tbitscore"
)


def _format_row(ref: str, qname: str, seqid: str,
                identity: float, distance: float, nident: int, aln_len: int,
                qs, qe, ss, se, bitscore) -> str:
    """Format one output row, handling unmapped sentinel values."""
    if aln_len == 0:   # unmapped
        return (f"{ref}\t{qname}\t{seqid}\t"
                f"0.00\t1.000000\t0\t0\t.\t.\t.\t.\t.")
    return (f"{ref}\t{qname}\t{seqid}\t"
            f"{identity:.2f}\t{distance:.6f}\t{nident}\t{aln_len}\t"
            f"{qs}\t{qe}\t{ss}\t{se}\t{bitscore:.1f}")


def write_report(blast_out: Path, genes: List[dict], ref_name: str,
                 query_name: str, out_path: Path,
                 all_hits: bool = False) -> Tuple[int, int]:
    """Parse blastn output and write a TSV report.

    Best-only mode (all_hits=False): one row per gene (best hit by nident).
    All-hits mode (all_hits=True): every HSP gets a row.
    Unmapped genes always get a row (identity=0, distance=1, nident=0,
    aln_length=0, '.' coords).

    Returns (n_hit_rows, n_total_genes).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if all_hits:
        return _write_all_hits(blast_out, genes, ref_name, query_name, out_path)
    else:
        return _write_best_only(blast_out, genes, ref_name, query_name, out_path)


def _write_best_only(blast_out, genes, ref_name, query_name, out_path):
    """Best-hit-per-gene mode.  Every gene gets a row."""
    best: Dict[str, Tuple[int, str, float, int, int, int, int, int]] = {}
    with blast_out.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            try:
                qseqid = f[0]
                nident = int(f[2])
                bitscore = float(f[8])
                hit = (nident, f[1], bitscore,
                       int(f[3]), int(f[4]), int(f[5]),
                       int(f[6]), int(f[7]))
            except (ValueError, IndexError):
                continue
            prev = best.get(qseqid)
            if prev is None or nident > prev[0] or (nident == prev[0] and bitscore > prev[2]):
                best[qseqid] = hit

    n_mapped = 0
    with out_path.open("w") as fh:
        fh.write(REPORT_HEADER + "\n")
        for g in genes:
            hit = best.get(g["seqid"])
            if hit is None:
                fh.write(_format_row(ref_name, query_name, g["seqid"],
                                     0.0, 1.0, 0, 0,
                                     ".", ".", ".", ".", 0.0) + "\n")
                continue
            nident, sseqid, bitscore, aln_len, qs, qe, ss, se = hit
            identity = 100.0 * nident / max(1, aln_len)
            distance = 1.0 - identity / 100.0
            s_start, s_end = (ss, se) if ss <= se else (se, ss)
            q_start, q_end = (qs, qe) if qs <= qe else (qe, qs)
            fh.write(_format_row(ref_name, query_name, g["seqid"],
                                 identity, distance, nident, aln_len,
                                 q_start, q_end, s_start, s_end, bitscore) + "\n")
            n_mapped += 1

    return n_mapped, len(genes)


def _write_all_hits(blast_out, genes, ref_name, query_name, out_path):
    """All-HSPs mode: every blastn hit gets a row.  Unmapped genes included."""
    gene_seqs = {g["seqid"] for g in genes}
    n_rows = 0
    seen_genes: set = set()

    with out_path.open("w") as fh:
        fh.write(REPORT_HEADER + "\n")
        with blast_out.open() as fh_in:
            for line in fh_in:
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                try:
                    qseqid = f[0]
                    nident = int(f[2])
                    bitscore = float(f[8])
                    aln_len = int(f[3])
                    qs, qe = int(f[4]), int(f[5])
                    ss, se = int(f[6]), int(f[7])
                except (ValueError, IndexError):
                    continue
                if qseqid not in gene_seqs:
                    continue
                seen_genes.add(qseqid)
                identity = 100.0 * nident / max(1, aln_len)
                distance = 1.0 - identity / 100.0
                s_start, s_end = (ss, se) if ss <= se else (se, ss)
                q_start, q_end = (qs, qe) if qs <= qe else (qe, qs)
                fh.write(_format_row(ref_name, query_name, qseqid,
                                     identity, distance, nident, aln_len,
                                     q_start, q_end, s_start, s_end, bitscore) + "\n")
                n_rows += 1

        for g in genes:
            if g["seqid"] not in seen_genes:
                fh.write(_format_row(ref_name, query_name, g["seqid"],
                                     0.0, 1.0, 0, 0,
                                     ".", ".", ".", ".", 0.0) + "\n")

    return n_rows, len(genes)


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="BLASTN all annotated genes from a GenBank file against a reference.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-g", "--gbk", required=True, type=Path, help="GenBank input file")
    p.add_argument("-r", "--reference", required=True, type=Path, help="Reference FASTA")
    p.add_argument("-o", "--out", required=True, type=Path, help="Output TSV report")
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--ref-name", default=None, help="Reference label [default: filename stem]")
    p.add_argument("--query-name", default=None, help="Query label [default: GBK accession]")
    p.add_argument("--genome", type=Path, default=None,
                   help="Optional genome FASTA to source gene sequences from (keyed by contig id)")
    p.add_argument("--all-hits", action="store_true",
                   help="Report ALL blastn HSPs (multiple rows per gene) instead of best-only")
    p.add_argument("--keep-tmp", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    for cmd in ("blastn", "makeblastdb"):
        if shutil.which(cmd) is None:
            raise SystemExit(f"{cmd} not found on PATH. Install NCBI BLAST+.")

    if not args.gbk.is_file():
        raise SystemExit(f"GBK file not found: {args.gbk}")
    if not args.reference.is_file():
        raise SystemExit(f"reference not found: {args.reference}")

    # ── Parse GBK ──
    print(f"reading {args.gbk}", file=sys.stderr)
    genes = extract_genes(args.gbk, genome_fasta=args.genome)
    if not genes:
        raise SystemExit("no gene features found in GBK")

    query_name = args.query_name or genes[0]["contig"]
    ref_name = args.ref_name or args.reference.stem

    # ── Work in temp directory ──
    tmpdir = tempfile.mkdtemp(prefix="blastn_genes_")
    try:
        gene_fa = Path(tmpdir) / "genes.fa"
        print(f"writing {len(genes)} genes -> {gene_fa}", file=sys.stderr)
        write_gene_fasta(genes, gene_fa)

        db = Path(tmpdir) / "ref_db"
        print(f"building BLAST database from {args.reference}", file=sys.stderr)
        makeblastdb(args.reference, db)

        blast_out = Path(tmpdir) / "genes.blastn6"
        mode = "all-hits" if args.all_hits else "best-only"
        print(f"running blastn ({mode}, high sensitivity) -> {blast_out}",
              file=sys.stderr)
        run_blastn(gene_fa, db, blast_out, args.threads, all_hits=args.all_hits)

        n_hit_rows, n_total = write_report(
            blast_out, genes, ref_name, query_name, args.out,
            all_hits=args.all_hits,
        )
        pct = 100.0 * n_hit_rows / max(1, n_total) if not args.all_hits else 0
        if args.all_hits:
            print(f"wrote {args.out}  ({n_hit_rows} hit rows, {n_total} genes)",
                  file=sys.stderr)
        else:
            print(f"wrote {args.out}  ({n_hit_rows}/{n_total} genes mapped, {pct:.1f}%)",
                  file=sys.stderr)

    finally:
        if args.keep_tmp:
            print(f"kept workdir: {tmpdir}", file=sys.stderr)
        else:
            shutil.rmtree(tmpdir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
