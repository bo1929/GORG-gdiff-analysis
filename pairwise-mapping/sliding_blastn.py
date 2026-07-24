#!/usr/bin/env python3
"""
sliding_blastn.py -- sliding-window BLASTn identity (query vs subject).

Edit CONFIG below. CLI overrides CONFIG.

Windows:
  identity = 100 * nident / window_length
  distance = 1 - identity / 100

--full-out (best HSP per contig):
  identity = 100 * nident / aligned_query_span
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import align_util as u

# =============================================================================
# CONFIG
# =============================================================================

DEFAULT_THREADS = 1
BLAST_TASK = "blastn"
BLAST_WORD_SIZE = 7  # smaller = more sensitive
BLAST_EVALUE = "1000"
BLAST_DUST = "no"
BLAST_MAX_TARGET_SEQS = 5
BLAST_MAX_HSPS = 5
BLAST_OUTFMT = (
    "6 qseqid sseqid nident length qlen qstart qend "
    "sstart send bitscore pident evalue"
)

# =============================================================================

Hit = Tuple[int, str, float, int, int, int, int]  # nident, sid, bits, qs, qe, ss, se


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip().upper())
        if name is not None:
            records.append((name, "".join(chunks)))
    if not records:
        raise SystemExit(f"no sequences found in {path}")
    return records


def windows_for_seq(seq_len: int, window: int, step: int) -> Iterable[Tuple[int, int]]:
    if seq_len < window:
        if seq_len > 0:
            yield 0, seq_len
        return
    start = 0
    while start + window <= seq_len:
        yield start, start + window
        start += step
    last_full = ((seq_len - window) // step) * step
    if last_full + window < seq_len:
        yield seq_len - window, seq_len


def write_window_fasta(
    records: List[Tuple[str, str]], window: int, step: int, out_fa: Path
) -> List[Tuple[str, str, int, int, int]]:
    meta: List[Tuple[str, str, int, int, int]] = []
    with out_fa.open("w") as oh:
        for contig, seq in records:
            for start0, end0 in windows_for_seq(len(seq), window, step):
                start1, end1 = start0 + 1, end0
                wid = f"{contig}:{start1}-{end1}"
                oh.write(f">{wid}\n{seq[start0:end0]}\n")
                meta.append((wid, contig, start1, end1, end0 - start0))
    return meta


def run_blastn(
    query_fa: Path,
    db: Path,
    out_path: Path,
    *,
    threads: int,
    max_target_seqs: int,
    max_hsps: int,
    word_size: int,
    evalue: str,
) -> None:
    u.run_cmd(
        [
            "blastn",
            "-task", BLAST_TASK,
            "-query", str(query_fa),
            "-db", str(db),
            "-word_size", str(word_size),
            "-dust", BLAST_DUST,
            "-evalue", evalue,
            "-max_target_seqs", str(max_target_seqs),
            "-max_hsps", str(max_hsps),
            "-num_threads", str(threads),
            "-outfmt", BLAST_OUTFMT,
        ],
        stdout_path=out_path,
    )


def best_hit_by_nident(blast_path: Path) -> Dict[str, Hit]:
    """Best HSP per query id (single subject or global best)."""
    best: Dict[str, Hit] = {}
    with blast_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qid, nident, bits = f[0], int(f[2]), float(f[9])
            hit: Hit = (nident, f[1], bits, int(f[5]), int(f[6]), int(f[7]), int(f[8]))
            prev = best.get(qid)
            if prev is None or nident > prev[0] or (nident == prev[0] and bits > prev[2]):
                best[qid] = hit
    return best


def best_hit_by_nident_per_subject(blast_path: Path) -> Dict[Tuple[str, str], Hit]:
    """Best HSP per (query id, subject genome id)."""
    best: Dict[Tuple[str, str], Hit] = {}
    with blast_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qid, sid = f[0], f[1]
            ssag = u.sag_of(sid)
            nident, bits = int(f[2]), float(f[9])
            key = (qid, ssag)
            hit: Hit = (nident, sid, bits, int(f[5]), int(f[6]), int(f[7]), int(f[8]))
            prev = best.get(key)
            if prev is None or nident > prev[0] or (nident == prev[0] and bits > prev[2]):
                best[key] = hit
    return best


def write_window_table(
    meta: List[Tuple[str, str, int, int, int]], best: Dict[str, Hit], out_path: Path
) -> None:
    with out_path.open("w") as oh:
        oh.write(
            "window_id\tcontig\tq_start\tq_end\twindow_len\t"
            "nident\tidentity\tdistance\tssag\tsseqid\tq_aln_start\tq_aln_end\t"
            "s_start\ts_end\tbitscore\n"
        )
        for wid, contig, start1, end1, win_len in meta:
            hit = best.get(wid)
            if hit is None:
                oh.write(
                    f"{wid}\t{contig}\t{start1}\t{end1}\t{win_len}\t"
                    f"0\t0.000000\t1.000000\t*\t*\t*\t*\t*\t*\t*\n"
                )
                continue
            nident, sid, bitscore, qs, qe, ss, se = hit
            identity, distance = u.identity_distance(nident, win_len)
            oh.write(
                f"{wid}\t{contig}\t{start1}\t{end1}\t{win_len}\t"
                f"{nident}\t{identity:.6f}\t{distance:.6f}\t"
                f"{u.sag_of(sid)}\t{sid}\t{qs}\t{qe}\t{ss}\t{se}\t{bitscore:.1f}\n"
            )


def write_window_table_by_subject(
    meta: List[Tuple[str, str, int, int, int]],
    best: Dict[Tuple[str, str], Hit],
    out_path: Path,
) -> int:
    """One row per (window, subject) with a hit; empty windows omitted."""
    by_window: Dict[str, List[str]] = {}
    for wid, ssag in best:
        by_window.setdefault(wid, []).append(ssag)
    for ssags in by_window.values():
        ssags.sort()

    n = 0
    with out_path.open("w") as oh:
        oh.write(
            "window_id\tcontig\tq_start\tq_end\twindow_len\t"
            "nident\tidentity\tdistance\tssag\tsseqid\tq_aln_start\tq_aln_end\t"
            "s_start\ts_end\tbitscore\n"
        )
        for wid, contig, start1, end1, win_len in meta:
            for ssag in by_window.get(wid, ()):
                nident, sid, bitscore, qs, qe, ss, se = best[(wid, ssag)]
                identity, distance = u.identity_distance(nident, win_len)
                oh.write(
                    f"{wid}\t{contig}\t{start1}\t{end1}\t{win_len}\t"
                    f"{nident}\t{identity:.6f}\t{distance:.6f}\t"
                    f"{ssag}\t{sid}\t{qs}\t{qe}\t{ss}\t{se}\t{bitscore:.1f}\n"
                )
                n += 1
    return n



def write_full_summary(
    query_records: List[Tuple[str, str]], best: Dict[str, Hit], out_path: Path
) -> None:
    with out_path.open("w") as oh:
        oh.write(
            "qseqid\tqlen\taln_qlen\tnident\tidentity\tdistance\tsseqid\t"
            "q_aln_start\tq_aln_end\ts_start\ts_end\tbitscore\n"
        )
        total_qlen = total_aln = total_nident = 0
        for qid, seq in query_records:
            qlen = len(seq)
            total_qlen += qlen
            hit = best.get(qid)
            if hit is None:
                oh.write(
                    f"{qid}\t{qlen}\t0\t0\t0.000000\t1.000000\t*\t*\t*\t*\t*\t*\n"
                )
                continue
            nident, sid, bitscore, qs, qe, ss, se = hit
            aln_qlen = abs(qe - qs) + 1
            total_aln += aln_qlen
            total_nident += nident
            identity, distance = u.identity_distance(nident, aln_qlen)
            oh.write(
                f"{qid}\t{qlen}\t{aln_qlen}\t{nident}\t{identity:.6f}\t"
                f"{distance:.6f}\t{sid}\t{qs}\t{qe}\t{ss}\t{se}\t{bitscore:.1f}\n"
            )
        if total_aln:
            g_ident, g_dist = u.identity_distance(total_nident, total_aln)
            oh.write(
                f"__ALL__\t{total_qlen}\t{total_aln}\t{total_nident}\t"
                f"{g_ident:.6f}\t{g_dist:.6f}\t*\t*\t*\t*\t*\t*\n"
            )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sliding-window BLASTn identities.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    u.add_io_args(p)
    p.add_argument("-W", "--window", type=int, required=True, help="window size (bp)")
    p.add_argument("--step", type=int, default=None, help="step size (default: W/2)")
    p.add_argument("--full-out", type=Path, default=None, help="full-query summary TSV")
    p.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS)
    p.add_argument("--word-size", type=int, default=BLAST_WORD_SIZE)
    p.add_argument("--evalue", default=BLAST_EVALUE)
    p.add_argument("--max-target-seqs", type=int, default=BLAST_MAX_TARGET_SEQS)
    p.add_argument("--max-hsps", type=int, default=BLAST_MAX_HSPS)
    p.add_argument(
        "--by-subject",
        action="store_true",
        help="keep best hit per (window, subject genome); for multi-genome subject DBs",
    )
    p.add_argument("--workdir", type=Path, default=None)
    p.add_argument("--keep-tmp", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    u.require_tools("blastn", "makeblastdb")
    if args.window <= 0:
        raise SystemExit("--window must be > 0")
    step = args.step if args.step is not None else max(1, args.window // 2)
    if step <= 0:
        raise SystemExit("--step must be > 0")
    u.check_fastas(args.query, args.subject)

    with u.managed_workdir(
        args.workdir, prefix="sliding_blastn_", keep_tmp=args.keep_tmp
    ) as work:
        query_records = read_fasta(args.query)
        db = work / "subject_db"
        print(f"building subject DB -> {db}", file=sys.stderr)
        u.run_cmd(
            [
                "makeblastdb",
                "-in", str(args.subject.resolve()),
                "-dbtype", "nucl",
                "-out", str(db),
                "-parse_seqids",
            ]
        )

        win_fa = work / "windows.fa"
        print(f"writing windows (W={args.window}, step={step}) -> {win_fa}", file=sys.stderr)
        meta = write_window_fasta(query_records, args.window, step, win_fa)
        print(f"  {len(meta)} windows", file=sys.stderr)

        blast_out = work / "windows.blastn6"
        print(f"running blastn (windows) -> {blast_out}", file=sys.stderr)
        run_blastn(
            win_fa, db, blast_out,
            threads=args.threads,
            max_target_seqs=args.max_target_seqs,
            max_hsps=args.max_hsps,
            word_size=args.word_size,
            evalue=args.evalue,
        )
        if args.by_subject:
            best_ps = best_hit_by_nident_per_subject(blast_out)
            n_win = len({w for w, _ in best_ps})
            print(f"  {n_win}/{len(meta)} windows with >=1 subject hit "
                  f"({len(best_ps)} window-subject pairs)", file=sys.stderr)
            n = write_window_table_by_subject(meta, best_ps, args.out)
            print(f"wrote {n} rows -> {args.out}", file=sys.stderr)
        else:
            best = best_hit_by_nident(blast_out)
            print(
                f"  {sum(1 for w, *_ in meta if w in best)}/{len(meta)} windows with a hit",
                file=sys.stderr,
            )
            write_window_table(meta, best, args.out)
            print(f"wrote {args.out}", file=sys.stderr)

        if args.full_out is not None:
            if args.by_subject:
                print("note: --full-out with --by-subject uses global best per contig",
                      file=sys.stderr)
            full_blast = work / "full.blastn6"
            print(f"running blastn (full query) -> {full_blast}", file=sys.stderr)
            run_blastn(
                args.query.resolve(), db, full_blast,
                threads=args.threads,
                max_target_seqs=args.max_target_seqs,
                max_hsps=args.max_hsps,
                word_size=args.word_size,
                evalue=args.evalue,
            )
            write_full_summary(query_records, best_hit_by_nident(full_blast), args.full_out)
            print(f"wrote {args.full_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
