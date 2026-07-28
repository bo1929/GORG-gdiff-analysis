#!/usr/bin/env python3
"""
blastn_blocks.py -- pairwise blocks via ultra-sensitive BLASTn.

Edit CONFIG below. CLI overrides CONFIG.

Default: full-query blastn HSPs as blocks (same TSV schema as mummer/minimap).
Optional -W/--window: sliding-window mode (window TSV; identity vs window length).

  identity = 100 * nident / aln_span
  distance = 1 - identity / 100

Block mode: aln_span = aligned query span. Coords 1-based inclusive.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import align_util as u

# =============================================================================
# CONFIG
# =============================================================================

DEFAULT_THREADS = 8
DEFAULT_MIN_ALN_LEN = 0
DEFAULT_MIN_IDENTITY = 0.0

# Default = ultra-sensitive (also applied by --sensitive)
BLAST_TASK = "blastn"
BLAST_DUST = "no"
BLAST_MAX_TARGET_SEQS = 5
BLAST_MAX_HSPS = 5
BLAST_OUTFMT = (
    "6 qseqid sseqid nident length qlen qstart qend "
    "sstart send bitscore pident evalue"
)

SENSITIVE_PRESET = {
    "word_size": 7,  # smaller = more sensitive (NCBI default is 11)
    "evalue": "1000",
}

# =============================================================================

TOOL = "blastn"
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


def parse_blast_blocks(blast_path: Path) -> Iterator[u.BlockRow]:
    with blast_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qid, sid = f[0], f[1]
            nident = int(f[2])
            qs, qe = int(f[5]), int(f[6])
            ss, se = int(f[7]), int(f[8])
            bits = float(f[9])
            s_strand = "+" if ss <= se else "-"
            s_start, s_end = (ss, se) if ss <= se else (se, ss)
            q_start, q_end = (qs, qe) if qs <= qe else (qe, qs)
            aln_span = q_end - q_start + 1
            yield (
                qid, q_start, q_end, "+",
                sid, s_start, s_end, s_strand,
                nident, aln_span, bits,
            )


def best_hit_by_nident(blast_path: Path) -> Dict[str, Hit]:
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


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Emit pairwise alignment blocks via ultra-sensitive BLASTn.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    u.add_io_args(p)
    u.add_block_filter_args(
        p,
        threads=DEFAULT_THREADS,
        min_aln_len=DEFAULT_MIN_ALN_LEN,
        min_identity=DEFAULT_MIN_IDENTITY,
    )
    p.add_argument(
        "-W",
        "--window",
        type=int,
        default=None,
        help="if set, sliding-window mode (window TSV instead of blocks)",
    )
    p.add_argument("--step", type=int, default=None, help="window step (default: W/2)")
    p.add_argument("--max-target-seqs", type=int, default=BLAST_MAX_TARGET_SEQS)
    p.add_argument("--max-hsps", type=int, default=BLAST_MAX_HSPS)
    p.add_argument(
        "--by-subject",
        action="store_true",
        help="sliding mode only: best hit per (window, subject genome)",
    )

    sens = p.add_argument_group("sensitivity (edit SENSITIVE_PRESET; default is ultra-sensitive)")
    s = SENSITIVE_PRESET
    sens.add_argument(
        "--sensitive",
        action="store_true",
        help=f"apply preset -word_size {s['word_size']} -evalue {s['evalue']} (already default)",
    )
    sens.add_argument("--word-size", type=int, default=SENSITIVE_PRESET["word_size"])
    sens.add_argument("--evalue", default=SENSITIVE_PRESET["evalue"])
    return p.parse_args(argv)


def resolve_blast_params(args: argparse.Namespace) -> Tuple[int, str]:
    if args.sensitive:
        return int(SENSITIVE_PRESET["word_size"]), str(SENSITIVE_PRESET["evalue"])
    return int(args.word_size), str(args.evalue)


def make_db(subject: Path, db: Path) -> None:
    print(f"building subject DB -> {db}", file=sys.stderr)
    u.run_cmd(
        [
            "makeblastdb",
            "-in", str(subject.resolve()),
            "-dbtype", "nucl",
            "-out", str(db),
            "-parse_seqids",
        ]
    )


def run_sliding(
    args: argparse.Namespace,
    work: Path,
    db: Path,
    word_size: int,
    evalue: str,
) -> int:
    assert args.window is not None
    if args.window <= 0:
        raise SystemExit("--window must be > 0")
    step = args.step if args.step is not None else max(1, args.window // 2)
    if step <= 0:
        raise SystemExit("--step must be > 0")

    query_records = read_fasta(args.query)
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
        word_size=word_size,
        evalue=evalue,
    )
    if args.by_subject:
        best_ps = best_hit_by_nident_per_subject(blast_out)
        n_win = len({w for w, _ in best_ps})
        print(
            f"  {n_win}/{len(meta)} windows with >=1 subject hit "
            f"({len(best_ps)} window-subject pairs)",
            file=sys.stderr,
        )
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
    return 0


def run_blocks(
    args: argparse.Namespace,
    work: Path,
    db: Path,
    word_size: int,
    evalue: str,
) -> int:
    if args.by_subject:
        print("note: --by-subject applies only with -W/--window; ignored", file=sys.stderr)

    blast_out = work / "aln.blastn6"
    print(f"running blastn -> {blast_out}", file=sys.stderr)
    run_blastn(
        args.query.resolve(), db, blast_out,
        threads=args.threads,
        max_target_seqs=args.max_target_seqs,
        max_hsps=args.max_hsps,
        word_size=word_size,
        evalue=evalue,
    )
    n = u.write_blocks(
        parse_blast_blocks(blast_out),
        args.out,
        tool=TOOL,
        min_aln_len=args.min_aln_len,
        min_identity=args.min_identity,
    )
    print(f"wrote {n} blocks -> {args.out}", file=sys.stderr)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    u.require_tools("blastn", "makeblastdb")
    u.check_fastas(args.query, args.subject)
    word_size, evalue = resolve_blast_params(args)
    print(
        f"blastn options: -task {BLAST_TASK} -word_size {word_size} "
        f"-evalue {evalue} -dust {BLAST_DUST}",
        file=sys.stderr,
    )

    with u.managed_workdir(
        args.workdir, prefix="blastn_blocks_", keep_tmp=args.keep_tmp
    ) as work:
        db = work / "subject_db"
        make_db(args.subject, db)
        if args.window is not None:
            return run_sliding(args, work, db, word_size, evalue)
        return run_blocks(args, work, db, word_size, evalue)


if __name__ == "__main__":
    raise SystemExit(main())
