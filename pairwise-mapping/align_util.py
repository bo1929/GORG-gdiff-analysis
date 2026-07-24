#!/usr/bin/env python3
"""Shared helpers for sliding_blastn / minimap_blocks / mummer_blocks.

Not a config file -- each script keeps its own CONFIG section.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple

# Block TSV schema used by minimap_blocks.py and mummer_blocks.py
BLOCK_HEADER = (
    "q_contig\tq_start\tq_end\tq_strand\t"
    "s_contig\ts_start\ts_end\ts_strand\t"
    "nident\taln_span\tidentity\tdistance\t"
    "tool\traw_score"
)

BlockRow = Tuple[
    str, int, int, str,  # q
    str, int, int, str,  # s
    int, int, object,  # nident, aln_span, raw_score
]


def identity_distance(nident: int, aln_span: int) -> Tuple[float, float]:
    """identity percent and distance = 1 - nident/aln_span."""
    if aln_span <= 0:
        return 0.0, 1.0
    frac = nident / aln_span
    return 100.0 * frac, 1.0 - frac


def require_tools(*names: str, hint: str = "") -> None:
    missing = [n for n in names if shutil.which(n) is None]
    if missing:
        msg = "required tool(s) not found on PATH: " + ", ".join(missing)
        if hint:
            msg += "\n" + hint
        raise SystemExit(msg)


def run_cmd(cmd: Sequence[str], *, stdout_path: Optional[Path] = None) -> None:
    if stdout_path is None:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    else:
        with stdout_path.open("w") as oh:
            proc = subprocess.run(
                cmd, check=False, stdout=oh, stderr=subprocess.PIPE, text=True
            )
    if proc.returncode != 0:
        err = proc.stderr or ""
        raise SystemExit(f"command failed ({proc.returncode}): {' '.join(cmd)}\n{err}")


def check_fastas(query: Path, subject: Path) -> None:
    for path, label in ((query, "query"), (subject, "subject")):
        if not path.is_file():
            raise SystemExit(f"{label} not found: {path}")


def add_io_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("-q", "--query", required=True, type=Path, help="query FASTA")
    p.add_argument("-s", "--subject", required=True, type=Path, help="subject FASTA")
    p.add_argument("-o", "--out", required=True, type=Path, help="output TSV")


def add_block_filter_args(
    p: argparse.ArgumentParser, *, threads: int, min_aln_len: int, min_identity: float
) -> None:
    p.add_argument("-t", "--threads", type=int, default=threads)
    p.add_argument("--min-aln-len", type=int, default=min_aln_len)
    p.add_argument(
        "--min-identity",
        type=float,
        default=min_identity,
        help="drop blocks below this identity (percent)",
    )
    p.add_argument("--workdir", type=Path, default=None)
    p.add_argument("--keep-tmp", action="store_true")


def write_blocks(
    blocks: Iterator[BlockRow],
    out_path: Path,
    *,
    tool: str,
    min_aln_len: int,
    min_identity: float,
) -> int:
    n_out = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as oh:
        oh.write(BLOCK_HEADER + "\n")
        for (
            q_contig,
            q_start,
            q_end,
            q_strand,
            s_contig,
            s_start,
            s_end,
            s_strand,
            nident,
            aln_span,
            raw_score,
        ) in blocks:
            if aln_span < min_aln_len:
                continue
            identity, distance = identity_distance(nident, aln_span)
            if identity < min_identity:
                continue
            oh.write(
                f"{q_contig}\t{q_start}\t{q_end}\t{q_strand}\t"
                f"{s_contig}\t{s_start}\t{s_end}\t{s_strand}\t"
                f"{nident}\t{aln_span}\t{identity:.6f}\t{distance:.6f}\t"
                f"{tool}\t{raw_score}\n"
            )
            n_out += 1
    return n_out


@contextmanager
def managed_workdir(path: Optional[Path], *, prefix: str, keep_tmp: bool):
    own = path is None
    work = Path(tempfile.mkdtemp(prefix=prefix)) if own else path
    assert work is not None
    if not own:
        work.mkdir(parents=True, exist_ok=True)
    try:
        yield work
    finally:
        if own and not keep_tmp:
            shutil.rmtree(work, ignore_errors=True)
        elif own and keep_tmp:
            print(f"kept workdir: {work}", file=sys.stderr)
