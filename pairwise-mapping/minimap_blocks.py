#!/usr/bin/env python3
"""
minimap_blocks.py -- pairwise blocks via minimap2 (distant-genome baseline).

Edit CONFIG below. CLI overrides CONFIG.

  identity = 100 * nident / aln_span
  distance = 1 - identity / 100

nident from cigar '=' (--eqx); aln_span = query block span. Coords 1-based inclusive.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import align_util as u

# =============================================================================
# CONFIG
# =============================================================================

DEFAULT_PRESET = "asm20"  # asm5 ~0.1%, asm10 ~1%, asm20 ~5% divergence
DEFAULT_THREADS = 8
DEFAULT_MIN_ALN_LEN = 0
DEFAULT_MIN_IDENTITY = 0.0

# Applied with --sensitive (on top of DEFAULT_PRESET)
SENSITIVE_PRESET = {"kmer": 15, "mini_window": 5, "min_chain_score": 10}

# =============================================================================

TOOL = "minimap2"
_CIGAR_OP = re.compile(r"(\d+)([=XIDMSH])")


@dataclass(frozen=True)
class MinimapOptions:
    kmer: Optional[int] = None
    mini_window: Optional[int] = None
    min_chain_score: Optional[int] = None

    def args(self) -> List[str]:
        out: List[str] = []
        if self.kmer is not None:
            out += ["-k", str(self.kmer)]
        if self.mini_window is not None:
            out += ["-w", str(self.mini_window)]
        if self.min_chain_score is not None:
            out += ["-m", str(self.min_chain_score)]
        return out


def nident_from_eqx_cigar(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR_OP.findall(cigar) if op == "=")


def parse_paf_blocks(paf_path: Path) -> Iterator[u.BlockRow]:
    with paf_path.open() as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise SystemExit(f"{paf_path}:{lineno}: expected >=12 PAF columns")
            cigar = next((t[5:] for t in f[12:] if t.startswith("cg:Z:")), None)
            if cigar is None:
                raise SystemExit(f"{paf_path}:{lineno}: missing cg:Z:; need -c --eqx")
            if "=" not in cigar and "X" not in cigar:
                raise SystemExit(f"{paf_path}:{lineno}: cigar lacks =/X; need --eqx")
            qs, qe = int(f[2]) + 1, int(f[3])
            ss, se = int(f[7]) + 1, int(f[8])
            yield (
                f[0], qs, qe, f[4],
                f[5], ss, se, "+",
                nident_from_eqx_cigar(cigar), qe - qs + 1, int(f[11]),
            )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Emit pairwise alignment blocks via minimap2.",
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
        "--preset",
        default=DEFAULT_PRESET,
        choices=("asm5", "asm10", "asm20"),
    )
    sens = p.add_argument_group("sensitivity (edit SENSITIVE_PRESET)")
    s = SENSITIVE_PRESET
    sens.add_argument(
        "--sensitive",
        action="store_true",
        help=f"apply preset -k {s['kmer']} -w {s['mini_window']} -m {s['min_chain_score']}",
    )
    sens.add_argument("-k", "--kmer", type=int, default=None)
    sens.add_argument("-w", "--mini-window", type=int, default=None)
    sens.add_argument("-m", "--min-chain-score", type=int, default=None)
    return p.parse_args(argv)


def resolve_opts(args: argparse.Namespace) -> MinimapOptions:
    opts = dict(SENSITIVE_PRESET) if args.sensitive else {}
    if args.kmer is not None:
        opts["kmer"] = args.kmer
    if args.mini_window is not None:
        opts["mini_window"] = args.mini_window
    if args.min_chain_score is not None:
        opts["min_chain_score"] = args.min_chain_score
    return MinimapOptions(**{k: opts.get(k) for k in ("kmer", "mini_window", "min_chain_score")})


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    u.require_tools("minimap2")
    u.check_fastas(args.query, args.subject)
    opts = resolve_opts(args)

    with u.managed_workdir(
        args.workdir, prefix="minimap_blocks_", keep_tmp=args.keep_tmp
    ) as work:
        paf = work / "aln.paf"
        cmd = [
            "minimap2", "-x", args.preset, *opts.args(),
            "-c", "--eqx", "-t", str(args.threads),
            str(args.subject.resolve()), str(args.query.resolve()),
        ]
        print(f"running {' '.join(cmd[:8])} ... -> {paf}", file=sys.stderr)
        u.run_cmd(cmd, stdout_path=paf)
        n = u.write_blocks(
            parse_paf_blocks(paf),
            args.out,
            tool=TOOL,
            min_aln_len=args.min_aln_len,
            min_identity=args.min_identity,
        )
        print(f"wrote {n} blocks -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
