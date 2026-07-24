#!/usr/bin/env python3
"""
mummer_blocks.py -- pairwise blocks via nucmer (ground-truth-style WGA).

Edit CONFIG below. CLI overrides CONFIG.

  nucmer -> delta-filter -> show-coords
  identity = 100 * nident / aln_span
  distance = 1 - identity / 100

nident = round(%IDY/100 * LEN2); aln_span = LEN2. Coords 1-based inclusive.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional

import align_util as u

# =============================================================================
# CONFIG
# =============================================================================

DEFAULT_THREADS = 8
DEFAULT_MIN_ALN_LEN = 0
DEFAULT_MIN_IDENTITY = 0.0

DELTA_FILTER_ARGS = ["-1"]  # one-to-one
SHOW_COORDS_ARGS = ["-T", "-H", "-r", "-l", "-c"]

SENSITIVE_PRESET = {
    "maxmatch": True,
    "mincluster": 25,  # -c
    "breaklen": 500,  # -b
    "maxgap": 200,  # -g
    "minmatch": 15,  # -l
}

# =============================================================================

TOOL = "nucmer"


@dataclass(frozen=True)
class NucmerOptions:
    maxmatch: bool = False
    mincluster: Optional[int] = None
    breaklen: Optional[int] = None
    maxgap: Optional[int] = None
    minmatch: Optional[int] = None

    def args(self) -> List[str]:
        out: List[str] = []
        if self.maxmatch:
            out.append("--maxmatch")
        if self.mincluster is not None:
            out += ["-c", str(self.mincluster)]
        if self.breaklen is not None:
            out += ["-b", str(self.breaklen)]
        if self.maxgap is not None:
            out += ["-g", str(self.maxgap)]
        if self.minmatch is not None:
            out += ["-l", str(self.minmatch)]
        return out


def nucmer_supports_threads() -> bool:
    proc = subprocess.run(
        ["nucmer", "-h"], check=False, capture_output=True, text=True
    )
    text = (proc.stdout or "") + (proc.stderr or "")
    return bool(re.search(r"(?m)^\s*-t\b|^\s*--threads\b", text))


def parse_show_coords(coords_path: Path) -> Iterator[u.BlockRow]:
    with coords_path.open() as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                int(f[0])
            except (ValueError, IndexError):
                continue
            if len(f) < 13:
                if len(f) == 12:
                    tags = f[11].split()
                    if len(tags) != 2:
                        raise SystemExit(f"{coords_path}:{lineno}: bad TAGS: {line}")
                    f = f[:11] + tags
                else:
                    raise SystemExit(
                        f"{coords_path}:{lineno}: expected >=13 cols, got {len(f)}"
                    )
            s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            len2, idy = int(f[5]), float(f[6])
            s_contig, q_contig = f[11], f[12]
            s_strand = "+" if s1 <= e1 else "-"
            q_strand = "+" if s2 <= e2 else "-"
            s_start, s_end = (s1, e1) if s1 <= e1 else (e1, s1)
            q_start, q_end = (s2, e2) if s2 <= e2 else (e2, s2)
            aln_span = len2 if len2 > 0 else q_end - q_start + 1
            yield (
                q_contig, q_start, q_end, q_strand,
                s_contig, s_start, s_end, s_strand,
                int(round(idy / 100.0 * aln_span)), aln_span, idy,
            )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Emit pairwise alignment blocks via nucmer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    u.add_io_args(p)
    u.add_block_filter_args(
        p,
        threads=DEFAULT_THREADS,
        min_aln_len=DEFAULT_MIN_ALN_LEN,
        min_identity=DEFAULT_MIN_IDENTITY,
    )
    sens = p.add_argument_group("sensitivity (edit SENSITIVE_PRESET)")
    s = SENSITIVE_PRESET
    sens.add_argument(
        "--sensitive",
        action="store_true",
        help=(
            f"apply preset --maxmatch -c {s['mincluster']} "
            f"-b {s['breaklen']} -g {s['maxgap']} -l {s['minmatch']}"
        ),
    )
    sens.add_argument("--maxmatch", action="store_true")
    sens.add_argument("--mincluster", "-c", type=int, default=None)
    sens.add_argument("--breaklen", "-b", type=int, default=None)
    sens.add_argument("--maxgap", "-g", type=int, default=None)
    sens.add_argument("--minmatch", "-l", type=int, default=None)
    return p.parse_args(argv)


def resolve_opts(args: argparse.Namespace) -> NucmerOptions:
    opts = dict(SENSITIVE_PRESET) if args.sensitive else {}
    if args.maxmatch:
        opts["maxmatch"] = True
    for key in ("mincluster", "breaklen", "maxgap", "minmatch"):
        val = getattr(args, key)
        if val is not None:
            opts[key] = val
    return NucmerOptions(
        maxmatch=bool(opts.get("maxmatch", False)),
        mincluster=opts.get("mincluster"),
        breaklen=opts.get("breaklen"),
        maxgap=opts.get("maxgap"),
        minmatch=opts.get("minmatch"),
    )


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    u.require_tools(
        "nucmer", "delta-filter", "show-coords", hint="Install with: brew install mummer"
    )
    u.check_fastas(args.query, args.subject)
    opts = resolve_opts(args)
    if opts.args():
        print(f"nucmer options: {' '.join(opts.args())}", file=sys.stderr)

    with u.managed_workdir(
        args.workdir, prefix="mummer_blocks_", keep_tmp=args.keep_tmp
    ) as work:
        prefix = work / "aln"
        cmd = ["nucmer"]
        if args.threads != 1 and nucmer_supports_threads():
            cmd += ["-t", str(args.threads)]
        elif args.threads != 1:
            print(
                f"note: this nucmer build has no -t; ignoring --threads {args.threads}",
                file=sys.stderr,
            )
        cmd += opts.args() + [
            "-p", str(prefix),
            str(args.subject.resolve()),
            str(args.query.resolve()),
        ]
        print(f"running nucmer -> {prefix}.delta", file=sys.stderr)
        u.run_cmd(cmd)
        delta = Path(str(prefix) + ".delta")
        if not delta.is_file():
            raise SystemExit(f"nucmer did not produce {delta}")

        filtered = work / "aln.1delta"
        print(f"delta-filter {' '.join(DELTA_FILTER_ARGS)} -> {filtered}", file=sys.stderr)
        u.run_cmd(["delta-filter", *DELTA_FILTER_ARGS, str(delta)], stdout_path=filtered)

        coords = work / "aln.coords"
        print(f"show-coords -> {coords}", file=sys.stderr)
        u.run_cmd(["show-coords", *SHOW_COORDS_ARGS, str(filtered)], stdout_path=coords)

        n = u.write_blocks(
            parse_show_coords(coords),
            args.out,
            tool=TOOL,
            min_aln_len=args.min_aln_len,
            min_identity=args.min_identity,
        )
        print(f"wrote {n} blocks -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
