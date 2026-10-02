#!/usr/bin/env python3
"""Write a samtools-style .fai index for a FASTA, without needing samtools.

wfmash requires a .fai for its target (and benefits from one for the query), and
samtools is not installed on this machine.  The .fai format is five tab-separated
columns per sequence:

    name  length  offset  linebases  linewidth

``offset`` is the byte offset of the first base, ``linebases`` the number of bases
per line, ``linewidth`` the number of bytes per line including the newline.  A
sequence whose length is not a multiple of ``linebases`` still gets indexed; only
the offset of the first base is needed to seek.

Usage:
    python3 make-fai.py rosina_chr15.fa sara_chr15_region.fa
    python3 make-fai.py --stdout some.fa
"""

from __future__ import annotations

import argparse
import os
import sys


def fai_records(path: str):
    """Yield (name, length, offset, linebases, linewidth) for each sequence."""
    name = None
    length = 0
    offset = 0
    linebases = 0
    linewidth = 0
    in_seq = False

    with open(path, "rb") as fh:
        while True:
            line_start = fh.tell()
            raw = fh.readline()
            if not raw:
                break
            if raw.startswith(b">"):
                if name is not None:
                    yield name, length, offset, linebases, linewidth
                name = raw[1:].split()[0].decode("ascii", "replace")
                length = 0
                in_seq = True
                # The first base line determines linebases/linewidth.
                offset = fh.tell()
                continue
            if not in_seq:
                raise SystemExit(f"{path}: sequence data before the first '>' header")
            body = raw.rstrip(b"\r\n")
            length += len(body)
            if linebases == 0:
                linebases = len(body)
                linewidth = len(raw)
            elif len(body) != linebases:
                # Ragged last line of a sequence is fine; only the first matters.
                pass

    if name is not None:
        yield name, length, offset, linebases, linewidth


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta", nargs="+", help="FASTA file(s) to index")
    p.add_argument("--stdout", action="store_true", help="print instead of writing <fasta>.fai")
    args = p.parse_args()

    for path in args.fasta:
        recs = list(fai_records(path))
        if not recs:
            print(f"{path}: no sequences found", file=sys.stderr)
            return 1
        lines = ["\t".join(str(x) for x in r) for r in recs]
        total = sum(r[1] for r in recs)
        if args.stdout:
            print("\n".join(lines))
            continue
        out = f"{path}.fai"
        if os.path.exists(out):
            print(f"{out} exists; leaving it alone", file=sys.stderr)
            continue
        with open(out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"{out}: {len(recs)} sequence(s), {total:,} bp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
