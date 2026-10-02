#!/usr/bin/env python3
"""STAGE 05 -- merge the per-pair tables into whole-run tables.

  05-collect.py --outdir DIR   # [./results]

IN
  OUTDIR/<pair>/genome.tsv          one row of genome-wide coverage per pair
  OUTDIR/<pair>/resources.tsv       one timed row per command
  OUTDIR/<pair>/genes-joined.tsv    one row per gene, from stage 04
  OUTDIR/plan.tsv                   optional; supplies the pair order

OUT
  OUTDIR/genome-all.tsv        one row per pair
  OUTDIR/resources-all.tsv     one timed row per command, tagged with its pair
  OUTDIR/genes-joined-all.tsv  every pair's per-gene rows
  stdout                       row counts

The gene table is by far the biggest -- one row per gene per pair, so ~1090 rows x 290 pairs for
the full design. It is therefore streamed a pair at a time rather than sorted in memory, and its
rows stay grouped by pair (in plan order when plan.tsv is present, otherwise by pair name) with
each pair's genes in genomic order. The two small tables are loaded and sorted.

Every merge checks that the files agree on their columns before combining them: a mismatched
header would otherwise line up values under the wrong names and fail silently.
"""

from __future__ import annotations

import argparse
import glob
import os


def rows_of(path):
    """Non-comment, non-blank lines of a table."""
    with open(path) as fh:
        return [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]


def merge(pattern, out, key_order, numeric=()):
    """Combine small per-pair tables, sorted, into one file."""
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"  no files matched {pattern}")
        return 0
    head, rows = None, []
    for path in files:
        lines = rows_of(path)
        if not lines:
            continue
        cols = lines[0].split("\t")
        if head is None:
            head = cols
        elif cols != head:
            raise SystemExit(f"{path}: columns differ from {files[0]}; refusing to combine "
                             f"misaligned tables")
        for line in lines[1:]:
            rows.append(line.split("\t"))
    if head is None:
        return 0
    order = {c: i for i, c in enumerate(head)}

    def key(r):
        # Numeric columns sort by value (divergence 0.18 must follow 0.05); text columns
        # lexicographically. Every element is a 3-tuple so the keys stay comparable.
        k = []
        for c in key_order:
            if c not in order:
                continue
            v = r[order[c]]
            if c in numeric:
                try:
                    k.append((0, float(v), ""))
                    continue
                except ValueError:
                    pass
            k.append((1, 0.0, v))
        return tuple(k)

    rows.sort(key=key)
    with open(out, "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"  {os.path.basename(out)}: {len(rows)} rows from {len(files)} pairs")
    return len(rows)


def concat(pattern, out, order=None):
    """Stream one big table from many files, writing the header once.

    `order` is a list of pair directory names; files are emitted in that order, with anything
    not listed after them. Rows keep their order within a file, so each pair's genes stay in
    genomic order.
    """
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"  no files matched {pattern}")
        return 0
    if order:
        rank = {v: i for i, v in enumerate(order)}
        files.sort(key=lambda p: (rank.get(os.path.basename(os.path.dirname(p)), len(rank)), p))
    head, n, used = None, 0, 0
    with open(out, "w") as fh:
        for path in files:
            lines = rows_of(path)
            if not lines:
                continue
            cols = lines[0].split("\t")
            if head is None:
                head = cols
                fh.write("\t".join(head) + "\n")
            elif cols != head:
                raise SystemExit(f"{path}: columns differ from {files[0]}; refusing to "
                                 f"concatenate misaligned tables")
            for line in lines[1:]:
                fh.write(line + "\n")
                n += 1
            used += 1
    print(f"  {os.path.basename(out)}: {n} rows from {used} pairs")
    return n


def plan_order(outdir):
    """Pair names in design order, from plan.tsv if stage 01 left one here."""
    path = os.path.join(outdir, "plan.tsv")
    if not os.path.isfile(path):
        return None
    out, head = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if head is None:
                head = f
                continue
            out.append(f[head.index("variant")])
    return out or None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="./results",
                   help="directory holding the per-pair subdirectories [./results]")
    args = p.parse_args()

    # `pair` is the variant directory name, written by stage 02's timeit.py calls.
    merge(os.path.join(args.outdir, "*", "genome.tsv"),
          os.path.join(args.outdir, "genome-all.tsv"),
          ["seed", "alpha", "level", "gnd"], numeric={"level", "gnd"})
    merge(os.path.join(args.outdir, "*", "resources.tsv"),
          os.path.join(args.outdir, "resources-all.tsv"), ["pair", "label"])
    concat(os.path.join(args.outdir, "*", "genes-joined.tsv"),
           os.path.join(args.outdir, "genes-joined-all.tsv"),
           order=plan_order(args.outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
