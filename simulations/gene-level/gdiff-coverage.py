#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import csv
import glob
import os


def merged_length(intervals):
    """Total length of the union of half-open [a, b) intervals."""
    total = 0
    cur_a = cur_b = None
    for a, b in sorted(intervals):
        if cur_b is None or a > cur_b:
            if cur_b is not None:
                total += cur_b - cur_a
            cur_a, cur_b = a, b
        elif b > cur_b:
            cur_b = b
    if cur_b is not None:
        total += cur_b - cur_a
    return total


def pair_coverage(path):
    """-> (query_bp, any_bp, finite_bp) for one pair directory, or None if it has no roll.tsv."""
    fai, roll = os.path.join(path, "qry.fa.fai"), os.path.join(path, "roll.tsv")
    if not os.path.isfile(roll) or not os.path.isfile(fai):
        return None
    query_bp = 0
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2 and f[1].isdigit():
                query_bp += int(f[1])

    windows = collections.defaultdict(list)
    finite = collections.defaultdict(list)
    with open(roll) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0] == "seq":
                continue
            a, b = int(f[1]) - 1, int(f[2])  # 1-based inclusive -> 0-based half-open
            windows[f[0]].append((a, b))
            if f[5] != "nan":
                finite[f[0]].append((a, b))

    any_bp = sum(merged_length(v) for v in windows.values())
    finite_bp = sum(merged_length(v) for v in finite.values())
    return query_bp, any_bp, finite_bp


def pct(part, whole):
    return 100.0 * part / whole if whole else float("nan")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="./results", help="results directory holding <pair>/roll.tsv [./results]")
    p.add_argument("--per-pair", action="store_true", help="also list every pair")
    args = p.parse_args()

    if not os.path.isdir(args.outdir):
        raise SystemExit(f"error: no such directory: {args.outdir}")

    # plan.tsv is optional: it only supplies the (alpha, level) grouping.
    cell_of = {}
    plan = os.path.join(args.outdir, "plan.tsv")
    if os.path.isfile(plan):
        with open(plan) as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                cell_of[r["variant"]] = (r["alpha"], r["level"])

    pairs = []
    skipped = []
    for d in sorted(glob.glob(os.path.join(args.outdir, "*"))):
        if not os.path.isdir(d):
            continue
        got = pair_coverage(d)
        if got is None:
            skipped.append(os.path.basename(d))
            continue
        query_bp, any_bp, finite_bp = got
        pairs.append((os.path.basename(d), query_bp, any_bp, finite_bp))

    if not pairs:
        raise SystemExit(f"error: no <pair>/roll.tsv under {args.outdir}")

    cells = collections.defaultdict(lambda: [0, 0, 0, 0])  # pairs, query, any, finite
    for name, query_bp, any_bp, finite_bp in pairs:
        c = cells[cell_of.get(name, ("?", "?"))]
        c[0] += 1
        c[1] += query_bp
        c[2] += any_bp
        c[3] += finite_bp

    print(f"directory {args.outdir}")
    print(f"pairs     {len(pairs)}" + (f"   (skipped {len(skipped)} without roll.tsv)" if skipped else ""))
    print()

    head = f"  {'alpha':>5} {'level':>7} {'pairs':>6} {'query bp':>14} " f"{'gdiff_cov':>10} {'nan_only':>9} {'no_window':>10}"
    print(head)
    for cell in sorted(cells, key=lambda k: (k[0], float(k[1]) if k[1] not in ("?",) else 0)):
        n, q, a, f = cells[cell]
        print(f"  {cell[0]:>5} {cell[1]:>7} {n:>6} {q:>14,} " f"{pct(f, q):>9.4f}% {pct(a - f, q):>8.4f}% {pct(q - a, q):>9.4f}%")

    q = sum(x[1] for x in pairs)
    a = sum(x[2] for x in pairs)
    f = sum(x[3] for x in pairs)
    print("  " + "-" * (len(head) - 2))
    print(f"  {'ALL':>5} {'':>7} {len(pairs):>6} {q:>14,} " f"{pct(f, q):>9.4f}% {pct(a - f, q):>8.4f}% {pct(q - a, q):>9.4f}%")

    if args.per_pair:
        print()
        print(f"  {'pair':<52} {'query bp':>12} {'gdiff_cov':>10} {'nan_only':>9} {'no_window':>10}")
        for name, query_bp, any_bp, finite_bp in pairs:
            print(
                f"  {name:<52} {query_bp:>12,} {pct(finite_bp, query_bp):>9.4f}% "
                f"{pct(any_bp - finite_bp, query_bp):>8.4f}% "
                f"{pct(query_bp - any_bp, query_bp):>9.4f}%"
            )

    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
