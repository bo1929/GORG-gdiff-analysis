#!/usr/bin/env python3
"""STAGE 06 -- gdiff and wfmash against the exact per-gene truth.

  06-evaluate.py                    # the whole design
  06-evaluate.py --outdir DIR       # defaults: --in DIR/*/genes-joined.tsv -o DIR/metrics.tsv

IN
  --in FILE...     genes-joined.tsv from stage 04; several files are pooled together
  --min-wf-cov F   wfmash coverage gate [0.5]

OUT
  -o FILE          metrics.tsv, one row per (alpha, level, method), with '#' provenance lines
  stdout           the same table

METHODS
  gdiff          every included gene with a finite d
  wfmash         included genes with wfmash coverage >= --min-wf-cov
  intersection   gdiff restricted to the genes wfmash also covers -- the subset both can see

Rows are pooled by (alpha, level), the *target* divergence. Each seed realises a slightly
different divergence, so pooling by the realised gnd would give every seed its own cell.
"""

from __future__ import annotations

import argparse
import glob
import math
import os


def read_tsv(paths):
    rows = []
    for path in paths:
        with open(path) as fh:
            head = None
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if not f or not f[0]:
                    continue
                if head is None:
                    head = f
                    continue
                rows.append(dict(zip(head, f)))
    return rows


def num(x):
    return float("nan") if x in ("", "NA", "nan") else float(x)


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return r

    return pearson(rank(xs), rank(ys))


FIELDS = [
    "alpha",
    "level",
    "gnd",
    "method",
    "n_pairs",
    "n_total",
    "n_included",
    "excl_genes_pct",
    "excl_bp_pct",
    "n",
    "mean_truth",
    "mean_est",
    "bias",
    "mae",
    "rmse",
    "pearson",
    "spearman",
]


def cell(r):
    """The design cell a gene row belongs to: (alpha, target level).

    The key is the *target* level, not the realised ``gnd``. Each seed lands on a slightly
    different divergence, so grouping by ``gnd`` would give every seed its own cell and quietly
    turn a 10-seed run into 10 single-genome runs.
    """
    return (r["alpha"], r["level"])


def stats(rows, method, min_cov):
    """(xs, ys) = (estimate, truth) for the genes this method can see."""
    xs, ys = [], []
    for r in rows:
        if int(r["included"]) != 1:
            continue
        if method in ("wfmash", "intersection") and num(r["wf_cov"]) < min_cov:
            continue
        # Explicit, not a substring test: `"gdiff" in "intersection"` is False, which silently
        # made `intersection` a duplicate of `wfmash` instead of gdiff on the shared gene set.
        est = num(r["wf_rate"] if method == "wfmash" else r["d_gd"])
        t = num(r["truth"])
        if math.isnan(est) or math.isnan(t):
            continue
        xs.append(est)
        ys.append(t)
    return xs, ys


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="./results", help="results directory [./results]")
    p.add_argument("--in", dest="inputs", nargs="+", help="per-gene tables [OUTDIR/*/genes-joined.tsv]")
    p.add_argument("--min-wf-cov", type=float, default=0.5, help="wfmash coverage gate [0.5]")
    p.add_argument("-o", "--out", help="output table [OUTDIR/metrics.tsv]")
    args = p.parse_args()

    if not args.inputs:
        args.inputs = sorted(glob.glob(os.path.join(args.outdir, "*", "genes-joined.tsv")))
        if not args.inputs:
            raise SystemExit(f"no */genes-joined.tsv under {args.outdir}; run stages 02-04 first")
    args.out = args.out or os.path.join(args.outdir, "metrics.tsv")

    rows = read_tsv(args.inputs)
    if not rows:
        raise SystemExit("no rows")
    groups = sorted({cell(r) for r in rows}, key=lambda k: (k[0], float(k[1])))

    out = []
    for alpha, level in groups:
        g = [r for r in rows if r["alpha"] == alpha and r["level"] == level]
        ntot = len(g)
        npairs = len({r.get("seed", "") for r in g})
        # the cell's realised divergence: mean over its pairs, not over its genes
        gnds = sorted({num(r["gnd"]) for r in g})
        mgnd = f"{sum(gnds)/len(gnds):.6f}" if gnds else "NA"
        nin = sum(1 for r in g if int(r["included"]) == 1)
        bp = sum(int(r["gene_len"]) for r in g)
        bpin = sum(int(r["gene_len"]) for r in g if int(r["included"]) == 1)
        for method in ("gdiff", "intersection", "wfmash"):
            xs, ys = stats(g, method, args.min_wf_cov)
            rec = {k: "NA" for k in FIELDS}
            rec.update(
                alpha=alpha,
                level=level,
                gnd=mgnd,
                method=method,
                n_pairs=npairs,
                n_total=ntot,
                n_included=nin,
                excl_genes_pct=f"{100*(ntot-nin)/ntot:.2f}",
                excl_bp_pct=f"{100*(bp-bpin)/bp:.2f}" if bp else "NA",
                n=len(xs),
            )
            if len(xs) >= 3:
                bias = sum(a - b for a, b in zip(xs, ys)) / len(xs)
                rec.update(
                    mean_truth=f"{sum(ys)/len(ys):.6f}",
                    mean_est=f"{sum(xs)/len(xs):.6f}",
                    bias=f"{bias:+.6f}",
                    mae=f"{sum(abs(a-b) for a, b in zip(xs, ys))/len(xs):.6f}",
                    rmse=f"{math.sqrt(sum((a-b)**2 for a, b in zip(xs, ys))/len(xs)):.6f}",
                    pearson=f"{pearson(xs, ys):+.4f}",
                    spearman=f"{spearman(xs, ys):+.4f}",
                )
            out.append(rec)

    with open(args.out, "w") as fh:
        fh.write(f"# inputs: {' '.join(args.inputs)}\n")
        fh.write(f"# wfmash coverage gate: {args.min_wf_cov}\n")
        fh.write("# bias = mean(estimate) - mean(truth); negative = less divergence than simulated\n")
        fh.write("\t".join(FIELDS) + "\n")
        for r in out:
            fh.write("\t".join(str(r[k]) for k in FIELDS) + "\n")

    w = [max(len(k), max(len(str(r[k])) for r in out)) for k in FIELDS]
    print("  ".join(k.ljust(x) for k, x in zip(FIELDS, w)))
    for r in out:
        print("  ".join(str(r[k]).ljust(x) for k, x in zip(FIELDS, w)))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
