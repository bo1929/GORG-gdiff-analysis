#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scipy import stats


def read_samples(path: str) -> np.ndarray:
    d_v = []
    lr_bg_v = []
    lr_ub_v = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                print(
                    f"  WARNING: {path}:{lineno} has only {len(cols)} columns, skipping",
                    file=sys.stderr,
                )
                continue
            try:
                d = float(cols[5])
                lr_bg = float(cols[6])
                lr_ub = float(cols[7])
            except ValueError:
                print(f"  WARNING: {path}:{lineno} could not be parsed, skipping", file=sys.stderr)
                continue
            if not np.isnan(d):
                d_v.append(d)
                lr_bg_v.append(lr_bg)
                lr_ub_v.append(lr_ub)
    return np.array(d_v), np.array(lr_bg_v), np.array(lr_ub_v)


def summarize_samples(
    d_v: np.ndarray, lr_bg_v: np.ndarray, lr_ub_v: np.ndarray
) -> tuple[float, float]:
    chisq_th = 10
    # d = np.median(d_v)
    # ani_pct = 100.0 * (1.0 - d)
    # return d, ani_pct
    # d_v = d_v[d_v < 0.25]
    # d_v = d_v[d_v < 0.25]
    x = len(d_v)
    d_median = np.median(d_v)
    # d_v = d_v[lr_ub_v > chisq_th]
    ix = np.argsort(d_v)
    ix_median = ix[d_v.shape[0] // 2]
    if np.median(lr_ub_v) < chisq_th:
        d = np.median(d_v)
        return d, 100.0 * (1.0 - d)

    # if lr_ub_v[ix_median] > chisq_th:
    #     d_v = d_v[lr_ub_v > chisq_th]
    # else:
    #     d = np.median(d_v)
    #     return d, 100.0 * (1.0 - d)

    # 2
    # if len(d_v) < (len(d_v) / 10):
    #     # raise ValueError(f"too few valid d_v values ({len(d_v)}) to fit a Gamma distribution")
    #     return d_median, 100.0 * (1.0 - d)
    # else:
    #     # Gamma MLE fit — scipy fits shape (a), loc, scale.
    #     # For a standard Gamma: mean = shape * scale + loc.
    #     # We fix loc = 0 so the support starts at 0.
    #     shape, loc, scale = stats.gamma.fit(d_v, floc=0, method="MLE")
    #     d = shape * scale  # = shape / rate, since rate = 1/scale
    #     # d = np.median(d_v)
    #     # d = np.mean(d_v)
    # 2.5
    # d = max(d_median, d)

    # 3
    shape, loc, scale = stats.gamma.fit(d_v, floc=0, method="MLE")
    d = shape * scale  # = shape / rate, since rate = 1/scale
    p_v = stats.gamma.cdf(d_v, shape, loc=loc, scale=scale)
    d_v = d_v[np.logical_or(p_v < 0.90, lr_ub_v > chisq_th)]
    # d_v = d_v[lr_ub_v > chisq_th]
    d = np.mean(d_v)
    # d = np.median(d_v)
    y = len(d_v)
    # print(y / x, file=sys.stderr)
    return d, 100.0 * (1.0 - d)


def parse_pair_from_filename(filename: str) -> tuple[str, str] | None:
    """
    Filenames look like  AG-325-K03__AG-905-E11.tsv  (double underscore
    separates the two genome names).  Return (genome_a, genome_b) or None.
    """
    stem = Path(filename).stem  # strip .tsv
    if "__" not in stem:
        return None
    parts = stem.split("__", 1)
    return parts[0], parts[1]


def main():
    parser = argparse.ArgumentParser(
        description="Gamma-distribution distance estimates from gdiff samples"
    )
    parser.add_argument(
        "sample_dir",
        help="Directory containing per-pair gdiff-sample TSV files "
        "(e.g. results/gdiff-samples/long-window/)",
    )
    parser.add_argument(
        "-m",
        "--method",
        default="gdiff-gamma",
        help="Method label written to the output 'method' column " "(default: gdiff-gamma)",
    )
    parser.add_argument(
        "-p",
        "--param-setup",
        default="",
        help="Parameter-setup label written to the output 'param_setup' column "
        "(e.g. 'k=27,w=35,-l=1000,n=200')",
    )
    parser.add_argument("-o", "--output", default=None, help="Output TSV path (default: stdout)")
    args = parser.parse_args()

    sample_dir = Path(args.sample_dir)
    if not sample_dir.is_dir():
        print(f"ERROR: {sample_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    tsv_files = sorted(sample_dir.glob("*.tsv"))
    if not tsv_files:
        print(f"ERROR: no .tsv files found in {sample_dir}", file=sys.stderr)
        sys.exit(1)

    # Write output
    out = open(args.output, "w") if args.output else sys.stdout
    try:
        out.write("method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct\n")
        skipped = 0
        for fpath in tsv_files:
            pair = parse_pair_from_filename(fpath.name)
            if pair is None:
                print(f"  SKIPPING: cannot parse pair name from {fpath.name}", file=sys.stderr)
                skipped += 1
                continue

            genome_a, genome_b = pair
            d_v, lr_bg_v, lr_ub_v = read_samples(str(fpath))

            if len(d_v) == 0:
                print(f"  SKIPPING: no valid d values in {fpath.name}", file=sys.stderr)
                skipped += 1
                continue

            try:
                distance, ani_pct = summarize_samples(d_v, lr_bg_v, lr_ub_v)
            except ValueError as exc:
                print(f"  SKIPPING {fpath.name}: {exc}", file=sys.stderr)
                skipped += 1
                continue

            out.write(
                f"{args.method}\t{args.param_setup}\t"
                f"{genome_a}\t{genome_b}\t"
                f"{distance:.8f}\t{ani_pct:.6f}\n"
            )

    finally:
        if args.output:
            out.close()

    if skipped:
        print(
            f"Done: {len(tsv_files) - skipped} pairs written, " f"{skipped} skipped.",
            file=sys.stderr,
        )
    else:
        print(f"Done: {len(tsv_files)} pairs written.", file=sys.stderr)


if __name__ == "__main__":
    main()
