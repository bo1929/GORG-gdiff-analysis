#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import minimize
from scipy.stats import gamma

try:
    from tqdm import tqdm
except ImportError:  # progress bar is optional
    def tqdm(iterable, *args, **kwargs):
        return iterable


def meabsdev_winsorize(sample, threshold=2.0):
    sample = np.array(sample)
    median = np.median(sample)
    meabsdev = np.median(np.abs(sample - median))
    if meabsdev == 0:
        return sample
    lower_bound = median - (threshold * meabsdev)
    upper_bound = median + (threshold * meabsdev)
    return np.clip(sample, lower_bound, upper_bound)


def _parse_float(s: str) -> float:
    s = s.strip()
    if s.lower() in ("nan", "na", ".", ""):
        return float("nan")
    return float(s)


def fit_gamma_cencored(sample, upper_bound=0.4):
    censored_mask = np.logical_or(sample >= upper_bound, np.isnan(sample))
    observed_data = np.where(censored_mask, upper_bound, sample)

    censample = stats.CensoredData.right_censored(observed_data, censored_mask)

    shape, loc, scale = stats.gamma.fit(censample, floc=0.0)
    return shape, loc, scale


def fit_gamma_quantiles(sample, q_probs, loc=0.0):
    sample = np.asarray(sample)
    probs = np.asarray(q_probs)
    target_q = np.percentile(sample, probs * 100)

    sample_mean = np.mean(sample) - loc
    sample_var = np.var(sample, ddof=1)

    init_scale = max(1e-6, sample_var / max(1e-6, sample_mean))
    init_shape = max(1e-6, sample_mean / init_scale)

    def objective(params):
        shape, scale = params
        if shape <= 0 or scale <= 0:
            return np.inf
        fitted_q = gamma.ppf(probs, a=shape, loc=loc, scale=scale)
        return np.sum((fitted_q - target_q) ** 2)

    result = minimize(objective, x0=[init_shape, init_scale], bounds=[(1e-5, None), (1e-5, None)], method="L-BFGS-B")

    shape, scale = result.x
    return shape, loc, scale


def drop_nan(d_v, lr_ub_v, lr_bg_v):
    mask = ~np.isnan(d_v)
    d_v = d_v[mask]
    lr_ub_v = lr_ub_v[mask]
    lr_bg_v = lr_bg_v[mask]
    return d_v, lr_ub_v, lr_bg_v


def summarize_samples(
    d_v: np.ndarray,
    lr_bg_v: np.ndarray,
    lr_ub_v: np.ndarray,
    chisq_th: float = 3.841,
    q_filter: float = 0.95,
    stat: str = "mean",
    upper_bound: float = 0.4,
) -> tuple[float, float]:
    if len(d_v) == 0:
        raise ValueError("no finite d values")
    # d_v = d_v[(lr_ub_v < chisq_th) | (lr_bg_v > chisq_th)]
    # d = np.median(d_v)
    # d = np.median(meabsdev_winsorize(d_v[~np.isnan(d_v)]))
    # return d, 100.0 * (1.0 - d)

    lr_ub_v = np.where(np.isnan(d_v), -1e-4, lr_ub_v)
    d_v = np.where(np.isnan(d_v), upper_bound + 1e-4, d_v)
    d_v = np.where(lr_ub_v < chisq_th, upper_bound + 1e-4, d_v)
    shape, loc, scale = fit_gamma_cencored(d_v, upper_bound)

    # Filtering
    p_v = stats.gamma.cdf(d_v, shape, loc=loc, scale=scale)
    d_v = d_v[p_v < q_filter]

    if len(d_v) == 0:
        raise ValueError("no d values left after filtering")

    # Second Gamma fit
    # shape, loc, scale = stats.gamma.fit(d_v, floc=0, method="MLE")
    shape, loc, scale = fit_gamma_cencored(d_v, upper_bound)
    # shape, loc, scale = fit_gamma_quantiles(d_v, [0.25, 0.5, 0.75], loc=0.0)
    if stat == "median":
        d = stats.gamma.median(a=shape, scale=scale)
        # d = np.median(meabsdev_winsorize(d_v[~np.isnan(d_v)]))
    elif stat == "mean":
        d = shape * scale
        # d = np.mean(d_v)  # * (1.0 - 1.0 / (3.0 * shape))
    else:
        raise ValueError("Given statistic is not defined!")
    return d, 100.0 * (1.0 - d)


def parse_pair_from_filename(filename: str) -> tuple[str, str] | None:
    stem = Path(filename).stem
    if "__" not in stem:
        return None
    a, b = stem.split("__", 1)
    return a, b


def read_raw_sample_row(cols: list[str], d_i: int, lr_bg_i: int, lr_ub_i: int | None) -> tuple[float, float, float] | None:
    if len(cols) <= max(d_i, lr_bg_i if lr_bg_i is not None else 0, lr_ub_i if lr_ub_i is not None else 0):
        return None
    try:
        d = _parse_float(cols[d_i])
        lr_bg = _parse_float(cols[lr_bg_i])
        lr_ub = _parse_float(cols[lr_ub_i])
    except ValueError:
        return None
    return d, lr_bg, lr_ub


def read_samples_file(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d_v, lr_bg_v, lr_ub_v = [], [], []
    with path.open() as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if lineno == 1 and cols[0].lower() in ("qid", "config"):
                continue
            parsed = read_raw_sample_row(cols, d_i=5, lr_bg_i=6, lr_ub_i=7)
            if parsed is None:
                continue
            d, lr_bg, lr_ub = parsed
            d_v.append(d)
            lr_bg_v.append(lr_bg)
            lr_ub_v.append(lr_ub)
    return np.asarray(d_v), np.asarray(lr_bg_v), np.asarray(lr_ub_v)


def read_concat_samples(path: Path) -> dict[tuple[str, str, str], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    buckets: dict[tuple[str, str, str], list[tuple[float, float, float]]] = defaultdict(list)
    with path.open() as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if lineno == 1 and cols[0].lower() in ("config", "method"):
                continue
            if len(cols) < 10:
                print(f"  WARNING: {path}:{lineno} has only {len(cols)} columns, skipping", file=sys.stderr)
                continue
            cfg, ga, gb = cols[0], cols[1], cols[2]
            parsed = read_raw_sample_row(cols, d_i=8, lr_bg_i=9, lr_ub_i=10)
            if parsed is None:
                continue
            buckets[(cfg, ga, gb)].append(parsed)

    out: dict[tuple[str, str, str], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for key, rows in buckets.items():
        d_v = np.asarray([r[0] for r in rows])
        lr_bg_v = np.asarray([r[1] for r in rows])
        lr_ub_v = np.asarray([r[2] for r in rows])
        out[key] = (d_v, lr_bg_v, lr_ub_v)
    return out


def write_row(out, method: str, param_setup: str, genome_a: str, genome_b: str, distance: float, ani_pct: float) -> None:
    out.write(f"{method}\t{param_setup}\t{genome_a}\t{genome_b}\t" f"{distance:.8f}\t{ani_pct:.6f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Gamma-distribution distance estimates from gdiff samples")
    parser.add_argument("input", help="Either a directory of per-pair sample TSVs ")
    parser.add_argument("-m", "--method", default="gdiff-gamma", help="Method label for the output 'method' column (default: gdiff-gamma)")
    parser.add_argument("-p", "--param-setup", default="", help="param_setup label; for concat input, defaults to each row's config")
    parser.add_argument("-o", "--output", default=None, help="Output TSV path (default: stdout)")
    parser.add_argument("--chisq", type=float, default=3.841, help="lr_ub chi-square threshold (default 3.841)")
    parser.add_argument("-q", "--quantile", type=float, default=0.95, help="gamma cdf quantile filter (default 0.95)")
    parser.add_argument("-u", "--upper", dest="upper_bound", type=float, default=0.4, help="censored-gamma upper bound (default 0.4)")
    parser.add_argument("--stat", default="mean", choices=["mean", "median"], help="final statistic of the fitted gamma (default: mean)")
    parser.add_argument(
        "--combine",
        default="none",
        choices=["none", "min", "max", "mean"],
        help="combine the two genome directions; 'none' reports both (default: none)",
    )
    args = parser.parse_args()

    inp = Path(args.input)
    if not inp.exists():
        print(f"ERROR: {inp} does not exist", file=sys.stderr)
        sys.exit(1)

    out = open(args.output, "w") if args.output else sys.stdout
    n_written = 0
    n_skipped = 0
    try:
        out.write("method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct\n")

        rows = []  # (param_setup, genome_a, genome_b, d_v, lr_bg_v, lr_ub_v)
        if inp.is_file():
            groups = read_concat_samples(inp)
            if not groups:
                print(f"ERROR: no usable sample rows in {inp}", file=sys.stderr)
                sys.exit(1)
            for (cfg, genome_a, genome_b), (d_v, lr_bg_v, lr_ub_v) in sorted(groups.items()):
                param = args.param_setup if args.param_setup else cfg
                rows.append((param, genome_a, genome_b, d_v, lr_bg_v, lr_ub_v))
        else:
            tsv_files = sorted(inp.glob("*.tsv"))
            if not tsv_files:
                print(f"ERROR: no .tsv files found in {inp}", file=sys.stderr)
                sys.exit(1)
            for fpath in tsv_files:
                pair = parse_pair_from_filename(fpath.name)
                if pair is None:
                    print(f"  SKIPPING: cannot parse pair name from {fpath.name}", file=sys.stderr)
                    n_skipped += 1
                    continue
                genome_a, genome_b = pair
                d_v, lr_bg_v, lr_ub_v = read_samples_file(fpath)
                if len(d_v) == 0:
                    print(f"  SKIPPING: no valid d values in {fpath.name}", file=sys.stderr)
                    n_skipped += 1
                    continue
                rows.append((args.param_setup, genome_a, genome_b, d_v, lr_bg_v, lr_ub_v))

        def one(d_v, lr_bg_v, lr_ub_v):
            return summarize_samples(d_v, lr_bg_v, lr_ub_v, args.chisq, args.quantile, args.stat, args.upper_bound)

        if args.combine == "none":
            for param, genome_a, genome_b, d_v, lr_bg_v, lr_ub_v in tqdm(rows, desc="pairs"):
                try:
                    distance, ani_pct = one(d_v, lr_bg_v, lr_ub_v)
                except ValueError as exc:
                    print(f"  SKIPPING {genome_a}__{genome_b}: {exc}", file=sys.stderr)
                    n_skipped += 1
                    continue
                write_row(out, args.method, param, genome_a, genome_b, distance, ani_pct)
                n_written += 1
        else:
            combine = {"min": min, "max": max, "mean": lambda xs: sum(xs) / len(xs)}[args.combine]
            pairs = defaultdict(list)
            for param, genome_a, genome_b, d_v, lr_bg_v, lr_ub_v in rows:
                key = (param,) + tuple(sorted((genome_a, genome_b)))
                pairs[key].append((genome_a, genome_b, d_v, lr_bg_v, lr_ub_v))
            for (param, ga, gb), dirs_all in tqdm(sorted(pairs.items()), desc="pairs"):
                distances = []
                for genome_a, genome_b, d_v, lr_bg_v, lr_ub_v in dirs_all:
                    try:
                        distance, _ = one(d_v, lr_bg_v, lr_ub_v)
                    except ValueError as exc:
                        print(f"  SKIPPING {genome_a}__{genome_b}: {exc}", file=sys.stderr)
                        n_skipped += 1
                        continue
                    distances.append(distance)
                if not distances:
                    continue
                distance = combine(distances)
                write_row(out, args.method, param, ga, gb, distance, 100.0 * (1.0 - distance))
                n_written += 1

    finally:
        if args.output:
            out.close()

    print(f"Done: {n_written} pairs n_written, {n_skipped} n_skipped.", file=sys.stderr)


if __name__ == "__main__":
    main()
