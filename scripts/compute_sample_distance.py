#!/usr/bin/env python3
"""Gamma-distribution distance estimates from gdiff sample windows.

Input (either):
  1) Directory of per-pair TSVs named <genome_a>__<genome_b>.tsv with columns
       qid start end strand reference d lr_bg [lr_ub]
  2) Single concatenated TSV (e.g. samples/gdiff-<cfg>.tsv) with columns
       config genome_a genome_b qid start end strand reference d lr_bg [lr_ub]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats


def _parse_float(s: str) -> float:
    s = s.strip()
    if s.lower() in ("nan", "na", ".", ""):
        return float("nan")
    return float(s)


def summarize_samples(
    d_v: np.ndarray, lr_bg_v: np.ndarray, lr_ub_v: np.ndarray
) -> tuple[float, float]:
    chisq_th = 10
    if len(d_v) == 0:
        raise ValueError("no finite d values")
    if np.median(lr_ub_v) < chisq_th:
        d = float(np.median(d_v))
        return d, 100.0 * (1.0 - d)

    shape, loc, scale = stats.gamma.fit(d_v, floc=0, method="MLE")
    d = shape * scale
    p_v = stats.gamma.cdf(d_v, shape, loc=loc, scale=scale)
    d_v = d_v[np.logical_or(p_v < 0.90, lr_ub_v > chisq_th)]
    if len(d_v) == 0:
        raise ValueError("no d values left after filtering")
    d = float(np.mean(d_v))
    return d, 100.0 * (1.0 - d)


def parse_pair_from_filename(filename: str) -> tuple[str, str] | None:
    stem = Path(filename).stem
    if "__" not in stem:
        return None
    a, b = stem.split("__", 1)
    return a, b


def read_raw_sample_row(
    cols: list[str], d_i: int, lr_bg_i: int, lr_ub_i: int | None
) -> tuple[float, float, float] | None:
    """Return (d, lr_bg, lr_ub) or None if unusable / unmapped."""
    if len(cols) <= max(d_i, lr_bg_i):
        return None
    try:
        d = _parse_float(cols[d_i])
        lr_bg = _parse_float(cols[lr_bg_i]) if len(cols) > lr_bg_i else float("nan")
        if lr_ub_i is not None and len(cols) > lr_ub_i:
            lr_ub = _parse_float(cols[lr_ub_i])
        else:
            lr_ub = float("nan")
    except ValueError:
        return None
    if np.isnan(d):
        return None
    if np.isnan(lr_ub):
        lr_ub = 0.0
    if np.isnan(lr_bg):
        lr_bg = 0.0
    return d, lr_bg, lr_ub


def read_samples_file(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-pair raw sample file: qid start end strand reference d lr_bg [lr_ub]."""
    d_v, lr_bg_v, lr_ub_v = [], [], []
    with path.open() as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            # skip accidental header
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


def read_concat_samples(
    path: Path,
) -> dict[tuple[str, str, str], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Concatenated samples TSV:
      config genome_a genome_b qid start end strand reference d lr_bg [lr_ub]
    Returns map (config, genome_a, genome_b) -> (d, lr_bg, lr_ub).
    """
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
                print(
                    f"  WARNING: {path}:{lineno} has only {len(cols)} columns, skipping",
                    file=sys.stderr,
                )
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


def write_row(
    out, method: str, param_setup: str, genome_a: str, genome_b: str, distance: float, ani_pct: float
) -> None:
    out.write(
        f"{method}\t{param_setup}\t{genome_a}\t{genome_b}\t"
        f"{distance:.8f}\t{ani_pct:.6f}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Gamma-distribution distance estimates from gdiff samples"
    )
    parser.add_argument(
        "input",
        help="Either a directory of per-pair sample TSVs "
        "(named <genome_a>__<genome_b>.tsv) or a single concatenated samples TSV "
        "(config genome_a genome_b qid start end strand reference d lr_bg [lr_ub])",
    )
    parser.add_argument(
        "-m",
        "--method",
        default="gdiff-gamma",
        help="Method label for the output 'method' column (default: gdiff-gamma)",
    )
    parser.add_argument(
        "-p",
        "--param-setup",
        default="",
        help="param_setup label; for concat input, defaults to each row's config "
        "when left empty",
    )
    parser.add_argument("-o", "--output", default=None, help="Output TSV path (default: stdout)")
    args = parser.parse_args()

    inp = Path(args.input)
    if not inp.exists():
        print(f"ERROR: {inp} does not exist", file=sys.stderr)
        sys.exit(1)

    out = open(args.output, "w") if args.output else sys.stdout
    written = 0
    skipped = 0
    try:
        out.write("method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct\n")

        if inp.is_file():
            groups = read_concat_samples(inp)
            if not groups:
                print(f"ERROR: no usable sample rows in {inp}", file=sys.stderr)
                sys.exit(1)
            for (cfg, genome_a, genome_b), (d_v, lr_bg_v, lr_ub_v) in sorted(groups.items()):
                param = args.param_setup if args.param_setup else cfg
                try:
                    distance, ani_pct = summarize_samples(d_v, lr_bg_v, lr_ub_v)
                except ValueError as exc:
                    print(
                        f"  SKIPPING {cfg} {genome_a}__{genome_b}: {exc}",
                        file=sys.stderr,
                    )
                    skipped += 1
                    continue
                write_row(out, args.method, param, genome_a, genome_b, distance, ani_pct)
                written += 1
        else:
            tsv_files = sorted(inp.glob("*.tsv"))
            if not tsv_files:
                print(f"ERROR: no .tsv files found in {inp}", file=sys.stderr)
                sys.exit(1)
            for fpath in tsv_files:
                pair = parse_pair_from_filename(fpath.name)
                if pair is None:
                    print(
                        f"  SKIPPING: cannot parse pair name from {fpath.name}",
                        file=sys.stderr,
                    )
                    skipped += 1
                    continue
                genome_a, genome_b = pair
                d_v, lr_bg_v, lr_ub_v = read_samples_file(fpath)
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
                write_row(
                    out, args.method, args.param_setup, genome_a, genome_b, distance, ani_pct
                )
                written += 1
    finally:
        if args.output:
            out.close()

    print(f"Done: {written} pairs written, {skipped} skipped.", file=sys.stderr)


if __name__ == "__main__":
    main()
