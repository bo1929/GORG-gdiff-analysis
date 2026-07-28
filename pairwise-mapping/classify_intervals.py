#!/usr/bin/env python3
"""
classify_intervals.py — Classify genomic intervals from pairwise mapping tools.

Takes a block/window TSV from blastn_blocks, minimap_blocks, mummer_blocks, or
window_blastn, fits a background model to the distance distribution, and adds
statistical columns + three-way classification:

  unmapped    — no alignment hit in this interval
  background  — distance within expected range
  significant — distance is an outlier (e.g. unusually high similarity → HGT)

Two statistical tests:
  robust-z        Non-parametric: median / MAD → normal p-value (default)
  gamma-quantiles Parametric: fit Gamma(k,θ) to median, Q₀.₃, Q₀.₁ → CDF p-value

An external background-distance file can be provided (--bg-dist) so that, e.g.,
blastn-window distances are used as the null for classifying minimap/nucmer blocks.

Requires: numpy, scipy
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import optimize, stats

# ═══════════════════════════════════════════════════════════════════════════════
# Column detection
# ═══════════════════════════════════════════════════════════════════════════════

# Each tool schema: canonical names for the columns we care about.
# All other columns are passed through unchanged.
SCHEMAS = {
    # align_util.py block TSV (minimap_blocks, mummer_blocks, blastn_blocks blocks mode)
    "blocks": {
        "contig": "q_contig",
        "start": "q_start",
        "end": "q_end",
        "ssag": "ssag",
        "nident": "nident",
        "identity": "identity",
        "distance": "distance",
    },
    # blastn_blocks.py -W / window_blastn.py window mode
    "windows": {
        "contig": "contig",
        "start": "q_start",
        "end": "q_end",
        "ssag": "ssag",
        "nident": "nident",
        "identity": "identity",
        "distance": "distance",
    },
}


def detect_schema(headers: List[str]) -> dict:
    """Return the schema dict that best matches *headers*, or the first that fits."""
    header_set = set(headers)
    for name, colmap in SCHEMAS.items():
        if all(c in header_set for c in colmap.values()):
            return colmap
    # Fallback: try to find columns heuristically
    fallback: Dict[str, str] = {}
    for key, guesses in {
        "start": ["q_start", "start", "window_start"],
        "end": ["q_end", "end", "window_end"],
        "ssag": ["ssag", "s_contig"],
        "nident": ["nident"],
        "identity": ["identity"],
        "distance": ["distance"],
    }.items():
        for g in guesses:
            if g in header_set:
                fallback[key] = g
                break
    return fallback


# ═══════════════════════════════════════════════════════════════════════════════
# I/O helpers
# ═══════════════════════════════════════════════════════════════════════════════

def load_background_dists(path: str) -> np.ndarray:
    """Load distances from a file (one per line, skip # comments and blanks)."""
    vals: List[float] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                vals.append(float(line))
            except ValueError:
                continue
    if not vals:
        raise SystemExit(f"no valid distances found in {path}")
    return np.array(vals, dtype=np.float64)


def read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    """Read a TSV file. Returns (headers, rows) where each row is a dict."""
    headers: List[str] = []
    rows: List[Dict[str, str]] = []
    with open(path) as fh:
        for lineno, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if lineno == 0:
                # Detect if first line is a header or data
                # Heuristic: if first field looks like a header (not numeric, no _NODE)
                if fields[0].startswith("#"):
                    # comment line — find the real header later
                    while fields[0].startswith("#"):
                        lineno += 1
                        line = next(fh, "").rstrip("\n")
                        if not line:
                            break
                        fields = line.split("\t")
                if fields and not fields[0].startswith("#"):
                    headers = fields
                continue
            if line.startswith("#"):
                continue
            if len(fields) == len(headers):
                rows.append(dict(zip(headers, fields)))
    if not headers:
        raise SystemExit(f"could not find header in {path}")
    return headers, rows


# ═══════════════════════════════════════════════════════════════════════════════
# Statistical tests
# ═══════════════════════════════════════════════════════════════════════════════

def robust_z_test(
    distances: np.ndarray,
    bg_dists: np.ndarray,
) -> Tuple[float, float, np.ndarray]:
    """Compute robust-z p-values.

    Returns (median, mad_scaled, p_values) where p_values[i] = Φ((dᵢ - μ)/σ)
    (one-tailed lower: tests for unusually LOW distance).
    """
    median = float(np.median(bg_dists))
    mad = float(np.median(np.abs(bg_dists - median)))
    mad_scaled = mad * 1.4826
    if mad_scaled < 1e-12:
        # degenerate: all background distances identical
        pvals = np.where(distances < median - 1e-12, 0.0,
                         np.where(distances > median + 1e-12, 1.0, 0.5))
        return median, mad_scaled, pvals
    z = (distances - median) / mad_scaled
    pvals = stats.norm.cdf(z)  # one-tailed lower
    return median, mad_scaled, pvals


def _gamma_quantile_loss(
    params: np.ndarray,
    target_quantiles: np.ndarray,
    probs: np.ndarray,
) -> float:
    """Sum of squared differences between Gamma quantiles and targets."""
    k, theta = params
    if k <= 0 or theta <= 0:
        return 1e12
    fitted = stats.gamma.ppf(probs, a=k, scale=theta)
    return float(np.sum((fitted - target_quantiles) ** 2))


def gamma_quantiles_test(
    distances: np.ndarray,
    bg_dists: np.ndarray,
) -> Tuple[float, float, np.ndarray]:
    """Fit Gamma(k,θ) to bg_dists quantiles (0.1, 0.3, 0.5).

    Returns (shape, scale, p_values) where p_values[i] = Gamma_CDF(dᵢ).
    """
    probs = np.array([0.1, 0.3, 0.5])
    target_q = np.quantile(bg_dists, probs)

    # Initial guess via method of moments
    m = float(np.mean(bg_dists))
    v = float(np.var(bg_dists))
    if v < 1e-12:
        v = 1e-12
    k0 = m ** 2 / v
    theta0 = v / m

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = optimize.minimize(
            _gamma_quantile_loss,
            x0=np.array([k0, theta0]),
            args=(target_q, probs),
            bounds=[(1e-6, None), (1e-6, None)],
            method="L-BFGS-B",
        )

    if result.success:
        k, theta = result.x
    else:
        # Fallback: use method-of-moments
        print(
            f"  note: gamma fit did not converge; using method-of-moments "
            f"(k={k0:.3f}, theta={theta0:.4f})",
            file=sys.stderr,
        )
        k, theta = k0, theta0

    if k <= 0:
        k = 0.5
    if theta <= 0:
        theta = np.mean(bg_dists) / k

    pvals = stats.gamma.cdf(distances, a=k, scale=theta)
    return float(k), float(theta), pvals


# ═══════════════════════════════════════════════════════════════════════════════
# Benjamini-Hochberg
# ═══════════════════════════════════════════════════════════════════════════════

def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Return BH FDR q-values.  Unmapped entries (p=1) are handled naturally."""
    n = len(pvals)
    order = np.argsort(pvals)
    qvals = np.ones(n)
    min_q = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        q = min(1.0, pvals[idx] * n / (i + 1))
        if q < min_q:
            min_q = q
        qvals[idx] = min_q
    return qvals


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Classify genomic intervals from pairwise mapping tool output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, type=Path, help="Input TSV")
    p.add_argument("-o", "--output", required=True, type=Path, help="Output TSV")
    p.add_argument(
        "--bg-dist", type=Path, default=None,
        help="External file of background distances (one per line)",
    )
    p.add_argument(
        "--method", choices=("robust-z", "gamma-quantiles"),
        default="robust-z",
        help="Statistical test method",
    )
    p.add_argument(
        "--alpha", type=float, default=0.05,
        help="FDR threshold for 'significant' category",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if not args.input.is_file():
        raise SystemExit(f"input not found: {args.input}")

    # ── Read input TSV ──
    print(f"reading {args.input}", file=sys.stderr)
    headers, rows = read_tsv(str(args.input))

    # ── Detect schema ──
    col = detect_schema(headers)
    if not col or "ssag" not in col or "nident" not in col:
        # Try harder: look for distance or identity
        pass
    if not col:
        print(
            "warning: could not detect known schema; "
            "expect columns: ssag/s_contig, nident, distance (or identity)",
            file=sys.stderr,
        )
        # Try to work with whatever we have
        col = {}
        for h in headers:
            if "dist" in h.lower():
                col["distance"] = h
            elif "ident" in h.lower():
                col.setdefault("identity", h)
            elif "nident" in h.lower():
                col["nident"] = h
            elif h in ("ssag", "s_contig"):
                col["ssag"] = h

    # Required columns
    ssag_col = col.get("ssag", "ssag")
    nident_col = col.get("nident", "nident")
    distance_col = col.get("distance", "distance")
    identity_col = col.get("identity", "identity")

    have_distance = distance_col in headers
    have_identity = identity_col in headers

    # ── Build per-group distance arrays ──
    # Group by ssag, extract distance for each row, mark unmapped
    groups: Dict[str, List[int]] = {}  # ssag -> list of row indices
    unmapped: List[bool] = [False] * len(rows)
    distances: List[float] = [1.0] * len(rows)

    for i, row in enumerate(rows):
        ssag = row.get(ssag_col, ".")
        groups.setdefault(ssag, []).append(i)

        nident_str = row.get(nident_col, "")
        try:
            nident = int(nident_str)
        except (ValueError, TypeError):
            nident = 0

        if nident <= 0 or (have_distance and float(row.get(distance_col, "1")) >= 1.0):
            unmapped[i] = True
            distances[i] = 1.0
        elif have_distance:
            try:
                distances[i] = float(row[distance_col])
            except (ValueError, KeyError):
                distances[i] = 1.0
                unmapped[i] = True
        elif have_identity:
            try:
                ident = float(row[identity_col])
                distances[i] = 1.0 - ident / 100.0
            except (ValueError, KeyError):
                distances[i] = 1.0
                unmapped[i] = True
        else:
            distances[i] = 1.0
            unmapped[i] = True

    n_unmapped = sum(unmapped)
    print(
        f"  {len(rows)} rows, {n_unmapped} unmapped "
        f"({100 * n_unmapped / max(1, len(rows)):.1f}%)",
        file=sys.stderr,
    )
    print(f"  {len(groups)} subject group(s)", file=sys.stderr)

    # ── Load or extract background distances ──
    if args.bg_dist:
        print(f"loading background distances from {args.bg_dist}", file=sys.stderr)
        bg_all = load_background_dists(str(args.bg_dist))
        print(f"  {len(bg_all)} distances", file=sys.stderr)
        # Use same background for all groups
        use_global_bg = True
    else:
        use_global_bg = False

    # ── Per-group classification ──
    p_values = np.ones(len(rows))
    q_values = np.ones(len(rows))
    test_stats = np.full(len(rows), np.nan)
    categories: List[str] = ["unmapped" if u else "background" for u in unmapped]
    directions: List[str] = ["." for _ in rows]

    out_headers = list(headers) + ["test_stat", "p_value", "q_value", "category", "direction"]

    for ssag, indices in groups.items():
        group_dists = np.array([distances[i] for i in indices])
        group_unmapped = np.array([unmapped[i] for i in indices])
        mapped_mask = ~group_unmapped

        if not np.any(mapped_mask):
            # all unmapped — everything stays at defaults
            continue

        mapped_dists = group_dists[mapped_mask]

        if use_global_bg:
            bg_dists = bg_all
        else:
            bg_dists = mapped_dists

        # Remove zero-distance self-matches from background
        bg_dists = bg_dists[bg_dists > 1e-9]
        if len(bg_dists) < 2:
            print(
                f"  [{ssag}] only {len(bg_dists)} non-zero bg distances; "
                f"all mapped → background",
                file=sys.stderr,
            )
            continue

        print(
            f"  [{ssag}] {np.sum(mapped_mask)} mapped, "
            f"{len(bg_dists)} bg distances, "
            f"method={args.method}",
            file=sys.stderr,
        )

        # ── Run statistical test ──
        if args.method == "gamma-quantiles":
            shape, scale, pv = gamma_quantiles_test(mapped_dists, bg_dists)
            ts = stats.gamma.cdf(mapped_dists, a=shape, scale=scale)
            print(
                f"         gamma k={shape:.4f}  theta={scale:.4f}",
                file=sys.stderr,
            )
        else:
            median, mad_s, pv = robust_z_test(mapped_dists, bg_dists)
            ts = (mapped_dists - median) / max(mad_s, 1e-12)
            print(
                f"         median={median:.4f}  MAD={mad_s:.4f}",
                file=sys.stderr,
            )

        # ── Q-values (Benjamini-Hochberg) ──
        qv = benjamini_hochberg(pv)

        # ── Write back per-row values ──
        mapped_indices = np.array(indices)[mapped_mask]
        for j, mi in enumerate(mapped_indices):
            p_values[mi] = float(pv[j])
            q_values[mi] = float(qv[j])
            test_stats[mi] = float(ts[j])

            if qv[j] <= args.alpha:
                categories[mi] = "significant"
                # Use the full (per-group or global) background median for direction
                bg_median = float(np.median(bg_dists))
                directions[mi] = (
                    "more_similar" if distances[mi] < bg_median else "more_dissimilar"
                )
            else:
                categories[mi] = "background"
                directions[mi] = "."

    # ── Write output ──
    print(f"writing {args.output}", file=sys.stderr)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write("\t".join(out_headers) + "\n")
        for i, row in enumerate(rows):
            vals = [row.get(h, ".") for h in headers]
            ts_str = f"{test_stats[i]:.6g}" if not np.isnan(test_stats[i]) else "."
            vals.extend([
                ts_str,
                f"{p_values[i]:.6g}",
                f"{q_values[i]:.6g}",
                categories[i],
                directions[i],
            ])
            fh.write("\t".join(vals) + "\n")

    # ── Summary ──
    n_sig = sum(1 for c in categories if c == "significant")
    n_bg = sum(1 for c in categories if c == "background")
    n_more_sim = sum(1 for d in directions if d == "more_similar")
    print(
        f"  result: {n_unmapped} unmapped, {n_bg} background, "
        f"{n_sig} significant ({n_more_sim} more_similar)",
        file=sys.stderr,
    )
    print("done.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
