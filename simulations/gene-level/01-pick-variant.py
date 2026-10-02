#!/usr/bin/env python3
"""STAGE 01 -- choose the variant matching a target divergence, and write the design plan.

  01-pick-variant.py                                   # whole design -> OUTDIR/plan.tsv
  01-pick-variant.py --seed S --alpha A --level 0.01   # one cell     -> variant on stdout

Selection is by MEASURED divergence -- substitutions/bp between the seed and mutated FASTAs --
not by the gnd_aad.txt label: at least one label is stale, and picking on it ran that seed at a
third of its target. Measurements are cached, so re-running the design is instant.

OUT  plan.tsv columns: seed, alpha, level, variant.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import plan  # noqa: E402


def read_seq(path):
    with open(path, "rb") as fh:
        return b"".join(l.strip() for l in fh if not l.startswith(b">"))


def realised(seed_dir, seed, variant):
    """Substitutions/bp between the seed and mutated genome. The simulation is substitution-only,
    so the two are the same length and a positional compare is exact."""
    ref = read_seq(os.path.join(seed_dir, f"{seed}_contigs.fasta"))
    qry = read_seq(os.path.join(seed_dir, variant, "mutated.fasta"))
    if len(ref) != len(qry):
        raise SystemExit(
            f"{variant}: length mismatch {len(ref)} vs {len(qry)}; the simulation is " f"substitution-only, so this variant cannot be scored this way"
        )
    return sum(a != b for a, b in zip(ref, qry)) / len(ref)


def load_cache(path):
    if not path or not os.path.isfile(path):
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 3 and f[0] != "variant":
                out[f[0]] = (f[1], float(f[2]))
    return out


def save_cache(path, cache):
    tmp = path + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(tmp, "w") as fh:
        fh.write("variant\tlabel\trealised\n")
        for v, (lab, real) in sorted(cache.items()):
            fh.write(f"{v}\t{lab}\t{real:.10g}\n")
    os.replace(tmp, path)


def candidates(simdir, seed, alpha, cache):
    """[(realised, label, variant)] for every candidate, filling cache. -> (rows, cache_fresh)."""
    seed_dir = os.path.join(simdir, "genomes", seed)
    pat = os.path.join(seed_dir, f"mutated_BLOSUM62_{seed}_a{alpha}_gnd*")
    variants = sorted(os.path.basename(d) for d in glob.glob(pat) if os.path.isfile(os.path.join(d, "mutated.fasta")))
    if not variants:
        raise SystemExit(f"no variants for seed={seed} alpha={alpha}")
    fresh = False
    rows = []
    for v in variants:
        ga = os.path.join(seed_dir, v, "gnd_aad.txt")
        label = ""
        if os.path.isfile(ga):
            with open(ga) as fh:
                label = fh.read().split()[0]
        if v in cache:
            label = label or cache[v][0]
            real = cache[v][1]
        else:
            real = realised(seed_dir, seed, v)
            cache[v] = (label, real)
            fresh = True
        rows.append((real, label, v))
    return rows, fresh


def choose(rows, target, tol):
    """-> (variant, realised, warning). Picks the measured divergence nearest the target."""
    real, label, best = min(rows, key=lambda r: abs(r[0] - target))
    warn = ""
    try:
        if label and abs(real / float(label) - 1) > tol:
            warn = f"   ** labelled {float(label):.5f} but realises {real:.5f}: the label is " f"stale, selection used the measured value"
    except ValueError:
        pass
    if abs(real - target) > 0.25 * target:
        warn += f"   ** nearest variant realises {real:.5f}, " f"{100 * abs(real / target - 1):.0f}% off the {target} target"
    return best, real, warn


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--simdir", default=os.path.dirname(HERE),
                   help="the simulations/ tree [this script's parent]")
    p.add_argument("--cache", default=None,
                   help="TSV of variant/label/realised to reuse "
                        "[OUTDIR/.realised.tsv when planning a design]")
    p.add_argument("--gnd-tol", type=float, default=0.01, help="warn when a candidate's label and sequences differ by this fraction")
    one = p.add_argument_group("one cell")
    one.add_argument("--seed")
    one.add_argument("--alpha")
    one.add_argument("--level", type=float, help="target divergence, e.g. 0.01")
    many = p.add_argument_group("a design")
    many.add_argument("--seeds", help="comma- or space-separated [every seed in genomes/]")
    many.add_argument("--alphas", default="5 22")
    many.add_argument("--levels", default="0.001 0.01 0.05 0.10 0.18")
    many.add_argument("--outdir", default="./results", help="writes OUTDIR/plan.tsv")
    args = p.parse_args()

    # No --seed/--level means: the whole design, with defaults. This is the "just run it" mode.
    if not args.seed and not args.seeds:
        gens = os.path.join(args.simdir, "genomes")
        args.seeds = " ".join(sorted(d for d in os.listdir(gens)
                                     if os.path.isdir(os.path.join(gens, d))))
        if not args.seeds:
            raise SystemExit(f"no seed directories under {gens}")

    # The cache defaults into the outdir, so a bare run caches and a re-run is instant.
    # Without it every run re-measures all ~1800 variants, which takes about a minute and a half.
    if args.seeds and not args.cache:
        args.cache = os.path.join(args.outdir, ".realised.tsv")

    cache = load_cache(args.cache)
    dirty = False

    if args.seeds:
        rows = []
        for seed in plan.split_list(args.seeds):
            for alpha in plan.split_list(args.alphas):
                cands, fresh = candidates(args.simdir, seed, alpha, cache)
                dirty |= fresh
                for level in plan.split_list(args.levels):
                    variant, real, warn = choose(cands, float(level), args.gnd_tol)
                    print(f"  {seed} a{alpha} level {level} -> {variant} " f"(realised {real:.5f}){warn}", file=sys.stderr)
                    rows.append(dict(seed=seed, alpha=alpha, level=level, variant=variant))
        if dirty and args.cache:
            save_cache(args.cache, cache)
        out = os.path.join(args.outdir, "plan.tsv")
        os.makedirs(args.outdir, exist_ok=True)
        plan.write_plan(out, rows)
        print(f"wrote {out} ({len(rows)} pairs)")
        return 0

    if not (args.seed and args.alpha and args.level is not None):
        raise SystemExit("give either --seed/--alpha/--level, or --seeds/--alphas/--levels/--outdir")
    cands, fresh = candidates(args.simdir, args.seed, args.alpha, cache)
    if fresh and args.cache:
        save_cache(args.cache, cache)
    variant, real, warn = choose(cands, args.level, args.gnd_tol)
    print(variant)
    print(f"  target {args.level}  ->  {variant}  realised {real:.5f}{warn}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
