#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

PRUNED_RE = re.compile(r"^mutated_p(\d+)_s(\d+)\.fasta$")
GND_DIR_RE = re.compile(r"^mutated_BLOSUM62_(?P<seed>[A-Z]+-.+)_a(?P<alpha>\d+)_gnd(?P<val>\d+)$")


def read_gnd_aad(path: Path) -> tuple[float, float] | None:
    if not path.exists():
        return None
    parts = path.read_text().strip().split()
    if len(parts) < 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None


def iter_seeds(genomes_dir: Path):
    return sorted(p for p in genomes_dir.iterdir() if p.is_dir())


def build_rows(genomes_dir: Path):
    for seed_dir in iter_seeds(genomes_dir):
        seed = seed_dir.name
        for variant_dir in sorted(seed_dir.glob("mutated_BLOSUM62_*")):
            m = GND_DIR_RE.match(variant_dir.name)
            if not m:
                continue
            gnd_aad = variant_dir / "gnd_aad.txt"
            gnd_aad_txt = read_gnd_aad(gnd_aad)
            if gnd_aad_txt is None:
                print(f"  WARNING: missing/malformed {gnd_aad}, skipping {variant_dir.name}", file=__import__("sys").stderr)
                continue
            gnd, aad = gnd_aad_txt
            true_ani = 100.0 * (1.0 - gnd)

            suffix = f"_a{m['alpha']}_gnd{m['val']}"
            level = f"{seed}{suffix}"

            yield seed, f"{level}", level, gnd, aad, true_ani, 1, "", ""

            for pruned in sorted(variant_dir.glob("mutated_p*_s*.fasta")):
                pm = PRUNED_RE.match(pruned.name)
                if not pm:
                    continue
                miss, s = int(pm.group(1)), int(pm.group(2))
                key = f"{level}_p{miss}_s{s}"
                yield seed, key, level, gnd, aad, true_ani, 0, str(miss), str(s)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes", type=Path, default=Path("genomes"), help="directory containing seed genome dirs")
    parser.add_argument("--metadata", type=Path, default=Path("metadata.tsv"))
    parser.add_argument("--pairs", type=Path, default=Path("pairs.tsv"))
    args = parser.parse_args()

    rows = list(build_rows(args.genomes))

    with args.metadata.open("w") as fh:
        fh.write("key\tseed\tlevel\tGND\tAAD\ttrue_ani\tcomplete\tmissing_percent\ts\n")
        for seed, key, level, gnd, aad, ani, complete, miss, s in rows:
            gnd_str = f"{gnd:.9f}" if complete == 1 else f"{gnd:.6f}"
            fh.write(f"{key}\t{seed}\t{level}\t{gnd_str}\t{aad:.9f}\t{ani:.6f}\t{complete}\t{miss}\t{s}\n")

    with args.pairs.open("w") as fh:
        for seed, key, level, gnd, aad, ani, complete, miss, s in rows:
            fh.write(f"{seed}\t{key}\t{ani:.6f}\n")

    n_seeds = len(set(r[0] for r in rows))
    n_complete = sum(1 for r in rows if r[6] == 1)
    n_pruned = sum(1 for r in rows if r[6] == 0)
    print(f"Wrote {len(rows)} rows ({n_complete} complete, {n_pruned} pruned) across {n_seeds} seeds")
    print(f"  -> {args.metadata}")
    print(f"  -> {args.pairs}")


if __name__ == "__main__":
    main()
