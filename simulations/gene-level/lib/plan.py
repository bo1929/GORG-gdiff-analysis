#!/usr/bin/env python3
"""The design plan, and the file layout it implies.

Stage 01 writes ``plan.tsv`` -- one row per pair in the design::

    seed    alpha   level   variant

Stages 02, 03 and 04 each accept either one case (``--variant``) or ``--plan`` to run the whole
design, so a single stage can be re-run over every pair without redoing the stages before it.

Keeping the paths here means the stages agree on the layout by construction.
"""

from __future__ import annotations

import os
import re

COLS = ["seed", "alpha", "level", "variant"]


def read_plan(path):
    """-> [dict(seed=..., alpha=..., level=..., variant=...)]."""
    rows = []
    head = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if head is None:
                head = f
                continue
            rows.append(dict(zip(head, f)))
    missing = [c for c in COLS if head and c not in head]
    if missing:
        raise SystemExit(f"{path}: plan is missing column(s) {', '.join(missing)}")
    return rows


def write_plan(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")


def genome_paths(simdir, seed, alpha, variant):
    """Inputs that live in the simulation tree."""
    sd = os.path.join(simdir, "genomes", seed)
    vd = os.path.join(sd, variant)
    return {
        "variant_dir": vd,
        "contigs": os.path.join(sd, f"{seed}_contigs.fasta"),
        "genes_faa": os.path.join(vd, "mutated_genes.faa"),
        "query_fa": os.path.join(vd, "mutated.fasta"),
        "weights": os.path.join(sd, f"g_weights_a{alpha}.npy"),
        "gnd_file": os.path.join(vd, "gnd_aad.txt"),
    }


def pair_paths(outdir, variant):
    """Inputs and outputs that live in the results tree."""
    d = os.path.join(outdir, variant)
    return {
        "dir": d,
        "ref_fa": os.path.join(d, "ref.fa"),
        "qry_fa": os.path.join(d, "qry.fa"),
        "qry_fai": os.path.join(d, "qry.fa.fai"),
        "aln": os.path.join(d, "aln.paf"),
        "sketch": os.path.join(d, "ref.gs"),
        "roll": os.path.join(d, "roll.tsv"),
        "resources": os.path.join(d, "resources.tsv"),
        "genes": os.path.join(d, "genes.tsv"),
        "joined": os.path.join(d, "genes-joined.tsv"),
        "genome": os.path.join(d, "genome.tsv"),
    }


def read_gnd(path):
    """The gnd_aad.txt label, or nan. Compared against the realised value, never trusted."""
    try:
        with open(path) as fh:
            return float(fh.read().split()[0])
    except (OSError, ValueError, IndexError):
        return float("nan")


def split_list(s):
    """'a,b c' -> ['a', 'b', 'c']."""
    return [x for x in s.replace(",", " ").split() if x]


VARIANT_RE = re.compile(r"^mutated_BLOSUM62_(.+)_a([0-9]+)_gnd[0-9]+$")


def parse_variant(variant):
    """'mutated_BLOSUM62_AG-359-G18_a22_gnd002' -> ('AG-359-G18', '22'), or (None, None)."""
    m = VARIANT_RE.match(variant)
    return (m.group(1), m.group(2)) if m else (None, None)
