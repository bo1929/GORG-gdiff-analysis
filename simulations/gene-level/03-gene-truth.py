#!/usr/bin/env python3
"""STAGE 03 -- exact per-gene truth from the sequences.

  03-gene-truth.py                 # the whole design
  03-gene-truth.py --variant V     # one pair

truth(gene) is the realised substitution fraction over the whole gene span: the positions where
the query base differs from the reference, divided by the gene length. The simulation is
substitution-only, so reference and query coordinates correspond 1:1 and no alignment is needed.
w_mean (the simulator's weight mean, rescaled to the realised genome-wide rate) is written for
provenance only: it is the *expected* per-gene rate, which fluctuates gene by gene around the
realised count, and at low divergence that fluctuation is large.

--prefix MUST match stage 02's --prefix, or no window will land inside any gene.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import plan  # noqa: E402

GENE_RE = re.compile(r"^>(\S+)\s*#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?\d+)")


def read_contigs(path):
    """[(name, length)] in FASTA order - the order the weight rows follow."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append([line[1:].split()[0], 0])
            elif out:
                out[-1][1] += len(line.strip())
    return [(n, l) for n, l in out]


def read_seq(path):
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


def read_weights(path, contigs):
    with open(path) as fh:
        raw = json.load(fh)
    if len(raw) != len(contigs):
        raise SystemExit(f"{path}: {len(raw)} weight rows but {len(contigs)} contigs")
    rows = []
    for (name, length), row in zip(contigs, raw):
        if len(row) != length:
            raise SystemExit(f"{path}: row for {name} has {len(row)} values, contig is {length} bp")
        rows.append(np.array([0.0 if x is None else float(x) for x in row]))
    return rows


def read_genes(path, contig_names):
    """[(gene, contig, start, end, strand)]; contig matched by longest name prefix."""
    by_len = sorted(contig_names, key=len, reverse=True)
    out = []
    with open(path) as fh:
        for line in fh:
            m = GENE_RE.match(line)
            if not m:
                continue
            gene, s, e, strand = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
            contig = next((c for c in by_len if gene.startswith(c + "_")), None)
            if contig is None:
                raise SystemExit(f"{path}: cannot map gene {gene!r} to a contig")
            out.append((gene, contig, s, e, strand))
    return out


def run_one(simdir, outdir, seed, alpha, variant, k, l, contig_prefix, gnd_tol):
    g = plan.genome_paths(simdir, seed, alpha, variant)
    for key in ("contigs", "genes_faa", "query_fa", "weights"):
        if not os.path.isfile(g[key]):
            raise SystemExit(f"error: no {key} at {g[key]}")
    out = plan.pair_paths(outdir, variant)["genes"]
    os.makedirs(os.path.dirname(out), exist_ok=True)

    contigs = read_contigs(g["contigs"])
    rows = read_weights(g["weights"], contigs)
    names = {n for n, _ in contigs}
    order = {n: i for i, (n, _) in enumerate(contigs)}

    # realised divergence, measured from the sequences
    ref = read_seq(g["contigs"])
    qry = read_seq(g["query_fa"])
    if len(ref) != len(qry):
        raise SystemExit(f"length mismatch: {g['contigs']} {len(ref)} bp vs {g['query_fa']} " f"{len(qry)} bp; the simulation is substitution-only")
    same = np.frombuffer(ref.encode(), dtype=np.uint8) == np.frombuffer(qry.encode(), dtype=np.uint8)
    realised = float(len(same) - np.count_nonzero(same)) / len(ref)

    total_pos = sum(len(r) for r in rows)
    mean_w = sum(float(r.sum()) for r in rows) / total_pos
    scale = realised / mean_w

    label = plan.read_gnd(g["gnd_file"])
    warn = ""
    if label == label and label > 0 and abs(realised / label - 1) > gnd_tol:
        warn = (
            f"   ** gnd_aad.txt says {label:.6g} but the sequences realise {realised:.6g} "
            f"(ratio {realised / label:.3f}) - using the realised value"
        )

    genes = read_genes(g["genes_faa"], names)
    # Per-gene truth is counted from the two sequences at the same coordinates, so the contigs
    # need an offset table (read_seq concatenates them in FASTA order).
    offsets = np.cumsum([0] + [length for _, length in contigs])
    ref_b = np.frombuffer(ref.encode(), dtype=np.uint8)
    qry_b = np.frombuffer(qry.encode(), dtype=np.uint8)
    lens, tr = [], []
    with open(out, "w") as fh:
        fh.write(f"# realised divergence measured from the sequences: {realised:.10g}\n")
        fh.write(f"# gnd_aad.txt: {label:.10g}   scale = realised/mean(w) = {scale:.10g}\n")
        fh.write("# truth = realised substitutions over the gene span, counted from ref vs query\n")
        fh.write("contig\tgene\tstart\tend\tstrand\tgene_len\tw_mean\ttruth\n")
        for gene, contig, s, e, strand in genes:
            i = order[contig]
            base = int(offsets[i])
            b1, b2 = base + s - 1, base + min(e, contigs[i][1])
            if b2 <= b1:
                continue
            w = rows[i][s - 1 : e]
            w_mean = float(w.mean())
            truth = float(np.count_nonzero(ref_b[b1:b2] != qry_b[b1:b2])) / (b2 - b1)
            lens.append(e - s + 1)
            tr.append(truth)
            fh.write(f"{contig_prefix}{contig}\t{gene}\t{s}\t{e}\t{strand}\t" f"{e - s + 1}\t{w_mean:.10g}\t{truth:.10g}\n")

    lens = np.array(lens)
    tr = np.array(tr)
    win_bp = l + k - 1
    print(f"--- {variant}")
    print(f"    realised {realised:.6f}   gnd_aad.txt {label:.6f}{warn}")
    print(f"    genes    {len(lens)}  ({lens.sum():,} bp, " f"{100 * lens.sum() / total_pos:.1f}% of genome)")
    print(f"    length   min={lens.min()} median={np.median(lens):.0f} max={lens.max()}")
    short = lens < win_bp
    print(
        f"    below the {win_bp} bp window: {short.sum()} genes "
        f"({100 * short.mean():.1f}% of genes, {100 * lens[short].sum() / lens.sum():.1f}% "
        f"of gene bp) - these cannot be scored by gdiff"
    )
    print(f"    wrote    {out}")
    return 0


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variant")
    p.add_argument("--plan", help="run every row of the design")
    p.add_argument("--outdir", default="./results")
    p.add_argument("--simdir", default=os.path.dirname(here))
    p.add_argument("--seed")
    p.add_argument("--alpha")
    p.add_argument("-k", "--kmer", type=int, default=23, help="k-mer length [23]")
    p.add_argument("-l", "--window", type=int, default=333, help="roll window in k-mers [333]")
    p.add_argument(
        "--prefix",
        default="Q_",
        help="prefix added to contig names; MUST match stage 02's --prefix, or no " "window will be found inside any gene [Q_]",
    )
    p.add_argument("--gnd-tol", type=float, default=0.01, help="warn when |realised/label - 1| exceeds this [0.01]")
    args = p.parse_args()
    opt = (args.kmer, args.window, args.prefix, args.gnd_tol)

    # Default job is the whole design; --variant switches to a single pair.
    if not args.variant and not args.plan:
        args.plan = os.path.join(args.outdir, "plan.tsv")

    if args.plan:
        for row in plan.read_plan(args.plan):
            run_one(args.simdir, args.outdir, row["seed"], row["alpha"], row["variant"], *opt)
        return 0

    seed, alpha = args.seed, args.alpha
    if not (seed and alpha):
        seed, alpha = plan.parse_variant(args.variant)
        if not seed:
            raise SystemExit(f"cannot parse seed/alpha from {args.variant}; pass --seed and --alpha")
    run_one(args.simdir, args.outdir, seed, alpha, args.variant, *opt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
