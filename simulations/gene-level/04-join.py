#!/usr/bin/env python3
"""STAGE 04 -- join the per-gene truth with both methods' per-gene measurements.

  04-join.py                        # the whole design
  04-join.py --variant V --level F  # one pair

  truth    mean simulated per-position divergence over the gene span, from stage 03
  d_gd     unweighted mean of the roll windows lying entirely inside the gene
  wf_rate  mismatches / (matches + mismatches) over the gene span, insertions in neither term
  wf_cov   aligned bp / gene length

Every gene is written, included or not, so the excluded portion is reported rather than dropped.
The `gnd` column is the realised divergence from stage 03, not the gnd_aad.txt label.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import plan  # noqa: E402

REALISED_RE = re.compile(r"measured from the sequences:\s*([\d.eE+-]+)")
LABEL_RE = re.compile(r"gnd_aad\.txt:\s*([\d.eE+-]+)")

COLS = [
    "seed",
    "alpha",
    "level",
    "gnd",
    "contig",
    "gene",
    "start",
    "end",
    "strand",
    "gene_len",
    "truth",
    "included",
    "reason",
    "d_gd",
    "gd_nwin",
    "gd_cov",
    "wf_rate",
    "wf_cov",
    "wf_aligned_bp",
]


def parse_cigar(cg):
    ops = []
    num = ""
    for ch in cg:
        if ch.isdigit():
            num += ch
            continue
        n = int(num or 0)
        num = ""
        if ch in "M=X":
            # wfmash emits only =/X/I/D, but M is an ambiguous match that may contain
            # mismatches, so treating it as a match would understate the rate. Fail loudly.
            if ch == "M":
                raise SystemExit("CIGAR uses 'M'; its mismatches are not recoverable -- " "re-run wfmash so the cg tag uses '=' and 'X'")
            ops.append(("match" if ch == "=" else "mismatch", n))
        elif ch == "I":
            ops.append(("ins", n))
        elif ch in "DN":
            ops.append(("del", 0))
        elif ch not in "SH":
            raise SystemExit(f"unrecognised CIGAR operator {ch!r}")
    return ops


def parse_paf(path):
    """{query_name: {'qlen': int, 'alns': [(qstart, qend, strand, ops)]}}."""
    out = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise SystemExit(f"{path}:{lineno}: expected >= 12 PAF columns")
            tags = {t[:2]: t[5:] for t in f[12:] if len(t) >= 5 and t[2] == ":"}
            if "cg" not in tags:
                raise SystemExit(f"{path}:{lineno}: no cg tag")
            rec = out.setdefault(f[0], {"qlen": int(f[1]), "alns": []})
            rec["alns"].append((int(f[2]), int(f[3]), f[4], parse_cigar(tags["cg"])))
    return out


def build_track(alns, qlen):
    """Per-position flags over the query: aligned at all, and matching."""
    is_aligned = np.zeros(qlen, dtype=np.uint8)
    is_match = np.zeros(qlen, dtype=np.uint8)
    for qs, qe, strand, ops in alns:
        qpos = qs if strand == "+" else qe
        for kind, qspan in ops:
            if qspan == 0:
                continue
            if strand == "+":
                a, b = qpos, qpos + qspan
                qpos += qspan
            else:
                a, b = qpos - qspan, qpos
                qpos -= qspan
            a, b = max(a, 0), min(b, qlen)
            if a < b and kind != "ins":
                is_aligned[a:b] = 1
                if kind == "match":
                    is_match[a:b] = 1
    return is_aligned, is_match


def num(tok):
    return float("nan") if tok in ("", "NA", "nan") else float(tok)


def fmt(x, spec=".10g"):
    return "NA" if isinstance(x, float) and math.isnan(x) else format(x, spec)


def read_genes_tsv(path):
    """-> (rows, realised, label). The '#' lines carry the divergence stage 03 measured."""
    genes, realised, label, head = [], None, None, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                m = REALISED_RE.search(line)
                if m:
                    realised = float(m.group(1))
                m = LABEL_RE.search(line)
                if m:
                    label = float(m.group(1))
                continue
            f = line.rstrip("\n").split("\t")
            if head is None:
                head = f
                continue
            if f and f[0]:
                genes.append(dict(zip(head, f)))
    if realised is None:
        raise SystemExit(f"{path}: no 'realised divergence' comment line; regenerate it with " f"stage 03 rather than supplying the value by hand")
    return genes, realised, label


def run_one(outdir, variant, seed, alpha, level, k, l):
    p = plan.pair_paths(outdir, variant)
    for key in ("genes", "aln", "roll"):
        if not os.path.isfile(p[key]):
            raise SystemExit(f"error: no {key} at {p[key]}")
    win_bp = l + k - 1

    genes, realised, label = read_genes_tsv(p["genes"])
    gnd = f"{realised:.10g}"
    if label is not None and label > 0 and abs(realised / label - 1) > 0.01:
        print(f"    note: gnd_aad.txt label {label:.6g} differs from the realised " f"{realised:.6g}; reporting the realised value")

    paf = parse_paf(p["aln"])
    roll = {}
    with open(p["roll"]) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if not f or f[0] == "seq":
                continue
            val = f[5] if len(f) == 6 and f[3] == "." else f[4] if len(f) > 4 else "NA"
            roll.setdefault(f[0], {})[int(f[1])] = num(val)

    tracks = {}

    def track(contig):
        if contig not in tracks:
            rec = paf.get(contig)
            if rec is None:
                tracks[contig] = None
            else:
                ia, im = build_track(rec["alns"], rec["qlen"])
                tracks[contig] = (
                    np.concatenate(([0], np.cumsum(ia, dtype=np.int64))),
                    np.concatenate(([0], np.cumsum(im, dtype=np.int64))),
                    rec["qlen"],
                )
        return tracks[contig]

    n_in = n_out = bp_in = bp_out = 0
    reasons = {}
    with open(p["joined"], "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for g in genes:
            contig = g["contig"]
            start, end = int(g["start"]), int(g["end"])
            glen = int(g["gene_len"])
            truth = float(g["truth"])

            # roll windows lying entirely inside the gene
            wins = sorted(s for s in roll.get(contig, {}) if start <= s and s + win_bp - 1 <= end)
            ds = [d for d in (roll[contig][s] for s in wins) if not math.isnan(d)]

            reason = ""
            if glen < win_bp:
                reason = "shorter_than_window"
            elif not wins:
                reason = "no_window_inside"
            elif not ds:
                reason = "gdiff_nan"

            if reason:
                d_gd, gd_cov, nwin, included = float("nan"), 0.0, len(wins), 0
                n_out += 1
                bp_out += glen
                reasons[reason] = reasons.get(reason, 0) + 1
            else:
                d_gd = sum(ds) / len(ds)
                gd_cov = min(1.0, (wins[-1] + win_bp - 1 - wins[0] + 1) / glen)
                nwin, included = len(wins), 1
                n_in += 1
                bp_in += glen

            tr = track(contig)
            if tr is None:
                wf_rate, wf_cov, wf_bp = float("nan"), 0.0, 0
            else:
                ac, mc, _ = tr
                b1, b2 = start - 1, min(end, len(ac) - 1)
                abp = int(ac[b2] - ac[b1])
                mbp = int(mc[b2] - mc[b1])
                wf_rate = (abp - mbp) / abp if abp else float("nan")
                wf_cov = abp / (b2 - b1) if b2 > b1 else 0.0
                wf_bp = abp

            fh.write(
                "\t".join(
                    [
                        seed,
                        alpha,
                        level or "-",
                        gnd,
                        contig,
                        g["gene"],
                        str(start),
                        str(end),
                        g["strand"],
                        str(glen),
                        fmt(truth),
                        str(included),
                        reason or ".",
                        fmt(d_gd),
                        str(nwin),
                        f"{gd_cov:.6f}",
                        fmt(wf_rate),
                        f"{wf_cov:.6f}",
                        str(wf_bp),
                    ]
                )
                + "\n"
            )

    tot, tot_bp = n_in + n_out, bp_in + bp_out
    total_q = sum(int(l_.split("\t")[1]) for l_ in open(p["qry_fai"]) if l_.strip())
    aligned_q = sum(int(t[0][-1]) for t in tracks.values() if t is not None)
    with open(p["genome"], "w") as fh:
        fh.write(
            "seed\talpha\tlevel\tgnd\tquery_bp\taligned_bp\twfmash_cov\tn_alignments\t" "genes_total\tgenes_included\tgene_bp\tgene_bp_included\n"
        )
        fh.write(
            f"{seed}\t{alpha}\t{level or '-'}\t{gnd}\t{total_q}\t{aligned_q}\t"
            f"{aligned_q / total_q:.6f}\t{sum(len(r['alns']) for r in paf.values())}\t"
            f"{tot}\t{n_in}\t{tot_bp}\t{bp_in}\n"
        )

    print(f"--- {variant}")
    print(
        f"    genes {tot} ({tot_bp:,} bp): {n_in} included ({100 * bp_in / tot_bp:.1f}% of "
        f"gene bp), {n_out} excluded"
        + ("  [" + ", ".join(f"{r} {c}" for r, c in sorted(reasons.items(), key=lambda kv: -kv[1])) + "]" if reasons else "")
    )
    print(f"    genome cov {100 * aligned_q / total_q:.2f}%   wrote {os.path.basename(p['joined'])}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variant")
    p.add_argument("--plan", help="run every row of the design")
    p.add_argument("--outdir", default="./results")
    p.add_argument("--seed", default="")
    p.add_argument("--alpha", default="")
    p.add_argument("--level", default="", help="target divergence; the stage-06 pooling key")
    p.add_argument("-k", "--kmer", type=int, default=23, help="k-mer length [23]")
    p.add_argument("-l", "--window", type=int, default=333, help="roll window in k-mers [333]")
    args = p.parse_args()

    # Default job is the whole design; --variant switches to a single pair.
    if not args.variant and not args.plan:
        args.plan = os.path.join(args.outdir, "plan.tsv")

    if args.plan:
        for row in plan.read_plan(args.plan):
            run_one(args.outdir, row["variant"], row["seed"], row["alpha"], row["level"], args.kmer, args.window)
        return 0

    seed, alpha = args.seed, args.alpha
    if not (seed and alpha):
        seed, alpha = plan.parse_variant(args.variant)
        if not seed:
            raise SystemExit(f"cannot parse seed/alpha from {args.variant}; pass --seed and --alpha")
    run_one(args.outdir, args.variant, seed, alpha, args.level, args.kmer, args.window)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
