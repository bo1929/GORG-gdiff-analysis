# gene-level

Benchmarks **gdiff** and **wfmash** against the exact simulated divergence of every gene in
`simulations/genomes`. The simulator stored a per-position mutation rate, so each gene's true
divergence is known exactly and no aligner-derived truth is required.

## Usage

Six stages, executed in order. Each defaults to the complete design; `--variant` restricts a
per-pair stage to one pair.

```bash
./01-pick-variant.py --outdir results       # plan.tsv, .realised.tsv
./02-run-tools.sh    --outdir results -t 1  # <variant>/: qry.fa, aln.paf, ref.gs, roll.tsv, resources.tsv
./03-gene-truth.py   --outdir results       # <variant>/genes.tsv
./04-join.py         --outdir results       # <variant>/genes-joined.tsv, genome.tsv
./05-collect.py      --outdir results       # genome-all.tsv, resources-all.tsv, genes-joined-all.tsv
./06-evaluate.py     --outdir results       # metrics.tsv
```

Defaults: every seed in `genomes/` (29), alphas `5 22`, levels `0.001 0.01 0.05 0.10 0.18`, outdir
`./results`. The parameters `-k/--kmer`, `-l/--window`, `-s/--step`, `-p/--min-identity` and
`-t/--threads` must agree between stages; stage 03's `--prefix` must equal stage 02's.

`-t` reaches both tools: `wfmash -t` and `gdiff --num-threads` (a global gdiff option, which
precedes the subcommand). gdiff is deterministic at any thread count; wfmash is not above one
thread, where repeats of the same pair differ. Reproducible results require `-t 1`.

## Definitions

| Column | Definition |
|---|---|
| `truth` | Mean simulated per-position divergence over the gene span. Exact. |
| `d_gd` | Unweighted mean of `d` over the roll windows contained entirely within the gene. |
| `wf_rate` | Mismatches / (matches + mismatches) over the gene span; insertions in neither term. |
| `wf_cov` | Aligned bp / gene length, so `wf_rate` is conditional on alignment. |

`-l 333 -k 23` is a 355 bp window at a 100 bp step; the median gene holds 6 windows. As adjacent
windows overlap by 255 bp and only contained windows are used, `d_gd` weights the centre of a gene
above its ends, spans a median 89.8% of it, and `gd_nwin` is not a sample size.

Every gene is written, scored or not. Two length rules exclude 15.04% of genes and 4.78% of gene bp
at every level; `gdiff_nan` occurs only above 10% divergence and removes high-divergence genes, so
the bias at 18% is a lower bound. Genes are pooled by target `level`, never by realised `gnd`.

## Results

29 seeds × 2 alphas × 5 levels = 290 pairs, `-t 1`, coverage gate 0.5; alpha 22 shown.

| Level | gdiff bias | wfmash bias | gdiff mae | wfmash mae |
|---|---|---|---|---|
| 0.001 | +0.000001 | +0.000001 | 0.000956 | 0.000824 |
| 0.05 | +0.002486 | −0.000006 | 0.009128 | 0.006130 |
| 0.18 | +0.011997 | −0.000761 | 0.021442 | 0.011309 |

gdiff overestimates increasingly with divergence, from +0.1% relative at 0.001 to +6.3% at 0.18.
wfmash is unbiased to within 6e-5 up to 5% divergence and attains the lower mae, rmse and Spearman
correlation at every level. Genome-wide wfmash coverage is 96.8% up to 5% divergence and 87.0% at
18%. The full table is in `results/metrics.tsv`.

## Resource benchmark

`resource-scaling.sh` times one shared index against `1x1` (baseline), `1xN` (query scaling) and
`Nx1` (reference scaling), each scaled case comprising N cells. `-n` sets the number of genomes and
cycles the seed pool if exceeded; `--level` assigns one divergence to every genome, whereas
`--levels` draws one per genome and `--rng-seed` makes that draw reproducible. Output is exactly two
files, the results TSV (`-o`) and the log (`--log`); intermediates are removed on exit.

```bash
./resource-scaling.sh -n 100 --levels "0.001,0.01,0.05,0.10,0.18" -o bench.tsv
```

## Requirements

`numpy`; `GDIFF` (default `../../bin/gdiff`) and `WFMASH` (default `bin/wfmash`), overridable by
environment variable. `lib/make-fai.py` writes the `.fai` files wfmash requires. The genome
directories are not written to.
