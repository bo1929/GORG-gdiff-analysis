# gene-level

Benchmarks **gdiff** and **wfmash** against the realised substitution fraction of every gene in
`simulations/genomes`. The simulation is substitution-only, so a gene's true divergence is
counted directly from the reference and query sequences: no aligner and no expected rate is
involved.

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
| `truth` | Realised substitutions / gene span, counted from the reference and query sequences. Exact, and per gene. |
| `d_gd` | Unweighted mean of `d` over the roll windows contained entirely within the gene. |
| `wf_rate` | Mismatches / (matches + mismatches) over the gene span; insertions in neither term. |
| `wf_cov` | Aligned bp / gene length, so `wf_rate` is conditional on alignment. |

`truth` is counted from the sequences, not from the simulator's weights. The weights give the
*expected* rate, and at low divergence the realised count departs from it by more than either
tool's error: a 1 kb gene at D = 0.001 holds about one substitution, so its realised rate is
quantised to 0, 1/len, 2/len, ... `genes.tsv` keeps `w_mean` (the expected rate, rescaled to the
realised genome-wide divergence) for provenance only. Genes with no substitution have
`truth = 0`, where relative-error measures (MAPE and the like) are undefined -- use the absolute
`bias`, `mae` and `rmse` columns.

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
| 0.001 | −0.000001 | +0.000000 | 0.000391 | 0.000000 |
| 0.05 | +0.002492 | −0.000002 | 0.006382 | 0.000002 |
| 0.18 | +0.011636 | −0.000743 | 0.017990 | 0.000861 |

Both columns are the *intersection*: genes with `wf_cov >= 0.5` and a finite gdiff estimate --
the `intersection` row of `metrics.tsv` -- so the two tools are scored on the same genes.

gdiff overestimates increasingly with divergence, from −0.1% relative at 0.001 to +4.5% at 0.05
and +6.4% at 0.18. wfmash now agrees with `truth` almost exactly -- mae 2e-6 at 5%, 0.000861 at
18% -- because both are counted from the same sequences: the simulated substitutions are the only
differences, and `wf_rate` counts them directly. The residual at 18% sits in the ~4% of genes
wfmash cannot align end to end (`wf_cov < 0.99`), where it misses substitutions in the unaligned
part. The informative column is therefore gdiff's. Genome-wide wfmash coverage is 96.8% up to 5%
divergence and 87.0% at 18%. The full table is in `results/metrics.tsv`.

`metrics.tsv` also carries a `coverage_pct` column. The denominator is every gene except the
`shorter_than_window` and `no_window_inside` cases, and a gene a method cannot score counts as a
miss: gdiff is uncovered when `d_gd` is NA (`gdiff_nan`, 718 gene rows, all at 10-18%), wfmash
when `wf_cov < 0.5`. Over the whole design that is 99.8% for gdiff and 95.3% for wfmash; at 0.18
the intersection drops to 88.3% (alpha 22), almost all of it wfmash's coverage, while gdiff stays
at 99.8%. The figures report the same shares in the scatter and above the bars.

## Figures

`plot-wfmash-gorg.R` reads `results/` and writes four PDFs; run it from this directory.

```bash
Rscript plot-wfmash-gorg.R
```

| file | what it shows |
|---|---|
| `S-mape-corr-wfmash.pdf` | per-gene estimate vs realised D, binned means with SE; R and coverage per tool |
| `S-mape-cov-wfmash.pdf` | MAE by gene-D bin for both tools; the labels above the bars are coverage, i.e. the share of eligible genes in that bin the method can score |
| `S-running-time-wfmash-nt8.pdf` | running time at 8 threads |
| `S-peak-memory-wfmash-nt8.pdf` | peak memory at 8 threads |

Both tools are scored on the same genes: `included == 1`, `wf_cov >= 0.5` and finite estimates for
both (the intersection). Coverage treats a gene as a miss when the method has no estimate, over a
denominator that removes only the length-filter and no-window-inside cases. The bar heights are
MAE; the labels above them are that coverage.

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
