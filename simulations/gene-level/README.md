# gene-level

Benchmark **gdiff** and **wfmash** against the exact simulated divergence of every gene, on the
`simulations/genomes` set. The simulator stored a per-position mutation rate, so each gene's true
divergence is known exactly and no aligner-based truth is needed.

## Pipeline

Six stages, run one at a time. Each defaults to the **whole design**; `--variant` narrows a
per-pair stage to a single pair.

```bash
./01-pick-variant.py --outdir results       # design -> results/plan.tsv
./02-run-tools.sh    --outdir results -t 1  # wfmash + gdiff, timed
./03-gene-truth.py   --outdir results       # exact per-gene truth
./04-join.py         --outdir results       # truth + both methods, per gene
./05-collect.py      --outdir results       # merge the per-pair tables
./06-evaluate.py     --outdir results       # metrics.tsv
```

| stage | writes |
|---|---|
| 01 pick-variant | `plan.tsv`, `.realised.tsv` (measured-divergence cache) |
| 02 run-tools | `<variant>/` — `ref.fa`, `qry.fa`, `*.fai`, `aln.paf`, `ref.gs`, `roll.tsv`, `resources.tsv`, logs |
| 03 gene-truth | `<variant>/genes.tsv` |
| 04 join | `<variant>/genes-joined.tsv`, `<variant>/genome.tsv` |
| 05 collect | `genome-all.tsv`, `resources-all.tsv`, `genes-joined-all.tsv` |
| 06 evaluate | `metrics.tsv` |

Defaults: every seed in `genomes/` (29 of them), alphas `5 22`, levels
`0.001 0.01 0.05 0.10 0.18`, outdir `./results`. Tool parameters must match between stages:
`-k/--kmer`, `-l/--window`, `-s/--step`, `-p/--min-identity`, `-t/--threads` — and stage 03's
`--prefix` must equal stage 02's, or no window lands inside any gene.

`05-collect.py` merges the per-pair tables into three whole-run files: `genome-all.tsv` (one row
per pair), `resources-all.tsv` (one row per timed command) and `genes-joined-all.tsv` (every pair's
per-gene rows — 371,660 rows for the full design, so it is streamed rather than sorted in memory).
All three refuse to combine files whose columns disagree.

**Threads.** `-t` reaches both tools: `wfmash -t` and `gdiff --num-threads` (a *global* gdiff
option, so it precedes the subcommand). They differ in reproducibility, which matters:

| | `-t 1` | `-t >1` |
|---|---|---|
| gdiff | deterministic | deterministic — identical `d` at `-t 1` and `-t 4` |
| wfmash | deterministic | **not reproducible** — repeats of the same pair differ |

wfmash's chaining races above one thread: one pair gave 54/54/54/54/54/54 alignment blocks at
`-t 1` but 55/55/56/55/56/56 at `-t 4`, and two whole-design runs at `-t 4` moved the gene-level
metrics. **Use `-t 1` for results you need to reproduce** — the table below was produced that way.
Raise `-t` only for timing; stage 02 warns when you do.

## Definitions

Per gene, three estimates of the gene's mean per-base substitution rate:

| | |
|---|---|
| `truth` | Exact. Mean simulated per-position divergence over the gene span. |
| `d_gd` | gdiff. Unweighted mean of the `d` of roll windows lying entirely inside the gene. |
| `wf_rate` | wfmash. mismatches / (matches + mismatches) over the gene span; insertions in neither term. |
| `wf_cov` | aligned bp / gene length, so `wf_rate` is conditional on wfmash having aligned it. |

`-l 333 -k 23` is a 355 bp window, `-s 100` starts one per 100 bp; a median gene holds 6 windows.

**Averaging windows is not averaging over the gene.** Windows overlap by 255 bp, so `d_gd` weights
the middle of a gene more than its ends, and only fully-contained windows count, so their union
covers a median 89.8% of it. On this data the effect is small (per seed, ~73% of genes carry a
single simulated rate), but per gene the averaging rule moves the value by up to 2.6e-3 — so it is
part of the definition. Windows are correlated: `gd_nwin` is not a sample size.

**Exclusions.** Every gene is written, scored or not. Two length rules drop 15.04% of genes but
only 4.78% of gene bp, identically at every level. A third, `gdiff_nan`, appears only above 10%
divergence — 651 of 37,166 genes (1.75%) at alpha 5 / 0.18, 64 at alpha 22 / 0.18 — and takes the
*high*-divergence genes, so the bias reported at 18% is a floor.

**Pooling.** By target `level`, never by realised `gnd`: each seed lands on a slightly different
divergence, so pooling by `gnd` would give every seed its own cell.

## Resource benchmark

`resource-scaling.sh` times **one shared index** against 1 and N references × 1 and N queries:

```bash
./resource-scaling.sh --seeds "A,B,..." --level 0.05 -o scaling.tsv
```

It writes exactly two files — the results TSV (`-o`, default `./resource-scaling.tsv`) and a log
(`--log`, default alongside the TSV). FASTAs, indexes and alignments go to a temporary directory
that is removed on exit. The TSV has one row per timed command:
`case  tool  step  references  queries  cells  threads  wall_s  user_s  sys_s  peak_rss_mb  exit`.

A *cell* is one query against one reference. Cost tracks cells, not genomes: 1 → N on a shared
index multiplies the cells by N², so dividing wall time by the number of genomes overstates the
growth. Compare `1x1 → 1xN` (query scaling) with `1x1 → Nx1` (reference scaling) before blaming the
shared index.

## Requirements

`numpy`, plus `GDIFF` (default `../../bin/gdiff`) and `WFMASH` (default `bin/wfmash`), both
overridable by environment variable. `lib/make-fai.py` writes the `.fai` files wfmash needs. The
genome directories are never written to: the seed FASTA is symlinked into the output directory.
