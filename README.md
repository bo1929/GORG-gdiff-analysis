# GORG-gdiff-analysis

## Layout

| Path | What it is |
|---|---|
| `bin/` | Bundled binaries + `gdiff`/`dashing2` arch dispatchers (osx  -  x86-64) |
| `methods/` | One shell wrapper per estimator; uniform `<genome_dir> <pairs.tsv> [outdir]` CLI |
| `pairwise-mapping/` | Regional block extractors: BLASTn, minimap2, nucmer -> common TSV schema |
| `simulations/` | Divergence simulation corpus, per-method runners, and the gene-level study |
| `resource-benchmarking/` | All-vs-all wall time / CPU / peak-RSS benchmark over a genome sample |
| `gene_coordinates/` | Per-gene gdiff distances and the BLASTn-vs-gdiff gene comparison |
| `blastn_blocks/` | BLASTn window blocks, all-vs-all and for selected pairs |
| `results/` | Published ANI tables (`ani-comparison/`), roll windows, and figure PDFs |
| `scripts/` | Small C helpers: sample reconciliation, summarisation, in-sample ANI, per-pair Hartigan dip (see `scripts/README.md`) |
| `contigs-gt80-complete/` | 830 genome assemblies, >80% completeness — the main input set |
| `dataset-GORG/` | Wider GORG source data: 12,715 Prokka tables, GBK/16S/ORF annotations |

## Data and fixtures

- `all_pairs.tsv`: 251,534 query/subject/ANI% rows (comment header).
- `anib-groundtruth.csv`: the same 251,534 pairs with the full ANIb truth
  columns: ANI, two-way alignment coverage and length, gene and ortholog
  counts, AAI, AAD and 16S divergence.
- `selected_pairs.tsv` (124 pairs, tab-separated with a comment header) and
  `selected_genomes-{ref,queries}.txt` (100/101 FASTA paths): the working
  subset used for the gene-level and regional comparisons.
- `contigs-gt80-complete/*.fasta`: genome FASTA keyed by strain id; the 828
  ids that also appear in `dataset-GORG` are the annotated core.
- `dataset-GORG/`: `tbl/` (12,715 `<id>_prokka-swissprot.tsv.xz`), `gbk/`
  (179 curated GenBank), `panspecies_gbk/` + `panspecies_blastn/` (828 genomes,
  16S vs ORF hits), `per_hit_mappings/` (197 per-genome hit tables).
- `simulations/genomes/<seed>/`: 29 seed genomes, each with the baseline
  contigs, gene calls, per-gene weights, and `mutated_BLOSUM62_*_a<alpha>_gnd*`
  variant directories at five divergence levels.
- `simulations/metadata.tsv`, `pairs.tsv`: realised GND/AAD/true ANI per
  variant and the variant - baseline pair list.

## Engines and wrappers

`methods/` runs each estimator over an explicit pair list and writes
`distances/<method>-<cfg>.tsv` in a shared schema
(`method, param_setup, genome_a, genome_b, …`): `gdiff-dist.sh`, `mash.sh`,
`skani.sh`, `dashing2.sh`, `fastani.sh`, `anib.sh` (pyani-plus), plus block-level
`blast-windows.sh`, `blast-genes.sh`, `minimap.sh`, `mummer.sh`. Shared caching
and conventions live in `methods/_lib.sh`.

`simulations/run_*.sh` are the variant-oriented equivalents, pairing every
mutated genome against its baseline; `iterate_variant_pairs.sh` and
`build_pairs.py` generate those pairs.

## Gene-level study

`simulations/gene-level/` is a self-contained six-stage pipeline (`01-pick-variant`
-> `06-evaluate`) with its own README, `lib/`, and bundled `wfmash`; it scores
gdiff roll windows and wfmash against the exact simulated per-gene divergence.
`resource-scaling.sh` times one shared index across `1x1`, `1xN` and `Nx1` cells.

## Analysis scripts

| Script | Role |
|---|---|
| `gorg-ani.R` | Main figure set: gdiff vs mash/skani/dashing2/fastANI vs ANIb truth |
| `plot-explore.R` | Regional gdiff-roll vs BLASTn block overlay for one pair |
| `simulations/eval-simulations.R` | Simulation accuracy/bias/variance figures |
| `simulations/cv_analysis.R` | Per-window coefficient-of-variation analysis |
| `simulations/gene-level/plot-wfmash-gorg.R` | Gene-level gdiff-vs-wfmash and resource plots |
| `gene_coordinates/compare-blastn.R` | gdiff vs BLASTn per-gene distance correlation |
| `resource-benchmarking/plot-resource.R` | Peak memory and running-time PDFs |
| `blastn_from_gbk.py` | GBK genes -> BLASTn -> block TSV (chain mode, sentinel rows) |

Outputs: `results/ani-comparison/*/distances/` (gzipped per-method tables),
`results/gdiff-roll/*.tsv.gz` (per-genome window distances),
`simulations/results/` and `simulations/results-ORFvANI/`, and the `S-`, `G-`,
`R-` prefixed PDFs in `results/`, `simulations/` and `resource-benchmarking/`.

## Running

Python deps: `pip install -r requirements.txt`; non-Python tools expected on
`PATH` are listed at the top of that file (mash, skani, minimap2, BLAST+,
nucmer, fastANI, pyani-plus). `GDIFF`/`WFMASH` are overridable by environment
variable. `commands.md` holds the original sketch/roll invocations. Scripts
resolve paths relative to themselves, so most can be run from their own
directory.
