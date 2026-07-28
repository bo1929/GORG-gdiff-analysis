# Pairwise regional distance scripts

| Script | Role |
|--------|------|
| [`blastn_blocks.py`](blastn_blocks.py) | BLASTn blocks (optional sliding `-W`) |
| [`minimap_blocks.py`](minimap_blocks.py) | minimap2 blocks |
| [`mummer_blocks.py`](mummer_blocks.py) | nucmer blocks |
| [`run_benchmark.sh`](run_benchmark.sh) | Batch: queries x refs |
| [`align_util.py`](align_util.py) | Shared helpers |
| [`../outlier-detection/classify_intervals.py`](../outlier-detection/classify_intervals.py) | Post-hoc outlier classification |

Edit **CONFIG** at the top of each Python tool. Coords are 1-based inclusive.
`distance = 1 - nident/denominator`. Block TSVs include `ssag`.
`blastn_blocks.py` defaults to ultra-sensitive (`-word_size 7 -evalue 1000`).

## Batch (20 queries x ~800 refs)

```text
gdiff    sketch refs once; map once per query
minimap  one combined refs FASTA; one run per query
blast    same; ultra-sensitive blastn blocks per query
mummer   pairwise in parallel (JOBS)
```

```bash
# queries.txt: one SAG id per line
METHODS=gdiff,minimap,blast,mummer \
SENSITIVE=1 JOBS=8 THREADS=4 \
GDIFF=/path/to/gdiff \
./run_benchmark.sh /path/to/contigs-panspecies queries.txt /path/to/benchmark_out
```

Optional: `REFS=refs.txt` `FORCE=1` `WINDOW=1000`
`BLAST_ARGS` (default `--max-target-seqs 5000`) / `BLAST_WINDOW=1000` (sliding)
`GDIFF_SKETCH_ARGS` / `GDIFF_MAP_ARGS` (default map: `-l $((WINDOW/2)) -d ...`).
Skips non-empty outputs unless `FORCE=1`.

## Single pair

```bash
./mummer_blocks.py  --sensitive -q query.fa -s subject.fa -o mummer.tsv -t 8
./minimap_blocks.py --sensitive -q query.fa -s subject.fa -o minimap.tsv -t 8
./blastn_blocks.py              -q query.fa -s subject.fa -o blastn.tsv -t 8
./blastn_blocks.py -W 1000      -q query.fa -s subject.fa -o windows.tsv -t 8
```

## Outlier detection

For HGT / atypical-region detection (classifying intervals as unmapped / background /
significant), use [`../outlier-detection/classify_intervals.py`](../outlier-detection/classify_intervals.py).
It adds statistical columns (p-value, q-value, category) to any block/window TSV.
