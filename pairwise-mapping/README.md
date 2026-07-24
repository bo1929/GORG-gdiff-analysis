# Pairwise regional distance scripts

Query-vs-subject regional similarity for later outlier / distance analysis
(not mapping benchmarks). Subject coords are metadata. Coords are 1-based inclusive.

| Script | Role | Tool |
|--------|------|------|
| [`sliding_blastn.py`](sliding_blastn.py) | Sliding-window BLASTn | blastn |
| [`minimap_blocks.py`](minimap_blocks.py) | WGA blocks (distant baseline) | minimap2 |
| [`mummer_blocks.py`](mummer_blocks.py) | WGA blocks (1-to-1 truth-style) | nucmer |

Shared helpers: [`align_util.py`](align_util.py). Each script keeps its own **CONFIG** section at the top (CLI overrides CONFIG).

```text
identity = 100 * nident / denominator
distance = 1 - identity / 100
```

| Script | nident | Denominator |
|--------|--------|-------------|
| BLAST windows | BLAST `nident` | window length `W` |
| BLAST `--full-out` | BLAST `nident` | aligned query span |
| minimap | cigar `=` count | query block span |
| mummer | `round(%IDY/100 * LEN2)` | LEN2 |

MUMmer `nident` is back-calculated from `%IDY` (not delta-counted). Minimap counts `=` exactly.

## Block TSV (minimap + mummer)

`q_contig q_start q_end q_strand s_contig s_start s_end s_strand nident aln_span identity distance tool raw_score`

## Usage

```bash
# edit CONFIG near top of each script, then:
./mummer_blocks.py  --sensitive -q query.fa -s subject.fa -o mummer.tsv  -t 8
./minimap_blocks.py --sensitive -q query.fa -s subject.fa -o minimap.tsv -t 8
./sliding_blastn.py -q query.fa -s subject.fa -W 1000 -o windows.tsv --full-out full.tsv -t 8
```

`--sensitive` applies each script's `SENSITIVE_PRESET`. MUMmer uses `DELTA_FILTER_ARGS=["-1"]` (1-to-1); MUMmer 3.x ignores `-t`.
