# Pairwise regional distance scripts

| Script | Role |
|--------|------|
| [`blastn_blocks.py`](blastn_blocks.py) | BLASTn blocks (optional sliding `-W`) |
| [`minimap_blocks.py`](minimap_blocks.py) | minimap2 blocks |
| [`mummer_blocks.py`](mummer_blocks.py) | nucmer blocks |
| [`align_util.py`](align_util.py) | Shared helpers |
| [`../outlier-detection/classify_intervals.py`](../outlier-detection/classify_intervals.py) | Post-hoc outlier classification |

Edit **CONFIG** at the top of each Python tool. Coords are 1-based inclusive.

```bash
./mummer_blocks.py  --sensitive -q query.fa -s subject.fa -o mummer.tsv -t 8
./minimap_blocks.py --sensitive -q query.fa -s subject.fa -o minimap.tsv -t 8
./blastn_blocks.py              -q query.fa -s subject.fa -o blastn.tsv -t 8
./blastn_blocks.py -W 1000      -q query.fa -s subject.fa -o windows.tsv -t 8
```
