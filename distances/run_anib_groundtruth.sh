#!/usr/bin/env bash
# Deprecated wrapper: ANIb pair-list driver lives in methods/anib.sh.
# For a full-directory all-vs-all ground-truth sweep, pass a pairs file that
# lists every genome pair you care about, or build pairs from the genome dir.
#
# usage: run_anib_groundtruth.sh <genome_dir> [outdir]
#        -> forwards to methods/anib.sh with an all-pairs list of the dir.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
GENOME_DIR="$(cd "${1:?Usage: $0 <genome_dir> [outdir] [suffix]}" && pwd)"
OUT="${2:-./anib_output}"
SUF="${3:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"

# Build all unordered pairs among genomes in GENOME_DIR matching SUF.
pairs="$OUT/all_pairs.tsv"
{
  find "$GENOME_DIR" -maxdepth 1 -type f -name "*$SUF" -print |
    sed "s|.*/||; s|${SUF}\$||" | sort
} > "$OUT/genome_ids.txt"
awk 'NR==FNR { a[NR]=$1; n=NR; next }
     { for (i=1;i<=n;i++) if (a[i] < $1) print a[i] "\t" $1 }
    ' "$OUT/genome_ids.txt" "$OUT/genome_ids.txt" > "$pairs"

echo "Forwarding $(wc -l < "$pairs" | tr -d ' ') pairs to methods/anib.sh" >&2
ONLY="${ONLY:-default}" FORCE="${FORCE:-0}" \
  bash "$HERE/../methods/anib.sh" "$GENOME_DIR" "$pairs" "$OUT" "$SUF"
