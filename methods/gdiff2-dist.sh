#!/usr/bin/env bash
# usage: gdiff2-dist.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF2 JOBS=8 FORCE=0 ONLY=all SAMPLES=1
# out: samples/gdiff2-<cfg>.tsv  (distance summaries are not written)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./output}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache gdiff2
SAMP_DIR="$OUT/samples"; mkdir -p "$SAMP_DIR"
GDIFF2="${GDIFF2:-../gdiff/gdiff2}"
JOBS="${JOBS:-${THREADS:-32}}"; FORCE="${FORCE:-1}"
SAMPLES="${SAMPLES:-1}"; ONLY="${ONLY:-all}"
[ -x "$GDIFF2" ] || { echo "set GDIFF2=/path/to/gdiff2" >&2; exit 1; }

CONFIGS=(
  "xyz|k=23,w=23,h=11,frac=0.33,-l=500,n=1000|-k 23 -h 11 -w 23 --frac 0.33|-l 500 --sample-size 1000"
  "abc|k=23,w=23,frac=0.5,-l=500,n=1000|-k 23 -w 23 --frac 0.5|-l 500 --sample-size 1000"
)
SAMPLES_HEADER=$'config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg\tlr_ub'
load_pairs

NG="$(wc -l < "$CACHE/genomes.txt" | tr -d ' ')"
[ "$NG" -ge 2 ] || { echo "fewer than 2 genomes in $PAIRS_FILE" >&2; exit 1; }

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args dist_args <<< "$c"
  want "$name" || continue
  samp_tsv="$SAMP_DIR/gdiff2-$name.tsv"
  [ "$FORCE" != 1 ] && [ -s "$samp_tsv" ] && { echo "$name: skip"; continue; }
  echo "$name [$setup] jobs=$JOBS" >&2
  skdir="$CACHE/$name"; mkdir -p "$skdir"
  parts="$(mktemp -d)"
  bundle="$skdir/all.g2"

  # Rebuild the bundle only when forced or when the genome set changed.
  rebuild=0
  [ "$FORCE" = 1 ] && rebuild=1
  [ -s "$bundle" ] || rebuild=1
  if [ -s "$skdir/genomes.txt" ] && ! cmp -s "$CACHE/genomes.txt" "$skdir/genomes.txt"; then rebuild=1; fi
  if [ ! -s "$skdir/genomes.txt" ]; then rebuild=1; fi

  if [ "$rebuild" -eq 0 ]; then
    echo "  sketch: cached" >&2
  else
    rm -f "$bundle"
    while read -r g; do printf '%s\t%s\n' "$g" "$(fa "$g")"; done < "$CACHE/genomes.txt" > "$parts/input.list"
    echo "  sketch2 [$NG genomes, threads=$JOBS] -> $bundle" >&2
    # shellcheck disable=SC2086
    "$GDIFF2" --num-threads "$JOBS" sketch2 --input-list "$parts/input.list" $sk_args $dist_args -o "$bundle" || {
      echo "sketch2 failed for $name" >&2; rm -rf "$parts"; exit 1; }
    cp "$CACHE/genomes.txt" "$skdir/genomes.txt"
  fi

  # --- single dist2: within-mode over the whole bundle, internal threads ---
  echo "  dist [$((NG * (NG - 1) / 2)) pairs, threads=$JOBS]" >&2

  # Canonical unordered pair keys (what the v1 script emitted from this pairs
  # file); the full all-vs-all matrix is exact, subsets get filtered out.
  awk 'NF>=2{a=$1;b=$2;if(a<b)k=a"|"b;else k=b"|"a;if(!(k in S)){S[k]=1;print k}}' \
    "$CACHE/pairs.tsv" | sort -u > "$parts/canon"
  if [ ! -s "$parts/canon" ]; then rm -rf "$parts"; continue; fi

  if [ "$SAMPLES" = 1 ]; then
    mkdir -p "$SAMP_DIR"
    echo "  samples -> $SAMP_DIR/gdiff2-$name.tsv" >&2
    # shellcheck disable=SC2086
    "$GDIFF2" --num-threads "$JOBS" dist2 "$bundle" --output-samples -o "$parts/samp.tsv" || true
    {
      echo "$SAMPLES_HEADER"
      awk -F'\t' -v cfg="$name" -v f="$parts/canon" '
        BEGIN { OFS="\t"; while ((getline l < f) > 0) fl[l] = 1; close(f) }
        {
          if ($1 == "ab") { a = $6; b = $7 } else { a = $7; b = $6 }
          k = (a < b) ? a "|" b : b "|" a
          if (!(k in fl)) next
          print cfg, a, b, $2, $3, $4, $5, $7, $8, $9, $10
        }' "$parts/samp.tsv"
    } >> "$SAMP_DIR/gdiff2-$name.tsv"
  fi

  rm -rf "$parts"
done

echo "done -> $SAMP_DIR" >&2
