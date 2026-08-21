#!/usr/bin/env bash
# gdiff detect over pairs.tsv. Cache = subject sketches only.
# usage: gdiff-detect.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF JOBS=8 FORCE=0 ONLY=all
# out: blocks/gdiff-detect/<cfg>/<q>__<s>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/gdiff-detect"; OUTDIR="$OUT/blocks/gdiff-detect"
mkdir -p "$CACHE" "$OUTDIR"
GDIFF="${GDIFF:-../gdiff/gdiff}"
JOBS="${JOBS:-${THREADS:-8}}"; FORCE="${FORCE:-0}"; ONLY="${ONLY:-all}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

CONFIGS=(
  "k27w35l500|k=27,w=35,l=500|-k 27 -w 35|-l 500"
  "k25w41l500|k=25,w=41,l=500|-k 25 -w 41|-l 500"
  "k27w35l1000|k=27,w=35,l=1000|-k 27 -w 35|-l 1000"
  "k27w35l500-pq|k=27,w=35,l=500,per-sequence|-k 27 -w 35|-l 500 --per-sequence"
)
load_pairs

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args det_args <<< "$c"
  want "$name" || continue
  echo "$name [$setup] jobs=$JOBS" >&2
  skdir="$CACHE/$name"; mkdir -p "$skdir" "$OUTDIR/$name"

  cut -f2 "$CACHE/pairs.tsv" | sort -u > "$skdir/subjects.txt"
  while read -r s; do
    [ "$FORCE" != 1 ] && [ -s "$skdir/$s.gdiff" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 sketch -o "$skdir/$s.gdiff" $sk_args -i "$(fa "$s")" >/dev/null 2>&1
  done < "$skdir/subjects.txt"

  n=0
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    out="$OUTDIR/$name/${q}__${s}.tsv"
    [ "$FORCE" != 1 ] && [ -s "$out" ] && continue
    # shellcheck disable=SC2086
    ( "$GDIFF" --num-threads 1 detect "$(fa "$q")" "$skdir/$s.gdiff" $det_args -o "$out" >/dev/null 2>&1 \
        || echo "  detect failed: $q vs $s" >&2 ) &
    n=$((n + 1)); throttle "$n"
  done < "$CACHE/pairs.tsv"
  wait
done
echo "done -> $OUTDIR" >&2
