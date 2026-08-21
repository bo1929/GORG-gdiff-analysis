#!/usr/bin/env bash
# gdiff map over pairs.tsv. Cache = subject sketches only.
# usage: gdiff-map.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF JOBS=8 FORCE=0 ONLY=default
# out: blocks/gdiff-map/<cfg>/<q>__<s>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/gdiff-map"; OUTDIR="$OUT/blocks/gdiff-map"
mkdir -p "$CACHE" "$OUTDIR"
GDIFF="${GDIFF:-../gidiff/gdiff}"
JOBS="${JOBS:-${THREADS:-8}}"; FORCE="${FORCE:-0}"; ONLY="${ONLY:-default}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

DISTS="0.0001 0.01 0.025 0.05 0.075 0.1 0.15 0.20"
CONFIGS=(
  "default|k=27,w=35,l=500|-k 27 -w 35|-l 500 -d $DISTS"
  "short-win|k=27,w=35,l=50|-k 27 -w 35|-l 50 -d $DISTS"
  "medium-win|k=27,w=35,l=200|-k 27 -w 35|-l 200 -d $DISTS"
  "fast|k=27,w=45,b=4,l=500|-k 27 -w 45|-b 4 -l 500 -d $DISTS"
)
load_pairs

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args map_args <<< "$c"
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
    ( "$GDIFF" --num-threads 1 map "$(fa "$q")" "$skdir/$s.gdiff" $map_args -o "$out" >/dev/null 2>&1 ) &
    n=$((n + 1)); throttle "$n"
  done < "$CACHE/pairs.tsv"
  wait
done
echo "done -> $OUTDIR" >&2
