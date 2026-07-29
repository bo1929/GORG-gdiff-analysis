#!/usr/bin/env bash
# gdiff map (distance-based intervals) for the pairs in <pairs.tsv>,
# one run per CONFIGS entry. map args need -l and -d (1 or 8 thresholds).
# usage: gdiff_map.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF=../gdiff/gdiff THREADS=8 FORCE=0 ONLY=default (name(s), or "all")
# writes: <outdir>/mapping/<cfg>/<q>__<s>.tsv, cache in <outdir>/cache/gdiff-map/
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/gdiff-map"; MAPS="$OUT/mapping"
mkdir -p "$CACHE" "$MAPS"
GDIFF="${GDIFF:-../gdiff/gdiff}"; THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

DISTS="0.01 0.025 0.05 0.075 0.1 0.125 0.15 0.175"
CONFIGS=(
  "default|k=27,w=35 -l=500|-k 27 -w 35|-l 500 -d $DISTS"
  "k29w47|k=29,w=47 -l=500|-k 29 -w 47|-l 500 -d $DISTS"
  "len1000|k=27,w=35 -l=1000|-k 27 -w 35|-l 1000 -d $DISTS"
  "bin4|k=27,w=35 -b=4,l=500|-k 27 -w 35|-b 4 -l 500 -d $DISTS"
)
ONLY="${ONLY:-default}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$CACHE/pairs.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args map_args <<< "$c"
  want "$name" || continue
  echo "$name [$setup]"
  mkdir -p "$CACHE/$name" "$MAPS/$name"
  cut -f2 "$CACHE/pairs.tsv" | sort -u | while read -r s; do
    sk="$CACHE/$name/$s.gdiff"
    [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" sketch -o "$sk" --num-threads "$THREADS" $sk_args -i "$(fa "$s")" >/dev/null
  done
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    out="$MAPS/$name/${q}__${s}.tsv"
    [ "$FORCE" != 1 ] && [ -s "$out" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" map -q "$(fa "$q")" -i "$CACHE/$name/$s.gdiff" \
      --num-threads "$THREADS" $map_args -o "$out" 2>/dev/null
  done < "$CACHE/pairs.tsv"
done
echo "done -> $MAPS"
