#!/usr/bin/env bash
# gdiff detect (two-sided outlier regions) for the pairs in <pairs.tsv>,
# one run per CONFIGS entry. detect args need -l (k-mer window / min length).
# usage: gdiff_detect.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF=../gdiff/gdiff THREADS=8 FORCE=0 ONLY=default (name(s), or "all")
# writes: <outdir>/detect/<cfg>/<q>__<s>.tsv, cache in <outdir>/cache/gdiff-detect/
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/gdiff-detect"; OUTDIR="$OUT/gdiff-detect"
mkdir -p "$CACHE" "$OUTDIR"
GDIFF="${GDIFF:-../gdiff/gdiff}"; THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

CONFIGS=(
  "k27w35l500|k=27,w=35,-l=500|-k 27 -w 35|-l 500"
  "k29w47l500|k=29,w=47,-l=500|-k 29 -w 47|-l 500"
  "k27w35l1000|k=27,w=35,-l=1000|-k 27 -w 35|-l 1000"
  "k27w35l500-pq|k=27,w=35,-l=500,per-query|-k 27 -w 35|-l 500 --fit-scope per-query"
)
ONLY="${ONLY:-all}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$CACHE/pairs.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args det_args <<< "$c"
  want "$name" || continue
  echo "$name [$setup]"
  mkdir -p "$CACHE/$name" "$OUTDIR/$name"
  cut -f2 "$CACHE/pairs.tsv" | sort -u | while read -r s; do
    sk="$CACHE/$name/$s.gdiff"
    [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" sketch -o "$sk" --num-threads "$THREADS" $sk_args -i "$(fa "$s")" >/dev/null
  done
  echo "!!"
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    out="$OUTDIR/$name/${q}__${s}.tsv"
    [ "$FORCE" != 1 ] && [ -s "$out" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" detect -q "$(fa "$q")" -i "$CACHE/$name/$s.gdiff" \
      --num-threads "$THREADS" $det_args -o "$out" 2>/dev/null
  done < "$CACHE/pairs.tsv"
done
echo "done -> $OUTDIR"
