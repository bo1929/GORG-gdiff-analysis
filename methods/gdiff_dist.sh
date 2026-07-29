#!/usr/bin/env bash
# gdiff dist ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# ANI = DIST_COL of the dist summary (13 = CAPPED_MEDIAN, median of windows
# with >=1 matching k-mer; 4 = raw mean). Samples saved per pair.
# usage: gdiff_dist.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF=../gdiff/gdiff THREADS=8 FORCE=0 ONLY=default DIST_COL=13
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./gdiff_dist_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT/cache" "$OUT/distances" "$OUT/samples"; OUT="$(cd "$OUT" && pwd)"
GDIFF="${GDIFF:-../gdiff/gdiff}"; THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"
DIST_COL="${DIST_COL:-13}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

CONFIGS=(
  "default|k=27,w=35 -l=500|-k 27 -w 35|-l 500"
  "k29w47|k=29,w=47 -l=500|-k 29 -w 47|-l 500"
  "len1000|k=27,w=35 -l=1000|-k 27 -w 35|-l 1000"
  "sensitive|k=27,w=35 -l=500,b=2,n=400|-k 27 -w 35|-l 500 -b 2 --sample-size 400"
)
ONLY="${ONLY:-default}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$OUT/cache/pairs.tsv"
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct'
echo "$HDR" > "$OUT/distances/all_gdiff_dist.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args dist_args <<< "$c"
  want "$name" || continue
  tsv="$OUT/distances/$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; tail -n+2 "$tsv" >> "$OUT/distances/all_gdiff_dist.tsv"; continue
  fi
  echo "$name [$setup]"
  mkdir -p "$OUT/cache/$name" "$OUT/samples/$name"
  cut -f2 "$OUT/cache/pairs.tsv" | sort -u | while read -r s; do
    sk="$OUT/cache/$name/$s.gdiff"
    [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" sketch -o "$sk" --num-threads "$THREADS" $sk_args -i "$(fa "$s")" >/dev/null
  done
  { echo "$HDR"
    while read -r q s _; do
      [ "$q" = "$s" ] && continue
      # shellcheck disable=SC2086
      "$GDIFF" dist -q "$(fa "$q")" -i "$OUT/cache/$name/$s.gdiff" \
        --num-threads "$THREADS" $dist_args \
        -o "$OUT/cache/$name/${q}__${s}.summary" \
        --samples-output "$OUT/samples/$name/${q}__${s}.tsv" 2>/dev/null
      awk -v q="$q" -v s="$s" -v setup="$setup" -v col="$DIST_COL" \
        'NF>=col && $col~/^[0-9.]/ {
           printf "gdiff_dist\t%s\t%s\t%s\t%s\t%.6f\n",setup,q,s,$col,(1-$col)*100; exit
         }' "$OUT/cache/$name/${q}__${s}.summary"
    done < "$OUT/cache/pairs.tsv"
  } > "$tsv"
  tail -n+2 "$tsv" >> "$OUT/distances/all_gdiff_dist.tsv"
done
echo "done -> $OUT/distances (samples: $OUT/samples)"
