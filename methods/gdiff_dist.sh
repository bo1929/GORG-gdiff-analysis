#!/usr/bin/env bash
# gdiff dist ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# ANI = DIST_COL of the dist summary (13 = CAPPED_MEDIAN, median of windows
# with >=1 matching k-mer; 4 = raw mean). SAMPLES=1 also writes per-pair
# sampled regions via --samples-output.
# usage: gdiff_dist.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF=../gdiff/gdiff THREADS=8 FORCE=0 ONLY=default DIST_COL=13 SAMPLES=0
# writes: <outdir>/gdiff-dist/{gdiff-<cfg>.tsv, all_gdiff.tsv, cache/, samples/}
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
MDIR="$OUT/gdiff-dist"
mkdir -p "$MDIR/cache"; OUT="$(cd "$OUT" && pwd)"; MDIR="$OUT/gdiff-dist"
GDIFF="${GDIFF:-../gdiff/gdiff}"; THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"
DIST_COL="${DIST_COL:-13}"; SAMPLES="${SAMPLES:-0}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

CONFIGS=(
  "default|k=27,w=35,-l=500,n=200|-k 27 -w 35|-l 500 --sample-size 200"
  "short-k|k=23,w=31,-l=500,n=200|-k 23 -w 47|-l 500 --sample-size 200"
  "long-window|k=27,w=35,-l=1000,n=200|-k 27 -w 35|-l 1000 --sample-size 200"
  "fast|k=27,w=43,-l=500,b=2,n=100|-k 27 -w 43|-l 500 -b 2 --sample-size 100"
)
ONLY="${ONLY:-default}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$MDIR/cache/pairs.tsv"
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct'

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args dist_args <<< "$c"
  want "$name" || continue
  tsv="$MDIR/gdiff-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; continue
  fi
  echo "$name [$setup]"
  mkdir -p "$MDIR/cache/$name"
  [ "$SAMPLES" = 1 ] && mkdir -p "$MDIR/samples/$name"
  cut -f2 "$MDIR/cache/pairs.tsv" | sort -u | while read -r s; do
    sk="$MDIR/cache/$name/$s.gdiff"
    [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" sketch -o "$sk" --num-threads "$THREADS" $sk_args -i "$(fa "$s")" >/dev/null
  done
  { echo "$HDR"
    while read -r q s _; do
      [ "$q" = "$s" ] && continue
      sum="$MDIR/cache/$name/${q}__${s}.summary"
      # shellcheck disable=SC2086
      set -- "$GDIFF" dist -q "$(fa "$q")" -i "$MDIR/cache/$name/$s.gdiff" \
        --num-threads "$THREADS" $dist_args -o "$sum"
      [ "$SAMPLES" = 1 ] && set -- "$@" --samples-output "$MDIR/samples/$name/${q}__${s}.tsv"
      "$@" 2>/dev/null
      awk -v q="$q" -v s="$s" -v setup="$setup" -v col="$DIST_COL" \
        'NF>=col && $col~/^[0-9.]/ {
           printf "gdiff_dist\t%s\t%s\t%s\t%s\t%.6f\n",setup,q,s,$col,(1-$col)*100; exit
         }' "$sum"
    done < "$MDIR/cache/pairs.tsv"
  } > "$tsv"
done

{ echo "$HDR"
  for c in "${CONFIGS[@]}"; do
    IFS='|' read -r name _ <<< "$c"
    if want "$name" && [ -s "$MDIR/gdiff-$name.tsv" ]; then tail -n+2 "$MDIR/gdiff-$name.tsv"; fi
  done
} > "$MDIR/all_gdiff.tsv"
echo "done -> $MDIR"
