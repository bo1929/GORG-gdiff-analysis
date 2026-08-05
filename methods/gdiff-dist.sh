#!/usr/bin/env bash
# gdiff dist ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# ANI = DIST_COL of the 12-col dist summary (4 = mean, 9 = Q50 median).
# SAMPLES=1 also writes per-pair sampled regions via --samples-output and
# concatenates them (prefixed with config + pair info) into
# <outdir>/samples/all_<cfg>.tsv and <outdir>/samples/all_samples.tsv.
# Concatenated columns: config, genome_a, genome_b, then the raw dist sample
# row (query_id, seq_len, start, end, strand, ref_id, dist, info, lr_bg,
# lr_ub); dist is NaN for unmapped windows (no k-mer hits).
# usage: gdiff_dist.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF=../gdiff/gdiff THREADS=8 FORCE=0 ONLY=default DIST_COL=9 SAMPLES=0
# writes: <outdir>/distances/{gdiff-<cfg>.tsv, all_gdiff.tsv},
#         <outdir>/samples/<cfg>/ + all_<cfg>.tsv + all_samples.tsv (SAMPLES=1),
#         cache in <outdir>/cache/gdiff-dist/
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/gdiff-dist"; OUTDIR="$OUT/distances"
mkdir -p "$CACHE" "$OUTDIR"
GDIFF="${GDIFF:-../gdiff/gdiff}"; THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"
DIST_COL="${DIST_COL:-9}"; SAMPLES="${SAMPLES:-0}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

CONFIGS=(
  "default|k=27,w=35,-l=500,n=200|-k 27 -w 35|-l 500 --sample-size 200"
  "short-k|k=23,w=31,-l=500,n=200|-k 23 -w 47|-l 500 --sample-size 200"
  "long-window|k=27,w=35,-l=1000,n=200|-k 27 -w 35|-l 1000 --sample-size 200"
  "gigantic-window|k=27,w=35,-l=5000,n=300|-k 27 -w 35|-l 5000 --sample-size 300"
  "full-scale|k=27,w=37,-l=10000,n=500|-k 27 -w 37|-l 10000 --sample-size 500"
  "fast|k=27,w=43,-l=500,b=2,n=100|-k 27 -w 43|-l 500 -b 2 --sample-size 100"
)
ONLY="${ONLY:-all}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$CACHE/pairs.tsv"
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct'
SAMPLES_HDR=$'config\tgenome_a\tgenome_b\tquery_id\tseq_len\tstart\tend\tstrand\tref_id\tdist\tinfo\tlr_bg\tlr_ub'

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args dist_args <<< "$c"
  want "$name" || continue
  tsv="$OUTDIR/gdiff-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; continue
  fi
  echo "$name [$setup]"
  mkdir -p "$CACHE/$name"
  concat="$OUT/samples/all_$name.tsv"
  if [ "$SAMPLES" = 1 ]; then
    mkdir -p "$OUT/samples/$name"
    echo "$SAMPLES_HDR" > "$concat"
  fi
  cut -f2 "$CACHE/pairs.tsv" | sort -u | while read -r s; do
    sk="$CACHE/$name/$s.gdiff"
    [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" sketch -o "$sk" --num-threads "$THREADS" $sk_args -i "$(fa "$s")" >/dev/null
  done
  { echo "$HDR"
    while read -r q s _; do
      [ "$q" = "$s" ] && continue
      sum="$CACHE/$name/${q}__${s}.summary"
      samples_tsv="$OUT/samples/$name/${q}__${s}.tsv"
      # shellcheck disable=SC2086
      set -- "$GDIFF" dist -q "$(fa "$q")" -i "$CACHE/$name/$s.gdiff" \
        --num-threads "$THREADS" $dist_args -o "$sum"
      [ "$SAMPLES" = 1 ] && set -- "$@" --samples-output "$samples_tsv"
      "$@" 2>/dev/null
      awk -v q="$q" -v s="$s" -v setup="$setup" -v col="$DIST_COL" \
        'NF>=col && $col~/^[0-9.]/ {
           printf "gdiff_dist\t%s\t%s\t%s\t%s\t%.6f\n",setup,q,s,$col,(1-$col)*100; exit
         }' "$sum"
      if [ "$SAMPLES" = 1 ] && [ -s "$samples_tsv" ]; then
        awk -v name="$name" -v q="$q" -v s="$s" \
          'NF>0 { print name "\t" q "\t" s "\t" $0 }' "$samples_tsv" >> "$concat"
      fi
    done < "$CACHE/pairs.tsv"
  } > "$tsv"
done

{ echo "$HDR"
  for c in "${CONFIGS[@]}"; do
    IFS='|' read -r name _ <<< "$c"
    if want "$name" && [ -s "$OUTDIR/gdiff-$name.tsv" ]; then tail -n+2 "$OUTDIR/gdiff-$name.tsv"; fi
  done
} > "$OUTDIR/all_gdiff.tsv"

if [ "$SAMPLES" = 1 ]; then
  { echo "$SAMPLES_HDR"
    for c in "${CONFIGS[@]}"; do
      IFS='|' read -r name _ <<< "$c"
      if want "$name" && [ -s "$OUT/samples/all_$name.tsv" ]; then tail -n+2 "$OUT/samples/all_$name.tsv"; fi
    done
  } > "$OUT/samples/all_samples.tsv"
  echo "done -> $OUTDIR and $OUT/samples"
else
  echo "done -> $OUTDIR"
fi
