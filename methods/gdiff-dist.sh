#!/usr/bin/env bash
# gdiff dist ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# Distance is SYMMETRIC: dist is run both A->B and B->A; the reported value
# is the average of the two directional D_MED estimates.
# ANI = DIST_COL of the 6-col dist summary
# (query_file, reference, N, D_MED, D_MED_FILT, N_REMOVED):
# 4 = D_MED (median MLE distance), 5 = D_MED_FILT (median after outlier
# removal). The new CLI has no mean column; D_MED replaces the old Q50.
# SAMPLES=1 additionally runs dist --output-samples (a separate run; sampling
# is seeded and therefore identical) and concatenates the per-window rows
# (prefixed with config + pair info) into <outdir>/gdiff-samples/all_<cfg>.tsv
# and <outdir>/gdiff-samples/all_samples.tsv.
# Concatenated columns: config, genome_a, genome_b, then the raw dist sample
# row (qid, start, end, strand, reference, d, lr_bg); d is NaN for unmapped
# windows (no k-mer hits).
#
# Parallelism: pairs (and sketches) run as background jobs with gdiff
# --num-threads 1. Within-call multithreading does not help single-genome
# sketch / single-ref dist; job-level parallelism does.
#
# usage: gdiff_dist.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF=../gdiff/gdiff JOBS=8 FORCE=0 ONLY=default DIST_COL=4 SAMPLES=0
#      THREADS is accepted as an alias for JOBS (legacy).
# writes: <outdir>/distances/{gdiff-<cfg>.tsv, all_gdiff.tsv},
#         <outdir>/gdiff-samples/<cfg>/ + all_<cfg>.tsv + all_samples.tsv (SAMPLES=1),
#         cache in <outdir>/cache/gdiff-dist/
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/gdiff-dist"; OUTDIR="$OUT/distances"
mkdir -p "$CACHE" "$OUTDIR"
GDIFF="${GDIFF:-../gdiff/gdiff}"
# Prefer JOBS; fall back to THREADS for older call sites, then 8.
JOBS="${JOBS:-${THREADS:-8}}"
FORCE="${FORCE:-0}"
DIST_COL="${DIST_COL:-4}"; SAMPLES="${SAMPLES:-1}"
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

# bash 3.2-compatible throttle: wait after every JOBS launches
throttle() {
  local n="$1"
  if (( n > 0 && n % JOBS == 0 )); then wait; fi
}

# Run one pair (both directions + optional samples) and write a one-line row file.
run_pair() {
  local name="$1" setup="$2" dist_args="$3" q="$4" s="$5"
  local sum_fwd sum_rev samples_fwd samples_rev row
  sum_fwd="$CACHE/$name/${q}__${s}.summary"
  sum_rev="$CACHE/$name/${s}__${q}.summary"
  samples_fwd="$OUT/gdiff-samples/$name/${q}__${s}.tsv"
  samples_rev="$OUT/gdiff-samples/$name/${s}__${q}.tsv"
  row="$CACHE/$name/rows/${q}__${s}.row"

  # shellcheck disable=SC2086
  "$GDIFF" --num-threads 1 dist "$(fa "$q")" "$CACHE/$name/$s.gdiff" \
    $dist_args -o "$sum_fwd" 2>/dev/null
  # shellcheck disable=SC2086
  "$GDIFF" --num-threads 1 dist "$(fa "$s")" "$CACHE/$name/$q.gdiff" \
    $dist_args -o "$sum_rev" 2>/dev/null
  if [ "$SAMPLES" = 1 ]; then
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 dist "$(fa "$q")" "$CACHE/$name/$s.gdiff" \
      $dist_args --output-samples -o "$samples_fwd" 2>/dev/null
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 dist "$(fa "$s")" "$CACHE/$name/$q.gdiff" \
      $dist_args --output-samples -o "$samples_rev" 2>/dev/null
  fi

  awk -v q="$q" -v s="$s" -v setup="$setup" -v col="$DIST_COL" \
      -v fwd="$sum_fwd" -v rev="$sum_rev" '
    BEGIN {
      d_fwd = ""; d_rev = ""
      while ((getline line < fwd) > 0) {
        n = split(line, a, "\t")
        if (n >= col && a[col] ~ /^[0-9.]/) { d_fwd = a[col]+0; break }
      }
      close(fwd)
      while ((getline line < rev) > 0) {
        n = split(line, a, "\t")
        if (n >= col && a[col] ~ /^[0-9.]/) { d_rev = a[col]+0; break }
      }
      close(rev)
      if (d_fwd == "" && d_rev == "") exit
      if (d_fwd == "") d = d_rev
      else if (d_rev == "") d = d_fwd
      else d = (d_fwd + d_rev) / 2
      printf "gdiff_dist\t%s\t%s\t%s\t%.9f\t%.6f\n", setup, q, s, d, (1-d)*100
    }' /dev/null > "$row"
}

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$CACHE/pairs.tsv"
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct'
SAMPLES_HDR=$'config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg'

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args dist_args <<< "$c"
  want "$name" || continue
  tsv="$OUTDIR/gdiff-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; continue
  fi
  echo "$name [$setup] jobs=$JOBS"
  mkdir -p "$CACHE/$name/rows"
  concat="$OUT/gdiff-samples/all_$name.tsv"
  if [ "$SAMPLES" = 1 ]; then
    mkdir -p "$OUT/gdiff-samples/$name"
    echo "$SAMPLES_HDR" > "$concat"
  fi

  # sketch all genomes that appear in either column (needed for both directions).
  # Use a file (not a pipe) so background jobs stay in this shell for `wait`.
  { cut -f1 "$CACHE/pairs.tsv"; cut -f2 "$CACHE/pairs.tsv"; } | sort -u > "$CACHE/$name/genomes.txt"
  n=0
  while read -r s; do
    sk="$CACHE/$name/$s.gdiff"
    [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 sketch -o "$sk" $sk_args -i "$(fa "$s")" >/dev/null &
    n=$((n + 1))
    throttle "$n"
  done < "$CACHE/$name/genomes.txt"
  wait

  # pairs as parallel jobs; each writes CACHE/$name/rows/<q>__<s>.row
  n=0
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    run_pair "$name" "$setup" "$dist_args" "$q" "$s" &
    n=$((n + 1))
    throttle "$n"
  done < "$CACHE/pairs.tsv"
  wait

  { echo "$HDR"
    while read -r q s _; do
      [ "$q" = "$s" ] && continue
      row="$CACHE/$name/rows/${q}__${s}.row"
      [ -s "$row" ] && cat "$row"
    done < "$CACHE/pairs.tsv"
  } > "$tsv"

  if [ "$SAMPLES" = 1 ]; then
    while read -r q s _; do
      [ "$q" = "$s" ] && continue
      samples_fwd="$OUT/gdiff-samples/$name/${q}__${s}.tsv"
      samples_rev="$OUT/gdiff-samples/$name/${s}__${q}.tsv"
      if [ -s "$samples_fwd" ]; then
        awk -v name="$name" -v q="$q" -v s="$s" \
          'NF>0 { print name "\t" q "\t" s "\t" $0 }' "$samples_fwd" >> "$concat"
      fi
      if [ -s "$samples_rev" ]; then
        awk -v name="$name" -v q="$s" -v s="$q" \
          'NF>0 { print name "\t" q "\t" s "\t" $0 }' "$samples_rev" >> "$concat"
      fi
    done < "$CACHE/pairs.tsv"
  fi
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
      if want "$name" && [ -s "$OUT/gdiff-samples/all_$name.tsv" ]; then tail -n+2 "$OUT/gdiff-samples/all_$name.tsv"; fi
    done
  } > "$OUT/gdiff-samples/all_samples.tsv"
  echo "done -> $OUTDIR and $OUT/gdiff-samples"
else
  echo "done -> $OUTDIR"
fi
