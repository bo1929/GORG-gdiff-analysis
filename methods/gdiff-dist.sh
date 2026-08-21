#!/usr/bin/env bash
# Symmetric gdiff dist ANI over pairs.tsv. Cache = sketches only.
# usage: gdiff-dist.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: GDIFF JOBS=8 FORCE=0 ONLY=all DIST_COL=4 SAMPLES=0
# out: distances/gdiff-<cfg>.tsv  [samples/gdiff-<cfg>.tsv]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/gdiff"; DIST_DIR="$OUT/distances"; SAMP_DIR="$OUT/samples"
mkdir -p "$CACHE" "$DIST_DIR"
GDIFF="${GDIFF:-../gdiff/gdiff}"
JOBS="${JOBS:-${THREADS:-8}}"; FORCE="${FORCE:-0}"
DIST_COL="${DIST_COL:-4}"; SAMPLES="${SAMPLES:-0}"; ONLY="${ONLY:-all}"
[ -x "$GDIFF" ] || { echo "set GDIFF=/path/to/gdiff" >&2; exit 1; }

CONFIGS=(
  "sensible-cfg|k=25,w=37,h=11,frac=0.1,-l=1000,n=200,b=6|-k 27 -h 11 -w 37 --frac 0.1|-l 1000 --sample-size 200 -b 6"
  "short-k|k=23,w=31,h=11,frac=0.5,-l=1000,n=200,b=4|-k 23 -h 11 -w 47 --frac 0.2|-l 1000 --sample-size 200 -b 4"
  "long-window|k=27,w=37,frac=0.5,-l=2000,n=200,b=4|-k 27 -w 37 --frac 0.2|-l 2000 --sample-size 200 -b 2"
  "gigantic-window|k=27,w=37,frac=0.5,-l=5000,n=200,b=4|-k 27 -w 37 --frac 0.2|-l 5000 --sample-size 200 -b 4"
  "full-scale|k=27,w=37,frac=0.5,-l=10000,n=500,b=2|-k 27 -w 37 --frac 0.2|-l 10000 --sample-size 500 -b 2"
  "fast|k=27,w=43,frac=0.2,-l=500,n=100,b=4|-k 27 -w 43 --frac 0.2|-l 500 --sample-size 100 -b 4"
)
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct'
SAMP_HDR=$'config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg'
load_pairs

# Print one distance line. Optional: append sample rows to $6.
run_pair() {
  local cfg="$1" setup="$2" dist_args="$3" q="$4" s="$5" skdir="$6" samp_out="${7:-}"
  local fwd rev d_fwd d_rev d samp
  fwd="$(mktemp)"; rev="$(mktemp)"
  # shellcheck disable=SC2086
  "$GDIFF" --num-threads 1 dist "$(fa "$q")" "$skdir/$s.gdiff" $dist_args -o "$fwd" >/dev/null 2>&1 || true
  # shellcheck disable=SC2086
  "$GDIFF" --num-threads 1 dist "$(fa "$s")" "$skdir/$q.gdiff" $dist_args -o "$rev" >/dev/null 2>&1 || true
  d_fwd=$(awk -v c="$DIST_COL" 'NF>=c && $c~/^[0-9.]/{print $c; exit}' "$fwd" || true)
  d_rev=$(awk -v c="$DIST_COL" 'NF>=c && $c~/^[0-9.]/{print $c; exit}' "$rev" || true)
  rm -f "$fwd" "$rev"
  [ -n "${d_fwd:-}" ] || [ -n "${d_rev:-}" ] || return 0
  if [ -z "${d_fwd:-}" ]; then d=$d_rev
  elif [ -z "${d_rev:-}" ]; then d=$d_fwd
  else d=$(awk -v a="$d_fwd" -v b="$d_rev" 'BEGIN{printf "%.9f",(a+b)/2}'); fi
  printf 'gdiff_dist\t%s\t%s\t%s\t%s\t%.6f\n' "$setup" "$q" "$s" "$d" \
    "$(awk -v d="$d" 'BEGIN{print (1-d)*100}')"
  if [ -n "$samp_out" ]; then
    samp="$(mktemp)"
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 dist "$(fa "$q")" "$skdir/$s.gdiff" \
      $dist_args --output-samples -o "$samp" >/dev/null 2>&1 || true
    awk -v cfg="$cfg" -v q="$q" -v s="$s" 'NF{print cfg"\t"q"\t"s"\t"$0}' "$samp" >> "$samp_out"
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 dist "$(fa "$s")" "$skdir/$q.gdiff" \
      $dist_args --output-samples -o "$samp" >/dev/null 2>&1 || true
    awk -v cfg="$cfg" -v q="$s" -v s="$q" 'NF{print cfg"\t"q"\t"s"\t"$0}' "$samp" >> "$samp_out"
    rm -f "$samp"
  fi
}

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args dist_args <<< "$c"
  want "$name" || continue
  tsv="$DIST_DIR/gdiff-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then echo "$name: skip"; continue; fi
  echo "$name [$setup] jobs=$JOBS" >&2
  skdir="$CACHE/$name"; mkdir -p "$skdir"
  parts="$(mktemp -d)"

  # --- sketch ---
  todo=0
  while read -r g; do
    [ "$FORCE" != 1 ] && [ -s "$skdir/$g.gdiff" ] && continue
    todo=$((todo + 1))
  done < "$CACHE/genomes.txt"
  if [ "$todo" -eq 0 ]; then
    echo "  sketch: cached" >&2
  else
    progress 0 "$todo" "sketch"
    n=0
    while read -r g; do
      [ "$FORCE" != 1 ] && [ -s "$skdir/$g.gdiff" ] && continue
      # shellcheck disable=SC2086
      ( "$GDIFF" --num-threads 1 sketch -o "$skdir/$g.gdiff" $sk_args -i "$(fa "$g")" >/dev/null 2>&1 ) &
      n=$((n + 1))
      if (( n % JOBS == 0 )); then wait; progress "$n" "$todo" "sketch"; fi
    done < "$CACHE/genomes.txt"
    wait; progress "$todo" "$todo" "sketch"
  fi

  # --- dist (shards: at most JOBS temp files at a time) ---
  np=$(awk '$1!=$2{n++} END{print n+0}' "$CACHE/pairs.tsv")
  echo "$HDR" > "$tsv"
  : > "$parts/out"
  if [ "$SAMPLES" = 1 ]; then
    mkdir -p "$SAMP_DIR"
    echo "$SAMP_HDR" > "$SAMP_DIR/gdiff-$name.tsv"
    : > "$parts/samp"
  fi
  progress 0 "$np" "dist"
  n=0; bi=0
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    sp=""
    if [ "$SAMPLES" = 1 ]; then sp="$parts/s.$bi"; : > "$sp"; fi
    (
      run_pair "$name" "$setup" "$dist_args" "$q" "$s" "$skdir" "$sp" >> "$parts/p.$bi"
    ) &
    n=$((n + 1)); bi=$((bi + 1))
    if (( n % JOBS == 0 )); then
      wait
      cat "$parts"/p.* >> "$parts/out" 2>/dev/null || true
      if [ "$SAMPLES" = 1 ]; then cat "$parts"/s.* >> "$parts/samp" 2>/dev/null || true; fi
      rm -f "$parts"/p.* "$parts"/s.*
      bi=0
      progress "$n" "$np" "dist"
    fi
  done < "$CACHE/pairs.tsv"
  wait
  cat "$parts"/p.* >> "$parts/out" 2>/dev/null || true
  if [ "$SAMPLES" = 1 ]; then cat "$parts"/s.* >> "$parts/samp" 2>/dev/null || true; fi
  progress "$np" "$np" "dist"
  cat "$parts/out" >> "$tsv"
  if [ "$SAMPLES" = 1 ] && [ -s "$parts/samp" ]; then
    cat "$parts/samp" >> "$SAMP_DIR/gdiff-$name.tsv"
  fi
  rm -rf "$parts"
  echo "  $np pairs -> $tsv" >&2
done

emit_all_distances gdiff "$HDR" "${CONFIGS[@]}"
echo "done -> $DIST_DIR" >&2
