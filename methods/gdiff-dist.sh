#!/usr/bin/env bash
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
  "sensible-cfg|k=27,w=35,-l=1000,n=200,|-k 27 -w 35|-l 1000 --sample-size 200"
  "short-k|k=23,w=31,-l=1000,n=200|-k 23 -w 47|-l 1000 --sample-size 200"
  "long-window|k=29,w=37,-l=2000,n=200,frac=0.5|-k 29 -w 37|-l 2000 --frac 0.5 --sample-size 200"
  "gigantic-window|k=29,w=37,-l=5000,n=200,frac=0.5|-k 29 -w 37|-l 5000 --frac 0.5 --sample-size 200"
  "full-scale|k=29,w=37,-l=10000,n=500|-k 29 -w 37|-l 10000 --sample-size 500"
  "fast|k=27,w=43,-l=500,b=4,n=100,frac=0.2|-k 27 -w 43|-l 500 -b 4 --frac 0.2 --sample-size 100"
)
ONLY="${ONLY:-all}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

# bash 3.2-compatible throttle: wait after every JOBS launches, refresh progress.
# $1=launched count, $2=total, $3=label, $4=progress-marker dir
throttle() {
  local n="$1" total="${2:-0}" label="${3:-}" pdir="${4:-}" done
  # leave the final batch for finish_progress so we do not double-print 100%
  if (( n > 0 && n % JOBS == 0 && n < total )); then
    wait
    if [ -n "$pdir" ] && [ "$total" -gt 0 ]; then
      done=$(find "$pdir" -type f 2>/dev/null | wc -l | tr -d ' ')
      [ -z "$done" ] && done=0
      draw_progress "$done" "$total" "$label"
    fi
  fi
}

# ASCII progress bar on stderr: label [=====>    ] 12/40 (30%)
draw_progress() {
  local done="$1" total="$2" label="$3"
  local width=30 filled empty i bar pct
  [ "$total" -gt 0 ] || return 0
  pct=$(( done * 100 / total ))
  filled=$(( done * width / total ))
  [ "$filled" -gt "$width" ] && filled=$width
  empty=$(( width - filled ))
  bar=""
  i=0; while [ "$i" -lt "$filled" ]; do bar="${bar}="; i=$((i + 1)); done
  if [ "$filled" -lt "$width" ] && [ "$done" -lt "$total" ]; then
    bar="${bar}>"
    empty=$(( empty - 1 ))
  fi
  i=0; while [ "$i" -lt "$empty" ]; do bar="${bar} "; i=$((i + 1)); done
  printf '\r%-10s [%s] %d/%d (%d%%)' "$label" "$bar" "$done" "$total" "$pct" >&2
}

finish_progress() {
  local pdir="$1" total="$2" label="$3" done
  wait
  done=$(find "$pdir" -type f 2>/dev/null | wc -l | tr -d ' ')
  [ -z "$done" ] && done=0
  draw_progress "$done" "$total" "$label"
  echo >&2
}

# Prefer GNU /usr/bin/time -v; fall back to POSIX -p (BSD/macOS).
if /usr/bin/time -v true >/dev/null 2>&1; then
  TIME_STYLE=v
else
  TIME_STYLE=p
fi

# Parse /usr/bin/time logs under $1; print avg wall (and user/sys).
# $2 = unit label (e.g. genome, pair).
# Stream files one-by-one (no glob / batched argv) to avoid ARG_MAX.
summarize_times() {
  local tdir="$1" unit="${2:-item}" awk_prog
  [ -d "$tdir" ] || return 0
  find "$tdir" -type f -name '*.time' -print -quit | grep -q . || return 0
  if [ "$TIME_STYLE" = v ]; then
    awk_prog='
      /Elapsed \(wall clock\) time/ {
        t = $NF
        n = split(t, a, ":")
        if (n == 3) sec = a[1]*3600 + a[2]*60 + a[3]
        else if (n == 2) sec = a[1]*60 + a[2]
        else sec = t + 0
        wall += sec; nw++
      }
      /^[ \t]*User time \(seconds\):/ { user += $NF; nu++ }
      /^[ \t]*System time \(seconds\):/ { sys += $NF; ns++ }
      END {
        if (nw < 1) exit
        printf "  time: avg wall %.3fs/%s", wall/nw, unit
        if (nu == nw) printf ", user %.3fs", user/nu
        if (ns == nw) printf ", sys %.3fs", sys/ns
        printf "  (%d %ss, /usr/bin/time -v)\n", nw, unit
      }'
  else
    awk_prog='
      /^real / { wall += $2; nw++ }
      /^user / { user += $2; nu++ }
      /^sys /  { sys  += $2; ns++ }
      END {
        if (nw < 1) exit
        printf "  time: avg wall %.3fs/%s", wall/nw, unit
        if (nu == nw) printf ", user %.3fs", user/nu
        if (ns == nw) printf ", sys %.3fs", sys/ns
        printf "  (%d %ss, /usr/bin/time -p; -v unavailable)\n", nw, unit
      }'
  fi
  find "$tdir" -type f -name '*.time' -print0 |
    while IFS= read -r -d '' f; do
      cat "$f"
    done |
    awk -v unit="$unit" "$awk_prog" >&2
}

run_sketch_work() {
  local name="$1" sk_args="$2" s="$3"
  local sk="$CACHE/$name/$s.gdiff"
  # shellcheck disable=SC2086
  "$GDIFF" --num-threads 1 sketch -o "$sk" $sk_args -i "$(fa "$s")" >/dev/null 2>&1
}

# Time one genome sketch; $4 = optional progress marker.
run_sketch() {
  local name="$1" sk_args="$2" s="$3" marker="${4:-}"
  local tlog="$CACHE/$name/times_sketch/${s}.time"
  mkdir -p "$CACHE/$name/times_sketch"
  export -f fa run_sketch_work
  export DIR SUF GDIFF CACHE
  if [ "$TIME_STYLE" = v ]; then
    /usr/bin/time -v -o "$tlog" bash -c 'run_sketch_work "$@"' _ \
      "$name" "$sk_args" "$s"
  else
    /usr/bin/time -p -o "$tlog" bash -c 'run_sketch_work "$@"' _ \
      "$name" "$sk_args" "$s"
  fi
  [ -n "$marker" ] && touch "$marker"
}

# Core pair work (timed by run_pair via /usr/bin/time).
run_pair_work() {
  local name="$1" setup="$2" dist_args="$3" q="$4" s="$5"
  local sum_fwd sum_rev samples_fwd samples_rev row
  sum_fwd="$CACHE/$name/${q}__${s}.summary"
  sum_rev="$CACHE/$name/${s}__${q}.summary"
  samples_fwd="$OUT/gdiff-samples/$name/${q}__${s}.tsv"
  samples_rev="$OUT/gdiff-samples/$name/${s}__${q}.tsv"
  row="$CACHE/$name/rows/${q}__${s}.row"

  # shellcheck disable=SC2086
  "$GDIFF" --num-threads 1 dist "$(fa "$q")" "$CACHE/$name/$s.gdiff" \
    $dist_args -o "$sum_fwd" >/dev/null 2>&1
  # shellcheck disable=SC2086
  "$GDIFF" --num-threads 1 dist "$(fa "$s")" "$CACHE/$name/$q.gdiff" \
    $dist_args -o "$sum_rev" >/dev/null 2>&1
  if [ "$SAMPLES" = 1 ]; then
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 dist "$(fa "$q")" "$CACHE/$name/$s.gdiff" \
      $dist_args --output-samples -o "$samples_fwd" >/dev/null 2>&1
    # shellcheck disable=SC2086
    "$GDIFF" --num-threads 1 dist "$(fa "$s")" "$CACHE/$name/$q.gdiff" \
      $dist_args --output-samples -o "$samples_rev" >/dev/null 2>&1
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

# Run one pair under /usr/bin/time; optional $6 = progress done-marker.
run_pair() {
  local name="$1" setup="$2" dist_args="$3" q="$4" s="$5" marker="${6:-}"
  local tlog="$CACHE/$name/times/${q}__${s}.time"
  mkdir -p "$CACHE/$name/times"
  # /usr/bin/time cannot invoke a shell function directly; re-enter this script's
  # functions via an exported-function child bash.
  export -f fa run_pair_work
  export DIR SUF GDIFF CACHE OUT SAMPLES DIST_COL
  if [ "$TIME_STYLE" = v ]; then
    /usr/bin/time -v -o "$tlog" bash -c 'run_pair_work "$@"' _ \
      "$name" "$setup" "$dist_args" "$q" "$s"
  else
    /usr/bin/time -p -o "$tlog" bash -c 'run_pair_work "$@"' _ \
      "$name" "$setup" "$dist_args" "$q" "$s"
  fi
  [ -n "$marker" ] && touch "$marker"
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
  pdir="$CACHE/$name/progress_sketch"
  rm -rf "$pdir" "$CACHE/$name/times_sketch"
  mkdir -p "$pdir" "$CACHE/$name/times_sketch"
  todo=0
  while read -r s; do
    sk="$CACHE/$name/$s.gdiff"
    [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
    todo=$((todo + 1))
  done < "$CACHE/$name/genomes.txt"
  if [ "$todo" -eq 0 ]; then
    echo "  sketch: all cached" >&2
  else
    draw_progress 0 "$todo" "  sketch"
    n=0
    while read -r s; do
      sk="$CACHE/$name/$s.gdiff"
      [ "$FORCE" != 1 ] && [ -s "$sk" ] && continue
      run_sketch "$name" "$sk_args" "$s" "$pdir/$s" &
      n=$((n + 1))
      throttle "$n" "$todo" "  sketch" "$pdir"
    done < "$CACHE/$name/genomes.txt"
    finish_progress "$pdir" "$todo" "  sketch"
    summarize_times "$CACHE/$name/times_sketch" "genome"
  fi

  # pairs as parallel jobs; each writes CACHE/$name/rows/<q>__<s>.row
  pdir="$CACHE/$name/progress_dist"
  rm -rf "$pdir" "$CACHE/$name/times"
  mkdir -p "$pdir" "$CACHE/$name/times"
  todo=$(awk '$1 != $2 { n++ } END { print n+0 }' "$CACHE/pairs.tsv")
  if [ "$todo" -eq 0 ]; then
    echo "  dist: nothing to do" >&2
  else
    draw_progress 0 "$todo" "  dist"
    n=0
    while read -r q s _; do
      [ "$q" = "$s" ] && continue
      run_pair "$name" "$setup" "$dist_args" "$q" "$s" "$pdir/${q}__${s}" &
      n=$((n + 1))
      throttle "$n" "$todo" "  dist" "$pdir"
    done < "$CACHE/pairs.tsv"
    finish_progress "$pdir" "$todo" "  dist"
    summarize_times "$CACHE/$name/times" "pair"
  fi

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
      sum_xy="$OUT/gdiff-samples/$name/${q}__${s}.tsv"
      sum_yx="$OUT/gdiff-samples/$name/${s}__${q}.tsv"
      if [ -s "$sum_xy" ]; then
        awk -v name="$name" -v q="$q" -v s="$s" \
          'NF>0 { print name "\t" q "\t" s "\t" $0 }' "$sum_xy" >> "$concat"
      fi
      if [ -s "$sum_yx" ]; then
        awk -v name="$name" -v q="$s" -v s="$q" \
          'NF>0 { print name "\t" q "\t" s "\t" $0 }' "$sum_yx" >> "$concat"
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
