#!/usr/bin/env bash
# Dashing2 (MinHash) ANI over pairs.tsv, one row per pair in its own
# orientation (genome_a = first id). Configs mirror simulations/run_dashing2.sh;
# containment/symmetric emit a fraction, ANI = 100*(1 + ln(frac)/k) as there.
# usage: dashing2.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: DASHING2 JOBS=8 THREADS=4 FORCE=0 ONLY=all CACHE_ROOT
# out: distances/dashing2-<cfg>.tsv  all_dashing2.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache dashing2
DIST_DIR="$OUT/distances"; mkdir -p "$DIST_DIR"
JOBS="${JOBS:-${THREADS:-8}}"; THREADS="${THREADS:-4}"
FORCE="${FORCE:-0}"; ONLY="${ONLY:-all}"

DASHING2="${DASHING2:-./dashing2-s512bw}"
[ -z "$DASHING2" ] && command -v dashing2 >/dev/null 2>&1 && DASHING2="$(command -v dashing2)"
[ -z "$DASHING2" ] && [ -x "$REPO_ROOT/dashing2-s512bw" ] && DASHING2="$REPO_ROOT/dashing2-s512bw"
[ -x "${DASHING2:-}" ] || { echo "set DASHING2=/path/to/dashing2" >&2; exit 1; }

CONFIGS=(
  "v1|mash,k=31,S=1024|--mash-distance -k 31 -S 1024"
  "v2|containment,k=32,S=1024|--containment"
  "v3|symmetric-containment,k=31,s=1024|--symmetric-containment -k 31 -S 1024"
  "v4|symmetric-containment,k=23,s=2048|--symmetric-containment -k 23 -S 2048"
  "v5|mash,k=21,S=5000|-k 21 -S 5000 --mash-distance"
)
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct'
load_pairs

# run_pair <setup> <flags> <q> <s> <dir>
# Emit one distance row on stdout; none when dashing2 reports '-', inf or NaN
# (low-sensitivity modes go blind on distant genomes). Containment is
# asymmetric, so -F is genome_b (ref) and -Q genome_a (query).
run_pair() {
  local setup=$1 flags=$2 q=$3 s=$4 dir=$5 k m=mash v
  printf '%s\n' "$(fa "$s")" > "$dir.rl"
  printf '%s\n' "$(fa "$q")" > "$dir.ql"
  case "$flags" in *containment*) m=containment;; esac
  k=$(printf '%s' "$flags" | sed -n 's/.*-k \([0-9]*\).*/\1/p'); k="${k:-32}"
  # shellcheck disable=SC2086
  v=$("$DASHING2" cmp -F "$dir.rl" -Q "$dir.ql" $flags -p "$THREADS" 2>/dev/null \
        | awk -F'\t' '!/^#/ { print $2; exit }')
  rm -f "$dir.rl" "$dir.ql"
  [ -n "$v" ] || return
  case "$v" in "-"|inf|NaN|nan) return;; esac
  awk -v q="$q" -v s="$s" -v setup="$setup" -v x="$v" -v k="$k" -v m="$m" '
    BEGIN {
      if (m == "containment") {
        if (x <= 0) exit
        ani = 100*(1 + log(x)/k)
        printf "dashing2\t%s\t%s\t%s\t%.6f\t%.6f\n", setup, q, s, 1-ani/100, ani
      } else
        printf "dashing2\t%s\t%s\t%s\t%.6f\t%.6f\n", setup, q, s, x, (1-x)*100
    }'
}

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$DIST_DIR/dashing2-$name.tsv"
  [ "$FORCE" != 1 ] && [ -s "$tsv" ] && { echo "$name: skip"; continue; }
  echo "$name [$setup] jobs=$JOBS" >&2
  parts=$(mktemp -d); : > "$parts/out"
  np=$(awk '$1!=$2{n++} END{print n+0}' "$CACHE/pairs.tsv")
  echo "$HDR" > "$tsv"
  progress 0 "$np" dashing2
  n=0; bi=0
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    ( run_pair "$setup" "$flags" "$q" "$s" "$parts/d.$bi" >> "$parts/p.$bi" ) &
    n=$((n+1)); bi=$((bi+1))
    if (( n % JOBS == 0 )); then
      wait
      cat "$parts"/p.* >> "$parts/out" 2>/dev/null || true
      rm -f "$parts"/p.* "$parts"/d.*; bi=0
      progress "$n" "$np" dashing2
    fi
  done < "$CACHE/pairs.tsv"
  wait
  cat "$parts"/p.* >> "$parts/out" 2>/dev/null || true
  progress "$np" "$np" dashing2
  cat "$parts/out" >> "$tsv"
  rm -rf "$parts"
  echo "  $(tail -n +2 "$tsv" | wc -l | tr -d ' ') rows -> $tsv" >&2
done
emit_all_distances dashing2 "$HDR" "${CONFIGS[@]}"
echo "done -> $DIST_DIR" >&2
