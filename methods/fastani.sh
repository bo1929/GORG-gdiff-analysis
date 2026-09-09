#!/usr/bin/env bash
# FastANI ANI over pairs.tsv. Asymmetric tool: each pair is emitted BOTH ways
# (a->b, b->a). SAMPLES=1 (default) also writes gdiff-schema per-window samples
# (fastANI --visualize) to samples/fastani-<cfg>.tsv.
# usage: fastani.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: FASTANI JOBS=8 THREADS=4 FORCE=0 ONLY=all SAMPLES=1 CACHE_ROOT
# out: distances/fastani-<cfg>.tsv  [samples/fastani-<cfg>.tsv]  all_fastani.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache fastani
DIST_DIR="$OUT/distances"; SAMP_DIR="$OUT/samples"; mkdir -p "$DIST_DIR"
JOBS="${JOBS:-${THREADS:-16}}"; THREADS="${THREADS:-2}"
FORCE="${FORCE:-1}"; ONLY="${ONLY:-all}"; SAMPLES="${SAMPLES:-1}"

FASTANI="${FASTANI:-}"
[ -z "$FASTANI" ] && command -v fastANI >/dev/null 2>&1 && FASTANI="$(command -v fastANI)"
[ -z "$FASTANI" ] && command -v micromamba >/dev/null 2>&1 \
  && FASTANI="$(micromamba run -n base which fastANI 2>/dev/null | tail -1)"
[ -x "${FASTANI:-}" ] || { echo "set FASTANI=/path/to/fastANI (micromamba install -n base -c bioconda fastani)" >&2; exit 1; }

CONFIGS=(
  "frag1000|frag=1000,minfrac=0.1|--fragLen 1000 --minFraction 0.1"
  "frag3000|frag=3000,minfrac=0.1|--fragLen 3000 --minFraction 0.1"
)
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct'
SAMP_HDR=$'config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg\tlr_ub'
load_pairs

# run_pair <cfg> <setup> <flags> <q> <s> <dir> <samp>
# Distance rows (a->b then b->a) on stdout; if <samp> non-empty, gdiff-style
# sample rows append to it. fastANI reports ANI only (af_ref/af_query = '.');
# visual cols 7,8 are 0-based window -> 1-based, d = 1 - ANI/100.
run_pair() {
  local cfg=$1 setup=$2 flags=$3 q=$4 s=$5 dir=$6 samp=$7
  local ab="$dir.ab" ba="$dir.ba" ani_ab ani_ba
  # shellcheck disable=SC2086
  "$FASTANI" -q "$(fa "$q")" -r "$(fa "$s")" -o "$ab" $flags -t "$THREADS" --visualize >/dev/null 2>&1 || return
  # shellcheck disable=SC2086
  "$FASTANI" -q "$(fa "$s")" -r "$(fa "$q")" -o "$ba" $flags -t "$THREADS" --visualize >/dev/null 2>&1 || return
  ani_ab=$(awk -F'\t' '$3+0==$3 {print $3; exit}' "$ab")
  ani_ba=$(awk -F'\t' '$3+0==$3 {print $3; exit}' "$ba")
  [ -n "$ani_ab" ] && printf 'fastani\t%s\t%s\t%s\t%s\t.\t.\n' "$setup" "$q" "$s" "$ani_ab"
  [ -n "$ani_ba" ] && printf 'fastani\t%s\t%s\t%s\t%s\t.\t.\n' "$setup" "$s" "$q" "$ani_ba"
  if [ "$samp" ]; then
    [ -e "$ab.visual" ] && awk -F'\t' -v c="$cfg" -v a="$q" -v b="$s" '
      $3+0==$3 { printf "%s\t%s\t%s\t.\t%d\t%d\t.\t%s\t%.6f\t.\t.\n", c, a, b, $7+1, $8+1, b, 1-$3/100 }' \
      "$ab.visual" >> "$samp"
    [ -e "$ba.visual" ] && awk -F'\t' -v c="$cfg" -v a="$s" -v b="$q" '
      $3+0==$3 { printf "%s\t%s\t%s\t.\t%d\t%d\t.\t%s\t%.6f\t.\t.\n", c, a, b, $7+1, $8+1, b, 1-$3/100 }' \
      "$ba.visual" >> "$samp"
  fi
  rm -f "$ab" "$ba" "$ab.visual" "$ba.visual"
}

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$DIST_DIR/fastani-$name.tsv"
  [ "$FORCE" != 1 ] && [ -s "$tsv" ] && { echo "$name: skip"; continue; }
  echo "$name [$setup] jobs=$JOBS" >&2
  parts=$(mktemp -d); : > "$parts/out"; : > "$parts/samp"
  [ "$SAMPLES" = 1 ] && { mkdir -p "$SAMP_DIR"; echo "$SAMP_HDR" > "$SAMP_DIR/fastani-$name.tsv"; }
  np=$(awk '$1!=$2{n++} END{print n+0}' "$CACHE/pairs.tsv")
  echo "$HDR" > "$tsv"
  progress 0 "$np" fastani
  n=0; bi=0
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    sp=""; [ "$SAMPLES" = 1 ] && sp="$parts/s.$bi"
    ( run_pair "$name" "$setup" "$flags" "$q" "$s" "$parts/d.$bi" "$sp" >> "$parts/p.$bi" ) &
    n=$((n+1)); bi=$((bi+1))
    if (( n % JOBS == 0 )); then
      wait
      cat "$parts"/p.* >> "$parts/out" 2>/dev/null || true
      cat "$parts"/s.* >> "$parts/samp" 2>/dev/null || true
      rm -f "$parts"/p.* "$parts"/s.* "$parts"/d.*; bi=0
      progress "$n" "$np" fastani
    fi
  done < "$CACHE/pairs.tsv"
  wait
  cat "$parts"/p.* >> "$parts/out" 2>/dev/null || true
  cat "$parts"/s.* >> "$parts/samp" 2>/dev/null || true
  progress "$np" "$np" fastani
  cat "$parts/out" >> "$tsv"
  [ "$SAMPLES" = 1 ] && cat "$parts/samp" >> "$SAMP_DIR/fastani-$name.tsv"
  rm -rf "$parts"
  echo "  $(tail -n +2 "$tsv" | wc -l | tr -d ' ') rows -> $tsv" >&2
done
emit_all_distances fastani "$HDR" "${CONFIGS[@]}"
echo "done -> $DIST_DIR" >&2
