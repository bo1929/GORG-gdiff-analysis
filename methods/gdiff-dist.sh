#!/usr/bin/env bash
# usage: gdiff-dist.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: gdiff JOBS=8 FORCE=0 ONLY=all SAMPLES=1
# out: distances/gdiff-<cfg>.tsv (+ all_gdiff.tsv) [samples/gdiff-<cfg>.tsv if SAMPLES=1]
#
# Target the bundled `gdiff` (bin/gdiff dispatches to gdiff-osx/-x86), v0.2.0.
# CLI mapping vs. the old gdiff2 (`sketch2`/`dist2`):
#   sketch --input-list <list> ... -o all.gdsk      # one record per genome
#   dist   all.gdsk [--output-samples] -o samples.tsv
# v0.2.0 summary header:
#   genome_a genome_b d d_median d_mean d_upper d_highest d_ab d_ba
#   n_ab n_ba n_ub n_na n_filtered   (`d` = reconciled, lr-filtered distance)
# v0.2.0 --output-samples header:
#   config genome_a genome_b seq start end strand direction d lr_ub
#   (genome_a/genome_b canonical, `direction` ab|ba = which genome was query)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./output}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache gdiff
SAMP_DIR="$OUT/samples"; mkdir -p "$SAMP_DIR"
DIST_DIR="$OUT/distances"; mkdir -p "$DIST_DIR"
# bin/gdiff dispatches to the bundled build for this OS/arch.
gdiff="${gdiff:-$REPO_ROOT/bin/gdiff}"
JOBS="${JOBS:-${THREADS:-32}}"; FORCE="${FORCE:-1}"
SAMPLES="${SAMPLES:-1}"; ONLY="${ONLY:-all}"
[ -x "$gdiff" ] || { echo "set gdiff=/path/to/gdiff" >&2; exit 1; }

CONFIGS=(
  "xyz|k=23,w=23,h=11,frac=0.33,-l=500,n=1000|-k 23 -h 11 -w 23 --frac 0.33 -l 500 --sample-size 1000|--hdist-th 4"
  # "abc|k=23,w=23,frac=0.5,-l=500,n=1000|-k 23 -w 23 --frac 0.5 -l 500 --sample-size 1000|--hdist-th 4"
  # "xyzs|k=23,w=23,h=11,frac=0.33,-l=500,n=1000|-k 23 -h 11 -w 23 --frac 0.33 -l 500 --sample-size 1000|--hdist-th 3"
  # "fgh|k=23,w=23,frac=0.50,-l=1000,n=1000|-k 23 -w 23 --frac 0.50 -l 1000 --sample-size 1000|--hdist-th 3"
  # "klm|k=23,w=23,frac=0.66,-l=250,n=1000|-k 23 -w 23 --frac 0.66 -l 250 --sample-size 1000|--hdist-th 3"
  # "abcs|k=23,w=23,frac=0.5,-l=500,n=1000|-k 23 -w 23 --frac 0.5 -l 500 --sample-size 1000|--hdist-th 3"
)
SAMPLES_HEADER=$'config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg\tlr_ub'
DIST_HEADER=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct'
load_pairs

NG="$(wc -l < "$CACHE/genomes.txt" | tr -d ' ')"
[ "$NG" -ge 2 ] || { echo "fewer than 2 genomes in $PAIRS_FILE" >&2; exit 1; }

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup sk_args dist_args <<< "$c"
  want "$name" || continue
  samp_tsv="$SAMP_DIR/gdiff-$name.tsv"
  [ "$FORCE" != 1 ] && [ -s "$samp_tsv" ] && { echo "$name: skip"; continue; }
  echo "$name [$setup] jobs=$JOBS" >&2
  skdir="$CACHE/$name"; mkdir -p "$skdir"
  parts="$(mktemp -d)"
  bundle="$skdir/all.gdsk"

  # Rebuild the bundle only when forced or when the genome set changed.
  rebuild=0
  [ "$FORCE" = 1 ] && rebuild=1
  [ -s "$bundle" ] || rebuild=1
  if [ -s "$skdir/genomes.txt" ] && ! cmp -s "$CACHE/genomes.txt" "$skdir/genomes.txt"; then rebuild=1; fi
  if [ ! -s "$skdir/genomes.txt" ]; then rebuild=1; fi

  if [ "$rebuild" -eq 0 ]; then
    echo "  sketch: cached" >&2
  else
    rm -f "$bundle"
    while read -r g; do printf '%s\t%s\n' "$g" "$(fa "$g")"; done < "$CACHE/genomes.txt" > "$parts/input.list"
    echo "  sketch [$NG genomes, threads=$JOBS] -> $bundle" >&2
    # shellcheck disable=SC2086
    "$gdiff" --num-threads "$JOBS" sketch --input-list "$parts/input.list" $sk_args -o "$bundle" || {
      echo "sketch failed for $name" >&2; rm -rf "$parts"; exit 1; }
    cp "$CACHE/genomes.txt" "$skdir/genomes.txt"
  fi

  # --- single dist: within-mode over the whole bundle, internal threads ---
  echo "  dist [$((NG * (NG - 1) / 2)) pairs, threads=$JOBS]" >&2

  # Canonical unordered pair keys (what the v1 script emitted from this pairs
  # file); the full all-vs-all matrix is exact, subsets get filtered out. The
  # distance arguments (--hdist-th ...) are passed through to the dist call.
  awk 'NF>=2{a=$1;b=$2;if(a<b)k=a"|"b;else k=b"|"a;if(!(k in S)){S[k]=1;print k}}' \
    "$CACHE/pairs.tsv" | sort -u > "$parts/canon"
  if [ ! -s "$parts/canon" ]; then rm -rf "$parts"; continue; fi

  # --- distances: one run without --output-samples, filtered to the pair list ---
  echo "  distances -> $DIST_DIR/gdiff-$name.tsv" >&2
  # shellcheck disable=SC2086
  "$gdiff" --num-threads "$JOBS" dist "$bundle" $dist_args > "$parts/sum.tsv" || true
  {
    echo "$DIST_HEADER"
    # Summary schema (v0.2.0): genome_a genome_b d d_median ... The reconciled,
    # lr-filtered per-pair distance is column `d`.
    awk -F'\t' -v setup="$setup" -v f="$parts/canon" '
      BEGIN { OFS="\t"; while ((getline l < f) > 0) fl[l] = 1; close(f) }
      /^#/ { next }
      $1 == "genome_a" { for (i = 1; i <= NF; i++) if ($i == "d") { col = i; break } ; next }
      col && NF >= col {
        a = $1; b = $2
        k = (a < b) ? a "|" b : b "|" a
        if (!(k in fl)) next
        printf "gdiff\t%s\t%s\t%s\t%s\t%.6g\n", setup, a, b, $col, (1 - $col) * 100
      }' "$parts/sum.tsv"
  } > "$DIST_DIR/gdiff-$name.tsv"

  if [ "$SAMPLES" = 1 ]; then
    mkdir -p "$SAMP_DIR"
    echo "  samples -> $SAMP_DIR/gdiff-$name.tsv" >&2
    # shellcheck disable=SC2086
    "$gdiff" --num-threads "$JOBS" dist "$bundle" $dist_args --output-samples -o "$parts/samp.tsv" || true
    {
      echo "$SAMPLES_HEADER"
      # v0.2.0 sample schema (header row): config genome_a genome_b seq start
      # end strand direction d lr_ub. genome_a/genome_b are canonical and
      # `direction` (ab|ba) says which was the query; emit one row per window
      # with genome_a = query so the pair key matches pairs.tsv. lr_bg is not
      # reported by v0.2.0 -> ".".
      awk -F'\t' -v cfg="$name" -v f="$parts/canon" '
        BEGIN { OFS="\t"; while ((getline l < f) > 0) fl[l] = 1; close(f) }
        !/^#/ && $1 != "config" && NF >= 10 {
          k = ($2 < $3) ? $2 "|" $3 : $3 "|" $2
          if (!(k in fl)) next
          if ($8 == "ba") { ga = $3; gb = $2; ref = $2 }   # query = genome_b
          else            { ga = $2; gb = $3; ref = $3 }   # query = genome_a
          print cfg, ga, gb, $4, $5, $6, $7, ref, $9, ".", $10
        }' "$parts/samp.tsv"
    } > "$samp_tsv"
  fi

  rm -rf "$parts"
done

emit_all_distances gdiff "$DIST_HEADER" "${CONFIGS[@]}"
echo "done -> $DIST_DIR" >&2
