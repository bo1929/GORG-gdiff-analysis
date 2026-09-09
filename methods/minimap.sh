#!/usr/bin/env bash
# minimap2 block TSVs over pairs.tsv. No persistent workdirs.
# usage: minimap.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: JOBS=8 THREADS=1 FORCE=0 ONLY=default
# out: blocks/minimap/<cfg>/<q>__<s>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache minimap
OUTDIR="$OUT/blocks/minimap"
mkdir -p "$OUTDIR"
JOBS="${JOBS:-${THREADS:-8}}"; FORCE="${FORCE:-0}"; ONLY="${ONLY:-default}"
command -v minimap2 >/dev/null || { echo "missing: minimap2" >&2; exit 1; }

CONFIGS=(
)
load_pairs
PY="$HERE/../pairwise-mapping/minimap_blocks.py"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup args <<< "$c"
  want "$name" || continue
  echo "$name [$setup] jobs=$JOBS" >&2
  mkdir -p "$OUTDIR/$name"
  n=0
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    out="$OUTDIR/$name/${q}__${s}.tsv"
    [ "$FORCE" != 1 ] && [ -s "$out" ] && continue
    # shellcheck disable=SC2086
    ( python3 "$PY" -q "$(fa "$q")" -s "$(fa "$s")" -o "$out" -t 1 $args >/dev/null 2>&1 ) &
    n=$((n + 1)); throttle "$n"
  done < "$CACHE/pairs.tsv"
  wait
done
echo "done -> $OUTDIR" >&2
