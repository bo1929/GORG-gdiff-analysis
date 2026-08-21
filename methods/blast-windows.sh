#!/usr/bin/env bash
# blastn window/block TSVs over pairs.tsv. No persistent workdirs.
# usage: blast-windows.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: JOBS=32 FORCE=0 ONLY=win300
# out: blocks/blast/<cfg>/<q>__<s>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache blast
OUTDIR="$OUT/blocks/blast"
mkdir -p "$OUTDIR"
JOBS="${JOBS:-${THREADS:-32}}"; FORCE="${FORCE:-0}"; ONLY="${ONLY:-win300}"
command -v blastn >/dev/null || { echo "missing: blastn" >&2; exit 1; }

CONFIGS=(
  "win300|W=300,by-subject|-W 300 --by-subject"
  "win1000|W=1000,by-subject|-W 1000 --by-subject"
  "maxt5000|max-target-seqs=5000|--max-target-seqs 5000"
  "fast|w=11,e=10|--word-size 11 --evalue 10"
)
load_pairs
PY="$HERE/../pairwise-mapping/blastn_blocks.py"

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
