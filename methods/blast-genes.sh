#!/usr/bin/env bash
# Gene-level blastn: GBK genes vs subject FASTA. No persistent workdirs.
# usage: blast-genes.sh <gbk_dir> <fasta_dir> <pairs.tsv> [outdir] [gbk_suf] [fa_suf]
# env: JOBS=32 FORCE=0 ONLY=default
# out: blocks/blast-genes/<cfg>/<q>__<s>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GBK_DIR="$(cd "${1:?usage: $0 <gbk_dir> <fasta_dir> <pairs.tsv> [outdir] [gbk_suf] [fa_suf]}" && pwd)"
FA_DIR="$(cd "${2:?}" && pwd)"
PAIRS_FILE="${3:?}"
OUT="${4:-./methods_out}"; GBK_SUF="${5:-.gbk}"; FA_SUF="${6:-_contigs.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache blast-genes
OUTDIR="$OUT/blocks/blast-genes"
mkdir -p "$OUTDIR"
JOBS="${JOBS:-${THREADS:-32}}"; FORCE="${FORCE:-0}"; ONLY="${ONLY:-default}"
command -v blastn >/dev/null || { echo "missing: blastn" >&2; exit 1; }

CONFIGS=(
)
grep -v '^#' "$PAIRS_FILE" | awk 'NF>=2' > "$CACHE/pairs.tsv"
PY="$HERE/../blastn_from_gbk.py"
gbk() { local f="$GBK_DIR/$1$GBK_SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
fa_s() { local f="$FA_DIR/$1$FA_SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }

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
    ( python3 "$PY" -g "$(gbk "$q")" -s "$(fa_s "$s")" -o "$out" -t 1 $args >/dev/null 2>&1 ) &
    n=$((n + 1)); throttle "$n"
  done < "$CACHE/pairs.tsv"
  wait
done
echo "done -> $OUTDIR" >&2
