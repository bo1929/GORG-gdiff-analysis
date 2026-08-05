#!/usr/bin/env bash
# Gene-level blastn blocks for the pairs in <pairs.tsv>: all annotated genes
# of each query (GBK) mapped against the subject genome FASTA, one run per
# CONFIGS entry. Args are passed through to blastn_from_gbk.py.
# usage: blast-genes.sh <gbk_dir> <fasta_dir> <pairs.tsv> [outdir] [gbk_suffix] [fasta_suffix]
#   gbk_suffix=.gbk  fasta_suffix=_contigs.fasta
# env: THREADS=32 FORCE=0 ONLY=default (name(s), or "all")
#   Pairs run in parallel as background jobs (up to THREADS at a time,
#   one blast thread each).
# writes: <outdir>/blast-gene-blocks/<cfg>/<q>__<s>.tsv,
#         cache in <outdir>/cache/blast-genes/
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
GBK_DIR="$(cd "${1:?usage: $0 <gbk_dir> <fasta_dir> <pairs.tsv> [outdir] [gbk_suffix] [fasta_suffix]}" && pwd)"
FA_DIR="$(cd "${2:?usage: $0 <gbk_dir> <fasta_dir> <pairs.tsv> [outdir] [gbk_suffix] [fasta_suffix]}" && pwd)"
PAIRS="${3:?usage: $0 <gbk_dir> <fasta_dir> <pairs.tsv> [outdir] [gbk_suffix] [fasta_suffix]}"
OUT="${4:-./methods_out}"; SUF="${5:-.gbk}"; FA_SUF="${6:-_contigs.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/blast-genes"; OUTDIR="$OUT/blast-gene-blocks"
mkdir -p "$CACHE" "$OUTDIR"
THREADS="${THREADS:-32}"; FORCE="${FORCE:-0}"
command -v blastn >/dev/null || { echo "missing: blastn" >&2; exit 1; }

CONFIGS=(
  "default|w=7,e=1000,chain|"
  "all-hits|w=7,e=1000,all|--all-hits"
)
ONLY="${ONLY:-default}"

gbk_for() { local f="$GBK_DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
fa_for()  { local f="$FA_DIR/$1$FA_SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want()    { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$CACHE/pairs.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup args <<< "$c"
  want "$name" || continue
  echo "$name [$setup]"
  mkdir -p "$OUTDIR/$name"
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    out="$OUTDIR/$name/${q}__${s}.tsv"
    [ "$FORCE" != 1 ] && [ -s "$out" ] && continue
    # shellcheck disable=SC2086
    python3 "$HERE/../blastn_from_gbk.py" \
      -g "$(gbk_for "$q")" -s "$(fa_for "$s")" -o "$out" -t 1 \
      --workdir "$CACHE/.work_${name}_${q}__${s}" $args &
    # Limit concurrency: at most THREADS background jobs at a time
    if [[ $(jobs -r -p | wc -l) -ge ${THREADS} ]]; then
        wait -n
    fi
  done < "$CACHE/pairs.tsv"
  wait
done
echo "done -> $OUTDIR"
