#!/usr/bin/env bash
# blastn block TSVs for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# Args are passed through to pairwise-mapping/blastn_blocks.py
# usage: blast.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=8 FORCE=0 ONLY=default (name(s), or "all")
# writes: <outdir>/blocks/blast/<cfg>/<q>__<s>.tsv, cache in <outdir>/cache/blast/
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/blast"; OUTDIR="$OUT/blast-blocks"
mkdir -p "$CACHE" "$OUTDIR"
THREADS="${THREADS:-32}"; FORCE="${FORCE:-0}"
command -v blastn >/dev/null || { echo "missing: blastn" >&2; exit 1; }

CONFIGS=(
  "win300|W=300,by-subject|-W 300 --by-subject"
  "win1000|W=1000,by-subject|-W 1000 --by-subject"
  "maxt5000|max-target-seqs=5000|--max-target-seqs 5000"
  "fast|w=11,e=10|--word-size 11 --evalue 10"
)
ONLY="${ONLY:-default}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

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
    python3 "$HERE/../pairwise-mapping/blastn_blocks.py" \
      -q "$(fa "$q")" -s "$(fa "$s")" -o "$out" -t 1 \
      --workdir "$CACHE/.work_${name}_${q}__${s}" $args &
    # Optional: Limit concurrency manually (e.g., maximum 4 background jobs)
    if [[ $(jobs -r -p | wc -l) -ge ${THREADS} ]]; then
        wait -n # Wait for any single background job to finish before continuing
    fi
  done < "$CACHE/pairs.tsv"
done
echo "done -> $OUTDIR"
