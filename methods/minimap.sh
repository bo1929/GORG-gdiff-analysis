#!/usr/bin/env bash
# minimap2 block TSVs for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# Args are passed through to pairwise-mapping/minimap_blocks.py.
# usage: minimap.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=8 FORCE=0 ONLY=default (name(s), or "all")
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./minimap_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT/blocks"; OUT="$(cd "$OUT" && pwd)"
THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"
command -v minimap2 >/dev/null || { echo "missing: minimap2" >&2; exit 1; }

CONFIGS=(
  "default|asm20|"
  "sensitive|k=15,w=5,m=10|--sensitive"
  "asm5|preset=asm5|--preset asm5"
  "minid85|sensitive,min-id=85|--sensitive --min-identity 85"
)
ONLY="${ONLY:-default}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$OUT/pairs.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup args <<< "$c"
  want "$name" || continue
  echo "$name [$setup]"
  mkdir -p "$OUT/blocks/$name"
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    out="$OUT/blocks/$name/${q}__${s}.tsv"
    [ "$FORCE" != 1 ] && [ -s "$out" ] && continue
    # shellcheck disable=SC2086
    python3 "$HERE/../pairwise-mapping/minimap_blocks.py" \
      -q "$(fa "$q")" -s "$(fa "$s")" -o "$out" -t "$THREADS" \
      --workdir "$OUT/blocks/$name/.work_${q}__${s}" $args
  done < "$OUT/pairs.tsv"
done
echo "done -> $OUT/blocks"
