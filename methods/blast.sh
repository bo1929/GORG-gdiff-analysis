#!/usr/bin/env bash
# blastn block TSVs for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# Args are passed through to pairwise-mapping/blastn_blocks.py
# (default is already ultra-sensitive: -word_size 7 -evalue 1000).
# usage: blast.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=8 FORCE=0 ONLY=default (name(s), or "all")
# writes: <outdir>/blast/{<cfg>/<q>__<s>.tsv, cache/}
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
MDIR="$OUT/blast"
mkdir -p "$MDIR/cache"; OUT="$(cd "$OUT" && pwd)"; MDIR="$OUT/blast"
THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"
command -v blastn >/dev/null || { echo "missing: blastn" >&2; exit 1; }

CONFIGS=(
  "default|w=7,e=1000|"
  "maxtargets|max-target-seqs=5000|--max-target-seqs 5000"
  "fast|w=11,e=10|--word-size 11 --evalue 10"
  "window1k|W=1000,by-subject|-W 1000 --by-subject"
)
ONLY="${ONLY:-default}"

fa()   { local f="$DIR/$1$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"; }
want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$MDIR/cache/pairs.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup args <<< "$c"
  want "$name" || continue
  echo "$name [$setup]"
  mkdir -p "$MDIR/$name"
  while read -r q s _; do
    [ "$q" = "$s" ] && continue
    out="$MDIR/$name/${q}__${s}.tsv"
    [ "$FORCE" != 1 ] && [ -s "$out" ] && continue
    # shellcheck disable=SC2086
    python3 "$HERE/../pairwise-mapping/blastn_blocks.py" \
      -q "$(fa "$q")" -s "$(fa "$s")" -o "$out" -t "$THREADS" \
      --workdir "$MDIR/cache/.work_${name}_${q}__${s}" $args
  done < "$MDIR/cache/pairs.tsv"
done
echo "done -> $MDIR"
