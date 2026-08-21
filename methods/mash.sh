#!/usr/bin/env bash
# Mash ANI over pairs.tsv (triangle + filter). Cache = sketches only.
# usage: mash.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=16 FORCE=0 ONLY=all
# out: distances/mash-<cfg>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache mash
DIST_DIR="$OUT/distances"
mkdir -p "$DIST_DIR"
THREADS="${THREADS:-16}"; JOBS="${JOBS:-$THREADS}"; FORCE="${FORCE:-0}"; ONLY="${ONLY:-all}"
command -v mash >/dev/null || { echo "missing: mash" >&2; exit 1; }

CONFIGS=(
  "default|k=21,s=1000|-k 21 -s 1000"
  "large-sketch|k=21,s=10000|-k 21 -s 10000"
  "long-k|k=31,s=10000|-k 31 -s 10000"
  "sensitive|k=16,s=10000|-k 16 -s 10000"
  "fsensitive|k=16,s=50000|-k 16 -s 50000"
)
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tp_value\tshared_hashes\tani_pct'
load_pairs
while read -r id; do fa "$id"; done < "$CACHE/genomes.txt" > "$CACHE/genomes.fa.list"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$DIST_DIR/mash-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then echo "$name: skip"; continue; fi
  echo "$name [$setup]" >&2
  # shellcheck disable=SC2086
  mash sketch -p "$THREADS" $flags -o "$CACHE/$name.msh" -l "$CACHE/genomes.fa.list"
  tmp="$(mktemp)"
  # shellcheck disable=SC2086
  mash triangle -p "$THREADS" $flags -E -l "$CACHE/genomes.fa.list" > "$tmp"
  { echo "$HDR"
    awk -v setup="$setup" -v sfx="$SUFFIX" 'BEGIN{FS=OFS="\t"; sl=length(sfx)}
      FNR==NR { P[$1 FS $2]=1; next }
      NF>=5 {
        a=$1; b=$2; sub(".*[/]","",a); sub(".*[/]","",b)
        if (sl && substr(a,length(a)-sl+1)==sfx) a=substr(a,1,length(a)-sl)
        if (sl && substr(b,length(b)-sl+1)==sfx) b=substr(b,1,length(b)-sl)
        if ((a FS b) in P)      print "mash",setup,a,b,$3,$4,$5,(1-$3)*100
        else if ((b FS a) in P) print "mash",setup,b,a,$3,$4,$5,(1-$3)*100
      }' "$CACHE/pairs.tsv" "$tmp"
  } > "$tsv"
  rm -f "$tmp"
done
emit_all_distances mash "$HDR" "${CONFIGS[@]}"
echo "done -> $DIST_DIR" >&2
