#!/usr/bin/env bash
# skani ANI over pairs.tsv (triangle + filter).
# usage: skani.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=16 FORCE=0 ONLY=all
# out: distances/skani-<cfg>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/skani"; DIST_DIR="$OUT/distances"
mkdir -p "$CACHE" "$DIST_DIR"
THREADS="${THREADS:-16}"; JOBS="${JOBS:-$THREADS}"; FORCE="${FORCE:-0}"; ONLY="${ONLY:-all}"
command -v skani >/dev/null || { echo "missing: skani" >&2; exit 1; }

CONFIGS=(
  "default|c=125,m=1000|-c 125 -m 1000"
  "fast|c=200,m=1000,fast|--fast"
  "slow|c=30,m=1000,slow|--slow"
  "sensitive|c=70,robust,min-af=5|-c 70 --robust --min-af 5"
  "fsensitive|c=30,m=100,robust,min-af=0|-c 30 -m 100 --robust --min-af 0"
)
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct'
load_pairs
while read -r id; do fa "$id"; done < "$CACHE/genomes.txt" > "$CACHE/genomes.fa.list"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$DIST_DIR/skani-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then echo "$name: skip"; continue; fi
  echo "$name [$setup]" >&2
  tmp="$(mktemp)"
  # shellcheck disable=SC2086
  skani triangle -t "$THREADS" -l "$CACHE/genomes.fa.list" -E $flags -o "$tmp"
  { echo "$HDR"
    awk -v setup="$setup" -v sfx="$SUFFIX" 'BEGIN{FS=OFS="\t"; sl=length(sfx)}
      FNR==NR { P[$1 FS $2]=1; next }
      FNR>1 && NF>=5 {
        a=$1; b=$2; sub(".*[/]","",a); sub(".*[/]","",b)
        if (sl && substr(a,length(a)-sl+1)==sfx) a=substr(a,1,length(a)-sl)
        if (sl && substr(b,length(b)-sl+1)==sfx) b=substr(b,1,length(b)-sl)
        if ((a FS b) in P)      print "skani",setup,a,b,$3,$4,$5
        else if ((b FS a) in P) print "skani",setup,b,a,$3,$4,$5
      }' "$CACHE/pairs.tsv" "$tmp"
  } > "$tsv"
  rm -f "$tmp"
done
emit_all_distances skani "$HDR" "${CONFIGS[@]}"
echo "done -> $DIST_DIR" >&2
