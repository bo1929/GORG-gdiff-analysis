#!/usr/bin/env bash
# skani ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# usage: skani.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=8 FORCE=0 ONLY=cfg1,cfg2
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./skani_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT/cache" "$OUT/distances"; OUT="$(cd "$OUT" && pwd)"
THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"

CONFIGS=(
  "default|c=125,m=1000|-c 125 -m 1000"
  "fast|c=200,m=1000,fast|--fast"
  "slow|c=30,m=1000,slow|--slow"
  "sensitive|c=70,robust,min-af=5|-c 70 --robust --min-af 5"
  "fsensitive|c=30,m=100,robust,min-af=0|-c 30 -m 100 --robust --min-af 0"
)
ONLY="${ONLY:-default}"

want() { [ -z "${ONLY:-}" ] || case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$OUT/cache/pairs.tsv"
cut -f1,2 "$OUT/cache/pairs.tsv" | tr '\t' '\n' | sort -u | while read -r id; do
  f="$DIR/$id$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"
done > "$OUT/cache/genomes.txt"
SUF_RE=$(printf '%s' "$SUF" | sed 's/[.[*\^$]/\\&/g')

HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct'
echo "$HDR" > "$OUT/distances/all_skani.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$OUT/distances/$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; tail -n+2 "$tsv" >> "$OUT/distances/all_skani.tsv"; continue
  fi
  echo "$name [$setup]"
  # shellcheck disable=SC2086
  skani triangle -t "$THREADS" -l "$OUT/cache/genomes.txt" -E $flags -o "$OUT/cache/$name.tmp"
  { echo "$HDR"
    awk -v setup="$setup" -v sfx="$SUF_RE" 'BEGIN{FS=OFS="\t"}
      FNR==NR { P[$1 FS $2]=1; next }
      FNR>1 && NF>=5 {
        a=$1; b=$2; sub(".*\/","",a); sub(".*\/","",b); sub(sfx"$","",a); sub(sfx"$","",b)
        if ((a FS b) in P)      print "skani",setup,a,b,$3,$4,$5
        else if ((b FS a) in P) print "skani",setup,b,a,$3,$4,$5
      }' "$OUT/cache/pairs.tsv" "$OUT/cache/$name.tmp"
  } > "$tsv"
  rm "$OUT/cache/$name.tmp"
  tail -n+2 "$tsv" >> "$OUT/distances/all_skani.tsv"
done
echo "done -> $OUT/distances"
