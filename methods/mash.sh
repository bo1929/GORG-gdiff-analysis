#!/usr/bin/env bash
# Mash ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# usage: mash.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=8 FORCE=0 ONLY=cfg1,cfg2
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./mash_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT/cache" "$OUT/distances"; OUT="$(cd "$OUT" && pwd)"
THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"

CONFIGS=(
  "default|21|1000"
  "sketch10k|21|10000"
  "k16_sketch10k|16|10000"
  "k31_sketch10k|31|10000"
  "maxsens|16|50000"
)

want() { [ -z "${ONLY:-}" ] || case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$OUT/cache/pairs.tsv"
cut -f1,2 "$OUT/cache/pairs.tsv" | tr '\t' '\n' | sort -u | while read -r id; do
  f="$DIR/$id$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"
done > "$OUT/cache/genomes.txt"
SUF_RE=$(printf '%s' "$SUF" | sed 's/[.[*\^$]/\\&/g')

HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tp_value\tshared_hashes\tani_pct'
echo "$HDR" > "$OUT/distances/all_mash.tsv"

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name k s <<< "$c"
  want "$name" || continue
  tsv="$OUT/distances/$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; tail -n+2 "$tsv" >> "$OUT/distances/all_mash.tsv"; continue
  fi
  echo "$name [k=$k,s=$s]"
  mash sketch -p "$THREADS" -k "$k" -s "$s" -o "$OUT/cache/$name.msh" -l "$OUT/cache/genomes.txt"
  mash triangle -p "$THREADS" -k "$k" -s "$s" -E -l "$OUT/cache/genomes.txt" > "$OUT/cache/$name.tmp"
  { echo "$HDR"
    awk -v setup="k=$k,s=$s" -v sfx="$SUF_RE" 'BEGIN{FS=OFS="\t"}
      FNR==NR { P[$1 FS $2]=1; next }
      NF>=5 {
        a=$1; b=$2; sub(".*\/","",a); sub(".*\/","",b); sub(sfx"$","",a); sub(sfx"$","",b)
        if ((a FS b) in P)      print "mash",setup,a,b,$3,$4,$5,(1-$3)*100
        else if ((b FS a) in P) print "mash",setup,b,a,$3,$4,$5,(1-$3)*100
      }' "$OUT/cache/pairs.tsv" "$OUT/cache/$name.tmp"
  } > "$tsv"
  rm "$OUT/cache/$name.tmp"
  tail -n+2 "$tsv" >> "$OUT/distances/all_mash.tsv"
done
echo "done -> $OUT/distances"
