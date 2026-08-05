#!/usr/bin/env bash
# skani ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# usage: skani.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=8 FORCE=0 ONLY=default (name(s), or "all")
# writes: <outdir>/distances/{skani-<cfg>.tsv, all_skani.tsv}, cache in <outdir>/cache/skani/
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
CACHE="$OUT/cache/skani"; DISTS="$OUT/distances"
mkdir -p "$CACHE" "$DISTS"
THREADS="${THREADS:-16}"; FORCE="${FORCE:-0}"

CONFIGS=(
  "default|c=125,m=1000|-c 125 -m 1000"
  "fast|c=200,m=1000,fast|--fast"
  "slow|c=30,m=1000,slow|--slow"
  "sensitive|c=70,robust,min-af=5|-c 70 --robust --min-af 5"
  "fsensitive|c=30,m=100,robust,min-af=0|-c 30 -m 100 --robust --min-af 0"
)
ONLY="${ONLY:-all}"

want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$CACHE/pairs.tsv"
cut -f1,2 "$CACHE/pairs.tsv" | tr '\t' '\n' | sort -u | while read -r id; do
  f="$DIR/$id$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"
done > "$CACHE/genomes.txt"

HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct'

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$DISTS/skani-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; continue
  fi
  echo "$name [$setup]"
  # shellcheck disable=SC2086
  skani triangle -t "$THREADS" -l "$CACHE/genomes.txt" -E $flags -o "$CACHE/$name.tmp"
  { echo "$HDR"
    awk -v setup="$setup" -v sfx="$SUF" 'BEGIN{FS=OFS="\t"; slen=length(sfx)}
      FNR==NR { P[$1 FS $2]=1; next }
      FNR>1 && NF>=5 {
        a=$1; b=$2; sub(".*[/]","",a); sub(".*[/]","",b)
        if (substr(a,length(a)-slen+1)==sfx) a=substr(a,1,length(a)-slen)
        if (substr(b,length(b)-slen+1)==sfx) b=substr(b,1,length(b)-slen)
        if ((a FS b) in P)      print "skani",setup,a,b,$3,$4,$5
        else if ((b FS a) in P) print "skani",setup,b,a,$3,$4,$5
      }' "$CACHE/pairs.tsv" "$CACHE/$name.tmp"
  } > "$tsv"
  rm "$CACHE/$name.tmp"
done

{ echo "$HDR"
  for c in "${CONFIGS[@]}"; do
    IFS='|' read -r name _ <<< "$c"
    if want "$name" && [ -s "$DISTS/skani-$name.tsv" ]; then tail -n+2 "$DISTS/skani-$name.tsv"; fi
  done
} > "$DISTS/all_skani.tsv"
echo "done -> $DISTS"
