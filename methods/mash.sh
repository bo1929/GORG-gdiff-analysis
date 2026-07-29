#!/usr/bin/env bash
# Mash ANI for the pairs in <pairs.tsv>, one run per CONFIGS entry.
# usage: mash.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS=8 FORCE=0 ONLY=default (name(s), or "all")
# writes: <outdir>/mash/{mash-<cfg>.tsv, all_mash.tsv, cache/}
set -euo pipefail
DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS="${2:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}"
OUT="${3:-./methods_out}"; SUF="${4:-.fasta}"
MDIR="$OUT/mash"
mkdir -p "$MDIR/cache"; OUT="$(cd "$OUT" && pwd)"; MDIR="$OUT/mash"
THREADS="${THREADS:-8}"; FORCE="${FORCE:-0}"

CONFIGS=(
  "default|k=21,s=1000|-k 21 -s 1000"
  "large-sketch|k=21,s=10000|-k 21 -s 10000"
  "long-k|k=31,s=10000|-k 31 -s 10000"
  "sensitive|k=16,s=10000|-k 16 -s 10000"
  "fsensitive|k=16,s=50000|-k 16 -s 50000"
)
ONLY="${ONLY:-default}"

want() { [ "$ONLY" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

grep -v '^#' "$PAIRS" | awk 'NF>=2' > "$MDIR/cache/pairs.tsv"
cut -f1,2 "$MDIR/cache/pairs.tsv" | tr '\t' '\n' | sort -u | while read -r id; do
  f="$DIR/$id$SUF"; [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }; echo "$f"
done > "$MDIR/cache/genomes.txt"

HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tp_value\tshared_hashes\tani_pct'

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$MDIR/mash-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then
    echo "$name: skip"; continue
  fi
  echo "$name [$setup]"
  # shellcheck disable=SC2086
  mash sketch -p "$THREADS" $flags -o "$MDIR/cache/$name.msh" -l "$MDIR/cache/genomes.txt"
  # shellcheck disable=SC2086
  mash triangle -p "$THREADS" $flags -E -l "$MDIR/cache/genomes.txt" > "$MDIR/cache/$name.tmp"
  { echo "$HDR"
    awk -v setup="$setup" -v sfx="$SUF" 'BEGIN{FS=OFS="\t"; slen=length(sfx)}
      FNR==NR { P[$1 FS $2]=1; next }
      NF>=5 {
        a=$1; b=$2; sub(".*[/]","",a); sub(".*[/]","",b)
        if (substr(a,length(a)-slen+1)==sfx) a=substr(a,1,length(a)-slen)
        if (substr(b,length(b)-slen+1)==sfx) b=substr(b,1,length(b)-slen)
        if ((a FS b) in P)      print "mash",setup,a,b,$3,$4,$5,(1-$3)*100
        else if ((b FS a) in P) print "mash",setup,b,a,$3,$4,$5,(1-$3)*100
      }' "$MDIR/cache/pairs.tsv" "$MDIR/cache/$name.tmp"
  } > "$tsv"
  rm "$MDIR/cache/$name.tmp"
done

{ echo "$HDR"
  for c in "${CONFIGS[@]}"; do
    IFS='|' read -r name _ <<< "$c"
    if want "$name" && [ -s "$MDIR/mash-$name.tsv" ]; then tail -n+2 "$MDIR/mash-$name.tsv"; fi
  done
} > "$MDIR/all_mash.tsv"
echo "done -> $MDIR"
