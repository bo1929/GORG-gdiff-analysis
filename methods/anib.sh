#!/usr/bin/env bash
# ANIb (pyani-plus) over genomes in pairs.tsv. All-vs-all on the subset, then
# filter to requested pairs. Slow (BLAST); prefer tiny pair lists for tests.
# usage: anib.sh <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
# env: THREADS unused (pyani local executor) FORCE=0 ONLY=default
#      CACHE_ROOT (default <repo>/.cache)
# out: distances/anib-<cfg>.tsv
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=_lib.sh
source "$HERE/_lib.sh"

GENOME_DIR="$(cd "${1:?usage: $0 <genome_dir> <pairs.tsv> [outdir] [suffix]}" && pwd)"
PAIRS_FILE="${2:?}"
OUT="${3:-./methods_out}"; SUFFIX="${4:-.fasta}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
use_cache anib
DIST_DIR="$OUT/distances"
mkdir -p "$DIST_DIR"
FORCE="${FORCE:-0}"; ONLY="${ONLY:-default}"
python3 -c "import pyani_plus" >/dev/null 2>&1 || { echo "missing: pyani_plus" >&2; exit 1; }
command -v blastn >/dev/null || { echo "missing: blastn" >&2; exit 1; }

CONFIGS=(
  "default|frag=1020|--fragsize 1020"
)
HDR=$'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct'
load_pairs

# Stage only genomes that appear in pairs (symlink farm for pyani).
stage_genomes() {
  local indir="$1"
  rm -rf "$indir"
  mkdir -p "$indir"
  while read -r id; do
    ln -s "$(fa "$id")" "$indir/${id}.fasta"
  done < "$CACHE/genomes.txt"
}

for c in "${CONFIGS[@]}"; do
  IFS='|' read -r name setup flags <<< "$c"
  want "$name" || continue
  tsv="$DIST_DIR/anib-$name.tsv"
  if [ "$FORCE" != 1 ] && [ -s "$tsv" ]; then echo "$name: skip"; continue; fi
  echo "$name [$setup]" >&2

  run_dir="$CACHE/$name"
  indir="$run_dir/indir"
  db="$run_dir/anib.db"
  export_dir="$run_dir/export"
  mkdir -p "$run_dir" "$export_dir"
  stage_genomes "$indir"
  ngen=$(wc -l < "$CACHE/genomes.txt" | tr -d ' ')
  echo "  $ngen genomes -> pyani-plus anib" >&2

  [ "$FORCE" = 1 ] && rm -f "$db"
  # shellcheck disable=SC2086
  python3 -m pyani_plus.public_cli anib "$indir" \
    -d "$db" --create-db --name "methods_anib_$name" \
    --log - $flags

  rm -f "$export_dir"/ANIb_run_*.tsv "$export_dir"/anib_run_*.tsv
  python3 -m pyani_plus.public_cli export-run \
    -d "$db" -o "$export_dir" --label stem --log -

  export_tsv=$(ls "$export_dir"/ANIb_run_*.tsv "$export_dir"/anib_run_*.tsv 2>/dev/null | head -1 || true)
  [ -n "$export_tsv" ] || { echo "export failed: no ANIb_run_*.tsv" >&2; exit 1; }

  # Long-form: #Query Subject Identity Query-Cov Subject-Cov ... (fractions).
  # Average both directed comparisons; emit one row per pairs.tsv entry as percent.
  { echo "$HDR"
    awk -v setup="$setup" 'BEGIN{FS=OFS="\t"}
      FNR==NR {
        if ($1 != $2) { order[++np]=$1 FS $2; want[$1 FS $2]=1 }
        next
      }
      FNR==1 { next }  # header (#Query ...)
      NF>=5 && $3+0==$3 {
        a=$1; b=$2; ani=$3+0; cq=$4+0; cs=($5=="NA"||$5==""? "": $5+0)
        key=a FS b; rkey=b FS a
        if (!(key in want) && !(rkey in want)) next
        # store under canonical pairs.tsv orientation when possible
        if (key in want) {
          sum_ani[key]+=ani; n_ani[key]++
          if (cq!="") { sum_aq[key]+=cq; n_aq[key]++ }
          if (cs!="") { sum_ar[key]+=cs; n_ar[key]++ }
        } else {
          # reverse: query_cov of (b,a) is af for b as query -> af_ref when orienting as (a,b)
          sum_ani[rkey]+=ani; n_ani[rkey]++
          if (cq!="") { sum_ar[rkey]+=cq; n_ar[rkey]++ }
          if (cs!="") { sum_aq[rkey]+=cs; n_aq[rkey]++ }
        }
      }
      END {
        for (i=1;i<=np;i++) {
          k=order[i]; if (!(k in n_ani)) continue
          split(k, ab, FS); a=ab[1]; b=ab[2]
          ani_pct = 100*sum_ani[k]/n_ani[k]
          aq = (n_aq[k] ? 100*sum_aq[k]/n_aq[k] : "NA")
          ar = (n_ar[k] ? 100*sum_ar[k]/n_ar[k] : "NA")
          printf "anib\t%s\t%s\t%s\t%.6f\t%s\t%s\n", setup, a, b, ani_pct, ar, aq
        }
      }' "$CACHE/pairs.tsv" "$export_tsv"
  } > "$tsv"
  echo "  $(tail -n+2 "$tsv" | wc -l | tr -d ' ') pairs -> $tsv" >&2
done

emit_all_distances anib "$HDR" "${CONFIGS[@]}"
echo "done -> $DIST_DIR" >&2
