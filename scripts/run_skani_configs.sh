#!/usr/bin/env bash
set -euo pipefail

GENOME_DIR="${1:?Usage: $0 <genome_dir> [output_root]}"
OUTPUT_ROOT="${2:-./skani_output}"

GENOME_DIR="$(cd "$GENOME_DIR" && pwd)"
OUTPUT_ROOT="$(mkdir -p "$OUTPUT_ROOT" && cd "$OUTPUT_ROOT" && pwd)"

DIST_DIR="$OUTPUT_ROOT/distances"
THREADS="${SKANI_THREADS:-8}"

mkdir -p "$DIST_DIR"

GENOME_LIST="$OUTPUT_ROOT/genome_list.txt"
find -L "$GENOME_DIR" -type f \
  \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \
     -o -name "*.fasta.gz" -o -name "*.fa.gz" -o -name "*.fna.gz" \) \
  | sort > "$GENOME_LIST"

NUM_GENOMES=$(wc -l < "$GENOME_LIST" | tr -d ' ')
if [[ "$NUM_GENOMES" -eq 0 ]]; then
  echo "ERROR: No FASTA files found in $GENOME_DIR" >&2
  exit 1
fi
echo "Found $NUM_GENOMES genomes in $GENOME_DIR"
echo "Output: $OUTPUT_ROOT"
echo ""

echo "First 3 genomes:"
head -3 "$GENOME_LIST" | while read -r f; do
  echo "  $f  ($(wc -c < "$f" 2>/dev/null || echo "MISSING") bytes)"
done
echo ""

FIRST_GENOME=$(head -1 "$GENOME_LIST")
echo "Testing skani can read first genome..."
skani sketch "$FIRST_GENOME" -o "$OUTPUT_ROOT/.skani_test" 2>&1 || {
  echo "ERROR: skani failed to sketch $FIRST_GENOME" >&2
  exit 1
}
rm -rf "$OUTPUT_ROOT/.skani_test"
echo "Pre-flight OK."
echo ""

CONFIGS=(
  "cfg01_default|c=125,m=1000|-c 125 -m 1000"
  "cfg02_fast|c=200,m=1000,fast|--fast"
  "cfg03_medium|c=70,m=1000,medium|--medium"
  "cfg04_slow|c=30,m=1000,slow|--slow"
  "cfg05_sensitive|c=30,m=200|-c 30 -m 200"
  "cfg06_robust|c=70,m=1000,robust|-c 70 --robust"
  "cfg07_median|c=70,m=1000,median|-c 70 --median"
  "cfg08_raw_mean|c=70,m=1000,no_learned_ani|-c 70 --no-learned-ani"
  "cfg09_low_af|c=70,m=1000,min_af=5|-c 70 --min-af 5"
  "cfg10_max_sensitivity|c=30,m=100,robust,min_af=5|-c 30 -m 100 --robust --min-af 5"
)

HEADER=$'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct'

SUMMARY="$OUTPUT_ROOT/skani_runs.csv"
echo "name,param_setup,skani_flags,dist_file,wall_seconds,num_pairs" > "$SUMMARY"

ALL_TSV="$DIST_DIR/all_skani.tsv"
echo "$HEADER" > "$ALL_TSV"

for entry in "${CONFIGS[@]}"; do
  IFS='|' read -r cfg_name param_setup skani_flags <<< "$entry"

  echo "--- $cfg_name  [$param_setup]"
  echo "  flags: $skani_flags"

  RAW_FILE="$DIST_DIR/${cfg_name}_raw.tmp"
  DIST_FILE="$DIST_DIR/${cfg_name}.tsv"

  START_TIME=$SECONDS

  skani triangle \
    -t "$THREADS" \
    -l "$GENOME_LIST" \
    -E \
    $skani_flags \
    -o "$RAW_FILE"

  {
    echo "$HEADER"
    awk -v method="skani" -v setup="$param_setup" \
      'BEGIN { FS="\t"; OFS="\t" }
       NR > 1 && NF >= 5 {
         gsub(/^.*\/|\.(fasta|fa|fna)(\.gz)?$/, "", $1);
         gsub(/^.*\/|\.(fasta|fa|fna)(\.gz)?$/, "", $2);
         print method, setup, $1, $2, $3, $4, $5
       }' "$RAW_FILE"
  } > "$DIST_FILE"

  rm -f "$RAW_FILE"

  ELAPSED=$((SECONDS - START_TIME))
  NUM_PAIRS=$(tail -n+2 "$DIST_FILE" | wc -l | tr -d ' ')

  tail -n+2 "$DIST_FILE" >> "$ALL_TSV"
  echo "$cfg_name,\"$param_setup\",\"$skani_flags\",$DIST_FILE,$ELAPSED,$NUM_PAIRS" >> "$SUMMARY"
  echo "  Done in ${ELAPSED}s  ->  $DIST_FILE  ($NUM_PAIRS pairs)"
  echo ""
done

echo "---"
echo "All skani runs complete."
echo "Per-config: $DIST_DIR/<name>.tsv"
echo "Combined:   $ALL_TSV"
echo "Summary:    $SUMMARY"
