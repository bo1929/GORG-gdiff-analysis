#!/usr/bin/env bash
set -euo pipefail

GENOME_DIR="${1:?Usage: $0 <genome_dir> [output_root]}"
OUTPUT_ROOT="${2:-./mash_output}"

GENOME_DIR="$(cd "$GENOME_DIR" && pwd)"
OUTPUT_ROOT="$(mkdir -p "$OUTPUT_ROOT" && cd "$OUTPUT_ROOT" && pwd)"

SKETCH_DIR="$OUTPUT_ROOT/sketches"
DIST_DIR="$OUTPUT_ROOT/distances"
THREADS="${MASH_THREADS:-8}"

mkdir -p "$SKETCH_DIR" "$DIST_DIR"

CONFIGS=(
  "cfg01_default|21|1000"
  "cfg02_sketch5k|21|5000"
  "cfg03_sketch10k|21|10000"
  "cfg04_sketch50k|21|50000"
  "cfg05_k31_default|31|1000"
  "cfg06_k31_sketch10k|31|10000"
  "cfg07_k16_default|16|1000"
  "cfg08_k16_sketch10k|16|10000"
  "cfg09_k25_sketch20k|25|20000"
  "cfg10_maxsens|16|50000"
)

GENOME_LIST="$OUTPUT_ROOT/genome_list.txt"
find -L "$GENOME_DIR" -type f \
  \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.fq" -o -name "*.fastq" \
     -o -name "*.fasta.gz" -o -name "*.fa.gz" -o -name "*.fna.gz" \
     -o -name "*.fq.gz" -o -name "*.fastq.gz" \) \
  | sort > "$GENOME_LIST"

NUM_GENOMES=$(wc -l < "$GENOME_LIST" | tr -d ' ')
if [[ "$NUM_GENOMES" -eq 0 ]]; then
  echo "ERROR: No FASTA/FASTQ files found in $GENOME_DIR" >&2
  exit 1
fi
echo "Found $NUM_GENOMES genomes in $GENOME_DIR"
echo "Output: $OUTPUT_ROOT"
echo ""

HEADER=$'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tp_value\tshared_hashes\tani_pct'

SUMMARY="$OUTPUT_ROOT/mash_runs.csv"
echo "name,k,s,sketch_file,dist_file,wall_seconds" > "$SUMMARY"

ALL_TSV="$DIST_DIR/all_mash.tsv"
echo "$HEADER" > "$ALL_TSV"

for entry in "${CONFIGS[@]}"; do
  IFS='|' read -r cfg_name k_val s_val <<< "$entry"
  param_setup="k=${k_val},s=${s_val}"

  echo "--- $cfg_name  [$param_setup]"

  SKETCH_FILE="$SKETCH_DIR/${cfg_name}.msh"
  RAW_FILE="$DIST_DIR/${cfg_name}_raw.tmp"
  DIST_FILE="$DIST_DIR/${cfg_name}.tsv"

  START_TIME=$SECONDS

  echo "  Sketching $NUM_GENOMES genomes..."
  mash sketch \
    -p "$THREADS" \
    -k "$k_val" \
    -s "$s_val" \
    -o "$SKETCH_FILE" \
    -l "$GENOME_LIST"

  echo "  Computing pairwise distances..."
  mash triangle \
    -p "$THREADS" \
    -k "$k_val" \
    -s "$s_val" \
    -E \
    -l "$GENOME_LIST" \
    > "$RAW_FILE"

  echo "  Formatting TSV..."
  {
    echo "$HEADER"
    awk -v method="mash" -v setup="$param_setup" \
      'BEGIN { FS="\t"; OFS="\t" }
       NF >= 5 {
         gsub(/^.*\/|\.(fasta|fa|fna|fq|fastq)(\.gz)?$/, "", $1);
         gsub(/^.*\/|\.(fasta|fa|fna|fq|fastq)(\.gz)?$/, "", $2);
         ani = (1 - $3) * 100;
         print method, setup, $1, $2, $3, $4, $5, ani
       }' "$RAW_FILE"
  } > "$DIST_FILE"

  rm -f "$RAW_FILE"

  ELAPSED=$((SECONDS - START_TIME))
  NUM_PAIRS=$(tail -n+2 "$DIST_FILE" | wc -l | tr -d ' ')

  tail -n+2 "$DIST_FILE" >> "$ALL_TSV"
  echo "$cfg_name,$k_val,$s_val,$SKETCH_FILE,$DIST_FILE,$ELAPSED" >> "$SUMMARY"
  echo "  Done in ${ELAPSED}s  ->  $DIST_FILE  ($NUM_PAIRS pairs)"
  echo ""
done

echo "---"
echo "All Mash runs complete."
echo "Per-config: $DIST_DIR/<name>.tsv"
echo "Combined:   $ALL_TSV"
echo "Summary:    $SUMMARY"
