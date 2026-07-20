#!/usr/bin/env bash
# ============================================================================
# run_mash_configs.sh — Test 10 Mash configurations for pairwise distance
# ============================================================================
# Usage:
#   ./run_mash_configs.sh <genome_dir> [output_root]
#
# Each config varies sketch size (-s) and/or k-mer size (-k).
# Every output row carries the method name and parameter setup string so that
# all configs can be concatenated and compared side by side.
#
# Output per config:
#   <output_root>/
#     sketches/<cfg>.msh
#     distances/<cfg>.tsv    ← uniform-column TSV
#     all_mash.tsv            ← concatenation of all config TSVs
# ============================================================================

set -euo pipefail

GENOME_DIR="${1:?Usage: $0 <genome_dir> [output_root]}"
OUTPUT_ROOT="${2:-./mash_output}"

GENOME_DIR="$(cd "$GENOME_DIR" && pwd)"
OUTPUT_ROOT="$(mkdir -p "$OUTPUT_ROOT" && cd "$OUTPUT_ROOT" && pwd)"

SKETCH_DIR="$OUTPUT_ROOT/sketches"
DIST_DIR="$OUTPUT_ROOT/distances"
THREADS="${MASH_THREADS:-8}"

mkdir -p "$SKETCH_DIR" "$DIST_DIR"

# ------------------------------------------------------------------
# Configuration definitions — indexed array of  "name|k|s"
# Key trade-offs:
#   - Larger k  → more specific, less sensitive at distance
#   - Smaller k → more sensitive for distant relatives, more noise
#   - Larger s  → better resolution / lower variance, more memory+CPU
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# Build genome list
# ------------------------------------------------------------------
GENOME_LIST="$OUTPUT_ROOT/genome_list.txt"
find "$GENOME_DIR" -type f \
  \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.fq" -o -name "*.fastq" -o -name "*.gbk" -o -name "*.gbff" \) \
  | sort > "$GENOME_LIST"

NUM_GENOMES=$(wc -l < "$GENOME_LIST" | tr -d ' ')
echo "Found $NUM_GENOMES genome files in $GENOME_DIR"
echo "Output root: $OUTPUT_ROOT"
echo ""

# Column header for the uniform TSV
HEADER="method\tparam_setup\tgenome_a\tgenome_b\tdistance\tp_value\tshared_hashes\tani_pct"

# ------------------------------------------------------------------
# Summary log (metadata about each run)
# ------------------------------------------------------------------
SUMMARY="$OUTPUT_ROOT/mash_runs.csv"
echo "name,k,s,sketch_file,dist_file,wall_seconds" > "$SUMMARY"

# ------------------------------------------------------------------
# Run each configuration
# ------------------------------------------------------------------
FIRST_TSV=true
ALL_TSV="$DIST_DIR/all_mash.tsv"
printf '%s\n' "$HEADER" > "$ALL_TSV"

for entry in "${CONFIGS[@]}"; do
  IFS='|' read -r cfg_name k_val s_val <<< "$entry"

  # Build the parameter setup identifier string (human-readable, one-line)
  param_setup="k=${k_val},s=${s_val}"

  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "▶ $cfg_name   [$param_setup]"
  echo ""

  SKETCH_FILE="$SKETCH_DIR/${cfg_name}.msh"
  RAW_FILE="$DIST_DIR/${cfg_name}_raw.tmp"
  DIST_FILE="$DIST_DIR/${cfg_name}.tsv"

  START_TIME=$SECONDS

  # --- Sketch ---
  echo "  Sketching $NUM_GENOMES genomes ..."
  mash sketch \
    -p "$THREADS" \
    -k "$k_val" \
    -s "$s_val" \
    -o "$SKETCH_FILE" \
    -l "$GENOME_LIST"

  # --- Pairwise triangle in edge-list mode ---
  # mash triangle reads the list and sketches on the fly unless given .msh files.
  # Using sketch files is faster but requires -k to match; we use raw FASTA for
  # simplicity and correctness.
  echo "  Computing pairwise distances (triangle -E) ..."
  mash triangle \
    -p "$THREADS" \
    -k "$k_val" \
    -s "$s_val" \
    -E \
    -l "$GENOME_LIST" \
    > "$RAW_FILE"

  # --- Transform into uniform TSV ---
  # mash triangle -E outputs tab-separated:
  #   seq1  seq2  distance  p-value  shared-hashes
  # We prepend method + param_setup, and append ANI = (1 - distance) * 100
  echo "  Formatting TSV ..."
  {
    printf '%s\n' "$HEADER"
    awk -v method="mash" -v setup="$param_setup" \
      'BEGIN { FS="\t"; OFS="\t" }
       NF >= 5 {
         gsub(/^.*\/|\.(fasta|fa|fna|fq|fastq|gbk|gbff)$/, "", $1);
         gsub(/^.*\/|\.(fasta|fa|fna|fq|fastq|gbk|gbff)$/, "", $2);
         ani = (1 - $3) * 100;
         print method, setup, $1, $2, $3, $4, $5, ani
       }' "$RAW_FILE"
  } > "$DIST_FILE"

  rm -f "$RAW_FILE"

  ELAPSED=$((SECONDS - START_TIME))
  NUM_PAIRS=$(tail -n+2 "$DIST_FILE" | wc -l | tr -d ' ')

  # Append to the all-in-one TSV (skip header)
  tail -n+2 "$DIST_FILE" >> "$ALL_TSV"

  echo "$cfg_name,$k_val,$s_val,$SKETCH_FILE,$DIST_FILE,$ELAPSED" >> "$SUMMARY"
  echo "  ✓ Done in ${ELAPSED}s  →  $DIST_FILE  ($NUM_PAIRS pairs)"
  echo ""
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "All Mash configurations complete."
echo ""
echo "Per-config TSV:     $DIST_DIR/<name>.tsv"
echo "Concatenated TSV:   $ALL_TSV"
echo "Run summary:        $SUMMARY"
echo ""
echo "Columns:  method, param_setup, genome_a, genome_b, distance,"
echo "          p_value, shared_hashes, ani_pct"
echo ""
echo "  distance = Mash distance  (0 = identical, ~0.05 = 95%% ANI, ...)"
echo "  ani_pct  = (1 - distance) × 100"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
