#!/usr/bin/env bash
# ============================================================================
# run_skani_configs.sh — Test 10 skani configurations for pairwise ANI
# ============================================================================
# Usage:
#   ./run_skani_configs.sh <genome_dir> [output_root]
#
# Each config varies compression (-c), marker compression (-m), presets,
# ANI estimator, and AF thresholds. Every output row carries the method
# name and parameter setup string so all configs can be concatenated and
# compared side by side.
#
# Output per config:
#   <output_root>/
#     distances/<cfg>.tsv    ← uniform-column TSV
#     all_skani.tsv           ← concatenation of all config TSVs
# ============================================================================

set -euo pipefail

GENOME_DIR="${1:?Usage: $0 <genome_dir> [output_root]}"
OUTPUT_ROOT="${2:-./skani_output}"

GENOME_DIR="$(cd "$GENOME_DIR" && pwd)"
OUTPUT_ROOT="$(mkdir -p "$OUTPUT_ROOT" && cd "$OUTPUT_ROOT" && pwd)"

DIST_DIR="$OUTPUT_ROOT/distances"
THREADS="${SKANI_THREADS:-8}"

mkdir -p "$DIST_DIR"

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

# ------------------------------------------------------------------
# Configuration definitions — indexed array of  "name|param_setup_string|skani_flags"
#
# The param_setup_string is the human-readable shorthand that goes into
# the TSV column. The skani_flags are the actual CLI arguments.
#
# Key parameters:
#   -c <C>       Compression factor. Lower = denser k-mers, more accurate
#                AF/ANI esp. for fragmented genomes. [default: 125]
#   -m <MARKER>  Marker k-mer compression for filtering. Lower = more
#                sensitive filtering. [default: 1000]
#   --min-af N   Only report pairs with > N% aligned fraction. [default: 15]
#   --robust     Trim 10/90% quantiles before estimating mean.
#   --median     Use median instead of mean ANI.
#   --no-learned-ani  Disable regression model; use raw mean ANI.
#   --fast       = -c 200
#   --medium     = -c 70
#   --slow       = -c 30
# ------------------------------------------------------------------
CONFIGS=(
  # 1. Default
  "cfg01_default|c=125,m=1000|-c 125 -m 1000"

  # 2. Fast preset
  "cfg02_fast|c=200,m=1000,fast|--fast"

  # 3. Medium preset
  "cfg03_medium|c=70,m=1000,medium|--medium"

  # 4. Slow preset
  "cfg04_slow|c=30,m=1000,slow|--slow"

  # 5. Low compression + sensitive markers
  "cfg05_sensitive|c=30,m=200|-c 30 -m 200"

  # 6. Robust ANI (trimmed mean) on medium compression
  "cfg06_robust|c=70,m=1000,robust|-c 70 --robust"

  # 7. Median ANI on medium compression
  "cfg07_median|c=70,m=1000,median|-c 70 --median"

  # 8. Raw mean ANI (no learned regression correction)
  "cfg08_raw_mean|c=70,m=1000,no_learned_ani|-c 70 --no-learned-ani"

  # 9. Low AF threshold (keep more distant pairs)
  "cfg09_low_af|c=70,m=1000,min_af=5|-c 70 --min-af 5"

  # 10. Maximum sensitivity — most permissive filtering + robust estimation
  "cfg10_max_sensitivity|c=30,m=100,robust,min_af=5|-c 30 -m 100 --robust --min-af 5"
)

# ------------------------------------------------------------------
# Column header for the uniform TSV
# ------------------------------------------------------------------
# skani triangle -E outputs tab-separated:
#   Ref  Query  ANI  AlignFracRef  AlignFracQuery  [RefNameLen  QueryNameLen]
# We keep the core columns and add method / param_setup.
HEADER="method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct"

# ------------------------------------------------------------------
# Summary log
# ------------------------------------------------------------------
SUMMARY="$OUTPUT_ROOT/skani_runs.csv"
echo "name,param_setup,skani_flags,dist_file,wall_seconds,num_pairs" > "$SUMMARY"

# ------------------------------------------------------------------
# All-in-one concatenated output
# ------------------------------------------------------------------
ALL_TSV="$DIST_DIR/all_skani.tsv"
printf '%s\n' "$HEADER" > "$ALL_TSV"

# ------------------------------------------------------------------
# Run each configuration
# ------------------------------------------------------------------
for entry in "${CONFIGS[@]}"; do
  IFS='|' read -r cfg_name param_setup skani_flags <<< "$entry"

  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "▶ $cfg_name   [$param_setup]"
  echo "  skani flags: $skani_flags"
  echo ""

  RAW_FILE="$DIST_DIR/${cfg_name}_raw.tmp"
  DIST_FILE="$DIST_DIR/${cfg_name}.tsv"

  START_TIME=$SECONDS

  # --- Run skani triangle (all-to-all, sparse/edge-list output) ---
  # skani does its own internal sketching; no separate sketch step needed.
  skani triangle \
    -t "$THREADS" \
    -l "$GENOME_LIST" \
    -E \
    $skani_flags \
    -o "$RAW_FILE"

  # --- Transform into uniform TSV ---
  echo "  Formatting TSV ..."
  {
    printf '%s\n' "$HEADER"
    awk -v method="skani" -v setup="$param_setup" \
      'BEGIN { FS="\t"; OFS="\t" }
       NF >= 5 {
         # Strip paths and extensions from genome names
         gsub(/^.*\/|\.(fasta|fa|fna|fq|fastq|gbk|gbff)$/, "", $1);
         gsub(/^.*\/|\.(fasta|fa|fna|fq|fastq|gbk|gbff)$/, "", $2);
         print method, setup, $1, $2, $3, $4, $5
       }' "$RAW_FILE"
  } > "$DIST_FILE"

  rm -f "$RAW_FILE"

  ELAPSED=$((SECONDS - START_TIME))
  NUM_PAIRS=$(tail -n+2 "$DIST_FILE" | wc -l | tr -d ' ')

  # Append to the all-in-one TSV (skip header)
  tail -n+2 "$DIST_FILE" >> "$ALL_TSV"

  echo "$cfg_name,\"$param_setup\",\"$skani_flags\",$DIST_FILE,$ELAPSED,$NUM_PAIRS" >> "$SUMMARY"
  echo "  ✓ Done in ${ELAPSED}s  →  $DIST_FILE  ($NUM_PAIRS pairs)"
  echo ""
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "All skani configurations complete."
echo ""
echo "Per-config TSV:     $DIST_DIR/<name>.tsv"
echo "Concatenated TSV:   $ALL_TSV"
echo "Run summary:        $SUMMARY"
echo ""
echo "Columns:  method, param_setup, genome_a, genome_b, ani_pct,"
echo "          af_ref_pct, af_query_pct"
echo ""
echo "  ani_pct      = estimated ANI (0–100)"
echo "  af_ref_pct   = fraction of reference genome aligned"
echo "  af_query_pct = fraction of query genome aligned"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
