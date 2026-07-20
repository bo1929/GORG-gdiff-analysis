#!/usr/bin/env bash
set -euo pipefail

GENOME_DIR="${1:?Usage: $0 <genome_dir> [output_root]}"
OUTPUT_ROOT="${2:-./anib_output}"

GENOME_DIR="$(cd "$GENOME_DIR" && pwd)"
OUTPUT_ROOT="$(mkdir -p "$OUTPUT_ROOT" && cd "$OUTPUT_ROOT" && pwd)"

DIST_DIR="$OUTPUT_ROOT/distances"
DB_FILE="$OUTPUT_ROOT/anib.db"
EXPORT_DIR="$OUTPUT_ROOT/pyani_export"
mkdir -p "$DIST_DIR" "$EXPORT_DIR"

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

echo "Running pyani-plus ANIb..."
python3 -m pyani_plus.public_cli anib \
  "$GENOME_DIR" \
  -d "$DB_FILE" \
  --create-db \
  --name "gorg_anib_groundtruth"

echo ""
echo "Exporting results..."
python3 -m pyani_plus.public_cli export-run \
  -d "$DB_FILE" \
  -o "$EXPORT_DIR"

EXPORT_TSV=$(ls "$EXPORT_DIR"/anib_run_*.tsv 2>/dev/null | head -1)
if [[ -z "$EXPORT_TSV" ]]; then
  echo "ERROR: No export file found in $EXPORT_DIR" >&2
  exit 1
fi
echo "Export: $EXPORT_TSV"

HEADER=$'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct'
DIST_FILE="$DIST_DIR/anib.tsv"

echo "Formatting TSV..."
{
  echo "$HEADER"
  awk 'BEGIN { FS="\t"; OFS="\t" }
       NR > 1 && NF >= 5 && $3 != "None" {
         print "ANIb", "frag=1020,mode=ANIb", $1, $2, $3, $4, $5
       }' "$EXPORT_TSV"
} > "$DIST_FILE"

NUM_PAIRS=$(tail -n+2 "$DIST_FILE" | wc -l | tr -d ' ')
echo "---"
echo "ANIb ground truth complete."
echo "Output: $DIST_FILE  ($NUM_PAIRS pairs)"
