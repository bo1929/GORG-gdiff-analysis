#!/usr/bin/env bash
# gdiff ANI+dist (and per-window samples) for every variant vs its baseline.
#
# Targets the bundled `gdiff` (bin/gdiff dispatches to gdiff-osx/-x86), v0.2.0.
#   sketch --input-list <name<TAB>path> ... -o base.gdsk   # baseline, once per seed
#   sketch --input-list <name<TAB>path> ... -o var.gdsk    # variant
#   dist   var.gdsk base.gdsk [--output-samples]           # both directions
# Both sketches must share k/w/h/-l/--sample-size (gdiff rejects mismatched
# pairs), so the variant sketch is built with the same SKETCH_ARGS.
#
# v0.2.0 outputs (header row on both):
#   summary : genome_a genome_b d d_median d_mean d_upper d_highest d_ab d_ba
#             n_ab n_ba n_ub n_na n_filtered       (genome_a/genome_b canonical;
#             `d` is the reconciled distance, lr-filtered via --lr-th/--min-portion)
#   samples : config genome_a genome_b seq start end strand direction d lr_ub
#             (genome_a/genome_b canonical, `direction` ab|ba = which was query)
# We normalise samples to the repo schema with genome_a = query:
#   config genome_a genome_b qid start end strand reference d lr_bg lr_ub
#
# Outputs written in the same format as
#   results/ani-comparison/gdiff-estimates-selected/{distances,samples}/
#   distances/gdiff-<cfg>.tsv  : method  param_setup  genome_a  genome_b  distance  ani_pct
#   samples/gdiff-<cfg>.tsv    : config  genome_a  genome_b  qid ... (see below)
#
# Config via env:
#   GDIFF        binary                 (default: ../bin/gdiff, arch dispatcher)
#   GENOMES      genome dir             (default: genomes)
#   SUFFIX       baseline suffix        (default: _contigs.fasta)
#   ONLY         comma-separated seed subset (default: all)
#   JOBS         parallel worker count   (default: 16)
#   SKETCH_ARGS  e.g. "-k 23 -w 23 --frac 0.5 -l 500 --sample-size 1000"
#                --frac subsampling ratio in [0,1]; default 1.0 (not given).
#                It is reported as a 0-100 percent in CONFIG/PARAM, e.g. frac100.
#   DIST_ARGS    e.g. "--hdist-th 4 [--lr-th 3.841] [--min-portion 0.66]"
#   CONFIG       short label (default: derived from SKETCH_ARGS)
#   OUTDIR       results dir (default: results/)   -> results/distances|samples/
#   NO_PROGRESS  set to 1 to disable the progress bar
#
# Usage:
#   SKETCH_ARGS="-k 23 -w 23 --frac 0.5 -l 500 --sample-size 1000" JOBS=8 ./run_gdiff.sh
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"

GDIFF="${GDIFF:-$REPO_ROOT/bin/gdiff}"
GENOMES="${GENOMES:-genomes}"
SUFFIX="${SUFFIX:-_contigs.fasta}"
ONLY="${ONLY:-}"
JOBS="${JOBS:-16}"
SKETCH_ARGS="${SKETCH_ARGS:--k 23 -w 23 -l 333 --frac 0.50 -h 9 --sample-size 1000}"
DIST_ARGS="${DIST_ARGS:---hdist-th 3}"
OUTDIR="${OUTDIR:-results/}"
NO_PROGRESS="${NO_PROGRESS:-0}"

# In-place progress on stderr; newline when done==total.
progress() {
    [ "$NO_PROGRESS" = 1 ] && return 0
    local done="$1" total="$2" width=40 pct filled empty i bar
    [ "$total" -gt 0 ] || return 0
    pct=$(( done * 100 / total ))
    filled=$(( done * width / total )); [ "$filled" -gt "$width" ] && filled=$width
    empty=$(( width - filled )); bar=""
    i=0; while [ "$i" -lt "$filled" ]; do bar="${bar}#"; i=$((i+1)); done
    i=0; while [ "$i" -lt "$empty" ]; do bar="${bar} "; i=$((i+1)); done
    printf '\rprogress [%s] %d/%d (%d%%)' "$bar" "$done" "$total" "$pct" >&2
    [ "$done" -ge "$total" ] && printf '\n' >&2
    return 0
}

k=$(echo "$SKETCH_ARGS" | grep -o '\-k [0-9]*' | awk '{print $2}')
w=$(echo "$SKETCH_ARGS" | grep -o '\-w [0-9]*' | awk '{print $2}')
l=$(echo "$SKETCH_ARGS" | grep -o '\-l [0-9]*' | awk '{print $2}')
h=$(echo "$SKETCH_ARGS" | grep -o '\-h [0-9]*' | awk '{print $2}')
n=$(echo "$SKETCH_ARGS" | grep -o '\-sample-size [0-9]*' | awk '{print $2}')
frac=$(echo "$SKETCH_ARGS" | grep -o '\-\-frac [0-9.]*' | awk '{print $2}') || true
frac="${frac:-1.0}"
delta=$(echo "$DIST_ARGS" | grep -o '\-hdist-th [0-9]*' | awk '{print $2}')
fracpct=$(awk -v f="$frac" 'BEGIN{ printf "%d", (f*100) + 0.5 }')
CONFIG="${CONFIG:-frac${fracpct}-k${k}-w${w}-h${h}-l${l}-n${n}-delta${delta}}"
PARAM="k=${k},w=${w},frac=${frac},-l=${l},n=${n},delta=${delta}"

DST="$OUTDIR/distances/gdiff-$CONFIG.tsv"
SAMP="$OUTDIR/samples/gdiff-$CONFIG.tsv"
mkdir -p "$OUTDIR"/distances "$OUTDIR"/samples

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/s"

# Only samples are saved as a per-window table; the per-pair summary is folded
# into DST (one row per variant, symmetric distance).
printf 'config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg\tlr_ub\n' > "$SAMP"
printf 'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct\n' > "$DST"

# process_variant <idx> <seed> <base_sketch> <var_fasta> <var_key>
# Runs the symmetric `gdiff dist (var, base)`. v0.2.0 --output-samples emits one
# header + window rows: config genome_a genome_b seq start end strand direction
# d lr_ub, with genome_a/genome_b canonically ordered and `direction` = query.
# We normalise to the repo schema with genome_a = query (genome_a = the variant
# in the ba rows): config genome_a genome_b qid start end strand reference d
# lr_bg lr_ub (lr_bg is not reported by v0.2.0 -> ".").
# The per-pair summary reports the reconciled distance in column `d`.
process_variant() {
    local idx="$1" seed="$2" bsk="$3" va="$4" vb="$5"
    local smp="$TMP/s/$idx" sum="$TMP/d/$idx" raw="$TMP/raw_$idx"

    printf '%s\t%s\n' "$vb" "$va" > "$TMP/q_$idx.list"
    "$GDIFF" sketch --input-list "$TMP/q_$idx.list" $SKETCH_ARGS -o "$TMP/q_$idx.gdsk" >/dev/null 2>&1
    rm -f "$TMP/q_$idx.list"

    # shellcheck disable=SC2086
    "$GDIFF" dist "$TMP/q_$idx.gdsk" "$bsk" $DIST_ARGS --output-samples > "$raw" 2>/dev/null
    awk -F'\t' -v c="$CONFIG" '
      NR > 1 && $1 != "config" && NF >= 10 {
        if ($8 == "ba") { ga = $3; gb = $2; ref = $2 }   # query = genome_b
        else            { ga = $2; gb = $3; ref = $3 }   # query = genome_a
        print c"\t"ga"\t"gb"\t"$4"\t"$5"\t"$6"\t"$7"\t"ref"\t"$9"\t.\t"$10
      }' "$raw" > "$smp"

    # shellcheck disable=SC2086
    "$GDIFF" dist "$TMP/q_$idx.gdsk" "$bsk" $DIST_ARGS > "$raw" 2>/dev/null
    awk -F'\t' -v p="$PARAM" -v a="$vb" -v b="$seed" '
      /^#/ { next }
      $1 == "genome_a" { for (i = 1; i <= NF; i++) if ($i == "d") { col = i; break } ; next }
      col && NF >= col { printf "gdiff\t%s\t%s\t%s\t%s\t%.6g\n", p, a, b, $col, (1 - $col) * 100; exit }
    ' "$raw" > "$sum"
    rm -f "$TMP/q_$idx.gdsk" "$raw"
}
export -f process_variant
export GDIFF SKETCH_ARGS DIST_ARGS PARAM CONFIG SUFFIX TMP GENOMES

# Count total variants (respecting ONLY) for the progress bar.
count_variants() {
    local seed n=0 d core vf f
    for seed in "$GENOMES"/*/; do
        seed=$(basename "$seed")
        [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
        if [ -n "$ONLY" ]; then case ",$ONLY," in *",$seed,"*) ;; *) continue ;; esac; fi
        for d in "$GENOMES/$seed"/mutated_BLOSUM62_*; do
            [ -d "$d" ] || continue; [ -f "$d/mutated.fasta" ] || continue
            n=$((n+1))
            for vf in "$d"/mutated_p*_s*.fasta; do [ -f "$vf" ] && n=$((n+1)); done
        done
    done
    echo "$n"
}

IDX=0
DONE=0
TOTAL=$(count_variants)
mkdir -p "$TMP/d"

for seed in "$GENOMES"/*/; do
    seed=$(basename "$seed")
    [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
    if [ -n "$ONLY" ]; then case ",$ONLY," in *",$seed,"*) ;; *) continue ;; esac; fi
    # echo "$seed" >&2
    # Name the baseline record after the seed so sample `reference` fields (and
    # the distance summary) carry the seed id instead of the FASTA basename.
    # shellcheck disable=SC2086
    printf '%s\t%s\n' "$seed" "$GENOMES/$seed/$seed$SUFFIX" > "$TMP/base.list"
    "$GDIFF" sketch --input-list "$TMP/base.list" $SKETCH_ARGS -o "$TMP/base.gdsk" >/dev/null 2>&1
    rm -f "$TMP/base.list"

    var_list="$TMP/vars_$seed"
    : > "$var_list"
    for d in "$GENOMES/$seed"/mutated_BLOSUM62_*; do
        [ -d "$d" ] || continue; [ -f "$d/mutated.fasta" ] || continue
        core="$(basename "$d")"; core="${core#mutated_BLOSUM62_}"
        for vf in "$d/mutated.fasta" "$d"/mutated_p*_s*.fasta; do
            [ -f "$vf" ] || continue
            vb="$core"
            if [[ "$(basename "$vf")" != mutated.fasta ]]; then
                f="$(basename "$vf")"; f="${f#mutated_}"; f="${f%.fasta}"
                vb="${core}_${f}"
            fi
            printf '%s\t%s\n' "$vf" "$vb" >> "$var_list"
        done
    done

    running=0
    while IFS= read -r line; do
        vf="${line%%	*}"; vb="${line##*	}"
        process_variant "$IDX" "$seed" "$TMP/base.gdsk" "$vf" "$vb" &
        IDX=$((IDX+1)); running=$((running+1))
        if (( running % JOBS == 0 )); then
            wait; running=0
            DONE="$IDX"; progress "$DONE" "$TOTAL"
        fi
    done < "$var_list"
    wait
    DONE="$IDX"; progress "$DONE" "$TOTAL"
 done

for i in $(seq 0 $((IDX-1))); do cat "$TMP/s/$i" 2>/dev/null; done >> "$SAMP"
for i in $(seq 0 $((IDX-1))); do cat "$TMP/d/$i" 2>/dev/null; done >> "$DST"

echo "wrote $SAMP"
echo "wrote $DST"
