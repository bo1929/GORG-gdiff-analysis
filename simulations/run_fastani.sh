#!/usr/bin/env bash
# FastANI on every variant vs its baseline genome.
# Outputs (both like gdiff) saved under results/ani-comparison/fastani.../:
#   distances/fastani-<cfg>.tsv : method param_setup genome_a genome_b ani_pct af_ref_pct af_query_pct
#   samples/fastani-<cfg>.tsv   : config genome_a genome_b qid start end strand reference d lr_bg lr_ub
# The sample rows are per-fragment windows from fastANI --visualize
# (query coords + window ANI -> distance d = 1 - ANI/100).
#
# fastANI is asymmetric (issue #36): ANI differs by query/ref choice. We run
# BOTH directions per pair and emit a->b and b->a rows concatenated, like gdiff:
#   a->b (variant query, baseline ref) and b->a (baseline query, variant ref).
# One invocation per (variant x direction) with --visualize; JOBS parallelizes
# variants.
#
# Config via env:
#   GENOMES      genome dir              (default: genomes)
#   SUFFIX       baseline suffix         (default: _contigs.fasta)
#   ONLY         comma-separated seed subset (default: all)
#   JOBS         parallel worker count   (default: 8)
#   FRAG         fragment length         (default: 1000)
#   MINFRAC      --minFraction           (default: 0.1)
#   THREADS      fastANI threads per run (default: 4)
#   CONFIG       short identifier label  (default: derived from FRAG)
#   OUTDIR       results dir (default: results/)
#   NO_PROGRESS  set to 1 to disable the progress bar
#
# Usage:
#   FRAG=1000 ./run_fastani.sh
#   ONLY=AG-359-G18 ./run_fastani.sh
set -euo pipefail

GENOMES="${GENOMES:-genomes}"
SUFFIX="${SUFFIX:-_contigs.fasta}"
ONLY="${ONLY:-}"
JOBS="${JOBS:-8}"
FRAG="${FRAG:-1000}"
MINFRAC="${MINFRAC:-0.1}"
THREADS="${THREADS:-4}"
CONFIG="${CONFIG:-frag${FRAG}}"
OUTDIR="${OUTDIR:-results/}"
NO_PROGRESS="${NO_PROGRESS:-0}"

c=$(echo "$FRAG" | tr -d ' ')
PARAMS="frag=${c},minfrac=${MINFRAC}"
DIST="$OUTDIR/distances/fastani-${CONFIG}.tsv"
SAMP="$OUTDIR/samples/fastani-${CONFIG}.tsv"
mkdir -p "$OUTDIR"/distances "$OUTDIR"/samples

# Resolve fastANI (available in micromamba base env; not on default PATH).
FASTANI="${FASTANI:-}"
[ -z "$FASTANI" ] && command -v fastANI >/dev/null 2>&1 && FASTANI="$(command -v fastANI)"
[ -z "$FASTANI" ] && command -v micromamba >/dev/null 2>&1 \
    && FASTANI="$(micromamba run -n base which fastANI 2>/dev/null | tail -1)"
if [ -z "$FASTANI" ] || [ ! -x "$FASTANI" ]; then
    echo "missing: fastANI (e.g. micromamba install -n base -c bioconda fastani)" >&2
    exit 1
fi
export FASTANI

# Count total variants (respecting ONLY) for the progress bar.
count_variants() {
    local seed n=0 d vf
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

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/e" "$TMP/s"

printf 'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct\n' > "$DIST"
printf 'config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg\tlr_ub\n' > "$SAMP"

# variant_key <path>: emit the metadata-style variant key from a fasta path.
variant_key() {
    local pth="$1" d f
    d=$(dirname "$pth"); f=$(basename "$pth")
    d="$(basename "$d")"; d="${d#mutated_BLOSUM62_}"
    if [ "$f" = "mutated.fasta" ]; then
        printf '%s\n' "$d"
    else
        f="${f#mutated_}"; f="${f%.fasta}"
        printf '%s_%s\n' "$d" "$f"
    fi
}
export -f variant_key

# process_variant <idx> <seed> <var_fasta> <ref_fasta> <var_key>
# fastANI is asymmetric: ANI(gquery=variant, rref=baseline) differs from
# ANI(gquery=baseline, rref=variant) (issue #36). We run BOTH directions and
# emit a->b and b->a rows (with per-window samples), like gdiff.
process_variant() {
    local idx="$1" seed="$2" vf="$3" ref="$4" vb="$5"
    if [[ "$(pwd)/$vf" == "$(pwd)/$ref" ]]; then return 0; fi
    local out_ab="$TMP/fa_${idx}_ab" out_ba="$TMP/fa_${idx}_ba"

    # a->b: variant as query, baseline as reference
    if ! "$FASTANI" -q "$(pwd)/$vf" -r "$(pwd)/$ref" -o "$out_ab" \
        --fragLen "$FRAG" --minFraction "$MINFRAC" -t "$THREADS" --visualize >/dev/null 2>&1; then
        echo "fastani failed for $vb (a->b)" >&2; return 0; fi
    # b->a: baseline as query, variant as reference
    if ! "$FASTANI" -q "$(pwd)/$ref" -r "$(pwd)/$vf" -o "$out_ba" \
        --fragLen "$FRAG" --minFraction "$MINFRAC" -t "$THREADS" --visualize >/dev/null 2>&1; then
        echo "fastani failed for $vb (b->a)" >&2; return 0; fi

    # distance summaries
    ani_ab=$(awk -F'\t' '{print $3}' "$out_ab")
    ani_ba=$(awk -F'\t' '{print $3}' "$out_ba")
    {
        printf 'fastani\t%s\t%s\t%s\t%s\t.\t.\n' "$PARAMS" "$vb" "$seed" "$ani_ab"
        printf 'fastani\t%s\t%s\t%s\t%s\t.\t.\n' "$PARAMS" "$seed" "$vb" "$ani_ba"
    } > "$TMP/e/$idx"

    # per-window samples for both directions
    {
        # a->b : genome_a=variant, genome_b=baseline, reference=baseline
        awk -F'\t' -v c="$CONFIG" -v a="$vb" -v b="$seed" '
            { d = 1 - $3/100
              printf "%s\t%s\t%s\t.\t%d\t%d\t.\t%s\t%.6f\t.\t.\n",
                     c, a, b, $7+1, $8+1, b, d }' "$out_ab.visual"
        # b->a : genome_a=baseline, genome_b=variant, reference=variant
        awk -F'\t' -v c="$CONFIG" -v a="$vb" -v b="$seed" '
            { d = 1 - $3/100
              printf "%s\t%s\t%s\t.\t%d\t%d\t.\t%s\t%.6f\t.\t.\n",
                     c, b, a, $7+1, $8+1, a, d }' "$out_ba.visual"
    } > "$TMP/s/$idx"
    rm -f "$out_ab" "$out_ba" "$out_ab.visual" "$out_ba.visual"
}
export -f process_variant
export GENOMES SUFFIX ONLY FRAG MINFRAC THREADS PARAMS CONFIG TMP

IDX=0
DONE=0
TOTAL=$(count_variants)

for seed in "$GENOMES"/*/; do
    seed=$(basename "$seed")
    [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
    if [ -n "$ONLY" ]; then case ",$ONLY," in *",$seed,"*) ;; *) continue ;; esac; fi
    # echo "$seed" >&2
    ref="$GENOMES/$seed/$seed$SUFFIX"

    running=0
    for d in "$GENOMES/$seed"/mutated_BLOSUM62_*; do
        [ -d "$d" ] || continue; [ -f "$d/mutated.fasta" ] || continue
        core="$(basename "$d")"; core="${core#mutated_BLOSUM62_}"
        for vf in "$d/mutated.fasta" "$d"/mutated_p*_s*.fasta; do
            [ -f "$vf" ] || continue
            vb="$core"
            if [[ "$(basename "$vf")" != mutated.fasta ]]; then
                x="$(basename "$vf")"; x="${x#mutated_}"; x="${x%.fasta}"
                vb="${core}_${x}"
            fi
            process_variant "$IDX" "$seed" "$vf" "$ref" "$vb" &
            IDX=$((IDX+1)); running=$((running+1))
            if (( running % JOBS == 0 )); then wait; running=0; DONE="$IDX"; progress "$DONE" "$TOTAL"; fi
        done
    done
    wait
    DONE="$IDX"; progress "$DONE" "$TOTAL"
done

for i in $(seq 0 $((IDX-1))); do cat "$TMP/e/$i" 2>/dev/null; done >> "$DIST"
for i in $(seq 0 $((IDX-1))); do cat "$TMP/s/$i" 2>/dev/null; done >> "$SAMP"

echo "wrote $DIST and $SAMP"
