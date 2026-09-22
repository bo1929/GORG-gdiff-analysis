#!/usr/bin/env bash
# skani ANI on every variant vs its baseline genome.
# Output saved in the same format as
#   results/ani-comparison/skani-estimates-selected/distances/skani-<cfg>.tsv
#   : method  param_setup  genome_a  genome_b  ani_pct  af_ref_pct  af_query_pct
#
# Config via env:
#   GENOMES      genome dir              (default: genomes)
#   SUFFIX       baseline suffix         (default: _contigs.fasta)
#   ONLY         comma-separated seed subset (default: all)
#   JOBS         parallel worker count   (default: 8)
#   ARGS         skani options           (default: "-c 30 -m 300 --slow")
#   CONFIG       short identifier label (default: derived from ARGS)
#   OUTDIR       results dir (default: results/)
#   NO_PROGRESS  set to 1 to disable the progress bar
#
# Usage:
#   ARGS="-c 200 -m 5000 --fast" CONFIG=fast ./run_skani.sh
# set -euo pipefail

GENOMES="${GENOMES:-genomes}"
SUFFIX="${SUFFIX:-_contigs.fasta}"
ONLY="${ONLY:-}"
JOBS="${JOBS:-8}"
# ARGS="${ARGS:---slow --min-af 0}"
ARGS="${ARGS:---robust --min-af 0 -c 70}"

c=$(echo "$ARGS" | grep -o '\-c [0-9]*' | awk '{print $2}')
m=$(echo "$ARGS" | grep -o '\-m [0-9]*' | awk '{print $2}')
PARAMS="c=${c},m=${m}"
CFGTMP=$(echo "$ARGS" | tr -d ' ' | sed 's/--/-/g')
CFGTMP="${CFGTMP#-}"
CONFIG="${CONFIG:-$CFGTMP}"
OUTDIR="${OUTDIR:-results//}"
NO_PROGRESS="${NO_PROGRESS:-0}"

DIST="$OUTDIR/distances/skani-${CONFIG}.tsv"
mkdir -p "$OUTDIR"/distances

command -v skani >/dev/null || { echo "missing: skani" >&2; exit 1; }

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
mkdir -p "$TMP/e"

printf 'method\tparam_setup\tgenome_a\tgenome_b\tani_pct\taf_ref_pct\taf_query_pct\n' > "$DIST"

# process_variant <idx> <seed> <ref_fasta> <var_fasta> <var_key>
process_variant() {
    local idx="$1" seed="$2" ref="$3" va="$4" vb="$5"
    local out ani af_ref af_qry
    # skani dist: Ref Query ANI Align_fraction_ref Align_fraction_query ...
    out="$(skani dist "$va" "$ref" $ARGS 2>/dev/null | tail -n1)"
    if [ -z "$out" ]; then
        printf 'skani\t%s\t%s\t%s\t.\t.\t.\n' "$PARAMS" "$vb" "$seed" > "$TMP/e/$idx"
        return
    fi
    read -r ani af_ref af_qry <<< "$(printf '%s' "$out" | awk '{print $3, $4, $5}')"
    printf 'skani\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$PARAMS" "$vb" "$seed" "$ani" "$af_ref" "$af_qry" > "$TMP/e/$idx"
}
export -f process_variant
export GENOMES SUFFIX ONLY ARGS PARAMS TMP

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

IDX=0
DONE=0
TOTAL=$(count_variants)

for seed in "$GENOMES"/*/; do
    seed=$(basename "$seed")
    [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
    if [ -n "$ONLY" ]; then case ",$ONLY," in *",$seed,"*) ;; *) continue ;; esac; fi
    # echo "$seed" >&2
    ref="$GENOMES/$seed/$seed$SUFFIX"

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
        process_variant "$IDX" "$seed" "$ref" "$vf" "$vb" &
        IDX=$((IDX+1)); running=$((running+1))
        if (( running % JOBS == 0 )); then
            wait; running=0
            DONE="$IDX"; progress "$DONE" "$TOTAL"
        fi
    done < "$var_list"
    wait
    DONE="$IDX"; progress "$DONE" "$TOTAL"
 done

for i in $(seq 0 $((IDX-1))); do cat "$TMP/e/$i" 2>/dev/null; done >> "$DIST"

echo "wrote $DIST"
