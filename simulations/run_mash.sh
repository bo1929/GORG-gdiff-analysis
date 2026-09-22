#!/usr/bin/env bash
# Mash ANI/dist on every variant vs its baseline genome.
# Output saved in the same format as
#   results/ani-comparison/mash-estimates-selected/distances/mash-<cfg>.tsv
#   : method  param_setup  genome_a  genome_b  distance  p_value  shared_hashes  ani_pct
#
# Config via env:
#   GENOMES      genome dir             (default: genomes)
#   SUFFIX       baseline suffix        (default: _contigs.fasta)
#   ONLY         comma-separated seed subset (default: all)
#   JOBS         parallel worker count   (default: 16)
#   k            mash kmer size          (default: 19)
#   s            mash sketch size        (default: 10000)
#   CONFIG       short label (default: derived from k and s)
#   OUTDIR       results dir (default: results/)
#   NO_PROGRESS  set to 1 to disable the progress bar
#
# Usage:
#   k=31 s=10000 CONFIG=long-k ./run_mash.sh
set -euo pipefail

GENOMES="${GENOMES:-genomes}"
SUFFIX="${SUFFIX:-_contigs.fasta}"
ONLY="${ONLY:-}"
JOBS="${JOBS:-16}"
k="${k:-19}"
s="${s:-10000}"
CONFIG="${CONFIG:-k$k-s$s}"
OUTDIR="${OUTDIR:-results/}"
NO_PROGRESS="${NO_PROGRESS:-0}"

DIST="$OUTDIR/distances/mash-$CONFIG.tsv"
mkdir -p "$OUTDIR"/distances

command -v mash >/dev/null || { echo "missing: mash" >&2; exit 1; }

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

printf 'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tp_value\tshared_hashes\tani_pct\n' > "$DIST"

# process_variant <idx> <seed> <base_msh> <var_fasta> <var_key>
process_variant() {
    local idx="$1" seed="$2" bmsh="$3" va="$4" vb="$5"
    local out d p shared ani
    # mash dist: base sketch + variant fasta -> ref query dist pval shared
    out="$(mash dist "$bmsh" "$va" 2>/dev/null | tail -n1)"
    [ -n "$out" ] || return
    read -r d p shared <<< "$(printf '%s' "$out" | awk -F'\t' '{print $3, $4, $5}')"
    ani=$(awk -v x="$d" 'BEGIN{printf "%.6f",(1-x)*100}')
    printf 'mash\tk=%s,s=%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$k" "$s" "$vb" "$seed" "$d" "$p" "$shared" "$ani" > "$TMP/e/$idx"
}
export -f process_variant
export GENOMES SUFFIX ONLY k s TMP

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

CACHE="${CACHE:-.cache/mash-sim}"
mkdir -p "$CACHE"

IDX=0
DONE=0
TOTAL=$(count_variants)

for seed in "$GENOMES"/*/; do
    seed=$(basename "$seed")
    [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
    if [ -n "$ONLY" ]; then case ",$ONLY," in *",$seed,"*) ;; *) continue ;; esac; fi
    # echo "$seed" >&2

    msh="$CACHE/${seed}.msh"
    if [ ! -f "$msh" ]; then
        mash sketch -k "$k" -s "$s" -o "${msh%.msh}" "$GENOMES/$seed/$seed$SUFFIX" >/dev/null 2>&1
    fi

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
        process_variant "$IDX" "$seed" "$msh" "$vf" "$vb" &
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
