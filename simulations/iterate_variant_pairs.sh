#!/usr/bin/env bash
# Iterate over each baseline genome and its variant -> baseline pairs.
#
# For every seed under genomes/, the baseline reference is the original
# <seed>_contigs.fasta. Each variant is paired against that baseline:
#   - complete variants:  every gnd/aad level dir (mutated.fasta)
#   - pruned variants:    every nested mutated_p<miss>_s<s>.fasta
#
# Prints one pair per line:  <variant_key>  <baseline_ref>
#
# Config via environment variables:
#   GENOMES   directory containing seed genome folders (default: genomes)
#   ONLY      comma-separated subset of seeds (default: all)
#   FORMAT    key | path   (default: key)
# Example:
#   ONLY=AG-359-G18 ./iterate_variant_pairs.sh
#   FORMAT=path  ./iterate_variant_pairs.sh
set -euo pipefail

GENOMES="${GENOMES:-genomes}"
SUFFIX="_contigs.fasta"        # baseline file suffix (minus seed prefix)
ONLY="${ONLY:-}"
FORMAT="${FORMAT:-key}"        # key | path

variants_of_seed() {
    local seed="$1" d p base
    for d in "$GENOMES/$seed"/mutated_BLOSUM62_*; do
        [ -d "$d" ] || continue
        [ -f "$d/mutated.fasta" ] || continue
        base=$(basename "$d")
        emit "$seed" "${base#mutated_BLOSUM62_}" "$d/mutated.fasta"
        for p in "$d"/mutated_p*_s*.fasta; do
            [ -f "$p" ] || continue
            emit "$seed" "${base#mutated_BLOSUM62_}_${p##*mutated_}" "$p"
        done
    done
}

emit() {
    local seed="$1" vkey="$2" vpath="$3"
    if [ "$FORMAT" = "path" ]; then
        printf '%s\t%s\n' "$vpath" "$GENOMES/$seed/$seed$SUFFIX"
    else
        printf '%s\t%s\n' "$vkey" "$seed$SUFFIX"
    fi
}

main() {
    local seed
    for seed in "$GENOMES"/*/; do
        seed=$(basename "$seed")
        [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
        if [ -n "$ONLY" ]; then
            case ",$ONLY," in
                *",$seed,"*) ;;
                *) continue ;;
            esac
        fi
        variants_of_seed "$seed"
    done
}

main
