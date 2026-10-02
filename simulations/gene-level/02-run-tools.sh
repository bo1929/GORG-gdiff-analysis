#!/usr/bin/env bash
# --- help ---
# STAGE 02 - run wfmash and gdiff sketch/roll over pairs, timed.
#
#   ./02-run-tools.sh                                  # the whole design
#   ./02-run-tools.sh --variant V                      # one pair
#
# IN   --plan FILE   design to run [OUTDIR/plan.tsv]   --variant V   one pair instead
#      --outdir DIR                results directory [./results]
#      --seed, --alpha             taken from the plan, or parsed from the variant name
#      -k, --kmer INT              k-mer length [23]
#      -l, --window INT            roll window in k-mers [333]
#      -s, --step INT              roll step in k-mers [100]
#      -p, --min-identity INT      wfmash minimum identity % [70]
#      -t, --threads INT           threads, BOTH tools [4]. gdiff's --num-threads is
#                                  global, so it precedes the subcommand. wfmash is NOT
#                                  reproducible at -t >1; gdiff is.
#      --prefix STR                query contig prefix [Q_]
#      --skip-existing             leave pairs whose roll.tsv already exists
# OUT  DIR/<variant>/  ref.fa qry.fa *.fai aln.paf ref.gs roll.tsv resources.tsv *.log
#
# Query ids are prefixed because wfmash's skip-self filter matches on sequence name.
# --- end help ---
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SIMDIR="${SIMDIR:-$(cd "$HERE/.." && pwd)}"
GDIFF="${GDIFF:-$(cd "$HERE/../.." && pwd)/bin/gdiff}"
WFMASH="${WFMASH:-$HERE/bin/wfmash}"
MAKEFAI="${MAKEFAI:-$HERE/lib/make-fai.py}"
TIMEIT="${TIMEIT:-$HERE/lib/timeit.py}"

VARIANT=""; PLAN=""; SEED=""; ALPHA=""; OUTDIR="./results"
K=23; L=333; S=100; PCT=70; THREADS=4; PREFIX="Q_"; SKIP=0

while [ $# -gt 0 ]; do
    case "$1" in
        --variant) VARIANT="$2"; shift 2 ;;
        --plan) PLAN="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --alpha) ALPHA="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        -k|--kmer) K="$2"; shift 2 ;;
        -l|--window) L="$2"; shift 2 ;;
        -s|--step) S="$2"; shift 2 ;;
        -p|--min-identity) PCT="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --skip-existing) SKIP=1; shift ;;
        -h|--help)
            # Check both markers first: a missing END marker makes sed run to EOF and print the
            # whole script, which is worse than printing nothing.
            if ! grep -qx '# --- help ---' "$0" || ! grep -qx '# --- end help ---' "$0"; then
                echo "error: $0 is missing its '# --- help ---' / '# --- end help ---' markers" >&2
                exit 1
            fi
            sed -n '/^# --- help ---$/,/^# --- end help ---$/p' "$0" | sed 's/^# \{0,1\}//; /^--- /d'
            exit 0
            ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

# wfmash's chaining races at -t >1: repeated runs of the SAME pair give different alignment
# counts, so its numbers are not reproducible. gdiff is deterministic at any thread count.
if [ "$THREADS" -gt 1 ]; then
    echo "note: wfmash is not reproducible at -t >1 (repeated runs differ). gdiff is." >&2
    echo "      Use -t 1 for results you need to reproduce; -t $THREADS is timing-only." >&2
fi

for tool in "$GDIFF" "$WFMASH"; do
    [ -x "$tool" ] || { echo "error: not executable: $tool" >&2; exit 1; }
    # A truncated (0-byte) file is still +x and runs as an empty script with exit 0, so it
    # fails silently rather than loudly. Check the size too.
    [ -s "$tool" ] || { echo "error: $tool is 0 bytes -- a truncated binary is not a tool" >&2; exit 1; }
done

run_one() {  # run_one <seed> <alpha> <variant>
    local seed="$1" alpha="$2" variant="$3"
    # declared before use: in `local a=1 b=$a`, $a is expanded before a is assigned.
    local seeddir ref vdir out
    seeddir="$SIMDIR/genomes/$seed"
    ref="$seeddir/${seed}_contigs.fasta"
    vdir="$seeddir/$variant"
    out="$OUTDIR/$variant"

    [ -f "$ref" ] || { echo "error: no seed FASTA at $ref" >&2; exit 1; }
    [ -f "$vdir/mutated.fasta" ] || { echo "error: no $vdir/mutated.fasta" >&2; exit 1; }

    if [ "$SKIP" = 1 ] && [ -s "$out/roll.tsv" ]; then
        echo "--- $variant : already run, skipping"
        return 0
    fi
    mkdir -p "$out"
    echo "--- $variant (seed $seed, alpha $alpha)"
    ln -sf "$ref" "$out/ref.fa"
    sed "s/^>/>${PREFIX}/" "$vdir/mutated.fasta" > "$out/qry.fa"
    # make-fai.py leaves an existing .fai alone, so clear them: both FASTAs are rewritten above
    # and a stale index would silently index the wrong sequence.
    rm -f "$out/ref.fa.fai" "$out/qry.fa.fai"
    python3 "$MAKEFAI" "$out/ref.fa" "$out/qry.fa" >/dev/null

    # wfmash writes its PAF to stdout, so a plain redirect captures it. timeit's own summary goes
    # to stderr, which is why no WFMASH_OUT-style hook is needed -- that would have tied this to
    # one patched build instead of any wfmash, on either platform.
    python3 "$TIMEIT" --tsv "$out/resources.tsv" --pair "$variant" \
        --label wfmash --threads "$THREADS" -- \
        "$WFMASH" -t "$THREADS" -p "$PCT" -B "$out" "$out/ref.fa" "$out/qry.fa" \
        > "$out/aln.paf" 2> "$out/wfmash.log"
    if [ ! -s "$out/aln.paf" ]; then
        echo "error: wfmash produced no alignments for $variant." >&2
        echo "       wfmash: $WFMASH" >&2
        sed 's/^/       /' "$out/wfmash.log" >&2
        exit 1
    fi
    # --num-threads is a GLOBAL gdiff option: it must precede the subcommand.
    python3 "$TIMEIT" --tsv "$out/resources.tsv" --pair "$variant" --label gdiff_sketch \
        --threads "$THREADS" -- \
        "$GDIFF" --num-threads "$THREADS" sketch -i "$out/ref.fa" -o "$out/ref.gs" \
        -k "$K" -w "$K" --frac 0.5 -l 0 > "$out/sketch.log" 2>&1
    python3 "$TIMEIT" --tsv "$out/resources.tsv" --pair "$variant" --label gdiff_roll \
        --threads "$THREADS" -- \
        "$GDIFF" --num-threads "$THREADS" roll "$out/qry.fa" "$out/ref.gs" -l "$L" -s "$S" \
        -o "$out/roll.tsv" > "$out/roll.log" 2>&1
    [ -f "$out/aln.paf" ] || {
        echo "error: wfmash exited 0 but wrote no $out/aln.paf (see $out/wfmash.log)" >&2
        exit 1
    }
    echo "    $(wc -l < "$out/aln.paf" | tr -d ' ') alignments, " \
         "$(($(wc -l < "$out/roll.tsv") - 3)) windows"
}

# Default job is the whole design; --variant switches to a single pair.
if [ -z "$VARIANT" ] && [ -z "$PLAN" ]; then
    PLAN="$OUTDIR/plan.tsv"
fi

if [ -n "$PLAN" ]; then
    [ -f "$PLAN" ] || { echo "error: no plan at $PLAN" >&2; exit 1; }
    n=0
    while IFS=$'\t' read -r p_seed p_alpha _p_level p_variant; do
        [ "$p_seed" = "seed" ] && continue       # header
        [ -n "${p_variant:-}" ] || continue
        run_one "$p_seed" "$p_alpha" "$p_variant"
        n=$((n + 1))
    done < "$PLAN"
    echo "stage 02: $n pairs"
    exit 0
fi

if [ -z "$SEED" ] || [ -z "$ALPHA" ]; then
    if [[ "$VARIANT" =~ ^mutated_BLOSUM62_(.+)_a([0-9]+)_gnd[0-9]+$ ]]; then
        [ -n "$SEED" ] || SEED="${BASH_REMATCH[1]}"
        [ -n "$ALPHA" ] || ALPHA="${BASH_REMATCH[2]}"
    else
        echo "error: cannot parse seed/alpha from $VARIANT; pass --seed and --alpha" >&2
        exit 1
    fi
fi
run_one "$SEED" "$ALPHA" "$VARIANT"
