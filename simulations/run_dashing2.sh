#!/usr/bin/env bash
# Dashing2 (MinHash k-mer distance) on every variant vs its baseline genome.
# Output saved in the same format as the other estimator scripts:
#   results/ani-comparison/dashing2-estimates-selected/distances/dashing2-<cfg>.tsv
#   : method  param_setup  genome_a  genome_b  distance  ani_pct
# Uses dashing2 cmp --mash-distance (Poisson model, ANI-compatible like mash:
# ani_pct = (1 - distance) * 100). One invocation per seed (-F baseline,
# -Q variants) rather than one process per pair; JOBS parallelizes seeds.
#
# MODE selects the parameter configuration (same set as methods/dashing2.sh):
#   v1  --mash-distance -k 31 -S 1024
#   v2  --containment                         (ANI = 100*(1+ln(C)/32))
#   v3  --symmetric-containment -k 31 -S 1024
#   v4  --symmetric-containment -k 23 -S 2048
#   v5  --mash-distance -k 21 -S 5000
#
# Config via env:
#   DASHING2     dashing2 binary          (default: ../bin/dashing2, arch
#                                         dispatcher, then PATH)
#   GENOMES      genome dir               (default: genomes)
#   SUFFIX       baseline suffix          (default: _contigs.fasta)
#   ONLY         comma-separated seed subset (default: all)
#   JOBS         parallel seed workers    (default: 8)
#   THREADS      dashing2 threads         (default: 4)
#   SKETCH_ARGS  extra dashing2 sketch/cmp args (default: "");
#                e.g. "-k 31" for custom kmer (defaults used if empty)
#   MODE         config preset: v1|v2|v3|v4|v5 (default v1)
#   CONFIG       short label for filename (default: derived from MODE)
#   OUTDIR       results dir (default: results/)
#   NO_PROGRESS  set to 1 to disable the progress bar
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"

DASHING2="${DASHING2:-}"
GENOMES="${GENOMES:-genomes}"
SUFFIX="${SUFFIX:-_contigs.fasta}"
ONLY="${ONLY:-}"
JOBS="${JOBS:-8}"
THREADS="${THREADS:-4}"
MODE="${MODE:-v4}"
SKETCH_ARGS="${SKETCH_ARGS:-}"
CONFIG="${CONFIG:-}"
OUTDIR="${OUTDIR:-results/}"
NO_PROGRESS="${NO_PROGRESS:-0}"

# Resolve dashing2 binary (bin/dashing2 arch dispatcher, then PATH).
DASHING2="${DASHING2:-}"
[ -z "$DASHING2" ] && [ -x "$REPO_ROOT/bin/dashing2" ] && DASHING2="$REPO_ROOT/bin/dashing2"
[ -z "$DASHING2" ] && command -v dashing2 >/dev/null 2>&1 && DASHING2="$(command -v dashing2)"
[ -x "${DASHING2:-}" ] || { echo "missing dashing2 binary (set DASHING2=)" >&2; exit 1; }

# Same parameter configurations as methods/dashing2.sh (v1-v5).
case "$MODE" in
  v1) DIST_FLAG="--mash-distance";         EXTRA="-k 31 -S 1024";      KMER=31; METRIC="mash";          CFG="v1"; ;;
  v2) DIST_FLAG="--containment";           EXTRA="";                   KMER=32; METRIC="containment";   CFG="v2"; ;;
  v3) DIST_FLAG="--symmetric-containment"; EXTRA="-k 31 -S 1024";      KMER=31; METRIC="containment";   CFG="v3"; ;;
  v4) DIST_FLAG="--symmetric-containment"; EXTRA="-k 23 -S 2048";      KMER=23; METRIC="containment";   CFG="v4"; ;;
  v5) DIST_FLAG="--mash-distance";         EXTRA="-k 21 -S 5000";      KMER=21; METRIC="mash";          CFG="v5"; ;;
  *) echo "unknown MODE=$MODE (v1|v2|v3|v4|v5)" >&2; exit 1 ;;
esac
[ -n "$SKETCH_ARGS" ] && EXTRA="$SKETCH_ARGS"
CONFIG="${CONFIG:-$CFG}"
PARAMS="${EXTRA:-dashing-default}"
DIST="$OUTDIR/distances/dashing2-${CONFIG:-$MODE}.tsv"
mkdir -p "$OUTDIR"/distances

# Count selected seeds for progress.
count_seeds() {
    local seed n=0
    for seed in "$GENOMES"/*/; do
        seed=$(basename "$seed")
        [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
        if [ -n "$ONLY" ]; then case ",$ONLY," in *",$seed,"*) ;; *) continue ;; esac; fi
        n=$((n+1))
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
mkdir -p "$TMP/e"

printf 'method\tparam_setup\tgenome_a\tgenome_b\tdistance\tani_pct\n' > "$DIST"

# run_dashing2_seed: build query list (variants) + ref list (baseline), run
# dashing2 cmp once, emit variant -> baseline rows to stdout.
run_dashing2_seed() {
    local seed="$1"
    local ql="$TMP/ql_$seed" rl="$TMP/rl_$seed" keylist="$TMP/keys_$seed" out="$TMP/cmp_$seed.txt"
    : > "$ql"; : > "$rl"; : > "$keylist"
    printf '%s\n' "$(pwd)/$GENOMES/$seed/$seed$SUFFIX" > "$rl"

    for d in "$GENOMES/$seed"/mutated_BLOSUM62_*; do
        [ -d "$d" ] || continue; [ -f "$d/mutated.fasta" ] || continue
        core="$(basename "$d")"; core="${core#mutated_BLOSUM62_}"
        printf '%s\n' "$(pwd)/$d/mutated.fasta" >> "$ql"
        printf '%s\n' "$core" >> "$keylist"
        for vf in "$d"/mutated_p*_s*.fasta; do
            [ -f "$vf" ] || continue
            f="$(basename "$vf")"; f="${f#mutated_}"; f="${f%.fasta}"
            printf '%s\n' "$(pwd)/$vf" >> "$ql"
            printf '%s\n' "${core}_${f}" >> "$keylist"
        done
    done
    [ -s "$ql" ] || return 1

    # rectangular matrix: one row per ref, tab-separated per query
    # shellcheck disable=SC2086
    "$DASHING2" cmp -F "$rl" -Q "$ql" $DIST_FLAG $EXTRA -p "$THREADS" 2>/dev/null \
        | grep -v '^[#Dashing2]' > "$out"

    # map: row 1 only (single ref = this seed's baseline); columns = queries in ql order
    local base="$seed$SUFFIX"
    if [ -s "$out" ]; then
        awk -F'\t' -v metric="$METRIC" -v params="$PARAMS" -v kk="$KMER" -v base="$base" -v seed="$seed" -v kf="$keylist" '
            BEGIN { ki=1; while ((getline k < kf) > 0) KEY[ki++] = k; close(kf) }
            NR==1 {
                for (i=2; i<=NF; i++) {
                    if ($i != "-" && $i != "NaN") {
                        v = $i + 0
                        if (metric == "containment") {
                            ani = 100 * (1 + log(v)/kk)
                            printf "dashing2\t%s\t%s\t%s\t%.6f\t%.6f\n", params, KEY[i-1], seed, 1-ani/100, ani
                        } else {
                            printf "dashing2\t%s\t%s\t%s\t%.6f\t%.6f\n", params, KEY[i-1], seed, v, (1-v)*100
                        }
                    }
                }
            }' "$out"
    fi
}
export -f run_dashing2_seed
export DASHING2 GENOMES SUFFIX ONLY THREADS DIST_FLAG EXTRA METRIC KMER PARAMS TMP

IDX=0
DONE=0
SEEDS=$(count_seeds)
running=0
for seed in "$GENOMES"/*/; do
    seed=$(basename "$seed")
    [ -f "$GENOMES/$seed/$seed$SUFFIX" ] || continue
    if [ -n "$ONLY" ]; then case ",$ONLY," in *",$seed,"*) ;; *) continue ;; esac; fi
    echo "$seed" >&2
    run_dashing2_seed "$seed" > "$TMP/e/$IDX" &
    IDX=$((IDX+1)); running=$((running+1))
    if (( running % JOBS == 0 )); then wait; running=0; DONE="$IDX"; progress "$DONE" "$SEEDS"; fi
done
wait
DONE="$IDX"; progress "$DONE" "$SEEDS"

for i in $(seq 0 $((IDX-1))); do cat "$TMP/e/$i" 2>/dev/null; done >> "$DIST"

echo "wrote $DIST"
