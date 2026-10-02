#!/usr/bin/env bash
# --- help ---
# RESOURCE BENCHMARK - gdiff and wfmash, 1 reference set against N, on ONE shared index.
#
#   ./resource-scaling.sh [options]
#
# Three cases, each timed against a single shared index:
#
#   1x1   1 cell      baseline
#   1xN   N cells     query scaling      (1 reference, N queries)
#   Nx1   N cells     reference scaling  (N references, 1 query)
#
# 1xN and Nx1 are both N cells, so comparing them shows which side drives the cost without the
# N^2 blow-up of an all-vs-all run. That is why there is no NxN case here.
#
# IN   --seeds LIST            seed pool to draw from [every seed in genomes/]
#      -n, --n INT             genomes to use [10]. If N exceeds the pool it is cycled, and the
#                              repeated copies get a c<K>_ contig prefix so names stay unique.
#      --alpha N               rate-heterogeneity alpha [22]
#      --level F               fixed target divergence for every genome [0.05]
#      --levels LIST           instead, draw each genome's level at random from this list, e.g.
#                              "0.001,0.01,0.05,0.10,0.18" -- a mix rather than one level. Each
#                              genome draws independently, so a cycled copy gets its own level.
#      --rng-seed INT          seed for that draw [1], so the mix is reproducible
#      -o, --out FILE          results TSV [./resource-scaling.tsv]
#      --log FILE              run log [<out>.log]
#      -k, --kmer INT          k-mer length [23]
#      -l, --window INT        roll window in k-mers [333]
#      -s, --step INT          roll step in k-mers [100]
#      -p, --min-identity INT  wfmash minimum identity % [70]
#      -t, --threads INT       threads, BOTH tools [4]
#
# OUT  exactly two files: the results TSV and the log. Everything else -- the FASTAs, the
#      indexes, the alignments -- goes to a temporary directory that is removed on exit.
#
#      results TSV columns:
#        case  tool  step  references  queries  cells  threads  wall_s  user_s  sys_s
#        peak_rss_mb  exit
#
# --num-threads is a GLOBAL gdiff option, so it precedes the subcommand. gdiff is deterministic
# at any thread count; wfmash is NOT at -t >1, so its alignment counts vary between runs.
# --- end help ---
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SIMDIR="${SIMDIR:-$(cd "$HERE/.." && pwd)}"
GDIFF="${GDIFF:-$(cd "$HERE/../.." && pwd)/bin/gdiff}"
WFMASH="${WFMASH:-$HERE/bin/wfmash}"
MAKEFAI="${MAKEFAI:-$HERE/lib/make-fai.py}"
TIMEIT="${TIMEIT:-$HERE/lib/timeit.py}"

SEEDS=""; NSET=""; ALPHA=22; LEVEL=0.05; LEVELS=""; RNG_SEED=1
OUT="./resource-scaling.tsv"; LOG=""
K=23; L=333; S=100; PCT=70; THREADS=4

while [ $# -gt 0 ]; do
    case "$1" in
        --seeds) SEEDS="${2//,/ }"; shift 2 ;;
        -n|--n) NSET="$2"; shift 2 ;;
        --alpha) ALPHA="$2"; shift 2 ;;
        --level) LEVEL="$2"; shift 2 ;;
        --levels) LEVELS="${2//,/ }"; shift 2 ;;
        --rng-seed) RNG_SEED="$2"; shift 2 ;;
        -o|--out) OUT="$2"; shift 2 ;;
        --log) LOG="$2"; shift 2 ;;
        -k|--kmer) K="$2"; shift 2 ;;
        -l|--window) L="$2"; shift 2 ;;
        -s|--step) S="$2"; shift 2 ;;
        -p|--min-identity) PCT="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -h|--help)
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

case "$THREADS" in
    ''|*[!0-9]*) echo "error: --threads must be a positive integer, got '$THREADS'" >&2; exit 1 ;;
esac
[ "$THREADS" -ge 1 ] || { echo "error: --threads must be >= 1" >&2; exit 1; }
if [ -n "$NSET" ]; then
    case "$NSET" in
        ''|*[!0-9]*) echo "error: --n must be a positive integer, got '$NSET'" >&2; exit 1 ;;
    esac
fi

for tool in "$GDIFF" "$WFMASH"; do
    [ -x "$tool" ] || { echo "error: not executable: $tool" >&2; exit 1; }
    # A truncated (0-byte) file is still +x and runs as an empty script with exit 0, so it fails
    # silently rather than loudly. Check the size too.
    [ -s "$tool" ] || { echo "error: $tool is 0 bytes -- a truncated binary is not a tool" >&2; exit 1; }
done

LOG="${LOG:-${OUT%.tsv}.log}"
mkdir -p "$(dirname "$OUT")" "$(dirname "$LOG")"

# Everything printed below goes to the log and to the terminal.
exec > >(tee "$LOG") 2>&1

# Scratch space: FASTAs, indexes, alignments. Removed on exit, so the only things left behind
# are the results TSV and the log.
TMP="$(mktemp -d "${TMPDIR:-/tmp}/resource-scaling.XXXXXX")"
cleanup() { [ -n "${TMP:-}" ] && [ -d "$TMP" ] && rm -rf "$TMP"; }
trap cleanup EXIT

RAW="$TMP/raw.tsv"
if [ -n "$LEVELS" ]; then
    echo "=== resource benchmark: 1 reference set vs N (alpha $ALPHA, mixed levels) ==="
else
    echo "=== resource benchmark: 1 reference set vs N (alpha $ALPHA, target level $LEVEL) ==="
fi
echo "results $OUT"
echo "log     $LOG"
echo "scratch $TMP (removed on exit)"
echo "threads $THREADS (wfmash -t, gdiff --num-threads)"

# ---------------------------------------------------------------- the pool of genomes
if [ -z "$SEEDS" ]; then
    SEEDS="$(cd "$SIMDIR/genomes" && ls -d */ 2>/dev/null | sed 's:/$::' | sort | tr '\n' ' ')"
fi
read -ra POOL <<< "$SEEDS"
POOLN=${#POOL[@]}
[ "$POOLN" -ge 1 ] || { echo "error: no seeds found under $SIMDIR/genomes" >&2; exit 1; }

N="${NSET:-10}"
[ "$N" -ge 2 ] || { echo "error: --n must be >= 2" >&2; exit 1; }
COPIES=$(( (N + POOLN - 1) / POOLN ))
echo "pool    $POOLN seed(s); using N=$N genome(s) in $COPIES pass(es) over the pool"
if [ "$COPIES" -gt 1 ]; then
    echo "        the pool is cycled: repeated copies get a c<K>_ contig prefix so names stay unique"
fi

# ---------------------------------------------------------------- levels, one per genome
# With --levels every genome draws its own target divergence, so the set is a mix; with --level
# they all get the same one. The draw is seeded, so a mix is reproducible.
if [ -n "$LEVELS" ]; then
    printf '%s\n' "$LEVELS" > "$TMP/levelpool.txt"
    python3 - "$TMP/levelpool.txt" "$RNG_SEED" "$N" > "$TMP/levels.txt" <<'PY'
import random, sys
pool = open(sys.argv[1]).read().split()
rng = random.Random(int(sys.argv[2]))
for _ in range(int(sys.argv[3])):
    print(rng.choice(pool))
PY
else
    : > "$TMP/levels.txt"
    for _ in $(seq 1 "$N"); do echo "$LEVEL" >> "$TMP/levels.txt"; done
fi
LEVEL_OF=()
while IFS= read -r lv; do LEVEL_OF+=("$lv"); done < "$TMP/levels.txt"

# ---------------------------------------------------------------- inputs
echo
echo "--- inputs"
if [ -n "$LEVELS" ]; then
    printf '    levels drawn from [%s] with rng seed %s: ' "$LEVELS" "$RNG_SEED"
    python3 -c '
import collections, sys
c = collections.Counter(sys.argv[1:])
print(", ".join(f"{k} x{v}" for k, v in sorted(c.items(), key=lambda kv: float(kv[0]))))
' "${LEVEL_OF[@]}"
fi

REF_FILES=()
for idx in $(seq 0 $((N - 1))); do
    k=$((idx / POOLN)); i=$((idx % POOLN))
    seed="${POOL[$i]}"; level="${LEVEL_OF[$idx]}"
    variant="$(python3 "$HERE/01-pick-variant.py" --simdir "$SIMDIR" --seed "$seed" \
        --alpha "$ALPHA" --level "$level" --cache "$TMP/.realised.tsv" 2>>"$TMP/pick.log")"
    ref="$SIMDIR/genomes/$seed/${seed}_contigs.fasta"
    qry="$SIMDIR/genomes/$seed/$variant/mutated.fasta"
    [ -f "$ref" ] || { echo "error: no $ref" >&2; exit 1; }
    [ -f "$qry" ] || { echo "error: no $qry" >&2; exit 1; }
    pfx=""; [ "$k" -gt 0 ] && pfx="c${k}_"
    sed "s/^>/>${pfx}/" "$ref" > "$TMP/ref.$idx.fa"
    sed "s/^>/>Q_${pfx}/" "$qry" > "$TMP/qry.$idx.fa"
    REF_FILES+=("$TMP/ref.$idx.fa")
    if [ "$N" -le 20 ] || [ $((idx % 10)) -eq 0 ] || [ "$idx" -eq $((N - 1)) ]; then
        echo "    $((idx + 1))/$N  level $level  $seed  ->  $variant${pfx:+  (copy $k)}"
    fi
done

# One reference and one query for the 1-side; all N for the N-side.
cp "$TMP/ref.0.fa" "$TMP/ref1.fa"
cp "$TMP/qry.0.fa" "$TMP/qry1.fa"
cat "${REF_FILES[@]}" > "$TMP/refN.fa"
: > "$TMP/qryN.fa"
for idx in $(seq 0 $((N - 1))); do cat "$TMP/qry.$idx.fa" >> "$TMP/qryN.fa"; done

# Duplicate contig names would silently merge references.
dupes=$(grep '^>' "$TMP/refN.fa" | sort | uniq -d | wc -l | tr -d ' ')
[ "$dupes" = 0 ] || { echo "error: $dupes duplicate contig names in refN.fa; refusing to combine" >&2; exit 1; }
python3 "$MAKEFAI" "$TMP/ref1.fa" "$TMP/qry1.fa" "$TMP/refN.fa" "$TMP/qryN.fa" >/dev/null
echo "        refN.fa $(du -h "$TMP/refN.fa" | cut -f1), qryN.fa $(du -h "$TMP/qryN.fa" | cut -f1)"

# ---------------------------------------------------------------- indexes, built once
# gdiff takes all N references in one `sketch -i` call (one container). wfmash has no single-shot
# index over several files, so the references are concatenated and indexed with -W; mapping then
# loads it back with -I.
echo
echo "--- indexes (built once, shared by every case)"
python3 "$TIMEIT" --tsv "$RAW" --pair index_1 --label gdiff_index --threads "$THREADS" -- \
    "$GDIFF" --num-threads "$THREADS" sketch -i "$TMP/ref1.fa" -o "$TMP/ref1.gs" \
    -k "$K" -w "$K" --frac 0.5 -l 0 > "$TMP/gdiff_index_1.log" 2>&1
python3 "$TIMEIT" --tsv "$RAW" --pair index_N --label gdiff_index --threads "$THREADS" -- \
    "$GDIFF" --num-threads "$THREADS" sketch -i "${REF_FILES[@]}" -o "$TMP/refN.gs" \
    -k "$K" -w "$K" --frac 0.5 -l 0 > "$TMP/gdiff_index_N.log" 2>&1
python3 "$TIMEIT" --tsv "$RAW" --pair index_1 --label wfmash_index --threads "$THREADS" -- \
    "$WFMASH" -W "$TMP/idx1" -t "$THREADS" -p "$PCT" -B "$TMP" "$TMP/ref1.fa" "$TMP/qry1.fa" \
    > "$TMP/wfmash_index_1.log" 2>&1
python3 "$TIMEIT" --tsv "$RAW" --pair index_N --label wfmash_index --threads "$THREADS" -- \
    "$WFMASH" -W "$TMP/idxN" -t "$THREADS" -p "$PCT" -B "$TMP" "$TMP/refN.fa" "$TMP/qryN.fa" \
    > "$TMP/wfmash_index_N.log" 2>&1

# ---------------------------------------------------------------- the three cases
echo
echo "--- mapping (queries x references)"
for cse in 1x1 1xN Nx1; do
    case "$cse" in
        1x1) ref=ref1; qry=qry1; gs=ref1.gs; idx=idx1 ;;
        1xN) ref=ref1; qry=qryN; gs=ref1.gs; idx=idx1 ;;
        Nx1) ref=refN; qry=qry1; gs=refN.gs; idx=idxN ;;
    esac
    python3 "$TIMEIT" --tsv "$RAW" --pair "$cse" --label gdiff_map --threads "$THREADS" -- \
        "$GDIFF" --num-threads "$THREADS" roll "$TMP/$qry.fa" "$TMP/$gs" -l "$L" -s "$S" \
        -o "$TMP/roll.$cse.tsv" > "$TMP/gdiff_map_$cse.log" 2>&1
    # WFMASH_OUT, not a redirect: timeit prints its own summary to stdout, which would otherwise
    # land in the PAF and be counted as an alignment.
    WFMASH_OUT="$TMP/aln.$cse.paf" python3 "$TIMEIT" --tsv "$RAW" --pair "$cse" \
        --label wfmash_map --threads "$THREADS" -- \
        "$WFMASH" -I "$TMP/$idx" -t "$THREADS" -p "$PCT" -B "$TMP" "$TMP/$ref.fa" "$TMP/$qry.fa" \
        > "$TMP/wfmash_map_$cse.log" 2>&1
    echo "    $cse  gdiff windows $(($(wc -l < "$TMP/roll.$cse.tsv") - 3))," \
         "wfmash alignments $(wc -l < "$TMP/aln.$cse.paf" | tr -d ' ')"
done

# ---------------------------------------------------------------- results
python3 - "$RAW" "$OUT" "$N" <<'PY'
import csv, sys

raw, out, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
CASES = {
    "index_1": ("index", 1, 0, 0),
    "index_N": ("index", n, 0, 0),
    "1x1": ("map", 1, 1, 1),
    "1xN": ("map", 1, n, n),
    "Nx1": ("map", n, 1, n),
}
FIELDS = ["case", "tool", "step", "references", "queries", "cells", "threads",
          "wall_s", "user_s", "sys_s", "peak_rss_mb", "exit"]

rows = []
with open(raw) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        _, refs, qrys, cells = CASES[r["pair"]]
        tool, step = r["label"].split("_", 1)
        rows.append({"case": r["pair"], "tool": tool, "step": step, "references": refs,
                     "queries": qrys, "cells": cells, "threads": r["threads"],
                     "wall_s": f"{float(r['wall_s']):.3f}",
                     "user_s": f"{float(r['user_s']):.3f}",
                     "sys_s": f"{float(r['sys_s']):.3f}",
                     "peak_rss_mb": f"{float(r['max_rss_mb']):.1f}",
                     "exit": r["exit"]})
order = {k: i for i, k in enumerate(CASES)}
rows.sort(key=lambda r: (order[r["case"]], r["tool"]))
with open(out, "w") as fh:
    fh.write("\t".join(FIELDS) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[f]) for f in FIELDS) + "\n")

get = {(r["case"], r["tool"]): r for r in rows}
wall = lambda case, tool: float(get[(case, tool)]["wall_s"])
rss = lambda case, tool: float(get[(case, tool)]["peak_rss_mb"])

print()
print("--- threads actually used")
for tool in ("gdiff", "wfmash"):
    ts = sorted({r["threads"] for r in rows if r["tool"] == tool})
    print(f"    {tool:<7} {', '.join(ts)}")
bad = [r for r in rows if r["exit"] != "0"]
if bad:
    print(f"    WARNING: {len(bad)} command(s) exited non-zero: "
          + ", ".join(f"{r['case']}/{r['tool']}" for r in bad))

print()
print("--- indexes, built once and shared")
print(f"    {'references':<12} {'gdiff':>18} {'wfmash':>18}")
for case, lab in (("index_1", "1"), ("index_N", str(n))):
    print(f"    {lab:<12} {wall(case, 'gdiff'):>9.2f} s {rss(case, 'gdiff'):>6.0f} MiB "
          f"{wall(case, 'wfmash'):>9.2f} s {rss(case, 'wfmash'):>6.0f} MiB")

print()
print("--- mapping (a cell is one query against one reference)")
print(f"    {'case':<5} {'cells':>6} | {'gdiff s':>9} {'s/cell':>8} | "
      f"{'wfmash s':>10} {'s/cell':>8}")
for case in ("1x1", "1xN", "Nx1"):
    c = CASES[case][3]
    g, w = wall(case, "gdiff"), wall(case, "wfmash")
    print(f"    {case:<5} {c:>6} | {g:>9.2f} {g / c:>8.4f} | {w:>10.2f} {w / c:>8.4f}")

print()
print(f"--- scaling against the 1x1 baseline, both at {n} cells")
for tool in ("gdiff", "wfmash"):
    b = wall("1x1", tool)
    print(f"    {tool:<7} baseline {b:7.2f} s | 1xN x{wall('1xN', tool) / b:<6.2f} "
          f"| Nx1 x{wall('Nx1', tool) / b:<6.2f}")

print()
print("--- reading it")
print("    1xN (1 reference, N queries) and Nx1 (N references, 1 query) are both N cells, so")
print("    the two columns separate query-side from reference-side cost: whichever is larger is")
print("    the side that dominates at this N. Per-cell cost falling with N is fixed")
print("    per-invocation cost being amortised, not a scaling penalty. Peak RSS is a per-process")
print("    high-water mark.")
print()
print(f"wrote {out}")
PY

echo
echo "--- results"
column -t "$OUT" 2>/dev/null || cat "$OUT"
