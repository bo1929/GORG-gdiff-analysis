#!/usr/bin/env bash
# All-by-all resource benchmark: one config per method, sketch then dist.
# Edit the defaults below, then: ./run_benchmarks.sh [method ...]
# Everything (sample, sketches, distances, timings) goes in output/.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

# ---- defaults (edit these) ------------------------------------------------
THREADS=16
SAMPLE_N=100
SEED=79
FORCE=1 # 1 = redraw sample and redo selected methods

PAIRS_TSV="$REPO/all_pairs.tsv"
GENOME_DIR="$REPO/contigs-gt80-complete"
GENOME_SUFFIX="_contigs.fasta"
OUT="$HERE/output"

DASHING2="$REPO/bin/dashing2"
# GDIFF="$REPO/bin/gdiff"
GDIFF="../../gdiff/gdiff"
MASH=mash
SKANI=skani
FASTANI=fastani # also accepts fastANI
ANIB=anib

DASHING2_ARGS="--symmetric-containment -k 23 -S 2048"
SKANI_ARGS="--slow"
GDIFF_SKETCH_ARGS="-k 23 -w 23 --frac 0.5 -l 500 --sample-size 1000"
GDIFF_DIST_ARGS="--hdist-th 3"
MASH_ARGS="-k 19 -s 10000"
FASTANI_ARGS="--fragLen 3000 --minFraction 0.1"

METHODS=(dashing2 skani mash gdiff fastani)
# --------------------------------------------------------------------------

GENOME_LIST="$OUT/genomes.list"
GDIFF_LIST="$OUT/genomes.gdiff.list"
SAMPLE_TSV="$OUT/sample_genomes.tsv"
SAMPLE_PAIRS="$OUT/sample_pairs.tsv"
RES="$OUT/resources.tsv"

if [ "$#" -gt 0 ]; then
  METHODS=("$@")
fi

case "$(uname -s)" in
  Darwin) TIME_BIN=/usr/bin/time; TIME_FLAGS=-l ;;
  *)      TIME_BIN=/usr/bin/time; TIME_FLAGS=-v ;;
esac

# ---- helpers --------------------------------------------------------------

sample_genomes() {
  mkdir -p "$OUT"
  [ -s "$PAIRS_TSV" ] || { echo "missing pairs file: $PAIRS_TSV" >&2; exit 1; }
  [ -d "$GENOME_DIR" ] || { echo "missing genome dir: $GENOME_DIR" >&2; exit 1; }

  local header="# n=$SAMPLE_N seed=$SEED"
  if [ "$FORCE" != 1 ] && [ -s "$SAMPLE_TSV" ] && [ "$(head -1 "$SAMPLE_TSV")" = "$header" ]; then
    echo "keeping $SAMPLE_TSV" >&2
  else
    echo "drawing $SAMPLE_N genomes from $PAIRS_TSV (seed=$SEED)" >&2
    python3 - "$PAIRS_TSV" "$GENOME_DIR" "$GENOME_SUFFIX" "$SAMPLE_N" "$SEED" \
      "$SAMPLE_TSV" "$SAMPLE_PAIRS" <<'PY'
import os, random, sys

pairs_path, genome_dir, suffix = sys.argv[1], sys.argv[2], sys.argv[3]
n, seed = int(sys.argv[4]), int(sys.argv[5])
out_ids, out_pairs = sys.argv[6], sys.argv[7]

ids = set()
with open(pairs_path) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            ids.add(f[0]); ids.add(f[1])
ids = sorted(i for i in ids if os.path.isfile(os.path.join(genome_dir, i + suffix)))
if len(ids) < n:
    sys.exit(f"only {len(ids)} genomes available in {genome_dir}, need {n}")

rng = random.Random(seed)
sample = sorted(rng.sample(ids, n))

with open(out_ids, "w") as fh:
    fh.write(f"# n={n} seed={seed}\n")
    for g in sample:
        fh.write(g + "\n")
with open(out_pairs, "w") as fh:
    fh.write(f"# all-by-all of {n} genomes (n={n} seed={seed})\n")
    fh.write("# query\tsubject\n")
    for i in range(n):
        for j in range(i + 1, n):
            fh.write(f"{sample[i]}\t{sample[j]}\n")
print(f"wrote {out_ids} ({n} genomes)", file=sys.stderr)
print(f"wrote {out_pairs} ({n*(n-1)//2} pairs)", file=sys.stderr)
PY
  fi

  : > "$GENOME_LIST"
  local id
  while read -r id; do
    [ -n "$id" ] || continue
    case "$id" in \#*) continue ;; esac
    printf '%s/%s%s\n' "$GENOME_DIR" "$id" "$GENOME_SUFFIX" >> "$GENOME_LIST"
  done < "$SAMPLE_TSV"

  local missing=0 f
  while read -r f; do
    [ -f "$f" ] || { echo "missing FASTA: $f" >&2; missing=1; }
  done < "$GENOME_LIST"
  [ "$missing" -eq 0 ] || exit 1

  awk -F/ -v OFS='\t' '{ n=split($NF, a, "."); print a[1], $0 }' "$GENOME_LIST" \
    > "$GDIFF_LIST"
}

# Parse /usr/bin/time -l (macOS) or -v (Linux) into wall user sys maxrss_mib.
parse_time() {
  awk '
    BEGIN { wall="."; user="."; sys="."; mb="." }
    /^[[:space:]]*[0-9]+([.][0-9]+)?[[:space:]]+real[[:space:]]/ {
      wall=$1; user=$3; sys=$5
    }
    /^[[:space:]]*[0-9]+([.][0-9]+)?[[:space:]]+maximum resident set size/ {
      mb = $1/1048576
    }
    /^[[:space:]]*Maximum resident set size \(kbytes\):/ { mb = $NF/1024 }
    /^[[:space:]]*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):/ {
      t=$NF; n=split(t, a, ":")
      wall = (n==3) ? a[1]*3600 + a[2]*60 + a[3] : a[1]*60 + a[2]
    }
    /^[[:space:]]*User time \(seconds\):/ { user=$NF }
    /^[[:space:]]*System time \(seconds\):/ { sys=$NF }
    END { printf "%s\t%s\t%s\t%s\n", wall, user, sys, mb }
  ' "$1"
}

# measure method config phase cmd...
measure() {
  local method="$1" cfg="$2" phase="$3"; shift 3
  local raw="$OUT/$method.$cfg.$phase.time" mrc=0 metrics
  # shellcheck disable=SC2086
  "$TIME_BIN" $TIME_FLAGS "$@" \
    > "$OUT/$method.$cfg.$phase.out" \
    2> "$raw" || mrc=$?
  metrics=$(parse_time "$raw" 2>/dev/null || true)
  [ -n "$metrics" ] || metrics=$'.\t.\t.\t.'
  printf '%s\t%s\t%s\t%s\t%s\n' "$method" "$cfg" "$phase" "$THREADS" "$metrics" >> "$RES"
  if [ "$mrc" -ne 0 ]; then
    echo "  [$method/$cfg/$phase] exited $mrc (see $raw)" >&2
  fi
  return 0
}

ensure_header() {
  if [ ! -s "$RES" ]; then
    printf 'method\tconfig\tphase\tthreads\twall_sec\tuser_sec\tsys_sec\tmaxrss_mib\n' > "$RES"
  fi
}

# Drop existing rows for a method so a re-run replaces them.
drop_method_rows() {
  local m="$1" tmp="$OUT/.resources.tmp"
  awk -F'\t' -v m="$m" 'NR==1 || $1!=m' "$RES" > "$tmp"
  mv "$tmp" "$RES"
}

# Skip a method when FORCE=0 and it already has a row in resources.tsv.
should_run() {
  local m="$1"
  if [ "$FORCE" = 1 ]; then
    drop_method_rows "$m"
    return 0
  fi
  if awk -F'\t' -v m="$m" 'NR>1 && $1==m { f=1 } END { exit !f }' "$RES"; then
    echo "  [$m] already in $RES (set FORCE=1 to redo)" >&2
    return 1
  fi
  return 0
}

# ---- methods --------------------------------------------------------------

run_dashing2() {
  local cfg=v4 sk="$OUT/dashing2.sketches" dist="$OUT/dashing2.v4.dist.tsv"
  rm -rf "$sk"; mkdir -p "$sk"
  echo "[dashing2/$cfg] sketch -> $sk" >&2
  # shellcheck disable=SC2086
  measure dashing2 "$cfg" sketch \
    "$DASHING2" sketch $DASHING2_ARGS -p "$THREADS" \
      --cache-sketches --outprefix "$sk" -F "$GENOME_LIST"
  echo "[dashing2/$cfg] dist ($NP pairs)" >&2
  # shellcheck disable=SC2086
  measure dashing2 "$cfg" dist \
    "$DASHING2" cmp $DASHING2_ARGS -p "$THREADS" \
      --cache-sketches --outprefix "$sk" \
      -F "$GENOME_LIST" -Q "$GENOME_LIST" --cmpout "$dist"
  echo "  -> $dist" >&2
}

run_skani() {
  local cfg=slow sk="$OUT/skani.sketches" sl="$OUT/skani.sketch.list"
  local dist="$OUT/skani.slow.dist.tsv"
  rm -rf "$sk"
  echo "[skani/$cfg] sketch -> $sk" >&2
  # --separate-sketches: triangle reads per-genome .sketch files.
  # shellcheck disable=SC2086
  measure skani "$cfg" sketch \
    "$SKANI" sketch -t "$THREADS" $SKANI_ARGS --separate-sketches \
      -l "$GENOME_LIST" -o "$sk"
  find "$sk" -name '*.sketch' | sort > "$sl"
  echo "[skani/$cfg] dist ($NP pairs)" >&2
  # shellcheck disable=SC2086
  measure skani "$cfg" dist \
    "$SKANI" triangle -t "$THREADS" $SKANI_ARGS --min-af 0 -l "$sl" -o "$dist"
  echo "  -> $dist" >&2
}

run_gdiff() {
  local cfg=abcs bundle="$OUT/gdiff.abcs.gdsk" dist="$OUT/gdiff.abcs.dist.tsv"
  rm -f "$bundle"
  echo "[gdiff/$cfg] sketch -> $bundle" >&2
  # shellcheck disable=SC2086
  measure gdiff "$cfg" sketch \
    "$GDIFF" --num-threads "$THREADS" sketch \
      --input-list "$GDIFF_LIST" $GDIFF_SKETCH_ARGS -o "$bundle"
  echo "[gdiff/$cfg] dist ($NP pairs)" >&2
  # shellcheck disable=SC2086
  measure gdiff "$cfg" dist \
    "$GDIFF" --num-threads "$THREADS" dist "$bundle" $GDIFF_DIST_ARGS -o "$dist"
  echo "  -> $dist" >&2
}

run_mash() {
  local cfg=sensitive msh="$OUT/mash.msh" dist="$OUT/mash.sensitive.dist.tsv"
  rm -f "$msh"
  echo "[mash/$cfg] sketch -> $msh" >&2
  # shellcheck disable=SC2086
  measure mash "$cfg" sketch \
    "$MASH" sketch -p "$THREADS" $MASH_ARGS -l -o "${msh%.msh}" "$GENOME_LIST"
  echo "[mash/$cfg] dist ($NP pairs)" >&2
  measure mash "$cfg" dist \
    sh -c "'$MASH' dist -p '$THREADS' '$msh' '$msh' > '$dist'"
  echo "  -> $dist" >&2
}

run_fastani() {
  local cfg=frag3000 dist="$OUT/fastani.frag3000.dist.tsv" bin="$FASTANI"
  if ! command -v "$bin" >/dev/null; then
    if [ "$bin" = "fastani" ] && command -v fastANI >/dev/null; then
      bin=fastANI
    else
      echo "fastani not found: $FASTANI (also tried fastANI)" >&2
      return 1
    fi
  fi
  echo "[fastani/$cfg] dist ($NP pairs)" >&2
  # shellcheck disable=SC2086
  measure fastani "$cfg" dist \
    "$bin" --ql "$GENOME_LIST" --rl "$GENOME_LIST" \
      -t "$THREADS" $FASTANI_ARGS -o "$dist"
  echo "  -> $dist" >&2
}

# ---- run ------------------------------------------------------------------

sample_genomes
ensure_header
NG=$(grep -cv '^#' "$GENOME_LIST" || true)
NP=$(( NG * (NG - 1) / 2 ))
echo "threads=$THREADS pairs=$NP results=$RES" >&2

rc=0
for m in "${METHODS[@]}"; do
  should_run "$m" || continue
  case "$m" in
    dashing2) run_dashing2 ;;
    skani)    run_skani ;;
    mash)     run_mash ;;
    gdiff)    run_gdiff ;;
    fastani)  run_fastani || rc=1 ;;
    *) echo "unknown method: $m (want: dashing2 skani mash gdiff fastani anib)" >&2; rc=1 ;;
  esac
done

echo >&2
echo "== $RES ==" >&2
column -t -s $'\t' "$RES" >&2 || cat "$RES" >&2
exit "$rc"
