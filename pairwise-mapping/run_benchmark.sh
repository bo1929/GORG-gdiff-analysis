#!/usr/bin/env bash
# Batch pairwise regional methods: each query vs references.
#
# Usage:
#   ./run_benchmark.sh <contig_dir> <queries.txt> [outdir]
#
# queries.txt: one SAG id per line
#
# Env (optional):
#   METHODS=gdiff,minimap,blast,mummer
#   REFS=refs.txt          # default: all genomes in contig_dir
#   JOBS=8 THREADS=4
#   SENSITIVE=1 FORCE=0
#   GDIFF=/path/to/gdiff WINDOW=1000 MAX_TARGET_SEQS=5000
#   GDIFF_SKETCH_ARGS=""   # extra args after: gdiff sketch -o <sketch>
#   GDIFF_MAP_ARGS="..."   # default: -l $((WINDOW/2)) -d <dists>
set -euo pipefail

CONTIG_DIR="${1:?Usage: $0 <contig_dir> <queries.txt> [outdir]}"
QUERIES_FILE="${2:?Usage: $0 <contig_dir> <queries.txt> [outdir]}"
OUTDIR="${3:-./benchmark_out}"

CONTIG_DIR="$(cd "$CONTIG_DIR" && pwd)"
OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"

METHODS="${METHODS:-gdiff,minimap,blast,mummer}"
JOBS="${JOBS:-8}"
THREADS="${THREADS:-4}"
SENSITIVE="${SENSITIVE:-0}"
FORCE="${FORCE:-0}"
WINDOW="${WINDOW:-1000}"
GDIFF_BIN="${GDIFF:-$(command -v gdiff || true)}"
GDIFF_SKETCH_ARGS="${GDIFF_SKETCH_ARGS:-}"
GDIFF_MAP_ARGS="${GDIFF_MAP_ARGS:--l $((WINDOW / 2)) -d 0.01 0.1 0.25 0.5 0.75 1.0 1.25 1.5}"
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-5000}"

mkdir -p "$OUTDIR"/cache "$OUTDIR"/gdiff "$OUTDIR"/minimap \
  "$OUTDIR"/blast "$OUTDIR"/mummer "$OUTDIR"/manifests

have() { case ",$METHODS," in *",$1,"*) return 0;; *) return 1;; esac; }
need() { command -v "$1" >/dev/null || { echo "missing: $1" >&2; exit 1; }; }
skip() { [ "$FORCE" != 1 ] && [ -s "$1" ]; }

fa_for() {
  # SAG id -> contig FASTA path
  f="$CONTIG_DIR/${1}_contigs.fasta"
  [ -f "$f" ] || { echo "missing FASTA for $1: $f" >&2; exit 1; }
  printf '%s\n' "$f"
}

# Load query ids
QUERIES=""
while read -r id _; do
  case "$id" in ""|\#*) continue;; esac
  QUERIES="$QUERIES $id"
done < "$QUERIES_FILE"

# Load ref ids
if [ -n "${REFS:-}" ]; then
  REFS_LIST=""
  while read -r id _; do
    case "$id" in ""|\#*) continue;; esac
    REFS_LIST="$REFS_LIST $id"
  done < "$REFS"
else
  REFS_LIST=""
  for f in "$CONTIG_DIR"/*_contigs.fasta; do
    [ -f "$f" ] || continue
    b=$(basename "$f" _contigs.fasta)
    REFS_LIST="$REFS_LIST $b"
  done
fi

# trim spaces
QUERIES=$(echo "$QUERIES" | xargs)
REFS_LIST=$(echo "$REFS_LIST" | xargs)
[ -n "$QUERIES" ] || { echo "no queries" >&2; exit 1; }
[ -n "$REFS_LIST" ] || { echo "no refs" >&2; exit 1; }

for q in $QUERIES; do fa_for "$q" >/dev/null; done
for r in $REFS_LIST; do fa_for "$r" >/dev/null; done

printf '%s\n' $QUERIES >"$OUTDIR/manifests/queries.txt"
printf '%s\n' $REFS_LIST >"$OUTDIR/manifests/refs.txt"
{
  printf 'query\tsubject\n'
  for q in $QUERIES; do
    for r in $REFS_LIST; do
      [ "$q" = "$r" ] && continue
      printf '%s\t%s\n' "$q" "$r"
    done
  done
} >"$OUTDIR/manifests/pairs.tsv"

nq=$(wc -l <"$OUTDIR/manifests/queries.txt" | tr -d ' ')
nr=$(wc -l <"$OUTDIR/manifests/refs.txt" | tr -d ' ')
echo "queries=$nq refs=$nr methods=$METHODS"
echo "outdir=$OUTDIR"

# Combined refs FASTA
REFS_FA="$OUTDIR/cache/refs_combined.fasta"
if ! skip "$REFS_FA"; then
  echo "writing $REFS_FA"
  : >"$REFS_FA"
  for r in $REFS_LIST; do
    cat "$(fa_for "$r")" >>"$REFS_FA"
    echo >>"$REFS_FA"
  done
fi

# gdiff sketch once
if have gdiff; then
  [ -n "$GDIFF_BIN" ] && [ -x "$GDIFF_BIN" ] || { echo "ERROR: set GDIFF=/path/to/gdiff" >&2; exit 1; }
  SKETCH="$OUTDIR/cache/refs.gdiff"
  if ! skip "$SKETCH"; then
    echo "gdiff sketch -> $SKETCH"
    # shellcheck disable=SC2086
    set -- "$GDIFF_BIN" sketch -o "$SKETCH" $GDIFF_SKETCH_ARGS
    for r in $REFS_LIST; do
      set -- "$@" -i "$(fa_for "$r")"
    done
    "$@"
  fi
fi

for q in $QUERIES; do
  QFA=$(fa_for "$q")
  echo "=== $q ==="

  if have gdiff; then
    out="$OUTDIR/gdiff/$q.tsv"
    if skip "$out"; then echo "  gdiff: skip"
    else
      echo "  gdiff"
      # shellcheck disable=SC2086
      "$GDIFF_BIN" map -q "$QFA" -i "$SKETCH" -o "$out" $GDIFF_MAP_ARGS
    fi
  fi

  if have minimap; then
    need minimap2
    out="$OUTDIR/minimap/$q.tsv"
    if skip "$out"; then echo "  minimap: skip"
    else
      echo "  minimap"
      set -- python3 "$HERE/minimap_blocks.py" -q "$QFA" -s "$REFS_FA" -o "$out" \
        -t "$THREADS" --workdir "$OUTDIR/minimap/.work_$q"
      [ "$SENSITIVE" = 1 ] && set -- "$@" --sensitive
      "$@"
    fi
  fi

  if have blast; then
    need blastn; need makeblastdb
    out="$OUTDIR/blast/$q.windows.tsv"
    if skip "$out"; then echo "  blast: skip"
    else
      echo "  blast"
      python3 "$HERE/sliding_blastn.py" -q "$QFA" -s "$REFS_FA" -o "$out" \
        -W "$WINDOW" -t "$THREADS" --by-subject --max-target-seqs "$MAX_TARGET_SEQS" \
        --workdir "$OUTDIR/blast/.work_$q"
    fi
  fi

  if have mummer; then
    need nucmer; need delta-filter; need show-coords
    mkdir -p "$OUTDIR/mummer/$q"
    joblist="$OUTDIR/mummer/$q/jobs.tsv"
    : >"$joblist"
    for r in $REFS_LIST; do
      [ "$q" = "$r" ] && continue
      out="$OUTDIR/mummer/$q/$r.tsv"
      skip "$out" && continue
      printf '%s\t%s\t%s\n' "$QFA" "$(fa_for "$r")" "$out" >>"$joblist"
    done
    nj=$(grep -c . "$joblist" 2>/dev/null || echo 0)
    echo "  mummer: $nj pairs (P=$JOBS)"
    [ "$nj" -eq 0 ] && continue
    export HERE SENSITIVE
    # shellcheck disable=SC2016
    cat "$joblist" | xargs -P "$JOBS" -n 3 bash -c '
      qfa="$1"; sfa="$2"; out="$3"
      work="$(dirname "$out")/.work_$(basename "$out" .tsv)"
      set -- python3 "$HERE/mummer_blocks.py" -q "$qfa" -s "$sfa" -o "$out" -t 1 --workdir "$work"
      [ "$SENSITIVE" = 1 ] && set -- "$@" --sensitive
      "$@"
    ' _
  fi
done

echo "done -> $OUTDIR"
