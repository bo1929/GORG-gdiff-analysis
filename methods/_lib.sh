# methods/_lib.sh — sourced by method drivers (not run directly).
# Conventions:
#   CLI:  <genome_dir> <pairs.tsv> [outdir] [suffix=.fasta]
#   env:  JOBS (pair parallelism, default 8)  THREADS (tool threads)
#         FORCE=0  ONLY=all|name[,name...]
#   out:  $OUT/distances/<method>-<cfg>.tsv  (+ all_<method>.tsv)
#         $OUT/blocks/<method>/<cfg>/<q>__<s>.tsv
#         $OUT/cache/<method>/   (pairs.tsv, genomes.txt, reusable sketches only)
# No per-pair cache files. Temps live under $TMPDIR and are removed.

want() { [ "${ONLY:-all}" = all ] && return 0; case ",$ONLY," in *",$1,"*) return 0;; esac; return 1; }

fa() {
  local f="${GENOME_DIR}/$1${SUFFIX}"
  [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }
  echo "$f"
}

# Batch-wait every JOBS launches (bash 3.2-safe; no wait -n).
throttle() { local n="$1"; if (( n > 0 && n % JOBS == 0 )); then wait; fi; }

# In-place progress on stderr. Call with done/total/label; newline when done==total.
progress() {
  local done="$1" total="$2" label="$3" width=28 pct filled empty i bar
  [ "$total" -gt 0 ] || return 0
  pct=$(( done * 100 / total ))
  filled=$(( done * width / total )); [ "$filled" -gt "$width" ] && filled=$width
  empty=$(( width - filled )); bar=""
  i=0; while [ "$i" -lt "$filled" ]; do bar="${bar}="; i=$((i+1)); done
  if [ "$done" -lt "$total" ] && [ "$filled" -lt "$width" ]; then bar="${bar}>"; empty=$((empty-1)); fi
  i=0; while [ "$i" -lt "$empty" ]; do bar="${bar} "; i=$((i+1)); done
  printf '\r%-8s [%s] %d/%d (%d%%)' "$label" "$bar" "$done" "$total" "$pct" >&2
  if [ "$done" -ge "$total" ]; then echo >&2; fi
  return 0
}

# Load pairs -> $CACHE/pairs.tsv; unique ids -> $CACHE/genomes.txt (ids only).
load_pairs() {
  mkdir -p "$CACHE"
  grep -v '^#' "$PAIRS_FILE" | awk 'NF>=2' > "$CACHE/pairs.tsv"
  cut -f1,2 "$CACHE/pairs.tsv" | tr '\t' '\n' | sort -u > "$CACHE/genomes.txt"
}

# Concat selected per-config distance TSVs into all_<method>.tsv
emit_all_distances() {
  local method="$1" hdr="$2"; shift 2
  local name tsv
  { echo "$hdr"
    for c in "$@"; do
      IFS='|' read -r name _ <<< "$c"
      tsv="$DIST_DIR/${method}-$name.tsv"
      if want "$name" && [ -s "$tsv" ]; then tail -n+2 "$tsv"; fi
    done
  } > "$DIST_DIR/all_${method}.tsv"
  return 0
}
