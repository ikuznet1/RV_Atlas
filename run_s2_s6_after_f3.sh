#!/usr/bin/env bash
## Wait for the running Figure_3.R to finish, then run S2 -> S6 sequentially.
## Self-contained: polls for the F3 Rscript process; no external trigger needed.

source /Users/ikuz/miniforge3/etc/profile.d/conda.sh
conda activate rv_atlas
cd /Users/ikuz/Documents/RV_Atlas
export R_MAX_VSIZE=200Gb

echo "=== GUARD: waiting for Figure_3.R to finish ($(date)) ==="
# Rscript exec's the R binary, so the live process is
#   .../bin/exec/R --no-echo --no-restore --file=.../Figure_3.R
# Match on --file=...Figure_3.R, NOT the literal string "Rscript".
while pgrep -fl -- "--file=.*Figure_3\.R" > /dev/null 2>&1; do
  sleep 60
done
echo "=== GUARD: Figure_3.R no longer running ($(date)) — starting S2/S6 ==="

run_one() {
  local script="$1"
  local log="/tmp/${script%.R}_run.log"
  echo ""; echo "----- $(date '+%H:%M:%S') START: $script -----"
  local t0=$(date +%s)
  Rscript "$script" > "$log" 2>&1
  local rc=$?
  local t1=$(date +%s)
  if [ $rc -eq 0 ]; then
    echo "----- $(date '+%H:%M:%S') OK    : $script  ($((t1-t0))s)"
  else
    echo "----- $(date '+%H:%M:%S') FAIL  : $script  ($((t1-t0))s, rc=$rc)"
    grep -E "^Error|Execution halted" "$log" | tail -3 | sed 's/^/    /'
  fi
  echo "  disk: $(df -h ~/Documents | awk 'NR==2{print $4}')"
}

for s in Supplementary_Figure_2.R Supplementary_Figure_6.R; do
  run_one "$s"
done

echo ""; echo "=== S2/S6 (post-F3) DONE: $(date) ==="
