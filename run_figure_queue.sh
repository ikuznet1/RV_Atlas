#!/usr/bin/env bash
## Sequential queue: Figure_3..Figure_8 then Supplementary_Figure_1..8.
## Each Rscript invocation has its own log file under /tmp/.
## Continues past failures so we capture all errors in one pass.

## (no `set -u` — conda's activate-gfortran script references unset GFORTRAN)
## (no `set -e` — we want to continue past Rscript failures)

source /Users/ikuz/miniforge3/etc/profile.d/conda.sh
conda activate rv_atlas
cd /Users/ikuz/Documents/RV_Atlas

start_global=$(date +%s)
echo "========================================================"
echo "FIGURE QUEUE START: $(date)"
echo "========================================================"

run_one() {
  local script="$1"
  local logfile="/tmp/${script%.R}_run.log"
  local tag="${script%.R}"
  echo ""
  echo "----- $(date '+%H:%M:%S') START: $script -----"
  echo "  log: $logfile"
  local t0=$(date +%s)
  Rscript "$script" > "$logfile" 2>&1
  local rc=$?
  local t1=$(date +%s)
  local dur=$((t1 - t0))
  if [ $rc -eq 0 ]; then
    echo "----- $(date '+%H:%M:%S') OK    : $script  (${dur}s, rc=$rc)"
  else
    echo "----- $(date '+%H:%M:%S') FAIL  : $script  (${dur}s, rc=$rc)"
    echo "  last error:"
    grep -E "^Error|Execution halted" "$logfile" | tail -3 | sed 's/^/    /'
  fi
  # Disk check (warn if low)
  local free=$(df -h ~/Documents | awk 'NR==2 {print $4}')
  echo "  disk free: $free"
}

# Phase 1: Figures 3..8
for n in 3 4 5 6 7 8; do
  run_one "Figure_${n}.R"
done

# Phase 2: Supplementary 1..8
for n in 1 2 3 4 5 6 7 8; do
  run_one "Supplementary_Figure_${n}.R"
done

end_global=$(date +%s)
total=$((end_global - start_global))
echo ""
echo "========================================================"
echo "FIGURE QUEUE END: $(date)  (total ${total}s)"
echo "========================================================"
