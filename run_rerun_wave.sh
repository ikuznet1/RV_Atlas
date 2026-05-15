#!/usr/bin/env bash
## Re-run wave for scripts that failed in the first queue and have now been fixed.
## F8 is intentionally skipped (OOM on 48 GB Mac; user okayed skipping).

source /Users/ikuz/miniforge3/etc/profile.d/conda.sh
conda activate rv_atlas
cd /Users/ikuz/Documents/RV_Atlas

export R_MAX_VSIZE=200Gb

start_global=$(date +%s)
echo "========================================================"
echo "RERUN WAVE START: $(date)"
echo "  scripts: F3, F5, F7, S2, S6, S8"
echo "  excluded: F8 (OOM expected, skip per instruction)"
echo "========================================================"

run_one() {
  local script="$1"
  local logfile="/tmp/${script%.R}_run.log"
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
    grep -E "^Error|Execution halted" "$logfile" | tail -3 | sed 's/^/    /'
  fi
  local free=$(df -h ~/Documents | awk 'NR==2 {print $4}')
  echo "  disk free: $free"
}

# Re-run order:
# 1. F5 (fast — should now find cm_subclust_new_new.rds in dependencies)
# 2. F7 (medium — assay.use harmony fix)
# 3. F3 (medium — instrumented; will give traceback if it still fails)
# 4. S2 (medium — parse error fixed)
# 5. S6 (long, 4hr last time — parses fine now)
# 6. S8 (medium — same harmony fix as F7)

for s in Figure_5.R Figure_7.R Figure_3.R Supplementary_Figure_2.R Supplementary_Figure_6.R Supplementary_Figure_8.R; do
  run_one "$s"
done

end_global=$(date +%s)
echo ""
echo "========================================================"
echo "RERUN WAVE END: $(date)  (total $((end_global - start_global))s)"
echo "========================================================"
