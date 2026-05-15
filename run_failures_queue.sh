#!/usr/bin/env bash
## Re-queue for scripts that failed in the first 2 waves and have new fixes:
##   F3 (parse-error: instrumented + safer .filter_to_panel)
##   F7 (harmony assay.use fix + TOM redirect)
##   S2 (Louvain cluster-8 guard)
##   S6 (parse fix; re-run)
##   S8 (ProjectModules group.by.vars dropped)
## F5 + F8 intentionally skipped (F5 just finished, F8 OOMs).

source /Users/ikuz/miniforge3/etc/profile.d/conda.sh
conda activate rv_atlas
cd /Users/ikuz/Documents/RV_Atlas

export R_MAX_VSIZE=200Gb

start_global=$(date +%s)
echo "========================================================"
echo "FAILURE-RECOVERY QUEUE START: $(date)"
echo "  scripts: F3 → F7 → S2 → S6 → S8"
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
  echo "  disk free: $(df -h ~/Documents | awk 'NR==2 {print $4}')"
}

for s in Figure_3.R Figure_7.R Supplementary_Figure_2.R Supplementary_Figure_6.R Supplementary_Figure_8.R; do
  run_one "$s"
done

end_global=$(date +%s)
echo ""
echo "========================================================"
echo "FAILURE-RECOVERY QUEUE END: $(date)  (total $((end_global - start_global))s)"
echo "========================================================"
