#!/usr/bin/env bash
## Final recovery: S2 (tile-block tryCatch) + S6 (GetAssayData layer= fix).
source /Users/ikuz/miniforge3/etc/profile.d/conda.sh
conda activate rv_atlas
cd /Users/ikuz/Documents/RV_Atlas
export R_MAX_VSIZE=200Gb

echo "=== S2+S6 RECOVERY START: $(date) ==="
for s in Supplementary_Figure_2.R Supplementary_Figure_6.R; do
  log="/tmp/${s%.R}_run.log"
  echo ""; echo "----- $(date '+%H:%M:%S') START: $s -----"
  t0=$(date +%s); Rscript "$s" > "$log" 2>&1; rc=$?; t1=$(date +%s)
  if [ $rc -eq 0 ]; then
    echo "----- $(date '+%H:%M:%S') OK    : $s  ($((t1-t0))s)"
  else
    echo "----- $(date '+%H:%M:%S') FAIL  : $s  ($((t1-t0))s, rc=$rc)"
    grep -E "^Error|Execution halted" "$log" | tail -3 | sed 's/^/    /'
  fi
  echo "  disk: $(df -h ~/Documents | awk 'NR==2{print $4}')"
done
echo ""; echo "=== S2+S6 RECOVERY END: $(date) ==="
