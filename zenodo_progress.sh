#!/usr/bin/env bash
# Monitor the Zenodo upload progress.
#   ./zenodo_progress.sh        one-shot snapshot
#   ./zenodo_progress.sh -w     live view, refreshes every 15s (Ctrl-C to quit)
LOG="/Users/ikuz/Documents/RV_Atlas/zenodo_upload.log"
TOTAL_GB=57.2
TOTAL_OBJ=45

snapshot() {
  local done gb run last cur lastmtime agemin sock pid
  done=$(grep -c '] OK ' "$LOG" 2>/dev/null)
  gb=$(grep -oE '  [0-9.]+ MB  ' "$LOG" 2>/dev/null | awk '{s+=$1} END{printf "%.1f", s/1024}')
  pid=$(pgrep -f zenodo_upload.py | tail -1)
  run=$([ -n "$pid" ] && echo "RUNNING (pid $pid)" || echo "STOPPED")
  pct=$(awk -v g="${gb:-0}" -v t="$TOTAL_GB" 'BEGIN{printf "%.0f", g/t*100}')
  # current in-flight object = the [n/45] in the most recent line, if not an OK
  last=$(grep -E '\] (OK|ERR)|Failed' "$LOG" 2>/dev/null | tail -1)
  cur=$(grep -oE '\[[0-9]+/45\] (OK|ERR) [^ ]+' "$LOG" 2>/dev/null | tail -1)
  lastmtime=$(stat -f %m "$LOG" 2>/dev/null)
  agemin=$(( ( $(date +%s) - ${lastmtime:-0} ) / 60 ))
  sock=""
  [ -n "$pid" ] && sock=$(lsof -nP -iTCP -a -p "$pid" 2>/dev/null | awk '/ESTABLISHED/{print $9; exit}')

  printf "Zenodo upload  %s\n" "$(date '+%H:%M:%S')"
  printf "  objects : %s / %s\n" "${done:-0}" "$TOTAL_OBJ"
  printf "  bytes   : %s GB / %s GB  (~%s%%)\n" "${gb:-0}" "$TOTAL_GB" "$pct"
  printf "  state   : %s\n" "$run"
  printf "  current : %s\n" "${cur:-n/a}"
  printf "  last log line %s min ago\n" "$agemin"
  [ -n "$sock" ] && printf "  live socket -> %s\n" "$sock"
  printf "  recent  : %s\n" "${last:-n/a}"
}

if [ "$1" = "-w" ]; then
  while true; do clear; snapshot; echo; echo "(refreshing every 15s; Ctrl-C to quit)"; sleep 15; done
else
  snapshot
fi
