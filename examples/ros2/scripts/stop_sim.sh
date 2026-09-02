#!/usr/bin/env bash
# stop_sim.sh : kill exactly the processes start_sim.sh recorded.
#
# Deliberately does not pattern-match on process names. A pkill for "gz" or "px4"
# happily matches the terminal you launched from, and a pkill that catches its own
# parent shell is an unpleasant way to find that out.
set -uo pipefail
RUNDIR="${RUNDIR:-/tmp/sheaf_sim}"
[[ -f "$RUNDIR/pids" ]] || { echo "nothing recorded in $RUNDIR/pids"; exit 0; }
while read -r pid; do
  [[ -z "${pid:-}" ]] && continue
  if kill -0 "$pid" 2>/dev/null; then
    echo "[stop_sim] terminating $pid"
    kill "$pid" 2>/dev/null || true
  fi
done < "$RUNDIR/pids"
sleep 2
while read -r pid; do
  [[ -z "${pid:-}" ]] && continue
  kill -0 "$pid" 2>/dev/null && kill -9 "$pid" 2>/dev/null || true
done < "$RUNDIR/pids"
: > "$RUNDIR/pids"
echo "[stop_sim] done"
