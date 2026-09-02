# lib.sh : preflight checks shared by smoke_test.sh and watch.sh. Source, do not run.

# Refuse to start when another run already owns the ZMQ endpoint.
#
# Every container uses --network host, so two runs share one port namespace AND,
# through the repo bind mount, one RUNDIR. The second brain dies with "Address
# already in use", its scripts then read the FIRST run's brain.log, see the
# "sheaf server" banner, conclude all is well, and let their robots talk to the
# other run's controller. The result looks like a plausible but failing
# experiment rather than a collision, which is a very expensive way to lose an
# afternoon. Fail loudly here instead.
require_free_zmq_port() {
  local endpoint="$1" port
  port="${endpoint##*:}"
  # Probe by connecting, not by parsing `ss`: the ros base image does not ship
  # iproute2, and a listing tool that is missing makes a grep-based check pass
  # silently -- which is exactly how a stale server got to fly someone else's
  # robots. bash's /dev/tcp needs nothing installed and cannot be absent.
  if (exec 3<>"/dev/tcp/127.0.0.1/${port}") 2>/dev/null; then
    exec 3>&- 3<&-
    cat >&2 <<MSG
ERROR: something is already listening on port ${port} (the sheaf ZMQ endpoint).

Another run is almost certainly still up. Because every container shares the
host network and this repo's .run directory, starting a second one silently
crosses the two: robots from one run get commands from the other's controller.

  docker ps --filter ancestor=sheaf-apark:humble
  docker rm -f <name>

Then try again.
MSG
    return 1
  fi
  return 0
}

# Wipe the previous run's artifacts so a stale brain.log can never be mistaken
# for this run's. Cheap, and it removes the other half of the collision above.
reset_rundir() {
  local rundir="$1"
  mkdir -p "$rundir"
  rm -f "$rundir"/brain.log "$rundir"/bridge.log "$rundir"/px4_*.log \
        "$rundir"/state.csv "$rundir"/gz.log "$rundir"/pids "$rundir"/*.pid
}
