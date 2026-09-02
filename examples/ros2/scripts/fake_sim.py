#!/usr/bin/env python3
"""ROS-free harness for the sheaf server.

Emulates the swarm with the same kinematics the bridge assumes -- velocity
integrators for the quadrotors, an offset-point unicycle for the ground robot --
and drives it through the real ZMQ protocol against a real `server.jl`. That
isolates the control law and the wire format from Gazebo, ROS 2, PX4 and MAVROS,
so a protocol or gain regression fails here in seconds instead of inside a
simulator.

    julia --project=examples/ros2 examples/ros2/server.jl &
    python3 examples/ros2/scripts/fake_sim.py --duration 60 --csv /tmp/run.csv
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import zmq

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "ros2_ws", "src", "sheaf_ros2")
)

from sheaf_ros2.protocol import AgentState, format_state, parse_control  # noqa: E402
from sheaf_ros2.transforms import world_velocity_to_unicycle  # noqa: E402


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(here, "..", "config", "swarm.json"))
    ap.add_argument("--duration", type=float, default=60.0, help="simulated seconds")
    ap.add_argument("--csv", default="", help="optional trajectory log")
    ap.add_argument("--timeout-ms", type=int, default=5000)
    # Hold an agent's `ok` flag low for the first N simulated seconds, e.g.
    # --not-ready 3:12. Off the default path on purpose: it exists to exercise
    # the server's synchronized-start gate and its present-but-not-active
    # obstacle rows, which need one agent to join late while the rest are up.
    ap.add_argument("--not-ready", action="append", default=[],
                    metavar="AGENT:SECONDS",
                    help="report ok=False for agent AGENT until SECONDS")
    args = ap.parse_args()

    not_ready = {}
    for spec in args.not_ready:
        who, _, secs = spec.partition(":")
        not_ready[int(who)] = float(secs)

    with open(args.config) as fh:
        cfg = json.load(fh)

    dt = 1.0 / float(cfg["rate_hz"])
    agents = cfg["agents"]
    n = len(agents)

    pos = [list(map(float, a["spawn"])) for a in agents]
    for i, a in enumerate(agents):
        if a["kind"] == "px4":
            pos[i][2] = float(a["z_ref"])  # assume the takeoff already happened
    vel = [[0.0, 0.0, 0.0] for _ in range(n)]
    yaw = [0.0] * n

    ctx = zmq.Context()
    sock = ctx.socket(zmq.REQ)
    sock.setsockopt(zmq.RCVTIMEO, args.timeout_ms)
    sock.setsockopt(zmq.LINGER, 0)
    sock.connect(cfg["zmq_endpoint"])

    log = open(args.csv, "w") if args.csv else None
    tnames = [t["name"] for t in cfg.get("targets", [])]
    if log:
        cols = (["t"]
                + [f"{a['name']}_{c}" for a in agents for c in ("x", "y", "z")]
                + [f"{n}_{c}" for n in tnames for c in ("x", "y", "z")])
        log.write(",".join(cols + ["cons_residual", "track_residual"]) + "\n")

    steps = int(args.duration / dt)
    t = 0.0
    # Peak, not first: agents join the active sheaf only after the server's
    # join debounce, so the first ticks report a vacuous 0.0 residual.
    peak_cons = peak_track = 0.0
    cons = track = float("nan")
    # Residual history over the last quarter of the run, to tell "settled at a
    # nonzero value" apart from "still moving".
    tail_cons, tail_track = [], []
    tail_from = int(steps * 0.75)
    # Safety-filter activity, straight from the server's diagnostics: how often
    # the CBF filter had to touch a command, how often it had to fall back to
    # raw repulsion, and the closest approach it saw in its own pair metrics
    # (air-ground pairs are horizontal distances, air-air are in R^3).
    cbf_ticks = cbf_peak = cbf_fallback = 0
    cbf_minsep = float("inf")
    cbf_zone_peak = cbf_obs_peak = 0
    cbf_zone_margin = float("inf")
    # When the server's synchronized-start gate opened, and who it was last
    # waiting on, so a run can show that nothing translated before it did.
    started_at = None
    waiting_last = []
    target_trace = []

    for step in range(steps):
        state = [
            AgentState(
                id=i,
                kind=agents[i]["kind"],
                p=pos[i],
                v=vel[i],
                yaw=yaw[i],
                ok=(t >= not_ready.get(i, 0.0)),
            )
            for i in range(n)
        ]
        sock.send(format_state(t, state))
        out = parse_control(sock.recv())

        cons = out.diag.get("consensus_residual", float("nan"))
        track = out.diag.get("tracking_residual", float("nan"))
        peak_cons = max(peak_cons, cons)
        peak_track = max(peak_track, track)
        n_cbf = int(out.diag.get("cbf_active", 0))
        cbf_ticks += 1 if n_cbf else 0
        cbf_peak = max(cbf_peak, n_cbf)
        cbf_fallback += int(out.diag.get("cbf_fallback", 0))
        sep = float(out.diag.get("min_separation", -1.0))
        if sep >= 0.0:
            cbf_minsep = min(cbf_minsep, sep)
        n_zone = int(out.diag.get("cbf_zone_constraints", 0))
        cbf_zone_peak = max(cbf_zone_peak, n_zone)
        cbf_obs_peak = max(cbf_obs_peak, int(out.diag.get("cbf_obstacle_rows", 0)))
        if n_zone:
            cbf_zone_margin = min(cbf_zone_margin,
                                  float(out.diag.get("min_zone_margin", -1.0)))
        if out.diag.get("started", True):
            if started_at is None:
                started_at = t
        else:
            waiting_last = list(out.diag.get("waiting_for", []))
        if step >= tail_from:
            tail_cons.append(cons)
            tail_track.append(track)
        # the targets' own motion, for the over-constrained test after the run
        if step % 5 == 0:
            target_trace.append([list(tp.p) for tp in out.targets])

        for c in out.agents:
            i = c.id
            spec = agents[i]
            if spec["kind"] == "diffdrive":
                v, w = world_velocity_to_unicycle(
                    c.u[0], c.u[1], yaw[i],
                    float(spec["offset"]), float(spec["v_max"]), float(spec["omega_max"]),
                )
                yaw[i] += w * dt
                pos[i][0] += v * math.cos(yaw[i]) * dt
                pos[i][1] += v * math.sin(yaw[i]) * dt
                vel[i] = [v * math.cos(yaw[i]), v * math.sin(yaw[i]), 0.0]
            else:
                vel[i] = list(c.u)
                for k in range(3):
                    pos[i][k] += c.u[k] * dt

        if log:
            row = ([f"{t:.3f}"]
                   + [f"{c:.4f}" for p in pos for c in p]
                   + [f"{c:.4f}" for tp in out.targets for c in tp.p])
            log.write(",".join(row + [f"{cons:.6f}", f"{track:.6f}"]) + "\n")
        t += dt

    if log:
        log.close()

    sock.send(format_state(t, [], stop=True))
    try:
        sock.recv()
    except zmq.Again:
        pass

    print(f"steps            {steps}  (dt={dt:.3f}s, {t:.1f}s simulated)")
    print(f"consensus resid  peak {peak_cons:.4f} -> {cons:.4f}")
    print(f"tracking resid   peak {peak_track:.4f} -> {track:.4f}")
    sep_txt = "n/a" if cbf_minsep == float("inf") else f"{cbf_minsep:.3f} m"
    print(f"cbf filter       active on {cbf_ticks}/{steps} ticks "
          f"(peak {cbf_peak} constraints, {cbf_fallback} fallbacks), "
          f"closest pair {sep_txt}")
    if cbf_zone_peak or cbf_obs_peak:
        margin_txt = ("n/a" if cbf_zone_margin == float("inf")
                      else f"{cbf_zone_margin:+.3f} m")
        print(f"cbf geometry     peak {cbf_zone_peak} keep-out rows "
              f"(worst zone margin {margin_txt}), "
              f"peak {cbf_obs_peak} rows against present-but-inactive agents")
    if started_at is None:
        print(f"start gate       never opened, still waiting for "
              f"{', '.join(waiting_last) or 'nobody'}")
    else:
        print(f"start gate       opened at t = {started_at:.2f}s"
              + (f" (last waiting for {', '.join(waiting_last)})"
                 if waiting_last else ""))
    for i, a in enumerate(agents):
        p = pos[i]
        print(f"  {a['name']:<8} {a['kind']:<10} p = [{p[0]:7.3f} {p[1]:7.3f} {p[2]:7.3f}]")

    ok = cons < peak_cons and cons < 1.0 and track < peak_track and track < 1.0
    if ok:
        print("RESULT converging")
        return 0

    # A scenario with no global section settles at the least-squares compromise
    # instead of at zero: the sheaf giving the right answer, not the controller
    # failing. Reported as settled rather than as a failure only when the config
    # actually asks for it AND the residual has stopped moving; anything else
    # nonzero is still a failure, so a real regression cannot hide here.
    #
    # Is the scenario over-constrained? Not a question about which fields are
    # present: a pin's own offset is a CONSTANT, so it can absorb a constant
    # difference between what two agents are asked to track and nothing more.
    # Two agents joined by a consensus edge can hold formation only if the
    # targets they follow MOVE the same way. So compare the motion directly,
    # with each trajectory's own mean removed, and call it a conflict when the
    # time-varying parts disagree.
    pinned = {}
    for p in cfg.get("pinning", []):
        pinned.setdefault(p["agent"], []).append(p["target"])

    def centred_mean_motion(ks):
        # mean position of targets `ks` over time, with its own average removed
        series = [[sum(fr[k][d] for k in ks) / len(ks) for d in range(3)]
                  for fr in target_trace]
        avg = [sum(f[d] for f in series) / len(series) for d in range(3)]
        return [[f[d] - avg[d] for d in range(3)] for f in series]

    over_constrained = False
    if target_trace:
        for i, j in cfg.get("edges", []):
            if not pinned.get(i) or not pinned.get(j):
                continue
            mi, mj = centred_mean_motion(pinned[i]), centred_mean_motion(pinned[j])
            rms = math.sqrt(sum((a[d] - b[d]) ** 2 for a, b in zip(mi, mj)
                                for d in range(3)) / len(mi))
            if rms > 0.05:
                over_constrained = True
                break
    settled = (len(tail_cons) > 10
               and max(tail_cons) - min(tail_cons) < 0.01 * max(1.0, cons)
               and max(tail_track) - min(tail_track) < 0.01 * max(1.0, track))
    if over_constrained and settled:
        print(f"RESULT settled at a nonzero residual (over-constrained by design: "
              f"consensus edges join agents that track different targets)")
        return 0
    print("RESULT NOT converging" + (" (residual still moving)" if not settled else ""))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
