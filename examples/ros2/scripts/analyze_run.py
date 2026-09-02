#!/usr/bin/env python3
"""Score a sheaf-swarm run from its state CSV, using the scenario's own structure.

Reads the same CSV that server.jl writes (--log / SHEAF_LOG) and that
scripts/fake_sim.py writes (--csv), so a Gazebo run and a pure-kinematics run are
scored by identical code; --baseline prints them side by side and the gap between
the columns is the physics.

The metrics come from swarm.json rather than being hardcoded, because "error"
means different things for different robots here. A quadrotor pinned through
PROJECT_Z is only ever asked to match the target's altitude; scoring its
horizontal distance to the target measures a constraint nobody imposed and
reports a perfectly behaved swarm as broken. So each pinning edge is scored
through its own restriction map, and each consensus edge through the formation
offset it actually encodes.

    analyze_run.py RUN.csv [--baseline FAKE.csv] [--config swarm.json] [--tail 0.2]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

MAPS = {"xy": (0, 1), "z": (2,), "identity": (0, 1, 2)}


def load(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"{path}: empty")
    return rows


def score(rows, cfg, tail_frac):
    agents = cfg["agents"]
    targets = cfg.get("targets", [])
    names = [a["name"] for a in agents]
    off = {a["name"]: a["formation"] for a in agents}
    tail = rows[int(len(rows) * (1.0 - tail_frac)):]

    def p(row, n):
        return [float(row[f"{n}_{c}"]) for c in "xyz"]

    have_targets = bool(targets) and f"{targets[0]['name']}_x" in rows[0]

    out = {"cons": sum(float(r["cons_residual"]) for r in tail) / len(tail),
           "track": sum(float(r["track_residual"]) for r in tail) / len(tail),
           "edge": {}, "pin": {}, "minsep": (math.inf, "", "")}

    for i, j in cfg.get("edges", []):
        a, b = names[i], names[j]
        errs = []
        for r in tail:
            pa, pb = p(r, a), p(r, b)
            d = [(pa[k] - off[a][k]) - (pb[k] - off[b][k]) for k in range(3)]
            errs.append(math.sqrt(sum(c * c for c in d)))
        out["edge"][f"{a}-{b}"] = sum(errs) / len(errs)

    if have_targets:
        for pin in cfg.get("pinning", []):
            a = names[pin["agent"]]
            tgt = targets[pin["target"]]["name"]
            comps = MAPS[pin["map"]]
            errs = []
            for r in tail:
                pa, pt = p(r, a), p(r, tgt)
                d = [(pa[k] - off[a][k]) - pt[k] for k in comps]
                errs.append(math.sqrt(sum(c * c for c in d)))
            out["pin"][f"{a}~{pin['map']}"] = sum(errs) / len(errs)

    for r in rows:
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                pa, pb = p(r, a), p(r, b)
                # (0,0,0) is the log's "no telemetry yet" placeholder.
                if pa == [0.0, 0.0, 0.0] or pb == [0.0, 0.0, 0.0]:
                    continue
                d = math.dist(pa, pb)
                if d < out["minsep"][0]:
                    out["minsep"] = (d, a, b)
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("run_csv")
    ap.add_argument("--baseline", default="")
    ap.add_argument("--config", default=os.path.join(here, "..", "config", "swarm.json"))
    ap.add_argument("--tail", type=float, default=0.2)
    ap.add_argument("--max-cons", type=float, default=1.0)
    ap.add_argument("--max-track", type=float, default=1.0)
    ap.add_argument("--max-edge", type=float, default=0.6)
    ap.add_argument("--max-pin", type=float, default=0.6)
    args = ap.parse_args()

    cfg = json.load(open(args.config))
    runs = [("run", score(load(args.run_csv), cfg, args.tail))]
    if args.baseline:
        runs.append(("baseline", score(load(args.baseline), cfg, args.tail)))

    w = 14
    hdr = "metric".ljust(30) + "".join(n.rjust(w) for n, _ in runs)
    print(hdr)
    print("-" * len(hdr))

    def row(label, key, sub=None):
        vals = []
        for _, sc in runs:
            v = sc[key] if sub is None else sc[key].get(sub, float("nan"))
            vals.append(f"{v:.3f}".rjust(w))
        print(label.ljust(30) + "".join(vals))

    row("consensus residual", "cons")
    row("tracking residual", "track")
    for k in sorted(runs[0][1]["edge"]):
        row(f"formation {k} [m]", "edge", k)
    for k in sorted(runs[0][1]["pin"]):
        row(f"track {k} [m]", "pin", k)
    for name, sc in runs:
        d, a, b = sc["minsep"]
        print(f"min separation ({name}): {d:.3f} m  ({a}/{b})")

    sc0 = runs[0][1]
    fails = []
    if sc0["cons"] >= args.max_cons:
        fails.append(f"consensus {sc0['cons']:.3f} >= {args.max_cons}")
    if sc0["track"] >= args.max_track:
        fails.append(f"tracking {sc0['track']:.3f} >= {args.max_track}")
    for k, v in sorted(sc0["edge"].items()):
        if v >= args.max_edge:
            fails.append(f"formation {k} {v:.3f} >= {args.max_edge}")
    for k, v in sorted(sc0["pin"].items()):
        if v >= args.max_pin:
            fails.append(f"track {k} {v:.3f} >= {args.max_pin}")
    print("RESULT", "ok" if not fails else "FAIL (" + "; ".join(fails) + ")")
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
