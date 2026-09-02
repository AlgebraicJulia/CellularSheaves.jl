#!/usr/bin/env python3
"""Plot the safety filter against its own absence on the same scenario.

Reads two runs of config/collide.json, one with SHEAF_CBF=0 and one with it on,
and draws the distance between the two quadrotors over time. The point of the
figure is the floor: without the filter the pair passes through each other, with
it the distance flattens onto d_min and stays there.

    python3 scripts/plot_cbf.py off.csv on.csv docs/img/cbf-separation.png
"""
import csv, math, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def series(path):
    rows = list(csv.DictReader(open(path)))
    t = [float(r["t"]) for r in rows]
    d = [math.dist([float(r["uav0_x"]), float(r["uav0_y"]), float(r["uav0_z"])],
                   [float(r["uav1_x"]), float(r["uav1_y"]), float(r["uav1_z"])])
         for r in rows]
    return t, d

def main(off_csv, on_csv, out, d_min=1.0):
    to, do = series(off_csv)
    tn, dn = series(on_csv)
    fig, ax = plt.subplots(figsize=(9, 4.2), dpi=150)
    ax.axhspan(0, d_min, color="#b3261e", alpha=0.07, zorder=0)
    ax.axhline(d_min, color="#b3261e", lw=1.2, ls="--", zorder=2,
               label=f"$d_{{\\min}}$ = {d_min:.1f} m")
    ax.plot(to, do, color="#b3261e", lw=1.8, label=f"filter off (min {min(do):.3f} m)")
    ax.plot(tn, dn, color="#2456a6", lw=2.0, label=f"filter on (min {min(dn):.3f} m)")
    ax.set_xlabel("time [s]"); ax.set_ylabel("distance between uav0 and uav1 [m]")
    ax.set_title("Control barrier filter on a scenario built to collide", loc="left",
                 fontsize=11, fontweight="600")
    ax.set_ylim(0, max(max(do), max(dn)) * 1.05); ax.set_xlim(min(to), max(to))
    ax.grid(alpha=0.25, lw=0.6); ax.legend(loc="upper right", framealpha=0.95, fontsize=9)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out)
    print(f"{out}: off min {min(do):.4f} m, on min {min(dn):.4f} m")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
