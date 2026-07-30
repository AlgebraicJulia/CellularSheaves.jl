#!/usr/bin/env python
"""Build the §1-2 training matrix for the distance-sensor formula discovery.

One row per (iteration, candidate alpha). Labels d_lo/d_hi are signed decades to the
window edges (window = state==0 candidates with ncraig <= mincraig+1). Grouped by PROBLEM
(basename with solver/tol/dial suffixes stripped) so no instance splits across CV folds.
"""
import glob, os, re, json
import numpy as np, pandas as pd

ORACLE = os.path.join(os.path.dirname(__file__), "..", "oracle2")
OUT = os.path.join(os.path.dirname(__file__), "dataset.pkl")
LOG0 = 1e-300

# --- feature slate (§2) ---
LOG_FEATS = [  # positive-continuous -> log10(x + 1e-300)
    "r0_p","r0_c","r0_w","r1_p","r1_c","r1_w","pres0","pres1","cres0","cres1","wres0","wres1",
    "p_res0_dual","p_res0_prim","c_res0_dual","c_res0_prim","w_res0_dual","w_res0_prim",
    "p_res_exit","c_res_exit","w_res_exit","force_tol","floor_tol","mu","alpha","rho",
    "hdiag_min","hdiag_max","bar_hdiag_med","Ldiag_min","Ldiag_max","sigma2min",
]
RAW_FEATS = ["cbase","pbase","cpass","ppass","ncraig","wbase","wpass","bar_hdiag_frac_mid"]
META_LOG = ["alpha_anchor"]   # + m,n raw, ng -> log10(1+ng)

def group_of(base):
    g = base
    g = re.sub(r"_(ipm|hsd)$", "", g)
    g = re.sub(r"_dial[0-9.eE+-]+$", "", g)
    g = re.sub(r"_1e-?[0-9]+$", "", g)
    return g

def fnum(s):
    try: return float(s)
    except: return np.nan

def build():
    rows = []
    files = sorted(glob.glob(os.path.join(ORACLE, "*.csv")))
    for fn in files:
        base = os.path.basename(fn)[:-4]
        df = pd.read_csv(fn, low_memory=False)
        if df.empty or "state" not in df: continue
        metap = fn[:-4] + ".meta.json"
        meta = json.load(open(metap)) if os.path.exists(metap) else {}
        is_hsd = 1 if base.endswith("_hsd") else 0
        grp = group_of(base)
        # per-iteration window edges
        for it, sub in df.groupby("iter"):
            succ = sub[sub["state"] == 0]
            if succ.empty: continue
            mincr = succ["ncraig"].min()
            win = succ[succ["ncraig"] <= mincr + 1]
            lo = np.log10(win["alpha"].min()); hi = np.log10(win["alpha"].max())
            for _, r in sub.iterrows():
                la = np.log10(r["alpha"])
                d_lo = la - lo
                d_hi = hi - la
                zone = "below" if la < lo - 1e-9 else ("above" if la > hi + 1e-9 else "in")
                row = {"problem": grp, "file": base, "is_hsd": is_hsd,
                       "d_lo": d_lo, "d_hi": d_hi, "zone": zone}
                for c in LOG_FEATS:
                    v = fnum(r[c]) if c in sub.columns and pd.notna(r.get(c, np.nan)) else np.nan
                    row["L_" + c] = np.log10(v + LOG0) if (np.isfinite(v) and v >= 0) else np.nan
                for c in RAW_FEATS:
                    row[c] = fnum(r[c]) if c in sub.columns else np.nan
                # stall flag
                stall = 0
                for sc in ("pstat","cstat","wstat"):
                    if sc in sub.columns and str(r.get(sc,"")).strip() == "REFINE_STALLED": stall = 1
                row["stall"] = stall
                # meta features (audit-flagged)
                row["m_meta"] = meta.get("m", np.nan); row["n_meta"] = meta.get("n", np.nan)
                row["Lng_meta"] = np.log10(1 + meta.get("ng", 0.0))
                aa = meta.get("alpha_anchor", np.nan)
                row["L_alpha_anchor"] = np.log10(aa + LOG0) if (aa is not None and np.isfinite(aa) and aa > 0) else np.nan
                rows.append(row)
    D = pd.DataFrame(rows)
    D.to_pickle(OUT)
    return D

if __name__ == "__main__":
    D = build()
    print(f"rows={len(D):,}  problems={D['problem'].nunique()}  files={D['file'].nunique()}")
    print(f"feature cols={sum(c.startswith('L_') or c in RAW_FEATS or c in ('stall','is_hsd','m_meta','n_meta','Lng_meta') for c in D.columns)}")
    print("\nzone balance:"); print(D["zone"].value_counts())
    print("\nd_lo:  ", D["d_lo"].describe()[["mean","std","min","max"]].round(2).to_dict())
    print("d_hi:  ", D["d_hi"].describe()[["mean","std","min","max"]].round(2).to_dict())
    print(f"\nhsd rows={D['is_hsd'].sum():,}  ipm rows={(D['is_hsd']==0).sum():,}")
    print("top families by #problems:")
    fam = D.drop_duplicates("problem")["problem"].str.split("_").str[0].value_counts()
    print(fam.head(10).to_dict())
