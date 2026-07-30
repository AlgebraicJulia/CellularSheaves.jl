#!/usr/bin/env python
"""Checkpointed EBM fitter. Each fit is a job that pickles its model the instant it finishes.
Queue: the 2 FULL fits first (shape functions), then the 10 CV fits (5 folds x 2 targets, rigorous
grouped-CV MAE). Resumable: rerun skips any job whose pickle exists; halting keeps completed fits.
Full 76k data, NaN kept native (EBM handles it). Prints shapes as soon as the full fits land."""
import os, pickle, time, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor

HERE = os.path.dirname(__file__)
CKPT = os.path.join(HERE, "ebm_ckpt"); os.makedirs(CKPT, exist_ok=True)
D = pd.read_pickle(os.path.join(HERE, "dataset.pkl"))
FP = ["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
CORE = [c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")+tuple(FP)]
X_ALL = D[CORE].to_numpy(float)
groups = D["problem"].values
folds = list(GroupKFold(5).split(X_ALL, D["d_lo"].to_numpy(float), groups))

def job(tag, target, fold):
    """fold=None -> full fit (for shapes); fold=k -> train on fold-k's train, predict its test."""
    ckf = os.path.join(CKPT, tag + ".pkl")
    if os.path.exists(ckf):
        return pickle.load(open(ckf, "rb"))
    y = D[target].to_numpy(float)
    t0 = time.time()
    if fold is None:
        m = ExplainableBoostingRegressor(interactions=0, random_state=0).fit(X_ALL, y)
        rec = {"tag": tag, "target": target, "fold": None, "model": m, "secs": time.time()-t0}
    else:
        tr, te = folds[fold]
        m = ExplainableBoostingRegressor(interactions=0, random_state=0).fit(X_ALL[tr], y[tr])
        rec = {"tag": tag, "target": target, "fold": fold, "model": m,
               "te": te, "pred": m.predict(X_ALL[te]), "secs": time.time()-t0}
    pickle.dump(rec, open(ckf, "wb"))
    print(f"  [done {rec['secs']:.0f}s] {tag}", flush=True)
    return rec

def shapes(rec, keys=("L_Ldiag_min","L_force_tol","L_mu","L_r0_c","L_p_res0_dual","L_hdiag_min",
                      "pbase","cbase","ncraig")):
    m = rec["model"]; g = m.explain_global(); names = list(m.term_names_)
    imp = sorted(zip(CORE, m.term_importances()), key=lambda kv: -kv[1])
    print(f"\n=== EBM shapes ({rec['target']}) === top terms:")
    for n, i in imp[:8]: print(f"    {i:.3f}  {n}")
    for k in keys:
        if k not in names: continue
        d = g.data(names.index(k)); xs = list(d.get("names") or []); ys = list(d.get("scores") or [])
        n = min(len(xs), len(ys))
        if not n: continue
        sel = np.linspace(0, n-1, min(7, n)).astype(int)
        fx = lambda v: (f"{float(v):.2g}" if _isnum(v) else str(v))
        print(f"    shape {k}: " + ", ".join(f"{fx(xs[j])}:{ys[j]:+.2f}" for j in sel))

def _isnum(v):
    try: float(v); return True
    except: return False

if __name__ == "__main__":
    import sys
    QUEUE = [("full_d_lo","d_lo",None), ("full_d_hi","d_hi",None)]        # shapes first
    if "--cv" in sys.argv:
        for t in ("d_lo","d_hi"):
            for k in range(5): QUEUE.append((f"cv_{t}_f{k}", t, k))       # rigorous MAE second
    for tag, target, fold in QUEUE:
        rec = job(tag, target, fold)
        if fold is None: shapes(rec)
    # assemble any complete CV MAEs
    for t in ("d_lo","d_hi"):
        recs = [pickle.load(open(os.path.join(CKPT, f"cv_{t}_f{k}.pkl"),"rb"))
                for k in range(5) if os.path.exists(os.path.join(CKPT, f"cv_{t}_f{k}.pkl"))]
        if len(recs) == 5:
            pred = np.full(len(D), np.nan)
            for r in recs: pred[r["te"]] = r["pred"]
            print(f"\nEBM grouped-CV MAE {t} = {mean_absolute_error(D[t], pred):.3f}")
