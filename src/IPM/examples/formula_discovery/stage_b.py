#!/usr/bin/env python
"""Stage B — capacity ladder + EBM shape functions (§4B).
Ladder: depth-4 vs depth-2 vs depth-1(additive) vs linear -> interaction vs 1-D content.
EBM (GAM): grouped-CV MAE + read per-feature shapes (counter staircases, continuous bends)."""
import os, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor

D = pd.read_pickle(os.path.join(os.path.dirname(__file__), "dataset.pkl"))
FINGERPRINT = ["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
CORE = [c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")+tuple(FINGERPRINT)]
groups = D["problem"].values
gkf = GroupKFold(n_splits=5)

def cv(model_fn, target, feats):
    X = D[feats].to_numpy(float); y = D[target].to_numpy(float); pred=np.full(len(y),np.nan)
    for tr,te in gkf.split(X,y,groups):
        m=model_fn(); m.fit(X[tr],y[tr]); pred[te]=m.predict(X[te])
    return mean_absolute_error(y,pred)

def ladder(target):
    hgb=lambda d:(lambda: HistGradientBoostingRegressor(max_depth=d,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0))
    lin=lambda: make_pipeline(SimpleImputer(),Ridge(alpha=1.0))
    print(f"  {target}: depth4={cv(hgb(4),target,CORE):.3f}  depth2={cv(hgb(2),target,CORE):.3f} "
          f" depth1(add)={cv(hgb(1),target,CORE):.3f}  linear={cv(lin,target,CORE):.3f}")

def ebm_fit(target, nsub=20000):
    idx=np.random.default_rng(0).choice(len(D), min(nsub,len(D)), replace=False)
    Ds=D.iloc[idx]
    X=np.nan_to_num(Ds[CORE].to_numpy(float), nan=-1000.0); y=Ds[target].to_numpy(float); g=Ds["problem"].values
    tr,te=next(GroupKFold(5).split(X,y,g))                    # one grouped split for a quick MAE
    m=ExplainableBoostingRegressor(interactions=0,random_state=0); m.fit(X[tr],y[tr])
    mae=mean_absolute_error(y[te],m.predict(X[te]))
    m=ExplainableBoostingRegressor(interactions=0,random_state=0); m.fit(X,y)  # full-subsample fit for shapes
    return mae, m, m.explain_global(), sorted(zip(CORE,m.term_importances()),key=lambda kv:-kv[1])

def _isnum(x):
    try: float(x); return True
    except: return False

if __name__=="__main__":
    import sys
    if "--ladder" in sys.argv:
        print("=== capacity ladder (grouped-CV MAE, core feats) ===")
        ladder("d_lo"); ladder("d_hi")
    for target in ("d_lo","d_hi"):
        mae,m,g,imp=ebm_fit(target)
        print(f"\n=== EBM ({target}): grouped-CV MAE={mae:.3f} ===  top terms by importance:")
        for name,im in imp[:8]: print(f"    {im:.3f}  {name}")
        # print shapes for a few key features
        keys = ["L_Ldiag_min","L_force_tol","L_mu","L_r0_c","pbase","ncraig","L_hdiag_min"]
        names = list(m.term_names_)
        for k in keys:
            if k not in names: continue
            i=names.index(k); d=g.data(i)
            xs=d.get('names'); ys=d.get('scores')
            if xs is None or ys is None: continue
            xs=list(xs); ys=list(ys)
            n=min(len(xs),len(ys));
            if n==0: continue
            sel=np.linspace(0,n-1,min(7,n)).astype(int)
            def fx(v): return f"{float(v):.2g}" if _isnum(v) else str(v)
            pts=", ".join(f"{fx(xs[j])}:{ys[j]:+.2f}" for j in sel)
            print(f"    shape {k}: {pts}")
