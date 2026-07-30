#!/usr/bin/env python
"""The aiming measurement: EBM interactions=15 = the PAIRWISE CEILING (best any mains+2D model can do).
Near depth-4 (~0.55/0.62) => gap is pairwise, closable by enumeration. Stalls midway => higher-order."""
import os, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
FP=["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
SURV=[c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")+tuple(FP)]
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
COMP={"u=la+s2":L("L_alpha")+L("L_sigma2min"),"khat=hdmx-s2":L("L_hdiag_max")-L("L_sigma2min"),
 "s2-ftol":L("L_sigma2min")-L("L_force_tol"),"margin":L("L_r0_c")-np.maximum(L("L_force_tol"),L("L_floor_tol")),
 "logN":L("L_r0_c")+L("L_alpha"),"tolratio":L("L_force_tol")-L("L_floor_tol")}
for k,v in COMP.items(): D[k]=v
SLATE=SURV+list(COMP); groups=D["problem"].values; gkf=GroupKFold(3)
for tgt in ("d_lo","d_hi"):
    X=D[SLATE].to_numpy(float); y=D[tgt].to_numpy(float); p=np.full(len(y),np.nan)
    for a,b in gkf.split(X,y,groups):
        m=ExplainableBoostingRegressor(interactions=15,outer_bags=2,max_rounds=3000,random_state=0).fit(X[a],y[a]); p[b]=m.predict(X[b])
    print(f"  EBM interactions=15 {tgt}: 3-fold grouped-CV MAE = {mean_absolute_error(y,p):.3f}  (pairwise ceiling)",flush=True)
print("DONE",flush=True)
