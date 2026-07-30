#!/usr/bin/env python
"""Stage D step 0 — pair nomination. EBM with interactions=12 (weak settings) on the Stage-A survivors;
ignore shapes, extract the FAST-ranked pairwise-interaction list = which products matter for PySR."""
import os, numpy as np, pandas as pd
from interpret.glassbox import ExplainableBoostingRegressor

HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl"))
SURV=["L_r0_c","L_r0_p","L_p_res0_dual","L_c_res0_dual","L_force_tol","L_floor_tol","L_mu","L_alpha",
      "L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_rho",
      "pbase","cbase","ncraig","cpass","ppass","stall","is_hsd","L_sigma2min"]
X=D[SURV].to_numpy(float)
for target in ("d_hi","d_lo"):
    m=ExplainableBoostingRegressor(interactions=12, outer_bags=4, max_rounds=3000, random_state=0)
    m.fit(X, D[target].to_numpy(float))
    imp=m.term_importances(); tf=m.term_features_
    pairs=[(imp[i], tf[i]) for i in range(len(tf)) if len(tf[i])==2]
    mains={SURV[tf[i][0]]: imp[i] for i in range(len(tf)) if len(tf[i])==1}
    pairs.sort(reverse=True)
    print(f"\n=== {target}: top FAST-ranked interaction PAIRS (importance) ===",flush=True)
    for im,(a,b) in pairs:
        print(f"   {im:.3f}   {SURV[a]} × {SURV[b]}",flush=True)
print("DONE",flush=True)
