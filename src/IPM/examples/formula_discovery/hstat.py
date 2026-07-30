#!/usr/bin/env python
"""Friedman's H-statistic on the depth-4 forest — asks the winning model which pairs it exploits
(sharper than FAST's data-level nominations). H2_jk = var(PD_jk - PD_j - PD_k)/var(PD_jk)."""
import os, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import partial_dependence
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
D["u"]=L("L_alpha")+L("L_sigma2min"); D["khat"]=L("L_hdiag_max")-L("L_sigma2min")
FP=["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
SLATE=[c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")+tuple(FP)]
samp=D.sample(8000,random_state=0)   # PD is an average; sample is fine and fast
Xs=samp[SLATE].to_numpy(float)
PAIRS={"d_hi":[("L_alpha","L_hdiag_min"),("L_hdiag_max","L_sigma2min"),("L_force_tol","L_sigma2min"),
               ("L_r0_c","L_Ldiag_max"),("L_alpha","L_Ldiag_min"),("u","L_Ldiag_min"),("L_force_tol","L_r0_c")],
       "d_lo":[("L_hdiag_min","L_hdiag_max"),("L_r0_p","L_sigma2min"),("L_alpha","L_sigma2min"),
               ("u","L_alpha"),("L_force_tol","L_floor_tol"),("L_hdiag_max","is_hsd")]}
def H2(m,jf,kf):
    j=SLATE.index(jf); k=SLATE.index(kf)
    pj=partial_dependence(m,Xs,[j],grid_resolution=12,kind="average")["average"][0]
    pk=partial_dependence(m,Xs,[k],grid_resolution=12,kind="average")["average"][0]
    pjk=partial_dependence(m,Xs,[(j,k)],grid_resolution=12,kind="average")["average"][0]
    pj=pj-pj.mean(); pk=pk-pk.mean(); pjk=pjk-pjk.mean()
    inter=pjk-pj[:,None]-pk[None,:]
    return (inter**2).sum()/max((pjk**2).sum(),1e-12)
for tgt in ("d_hi","d_lo"):
    y=D[tgt].to_numpy(float)
    m=HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0).fit(D[SLATE].to_numpy(float),y)
    print(f"\n=== {tgt}: Friedman H² (fraction of pair-PD variance that is interaction) ===",flush=True)
    for a,b in sorted(PAIRS[tgt],key=lambda p:-H2(m,*p)):
        print(f"   H²={H2(m,a,b):.3f}   {a} × {b}",flush=True)
print("DONE",flush=True)
