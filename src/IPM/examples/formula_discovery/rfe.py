#!/usr/bin/env python
"""Recursive feature elimination on the depth-4 forest. Start with the full slate (survivors +
composites + counters + flags), repeatedly drop the least-important feature(s) by permutation importance,
track grouped-holdout MAE, and report the knee = smallest set within ~0.03 of the full-slate MAE."""
import numpy as np, pandas as pd, os
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
D["u"]=L("L_alpha")+L("L_sigma2min"); D["margin"]=L("L_r0_c")-np.maximum(L("L_force_tol"),L("L_floor_tol"))
D["tolratio"]=L("L_force_tol")-L("L_floor_tol"); D["khat"]=L("L_hdiag_max")-L("L_sigma2min")
CONT=["L_alpha","L_r0_c","L_r0_p","L_p_res0_dual","L_c_res0_dual","L_force_tol","L_floor_tol","L_mu",
      "L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_sigma2min","L_rho",
      "u","margin","tolratio","khat"]
CNT=["pbase","cbase","ncraig","cpass","ppass","stall","is_hsd"]
FULL=CONT+CNT
for c in FULL: D[c]=D[c].fillna(D[c].median() if c in CONT else 0)
# fixed grouped 80/20 split by problem
rng=np.random.default_rng(0); probs=D["problem"].unique().copy(); rng.shuffle(probs)
trp=set(probs[:int(0.8*len(probs))]); trm=D["problem"].isin(trp).to_numpy()
def fitmae(feats):
    Xtr=D.loc[trm,feats].to_numpy(float); ytr=D.loc[trm,"d_lo"].to_numpy(float)  # placeholder
    return None
def run(tgt):
    y=D[tgt].to_numpy(float); ytr,yte=y[trm],y[~trm]
    feats=list(FULL); curve=[]
    while len(feats)>=3:
        Xtr=D.loc[trm,feats].to_numpy(float); Xte=D.loc[~trm,feats].to_numpy(float)
        m=HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0).fit(Xtr,ytr)
        mae=mean_absolute_error(yte,m.predict(Xte))
        r=permutation_importance(m,Xte,yte,n_repeats=3,random_state=0,scoring="neg_mean_absolute_error")
        curve.append((len(feats),mae,list(feats)))
        k=3 if len(feats)>15 else 1
        drop=[feats[i] for i in np.argsort(r.importances_mean)[:k]]
        feats=[f for f in feats if f not in drop]
    full_mae=curve[0][1]; knee=min((c for c in curve if c[1]<=full_mae+0.03), key=lambda c:c[0])
    print(f"\n=== {tgt}: RFE curve (n features -> holdout MAE) ===",flush=True)
    for n,mae,_ in curve: print(f"   n={n:2d}: MAE={mae:.3f}{'  <-- knee' if n==knee[0] else ''}",flush=True)
    print(f"   ESSENTIAL SET ({knee[0]} feats, MAE={knee[1]:.3f}, full={full_mae:.3f}): {knee[2]}",flush=True)
for tgt in ("d_lo","d_hi"): run(tgt)
print("DONE",flush=True)
