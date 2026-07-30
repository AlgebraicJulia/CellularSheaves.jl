#!/usr/bin/env python
"""Fresh RFE on the PURE-IPM dataset (is_hsd dropped, IPM rows only), from the full feature slate down.
Backward elimination by permutation importance; track grouped holdout MAE; report the surviving set.
Question: do we recover the same essential five, or were they an artifact of the mixed-data derivation?
(The 5 were distilled from 8, which were generated on the MIXED ipm+hsd dataset.)"""
import numpy as np, pandas as pd, os
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
D=D[D["is_hsd"]==0].copy().reset_index(drop=True)          # PURE IPM
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
# composites (same defs as the mixed-data RFE, so the comparison is apples-to-apples)
D["u"]=L("L_alpha")+L("L_sigma2min"); D["margin"]=np.maximum(L("L_r0_c"),L("L_r0_p"))-np.maximum(L("L_force_tol"),L("L_floor_tol"))
D["tolratio"]=L("L_force_tol")-L("L_floor_tol"); D["khat"]=L("L_hdiag_max")-L("L_sigma2min")
CONT=["L_alpha","L_r0_c","L_r0_p","L_p_res0_dual","L_c_res0_dual","L_force_tol","L_floor_tol","L_mu",
      "L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_sigma2min","L_rho",
      "u","margin","tolratio","khat"]
CNT=["pbase","cbase","ncraig","cpass","ppass","stall"]     # is_hsd DROPPED
FULL=CONT+CNT
for c in FULL: D[c]=D[c].fillna(D[c].median() if c in CONT else 0)
# fixed grouped 80/20 by problem (same protocol that produced the original 8->5)
rng=np.random.default_rng(0); probs=D["problem"].unique().copy(); rng.shuffle(probs)
trp=set(probs[:int(0.8*len(probs))]); trm=D["problem"].isin(trp).to_numpy()
y=D["d_lo"].to_numpy(float); ytr,yte=y[trm],y[~trm]
print(f"PURE-IPM floor RFE: {len(D)} rows, {D['problem'].nunique()} problems, {len(FULL)} start features (is_hsd dropped)",flush=True)
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
print("\n=== PURE-IPM floor RFE curve (n features -> holdout MAE) ===",flush=True)
for n,mae,_ in curve: print(f"   n={n:2d}: MAE={mae:.3f}{'  <-- knee (+0.03)' if n==knee[0] else ''}",flush=True)
print(f"\nKNEE set ({knee[0]} feats, MAE={knee[1]:.3f}, full={full_mae:.3f}):\n   {knee[2]}",flush=True)
# what survives at exactly 5, for direct comparison to the mixed-derived five
five=[c for c in curve if c[0]==5]
if five: print(f"\nAt n=5 (MAE={five[0][1]:.3f}):  {five[0][2]}",flush=True)
ORIG5={"L_alpha","L_c_res0_dual","L_hdiag_max","L_bar_hdiag_med","margin"}
for n,mae,fs in curve:
    if n<=8:
        ov=ORIG5 & set(fs); print(f"   overlap with mixed-derived 5 at n={n}: {len(ov)}/5  {sorted(ov)}",flush=True)
print("DONE",flush=True)
