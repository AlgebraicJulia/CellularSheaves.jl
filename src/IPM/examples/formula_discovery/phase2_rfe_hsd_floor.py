#!/usr/bin/env python
"""RFE redone HONESTLY: canonical 21-family 5-fold CV (canonical_folds.json), not the leaky problem-level
80/20. IPM rows only, is_hsd dropped, full 25-feature slate. Backward elimination by permutation
importance averaged over the 5 canonical folds. Question: do the same five survive on honest folds?
Compare to the leaky-protocol five {L_alpha, L_c_res0_dual, L_hdiag_max, L_bar_hdiag_med, margin}."""
import os, re, json, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); D=D[D["is_hsd"]==1].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in D["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam])
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
D["u"]=L("L_alpha")+L("L_sigma2min"); D["margin"]=np.maximum(L("L_r0_c"),L("L_r0_p"))-np.maximum(L("L_force_tol"),L("L_floor_tol"))
D["tolratio"]=L("L_force_tol")-L("L_floor_tol"); D["khat"]=L("L_hdiag_max")-L("L_sigma2min")
CONT=["L_alpha","L_r0_c","L_r0_p","L_p_res0_dual","L_c_res0_dual","L_force_tol","L_floor_tol","L_mu",
      "L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_sigma2min","L_rho",
      "u","margin","tolratio","khat"]
CNT=["pbase","cbase","ncraig","cpass","ppass","stall"]; FULL=CONT+CNT
for c in FULL: D[c]=D[c].fillna(D[c].median() if c in CONT else 0)
y=D["d_lo"].to_numpy(float)
mk=lambda: HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0)
def cv_and_importance(feats):
    Xall=D[feats].to_numpy(float); p=np.full(len(y),np.nan); imp=np.zeros(len(feats))
    for j in range(CF["k"]):
        tr=foldid!=j; te=foldid==j; m=mk().fit(Xall[tr],y[tr]); p[te]=m.predict(Xall[te])
        r=permutation_importance(m,Xall[te],y[te],n_repeats=2,random_state=0,scoring="neg_mean_absolute_error")
        imp+=r.importances_mean*te.sum()
    return mean_absolute_error(y,p),imp/len(y)
print(f"CANONICAL RFE: {len(D)} IPM rows, {CF['n_families']} families, {len(FULL)} start features",flush=True)
feats=list(FULL); curve=[]
while len(feats)>=3:
    mae,imp=cv_and_importance(feats); curve.append((len(feats),mae,list(feats)))
    print(f"   n={len(feats):2d}: canon-CV MAE={mae:.3f}",flush=True)
    k=3 if len(feats)>15 else 1; drop=[feats[i] for i in np.argsort(imp)[:k]]; feats=[f for f in feats if f not in drop]
full_mae=curve[0][1]; knee=min((c for c in curve if c[1]<=full_mae+0.03),key=lambda c:c[0])
print(f"\nKNEE ({knee[0]} feats, MAE={knee[1]:.3f}, full={full_mae:.3f}):\n   {knee[2]}",flush=True)
ORIG5={"L_alpha","L_c_res0_dual","L_hdiag_max","L_bar_hdiag_med","margin"}
print("\noverlap with leaky-protocol five, by size:",flush=True)
for n,mae,fs in curve:
    if n<=9: ov=ORIG5&set(fs); print(f"   n={n}: {len(ov)}/5 present  {sorted(ov)}   MAE={mae:.3f}",flush=True)
json.dump([(n,round(m,4),fs) for n,m,fs in curve],open(os.path.join(HERE,"phase2_rfe_hsd_floor.json"),"w"),indent=1)
print("DONE",flush=True)
