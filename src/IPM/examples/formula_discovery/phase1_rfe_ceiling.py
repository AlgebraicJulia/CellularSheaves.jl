#!/usr/bin/env python
"""IPM CEILING (d_hi) — honest RFE on canonical folds, full 25-feature slate, IPM rows only.
Same protocol as the floor RFE (phase1_rfe_canonical.py) but target = d_hi. The ceiling is expected to
be the conditioning / lambda_min-shrink story (Ldiag_min, hdiag_min, dual residual) rather than the
floor's margin law -- let the RFE find its essential set fresh."""
import os, re, json, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); D=D[D["is_hsd"]==0].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in D["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam])
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
D["u"]=L("L_alpha")+L("L_sigma2min"); D["margin"]=np.maximum(L("L_r0_c"),L("L_r0_p"))-np.maximum(L("L_force_tol"),L("L_floor_tol"))
D["tolratio"]=L("L_force_tol")-L("L_floor_tol"); D["khat"]=L("L_hdiag_max")-L("L_sigma2min")
D["ceilmargin"]=L("L_Ldiag_min")+L("L_alpha")   # lambda_min(F)-side ceiling proxy from alpha-window theory
CONT=["L_alpha","L_r0_c","L_r0_p","L_p_res0_dual","L_c_res0_dual","L_force_tol","L_floor_tol","L_mu",
      "L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_sigma2min","L_rho",
      "u","margin","tolratio","khat","ceilmargin"]
CNT=["pbase","cbase","ncraig","cpass","ppass","stall"]; FULL=CONT+CNT
for c in FULL: D[c]=D[c].fillna(D[c].median() if c in CONT else 0)
y=D["d_hi"].to_numpy(float)
mk=lambda: HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0)
def cv_and_imp(feats):
    Xall=D[feats].to_numpy(float); p=np.full(len(y),np.nan); imp=np.zeros(len(feats))
    for j in range(CF["k"]):
        tr=foldid!=j; te=foldid==j; m=mk().fit(Xall[tr],y[tr]); p[te]=m.predict(Xall[te])
        r=permutation_importance(m,Xall[te],y[te],n_repeats=2,random_state=0,scoring="neg_mean_absolute_error")
        imp+=r.importances_mean*te.sum()
    return mean_absolute_error(y,p),imp/len(y)
print(f"CEILING RFE (d_hi): {len(D)} IPM rows, {CF['n_families']} families, {len(FULL)} start features",flush=True)
feats=list(FULL); curve=[]
while len(feats)>=3:
    mae,imp=cv_and_imp(feats); curve.append((len(feats),mae,list(feats)))
    print(f"   n={len(feats):2d}: canon-CV MAE={mae:.3f}   least-imp drop next: {[feats[i] for i in np.argsort(imp)[:(3 if len(feats)>15 else 1)]]}",flush=True)
    k=3 if len(feats)>15 else 1; feats=[f for i,f in enumerate(feats) if i not in np.argsort(imp)[:k]]
full_mae=curve[0][1]; knee=min((c for c in curve if c[1]<=full_mae+0.03),key=lambda c:c[0])
print(f"\nKNEE ({knee[0]} feats, MAE={knee[1]:.3f}, full={full_mae:.3f}):\n   {knee[2]}",flush=True)
for n,mae,fs in curve:
    if n<=10: print(f"   n={n}: MAE={mae:.3f}  {sorted(fs)}",flush=True)
json.dump([(n,round(m,4),fs) for n,m,fs in curve],open(os.path.join(HERE,"phase1_rfe_ceiling.json"),"w"),indent=1)
print("DONE",flush=True)
