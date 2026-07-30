#!/usr/bin/env python
"""Pickle the 5-feature IPM-only floor forest (the forest(5) reference) so it can be reused without
retraining. Saves the full-data model + the exact family-fold CV MAE + feature construction recipe."""
import os, pickle, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CK=os.path.join(HERE,"ebm_ckpt")
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
FEAT={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
      "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
ORDER=["la","cd","hx","bm","mg"]
X=np.column_stack([FEAT[f] for f in ORDER]); y=I["d_lo"].to_numpy(float)
fam=I["problem"].str.split("_").str[0].to_numpy()
def family_folds(k,seed):
    fams=np.array(sorted(np.unique(fam))); rng=np.random.default_rng(seed); rng.shuffle(fams)
    sizes={f:int((fam==f).sum()) for f in fams}; load=[0]*k; assign={}
    for f in sorted(fams,key=lambda f:-sizes[f]):
        j=int(np.argmin(load)); assign[f]=j; load[j]+=sizes[f]
    a=np.array([assign[f] for f in fam]); return [(a!=j,a==j) for j in range(k)]
mk=lambda: HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0)
p=np.full(len(y),np.nan)
for tr,te in family_folds(5,0): p[te]=mk().fit(X[tr],y[tr]).predict(X[te])
cvmae=float(mean_absolute_error(y,p))
model=mk().fit(X,y)   # full-data model for reuse
RECIPE={"la":"log10(alpha)","cd":"log10(c_res0_dual)","hx":"log10(hdiag_max)","bm":"log10(bar_hdiag_med)",
        "mg":"log10(max(r0_c,r0_p)) - log10(max(force_tol,floor_tol))"}
out=os.path.join(CK,"forest5_ipm_d_lo.pkl")
pickle.dump({"model":model,"features":ORDER,"recipe":RECIPE,"target":"d_lo","solver":"ipm",
             "protocol":"family-5fold GroupKFold(seed0), depth4/it300/lr0.06/l2_1.0","family_cv_mae":cvmae,
             "n_rows":int(len(y)),"n_families":int(len(np.unique(fam)))},open(out,"wb"))
print(f"saved {out}")
print(f"  features (in order): {ORDER}")
print(f"  recipe: {RECIPE}")
print(f"  family-5fold CV MAE = {cvmae:.4f}  (this is the forest(5) reference)")
print("DONE")
