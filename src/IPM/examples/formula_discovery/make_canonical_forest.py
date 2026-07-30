#!/usr/bin/env python
"""STEP 2: train & push THE canonical reference object -- ipm_floor_forest5.pkl -- from this pipeline,
scored on the committed canonical_folds.json. Full provenance travels with it. Deprecate all other
forest scaffolding so there is exactly one answer to 'which object'."""
import os, json, pickle, numpy as np, pandas as pd, sklearn
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CK=os.path.join(HERE,"ebm_ckpt")
CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
import re
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in I["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam])
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
FEATDEF={"la":"log10(alpha)","cd":"log10(c_res0_dual)","hx":"log10(hdiag_max)","bm":"log10(bar_hdiag_med)",
         "mg":"log10(max(r0_c,r0_p)) - log10(max(force_tol,floor_tol))"}
FEAT={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
      "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
ORDER=["la","cd","hx","bm","mg"]; X=np.column_stack([FEAT[f] for f in ORDER]); y=I["d_lo"].to_numpy(float)
HP=dict(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0)
mk=lambda: HistGradientBoostingRegressor(**HP)
insample=float(mean_absolute_error(y,mk().fit(X,y).predict(X)))
p=np.full(len(y),np.nan)
for j in range(CF["k"]):
    tr=foldid!=j; p[foldid==j]=mk().fit(X[tr],y[tr]).predict(X[foldid==j])
famcv=float(mean_absolute_error(y,p))
model=mk().fit(X,y)
prov={"model":model,"features":ORDER,"feature_recipe":FEATDEF,"target":"d_lo","solver":"ipm",
      "row_filter":"is_hsd==0 (IPM rows only)","hyperparams":HP,"seed":0,"sklearn":sklearn.__version__,
      "folds_artifact":"canonical_folds.json","folds_rule":CF["rule"],"n_families":CF["n_families"],
      "n_rows":int(len(y)),"handshake":{"insample_mae":round(insample,4),"canonical_familycv_mae":round(famcv,4)}}
out=os.path.join(CK,"ipm_floor_forest5.pkl"); pickle.dump(prov,open(out,"wb"))
print(f"WROTE canonical reference: {out}")
print(f"  features {ORDER}  recipe.margin = {FEATDEF['mg']}")
print(f"  hyperparams {HP}  sklearn {sklearn.__version__}")
print(f"  handshake: in-sample={insample:.4f}  canonical-family-CV={famcv:.4f}")
# deprecate the container scaffolding + my earlier wrong-folds pickle
DEP=["forest_d_lo.pkl","forest_d_hi.pkl","forest5_ipm_d_lo.pkl"]
notes=[]
for f in DEP:
    p2=os.path.join(CK,f)
    if os.path.exists(p2):
        os.rename(p2,p2+".DEPRECATED"); notes.append(f+" -> "+f+".DEPRECATED")
open(os.path.join(CK,"DEPRECATED.txt"),"w").write(
    "Deprecated forest scaffolding. THE canonical IPM-floor reference is ipm_floor_forest5.pkl\n"
    "scored on canonical_folds.json. Do not adopt the *.DEPRECATED pickles:\n  - forest_d_lo/hi.pkl: "
    "16-feature MIXED-data forests trained container-side (superseded)\n  - forest5_ipm_d_lo.pkl: "
    "earlier 5-feature forest on non-canonical (31-family) folds (superseded)\n")
print("  deprecated:",notes)
print("DONE")
