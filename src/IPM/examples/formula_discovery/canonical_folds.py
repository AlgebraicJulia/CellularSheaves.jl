#!/usr/bin/env python
"""Canonicalize the evaluation folds as a committed ARTIFACT (not a recipe).
Canonical family rule: token = problem.split('_')[0]; if it matches X<digits><suffix>, strip the
suffix -> X-base (X09bt->X09, X04b->X04). All other tokens (e01, 06, 13, SOS) kept verbatim.
This collapses dial/twin variants so a family never straddles folds via its twin. Expect 21 families.
Writes canonical_folds.json {rule, family->fold map, per-fold family lists}. Verifies the handshake:
forest(5) in-sample ~0.46, family-CV ~1.07 on these folds."""
import os, re, json, numpy as np, pandas as pd, sklearn
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
def canon_family(p):
    tok=str(p).split("_")[0]; m=re.match(r"^(X\d+)",tok)
    return m.group(1) if m else tok
fam=np.array([canon_family(p) for p in I["problem"]])
fams=sorted(np.unique(fam)); print(f"canonical families: {len(fams)}\n  {fams}")

# balanced deterministic fold assignment (largest-family-first -> least-loaded fold). No RNG => reproducible.
K=5; sizes={f:int((fam==f).sum()) for f in fams}; load=[0]*K; assign={}
for f in sorted(fams,key=lambda f:(-sizes[f],f)):
    j=int(np.argmin(load)); assign[f]=j; load[j]+=sizes[f]
folds=[(np.array([assign[x] for x in fam])!=j, np.array([assign[x] for x in fam])==j) for j in range(K)]
print(f"per-fold row loads: {load}")

Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
FEAT={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
      "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
X5=np.column_stack([FEAT[f] for f in ["la","cd","hx","bm","mg"]]); y=I["d_lo"].to_numpy(float)
mk=lambda: HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0)
insample=mean_absolute_error(y,mk().fit(X5,y).predict(X5))
p=np.full(len(y),np.nan)
for tr,te in folds: p[te]=mk().fit(X5[tr],y[tr]).predict(X5[te])
famcv=mean_absolute_error(y,p)
print(f"\n=== HANDSHAKE (forest5 on canonical folds) ===")
print(f"  in-sample MAE = {insample:.3f}   (expect ~0.46)")
print(f"  family-CV MAE = {famcv:.3f}   (expect ~1.07)")

out={"rule":"family=problem.split('_')[0]; if re.match('^(X\\\\d+)') strip suffix to X-base; else verbatim",
     "n_families":len(fams),"k":K,"assignment":"balanced largest-first least-loaded, no RNG (deterministic)",
     "family_to_fold":{f:int(assign[f]) for f in fams},
     "fold_families":{str(j):[f for f in fams if assign[f]==j] for j in range(K)},
     "handshake":{"forest5_insample":round(float(insample),4),"forest5_familycv":round(float(famcv),4),
                  "sklearn":sklearn.__version__}}
json.dump(out,open(os.path.join(HERE,"canonical_folds.json"),"w"),indent=1)
print("\nwrote canonical_folds.json  DONE")
