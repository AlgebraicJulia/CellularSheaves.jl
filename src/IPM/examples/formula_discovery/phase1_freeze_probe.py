#!/usr/bin/env python
"""Freeze probe on canonical folds: after the nested run showed only {la,mg} are 5/5-stable and the
greedy L0 soup loses to the champion, evaluate FIXED supports honestly (structure fixed -> only coefs
refit per fold, no support leakage). Report all-families and ex-X04 side by side, per-fold, for:
  - champion (7-term off-the-dome)
  - {la, mg} only  (the 5/5 stable core)
  - stable-core (la,mg + the 3/5 terms at size 12)
  - forest(5) reference"""
import os, re, json, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in I["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam])
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
RAW={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
     "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
MED={f:float(np.median(RAW[f])) for f in RAW}; C={f:RAW[f]-MED[f] for f in RAW}  # centered
y=I["d_lo"].to_numpy(float); exX04=~(fam=="X04")
la,cd,hx,bm,mg=C["la"],C["cd"],C["hx"],C["bm"],C["mg"]
def h(f,rk): return np.maximum(C[f]-(rk-MED[f]),0.0)   # relu(raw f - rk)
TERMS={ "la":la,"mg":mg,"cd*hx*mg":cd*hx*mg,"hx*bm":hx*bm,"la*cd*bm":la*cd*bm,
        "relu(cd-10.52)*bm":h("cd",-10.52)*bm,"relu(cd-10.52)*la":h("cd",-10.52)*la }
def cvcols(cols):
    p=np.full(len(y),np.nan)
    for j in range(CF["k"]):
        tr=foldid!=j; p[foldid==j]=LinearRegression().fit(cols[tr],y[tr]).predict(cols[foldid==j])
    return p
def report(name,cols):
    p=cvcols(cols); allm=mean_absolute_error(y,p); exm=mean_absolute_error(y[exX04],p[exX04])
    perf=[round(float(mean_absolute_error(y[foldid==j],p[foldid==j])),2) for j in range(CF["k"])]
    print(f"  {name:24} all={allm:.3f}  ex-X04={exm:.3f}   per-fold {perf}")
    return allm,exm
print("=== FREEZE PROBE (canonical folds; fixed structure, coefs refit per fold) ===")
for j in range(CF["k"]): print(f"  fold {j} test families: {CF['fold_families'][str(j)]}")
# champion
CX=np.column_stack([RAW['la'],RAW['mg'],np.maximum(RAW['bm']-2.3,0),np.maximum(RAW['hx']-2.8,0),
     np.maximum(RAW['cd']+3.6,0),np.maximum(RAW['cd']+3.6,0)*RAW['mg'],np.maximum(RAW['cd']+10.7,0)*np.maximum(RAW['mg']+1,0)])
report("champion (7)",CX)
report("{la,mg} core (2)",np.column_stack([TERMS['la'],TERMS['mg']]))
report("stable-core (7: 5/5+3/5)",np.column_stack([TERMS[t] for t in TERMS]))
# forest5 for reference (all + ex-X04)
Xf=np.column_stack([RAW['la'],RAW['cd'],RAW['hx'],RAW['bm'],RAW['mg']])
pf=np.full(len(y),np.nan)
for j in range(CF["k"]):
    tr=foldid!=j; pf[foldid==j]=HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0).fit(Xf[tr],y[tr]).predict(Xf[foldid==j])
print(f"  {'forest(5) ref':24} all={mean_absolute_error(y,pf):.3f}  ex-X04={mean_absolute_error(y[exX04],pf[exX04]):.3f}")
print("DONE")
