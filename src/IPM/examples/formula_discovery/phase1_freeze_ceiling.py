#!/usr/bin/env python
"""Audit battery + drop-test for the FROZEN IPM ceiling (d_hi) stable-core (11 terms, >=3/5).
Zone / log-alpha flatness, per-family table, drop-test each term (to see if 11 trims down),
all vs ex-X04. Saves phase1_frozen_ceiling.json."""
import os, re, json, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in I["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam]); y=I["d_hi"].to_numpy(float); zone=I["zone"].to_numpy()
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
RAWc={"Ld":Lc("L_Ldiag_min"),"cd":Lc("L_c_res0_dual"),"ft":Lc("L_force_tol"),"pd":Lc("L_p_res0_dual"),"la":Lc("L_alpha")}
MED={f:float(np.median(RAWc[f])) for f in RAWc}; C={f:RAWc[f]-MED[f] for f in RAWc}
def hi(f,rk): return np.maximum(C[f]-(rk-MED[f]),0.0)
Ld,cd,ft,pd,la=C["Ld"],C["cd"],C["ft"],C["pd"],C["la"]; la_raw=RAWc["la"]; exX04=~(fam=="X04")
TERMS={"relu(ft+9.86)":hi("ft",-9.86),"Ld":Ld,"relu(pd+10.21)":hi("pd",-10.21),"pd":pd,
       "relu(Ld+4.47)":hi("Ld",-4.47),"la":la,"Ld*la":Ld*la,"Ld*ft":Ld*ft,"ft*pd":ft*pd,
       "relu(ft+4.96)":hi("ft",-4.96),"ft*la":ft*la}
NAMES=list(TERMS); X=np.column_stack([TERMS[t] for t in NAMES])
def cvpred(cols):
    p=np.full(len(y),np.nan)
    for j in range(CF["k"]):
        tr=foldid!=j; p[foldid==j]=LinearRegression().fit(cols[tr],y[tr]).predict(cols[foldid==j])
    return p
p=cvpred(X); base=mean_absolute_error(y,p); full=LinearRegression().fit(X,y)
print(f"=== FROZEN ceiling (11 terms) canonical CV: all={base:.3f}  ex-X04={mean_absolute_error(y[exX04],p[exX04]):.3f}  (forest 1.080) ===")
print("=== AUDIT 1a: MAE by zone ===")
for z in ("below","in","above"):
    m=zone==z
    if m.any(): print(f"   {z:6}: MAE={mean_absolute_error(y[m],p[m]):.3f}  bias={np.mean(p[m]-y[m]):+.3f}  n={m.sum()}")
print("=== AUDIT 1b: MAE by log-alpha bin ===")
qs=np.quantile(la_raw,[0,.2,.4,.6,.8,1.0])
for i in range(5):
    lo,ha=qs[i],qs[i+1]; m=(la_raw>=lo)&(la_raw<=ha if i==4 else la_raw<ha)
    if m.any(): print(f"   logα∈[{lo:4.1f},{ha:4.1f}): MAE={mean_absolute_error(y[m],p[m]):.3f}  bias={np.mean(p[m]-y[m]):+.3f}  n={m.sum()}")
print("=== AUDIT 2: per-family MAE ===")
famtab={}
for f in sorted(np.unique(fam)):
    m=fam==f; famtab[f]=round(float(mean_absolute_error(y[m],p[m])),3); print(f"   {f:6}: MAE={famtab[f]:.3f}  n={m.sum()}")
print("=== AUDIT 3: drop-test (canonical-CV MAE with each term removed) ===")
drop={}
for i in range(len(NAMES)):
    keep=[k for k in range(len(NAMES)) if k!=i]; mj=mean_absolute_error(y,cvpred(X[:,keep])); drop[NAMES[i]]=round(mj-base,3)
    print(f"   drop {NAMES[i]:20}: MAE={mj:.3f}  Δ={mj-base:+.3f}")
frozen={"target":"d_hi","solver":"ipm","cv_all":round(base,3),"cv_exX04":round(float(mean_absolute_error(y[exX04],p[exX04])),3),
        "forest_ref":1.080,"intercept":round(float(full.intercept_),4),
        "terms":{NAMES[i]:{"coef":round(float(full.coef_[i]),4),"drop_delta":drop[NAMES[i]]} for i in range(len(NAMES))},
        "medians":{k:round(v,3) for k,v in MED.items()},"folds":"canonical_folds.json",
        "audits":{"by_zone":{z:round(float(mean_absolute_error(y[zone==z],p[zone==z])),3) for z in ("below","in","above") if (zone==z).any()},"per_family":famtab}}
json.dump(frozen,open(os.path.join(HERE,"phase1_frozen_ceiling.json"),"w"),indent=1)
print("\nsaved phase1_frozen_ceiling.json  DONE")
