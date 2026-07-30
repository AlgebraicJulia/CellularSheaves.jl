#!/usr/bin/env python
"""FREEZE the IPM/floor formula = the stability-selected 7-term core (5/5 + 3/5 terms from the nested
run), fixed structure. Fit coefficients on all IPM data; report per-fold coefficient ranges; run the
mandatory audit battery on canonical folds (zone / log-alpha flatness, per-family table, drop-test each
term). Save phase1_frozen_formula.json for the coworker."""
import os, re, json, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in I["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam])
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
RAW={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
     "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
MED={f:float(np.median(RAW[f])) for f in RAW}; C={f:RAW[f]-MED[f] for f in RAW}
y=I["d_lo"].to_numpy(float); zone=I["zone"].to_numpy(); la_raw=RAW["la"]
def h(f,rk): return np.maximum(C[f]-(rk-MED[f]),0.0)
la,cd,hx,bm,mg=C["la"],C["cd"],C["hx"],C["bm"],C["mg"]
TERMS={"la":la,"mg":mg,"cd*hx*mg":cd*hx*mg,"hx*bm":hx*bm,"la*cd*bm":la*cd*bm,
       "relu(cd-10.52)*bm":h("cd",-10.52)*bm,"relu(cd-10.52)*la":h("cd",-10.52)*la}
NAMES=list(TERMS); X=np.column_stack([TERMS[t] for t in NAMES])
def cvpred(cols):
    p=np.full(len(y),np.nan); coefs=[]
    for j in range(CF["k"]):
        tr=foldid!=j; m=LinearRegression().fit(cols[tr],y[tr]); p[foldid==j]=m.predict(cols[foldid==j]); coefs.append([m.intercept_]+list(m.coef_))
    return p,np.array(coefs)
p,coefs=cvpred(X); full=LinearRegression().fit(X,y)
exX04=~(fam=="X04")
print("=== FROZEN IPM/floor formula (7 terms, centered; medians below) ===")
print(f"  medians for centering: {({k:round(v,2) for k,v in MED.items()})}")
print(f"  canonical-fold CV MAE: all={mean_absolute_error(y,p):.3f}  ex-X04={mean_absolute_error(y[exX04],p[exX04]):.3f}   (forest5 1.109/0.809, champion 1.320/1.048)")
print(f"\n  d_lo = {full.intercept_:+.3f}")
order=np.argsort(-np.abs(full.coef_))
frozen={"target":"d_lo","solver":"ipm","terms":[],"intercept_alldata":round(float(full.intercept_),4),
        "medians":{k:round(v,4) for k,v in MED.items()},"folds":"canonical_folds.json",
        "cv_mae_all":round(float(mean_absolute_error(y,p)),4),"cv_mae_exX04":round(float(mean_absolute_error(y[exX04],p[exX04])),4)}
for i in order:
    cr=coefs[:,i+1]; act=(TERMS[NAMES[i]]!=0).mean()
    print(f"       {full.coef_[i]:+.3f} · {NAMES[i]:20}  [fold {cr.min():+.3f},{cr.max():+.3f}]  active {act*100:.0f}%")
    frozen["terms"].append({"term":NAMES[i],"coef_alldata":round(float(full.coef_[i]),4),
        "fold_range":[round(float(cr.min()),4),round(float(cr.max()),4)],"activity":round(float(act),3)})
print("  (centered: la=logα−9.0, cd=log c_res0_dual+7.27, hx=log hdiag_max−1.63, bm=log bar_hdiag_med+0.87, mg=margin+3.86)")

print("\n=== AUDIT 1a: MAE by zone ===")
for z in ("below","in","above"):
    m=zone==z
    if m.any(): print(f"   {z:6}: MAE={mean_absolute_error(y[m],p[m]):.3f}  bias={np.mean(p[m]-y[m]):+.3f}  n={m.sum()}")
print("=== AUDIT 1b: MAE by log-alpha bin ===")
qs=np.quantile(la_raw,[0,.2,.4,.6,.8,1.0])
for i in range(5):
    lo,hi=qs[i],qs[i+1]; m=(la_raw>=lo)&(la_raw<=hi if i==4 else la_raw<hi)
    if m.any(): print(f"   logα∈[{lo:4.1f},{hi:4.1f}): MAE={mean_absolute_error(y[m],p[m]):.3f}  bias={np.mean(p[m]-y[m]):+.3f}  n={m.sum()}")
print("=== AUDIT 2: per-family MAE ===")
famtab={}
for f in sorted(np.unique(fam)):
    m=fam==f; e=float(mean_absolute_error(y[m],p[m])); famtab[f]=round(e,3)
    print(f"   {f:6}: MAE={e:.3f}  n={m.sum()}")
print("=== AUDIT 3: drop-test (canonical-CV MAE with each term removed) ===")
base=mean_absolute_error(y,p); drop={}
for i in range(len(NAMES)):
    keep=[k for k in range(len(NAMES)) if k!=i]; pj,_=cvpred(X[:,keep]); mj=mean_absolute_error(y,pj)
    drop[NAMES[i]]=round(mj-base,3); print(f"   drop {NAMES[i]:20}: MAE={mj:.3f}  Δ={mj-base:+.3f}")
frozen["audits"]={"by_zone":{z:round(float(mean_absolute_error(y[zone==z],p[zone==z])),3) for z in ("below","in","above") if (zone==z).any()},
                  "per_family":famtab,"drop_test_delta":drop,"all_vs_exX04":[round(float(mean_absolute_error(y,p)),3),round(float(mean_absolute_error(y[exX04],p[exX04])),3)]}
json.dump(frozen,open(os.path.join(HERE,"phase1_frozen_formula.json"),"w"),indent=1)
print("\nsaved phase1_frozen_formula.json  DONE")
