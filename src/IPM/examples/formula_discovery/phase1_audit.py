#!/usr/bin/env python
"""Phase-1 audits for the candidate IPM/floor formula (MAE 1.27, forest(5)=1.07).
Reproduce family-CV MAE, then run the mandatory battery:
  (1) error by zone and by log-alpha bin   (flat or explain)
  (2) derivative identity  d(d_lo)/d(log alpha) -> should be ~+1 decade/decade
  (3) per-family error table + all-families vs ex-X04
  (4) drop-test on every term (esp. the two small cd-interaction terms)
  (5) coefficient ranges across the 5 reporting folds
Family = problem.split('_')[0]. Selection never touches the reporting folds (structure is fixed)."""
import os, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
I=D[D["is_hsd"]==0].copy()
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
la=Lc("L_alpha")
mg=np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))
cd=Lc("L_c_res0_dual"); hx=Lc("L_hdiag_max"); bm=Lc("L_bar_hdiag_med")
y=I["d_lo"].to_numpy(float); zone=I["zone"].to_numpy()
fam=I["problem"].str.split("_").str[0].to_numpy()

# ---- candidate basis (the 7 non-constant terms, in formula order) ----
def basis(la,mg,cd,hx,bm):
    return np.column_stack([
        la,                                   # T1  log alpha
        mg,                                   # T2  margin
        np.maximum(bm-2.3,0.0),               # T3  (bm-2.3)+
        np.maximum(hx-2.8,0.0),               # T4  (hx-2.8)+
        np.maximum(cd+3.6,0.0),               # T5  (cd+3.6)+
        np.maximum(cd+3.6,0.0)*mg,            # T6  (cd+3.6)+ * margin
        np.maximum(cd+10.7,0.0)*np.maximum(mg+1.0,0.0),  # T7 (cd+10.7)+*(mg+1)+
    ])
NAMES=["log a","margin","(bm-2.3)+","(hx-2.8)+","(cd+3.6)+","(cd+3.6)+*mg","(cd+10.7)+*(mg+1)+"]
CAND=np.array([0.49,-0.57,-0.77,0.54,0.16,-0.03,0.10]); CAND0=-1.50
X=basis(la,mg,cd,hx,bm)

def famcv(Xb,y,refit=True):
    gkf=GroupKFold(5); p=np.full(len(y),np.nan); coefs=[]
    for a,b in gkf.split(Xb,y,fam):
        m=LinearRegression().fit(Xb[a],y[a]); p[b]=m.predict(Xb[b]); coefs.append((m.intercept_,*m.coef_))
    return p,np.array(coefs)

print(f"IPM/floor rows={len(y)}  families={len(np.unique(fam))}  problems={I['problem'].nunique()}")
# literal formula (fixed published coefs, no refit)
plit=X@CAND+CAND0
print(f"\n[repro] literal formula (published coefs, no refit): MAE={mean_absolute_error(y,plit):.3f}")
p,coefs=famcv(X,y)
print(f"[repro] family-5fold CV (structure fixed, per-fold OLS refit): MAE={mean_absolute_error(y,p):.3f}   (brief: 1.27)")

# ---- AUDIT 1: error by zone, by log-alpha bin ----
print("\n=== AUDIT 1a: |error| by zone ===")
for z in ("below","in","above"):
    mk=zone==z
    if mk.any(): print(f"   {z:6}: MAE={mean_absolute_error(y[mk],p[mk]):.3f}  bias={np.mean(p[mk]-y[mk]):+.3f}  n={mk.sum()}")
print("=== AUDIT 1b: |error| by log-alpha bin ===")
qs=np.quantile(la,[0,.2,.4,.6,.8,1.0])
for i in range(5):
    lo,hi=qs[i],qs[i+1]; mk=(la>=lo)&(la<=hi if i==4 else la<hi)
    if mk.any(): print(f"   log a in [{lo:5.1f},{hi:5.1f}): MAE={mean_absolute_error(y[mk],p[mk]):.3f}  bias={np.mean(p[mk]-y[mk]):+.3f}  n={mk.sum()}")

# ---- AUDIT 2: derivative identity d(d_lo)/d(log alpha) ----
# within-iteration the oracle sweeps alpha on deepcopies (state=0): r0 co-moves, so measure
# empirical d(margin)/d(log a) per (problem,iter) group, then total slope = 0.49 - 0.57*dmargin/dloga.
print("\n=== AUDIT 2: derivative identity  d(d_lo)/d(log alpha)  (expect ~ +1) ===")
grp=I.groupby(["problem","iter"]) if "iter" in I.columns else None
if grp is not None:
    slopes=[]
    for _,idx in grp.groups.items():
        sub=I.loc[idx]; a=Lc("L_alpha")[I.index.get_indexer(idx)]
        if len(np.unique(a))<3: continue
        m2=np.maximum(sub["L_r0_c"].fillna(sub["L_r0_c"].median()).to_numpy(float),
                      sub["L_r0_p"].fillna(sub["L_r0_p"].median()).to_numpy(float))-\
           np.maximum(sub["L_force_tol"].fillna(0).to_numpy(float),sub["L_floor_tol"].fillna(0).to_numpy(float))
        A=np.polyfit(a,m2,1)[0]; slopes.append(A)
    slopes=np.array(slopes)
    dmarg=np.median(slopes)
    print(f"   empirical d(margin)/d(log a) [median over iters] = {dmarg:+.3f}   (floor law r0~1/a -> expect ~ -1)")
    print(f"   formula total slope  d(d_lo)/d(log a) = +0.49 + (-0.57)*({dmarg:+.2f}) = {0.49-0.57*dmarg:+.3f}   (expect ~ +1)")
else:
    print("   [no per-iter alpha sweep column found -- skipping]")

# ---- AUDIT 3: per-family error, all vs ex-X04 ----
print("\n=== AUDIT 3: per-family MAE ===")
for f in sorted(np.unique(fam)):
    mk=fam==f
    print(f"   {f:12}: MAE={mean_absolute_error(y[mk],p[mk]):.3f}  n={mk.sum():5d}")
exX04=~np.char.startswith(fam.astype(str),"X04")
print(f"\n   ALL families : MAE={mean_absolute_error(y,p):.3f}")
print(f"   ex-X04       : MAE={mean_absolute_error(y[exX04],p[exX04]):.3f}")
print(f"   X04-only     : MAE={mean_absolute_error(y[~exX04],p[~exX04]):.3f}  n={(~exX04).sum()}")

# ---- AUDIT 4: drop-test every term ----
print("\n=== AUDIT 4: drop-test (family-CV MAE with each term removed; delta vs full) ===")
base,_=famcv(X,y); baseMAE=mean_absolute_error(y,base)
for j in range(X.shape[1]):
    keep=[k for k in range(X.shape[1]) if k!=j]
    pj,_=famcv(X[:,keep],y); mj=mean_absolute_error(y,pj)
    print(f"   drop {NAMES[j]:20}: MAE={mj:.3f}  delta={mj-baseMAE:+.3f}")

# ---- AUDIT 5: coefficient ranges across folds ----
print("\n=== AUDIT 5: coefficient ranges across the 5 folds (const + 7 terms) ===")
lbl=["const"]+NAMES
for j in range(coefs.shape[1]):
    c=coefs[:,j]; print(f"   {lbl[j]:20}: [{c.min():+.3f}, {c.max():+.3f}]  mean={c.mean():+.3f}  published={([CAND0]+list(CAND))[j]:+.3f}")
print("\nDONE")
