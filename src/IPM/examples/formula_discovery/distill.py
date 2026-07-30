#!/usr/bin/env python
"""Stage-D closeout — distill the additive-EBM (a GAM) into a closed-form formula. Since we proved the
structure is additive (no interactions), approximate each EBM curve with dense hinges (piecewise-linear)
+ composites + counter staircase. Lasso->sparsity target->OLS->grouped CV. Also pickles the depth-4
forests for reuse. Prints the two final formulas."""
import os, pickle, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import lasso_path, LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); CK=os.path.join(HERE,"ebm_ckpt"); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
groups=D["problem"].values; gkf=GroupKFold(5); zone=D["zone"].to_numpy()
Lc=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
# base features (survivors + composites)
B={n:Lc(n) for n in ["L_alpha","L_r0_c","L_p_res0_dual","L_force_tol","L_floor_tol","L_mu",
    "L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_sigma2min"]}
B["u"]=Lc("L_alpha")+Lc("L_sigma2min"); B["margin"]=Lc("L_r0_c")-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))
B["tolratio"]=Lc("L_force_tol")-Lc("L_floor_tol"); B["khat"]=Lc("L_hdiag_max")-Lc("L_sigma2min")
SLATE=list(B)

def library():
    T={}
    for n,x in B.items():
        T[n]=x                                                   # linear term
        for q in (0.2,0.4,0.6,0.8):                              # dense hinges -> piecewise-linear
            k=np.quantile(x,q); T[f"({n}-{k:.1f})+"]=np.maximum(x-k,0.0)
    for c in ["pbase","cbase","ncraig"]:                         # counter staircase
        v=D[c].fillna(0).to_numpy(float)
        for kk in (2,3,5,8,12,17): T[f"[{c}>={kk}]"]=(v>=kk).astype(float)
    return pd.DataFrame(T)
Xdf=library(); names=list(Xdf.columns); X=Xdf.to_numpy(float); Xs=(X-X.mean(0))/(X.std(0)+1e-12)
def cvmae(supp,y):
    Xr=X[:,supp]; p=np.full(len(y),np.nan)
    for a,b in gkf.split(Xr,y,groups): p[b]=LinearRegression().fit(Xr[a],y[a]).predict(Xr[b])
    return mean_absolute_error(y,p),{z:mean_absolute_error(y[zone==z],p[zone==z]) for z in ("in","above")}

print(f"distillation library: {len(names)} terms (additive: linear + dense hinges + counters)",flush=True)
for tgt in ("d_lo","d_hi"):
    y=D[tgt].to_numpy(float)
    alphas,coefs,_=lasso_path(Xs,y,alphas=140,max_iter=30000)
    supps=(np.abs(coefs)>1e-8).sum(0)
    print(f"\n=== {tgt}: distillation frontier ===",flush=True)
    picks={}
    for nt in (10,15,18,22,30,45):
        j=int(np.argmin(np.abs(supps-nt))); supp=[i for i in range(len(names)) if abs(coefs[i,j])>1e-8]
        mae,z=cvmae(supp,y); picks[len(supp)]=(mae,supp)
        print(f"  n={len(supp):2d}: CV-MAE={mae:.3f}  [in {z['in']:.2f} above {z['above']:.2f}]",flush=True)
    # final = the ~15-18 term pick with best MAE
    n15=min(picks, key=lambda k:(picks[k][0]))
    mae,supp=picks[n15]; mf=LinearRegression().fit(X[:,supp],y)
    terms=sorted(zip([names[i] for i in supp],mf.coef_),key=lambda kv:-abs(kv[1]))
    print(f"  FINAL {tgt} ({len(supp)} terms) MAE={mae:.3f}:  {tgt} ≈ {mf.intercept_:+.2f}",flush=True)
    for n,c in terms: print(f"      {c:+.3f} · {n}",flush=True)
    # save the depth-4 forest for coworker reuse
    fr=HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0).fit(np.column_stack([B[n] for n in SLATE]),y)
    pickle.dump({"model":fr,"features":SLATE,"target":tgt},open(os.path.join(CK,f"forest_{tgt}.pkl"),"wb"))
print("\n(saved forest_d_lo.pkl, forest_d_hi.pkl to ebm_ckpt/)  DONE",flush=True)
