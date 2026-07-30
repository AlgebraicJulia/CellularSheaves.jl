#!/usr/bin/env python
"""Product-SINDy — brute force in a good basis. Stage-C library + ALL pairwise products of the
continuous survivors+composites + physics hinge-products (threshold × conditioning) + counter indicators.
Lasso path -> sparsity target -> OLS refit -> grouped 5-fold CV. Prints the frontier + best formula."""
import os, itertools, numpy as np, pandas as pd
from sklearn.linear_model import lasso_path, LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
groups=D["problem"].values; gkf=GroupKFold(5); zone=D["zone"].to_numpy()
Lc=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
# continuous base = survivors + composites
BASE={n:Lc(n) for n in ["L_alpha","L_r0_c","L_r0_p","L_p_res0_dual","L_c_res0_dual","L_force_tol",
      "L_floor_tol","L_mu","L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_sigma2min"]}
BASE["u"]=Lc("L_alpha")+Lc("L_sigma2min"); BASE["margin"]=Lc("L_r0_c")-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))
BASE["logN"]=Lc("L_r0_c")+Lc("L_alpha"); BASE["tolratio"]=Lc("L_force_tol")-Lc("L_floor_tol")
BASE["khat"]=Lc("L_hdiag_max")-Lc("L_sigma2min")
def library():
    T=dict(BASE)
    for a,b in itertools.combinations(BASE,2): T[f"{a}*{b}"]=BASE[a]*BASE[b]          # all pairwise products
    for thr,cond in [("u","L_Ldiag_min"),("u","L_hdiag_min"),("margin","L_alpha"),
                     ("u","L_alpha"),("margin","L_Ldiag_min")]:                        # hinge-products: threshold x conditioning
        T[f"({thr})+*{cond}"]=np.maximum(BASE[thr],0.0)*BASE[cond]
    for c in ["pbase","cbase","ncraig"]:                                              # counter indicators
        x=D[c].fillna(0).to_numpy(float)
        for kk in (2,3,5,8,12,17): T[f"[{c}>={kk}]"]=(x>=kk).astype(float)
    return pd.DataFrame(T)
Xdf=library(); names=list(Xdf.columns); X=Xdf.to_numpy(float)
Xs=(X-X.mean(0))/(X.std(0)+1e-12)
def cvmae(supp,y):
    Xr=X[:,supp]; p=np.full(len(y),np.nan)
    for a,b in gkf.split(Xr,y,groups): p[b]=LinearRegression().fit(Xr[a],y[a]).predict(Xr[b])
    return mean_absolute_error(y,p),{z:mean_absolute_error(y[zone==z],p[zone==z]) for z in ("in","above")}
print(f"library: {len(names)} terms (base {len(BASE)} + products)",flush=True)
for tgt in ("d_lo","d_hi"):
    y=D[tgt].to_numpy(float)
    alphas,coefs,_=lasso_path(Xs,y,n_alphas=120,max_iter=30000)
    supps=(np.abs(coefs)>1e-8).sum(0)
    print(f"\n=== {tgt}: product-SINDy frontier ===",flush=True)
    best=None
    for nt in (8,12,16,20):
        j=int(np.argmin(np.abs(supps-nt))); supp=[i for i in range(len(names)) if abs(coefs[i,j])>1e-8]
        mae,z=cvmae(supp,y)
        if best is None or mae<best[0]: best=(mae,supp)
        print(f"  n={len(supp):2d}: CV-MAE={mae:.3f}  [in {z['in']:.2f} above {z['above']:.2f}]",flush=True)
    mf=LinearRegression().fit(X[:,best[1]],y)
    terms=sorted(zip([names[i] for i in best[1]],mf.coef_),key=lambda kv:-abs(kv[1]))
    print(f"  BEST {tgt} MAE={best[0]:.3f} ≈ {mf.intercept_:+.2f}",flush=True)
    for n,c in terms[:12]: print(f"      {c:+.3f} · {n}",flush=True)
print("DONE",flush=True)
