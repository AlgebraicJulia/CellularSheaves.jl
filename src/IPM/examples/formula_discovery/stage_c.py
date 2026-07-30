#!/usr/bin/env python
"""Stage C — SINDy-style sparse formula (§4C).
Term library: curated survivors + physics composites (margin, logN, tol-ratio, ceiling-margin) +
hinges at quantile knots + counter indicators/ramps + a few pairwise products. Lasso (grouped-CV
alpha) -> OLS refit on the selected support -> grouped-CV MAE. Prints the readable formula."""
import os, numpy as np, pandas as pd
from sklearn.linear_model import lasso_path, LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error

D = pd.read_pickle(os.path.join(os.path.dirname(__file__), "dataset.pkl")).copy()
groups = D["problem"].values
gkf = GroupKFold(n_splits=5)

# curated survivors (Stage A): residual reps, tolerances, structure, counters
SURV = ["L_r0_c","L_r0_p","L_p_res0_dual","L_c_res0_dual","L_force_tol","L_floor_tol","L_mu","L_alpha",
        "L_Ldiag_min","L_Ldiag_max","L_hdiag_min","L_hdiag_max","L_bar_hdiag_med","L_rho",
        "pbase","cbase","ncraig","cpass","ppass","stall","is_hsd"]
CONT = [c for c in SURV if c.startswith("L_")]

def library(df):
    T = {}
    for c in SURV: T[c] = df[c].fillna(df[c].median()).to_numpy(float)
    L = lambda c: df[c].fillna(df[c].median()).to_numpy(float)
    # physics composites
    tolmax = np.maximum(L("L_force_tol"), L("L_floor_tol"))
    T["margin=r0-maxtol"]      = L("L_r0_c") - tolmax
    T["logN=r0+alpha"]         = L("L_r0_c") + L("L_alpha")
    T["tolratio=force-floor"]  = L("L_force_tol") - L("L_floor_tol")
    T["ceilmargin=Ldmin+alpha"]= L("L_Ldiag_min") + L("L_alpha")     # ~ log(alpha*lambda_min(F)) ceiling scale
    T["kappa=Ldmax-Ldmin"]     = L("L_Ldiag_max") - L("L_Ldiag_min") # log kappa(F) proxy
    T["r0margin_ceil=r0-Ldmin"]= L("L_r0_c") - L("L_Ldiag_min")
    # hinges at quantile knots for the strongest continuous terms
    for c in ["L_r0_c","L_force_tol","L_mu","L_Ldiag_min","L_hdiag_min","margin=r0-maxtol"]:
        x = T[c] if c in T else L(c)
        for q in (0.25,0.5,0.75):
            k = np.quantile(x, q); T[f"({c}-{k:.1f})+"] = np.maximum(x-k,0.0)
    # hinges at EBM-indicated bends (Stage B shapes)
    EBMK = {"L_force_tol":[-5.0], "L_Ldiag_min":[-2.0], "L_hdiag_min":[-2.0], "L_p_res0_dual":[-5.0]}
    for c,ks in EBMK.items():
        x = T[c] if c in T else L(c)
        for k in ks: T[f"({c}-{k:.1f})+ebm"] = np.maximum(x-k,0.0)
    # counter indicators + ramps (extended to cover the pbase staircase saturating ~17)
    for c in ["pbase","cbase","ncraig"]:
        x = df[c].fillna(0).to_numpy(float)
        for kk in (2,3,5,8,12,17): T[f"[{c}>={kk}]"] = (x>=kk).astype(float)
        T[f"({c}-2)+"] = np.maximum(x-2,0.0)
    return pd.DataFrame(T)

def cvmae(Xr, y):
    pred = np.full(len(y), np.nan)
    for tr,te in gkf.split(Xr,y,groups):
        pred[te] = LinearRegression().fit(Xr[tr],y[tr]).predict(Xr[te])
    z = {zz: mean_absolute_error(y[D.zone==zz], pred[D.zone==zz]) for zz in ("below","in","above")}
    return mean_absolute_error(y,pred), z

def fit(target, nterms=8):
    Xdf = library(D); names = list(Xdf.columns)
    X = Xdf.to_numpy(float); y = D[target].to_numpy(float)
    Xs = (X - X.mean(0)) / (X.std(0)+1e-12)
    # walk the lasso path; pick the alpha whose support size is closest to nterms
    alphas, coefs, _ = lasso_path(Xs, y, n_alphas=100, max_iter=20000)
    supps = (np.abs(coefs) > 1e-8).sum(0)                    # support size per alpha
    j = int(np.argmin(np.abs(supps - nterms)))
    supp = [i for i in range(len(names)) if abs(coefs[i,j]) > 1e-8]
    Xr = X[:, supp]
    mae, z = cvmae(Xr, y)
    mf = LinearRegression().fit(Xr, y)                        # full-fit coefs for display
    terms = sorted(zip([names[i] for i in supp], mf.coef_), key=lambda kv: -abs(kv[1]))
    return mae, z, mf.intercept_, terms, len(supp)

if __name__=="__main__":
    for target,ops in (("d_lo",("below","in")),("d_hi",("in","above"))):
        mae,z,b0,terms,ns=fit(target)
        opm=np.nan
        print(f"\n===== {target}: sparse formula ({ns} terms), grouped-CV MAE={mae:.3f}  "
              f"[below={z['below']:.2f} in={z['in']:.2f} above={z['above']:.2f}] =====")
        print(f"  {target} ≈ {b0:+.2f}")
        for name,co in terms[:14]:
            print(f"      {co:+.3f} · {name}")
