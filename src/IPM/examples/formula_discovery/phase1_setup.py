#!/usr/bin/env python
"""Phase-1 setup: IPM/floor. Install FAMILY folds, calibrate references (forest5/forest46/lasso-v1),
extract knots from the fast-EBM shapes on the essential five. Family = problem.split('_')[0]."""
import os, itertools, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import lasso_path, LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
D["margin"]=L("L_r0_c")-np.maximum(L("L_force_tol"),L("L_floor_tol"))
I=D[D["is_hsd"]==0].copy()
FIVE=["L_alpha","L_c_res0_dual","L_hdiag_max","L_bar_hdiag_med","margin"]
SH={"L_alpha":"la","L_c_res0_dual":"cd","L_hdiag_max":"hx","L_bar_hdiag_med":"bm","margin":"mg"}
for c in FIVE: I[c]=I[c].fillna(I[c].median())
fam=I["problem"].str.split("_").str[0].to_numpy()
print(f"IPM/floor: {len(I)} rows, {len(np.unique(fam))} families, {I['problem'].nunique()} problems")
y=I["d_lo"].to_numpy(float); gkf=GroupKFold(5)
def fam_cv(X):
    p=np.full(len(y),np.nan)
    for a,b in gkf.split(X,y,fam):
        p[b]=HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0).fit(X[a],y[a]).predict(X[b])
    return mean_absolute_error(y,p)
FULL=[c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone","is_hsd","m_meta","n_meta","Lng_meta","L_alpha_anchor")]
for c in FULL:
    if c not in I: I[c]=D.loc[I.index,c]
    I[c]=I[c].fillna(I[c].median() if I[c].dtype!=object else 0)
print("=== family-fold references ===")
print(f"  forest(5)  = {fam_cv(I[FIVE].to_numpy(float)):.3f}   (brief: 1.07)")
print(f"  forest(46) = {fam_cv(I[FULL].to_numpy(float)):.3f}   (brief: 1.32)")
# lasso v1: order-<=3 basis, blind 33/66 knots
def order3_basis(knots):
    B={f:I[f].to_numpy(float) for f in FIVE}; T=dict(B)
    for f in FIVE:
        for k in knots[f]: T[f"({SH[f]}-{k:.1f})+"]=np.maximum(B[f]-k,0.0)
    for a,b in itertools.combinations(FIVE,2): T[f"{SH[a]}*{SH[b]}"]=B[a]*B[b]
    for a,b,c in itertools.combinations(FIVE,3): T[f"{SH[a]}*{SH[b]}*{SH[c]}"]=B[a]*B[b]*B[c]
    return pd.DataFrame(T)
qk={f:[np.quantile(I[f],q) for q in (0.33,0.66)] for f in FIVE}
Xdf=order3_basis(qk); X=Xdf.to_numpy(float); Xs=(X-X.mean(0))/(X.std(0)+1e-12)
al,co,_=lasso_path(Xs,y,alphas=120,max_iter=30000); sp=(np.abs(co)>1e-8).sum(0)
j=int(np.argmin(np.abs(sp-25))); supp=[i for i in range(X.shape[1]) if abs(co[i,j])>1e-8]
p=np.full(len(y),np.nan)
for a,b in gkf.split(X,y,fam): p[b]=LinearRegression().fit(X[a][:,supp] if False else X[np.ix_(a,supp)],y[a]).predict(X[np.ix_(b,supp)])
print(f"  lasso-v1 ({len(supp)}t, blind knots) = {mean_absolute_error(y,p):.3f}   (brief: 1.22)")
# EBM knots (bend points) on the five
m=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(I[FIVE].to_numpy(float),y)
print("\n=== EBM shape bends -> knots (per feature) ===")
for i,f in enumerate(FIVE):
    cuts=np.asarray(m.bins_[i][0]).ravel(); sc=np.asarray(m.term_scores_[i]).ravel(); body=sc[1:1+len(cuts)+1]
    # bends = where slope sign changes or biggest 2nd-difference
    if len(body)>3:
        d2=np.abs(np.diff(body,2)); bidx=np.argsort(-d2)[:2]+1
        kn=sorted(cuts[min(bi,len(cuts)-1)] for bi in bidx)
    else: kn=[np.quantile(I[f],0.5)]
    print(f"  {f:16} knots ~ {[round(float(k),1) for k in kn]}")
print("DONE")
