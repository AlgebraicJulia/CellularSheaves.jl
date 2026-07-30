#!/usr/bin/env python
"""RuleFit — harvest the gradient-boosted forest's own decision rules (axis-aligned threshold
conjunctions) as binary features + linear terms, lasso-select, grouped-CV. Captures the threshold-
interaction structure trees use, which H²-per-pair can miss. If it beats the additive-composite EBM
(0.74/0.70), the forest's edge is threshold-rules; if not, the gap is irreducibly diffuse."""
import os, numpy as np, pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import lasso_path, LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
D["u"]=L("L_alpha")+L("L_sigma2min"); D["khat"]=L("L_hdiag_max")-L("L_sigma2min")
D["margin"]=L("L_r0_c")-np.maximum(L("L_force_tol"),L("L_floor_tol")); D["tolratio"]=L("L_force_tol")-L("L_floor_tol")
FP=["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
SLATE=[c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")+tuple(FP)]
Xslate=D[SLATE].fillna(D[SLATE].median()).to_numpy(float)
groups=D["problem"].values; gkf=GroupKFold(5); zone=D["zone"].to_numpy()

def extract_rules(gbr, X):
    """each internal node's root path = a rule; return dict name->binary column."""
    R={}
    for est in gbr.estimators_.ravel():
        t=est.tree_
        def rec(node, mask, cond):
            if t.children_left[node]==-1: return
            f=t.feature[node]; thr=t.threshold[node]
            lm=mask & (X[:,f]<=thr); rm=mask & (X[:,f]>thr)
            for m,c in ((lm,cond+[f"{SLATE[f]}<={thr:.2f}"]),(rm,cond+[f"{SLATE[f]}>{thr:.2f}"])):
                if len(c)>=2:            # only conjunctions (rules), skip single splits
                    key=" & ".join(sorted(c));
                    if key not in R and 0.02<m.mean()<0.98: R[key]=m.astype(float)
                rec(t.children_left[node] if m is lm else t.children_right[node], m, c)
        rec(0, np.ones(len(X),bool), [])
    return R

def cvmae(Xr,y):
    p=np.full(len(y),np.nan)
    for a,b in gkf.split(Xr,y,groups): p[b]=LinearRegression().fit(Xr[a],y[a]).predict(Xr[b])
    return mean_absolute_error(y,p),{z:mean_absolute_error(y[zone==z],p[zone==z]) for z in ("in","above")}

for tgt in ("d_lo","d_hi"):
    y=D[tgt].to_numpy(float)
    gbr=GradientBoostingRegressor(max_depth=3,n_estimators=80,learning_rate=0.1,subsample=0.6,random_state=0).fit(Xslate,y)
    rules=extract_rules(gbr,Xslate)
    names=SLATE+list(rules)
    X=np.column_stack([Xslate]+[rules[k] for k in rules]) if rules else Xslate
    Xs=(X-X.mean(0))/(X.std(0)+1e-12)
    alphas,coefs,_=lasso_path(Xs,y,alphas=100,max_iter=20000)
    supps=(np.abs(coefs)>1e-8).sum(0)
    print(f"\n=== {tgt}: RuleFit ({len(rules)} rules + {len(SLATE)} linear) ===",flush=True)
    best=None
    for nt in (10,15,20,30):
        j=int(np.argmin(np.abs(supps-nt))); supp=[i for i in range(len(names)) if abs(coefs[i,j])>1e-8]
        mae,z=cvmae(X[:,supp],y)
        if best is None or mae<best[0]: best=(mae,supp)
        print(f"  n={len(supp):2d}: CV-MAE={mae:.3f}  [in {z['in']:.2f} above {z['above']:.2f}]",flush=True)
    mf=LinearRegression().fit(X[:,best[1]],y)
    terms=sorted(zip([names[i] for i in best[1]],mf.coef_),key=lambda kv:-abs(kv[1]))
    print(f"  BEST {tgt} MAE={best[0]:.3f}:",flush=True)
    for n,c in terms[:10]: print(f"      {c:+.3f} · {n}",flush=True)
print("DONE",flush=True)
