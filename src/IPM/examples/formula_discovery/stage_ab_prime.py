#!/usr/bin/env python
"""Stage A'/B' — composite-augmented rerun (before Stage D). Slate = survivors + materialized composites
(u=logα+logσ²min, κ̂=hdmax−σ²min, σ²min−ftol, margin, logN, tolratio, top-3 FAST pairs/target as products).
(1) depth-4 forest CV; (2) capacity ladder (headline: depth4→additive gap); (3) EBM shapes, esp. u.
Decision: additive within ~0.05 of depth-4 ⇒ enrich Stage C, Stage D confirm-only; gap ≥0.15 ⇒ Stage D full."""
import os, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor

HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
FP=["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
SURV=[c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")+tuple(FP)]
L=lambda c: D[c].fillna(D[c].median()).to_numpy(float)
COMP={
 "u=la+s2":       L("L_alpha")+L("L_sigma2min"),
 "khat=hdmx-s2":  L("L_hdiag_max")-L("L_sigma2min"),
 "s2-ftol":       L("L_sigma2min")-L("L_force_tol"),
 "margin":        L("L_r0_c")-np.maximum(L("L_force_tol"),L("L_floor_tol")),
 "logN":          L("L_r0_c")+L("L_alpha"),
 "tolratio":      L("L_force_tol")-L("L_floor_tol"),
 "p_hdmx*s2":     L("L_hdiag_max")*L("L_sigma2min"),   # d_hi pairs
 "p_la*hdmn":     L("L_alpha")*L("L_hdiag_min"),
 "p_ft*s2":       L("L_force_tol")*L("L_sigma2min"),
 "p_hdmn*hdmx":   L("L_hdiag_min")*L("L_hdiag_max"),   # d_lo pairs
 "p_r0p*s2":      L("L_r0_p")*L("L_sigma2min"),
 "p_hdmx*ishsd":  L("L_hdiag_max")*L("is_hsd"),
}
for k,v in COMP.items(): D[k]=v
SLATE=SURV+list(COMP)
groups=D["problem"].values; gkf=GroupKFold(5)

def cv(mk,target,feats):
    X=D[feats].to_numpy(float); y=D[target].to_numpy(float); p=np.full(len(y),np.nan)
    for a,b in gkf.split(X,y,groups):
        m=mk(); m.fit(X[a],y[a]); p[b]=m.predict(X[b])
    return mean_absolute_error(y,p)
hgb=lambda d:(lambda:HistGradientBoostingRegressor(max_depth=d,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0))
lin=lambda:make_pipeline(SimpleImputer(),Ridge(1.0))

print("=== A'/B' composite-augmented (grouped-CV MAE) ===",flush=True)
for tgt in ("d_lo","d_hi"):
    d4=cv(hgb(4),tgt,SLATE); d1=cv(hgb(1),tgt,SLATE); li=cv(lin,tgt,SLATE)
    print(f"  {tgt}: depth4={d4:.3f}  depth1(add)={d1:.3f}  linear={li:.3f}  | gap(d4->add)={d1-d4:.3f}",flush=True)

# EBM fast + u shape
for tgt in ("d_hi","d_lo"):
    X=D[SLATE].to_numpy(float); y=D[tgt].to_numpy(float)
    m=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(X,y)
    p=np.full(len(y),np.nan)
    for a,b in gkf.split(X,y,groups):
        mm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(X[a],y[a]); p[b]=mm.predict(X[b])
    imp=sorted(zip(SLATE,m.term_importances()),key=lambda kv:-kv[1])
    print(f"\n=== EBM {tgt}: CV-MAE={mean_absolute_error(y,p):.3f}  top: "+", ".join(f"{n}({i:.2f})" for n,i in imp[:7]),flush=True)
    fi=SLATE.index("u=la+s2"); tf=m.term_features_
    ti=[k for k in range(len(tf)) if len(tf[k])==1 and tf[k][0]==fi][0]
    cuts=np.asarray(m.bins_[ti][0]).ravel(); sc=np.asarray(m.term_scores_[ti]).ravel(); body=sc[1:1+len(cuts)+1]
    sel=np.linspace(0,len(body)-1,min(9,len(body))).astype(int)
    xl=lambda j:(f"<{cuts[0]:.1f}" if j==0 else (f">{cuts[-1]:.1f}" if j>=len(cuts) else f"{cuts[j-1]:.1f}"))
    print(f"    shape u (=logα+logσ²min; theory bends u≈0 and u≈−log(ε·κ̂)): "+"  ".join(f"[{xl(j)}]{body[j]:+.2f}" for j in sel),flush=True)
print("DONE",flush=True)
