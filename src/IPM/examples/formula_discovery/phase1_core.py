#!/usr/bin/env python
"""PHASE 1 core: IPM/floor formula to a frozen candidate.
  0. family folds + calibration (forest5, forest46)
  1. knots from fast-EBM bend points on the five
  2. L0 greedy-forward + swap over order-<=3 basis; frontier {6,10,15,25}
     selection scored on INNER family folds (seed 1); reporting on OUTER family folds (seed 0)
  4. test the off-the-dome candidate (literal + structure-refit) side by side
Family = problem.split('_')[0]. All logs base-10 (dataset already log-transformed as L_*)."""
import os, itertools, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
np.set_printoptions(suppress=True)
HERE=os.path.dirname(__file__); D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
FEAT={"la":Lc("L_alpha"),
      "cd":Lc("L_c_res0_dual"),
      "hx":Lc("L_hdiag_max"),
      "bm":Lc("L_bar_hdiag_med"),
      "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
ORDER=["la","cd","hx","bm","mg"]
y=I["d_lo"].to_numpy(float); zone=I["zone"].to_numpy()
fam=I["problem"].str.split("_").str[0].to_numpy()

def family_folds(k,seed):
    """balanced folds holding out whole families; seed permutes assignment."""
    fams=np.array(sorted(np.unique(fam))); rng=np.random.default_rng(seed); rng.shuffle(fams)
    sizes={f:(fam==f).sum() for f in fams}; load=[0]*k; assign={}
    for f in sorted(fams,key=lambda f:-sizes[f]):
        j=int(np.argmin(load)); assign[f]=j; load[j]+=sizes[f]
    folds=[]
    for j in range(k):
        te=np.array([assign[f] for f in fam])==j
        folds.append((~te,te))
    return folds
OUT=family_folds(5,0); INN=family_folds(4,1)
def cv_mae(Xcols, folds):
    p=np.full(len(y),np.nan); coefs=[]
    for tr,te in folds:
        m=LinearRegression().fit(Xcols[tr],y[tr]); p[te]=m.predict(Xcols[te]); coefs.append((m.intercept_,*m.coef_))
    return mean_absolute_error(y,p),p,np.array(coefs)

# ---- 0. calibration ----
def forest_cv(X):
    p=np.full(len(y),np.nan)
    for tr,te in OUT:
        p[te]=HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,
              l2_regularization=1.0,random_state=0).fit(X[tr],y[tr]).predict(X[te])
    return mean_absolute_error(y,p)
X5=np.column_stack([FEAT[f] for f in ORDER])
print(f"IPM/floor: {len(y)} rows, {len(np.unique(fam))} families")
print(f"[calib] forest(5)  = {forest_cv(X5):.3f}   (brief 1.07)")
FULLC=[c for c in I.columns if c not in ("problem","file","is_hsd","d_lo","d_hi","zone","m_meta","n_meta","Lng_meta","L_alpha_anchor")]
Xfull=I[FULLC].fillna(I[FULLC].median()).to_numpy(float)
print(f"[calib] forest({len(FULLC)}) = {forest_cv(Xfull):.3f}   (brief 1.32)")

# ---- 1. knots from EBM bends (interior, separated) ----
print("\n[knots] fitting fast EBM on the five ...",flush=True)
ebm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(X5,y)
KNOTS={}
for i,f in enumerate(ORDER):
    cuts=np.asarray(ebm.bins_[i][0]).ravel()
    sc=np.asarray(ebm.term_scores_[i]).ravel()
    n=min(len(cuts),len(sc)-2);
    lo,hi=np.quantile(FEAT[f],[0.10,0.90]); minsep=(hi-lo)/5.0
    if n>=4:
        body=sc[1:1+n]; d2=np.abs(np.diff(body,2))  # |curvature| at interior cuts
        cand=[(d2[j-1],float(cuts[j])) for j in range(1,n-1) if lo<=cuts[j]<=hi]
        cand.sort(key=lambda t:-t[0]); kn=[]
        for _,c in cand:
            if all(abs(c-k)>=minsep for k in kn): kn.append(c)
            if len(kn)==2: break
        kn=sorted(kn) if kn else [float(np.quantile(FEAT[f],0.5))]
    else:
        kn=[float(np.quantile(FEAT[f],0.5))]
    KNOTS[f]=kn
    print(f"   {f}: knots {[round(k,2) for k in kn]}  (bulk [{lo:.1f},{hi:.1f}])",flush=True)

# ---- build order-<=3 basis with EBM knots ----
def build_basis():
    T={}; nm=[]
    for f in ORDER: T[f]=FEAT[f]; nm.append(f)
    HIN={}
    for f in ORDER:
        for k in KNOTS[f]:
            key=f"({f}{k:+.2f})+"; T[key]=np.maximum(FEAT[f]-k,0.0); nm.append(key); HIN[key]=f
    for a,b in itertools.combinations(ORDER,2):
        key=f"{a}*{b}"; T[key]=FEAT[a]*FEAT[b]; nm.append(key)
    for a,b,c in itertools.combinations(ORDER,3):
        key=f"{a}*{b}*{c}"; T[key]=FEAT[a]*FEAT[b]*FEAT[c]; nm.append(key)
    for hk,hf in HIN.items():                    # hinge x linear (other feats)
        for g in ORDER:
            if g==hf: continue
            key=f"{hk}*{g}"; T[key]=T[hk]*FEAT[g]; nm.append(key)
    M=np.column_stack([T[n] for n in nm]); return nm,M
NAMES,B=build_basis()
print(f"\n[basis] {len(NAMES)} candidate terms (order<=3, EBM knots)",flush=True)

# ---- 2. L0 greedy-forward + swap, scored on INNER folds ----
def score(cols):
    if not cols: return 9.9
    return cv_mae(B[:,cols],INN)[0]
def greedy_swap(target):
    sel=[]
    while len(sel)<target:
        best=None
        for j in range(len(NAMES)):
            if j in sel: continue
            m=score(sel+[j])
            if best is None or m<best[0]: best=(m,j)
        sel.append(best[1])
    improved=True
    while improved:
        improved=False; cur=score(sel)
        for si in list(sel):
            for j in range(len(NAMES)):
                if j in sel: continue
                cand=[x for x in sel if x!=si]+[j]; m=score(cand)
                if m<cur-1e-4: sel=cand; cur=m; improved=True; break
            if improved: break
    return sel
print("\n[L0] frontier (selection=inner folds, report=outer folds):",flush=True)
FRONT={}
for tgt in (6,10,15,25):
    sel=greedy_swap(tgt)
    rep,p,coefs=cv_mae(B[:,sel],OUT)
    FRONT[tgt]=(sel,rep,coefs)
    print(f"   n={tgt:2d}: outer-CV MAE={rep:.3f}   inner={score(sel):.3f}",flush=True)

# ---- 4. off-the-dome candidate ----
la,cd,hx,bm,mg=FEAT["la"],FEAT["cd"],FEAT["hx"],FEAT["bm"],FEAT["mg"]
CX=np.column_stack([la,mg,np.maximum(bm-2.3,0),np.maximum(hx-2.8,0),np.maximum(cd+3.6,0),
                    np.maximum(cd+3.6,0)*mg,np.maximum(cd+10.7,0)*np.maximum(mg+1,0)])
Ccoef=np.array([0.49,-0.57,-0.77,0.54,0.16,-0.03,0.10]); C0=-1.50
lit=mean_absolute_error(y,CX@Ccoef+C0)
refit,_,_=cv_mae(CX,OUT)
print(f"\n[candidate off-the-dome] literal(published coefs)={lit:.3f}   structure-refit outer-CV={refit:.3f}   (you guessed 1.27)")

# save frontier for the audit stage
import pickle
pickle.dump({"NAMES":NAMES,"FRONT":{k:v[0] for k,v in FRONT.items()},"KNOTS":KNOTS,
             "cand_lit":lit,"cand_refit":refit},open(os.path.join(HERE,"phase1_frontier.pkl"),"wb"))
print("\nsaved phase1_frontier.pkl  DONE",flush=True)
