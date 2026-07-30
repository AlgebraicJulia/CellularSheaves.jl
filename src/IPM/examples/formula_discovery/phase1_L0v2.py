#!/usr/bin/env python
"""PHASE 1 L0 v2 -- centered parametrization + unambiguous naming + champion head-to-head.
Fixes from review:
  #2 CENTER continuous vars at operating medians before forming products -> main effects keep physical
     slopes, interactions encode deviations, coefficients readable & fold-stable.
  #3 hinge naming: relu(f - knot); knot printed as its true location; activity % reported.
  #1 single support (selected on INNER folds) refit per OUTER fold -> support stable by construction;
     per-fold coefficient ranges reported. Champion (off-the-dome 7-term) scored on the SAME outer folds.
  #4 frontier {8,10,15}; tail-trim analysis on the largest.
Incremental + halt-safe: each milestone written to phase1_L0v2_results.json before continuing."""
import os, json, pickle, itertools, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
HERE=os.path.dirname(__file__)
RESJSON=os.path.join(HERE,"phase1_L0v2_results.json"); RESPKL=os.path.join(HERE,"phase1_L0v2_results.pkl")
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
RAW={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
     "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
ORDER=["la","cd","hx","bm","mg"]
MED={f:float(np.median(RAW[f])) for f in ORDER}              # operating medians
FEAT={f:RAW[f]-MED[f] for f in ORDER}                        # CENTERED base features
y=I["d_lo"].to_numpy(float); fam=I["problem"].str.split("_").str[0].to_numpy()
X5=np.column_stack([FEAT[f] for f in ORDER])

def family_folds(k,seed):
    fams=np.array(sorted(np.unique(fam))); rng=np.random.default_rng(seed); rng.shuffle(fams)
    sizes={f:int((fam==f).sum()) for f in fams}; load=[0]*k; assign={}
    for f in sorted(fams,key=lambda f:-sizes[f]):
        j=int(np.argmin(load)); assign[f]=j; load[j]+=sizes[f]
    a=np.array([assign[f] for f in fam]); return [(a!=j,a==j) for j in range(k)]
OUT=family_folds(5,0); INN=family_folds(4,1)
def cv(cols,folds):
    p=np.full(len(y),np.nan); coefs=[]
    for tr,te in folds:
        m=LinearRegression().fit(cols[tr],y[tr]); p[te]=m.predict(cols[te]); coefs.append([m.intercept_]+list(m.coef_))
    return mean_absolute_error(y,p),coefs
def score(sel): return cv(B[:,sel],INN)[0] if sel else 9.9

# ---- knots (interior, separated) on CENTERED features; store raw knot loc for naming ----
print("[knots] EBM (centered) ...",flush=True)
ebm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(X5,y)
KNOTS={}   # centered knot values
for i,f in enumerate(ORDER):
    cuts=np.asarray(ebm.bins_[i][0]).ravel(); sc=np.asarray(ebm.term_scores_[i]).ravel()
    n=min(len(cuts),len(sc)-2); lo,hi=np.quantile(FEAT[f],[0.10,0.90]); minsep=(hi-lo)/5.0
    if n>=4:
        body=sc[1:1+n]; d2=np.abs(np.diff(body,2))
        cand=sorted([(d2[j-1],float(cuts[j])) for j in range(1,n-1) if lo<=cuts[j]<=hi],key=lambda t:-t[0])
        kn=[]
        for _,c in cand:
            if all(abs(c-k)>=minsep for k in kn): kn.append(c)
            if len(kn)==2: break
        kn=sorted(kn) if kn else [0.0]
    else: kn=[0.0]
    KNOTS[f]=kn
    print(f"   {f}: raw knots {[round(k+MED[f],2) for k in kn]}",flush=True)

# ---- centered order<=3 basis; unambiguous names relu(f-rawknot) ----
T={}; NAMES=[]
for f in ORDER: T[f]=FEAT[f]; NAMES.append(f)         # centered linear, name = f  (means f-med)
HIN={}
for f in ORDER:
    for k in KNOTS[f]:
        rk=k+MED[f]; key=f"relu({f}-{rk:+.2f})"; T[key]=np.maximum(FEAT[f]-k,0.0); NAMES.append(key); HIN[key]=f
for a,b in itertools.combinations(ORDER,2): key=f"{a}*{b}"; T[key]=FEAT[a]*FEAT[b]; NAMES.append(key)
for a,b,c in itertools.combinations(ORDER,3): key=f"{a}*{b}*{c}"; T[key]=FEAT[a]*FEAT[b]*FEAT[c]; NAMES.append(key)
for hk,hf in HIN.items():
    for g in ORDER:
        if g!=hf: key=f"{hk}*{g}"; T[key]=T[hk]*FEAT[g]; NAMES.append(key)
B=np.column_stack([T[n] for n in NAMES])
ACT={n:float((T[n]!=0).mean()) for n in NAMES}
print(f"[basis] {len(NAMES)} centered terms  (linear names mean f-median; medians {({k:round(v,2) for k,v in MED.items()})})",flush=True)

def swap_refine(sel):
    sel=list(sel); improved=True
    while improved:
        improved=False; cur=score(sel)
        for si in list(sel):
            for j in range(len(NAMES)):
                if j in sel: continue
                cand=[x for x in sel if x!=si]+[j]
                if score(cand)<cur-1e-4: sel=cand; cur=score(cand); improved=True; break
            if improved: break
    return sel

MILE=[8,10,15]; results=[]
def save(): json.dump(results,open(RESJSON,"w"),indent=1); pickle.dump({"NAMES":NAMES,"KNOTS":KNOTS,"MED":MED,"ACT":ACT,"results":results},open(RESPKL,"wb"))

# ---- champion (off-the-dome) on the SAME outer folds ----
la,cd,hx,bm,mg=RAW["la"],RAW["cd"],RAW["hx"],RAW["bm"],RAW["mg"]
CX=np.column_stack([la,mg,np.maximum(bm-2.3,0),np.maximum(hx-2.8,0),np.maximum(cd+3.6,0),
                    np.maximum(cd+3.6,0)*mg,np.maximum(cd+10.7,0)*np.maximum(mg+1,0)])
champ_lit=mean_absolute_error(y,CX@np.array([0.49,-0.57,-0.77,0.54,0.16,-0.03,0.10])-1.50)
champ_refit,_=cv(CX,OUT)
print(f"\n[champion off-the-dome, 7 terms] literal={champ_lit:.3f}  structure-refit on MY outer folds={champ_refit:.3f}",flush=True)
results.append({"champion":{"terms":7,"literal":round(champ_lit,3),"refit_my_folds":round(champ_refit,3)}}); save()

print("\n[L0 v2] greedy-forward (once) + per-milestone swap+save:",flush=True)
sel=[]
for k in range(1,max(MILE)+1):
    best=None
    for j in range(len(NAMES)):
        if j in sel: continue
        m=score(sel+[j])
        if best is None or m<best[0]: best=(m,j)
    sel.append(best[1])
    if k in MILE:
        refined=swap_refine(sel); rep,coefs=cv(B[:,refined],OUT)
        terms=[NAMES[i] for i in refined]; fc=np.array(coefs)
        cmean=fc.mean(0); cmin=fc.min(0); cmax=fc.max(0)
        results.append({"n":k,"outer_cv_mae":round(rep,4),"gap_forest5":round(rep-0.936,4),
            "terms":terms,"coef_mean":[round(float(x),4) for x in cmean],
            "coef_range":[[round(float(a),4),round(float(b),4)] for a,b in zip(cmin,cmax)],
            "activity":[round(ACT[t],3) for t in terms]}); save()
        print(f"   n={k:2d}: outer-CV MAE={rep:.3f}  (gap forest5 {rep-0.936:+.3f}, vs champion {champ_refit:.3f})  [saved]",flush=True)
        # tail: how many terms carry |mean coef| < 0.02
        dead=[terms[i] for i in range(len(terms)) if abs(cmean[i+1])<0.02]
        print(f"        {len(dead)} terms with |coef|<0.02 (tail): {dead}",flush=True)
print("\nDONE",flush=True)
