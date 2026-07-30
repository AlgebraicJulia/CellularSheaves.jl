#!/usr/bin/env python
"""PHASE 1 L0 selection, INCREMENTAL + halt-safe.
Greedy-forward is run ONCE to 25 (nested prefixes give the frontier -- no 4x re-search).
After EACH milestone {6,10,15,25}: swap-refine a copy, score on OUTER folds, write results to disk.
So halting mid-run never loses completed milestones. Results -> phase1_L0_results.json (+ .pkl w/ coefs).
Selection scored on INNER family folds (seed 1); reporting on OUTER family folds (seed 0)."""
import os, json, pickle, itertools, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
HERE=os.path.dirname(__file__)
RESJSON=os.path.join(HERE,"phase1_L0_results.json"); RESPKL=os.path.join(HERE,"phase1_L0_results.pkl")
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
FEAT={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
      "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
ORDER=["la","cd","hx","bm","mg"]
y=I["d_lo"].to_numpy(float); fam=I["problem"].str.split("_").str[0].to_numpy()
X5=np.column_stack([FEAT[f] for f in ORDER])

def family_folds(k,seed):
    fams=np.array(sorted(np.unique(fam))); rng=np.random.default_rng(seed); rng.shuffle(fams)
    sizes={f:int((fam==f).sum()) for f in fams}; load=[0]*k; assign={}
    for f in sorted(fams,key=lambda f:-sizes[f]):
        j=int(np.argmin(load)); assign[f]=j; load[j]+=sizes[f]
    a=np.array([assign[f] for f in fam])
    return [(a!=j,a==j) for j in range(k)]
OUT=family_folds(5,0); INN=family_folds(4,1)
def cv(cols,folds):
    p=np.full(len(y),np.nan); coefs=[]
    for tr,te in folds:
        m=LinearRegression().fit(cols[tr],y[tr]); p[te]=m.predict(cols[te]); coefs.append([m.intercept_]+list(m.coef_))
    return mean_absolute_error(y,p),coefs
def score(sel): return cv(B[:,sel],INN)[0] if sel else 9.9

# ---- knots (interior, separated bends) ----
print("[knots] EBM ...",flush=True)
ebm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(X5,y)
KNOTS={}
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
        kn=sorted(kn) if kn else [float(np.quantile(FEAT[f],0.5))]
    else: kn=[float(np.quantile(FEAT[f],0.5))]
    KNOTS[f]=kn; print(f"   {f}: {[round(k,2) for k in kn]}",flush=True)

# ---- order<=3 basis with EBM knots ----
T={}; NAMES=[]
for f in ORDER: T[f]=FEAT[f]; NAMES.append(f)
HIN={}
for f in ORDER:
    for k in KNOTS[f]:
        key=f"({f}{k:+.2f})+"; T[key]=np.maximum(FEAT[f]-k,0.0); NAMES.append(key); HIN[key]=f
for a,b in itertools.combinations(ORDER,2): key=f"{a}*{b}"; T[key]=FEAT[a]*FEAT[b]; NAMES.append(key)
for a,b,c in itertools.combinations(ORDER,3): key=f"{a}*{b}*{c}"; T[key]=FEAT[a]*FEAT[b]*FEAT[c]; NAMES.append(key)
for hk,hf in HIN.items():
    for g in ORDER:
        if g!=hf: key=f"{hk}*{g}"; T[key]=T[hk]*FEAT[g]; NAMES.append(key)
B=np.column_stack([T[n] for n in NAMES])
print(f"[basis] {len(NAMES)} terms",flush=True)

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

MILE=[6,10,15,25]; results=[]
def save():
    json.dump(results,open(RESJSON,"w"),indent=1)
    pickle.dump({"NAMES":NAMES,"KNOTS":KNOTS,"results":results},open(RESPKL,"wb"))
print("\n[L0] greedy-forward (once) with per-milestone swap+save:",flush=True)
sel=[]  # greedy forward path (unperturbed)
for k in range(1,max(MILE)+1):
    best=None
    for j in range(len(NAMES)):
        if j in sel: continue
        m=score(sel+[j])
        if best is None or m<best[0]: best=(m,j)
    sel.append(best[1])
    if k in MILE:
        refined=swap_refine(sel)
        rep,coefs=cv(B[:,refined],OUT)
        terms=[NAMES[i] for i in refined]
        results.append({"n":k,"outer_cv_mae":round(rep,4),"inner_mae":round(score(refined),4),
                        "terms":terms,"fold_coefs":coefs,"gap_to_forest5":round(rep-0.936,4)})
        save()
        print(f"   n={k:2d}: outer-CV MAE={rep:.3f}  (gap to forest5 {rep-0.936:+.3f})  [saved]",flush=True)
        print(f"         terms: {terms}",flush=True)

# candidate off-the-dome (cheap, at the end)
la,cd,hx,bm,mg=FEAT["la"],FEAT["cd"],FEAT["hx"],FEAT["bm"],FEAT["mg"]
CX=np.column_stack([la,mg,np.maximum(bm-2.3,0),np.maximum(hx-2.8,0),np.maximum(cd+3.6,0),
                    np.maximum(cd+3.6,0)*mg,np.maximum(cd+10.7,0)*np.maximum(mg+1,0)])
Cc=np.array([0.49,-0.57,-0.77,0.54,0.16,-0.03,0.10]); C0=-1.50
lit=mean_absolute_error(y,CX@Cc+C0); refit,_=cv(CX,OUT)
cand={"literal":round(lit,3),"structure_refit_outer_cv":round(refit,3),"you_guessed":1.27}
results.append({"candidate_off_the_dome":cand}); save()
print(f"\n[candidate] literal={lit:.3f}  structure-refit outer-CV={refit:.3f}  (you guessed 1.27)  [saved]",flush=True)
print("DONE",flush=True)
