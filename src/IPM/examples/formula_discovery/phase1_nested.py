#!/usr/bin/env python
"""STEP 3: NESTED L0 on canonical folds -- closes support-level leakage + doubles as stability audit.
For each of the 5 canonical OUTER folds: select the term support using ONLY that fold's TRAINING
families (inner 4-fold family CV within train), refit coefficients on outer-train, score outer-test.
Report: honest outer-CV MAE at sizes {8,10,12}; per-term SELECTION COUNT across the 5 folds
(5/5 = the formula; 2/5 = weather). Champion rescored on the same canonical folds.
Greedy-forward (swap omitted in nested for tractability -- noted). Centered basis, relu(f-knot) names.
Halt-safe: writes phase1_nested_results.json after each outer fold."""
import os, re, json, itertools, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
HERE=os.path.dirname(__file__); RES=os.path.join(HERE,"phase1_nested_results.json")
CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in I["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam])
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
RAW={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
     "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
ORDER=["la","cd","hx","bm","mg"]; MED={f:float(np.median(RAW[f])) for f in ORDER}
FEAT={f:RAW[f]-MED[f] for f in ORDER}; y=I["d_lo"].to_numpy(float)
X5=np.column_stack([FEAT[f] for f in ORDER])

# knots (centered) via EBM on full data (knot LOCATIONS are a basis choice, fixed; only term SELECTION is nested)
ebm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(X5,y)
KNOTS={}
for i,f in enumerate(ORDER):
    cuts=np.asarray(ebm.bins_[i][0]).ravel(); sc=np.asarray(ebm.term_scores_[i]).ravel()
    n=min(len(cuts),len(sc)-2); lo,hi=np.quantile(FEAT[f],[0.10,0.90]); minsep=(hi-lo)/5.0
    if n>=4:
        body=sc[1:1+n]; d2=np.abs(np.diff(body,2))
        cand=sorted([(d2[j-1],float(cuts[j])) for j in range(1,n-1) if lo<=cuts[j]<=hi],key=lambda t:-t[0]); kn=[]
        for _,c in cand:
            if all(abs(c-k)>=minsep for k in kn): kn.append(c)
            if len(kn)==2: break
        kn=sorted(kn) if kn else [0.0]
    else: kn=[0.0]
    KNOTS[f]=kn
# centered basis
T={}; NAMES=[]
for f in ORDER: T[f]=FEAT[f]; NAMES.append(f)
HIN={}
for f in ORDER:
    for k in KNOTS[f]:
        rk=k+MED[f]; key=f"relu({f}-{rk:+.2f})"; T[key]=np.maximum(FEAT[f]-k,0.0); NAMES.append(key); HIN[key]=f
for a,b in itertools.combinations(ORDER,2): key=f"{a}*{b}"; T[key]=FEAT[a]*FEAT[b]; NAMES.append(key)
for a,b,c in itertools.combinations(ORDER,3): key=f"{a}*{b}*{c}"; T[key]=FEAT[a]*FEAT[b]*FEAT[c]; NAMES.append(key)
for hk,hf in HIN.items():
    for g in ORDER:
        if g!=hf: key=f"{hk}*{g}"; T[key]=T[hk]*FEAT[g]; NAMES.append(key)
B=np.column_stack([T[n] for n in NAMES]); ACT={n:float((T[n]!=0).mean()) for n in NAMES}
print(f"[nested] {len(NAMES)} centered terms, canonical folds ({CF['n_families']} families)",flush=True)

def inner_folds(trainmask):
    tf=sorted(np.unique(fam[trainmask])); load=[0]*4; asg={}
    for f in sorted(tf,key=lambda f:-int((fam==f).sum())):
        j=int(np.argmin(load)); asg[f]=j; load[j]+=int((fam==f).sum())
    ff=np.array([asg.get(x,-1) for x in fam])
    return [((trainmask)&(ff!=j),(trainmask)&(ff==j)) for j in range(4)]
def cvmae(cols,folds):
    p=np.full(len(y),np.nan); cnt=0
    for tr,te in folds:
        if te.sum()==0: continue
        p[te]=LinearRegression().fit(cols[tr],y[tr]).predict(cols[te])
    m=~np.isnan(p); return mean_absolute_error(y[m],p[m])
def greedy(trainmask,target):
    inf=inner_folds(trainmask); sel=[]
    while len(sel)<target:
        best=None
        for j in range(len(NAMES)):
            if j in sel: continue
            m=cvmae(B[:,sel+[j]],inf)
            if best is None or m<best[0]: best=(m,j)
        sel.append(best[1])
    return sel

SIZES=[8,10,12]; MAXT=max(SIZES)
pred={s:np.full(len(y),np.nan) for s in SIZES}; supports={s:[] for s in SIZES}
results={"note":"nested; greedy-forward, swap omitted; canonical folds","sizes":SIZES,"folds":[]}
def save(): json.dump(results,open(RES,"w"),indent=1)
for o in range(CF["k"]):
    tr=foldid!=o; te=foldid==o
    path=greedy(tr,MAXT)                       # nested selection on train only
    foldrec={"outer":o,"test_families":sorted(np.unique(fam[te]).tolist())}
    for s in SIZES:
        sup=path[:s]; supports[s].append([NAMES[i] for i in sup])
        m=LinearRegression().fit(B[np.ix_(tr,sup)],y[tr]); pred[s][te]=m.predict(B[np.ix_(te,sup)])
        foldrec[f"support_{s}"]=[NAMES[i] for i in sup]
        foldrec[f"test_mae_{s}"]=round(float(mean_absolute_error(y[te],pred[s][te])),4)
    results["folds"].append(foldrec); save()
    print(f"  outer fold {o} done (test fams {foldrec['test_families']}): "
          f"mae8={foldrec['test_mae_8']} mae10={foldrec['test_mae_10']} mae12={foldrec['test_mae_12']}",flush=True)

# aggregate
agg={}
for s in SIZES:
    mae=float(mean_absolute_error(y,pred[s]))
    from collections import Counter
    cnt=Counter(t for sup in supports[s] for t in sup)
    agg[s]={"outer_cv_mae":round(mae,4),"gap_forest5":round(mae-1.1087,4),
            "selection_counts":sorted(cnt.items(),key=lambda kv:(-kv[1],kv[0]))}
    print(f"\n=== size {s}: nested outer-CV MAE={mae:.3f} (gap to forest5 {mae-1.1087:+.3f}) ===",flush=True)
    for t,c in agg[s]["selection_counts"]:
        tag="FORMULA" if c==5 else ("weather" if c<=2 else "")
        print(f"   {c}/5  {t:26} (active {ACT[t]*100:.0f}%)  {tag}",flush=True)
# champion on canonical folds
la,cd,hx,bm,mg=RAW["la"],RAW["cd"],RAW["hx"],RAW["bm"],RAW["mg"]
CX=np.column_stack([la,mg,np.maximum(bm-2.3,0),np.maximum(hx-2.8,0),np.maximum(cd+3.6,0),
                    np.maximum(cd+3.6,0)*mg,np.maximum(cd+10.7,0)*np.maximum(mg+1,0)])
pc=np.full(len(y),np.nan)
for o in range(CF["k"]):
    tr=foldid!=o; te=foldid==o; pc[te]=LinearRegression().fit(CX[tr],y[tr]).predict(CX[te])
champ=float(mean_absolute_error(y,pc))
agg["champion_canonical_cv"]=round(champ,4)
results["aggregate"]=agg; save()
print(f"\n[champion off-the-dome, 7 terms] canonical-fold CV MAE = {champ:.3f}",flush=True)
print("DONE",flush=True)
