#!/usr/bin/env python
"""IPM CEILING (d_hi) L0 -- nested, leak-free, canonical folds. Same recipe as the floor freeze.
Features (from ceiling RFE): Ld=log Ldiag_min, cd=log c_res0_dual, ft=log force_tol, pd=log p_res0_dual,
la=log alpha (all centered continuous) + cp=cpass (COUNTER -> staircase indicators).
Basis order<=3: centered linear + EBM-knot hinges + cpass indicators + pairwise + triples.
Reports forest reference, nested frontier {8,10,12}, selection stability, then freezes the >=3/5 core."""
import os, re, json, itertools, numpy as np, pandas as pd
from collections import Counter
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
HERE=os.path.dirname(__file__); RES=os.path.join(HERE,"phase1_L0_ceiling_results.json")
CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); I=D[D["is_hsd"]==0].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in I["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam]); y=I["d_hi"].to_numpy(float); zone=I["zone"].to_numpy()
Lc=lambda c: I[c].fillna(I[c].median()).to_numpy(float)
RAWc={"Ld":Lc("L_Ldiag_min"),"cd":Lc("L_c_res0_dual"),"ft":Lc("L_force_tol"),"pd":Lc("L_p_res0_dual"),"la":Lc("L_alpha")}
CONT=["Ld","cd","ft","pd","la"]; MED={f:float(np.median(RAWc[f])) for f in CONT}; C={f:RAWc[f]-MED[f] for f in CONT}
cp=I["cpass"].fillna(0).to_numpy(float); cpmed=float(np.median(cp)); C["cp"]=cp-cpmed
ORDER=CONT+["cp"]; Xc=np.column_stack([C[f] for f in ORDER]); exX04=~(fam=="X04")
mk=lambda: HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0)
def forest_ref(feats):
    Xf=np.column_stack([RAWc.get(f,None) if f in RAWc else (cp if f=="cp" else None) for f in feats])
    p=np.full(len(y),np.nan)
    for j in range(CF["k"]):
        tr=foldid!=j; p[foldid==j]=mk().fit(Xf[tr],y[tr]).predict(Xf[foldid==j])
    return mean_absolute_error(y,p)
print(f"[ceiling] forest ref: 4-knee {{Ld,cd,ft,cp}}={forest_ref(['Ld','cd','ft','cp']):.3f}  6-feat={forest_ref(ORDER):.3f}",flush=True)
def cv(cols,split):
    p=np.full(len(y),np.nan)
    for tr,te in split:
        if te.sum()==0: continue
        p[te]=LinearRegression().fit(cols[tr],y[tr]).predict(cols[te])
    m=~np.isnan(p); return mean_absolute_error(y[m],p[m])
print("[knots] EBM ...",flush=True)
ebm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(Xc,y)
KN={}
for i,f in enumerate(CONT):
    cuts=np.asarray(ebm.bins_[i][0]).ravel(); sc=np.asarray(ebm.term_scores_[i]).ravel()
    n=min(len(cuts),len(sc)-2); lo,hi=np.quantile(C[f],[0.10,0.90]); minsep=(hi-lo)/5.0
    if n>=4:
        body=sc[1:1+n]; d2=np.abs(np.diff(body,2))
        cand=sorted([(d2[j-1],float(cuts[j])) for j in range(1,n-1) if lo<=cuts[j]<=hi],key=lambda t:-t[0]); kn=[]
        for _,c in cand:
            if all(abs(c-k)>=minsep for k in kn): kn.append(c)
            if len(kn)==2: break
        kn=sorted(kn) if kn else [0.0]
    else: kn=[0.0]
    KN[f]=kn
T={}; NAMES=[]
for f in ORDER: T[f]=C[f]; NAMES.append(f)
for f in CONT:
    for k in KN[f]:
        rk=k+MED[f]; key=f"relu({f}-{rk:+.2f})"; T[key]=np.maximum(C[f]-k,0.0); NAMES.append(key)
for thr in (2,3,5):
    key=f"[cpass>={thr}]"; T[key]=(cp>=thr).astype(float); NAMES.append(key)
for a,b in itertools.combinations(ORDER,2): key=f"{a}*{b}"; T[key]=C[a]*C[b]; NAMES.append(key)
for a,b,c in itertools.combinations(ORDER,3): key=f"{a}*{b}*{c}"; T[key]=C[a]*C[b]*C[c]; NAMES.append(key)
B=np.column_stack([T[n] for n in NAMES]); ACT={n:float((T[n]!=0).mean()) for n in NAMES}
print(f"[basis] {len(NAMES)} terms",flush=True)
def inner_folds(trainmask):
    tf=sorted(np.unique(fam[trainmask])); load=[0]*4; asg={}
    for f in sorted(tf,key=lambda f:-int((fam==f).sum())):
        j=int(np.argmin(load)); asg[f]=j; load[j]+=int((fam==f).sum())
    ff=np.array([asg.get(x,-1) for x in fam]); return [((trainmask)&(ff!=j),(trainmask)&(ff==j)) for j in range(4)]
def greedy(trainmask,target):
    inf=inner_folds(trainmask); sel=[]
    while len(sel)<target:
        best=None
        for j in range(len(NAMES)):
            if j in sel: continue
            m=cv(B[:,sel+[j]],inf)
            if best is None or m<best[0]: best=(m,j)
        sel.append(best[1])
    return sel
SIZES=[8,10,12]; pred={s:np.full(len(y),np.nan) for s in SIZES}; supports={s:[] for s in SIZES}
results={"note":"nested ceiling L0 d_hi, canonical folds","sizes":SIZES,"folds":[]}
def save(): json.dump(results,open(RES,"w"),indent=1)
for o in range(CF["k"]):
    tr=foldid!=o; te=foldid==o; path=greedy(tr,12); rec={"outer":o,"test_families":sorted(np.unique(fam[te]).tolist())}
    for s in SIZES:
        sup=path[:s]; supports[s].append([NAMES[i] for i in sup])
        m=LinearRegression().fit(B[np.ix_(tr,sup)],y[tr]); pred[s][te]=m.predict(B[np.ix_(te,sup)])
        rec[f"support_{s}"]=[NAMES[i] for i in sup]; rec[f"test_mae_{s}"]=round(float(mean_absolute_error(y[te],pred[s][te])),4)
    results["folds"].append(rec); save()
    print(f"  outer {o} ({rec['test_families']}): mae8={rec['test_mae_8']} mae10={rec['test_mae_10']} mae12={rec['test_mae_12']}",flush=True)
agg={}
for s in SIZES:
    mae=float(mean_absolute_error(y,pred[s])); cnt=Counter(t for sup in supports[s] for t in sup)
    agg[s]={"nested_cv_mae":round(mae,4),"selection_counts":sorted(cnt.items(),key=lambda kv:(-kv[1],kv[0]))}
    print(f"\n=== size {s}: nested outer-CV MAE={mae:.3f} ===",flush=True)
    for t,c in agg[s]["selection_counts"]:
        print(f"   {c}/5  {t:22} (act {ACT[t]*100:.0f}%){'  FORMULA' if c==5 else ('  weather' if c<=2 else '')}",flush=True)
stable=[t for t,c in agg[12]["selection_counts"] if c>=3]; Sidx=[NAMES.index(t) for t in stable]
pfix=np.full(len(y),np.nan); coefs=[]
for o in range(CF["k"]):
    trm=foldid!=o; tem=foldid==o; m=LinearRegression().fit(B[np.ix_(trm,Sidx)],y[trm]); pfix[tem]=m.predict(B[np.ix_(tem,Sidx)]); coefs.append([m.intercept_]+list(m.coef_))
allm=float(mean_absolute_error(y,pfix)); exm=float(mean_absolute_error(y[exX04],pfix[exX04])); coefs=np.array(coefs)
full=LinearRegression().fit(B[:,Sidx],y)
print(f"\n=== FROZEN ceiling stable-core ({len(stable)} terms, >=3/5) canonical CV: all={allm:.3f}  ex-X04={exm:.3f} ===",flush=True)
print(f"   d_hi = {full.intercept_:+.3f}",flush=True)
for r,i in sorted(enumerate(Sidx),key=lambda ri:-abs(full.coef_[ri[0]])):
    print(f"       {full.coef_[r]:+.3f} · {NAMES[i]:22} [fold {coefs[:,r+1].min():+.3f},{coefs[:,r+1].max():+.3f}] act {ACT[NAMES[i]]*100:.0f}%",flush=True)
results["frozen_core"]={"terms":stable,"cv_all":round(allm,3),"cv_exX04":round(exm,3),"intercept":round(float(full.intercept_),4),
    "coefs":{NAMES[i]:round(float(full.coef_[r]),4) for r,i in enumerate(Sidx)},"medians":{**{k:round(v,3) for k,v in MED.items()},"cp":cpmed}}
results["aggregate"]=agg; save()
print("\nsaved phase1_L0_ceiling_results.json  DONE",flush=True)
