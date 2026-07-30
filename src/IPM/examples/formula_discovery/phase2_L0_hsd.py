#!/usr/bin/env python
"""PHASE 2 rung (iv): free HSD nested L0 for one sensor. Usage: phase2_L0_hsd.py {floor|ceiling}.
HSD rows, canonical folds, features from the HSD RFE, nested (leak-free) greedy L0, freeze >=3/5 core,
drop-test trim, audits. Mirrors the IPM pipeline exactly."""
import os, re, sys, json, itertools, numpy as np, pandas as pd
from collections import Counter
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
SENSOR=sys.argv[1]; HERE=os.path.dirname(__file__); CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy(); H=D[D["is_hsd"]==1].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in H["problem"]])
foldid=np.array([CF["family_to_fold"][f] for f in fam]); zone=H["zone"].to_numpy()
Lc=lambda c: H[c].fillna(H[c].median()).to_numpy(float)
mg=np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))
if SENSOR=="floor":
    TGT="d_lo"; CONTdef={"la":"L_alpha","ft":"L_force_tol","pd":"L_p_res0_dual","Ld":"L_Ldiag_min"}
    CNT=("pbase","pb"); REF=["L_alpha","L_force_tol","L_p_res0_dual"]; REFmg=True
else:
    TGT="d_hi"; CONTdef={"Ld":"L_Ldiag_min","la":"L_alpha","cd":"L_c_res0_dual","ft":"L_force_tol","r0p":"L_r0_p"}
    CNT=("ppass","pp"); REF=["L_Ldiag_min","L_alpha","L_c_res0_dual","L_force_tol"]; REFmg=False
y=H[TGT].to_numpy(float); exX04=~(fam=="X04")
RAWc={k:Lc(v) for k,v in CONTdef.items()}
if SENSOR=="floor": RAWc["mg"]=mg
CONT=list(RAWc); MED={f:float(np.median(RAWc[f])) for f in CONT}; C={f:RAWc[f]-MED[f] for f in CONT}
cnt=H[CNT[0]].fillna(0).to_numpy(float); cntmed=float(np.median(cnt)); C[CNT[1]]=cnt-cntmed
ORDER=CONT+[CNT[1]]; Xc=np.column_stack([C[f] for f in ORDER])
mk=lambda: HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0)
def forestCV(X):
    p=np.full(len(y),np.nan)
    for j in range(CF["k"]): tr=foldid!=j; p[foldid==j]=mk().fit(X[tr],y[tr]).predict(X[foldid==j])
    return mean_absolute_error(y,p)
Xref=np.column_stack([Lc(c) for c in REF]+([mg] if REFmg else [])+[cnt])
print(f"[HSD {SENSOR}] {len(y)} rows. forest ref (essential+counter) = {forestCV(Xref):.3f}",flush=True)
def cv(cols,split):
    p=np.full(len(y),np.nan)
    for tr,te in split:
        if te.sum()==0: continue
        p[te]=LinearRegression().fit(cols[tr],y[tr]).predict(cols[te])
    m=~np.isnan(p); return mean_absolute_error(y[m],p[m])
ebm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=3000,random_state=0).fit(Xc,y)
KN={}
for i,f in enumerate(CONT):
    cuts=np.asarray(ebm.bins_[i][0]).ravel(); sc=np.asarray(ebm.term_scores_[i]).ravel()
    n=min(len(cuts),len(sc)-2); lo,hi=np.quantile(C[f],[0.10,0.90]); minsep=(hi-lo)/5.0
    if n>=4:
        body=sc[1:1+n]; d2=np.abs(np.diff(body,2))
        cand=sorted([(d2[j-1],float(cuts[j])) for j in range(1,n-1) if lo<=cuts[j]<=hi],key=lambda t:-t[0]); kn=[]
        for _,cc in cand:
            if all(abs(cc-k)>=minsep for k in kn): kn.append(cc)
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
    key=f"[{CNT[0]}>={thr}]"; T[key]=(cnt>=thr).astype(float); NAMES.append(key)
for a,b in itertools.combinations(ORDER,2): T[f"{a}*{b}"]=C[a]*C[b]; NAMES.append(f"{a}*{b}")
for a,b,c in itertools.combinations(ORDER,3): T[f"{a}*{b}*{c}"]=C[a]*C[b]*C[c]; NAMES.append(f"{a}*{b}*{c}")
B=np.column_stack([T[n] for n in NAMES]); ACT={n:float((T[n]!=0).mean()) for n in NAMES}
print(f"[basis] {len(NAMES)} terms over {ORDER}",flush=True)
def inner_folds(tm):
    tf=sorted(np.unique(fam[tm])); load=[0]*4; asg={}
    for f in sorted(tf,key=lambda f:-int((fam==f).sum())): j=int(np.argmin(load)); asg[f]=j; load[j]+=int((fam==f).sum())
    ff=np.array([asg.get(x,-1) for x in fam]); return [((tm)&(ff!=j),(tm)&(ff==j)) for j in range(4)]
def greedy(tm,target):
    inf=inner_folds(tm); sel=[]
    while len(sel)<target:
        best=None
        for j in range(len(NAMES)):
            if j in sel: continue
            m=cv(B[:,sel+[j]],inf)
            if best is None or m<best[0]: best=(m,j)
        sel.append(best[1])
    return sel
SIZES=[8,10,12]; pred={s:np.full(len(y),np.nan) for s in SIZES}; supports={s:[] for s in SIZES}
for o in range(CF["k"]):
    tr=foldid!=o; te=foldid==o; path=greedy(tr,12)
    for s in SIZES:
        sup=path[:s]; supports[s].append([NAMES[i] for i in sup]); pred[s][te]=LinearRegression().fit(B[np.ix_(tr,sup)],y[tr]).predict(B[np.ix_(te,sup)])
    print(f"  outer {o} ({sorted(np.unique(fam[te]))[:3]}...): mae10={mean_absolute_error(y[te],pred[10][te]):.3f}",flush=True)
for s in SIZES:
    cnt2=Counter(t for sup in supports[s] for t in sup)
    print(f"\n=== HSD {SENSOR} size {s}: nested CV MAE={mean_absolute_error(y,pred[s]):.3f} ===",flush=True)
    for t,c in sorted(cnt2.items(),key=lambda kv:(-kv[1],kv[0])):
        if c>=2: print(f"   {c}/5  {t:20} (act {ACT[t]*100:.0f}%){'  FORMULA' if c==5 else ''}",flush=True)
# freeze >=3/5 at size 12, drop-test trim
agg12=Counter(t for sup in supports[12] for t in sup); stable=[t for t,c in agg12.items() if c>=3]
Sidx=[NAMES.index(t) for t in stable]
def cvp(idx):
    p=np.full(len(y),np.nan)
    for o in range(CF["k"]): tr=foldid!=o; p[foldid==o]=LinearRegression().fit(B[np.ix_(tr,idx)],y[tr]).predict(B[np.ix_(foldid==o,idx)])
    return p
p=cvp(Sidx); base=mean_absolute_error(y,p)
# iteratively drop terms whose removal doesn't hurt (>= -0.002)
keep=list(Sidx); improved=True
while improved and len(keep)>3:
    improved=False
    for i in list(keep):
        d=mean_absolute_error(y,cvp([k for k in keep if k!=i]))-mean_absolute_error(y,cvp(keep))
        if d<0.003: keep=[k for k in keep if k!=i]; improved=True; break
p=cvp(keep); full=LinearRegression().fit(B[:,keep],y)
allm=mean_absolute_error(y,p); exm=mean_absolute_error(y[exX04],p[exX04])
print(f"\n=== FROZEN HSD {SENSOR} ({len(keep)} terms, drop-test trimmed) canonical CV: all={allm:.3f} ex-X04={exm:.3f} ===",flush=True)
print(f"   {TGT} = {full.intercept_:+.3f}",flush=True)
for r,i in sorted(enumerate(keep),key=lambda ri:-abs(full.coef_[ri[0]])):
    print(f"       {full.coef_[r]:+.3f} · {NAMES[i]:20} act {ACT[NAMES[i]]*100:.0f}%",flush=True)
out={"sensor":SENSOR,"target":TGT,"solver":"hsd","cv_all":round(float(allm),3),"cv_exX04":round(float(exm),3),
     "intercept":round(float(full.intercept_),4),"terms":{NAMES[i]:round(float(full.coef_[r]),4) for r,i in enumerate(keep)},
     "medians":{**{k:round(v,3) for k,v in MED.items()},CNT[1]:cntmed},"folds":"canonical_folds.json",
     "by_zone":{z:round(float(mean_absolute_error(y[zone==z],p[zone==z])),3) for z in ("below","in","above") if (zone==z).any()}}
json.dump(out,open(os.path.join(HERE,f"phase2_frozen_hsd_{SENSOR}.json"),"w"),indent=1)
print(f"\nsaved phase2_frozen_hsd_{SENSOR}.json  DONE",flush=True)
