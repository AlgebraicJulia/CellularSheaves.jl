#!/usr/bin/env python
"""PHASE 2 -- HSD transfer ladder for both sensors (floor d_lo, ceiling d_hi), canonical folds (same
family->fold map; applied to HSD rows). Rungs:
  (i)   IPM formula VERBATIM (IPM structure + IPM knots + IPM coefs + IPM medians) scored on HSD -> transfer floor
  (ii)  same structure (IPM knots), COEFFICIENTS refit on HSD (canonical CV)
  (iii) COEFFICIENTS + KNOTS refit (knots re-derived from HSD EBM), canonical CV
  ref   HSD forest on the sensor's essential features (achievability bound; caps what free-rung (iv) can do)
Gaps between rungs decompose solver differences: (i)->(ii)=levels/offsets, (ii)->(iii)=thresholds,
(iii)->ref=structure (=> whether (iv) free redo is worth running)."""
import os, re, json, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
from interpret.glassbox import ExplainableBoostingRegressor
HERE=os.path.dirname(__file__); CF=json.load(open(os.path.join(HERE,"canonical_folds.json")))
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl")).copy()
H=D[D["is_hsd"]==1].copy().reset_index(drop=True)
fam=np.array([(re.match(r"^(X\d+)",str(p).split("_")[0]).group(1) if re.match(r"^(X\d+)",str(p).split("_")[0]) else str(p).split("_")[0]) for p in H["problem"]])
foldid=np.array([CF["family_to_fold"].get(f,-1) for f in fam])
print(f"HSD rows: {len(H)}  families present: {sorted(set(fam))}  (unmapped: {sorted(set(fam[foldid<0]))})")
Lc=lambda c: H[c].fillna(H[c].median()).to_numpy(float)
def canonCV(X,yv):
    p=np.full(len(yv),np.nan)
    for j in range(CF["k"]):
        tr=foldid!=j; te=foldid==j
        if te.sum()==0: continue
        p[te]=LinearRegression().fit(X[tr],yv[tr]).predict(X[te])
    m=~np.isnan(p); return mean_absolute_error(yv[m],p[m])
def forestCV(X,yv):
    p=np.full(len(yv),np.nan)
    for j in range(CF["k"]):
        tr=foldid!=j; te=foldid==j
        if te.sum()==0: continue
        p[te]=HistGradientBoostingRegressor(max_depth=4,max_iter=300,learning_rate=0.06,l2_regularization=1.0,random_state=0).fit(X[tr],yv[tr]).predict(X[te])
    m=~np.isnan(p); return mean_absolute_error(yv[m],p[m])
def hsd_knot(feat,fallback):
    x=Lc(feat); ebm=ExplainableBoostingRegressor(interactions=0,outer_bags=4,max_rounds=1500,random_state=0).fit( x.reshape(-1,1), yv )
    cuts=np.asarray(ebm.bins_[0][0]).ravel(); sc=np.asarray(ebm.term_scores_[0]).ravel()
    n=min(len(cuts),len(sc)-2); lo,hi=np.quantile(x,[0.10,0.90])
    if n>=4:
        body=sc[1:1+n]; d2=np.abs(np.diff(body,2))
        cand=[(d2[j-1],float(cuts[j])) for j in range(1,n-1) if lo<=cuts[j]<=hi]
        if cand: return sorted(cand,key=lambda t:-t[0])[0][1]
    return fallback
LAD={}

# ---------- FLOOR ----------
yv=H["d_lo"].to_numpy(float)
FF=json.load(open(os.path.join(HERE,"phase1_frozen_formula.json"))); Mf=FF["medians"]
raw={"la":Lc("L_alpha"),"cd":Lc("L_c_res0_dual"),"hx":Lc("L_hdiag_max"),"bm":Lc("L_bar_hdiag_med"),
     "mg":np.maximum(Lc("L_r0_c"),Lc("L_r0_p"))-np.maximum(Lc("L_force_tol"),Lc("L_floor_tol"))}
def floor_terms(knot_cd, med):
    C={f:raw[f]-med[f] for f in raw}; H_=np.maximum(raw["cd"]-knot_cd,0.0)
    return {"la":C["la"],"mg":C["mg"],f"reluCd*bm":H_*C["bm"],f"reluCd*la":H_*C["la"],
            "hx*bm":C["hx"]*C["bm"],"cd*hx*mg":C["cd"]*C["hx"]*C["mg"],"la*cd*bm":C["la"]*C["cd"]*C["bm"]}
# (i) verbatim: IPM coefs + IPM medians + IPM knot -10.52
Tv=floor_terms(-10.52,Mf); order=["la","mg","reluCd*bm","reluCd*la","hx*bm","cd*hx*mg","la*cd*bm"]
coef=[t["coef_alldata"] for t in FF["terms"]]; nm=[t["term"] for t in FF["terms"]]
cmap={"la":"la","mg":"mg","relu(cd-10.52)*bm":"reluCd*bm","relu(cd-10.52)*la":"reluCd*la","hx*bm":"hx*bm","cd*hx*mg":"cd*hx*mg","la*cd*bm":"la*cd*bm"}
b={cmap[n]:c for n,c in zip(nm,coef)}
pv=FF["intercept_alldata"]+sum(b[k]*Tv[k] for k in order)
i_=mean_absolute_error(yv,pv)
# (ii) coef refit, IPM knots
Xii=np.column_stack([Tv[k] for k in order]); ii_=canonCV(Xii,yv)
# (iii) refit knot cd from HSD EBM + coefs
kcd=hsd_knot("L_c_res0_dual",-10.52); Tiii=floor_terms(kcd,Mf); Xiii=np.column_stack([Tiii[k] for k in order]); iii_=canonCV(Xiii,yv)
# forest ref (5 feats)
Xf=np.column_stack([raw["la"],raw["cd"],raw["hx"],raw["bm"],raw["mg"]]); ref_=forestCV(Xf,yv)
LAD["floor"]={"i_verbatim":round(i_,3),"ii_coef":round(ii_,3),"iii_knot":round(iii_,3),"forest_ref":round(ref_,3),"hsd_cd_knot":round(kcd,2),"ipm_cd_knot":-10.52}
print(f"\nFLOOR d_lo ladder (HSD): (i)verbatim={i_:.3f}  (ii)coef={ii_:.3f}  (iii)knot={iii_:.3f}  forest_ref={ref_:.3f}   [IPM-on-IPM was 1.204/forest1.109]")
print(f"   hsd cd-knot {kcd:.2f} vs ipm -10.52")

# ---------- CEILING ----------
yv=H["d_hi"].to_numpy(float)
FC=json.load(open(os.path.join(HERE,"phase1_frozen_ceiling.json"))); Mc=FC["medians"]
rawc={"Ld":Lc("L_Ldiag_min"),"ft":Lc("L_force_tol"),"pd":Lc("L_p_res0_dual"),"la":Lc("L_alpha")}
def ceil_terms(kft,kpd,med):
    C={f:rawc[f]-med[f] for f in rawc}; Hft=np.maximum(rawc["ft"]-kft,0.0); Hpd=np.maximum(rawc["pd"]-kpd,0.0)
    return {"reluFt":Hft,"Ld":C["Ld"],"pd":C["pd"],"reluPd":Hpd,"Ld*la":C["Ld"]*C["la"],
            "Ld*ft":C["Ld"]*C["ft"],"ft*la":C["ft"]*C["la"],"ft*pd":C["ft"]*C["pd"]}
co=FC["coefs"]; cmapc={"relu(ft+9.86)":"reluFt","Ld":"Ld","pd":"pd","relu(pd+10.21)":"reluPd","Ld*la":"Ld*la","Ld*ft":"Ld*ft","ft*la":"ft*la","ft*pd":"ft*pd"}
bc={cmapc[k]:v for k,v in co.items()}; ordc=list(bc)
Tv=ceil_terms(-9.86,-10.21,Mc)
pv=FC["intercept"]+sum(bc[k]*Tv[k] for k in ordc); i_=mean_absolute_error(yv,pv)
Xii=np.column_stack([Tv[k] for k in ordc]); ii_=canonCV(Xii,yv)
kft=hsd_knot("L_force_tol",-9.86); kpd=hsd_knot("L_p_res0_dual",-10.21)
Tiii=ceil_terms(kft,kpd,Mc); Xiii=np.column_stack([Tiii[k] for k in ordc]); iii_=canonCV(Xiii,yv)
Xf=np.column_stack([rawc["Ld"],Lc("L_c_res0_dual"),rawc["ft"],H["cpass"].fillna(0).to_numpy(float)]); ref_=forestCV(Xf,yv)
LAD["ceiling"]={"i_verbatim":round(i_,3),"ii_coef":round(ii_,3),"iii_knot":round(iii_,3),"forest_ref":round(ref_,3),"hsd_ft_knot":round(kft,2),"hsd_pd_knot":round(kpd,2)}
print(f"\nCEILING d_hi ladder (HSD): (i)verbatim={i_:.3f}  (ii)coef={ii_:.3f}  (iii)knot={iii_:.3f}  forest_ref={ref_:.3f}   [IPM-on-IPM was 0.933/forest1.080]")
print(f"   hsd ft-knot {kft:.2f} (ipm -9.86), pd-knot {kpd:.2f} (ipm -10.21)")
json.dump(LAD,open(os.path.join(HERE,"phase2_ladder.json"),"w"),indent=1)
print("\nsaved phase2_ladder.json  DONE")
