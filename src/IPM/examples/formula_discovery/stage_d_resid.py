#!/usr/bin/env python
"""Stage D residual track (§ brief item 2). f_C = Stage-C 8-term formula; PySR fits r = y - f_C so it
only hunts the missing interaction. Assisted inputs = survivors + physics composites + counter indicators.
Export FULL Pareto front; CV-score (f_C + beta*expr) under grouped 5-fold (refit outer coef per fold).
X04 policy: report CV MAE with AND without the X04* family. Usage: stage_d_resid.py <target> [niter]"""
import sys, os, numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error
sys.path.insert(0, os.path.dirname(__file__))
import stage_c as SC

HERE=os.path.dirname(__file__); D=SC.D
target = sys.argv[1] if len(sys.argv)>1 else "d_hi"
NITER  = int(sys.argv[2]) if len(sys.argv)>2 else 150
groups = D["problem"].values; y = D[target].to_numpy(float)

# --- reconstruct Stage-C 8-term f_C on all rows ---
Lib = SC.library(D)
mae_c, z_c, b0, terms, ns = SC.fit(target, nterms=8)
fC = np.full(len(D), b0, float)
for name,co in terms: fC += co*Lib[name].to_numpy(float)
r = y - fC
print(f"target={target}: Stage-C f_C 8-term CV-MAE={mae_c:.3f}; residual std={r.std():.2f}",flush=True)

# --- assisted inputs for PySR (short names) ---
def col(c): return D[c].fillna(D[c].median()).to_numpy(float)
FEAT = {
 "la":col("L_alpha"), "r0":col("L_r0_c"), "ft":col("L_force_tol"), "fl":col("L_floor_tol"),
 "mu":col("L_mu"), "Ldmn":col("L_Ldiag_min"), "Ldmx":col("L_Ldiag_max"),
 "hdmn":col("L_hdiag_min"), "hdmx":col("L_hdiag_max"), "pdual":col("L_p_res0_dual"),
 "s2":col("L_sigma2min"),
 "margin":col("L_r0_c")-np.maximum(col("L_force_tol"),col("L_floor_tol")),
 "ceilm":col("L_Ldiag_min")+col("L_alpha"),
 "pb3":(D["pbase"].fillna(0)>=3).to_numpy(float),
}
names=list(FEAT); Xall=np.column_stack([FEAT[n] for n in names])

# --- thin per problem, grouped train for PySR ---
rng=np.random.default_rng(0); keep=[]
for p,idx in D.groupby("problem").groups.items():
    idx=np.array(idx); keep.extend(idx if len(idx)<=120 else rng.choice(idx,120,replace=False))
keep=np.array(keep); probs=D.loc[keep,"problem"].unique().copy(); rng.shuffle(probs)
trp=set(probs[:int(0.85*len(probs))])
trm=np.array([g in trp for g in D.loc[keep,"problem"]])
tr_idx=keep[trm]
Xtr=Xall[tr_idx]; rtr=r[tr_idx]

from pysr import PySRRegressor
model=PySRRegressor(niterations=NITER, maxsize=24, binary_operators=["+","-","*","/"],
    unary_operators=[], populations=24, parallelism="serial", random_state=0, deterministic=True,
    progress=False, verbosity=0, output_directory=os.path.join(HERE,"pysr_out"))
model.fit(Xtr, rtr, variable_names=names)

# --- CV-score (f_C + beta*expr) under grouped 5-fold; refit beta per fold ---
gkf=GroupKFold(5); zone=D["zone"].to_numpy(); notX04=~D["problem"].astype(str).str.startswith("X04").to_numpy()
def cvscore(expr_pred):
    pred=np.full(len(y),np.nan)
    F=np.column_stack([fC, expr_pred])
    for a,b in gkf.split(F,y,groups):
        pred[b]=LinearRegression().fit(F[a],y[a]).predict(F[b])
    return (mean_absolute_error(y,pred),
            mean_absolute_error(y[notX04],pred[notX04]),
            {zz:mean_absolute_error(y[zone==zz],pred[zone==zz]) for zz in ("in","above")})
print("\n=== Pareto front: complexity | CV-MAE(C+expr) | noX04 | expr ===",flush=True)
best=None
for _,row in model.equations_.iterrows():
    c=int(row["complexity"])
    try: ep=model.predict(Xall, index=c)
    except Exception: continue
    mae,maeN,z=cvscore(ep)
    tag="  <<<" if (best is None or mae<best[0]) else ""
    if best is None or mae<best[0]: best=(mae,c,row["equation"])
    print(f"  c={c:2d}  {mae:.3f}  (noX04 {maeN:.3f})  [in {z['in']:.2f} above {z['above']:.2f}]  {row['equation']}{tag}",flush=True)
print(f"\nBEST by CV: c={best[1]} MAE={best[0]:.3f}  |  Stage-C alone was {mae_c:.3f}",flush=True)
print("DONE",flush=True)
