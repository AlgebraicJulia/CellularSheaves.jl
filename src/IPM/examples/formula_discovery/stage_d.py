#!/usr/bin/env python
"""Stage D — PySR symbolic regression (§4D). Binary ops +-*/ (interactions PySR can form that the
additive Stage C cannot). Survivor features + counter indicators (SR is poor at integer thresholds).
Diverse subsample (thin per problem so all 96 appear), problem-grouped train/test. Exports the Pareto
front and CV-scores each equation on held-out problems.  Usage: stage_d.py <target> [niter] [maxsize]"""
import sys, os, numpy as np, pandas as pd
from sklearn.metrics import mean_absolute_error

HERE=os.path.dirname(__file__)
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl"))
target = sys.argv[1] if len(sys.argv)>1 else "d_hi"
NITER  = int(sys.argv[2]) if len(sys.argv)>2 else 80
MAXSZ  = int(sys.argv[3]) if len(sys.argv)>3 else 22

# compact survivor set + counter indicators (short names for readable equations)
D["pb3"]=(D["pbase"].fillna(0)>=3).astype(float); D["pb8"]=(D["pbase"].fillna(0)>=8).astype(float)
D["nc5"]=(D["ncraig"].fillna(0)>=5).astype(float)
FEATS=["L_alpha","L_r0_c","L_force_tol","L_floor_tol","L_mu","L_Ldiag_min","L_hdiag_min",
       "L_p_res0_dual","pb3","pb8","nc5"]
short={"L_alpha":"la","L_r0_c":"r0","L_force_tol":"ft","L_floor_tol":"fl","L_mu":"mu",
       "L_Ldiag_min":"Ldmn","L_hdiag_min":"hdmn","L_p_res0_dual":"pdual","pb3":"pb3","pb8":"pb8","nc5":"nc5"}

# thin per problem to ~120 rows so all 96 problems appear; grouped train/test split
rng=np.random.default_rng(0)
keep=[]
for p,idx in D.groupby("problem").groups.items():
    idx=np.array(idx); keep.extend(idx if len(idx)<=120 else rng.choice(idx,120,replace=False))
Dt=D.loc[keep].copy()
probs=Dt["problem"].unique().copy(); rng.shuffle(probs)
ntr=int(0.8*len(probs)); trp=set(probs[:ntr])
tr=Dt[Dt["problem"].isin(trp)]; te=Dt[~Dt["problem"].isin(trp)]
Xtr=tr[FEATS].fillna(tr[FEATS].median()).to_numpy(float); ytr=tr[target].to_numpy(float)
Xte=te[FEATS].fillna(tr[FEATS].median()).to_numpy(float); yte=te[target].to_numpy(float)
print(f"target={target}: thinned {len(Dt)} rows; train {len(tr)}/{len(trp)}p  test {len(te)}/{te['problem'].nunique()}p  niter={NITER} maxsize={MAXSZ}",flush=True)

from pysr import PySRRegressor
model=PySRRegressor(
    niterations=NITER, maxsize=MAXSZ, binary_operators=["+","-","*","/"], unary_operators=[],
    populations=24, model_selection="best", random_state=0, deterministic=True, parallelism="serial",
    variable_names=[short[f] for f in FEATS], progress=False, verbosity=0,
    output_directory=os.path.join(HERE,"pysr_out"),
)
model.fit(Xtr,ytr)
print("\n=== Pareto front (complexity | held-out grouped MAE | equation) ===",flush=True)
eqs=model.equations_
for _,row in eqs.iterrows():
    try: pred=model.predict(Xte, index=int(row["complexity"]))  # predict with this eq
    except Exception: continue
    mae=mean_absolute_error(yte,pred)
    print(f"  c={int(row['complexity']):2d}  MAE={mae:.3f}   {row['equation']}",flush=True)
print("\nBEST (model_selection):", model.sympy(),flush=True)
print("DONE",flush=True)
