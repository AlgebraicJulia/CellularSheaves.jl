#!/usr/bin/env python
"""Fit + SAVE an EBM per target on a grouped subsample, then analyze.
Usage: ebm_fit_save.py <n_rows> <outer_bags> <max_rounds(0=default)> <tag>
Grouped subsample: whole problems until ~n_rows -> train; the rest -> held-out test (clean, no
group overlap). Pickles each fit to ebm_ckpt/<tag>_<target>.pkl the instant it finishes (resumable)."""
import sys, time, os, pickle, numpy as np, pandas as pd
from interpret.glassbox import ExplainableBoostingRegressor
from sklearn.metrics import mean_absolute_error

HERE=os.path.dirname(__file__); CK=os.path.join(HERE,"ebm_ckpt"); os.makedirs(CK,exist_ok=True)
D=pd.read_pickle(os.path.join(HERE,"dataset.pkl"))
FP=["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
CORE=[c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")+tuple(FP)]

N=int(sys.argv[1]); BAGS=int(sys.argv[2]); ROUNDS=int(sys.argv[3]); TAG=sys.argv[4]

# grouped subsample: shuffle problems, take whole problems until >= N rows
rng=np.random.default_rng(0); probs=D["problem"].unique().copy(); rng.shuffle(probs)
tr_probs=[]; nrows=0
for p in probs:
    tr_probs.append(p); nrows+=(D["problem"]==p).sum()
    if nrows>=N: break
tr_mask=D["problem"].isin(tr_probs).to_numpy(); te_mask=~tr_mask
Xtr=D[tr_mask][CORE].to_numpy(float); Xte=D[te_mask][CORE].to_numpy(float)
print(f"tag={TAG}: train {tr_mask.sum()} rows / {len(tr_probs)} problems  |  test {te_mask.sum()} rows / "
      f"{D[te_mask]['problem'].nunique()} problems  | bags={BAGS} rounds={ROUNDS or 'default'}",flush=True)

def ebm():
    kw=dict(interactions=0, outer_bags=BAGS, random_state=0)
    if ROUNDS>0: kw["max_rounds"]=ROUNDS
    return ExplainableBoostingRegressor(**kw)

def _isnum(v):
    try: float(v); return True
    except: return False

for target in ("d_lo","d_hi"):
    ckf=os.path.join(CK,f"{TAG}_{target}.pkl")
    if os.path.exists(ckf):
        rec=pickle.load(open(ckf,"rb")); m=rec["model"]
    else:
        y=D[tr_mask][target].to_numpy(float)
        t=time.time(); m=ebm().fit(Xtr,y); dt=time.time()-t
        rec={"tag":TAG,"target":target,"model":m,"tr_probs":list(tr_probs),"secs":dt}
        pickle.dump(rec,open(ckf,"wb")); print(f"  [fit {dt:.1f}s -> {os.path.basename(ckf)}]",flush=True)
    # held-out MAE
    yte=D[te_mask][target].to_numpy(float); pred=m.predict(Xte); zt=D[te_mask]["zone"].to_numpy()
    mae=mean_absolute_error(yte,pred)
    z={zz:mean_absolute_error(yte[zt==zz],pred[zt==zz]) for zz in ("below","in","above") if (zt==zz).any()}
    print(f"  {target}: held-out MAE={mae:.3f}  [below={z.get('below',float('nan')):.2f} "
          f"in={z.get('in',float('nan')):.2f} above={z.get('above',float('nan')):.2f}]",flush=True)
    # shapes + importances
    g=m.explain_global(); names=list(m.term_names_)
    imp=sorted(zip(CORE,m.term_importances()),key=lambda kv:-kv[1])
    print(f"    top terms: "+", ".join(f"{n}({i:.2f})" for n,i in imp[:6]),flush=True)
    for k in ("L_Ldiag_min","L_force_tol","L_mu","L_r0_c","L_p_res0_dual","L_hdiag_min","pbase","ncraig"):
        if k not in names: continue
        d=g.data(names.index(k)); xs=list(d.get("names") or []); ys=list(d.get("scores") or [])
        nn=min(len(xs),len(ys))
        if not nn: continue
        sel=np.linspace(0,nn-1,min(7,nn)).astype(int)
        fx=lambda v:(f"{float(v):.2g}" if _isnum(v) else str(v))
        print(f"    shape {k}: "+", ".join(f"{fx(xs[j])}:{ys[j]:+.2f}" for j in sel),flush=True)
print("DONE",flush=True)
