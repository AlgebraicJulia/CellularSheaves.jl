#!/usr/bin/env python
"""Stage A — forest reference + cluster-permutation feature discovery (§4).
Grouped 5-fold CV by problem. HistGradientBoosting depth 4. MAE in decades, pooled + by zone.
Run with and without fingerprint meta features (§5). Cluster-level permutation importance."""
import os, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error

D = pd.read_pickle(os.path.join(os.path.dirname(__file__), "dataset.pkl"))
FINGERPRINT = ["m_meta","n_meta","Lng_meta","L_alpha_anchor"]
FEATS = [c for c in D.columns if c not in ("problem","file","d_lo","d_hi","zone")]
CORE = [c for c in FEATS if c not in FINGERPRINT]
groups = D["problem"].values
gkf = GroupKFold(n_splits=5)

def cv_forest(target, feats, seed=0):
    X = D[feats].to_numpy(float); y = D[target].to_numpy(float)
    pred = np.full(len(y), np.nan)
    for tr, te in gkf.split(X, y, groups):
        m = HistGradientBoostingRegressor(max_depth=4, max_iter=300, learning_rate=0.06,
                                          random_state=seed, l2_regularization=1.0)
        m.fit(X[tr], y[tr]); pred[te] = m.predict(X[te])
    mae = mean_absolute_error(y, pred)
    zmae = {z: mean_absolute_error(y[D.zone==z], pred[D.zone==z]) for z in ("below","in","above")}
    return mae, zmae, pred

def report(target, opslice):
    print(f"\n===== {target}  (operational slice: {'+'.join(opslice)}) =====")
    for tag, feats in (("full (w/ fingerprints)", FEATS), ("core (no fingerprints)", CORE)):
        mae, z, _ = cv_forest(target, feats)
        ops = mean_absolute_error(D[D.zone.isin(opslice)][target],
                                  cv_forest(target, feats)[2][D.zone.isin(opslice)])
        print(f"  {tag:24}: pooled MAE={mae:.3f}  | below={z['below']:.2f} in={z['in']:.2f} above={z['above']:.2f}  | op={ops:.3f}")

def clusters(feats, thr=0.8):
    X = D[feats].to_numpy(float)
    C = pd.DataFrame(X, columns=feats).corr().abs().fillna(0).to_numpy()
    parent = list(range(len(feats)))
    def find(a):
        while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
        return a
    for i in range(len(feats)):
        for j in range(i+1,len(feats)):
            if C[i,j]>=thr: parent[find(i)]=find(j)
    cl={}
    for i,f in enumerate(feats): cl.setdefault(find(i),[]).append(f)
    return list(cl.values())

def cluster_perm(target, feats, seed=0):
    X = D[feats].to_numpy(float); y = D[target].to_numpy(float)
    cls = clusters(feats)
    dmg = {}
    for tr, te in gkf.split(X, y, groups):
        m = HistGradientBoostingRegressor(max_depth=4, max_iter=300, learning_rate=0.06,
                                          random_state=seed, l2_regularization=1.0)
        m.fit(X[tr], y[tr])
        base = mean_absolute_error(y[te], m.predict(X[te]))
        rng = np.random.default_rng(seed)
        for cl in cls:
            idx=[feats.index(f) for f in cl]
            Xp=X[te].copy()
            perm=rng.permutation(len(te))
            for k in idx: Xp[:,k]=X[te][perm,k]
            d = mean_absolute_error(y[te], m.predict(Xp)) - base
            key=" + ".join(cl)
            dmg[key]=dmg.get(key,0)+d/5
    return sorted(dmg.items(), key=lambda kv:-kv[1])

if __name__=="__main__":
    print(f"rows={len(D):,} problems={D['problem'].nunique()} features: full={len(FEATS)} core={len(CORE)}")
    report("d_lo", ("below","in"))
    report("d_hi", ("in","above"))
    for target in ("d_lo","d_hi"):
        print(f"\n--- cluster-permutation importance (core feats), {target} ---")
        for key,d in cluster_perm(target, CORE)[:10]:
            print(f"  +{d:.3f}  {key[:90]}")
