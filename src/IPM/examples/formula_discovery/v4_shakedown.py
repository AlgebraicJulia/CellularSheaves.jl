#!/usr/bin/env python
"""Logging batch v4 shakedown analysis. Reads the shakedown CSVs and runs:
  ACCEPTANCE: coverage of new cols on native solver; wbase nondegeneracy; tau/kappa finite/sane;
              sigma2max>=sigma2min rowwise; recovery-guard (neg sentinel) rate.
  KMAX DIAG:  within-sweep slope of log sigma2min vs log alpha per kmax {10,20,40}, split low/high alpha.
  u_max PREVIEW: is the in-window ceiling at u_max = log alpha + log sigma2max ~ 15.7 and problem-indep?
  WALL PREVIEW:  rho_w = (r1_w/r0_w)^(1/wbase) well-defined + informative on HSD?"""
import os, glob, numpy as np, pandas as pd
SH=os.environ.get("SHDIR","/Users/richardsamuelson/.claude/jobs/7f05e925/tmp/shakedown")
def load(pat):
    fs=glob.glob(os.path.join(SH,pat)); return {os.path.basename(f).split('.')[0]:pd.read_csv(f) for f in fs}
NEW=["sigma2max","ritz_beta","ritz_gap"]; HSDNEW=["tau","kappa","wbase","wpass"]
main={k:v for k,v in load("*_1e-8_*.csv").items() if not k.startswith("kmaxdiag")}
print("#"*70,"\nACCEPTANCE CHECKS\n"+"#"*70)
for name,d in sorted(main.items()):
    hsd="hsd" in name; print(f"\n--- {name}  ({len(d)} rows, {'HSD' if hsd else 'IPM'}) ---")
    for c in NEW:
        cov=d[c].notna().mean(); print(f"   cov {c:10}= {cov*100:5.1f}%  {'OK' if cov>=0.95 else 'LOW!'}")
    if hsd:
        for c in HSDNEW:
            cov=d[c].notna().mean(); print(f"   cov {c:10}= {cov*100:5.1f}%  {'OK' if cov>=0.95 else 'LOW!'}")
        wb=d["wbase"].dropna()
        print(f"   wbase distribution: min={wb.min():.0f} max={wb.max():.0f} nunique={wb.nunique()} "
              f"{'NONDEGEN OK' if wb.nunique()>1 else 'DEGENERATE (constant!) -- w-channel carries nothing'}")
        for c in ("tau","kappa"):
            v=d[c].dropna(); fin=np.isfinite(v).mean() if len(v) else 0
            print(f"   {c}: finite={fin*100:.0f}% range[{v.min():.2e},{v.max():.2e}] log10-range[{np.log10(np.abs(v)+1e-300).min():.1f},{np.log10(np.abs(v)+1e-300).max():.1f}]")
    # ordering + sentinel
    both=d[["sigma2min","sigma2max"]].dropna(); pos=both[both.sigma2max>0]
    ok=(pos.sigma2max>=pos.sigma2min).mean() if len(pos) else float('nan')
    print(f"   sigma2max>=sigma2min (pos rows): {ok*100:.1f}%  {'OK' if ok>=0.999 else 'VIOLATIONS!'}")
    print(f"   recovery-guard sentinel (s2max<0) rate: {(d.sigma2max<0).mean()*100:.1f}%")

print("\n"+"#"*70,"\nKMAX DIAGNOSTIC (slope of log10 sigma2min vs log10 alpha, per iteration)\n"+"#"*70)
km={k:pd.read_csv(os.path.join(SH,f"kmaxdiag_e02_k{k}.csv")) for k in (10,20,40) if os.path.exists(os.path.join(SH,f"kmaxdiag_e02_k{k}.csv"))}
for k,d in km.items():
    d=d[np.isfinite(d.sigma2min)&(d.sigma2min>0)].copy(); d["la"]=np.log10(d.alpha); d["ls"]=np.log10(d.sigma2min)
    med=d.la.median(); slopes_lo=[]; slopes_hi=[]; slopes_all=[]
    for it,g in d.groupby("iter"):
        if len(g)>=3: slopes_all.append(np.polyfit(g.la,g.ls,1)[0])
        glo=g[g.la<=med]; ghi=g[g.la>med]
        if len(glo)>=3: slopes_lo.append(np.polyfit(glo.la,glo.ls,1)[0])
        if len(ghi)>=3: slopes_hi.append(np.polyfit(ghi.la,ghi.ls,1)[0])
    print(f"  kmax={k}: slope all={np.median(slopes_all):+.3f}  low-α={np.nanmedian(slopes_lo) if slopes_lo else float('nan'):+.3f}  "
          f"high-α={np.nanmedian(slopes_hi) if slopes_hi else float('nan'):+.3f}  (n_iters={len(slopes_all)})")
print("  interpretation: slope FLATTENS with kmax at low-α => truncation bias; STAYS broken at high-α => 1/(1-mu) recovery amplification")

print("\n"+"#"*70,"\nu_max PREVIEW (in-window ceiling at u_max = log alpha + log sigma2max ~ 15.7?)\n"+"#"*70)
for name,d in sorted(main.items()):
    d=d.copy()
    # in-window ceiling per iter = max alpha with state==0; sigma2max ~ const per problem
    ceil=[]
    for it,g in d.groupby("iter"):
        win=g[(g.state==0)&(g.sigma2max>0)]
        if len(win):
            ac=win.alpha.max(); s2=win.loc[win.alpha.idxmax(),"sigma2max"]
            if s2>0: ceil.append(np.log10(ac)+np.log10(s2))
    if ceil: print(f"   {name:22}: u_max@ceiling median={np.median(ceil):.2f}  (pred 15.7)  [n_iters={len(ceil)}]")

print("\n"+"#"*70,"\nWALL PREVIEW (rho_w = (r1_w/r0_w)^(1/wbase) on HSD)\n"+"#"*70)
for name,d in sorted(main.items()):
    if "hsd" not in name: continue
    d=d.dropna(subset=["r1_w","r0_w","wbase"]).copy(); d=d[d.wbase>0]
    if not len(d): print(f"   {name}: no valid w-rows"); continue
    rho=np.power(np.clip(d.r1_w/d.r0_w,1e-300,None),1.0/d.wbase)
    print(f"   {name:22}: rho_w median={np.median(rho):.3e} range[{rho.min():.2e},{rho.max():.2e}] "
          f"nondegen={'yes' if rho.nunique()>5 else 'NO'}  (n={len(d)})")
print("\nDONE")
