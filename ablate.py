import json, sys, time, numpy as np, core

def cov_for(N, ns, vary, W=384, H=240, ref=0):
    """vary: set of draw names that follow the seed; the rest are frozen at `ref`."""
    L0,U0,T0,P0,S0 = core.draws(N, ref)
    out=[]
    for k in range(ns):
        L,U,T,P,S = core.draws(N, 20000+37*k+N)
        lam = L if 'lam' in vary else L0
        u   = U if 'u'   in vary else U0
        th  = T if 'th'  in vary else T0
        phi = P if 'phi' in vary else P0
        psi = S if 'psi' in vary else S0
        out.append(core.std_from(lam,u,th,phi,psi,W,H))
    v=np.array(out)
    return v.std(ddof=1)/v.mean(), v.mean()

N  = int(sys.argv[1]); NS = int(sys.argv[2])
ALL = {'lam','u','th','phi','psi'}
conds = [("all vary", ALL)] + [("only "+x, {x}) for x in ["lam","u","th","phi","psi"]]
res={}
t0=time.time()
for name, vary in conds:
    c,m = cov_for(N, NS, vary)
    res[name]={"cov":c,"mean":m}
    print("N=%-4d %-10s cov=%.5f  mean=%.5f  [%.0fs]"%(N,name,c,m,time.time()-t0), flush=True)
json.dump(res, open("ablation_N%d.json"%N,"w"), indent=1)
