import json, math, glob, numpy as np
rng = np.random.default_rng(777); B = 20000
wa = json.load(open("wave_analysis.json")); A, dA = wa["A"], wa["A_err"]
NW = 384; AN = A/NW

cells = {}
for f in glob.glob("cell_*.json"):
    d = json.load(open(f)); cells[(d["W"], d["H"])] = np.array(d["dT"])

GEOM = {"both":   [(192,120),(384,240),(768,480)],
        "width":  [(192,240),(384,240),(768,240)],
        "height": [(384,120),(384,240),(384,480)]}
PRED = {"both": -1.00, "width": -0.50, "height": -0.50}
XLAB = {"both": "L (both dims)", "width": "W (H fixed 240)", "height": "H (W fixed 384)"}

def floor_boot(v):
    """Bootstrap the floor sqrt(CoV^2 - A/N), propagating the uncertainty in A too."""
    n = len(v); idx = rng.integers(0, n, size=(B, n)); s = v[idx]
    cov_b = s.std(axis=1, ddof=1)/s.mean(axis=1)
    A_b = rng.normal(A, dA, B)
    var = cov_b**2 - A_b/NW
    var = var[var > 0]
    f = np.sqrt(var)
    c = v.std(ddof=1)/v.mean()
    return math.sqrt(c**2 - AN), float(f.std(ddof=1)), c

out = {}
for g, sizes in GEOM.items():
    xs, fs, es, cs = [], [], [], []
    for (W,H) in sizes:
        v = cells[(W,H)]
        f, ef, c = floor_boot(v)
        x = W if g in ("both","width") else H
        xs.append(x); fs.append(f); es.append(ef); cs.append(c)
    xs=np.array(xs,float); fs=np.array(fs); es=np.array(es)
    lx, lf = np.log(xs), np.log(fs); w=(fs/es)**2
    Wm=w.sum(); xb=(w*lx).sum()/Wm; yb=(w*lf).sum()/Wm
    sl=(w*(lx-xb)*(lf-yb)).sum()/(w*(lx-xb)**2).sum()
    se=1.0/math.sqrt((w*(lx-xb)**2).sum())
    out[g]=dict(x=xs.tolist(), cov=cs, floor=fs.tolist(), floor_err=es.tolist(),
                slope=sl, slope_err=se, predicted=PRED[g],
                sigma_from_pred=abs(sl-PRED[g])/se,
                sizes=[f"{a}x{b}" for a,b in sizes])
    print("\n%-7s  %s   predicted slope %.2f" % (g, XLAB[g], PRED[g]))
    for (W,H),c,f,e in zip(sizes,cs,fs,es):
        print("   %4dx%-4d  n=%3d  CoV=%.5f  floor=%.5f +- %.5f" % (W,H,len(cells[(W,H)]),c,f,e))
    r1=fs[0]/fs[1]; r2=fs[1]/fs[2]
    d1=r1*math.hypot(es[0]/fs[0],es[1]/fs[1]); d2=r2*math.hypot(es[1]/fs[1],es[2]/fs[2])
    exp_r = 2.0 if PRED[g]==-1.0 else math.sqrt(2.0)
    print("   ratios %.3f +- %.3f and %.3f +- %.3f   (expected %.3f)" % (r1,d1,r2,d2,exp_r))
    print("   fitted slope %.3f +- %.3f  -> %.1f sigma from %.2f, %.1f sigma from %.2f"
          % (sl,se,abs(sl-PRED[g])/se,PRED[g],abs(sl-(-1.0 if PRED[g]==-0.5 else -0.5))/se,
             (-1.0 if PRED[g]==-0.5 else -0.5)))
json.dump(out, open("width_analysis.json","w"), indent=1)
print("\nwrote width_analysis.json")
