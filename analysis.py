"""Bootstrap uncertainties on the coefficient of variation, a direct goodness-of-fit
test of the predicted law, the two-component refit, and the three-geometry scaling test."""
import json, math, numpy as np

rng = np.random.default_rng(12345)
B = 20000

def boot_cov(v, B=B):
    """Bootstrap the coefficient of variation: value, sd, and the sd/mean covariance."""
    n = len(v)
    idx = rng.integers(0, n, size=(B, n))
    s = v[idx]
    m  = s.mean(axis=1)
    sd = s.std(axis=1, ddof=1)
    cov_b = sd / m
    return (v.std(ddof=1)/v.mean(), cov_b.std(ddof=1),
            float(np.corrcoef(sd, m)[0,1]), float(np.cov(sd, m)[0,1]))

raw = json.load(open("raw_wave_dT.json"))
N = np.array(sorted(int(k) for k in raw))
rows=[]
print("  N     n    CoV      analytic err  bootstrap err   corr(sd,mean)   skew   exc.kurt")
for n_ in N:
    v = np.array(raw[str(n_)])
    c, be, corr, cv = boot_cov(v)
    ana = c/math.sqrt(2*len(v)-1)
    z = (v-v.mean())/v.std(ddof=1)
    skew = float((z**3).mean()); kurt = float((z**4).mean()-3.0)
    rows.append(dict(N=int(n_), n=len(v), cov=c, ana=ana, boot=be, corr=corr,
                     skew=skew, kurt=kurt))
    print("%4d %5d  %.5f   %.5f      %.5f       %+.3f        %+.3f  %+.3f"
          % (n_, len(v), c, ana, be, corr, skew, kurt))

y  = np.array([r["cov"] for r in rows])
eb = np.array([r["boot"] for r in rows])
ea = np.array([r["ana"] for r in rows])

def wls_loglog(N, y, e):
    lx, ly = np.log(N), np.log(y)
    w = (y/e)**2
    Wm=w.sum(); xb=(w*lx).sum()/Wm; yb=(w*ly).sum()/Wm
    sl=(w*(lx-xb)*(ly-yb)).sum()/(w*(lx-xb)**2).sum()
    se=1.0/math.sqrt((w*(lx-xb)**2).sum())
    return sl, se, math.exp(yb-sl*xb)

for lbl, e in (("analytic weights", ea), ("bootstrap weights", eb)):
    sl, se, amp = wls_loglog(N, y, e)
    print("\n%s: exponent %.4f +- %.4f  (%.1f sigma from -0.5), amplitude %.4f"
          % (lbl, sl, se, abs(sl+0.5)/se, amp))

# ---- direct GoF test of the FIXED predicted law 0.806 * N^-1/2 -------------
pred = 0.806 * N**-0.5
chi2_obs = float((((y-pred)/eb)**2).sum())
# null distribution by simulation: resample each point from its bootstrap sd about pred
sim = rng.normal(pred, eb, size=(200000, len(N)))
chi2_null = (((sim-pred)/eb)**2).sum(axis=1)
p = float((chi2_null >= chi2_obs).mean())
print("\nGoF against the fixed 0.806 N^-1/2 law: chi2 = %.1f on %d points, p = %.3g"
      % (chi2_obs, len(N), p))

# ---- two-component refit with bootstrap weights ---------------------------
def fit_two(N, y, e):
    yy = y**2; dy = 2*y*e; w = 1/dy**2
    X = np.column_stack([1.0/N, np.ones_like(N, dtype=float)])
    XtW = X.T*w; cov = np.linalg.inv(XtW@X); beta = cov@(XtW@yy)
    A,Bc = beta; dA,dB = math.sqrt(cov[0,0]), math.sqrt(cov[1,1])
    r = (yy - X@beta)/dy
    return A,dA,Bc,dB,float((r**2).sum()), len(N)-2, X@beta

for lbl, e in (("analytic", ea), ("bootstrap", eb)):
    A,dA,Bc,dB,chi2,dof,mod = fit_two(N,y,e)
    fl=math.sqrt(Bc); dfl=dB/(2*fl)
    print("%s weights: A=%.4f+-%.4f (%.1f sigma from 0.300)  floor=%.4f+-%.4f  chi2=%.2f/%d"
          % (lbl,A,dA,abs(A-0.300)/dA,fl,dfl,chi2,dof))

A,dA,Bc,dB,chi2,dof,mod = fit_two(N,y,eb)
json.dump(dict(rows=rows, A=A, A_err=dA, B=Bc, B_err=dB, floor=math.sqrt(Bc),
               floor_err=dB/(2*math.sqrt(Bc)), chi2=chi2, dof=dof,
               gof_chi2=chi2_obs, gof_p=p,
               model=[float(math.sqrt(m)) for m in mod]),
          open("wave_analysis.json","w"), indent=1)
print("\nwrote wave_analysis.json")
