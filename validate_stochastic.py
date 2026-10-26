"""Validation of the stochastic machinery, which is what Major 8 asks for.
No scipy: Bessel J0 and the KS p-value are computed from first principles."""
import json, math, numpy as np, core

rng = np.random.default_rng(2024)
res = {}

def ks_uniform(x, lo, hi):
    """KS statistic and asymptotic p against Uniform(lo,hi)."""
    u = np.sort((np.asarray(x)-lo)/(hi-lo)); n = len(u)
    i = np.arange(1, n+1)
    d = max(np.max(i/n - u), np.max(u - (i-1)/n))
    lam = (math.sqrt(n) + 0.12 + 0.11/math.sqrt(n)) * d
    p = 2*sum((-1)**(k-1)*math.exp(-2*k*k*lam*lam) for k in range(1, 101))
    return float(d), float(min(max(p, 0.0), 1.0))

# ---- 1. the five draws, pooled over 2000 realisations at N=48 --------------
NR, NW = 2000, 48
L=[];U=[];T=[];P=[];S=[]
for k in range(NR):
    l,u,t,p,s = core.draws(NW, 20000+37*k+NW)
    L.append(l);U.append(u);T.append(t);P.append(p);S.append(s)
L=np.concatenate(L);U=np.concatenate(U);T=np.concatenate(T);P=np.concatenate(P);S=np.concatenate(S)
d,pv = ks_uniform(np.log(L), math.log(18.0), math.log(260.0)); res["wavelength_log_uniform"]={"D":d,"p":pv,"n":len(L)}
print("log-wavelength  ~ U(ln18, ln260)   D=%.4f  p=%.3f  n=%d" % (d,pv,len(L)))
for nm,arr,lo,hi in [("amplitude factor",U,0.6,1.4),("direction",T,0,2*math.pi),
                     ("spatial phase",P,0,2*math.pi),("temporal phase",S,0,2*math.pi)]:
    d,pv = ks_uniform(arr, lo, hi); res[nm.replace(" ","_")]={"D":d,"p":pv,"n":len(arr)}
    print("%-16s~ U(%.2f, %.2f)%s D=%.4f  p=%.3f  n=%d" % (nm,lo,hi," "*6,d,pv,len(arr)))

# ---- 2. the normalisation and the A coefficient it should produce ----------
M = 4_000_000
g = np.random.default_rng(7)
lam = np.exp(g.uniform(np.log(18.0), np.log(260.0), M))
env = np.exp(-((np.log(lam/70.0))/0.75)**2)
u   = g.uniform(0.6, 1.4, M)
a2  = (env*u)**2
q_mc = float((a2**2).mean()/a2.mean()**2)
# Gauss-Legendre in ln(lambda), analytic over u
xs, ws = np.polynomial.legendre.leggauss(400)
lo, hi = math.log(18.0), math.log(260.0)
ll = 0.5*(hi-lo)*xs + 0.5*(hi+lo); wq = 0.5*(hi-lo)*ws/(hi-lo)
e2 = np.exp(-2*((ll-math.log(70.0))/0.75)**2)
Eu2, Eu4 = (1.4**3-0.6**3)/(3*0.8), (1.4**5-0.6**5)/(5*0.8)
Ea2 = float((wq*e2).sum())*Eu2
Ea4 = float((wq*e2**2).sum())*Eu4
q_gl = Ea4/Ea2**2
res["q"] = {"monte_carlo": q_mc, "gauss_legendre": q_gl,
            "A_pred": q_gl/8.0, "r": 1.5*q_gl-1.0, "naive_coeff": math.sqrt(1.5*q_gl-1.0)/2}
print("\nq = E[a^4]/E[a^2]^2 : Monte Carlo %.4f   Gauss-Legendre %.4f" % (q_mc, q_gl))
print("  -> A = q/8 = %.4f ;  r = 3q/2 - 1 = %.4f ;  naive sqrt(r)/2 = %.4f"
      % (q_gl/8, 1.5*q_gl-1, math.sqrt(1.5*q_gl-1)/2))
wa = json.load(open("wave_analysis.json"))
print("  measured A (bootstrap weights) = %.4f +- %.4f  -> %.1f sigma from prediction"
      % (wa["A"], wa["A_err"], abs(wa["A"]-q_gl/8)/wa["A_err"]))

# ---- 3. two-point covariance, analytic vs empirical ------------------------
def J0(z):
    th = np.linspace(0, math.pi, 2001)
    return np.trapezoid(np.cos(np.outer(np.asarray(z), np.sin(th))), th, axis=1)/math.pi
W=H=256; NWv=96; t=5.0; NRc=120
lags = np.arange(0, 41)
emp = np.zeros(len(lags)); acc_var=0.0
for k in range(NRc):
    l,u,th_,p,s = core.draws(NWv, 5000+k)
    amp = np.exp(-((np.log(l/70.0))/0.75)**2)*u
    kk = 2*np.pi/l; kx,ky = kk*np.cos(th_), kk*np.sin(th_)
    om = (2*np.pi/12.0)*(70.0/l)
    c = amp*np.cos(om*t+s); norm = 1.0/(2.8*math.sqrt(float((amp**2).sum())/4.0))
    yg,xg = np.mgrid[0:H,0:W].astype(np.float32)
    f = np.zeros((H,W))
    for i in range(NWv):
        f += c[i]*np.sin(kx[i]*xg + ky[i]*yg + p[i])
    f *= norm; f -= f.mean()
    acc_var += float((f*f).mean())
    for j,L_ in enumerate(lags):
        emp[j] += float((f[:, :W-L_]*f[:, L_:]).mean()) if L_>0 else float((f*f).mean())
emp/=NRc
# analytic: C(r) = norm^2 * sum_j c_j^2/2 * J0(k_j r), ensemble-averaged over draws
ana = np.zeros(len(lags)); 
for k in range(NRc):
    l,u,th_,p,s = core.draws(NWv, 5000+k)
    amp = np.exp(-((np.log(l/70.0))/0.75)**2)*u
    kk = 2*np.pi/l; om=(2*np.pi/12.0)*(70.0/l)
    c = amp*np.cos(om*t+s); norm = 1.0/(2.8*math.sqrt(float((amp**2).sum())/4.0))
    ana += norm**2*0.5*(c[:,None]**2*J0(np.outer(kk,lags).ravel()).reshape(len(kk),len(lags))).sum(axis=0)
ana/=NRc
rel = np.abs(emp-ana)/ana[0]
res["two_point_covariance"] = {"lags": lags.tolist(), "empirical": emp.tolist(),
                               "analytic": ana.tolist(), "max_abs_dev_rel_to_C0": float(rel.max())}
print("\ntwo-point covariance, %d realisations on a %dx%d map, lags 0-40 px" % (NRc,W,H))
for j in (0,2,5,10,20,40):
    print("   lag %2d  empirical %+.5f   analytic %+.5f   diff/C(0) %+.4f" % (lags[j],emp[j],ana[j],(emp[j]-ana[j])/ana[0]))
print("   worst deviation over all lags, relative to C(0): %.4f" % rel.max())

json.dump(res, open("validation_stochastic.json","w"), indent=1)
print("\nwrote validation_stochastic.json")
