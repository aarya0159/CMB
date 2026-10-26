"""Field evaluation with explicit control of each random draw, for the ablation."""
import math
import numpy as np
PEAK = 70.0

def draws(n, seed):
    """The five independent random draws, in the simulator's exact RNG order."""
    rng = np.random.default_rng(seed)
    lam = np.exp(rng.uniform(np.log(18.0), np.log(260.0), n))
    u   = rng.uniform(0.6, 1.4, n)          # amplitude random factor
    th  = rng.uniform(0.0, 2.0*np.pi, n)    # direction
    phi = rng.uniform(0.0, 2.0*np.pi, n)    # spatial phase
    psi = rng.uniform(0.0, 2.0*np.pi, n)    # temporal phase
    return lam, u, th, phi, psi

def std_from(lam, u, th, phi, psi, W, H, t=5.0, block=32):
    n = len(lam)
    amp = np.exp(-((np.log(lam/PEAK))/0.75)**2) * u
    k = 2.0*np.pi/lam
    kx, ky = k*np.cos(th), k*np.sin(th)
    omega = (2.0*np.pi/12.0)*(PEAK/lam)
    coeff = (amp*np.cos(omega*t + psi)).astype(np.float32)
    norm = 1.0/(2.8*math.sqrt(float(np.sum(amp**2))/4.0))
    ygrid, xgrid = np.mgrid[0:H, 0:W].astype(np.float32)
    acc = np.zeros((H, W), dtype=np.float64)
    for s in range(0, n, block):
        e = min(s+block, n)
        arg = (kx[s:e,None,None]*xgrid + ky[s:e,None,None]*ygrid + phi[s:e,None,None])
        acc += np.tensordot(coeff[s:e], np.sin(arg, dtype=np.float32), axes=1)
        del arg
    acc *= norm
    return float(acc.std())

def field_std(W, H, n, seed, t=5.0, block=32):
    return std_from(*draws(n, seed), W=W, H=H, t=t, block=block)
