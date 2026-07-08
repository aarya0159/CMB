"""Cosmic-variance scaling campaign on the ACTUAL simulator engine.
Sweeps the number of waves N and measures seed-to-seed scatter of the
map rms, testing the predicted sigma/mean = (sqrt(r)/2) N^{-1/2} law."""
import sys, types, json, math
import numpy as np

pg = types.ModuleType("pygame")
class _Rect:
    def __init__(self, *a):
        if len(a) == 4: self.x, self.y, self.w, self.h = a
pg.Rect = _Rect
sys.modules["pygame"] = pg
import cmb_simulator as sim

W, H, TFIX = sim.FIELD_W, sim.FIELD_H, 5.0
PLAN = [(6, 400), (12, 400), (24, 400), (48, 300),
        (96, 200), (192, 120), (384, 80)]

# ---------- theory constant r from the code's own amplitude distributions ----
rng = np.random.default_rng(7)
M = 4_000_000
lam = np.exp(rng.uniform(np.log(18.0), np.log(260.0), M))
env = np.exp(-((np.log(lam / 70.0)) / 0.75) ** 2)
u = rng.uniform(0.6, 1.4, M)
th = rng.uniform(0, 2 * np.pi, M)
s = (env * u) ** 2 * np.cos(th) ** 2
r = s.var() / s.mean() ** 2
print(f"theory: r = {r:.4f}  -> predicted sigma/mean = {math.sqrt(r)/2:.4f} * N^-1/2", flush=True)

# ---------- measurement ----------
results = {}
hists = {}
for N, nseeds in PLAN:
    stds = []
    for k in range(nseeds):
        fld = sim.CMBField(W, H, n_waves=N, seed=20_000 + 37 * k + N)
        stds.append(float(fld.evaluate(TFIX).std()))
    stds = np.array(stds)
    frac = stds.std(ddof=1) / stds.mean()
    err = frac / math.sqrt(2 * (nseeds - 1))
    results[N] = {"nseeds": nseeds, "mean_std": float(stds.mean()),
                  "frac_scatter": float(frac), "frac_err": float(err)}
    if N in (12, 48, 192):
        hists[N] = (3000.0 * 0.0507 * stds).tolist()   # in kelvin at screenshot settings
    print(f"N={N:4d}  seeds={nseeds}  <std>={stds.mean():.4f}  "
          f"sigma/mean={frac:.4f} +- {err:.4f}", flush=True)

# ---------- weighted log-log fit of the exponent ----------
Ns = np.array(list(results))
y = np.array([results[n]["frac_scatter"] for n in Ns])
ey = np.array([results[n]["frac_err"] for n in Ns])
lx, ly, w = np.log(Ns), np.log(y), (y / ey) ** 2   # weights for log-space
Wm = np.sum(w); xb = np.sum(w * lx) / Wm; yb = np.sum(w * ly) / Wm
slope = np.sum(w * (lx - xb) * (ly - yb)) / np.sum(w * (lx - xb) ** 2)
inter = yb - slope * xb
slope_err = 1.0 / math.sqrt(np.sum(w * (lx - xb) ** 2))
amp_fit = math.exp(inter)
print(f"\nFIT: exponent = {slope:.4f} +- {slope_err:.4f}   (prediction: -0.5)")
print(f"FIT: amplitude = {amp_fit:.4f}   (idealized theory: {math.sqrt(r)/2:.4f})")

json.dump({"r": r, "pred_amp": math.sqrt(r)/2, "results": results,
           "fit": {"slope": slope, "slope_err": slope_err, "amp": amp_fit},
           "hists": hists},
          open("variance_results.json", "w"))
print("saved variance_results.json")
