"""Trial runner: executes the ACTUAL simulator engine (cmb_simulator.py,
unmodified) and records telemetry exactly as the app computes it
(lines 311-324 of the source). Only the pygame *graphics* import is stubbed;
all physics code paths are the author's own."""
import sys, types, math, json
import numpy as np

# stub pygame (graphics only; never used by CMBField / age table / telemetry)
pg = types.ModuleType("pygame")
class _Rect:
    def __init__(self, *a):
        if len(a) == 4: self.x, self.y, self.w, self.h = a
        elif len(a) == 1: self.x, self.y, self.w, self.h = a[0]
pg.Rect = _Rect
sys.modules["pygame"] = pg

import cmb_simulator as sim

NU = sim.NU_PEAK_PER_K   # 5.878925757e10 Hz/K
WB = sim.WIEN_B          # 2.897771955e-3 m K
T0 = sim.T_TODAY

def telemetry(t_base, amp, field, t):
    """Exactly the app's telemetry block (main loop lines 315-324)."""
    f = field.evaluate(t)
    t_map = t_base * (1.0 + amp * f)
    return dict(
        t_avg=float(t_map.mean()),
        dt_rms=float(t_map.std()),
        nu_peak=NU * float(t_map.mean()),
        lam_peak=WB / float(t_map.mean()),
        redshift=t_base / T0 - 1.0,
        age_yr=sim.age_at_redshift(t_base / T0 - 1.0),
    )

W, H = sim.FIELD_W, sim.FIELD_H   # 384 x 240, the app's real grid
out = {}

# ---------- Experiment 1: Wien linearity sweep ----------
fld = sim.CMBField(W, H, seed=11)
rows = []
for T in [3000.0, 1500.0, 800.0, 400.0, 200.0, 100.0, 30.0, 10.0, 2.72548]:
    m = telemetry(T, 0.02, fld, 5.0)
    rows.append((T, m["t_avg"], m["nu_peak"]))
slopes = [nu / ta for (_, ta, nu) in rows]
out["E1"] = {"rows": rows, "slope_GHz_per_K": [s/1e9 for s in slopes]}

# ---------- Experiment 2: amplitude sweep at 3000 K ----------
rows2 = []
for amp in [0.3, 0.1, 0.0507, 0.02, 0.01, 0.001, 1e-4]:
    m = telemetry(3000.0, amp, fld, 5.0)
    gain = math.sqrt(amp / 0.02)          # the app's visual contrast gain
    rows2.append((amp, m["dt_rms"], gain))
out["E2"] = rows2

# ---------- Experiment 3: cooling / redshift ride ----------
rows3 = []
for T in [3000.0, 1500.0, 300.0, 30.0, 2.72548]:
    m = telemetry(T, 0.02, fld, 5.0)
    rows3.append((T, m["redshift"], m["age_yr"], m["nu_peak"], m["lam_peak"],
                  m["nu_peak"] / (1.0 + m["redshift"])))
out["E3"] = rows3

# ---------- Experiment 4: dT rms pulsing at screenshot settings ----------
samples = []
for t in np.arange(0.0, 24.0, 0.25):      # two full ~12 s standing-wave periods
    m = telemetry(3000.0, 0.0507, fld, float(t))
    samples.append(m["dt_rms"])
samples = np.array(samples)
out["E4"] = {"amp": 0.0507, "T": 3000.0,
             "dt_mean": float(samples.mean()),
             "dt_min": float(samples.min()), "dt_max": float(samples.max()),
             "spread_THz_mean": 58.7893e9 * samples.mean() / 1e12,
             "pct_of_peak": 100 * samples.mean() * 0.02 / 0.02 * 0.0507 / 0.0507 * (58.7893e9*samples.mean()) / (NU*3000.0)}

# ---------- Experiment 5: cosmic variance across seeds ----------
vals = []
for s in range(1, 13):
    f2 = sim.CMBField(W, H, seed=s)
    m = telemetry(3000.0, 0.0507, f2, 5.0)
    vals.append(m["dt_rms"])
vals = np.array(vals)
out["E5"] = {"n": len(vals), "mean": float(vals.mean()),
             "std": float(vals.std(ddof=1)),
             "min": float(vals.min()), "max": float(vals.max())}

# ---------- Validation table at recombination default ----------
m = telemetry(3000.0, 0.0507, fld, 5.0)
out["VAL"] = {**m,
  "pred_nu": NU*m["t_avg"], "pred_lam": WB/m["t_avg"],
  "z_pred": 3000.0/T0-1.0}

# ---------- time-averaged field std (the normalization the map realizes) ----------
stds = [float(fld.evaluate(float(t)).std()) for t in np.arange(0, 24, 0.25)]
out["FIELD_STD"] = {"mean": float(np.mean(stds)), "min": float(np.min(stds)),
                    "max": float(np.max(stds))}

print(json.dumps(out, indent=1, default=str))
