# Finite Map Variance in a Simulated Cosmic Microwave Background

Code and raw outputs supporting every number, table and figure in *Finite Map Variance in a
Simulated Cosmic Microwave Background: Testing a Random-Walk Scaling Law Across 1,900 Simulated
Universes*.

Aarya Pandey and Vijaya Geet Sudhams Bathamsetty, Global Indian International School, Singapore.

---

## Contents

| File | Produces |
|---|---|
| `cmb_simulator.py` | The simulator engine, unmodified. Every script below imports it rather than reimplementing its physics. |
| `index.html` | The browser version of the simulator. |
| `core.py` | A memory-lean field evaluator, written so the largest maps are affordable. Agrees with `cmb_simulator.CMBField` to 3 parts in 10⁷. |
| `variance_campaign.py` | Table 1 and the measured exponent. Writes `variance_results.json`. |
| `width_campaign.py` | The original map-width campaign. Writes `width_results.json`. |
| `cell.py` | One map size of the extended map-size campaign, resumable. Writes `cell_<W>x<H>_<n>.json`. |
| `analysis.py` | Bootstrap uncertainties, the goodness-of-fit test, the two-component fit. Writes `wave_analysis.json`. |
| `width_analysis.py` | The three-geometry test. Writes `width_analysis.json`. |
| `ablate.py` | Which random draw drives the scatter. Writes `ablation_N<N>.json`. |
| `validate_stochastic.py` | Distribution, covariance and normalisation tests. Writes `validation_stochastic.json`. |
| `validation_trials.py` | Table 2, the telemetry comparison. |
| `fit_two_component.py` | The earlier two-component fit and floor extraction. Writes `fit_results.json`. |
| `make_figures.py` | Figure 1, all three panels. Writes `figure1_revised.png` and `.pdf`. |
| `raw_wave_dT.json`, `cell_*.json` | The individual ΔT value from every realisation in both campaigns. |

## Running

Python 3 with NumPy and Matplotlib. SciPy is not required: the Bessel function and the
Kolmogorov-Smirnov statistic used in the validation tests are implemented directly.
`pygame` is not required either, because the scripts stub the graphics import so the
simulator's physics can be driven without a display.

```
python variance_campaign.py       # Table 1              -> variance_results.json
python cell.py 192 120 400        # one map size, repeat for each size below
python analysis.py                # bootstrap + GoF      -> wave_analysis.json
python width_analysis.py          # three geometries     -> width_analysis.json
python ablate.py 12 300           # and  ablate.py 384 150
python validate_stochastic.py     # stochastic tests     -> validation_stochastic.json
python make_figures.py            # Figure 1             -> figure1_revised.png
```

The seven map sizes are 192×120, 192×240, 384×120, 384×240, 384×480, 768×240 and 768×480, at
400, 400, 400, 300, 200, 200 and 200 realisations respectively. `cell.py` takes a time budget as
its fourth argument and resumes where it left off, so a large map can be run in several passes.

`analysis.py` must run before `width_analysis.py`, which needs the fitted A.

## Reproducibility

Every realisation is generated from a deterministic integer seed, so these numbers come out
identically on any machine.

- **Wave-count campaign.** Seed `20000 + 37k + N` for the kth realisation at wave count N. Seven
  wave counts, 1,900 realisations, on the 384×240 grid, each field evaluated at a fixed instant
  t = 5.0 s in the standing-wave cycle.
- **Map-size campaign.** Seed `k + 1`, the same sequence at every map size, wave count held at
  N = 384, same fixed instant.
- **Validation.** 480 comparisons: five baseline temperatures from 3,000 K to 2.72548 K, twelve
  seeds and eight instants in the cycle.

## Provenance of the map-size campaign

The three-width campaign reported in the submitted manuscript was run as an inline command rather
than a saved script, and its seed list was not retained. It cannot be reproduced exactly, and a
search across roughly 26.7 million candidate seed sequences, tested against the three surviving
double-precision outputs, did not recover it.

`cell.py` and `width_analysis.py` regenerate that campaign from the documented seed rule above and
extend it from three map sizes to seven. Only the seed sequence is new. The regenerated campaign
also uses more realisations than the original, because at the largest map size with 30 realisations
the measured scatter sits only slightly above `sqrt(A/N)`, so the subtraction removes most of it
and the surviving floor is governed by which seeds happened to be drawn.

The wave-count campaign survived as a saved script and reproduces Table 1, the measured exponent,
the two-component fit, the floor and the chi-squared exactly.

---

## What the revision found

### 1. The deposited code reproduces the submitted numbers exactly

| N | submitted | re-run |
|---|---|---|
| 6 | 0.233 | 0.23326 |
| 12 | 0.185 | 0.18521 |
| 24 | 0.128 | 0.12772 |
| 48 | 0.097 | 0.09722 |
| 96 | 0.080 | 0.08004 |
| 192 | 0.070 | 0.07029 |
| 384 | 0.073 | 0.07341 |

The map-width campaign reproduces as well: 0.11150, 0.06367 and 0.03858 against the submitted
0.111501, 0.063670 and 0.038580.

### 2. Bootstrap uncertainties

The per-seed ΔT values are retained and the coefficient of variation bootstrapped at each N, with
20,000 resamples.

| N | n | CoV | analytic err | bootstrap err | corr(sd, mean) | skew |
|---|---|---|---|---|---|---|
| 6 | 400 | 0.23326 | 0.00825 | 0.00997 | −0.414 | −0.582 |
| 12 | 400 | 0.18521 | 0.00655 | 0.00743 | −0.284 | −0.410 |
| 24 | 400 | 0.12772 | 0.00452 | 0.00461 | −0.229 | −0.313 |
| 48 | 300 | 0.09722 | 0.00397 | 0.00378 | −0.053 | −0.068 |
| 96 | 200 | 0.08004 | 0.00401 | 0.00435 | +0.204 | +0.309 |
| 192 | 120 | 0.07029 | 0.00455 | 0.00425 | +0.454 | +0.619 |
| 384 | 80 | 0.07341 | 0.00582 | 0.00502 | −0.339 | −0.385 |

Two things the analytic formula misses and the bootstrap does not. The ΔT distributions are not
normal at small N, with skew −0.58 at N = 6 falling to near zero by N = 48. And the numerator and
denominator are correlated, with a sign that is not even fixed: corr(sd, mean) runs from −0.41 at
N = 6 to +0.45 at N = 192.

With bootstrap weights the measured exponent is −0.3288 ± 0.0143, which is 12.0σ from the predicted
−0.5 rather than 11.2σ. The two-component fit gives A = 0.3156 ± 0.0188 against the predicted
0.300, a difference of 0.8σ, with floor 0.0589 ± 0.0035 and χ² = 6.11 on 5 degrees of freedom.

### 3. Direct test of the predicted law

Fitting a free exponent and quoting its distance from −0.5 is the weaker test, because the paper's
own case is that the data bend and a single power law is the wrong shape. Tested directly against
the fixed 0.806·N^(−1/2) law, χ² = 272.1 on seven points. The null distribution was built by
simulation, 200,000 draws, with zero occurrences at or above the observed value, so p < 5 × 10⁻⁶.

### 4. Three geometries

Same seed rule, fixed N = 384.

| geometry | maps | predicted slope | measured slope |
|---|---|---|---|
| both dimensions | 192×120, 384×240, 768×480 | −1.00 | **−0.983 ± 0.071** |
| width only, H = 240 | 192×240, 384×240, 768×240 | −0.50 | **−0.556 ± 0.059** |
| height only, W = 384 | 384×120, 384×240, 384×480 | −0.50 | **−0.494 ± 0.059** |

Each family follows the law its own geometry implies, and each sits more than six standard
deviations from the other law. A fixed-height sweep really does give one over the square root of
width; the manuscript's inverse-width law holds for the geometry the manuscript ran, in which both
dimensions scale together.

### 5. The floor is set by map area alone

| W | H | area | aspect | floor | floor × √area |
|---|---|---|---|---|---|
| 192 | 120 | 23,040 | 1.60 | 0.10775 ± 0.00421 | 16.36 |
| 192 | 240 | 46,080 | 0.80 | 0.08159 ± 0.00315 | 17.51 |
| 384 | 120 | 46,080 | 3.20 | 0.07569 ± 0.00313 | 16.25 |
| 384 | 240 | 92,160 | 1.60 | 0.05685 ± 0.00312 | 17.26 |
| 384 | 480 | 184,320 | 0.80 | 0.03711 ± 0.00278 | 15.93 |
| 768 | 240 | 184,320 | 3.20 | 0.03716 ± 0.00302 | 15.95 |
| 768 | 480 | 368,640 | 1.60 | 0.02582 ± 0.00300 | 15.67 |

The floor falls as area to the power −0.509 ± 0.027 against a predicted −0.500, which is 0.3σ,
across a sixteenfold range in area and three aspect ratios. The two maps at equal area and very
different shape, 384×480 at aspect 0.80 and 768×240 at aspect 3.20, give floors agreeing to 0.1 per
cent. The floor does not depend on width, height or shape. It depends on area.

**Where the law comes from.** A map holds roughly `N_eff = area / λ_c²` independent coherence
cells. A variance estimated from `N_eff` independent samples carries a relative uncertainty of
order `√(2/N_eff)`, so its square root, which is what ΔT is, carries about `1/√(2 N_eff)`. Hence

```
floor = λ_c / √(2 · area)
```

with no dependence on shape. Scaling both dimensions makes area grow as L² and gives 1/L; scaling
one dimension makes area grow as L and gives 1/√L. Both measured slopes follow from the same
expression.

Calibrated against the measurement, `floor × √area = 16.39`, implying an effective coherence length
λ_c = 23.2 pixels. Wavelengths are drawn log-uniformly on 18 to 260 pixels with the amplitude
envelope peaked at 70, so a coherence length of 23 pixels sits inside the sampled range at the
short end. The exponent is predicted from first principles and confirmed at 0.3σ; the coefficient
is calibrated from the data, not predicted.

### 6. Which randomness drives the scatter

Each of the five random draws varied alone, with the other four held fixed.

| condition | CoV at N = 12 | CoV at N = 384 |
|---|---|---|
| all vary | 0.18152 | 0.07061 |
| wavelength only | 0.18142 | 0.05879 |
| amplitude factor only | 0.04544 | 0.02391 |
| direction only | 0.06411 | 0.05810 |
| spatial phase only | 0.04184 | 0.05900 |
| temporal phase only | 0.20389 | 0.06666 |

The 0.806 prediction accounts for wavelength, amplitude factor and temporal phase, and does not
account for direction or spatial phase. That is visible in the derivation itself:
`variance_campaign.py` computes `r` from `(env·u)² · cos²(th)`, where the envelope carries the
wavelength dependence, `u` is the amplitude factor and the cosine stands in for the temporal phase.
Direction and spatial phase enter only through the spatial basis.

### 7. Validation of the stochastic machinery

`validation_trials.py` compares `NU_PEAK_PER_K · t_avg` against `NU_PER_K_IND · t_avg`, the same
measured map mean on both sides, so it tests the constants and never touches the field. The tests
in `validate_stochastic.py` do test the field.

The five draws, pooled over 2,000 realisations at 96,000 samples each, Kolmogorov-Smirnov against
the intended distribution:

| draw | intended | D | p |
|---|---|---|---|
| log-wavelength | U(ln 18, ln 260) | 0.0018 | 0.92 |
| amplitude factor | U(0.6, 1.4) | 0.0033 | 0.26 |
| direction | U(0, 2π) | 0.0031 | 0.31 |
| spatial phase | U(0, 2π) | 0.0025 | 0.59 |
| temporal phase | U(0, 2π) | 0.0032 | 0.27 |

The normalisation reproduces the coefficient it should: q = E[a⁴]/E[a²]² = 2.4009 by 400-point
Gauss-Legendre quadrature and 2.4011 by four-million-sample Monte Carlo, giving A = q/8 = 0.3001
and a naive coefficient of 0.8064, both confirming the manuscript's 0.300 and 0.806 independently.

For a sum of randomly phased modes the ensemble two-point covariance is
`C(r) = norm² · Σⱼ (cⱼ²/2) · J₀(kⱼ r)`. Measured over 120 realisations on a 256×256 map, the
empirical radial covariance matches that expression to within 1.4 per cent of C(0) at separations
up to 10 pixels, degrading to 4.1 per cent at 40 pixels. That degradation grows with separation,
which is where a finite map runs out of independent separations, so the one place the analytic
prediction fails is the finite-map effect this paper is about.

### 8. Figure

`figure1_revised.png` and `.pdf`, authored at the width it is placed. Panel (a) is scatter against
wave count with bootstrap error bars; panel (b) is floor against area for all seven maps with the
equal-area pair labelled; panel (c) is the three geometries each against its own scaled dimension.

Panel (b) replaces the submitted Figure 1b, which plotted three points against width with a dashed
line normalised at the first point, so part of the agreement was built in. The new panel plots
seven points against area with a one-parameter law and no normalisation freedom.

`figure1_as_submitted.png` is the original figure, retained for comparison.

## License

MIT. See `LICENSE`.
