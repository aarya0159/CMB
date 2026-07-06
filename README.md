# CMB Blackbody Simulator

An interactive, single-file simulator of the cosmic microwave background (CMB): the blackbody glow left over from the Big Bang, together with the tiny temperature ripples (anisotropies) that mark the seeds of galaxies and dark matter structure.

Built with plain **JavaScript and HTML5 Canvas**. No libraries, no build step, no installation: the entire simulator is one `index.html` file that runs in any modern browser.

This simulator accompanies our research paper, *"The Universe as a Blackbody: Exploring the Cosmic Microwave Background and Its Tiny Temperature Ripples with an Interactive Simulator"* (Geet et al., 2026). Every experiment in Section 4 of the paper can be reproduced here.

## Run it

**Option 1 (instant):** download `index.html` and double-click it. It opens in your browser.

**Option 2 (live on the web):** this repository can be served with GitHub Pages. In the repo settings, go to *Pages*, choose *Deploy from a branch*, select `main` and the root folder, and save. The simulator will be live at `https://<username>.github.io/cmb-simulator/`.

## What you are looking at

The animated map is a smooth random temperature field: red blobs are directions slightly hotter than average, blue blobs slightly colder, the same colour convention as real CMB maps. The panel below recalculates the blackbody physics live for whatever universe you dial in.

## Controls

| Control | Keyboard | What it does |
|---|---|---|
| Universe Age / Cooling slider | LEFT / RIGHT | Sets the average temperature, from 3,000 K (the moment the CMB was released) down to 2.72548 K (today). Cooling the universe is the same as moving forward in time. |
| Density Perturbation Amplitude slider | UP / DOWN | Sets the ripple strength ΔT/T, from 10⁻⁶ to 10⁻¹. The real universe's value is about 10⁻⁵, which is invisibly small; try it. |
| Pause | SPACE | Freezes the field for inspection. |
| New Seed | R | Draws a fresh random universe with the same physics (cosmic variance). |

## Telemetry

| Readout | Formula behind it |
|---|---|
| Average Temperature | Mean of the rendered temperature map |
| Peak Frequency (Planck) | ν_peak = 58.7893 GHz/K × T (frequency form of Wien's law, x* = 2.8214) |
| Peak Wavelength (Wien) | λ_peak = 2.897772 mm·K / T (wavelength form; note it points at a different part of the spectrum than the frequency peak, which is not a bug) |
| ΔT (rms anisotropy) | Standard deviation of the rendered map |
| Redshift z / Age | z = T/2.72548 − 1, and age ≈ 13.8 Gyr / (1+z)^{3/2} |

## Suggested experiments (from the paper)

1. **Wien linearity:** step the temperature down and record (T, peak frequency) pairs; the ratio is always 58.79 GHz per K.
2. **One part in 100,000:** lower the amplitude toward 10⁻⁵ and watch the map turn perfectly uniform while the ΔT readout stays nonzero.
3. **Ride the redshift:** cool from 3,000 K to today and watch a near-infrared glow (966 nm) stretch into millimetre microwaves.
4. **Peak shift:** multiply the ΔT readout by 58.79 GHz/K to get the peak-frequency spread between hot and cold blobs; scaled to the real universe's 10 µK this gives 588 kHz.
5. **Cosmic variance:** press R repeatedly and watch the measured ΔT wobble between random universes.

## Physics constants used

T₀ = 2.72548 K (Fixsen 2009) · Wien frequency constant 58.7893 GHz/K · Wien wavelength constant 2.897772 mm·K · recombination at ≈ 3,000 K, z ≈ 1,100, age ≈ 378,000 years.

## License

MIT. See `LICENSE`.
