# Hydrology

Four homework pieces, all published with GitHub Pages:

| | Page | Source |
|---|---|---|
| HW 1 · Darcy flux mass balance | https://owendinsmore.github.io/hydrology/ | `darcy_flux_mass_balance.ipynb` |
| HW 2 · Pore-ensemble SWRC simulator | https://owendinsmore.github.io/hydrology/swrc/ | `swrc/` |
| HW 3 · Brooks–Corey fit and Poiseuille conductivity | https://owendinsmore.github.io/hydrology/conductivity/ | `conductivity/` |
| HW 4 · Steady moisture profiles above a water table | https://owendinsmore.github.io/hydrology/profiles/ | `profiles/` |
| Research · Layered and branching pore networks | https://owendinsmore.github.io/hydrology/networks/ | `networks/` |

## Research · Layered and branching pore networks

A static web app (`networks/`, no build step) that builds five pore networks on the same soil and cylinder,
drains them, and solves for conductivity in two directions. Each network is a graph of bonds that hold the
water and set the capillary barriers. Nodes are zero-volume junctions.

| Model | Idea | Closest prior work |
|---|---|---|
| Straight bundle | HW 2 and 3 baseline | Purcell; K = ρg/8μ · Σr⁴/R² |
| 1 · Stacked plates | Plates one grain thick with random holes, joined only where holes overlap (Sampson orifice r³/3μ) | Childs and Collis-George (1950), Mualem (1976) |
| 2 · Gapped plates | Plates separated by a rough gap (b³/12μ, drains at 2σ/ρgb); optional side holes drilled through the glass | Guido's suggestion; b → 0 gives model 1 |
| 3 · Branching tree | Vertical trunks with side branches bifurcating by Murray's law; tip fusion closes loops | Bethe-lattice percolation, fractal tree networks |
| Pore-throat lattice | Cubic lattice thinned to coordination z | Fatt (1956) pore-network models |

- Pore radii come from d₅₀: r̄ = 0.10 d₅₀ and layer thickness = d₅₀ by default, with r̄ and the gap aperture
  scaled by the void ratio φ/(1−φ). φ is the total porosity (holes, gaps, side holes, branches); the sample is
  resized until the total matches.
- Drainage is quasi-static invasion from the top face. With "air must reach a pore" on, a bond drains only when
  its entry radius exceeds r* = 2σ/ρg|h| and air can reach it through drained bonds.
- A drained bond keeps 0.395 (r*/r)² of its water in its corners and κ (r*/r)⁴ of its conductance (after
  Tuller and Or 2001). κ = 1.5 × 10⁻⁴ is fitted to measured sands, see below. κ = 0 gives round tubes, which
  lose all conductivity while S is still 0.3 to 0.5.
- K_v comes from Kirchhoff's law with the top face at unit head and the bottom at zero. K_h uses two 80° rim
  electrodes, scaled by the disk's own shape factor from a finite-difference solve, so an isotropic lattice gives
  K_h/K_v ≈ 1.
- Each run reports K_sat against Kozeny–Carman, (ρg/μ) d₅₀²/180 · φ³/(1−φ)², and Mualem's ℓ fitted to the
  network, ln K_r − 2 ln M(S) = ℓ ln S, where M(S) comes from the network's own θ(h).
- The solver is conjugate gradients with an IC(0) preconditioner, because conductances span about ten decades.
- `node --test tests/networks.test.mjs` checks the bundle against the HW 3 formula, flow conservation,
  isotropy of a uniform lattice, monotone drainage, total porosity, corner flow, and the calibration below.

The gapped-plate defaults (r̄ = 0.10 d₅₀, b = 0.26 d₅₀) were calibrated at φ = 0.35 to two targets:
K_h/K_v ≈ 1 and K_sat near Kozeny–Carman. With void-ratio scaling they then hold from φ = 0.25 to 0.45
(K/KC 1.45 to 1.56, K_h/K_v 0.91 to 1.05) and at any d₅₀.

Default results (fine sand, mean ± s.d. over five seeds):

| Model | K_sat / Kozeny–Carman | K_h/K_v |
|---|---:|---:|
| Straight bundle | 2.9 ± 0.4 | 0 |
| 1 · Stacked plates | 0.24 ± 0.03 | 0.014 ± 0.005 |
| 2 · Gapped plates | 1.49 ± 0.02 | 0.98 ± 0.03 |
| 3 · Branching tree | 4.4 ± 0.9 | 0.011 ± 0.004 |
| Pore-throat lattice | 0.039 ± 0.002 | 1.20 ± 0.08 |

The fitted Mualem ℓ is not a stable property of any model here. For gapped plates it moves from −0.77 at
κ = 0.05 to 1.42 at κ = 1.5 × 10⁻⁴, so it is reported on the page but should not be quoted as a result.

### Test against measured data

`node tools/validate.mjs <DataFiles_v2> [cache.json] [κ]` tests the gapped plates against Hohenbrink et al.
(2023), doi:10.5880/fidgeo.2023.012 (CC BY 4.0, not stored here): 69 sand and loamy-sand samples with
evaporation-method K(h), 38 of them with a measured K_sat. The network sees only each sample's retention
curve. Its pore-size spread and one length scale are fitted to θ(h), then K(h) is predicted and compared
with the measured points. Errors are medians, in decades of log₁₀ K.

| Prediction | K data used | Absolute error | Shape error (best shift removed) |
|---|---|---:|---:|
| Gapped plates, κ = 0.05 (first guess) | none | 2.43 | 0.14 |
| Gapped plates, κ = 1.5 × 10⁻⁴ | one global κ | 0.49 | 0.23 |
| Mualem–van Genuchten, ℓ = 0.5 | measured K_sat per sample | 0.86 | 0.38 |
| van Genuchten–Mualem, K_s and ℓ fitted | all K points per sample | 0.16 | – |

- The measured points sit at a median saturation of 0.12, so this tests the dry end, where corner flow sets
  the level. The bias moves one decade per decade of κ, and κ was chosen to remove it. It is fitted on the
  same samples, so 0.49 is not an out-of-sample error.
- By lab at κ = 1.5 × 10⁻⁴: TUBS (49 samples) 0.36, KIT (6) 0.44, UFZ (13) 1.14.
- Network K_sat is a median 0.39 decades (2.5×) above measured for the 38 samples that have one.
- Not tested on sands: saturations between about 0.3 and 1, and horizontal conductivity.

**Held-out check on finer textures.** Passing a texture list as the fifth argument runs other classes with κ
frozen. On all 559 usable samples (σ grid widened to 1.4), the model fails outside sands:

| Texture | Samples | Retention fit RMSE (S) | Absolute error | Shape error | Mualem + measured K_sat |
|---|---:|---:|---:|---:|---:|
| Sand | 46 | 0.06 | 0.61 | 0.24 | 1.29 |
| Loamy sand | 24 | 0.11 | 0.85 | 0.29 | 0.58 |
| Sandy loam | 131 | 0.17 | 1.69 | 0.68 | 0.86 |
| Silt loam | 141 | 0.13 | 2.48 | 0.61 | 1.76 |
| Loam | 44 | 0.16 | 2.36 | 0.64 | 1.69 |
| Clay | 26 | 0.09 | 3.52 | 0.61 | 1.70 |

- The failure starts in the retention fit: a single unimodal set of capillary holes cannot reproduce the
  retention curves of loams and clays. The K error follows from that, so these numbers do not test κ itself.
- The network's K_sat is a median 0.79 decades below measured across 402 samples, consistent with structural
  macropores the model does not have.
- Adding σ = 1.4 to the library moved the sand results (Sa 0.61, LoSa 0.85, against 0.49 for both together
  before), so the sand error is sensitive to the fitting grid at about ±0.1 to 0.3 decades.

**Film water.** Clay and silt surfaces hold water as adsorbed films. The model now adds a non-capillary store
θ_nc beside the network (`film` parameter). Its retention follows the PDI model (Peters 2013; Iden and Durner
2014): full up to the network's air entry, then log-linear to zero at pF 6.8. Its flow follows Tuller and Or
(2001): van der Waals film thickness δ = (A/6πρgh)^⅓ with A = 6 × 10⁻²⁰ J, giving
K_film = ρg/9μ · θ_nc(h) · δ². θ_nc is the only new parameter and is fitted to retention, never to K.
The capillary porosity is then θ_s − θ_nc. Same 11 classes, κ frozen, θ_nc fitted per sample
(524 samples; 35 dropped because their capillary porosity fell below the library's 0.1):

| Texture | Median θ_nc | θ fit, no film → film | Absolute error, no film → film | Shape error | Mualem + measured K_sat |
|---|---:|---:|---:|---:|---:|
| Sand | 0.09 | 0.019 → 0.010 | 0.62 → 0.60 | 0.25 | 1.59 |
| Loamy sand | 0.15 | 0.038 → 0.010 | 0.85 → 0.72 | 0.33 | 0.76 |
| Sandy loam | 0.26 | 0.062 → 0.011 | 1.69 → 0.88 | 0.46 | 0.78 |
| Silt loam | 0.35 | 0.057 → 0.012 | 2.50 → 1.62 | 0.74 | 1.80 |
| Loam | 0.35 | 0.063 → 0.011 | 2.36 → 1.74 | 0.81 | 2.11 |
| Clay | 0.43 | 0.050 → 0.010 | 3.60 → 2.75 | 0.69 | 1.74 |
| All | | 0.060 → 0.012 | 2.33 → 1.35 | 0.57 | 1.34 |

- Films fix retention everywhere (median error 0.06 → 0.012 in θ) and cut the conductivity error by about
  one decade in fine soils.
- Overall the network now ties Mualem scaled to measured K_sat (1.35 against 1.34) while using no K data.
  It beats it in loam and silt loam and loses in clay and sandy loam.
- The shape is still poor outside sands: 0.57 decades against Mualem's 0.15. Measured points in fine soils sit
  at a median saturation of 0.67, the wet range, where structure and macropores set the slope. Film flow is
  negligible there (about 10⁻¹¹ m/s), so it changes the level through the retention fit, not by carrying water.
- The κ that would remove the bias still spans 2 × 10⁻⁷ (clay) to 9 × 10⁻⁴ (silty clay), so the corner
  constant does not transfer between textures.

**Macropores (tested, not adopted).** A structural pore set was added in parallel: straight, randomly oriented
tubes with a lognormal radius by volume (σ = 0.5), so K_ma = ρg/24μ · w · E_v[r²; r < r*]
(`macroSaturation`, `macroK`; `node tools/validate.mjs … macro`). Volume w and median radius r_m are fitted to
retention with θ_nc. It made no difference outside sands: w came out at 0.003–0.015, absolute error 1.37 and
shape error 0.57, both unchanged. The reason is the sign of the error. Against log h across pF 2–3, the shape
residual has a median slope of −1 to −2.6 decades per decade in every texture except clay, so the network's K
falls *too fast* with drying. Macropores add K only at the wet end, which steepens the curve further.

**Corner exponent.** The steep part is the corner conductance (r*/r)⁴ ∝ h⁻⁴ in drained pores
(`cornerExp`, default 4; `CORNER_EXP=n node tools/validate.mjs …`). The library was rebuilt at other exponents,
with κ held at 1.5 × 10⁻⁴, and run on all textures:

| Corner exponent | 4 (triangle) | 3 | 2 |
|---|---:|---:|---:|
| Shape error, all | 0.57 | 0.50 | 0.44 |
| Shape error, sandy loam | 0.46 | 0.34 | 0.20 |
| Shape error, loamy sand | 0.33 | 0.21 | 0.15 |
| Shape error, silt loam | 0.74 | 0.66 | 0.58 |
| Absolute error, all (κ not re-tuned) | 1.35 | 1.29 | 1.61 |

Re-tuning κ for exponent 2 on the sands alone gives κ = 2 × 10⁻⁶. On the sands that gives median bias 0.06,
absolute error 0.60 and shape error 0.14, against Mualem's 0.86. On the held-out textures it is worse than
exponent 4: absolute error 1.69 against 1.35 overall. Sandy loam keeps a good shape (0.22) but sits
0.8 decades low. Loam and silt loam shape errors double to about 1.0. With a smaller κ, the curve drops
sharply where the capillary network disconnects. The κ that would remove the bias spans 3 × 10⁻⁹ (clay) to
7 × 10⁻⁵ (silty clay), wider than at exponent 4.

A shallower exponent fixes the shape in coarse-to-medium soils. Fine soils stay near 0.5–0.6 even once their
slope residual is small, so what remains there is curvature, not slope. The default stays at 4 because a
lower exponent has no geometric basis yet. One candidate is rough-walled or partly filled corners
(Tuller & Or 2001).

## HW 3 · Conductivity

A one-screen static web app (`conductivity/`, no build step, built on `swrc/styles.css`) that extends HW 2 to K(Θ_E).

- Draws 10⁶ pores from an exponential with r̄ = 0.01 mm and removes everything above 3r̄, in a sample of radius
  2 cm. That leaves about 950,000 pores and φ ≈ 0.288.
- h_A = 2σ/(ρ g r_max) ≈ 490 mm. λ and c are set by sliders, with a button that snaps to the least-squares fit.
- K = (ρg/8μ) Σ r⁴ / R² over the pores still full at each h, which gives K_sat ≈ 1.35 × 10⁻⁴ m/s.
- ρg sits in the numerator of Q, per the corrected eq. 2 (ΔP + ρgΔz = ρgΔH). It sets K_sat but not K/K_sat or c.

| Fit range | λ | c | 1 + 2/λ |
|---|---:|---:|---:|
| Θ_E ≥ 0.1 | 1.71 | 2.02 | 2.17 |
| Θ_E ≥ 0.01 | 1.99 | 1.91 | 2.00 |
| Θ_E ≥ 0.001 | 2.17 | 1.85 | 1.92 |
| Θ_E ≥ 0.0001 | 2.30 | 1.82 | 1.87 |

## HW 2 · Pore-ensemble SWRC

A static web app (`swrc/`, no build step) that builds a soil water retention curve from a bundle of
straight capillaries drilled through a cylinder.

- Sample N pore radii from an exponential, gamma, or log-normal distribution with mean r̄, optionally
  truncated so no pore is larger than r_max.
- Capillary rise sets the head at which each pore drains: h = −2σ cos γ / (ρ_w g r). With the default
  fluid constants that is h ≈ −14.7 / r (h in mm, r in mm).
- Porosity is the summed pore volume over the cylinder volume. Length cancels, so φ = Σ r_i² / R².
- At a given head every pore with r > r* = −2σ cos γ / (ρ_w g h) has drained; θ(h) is what remains.
- The page plots θ against |h| on a log axis, shows porosity, θ, saturation, r*, and drained-pore
  counts for the current head, and renders a rotatable 3D cylinder in which a sample of pores is
  coloured by whether it still holds water. Buttons re-run the ensemble, draw a new seed, animate a
  drainage sweep, reset the defaults, and export the curve as CSV.

With the sketch's defaults (10⁶ pores, r̄ = 0.01 mm, D = 20 mm) the summed pore area is about twice the
cylinder cross-section, so the reported porosity is about 2. The page flags this rather than hiding it.

## HW 1 · Darcy flux mass balance

A single Jupyter notebook that visualizes a 2D Darcy flux field and computes the water balance of a
1 cm control volume three different ways.

## Problem

Flux field (cm/hr), with no z-component and no z-dependence:

```
q_x(x, y) = 3x^2 + 5xy + 7y^2
q_y(x, y) = -2x^3 + 4x^2 y^3 - 2y
```

Continuity with no sources or sinks:

```
dθ/dt = -div(q) = -(dq_x/dx + dq_y/dy)
dq_x/dx = 6x + 5y
dq_y/dy = 12x^2 y^2 - 2
```

Control volume: x ∈ [2, 3], y ∈ [3, 4], dz = 1 cm, so V = 1 cm³ and every face area is 1 cm².
Flow is +x and +y everywhere on the box. The left and bottom faces are inflow; the right and top faces are outflow.

## Answers

| Method | dθ/dt (1/hr) | Q_in (cm³/hr) | Q_out (cm³/hr) |
|---|---:|---:|---:|
| 1. Point divergence at (2.5, 3.5) | −949.25 | — | — |
| 2. Face-center flux (midpoint rule) | −955.50 | 770.50 | 1726.00 |
| 3. Exact face integrals | −967.8333 | 778.8333 | 1746.6667 |

Exact face flows (cm³/hr): left 133.3333, bottom 645.5, right 165.8333, top 1580.8333.
The volume integral of div(q) is 967.8333, matching the net outflow as the divergence theorem requires.

All three methods share the same x-contribution (−32.50 1/hr) because dq_x/dx is linear. All of the disagreement
comes from the y-direction, where the 4x²y³ term makes q_y grow about 2.4× from the bottom face to the top.
As the box shrinks, the face-center estimate converges to the point divergence as −949.25 − 25h², where h is the half-width.

The notebook asserts every number above, so a successful run is also a verification.

## Contents

- `darcy_flux_mass_balance.ipynb`: the analysis, with three figures (the control volume in xyz space, the q_y surface,
  and streamlines over log10|q|). The two 3D views are drawn as SVG, so they need no WebGL.
- `requirements.txt`: exact dependency pins, including [RiskPlot](https://pypi.org/project/riskplot/) for the
  flux-magnitude contour behind the streamlines
- `.github/workflows/pages.yml`: executes the notebook with `nbconvert --execute` and publishes the HTML to GitHub Pages

## Run locally

Python 3.12:

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
jupyter nbconvert --to html --execute darcy_flux_mass_balance.ipynb --output-dir _site --output index.html
```
