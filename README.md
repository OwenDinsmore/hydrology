# Darcy flux mass balance

A single Jupyter notebook that visualizes a 2D Darcy flux field and computes the water balance of a
1 cm control volume three different ways.

**Rendered notebook:** https://owendinsmore.github.io/hydrology/

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

- `darcy_flux_mass_balance.ipynb`: the analysis (interactive 3D view of the control volume, flux matrices, flow plots,
  three methods, mass balance, convergence)
- `requirements.txt`: exact dependency pins, including [RiskPlot](https://pypi.org/project/riskplot/) for the heatmaps,
  field contours, interactive surface, and budget waterfall
- `.github/workflows/pages.yml`: executes the notebook with `nbconvert --execute` and publishes the HTML to GitHub Pages

## Run locally

Python 3.12:

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
jupyter nbconvert --to html --execute darcy_flux_mass_balance.ipynb --output-dir _site --output index.html
```
