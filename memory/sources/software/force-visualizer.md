---
title: Analog four-force visualizer
type: source
status: current
sources:
  - software/force_visualizer/README.md
  - software/force_visualizer/notes/build_spec.md
  - software/force_visualizer/app.py
  - software/force_visualizer/render_all.py
  - software/force_visualizer/report.py
  - software/force_visualizer/output/verification_report.txt
last_updated: 2026-09-03
---

# Analog four-force visualizer

## Purpose

This package presents a calibrated 2D visualization of the toy model's
gravity, light, electric, and magnetic effective laws. It has an interactive
Matplotlib application, a headless four-GIF renderer, and a deterministic
numeric report. It does not solve the full nonlinear PDE or demonstrate that
the four forces emerge from one microscopic model.

Sources:

- `software/force_visualizer/README.md` — opening scope statement
- `software/force_visualizer/notes/build_spec.md` — headings `What this is (and is NOT)` and `Out of scope / do-not`

## What it shows

- Gravity: a Newtonian/1PN orbit and inward medium-flow tracers, with an
  optional calibrated 2.5PN radiation-reaction benchmark.
- Light: two transverse shear waves, lensing, and a deliberately exposed stray
  longitudinal component.
- Charge: sign-dependent Coulomb motion and field vectors, with an optional
  Yukawa partner correction.
- Magnetism: moving end-on current sources, circulating fields, and attraction
  or repulsion for parallel or antiparallel currents, including an exposed
  scalar admixture.

The live app supplies sliders and a departures checkbox for each sector.
`render_all.py` uses the noninteractive Agg backend and writes `gravity.gif`,
`light.gif`, `charge.gif`, and `magnetism.gif`; hiding departures is explicitly
counterfactual.

Sources:

- `software/force_visualizer/README.md` — interactive and field-visualization paragraphs
- `software/force_visualizer/app.py` — classes `_GravityView`, `_LightView`, `_ChargeView`, `_MagnetismView`, and `ForceVisualizerApp`
- `software/force_visualizer/render_all.py` — function `render_all`

## Numeric verification

`report.py` samples the rendering-independent physics core. The committed
report records 18/18 checks passing across 0PN energy conservation, 1PN
precession, gravity-field direction, light dispersion and polarization,
lensing, Coulomb falloff and signs, and magnetic falloff, signs, scalar ratio,
and circulation.

Those checks show that the implementation follows its chosen effective laws.
They do not turn calibrated magnitudes, cone choices, or source laws into
predictions of the parent toy model. In particular, the native 2.5PN
normalization remains blocked and its optional visualization is a standard
calibrated benchmark.

Sources:

- `software/force_visualizer/report.py` — functions `_add_gravity`, `_add_light`, `_add_charge`, `_add_magnetism`, and `build_report`
- `software/force_visualizer/output/verification_report.txt` — final summary and section `CHARACTERIZED DEPARTURES`

## Characterized departures

The default display keeps four limitations visible: gravity return residuals,
the light sector's longitudinal contaminant, a charge-sector Yukawa partner,
and magnetic scalar admixture/preferred-frame sensitivity. The build
specification is design and acceptance context; its `DRAFT` label should not
be mistaken for the status of the implemented application and stored report.

## Related memory

- `memory/scripts/force-visualizer.md` catalogs the three entry points.
