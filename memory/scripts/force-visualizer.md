---
title: Force visualizer script catalog
type: script_catalog
status: current
sources:
  - software/force_visualizer/app.py
  - software/force_visualizer/render_all.py
  - software/force_visualizer/report.py
  - software/force_visualizer/output/verification_report.txt
last_updated: 2026-09-03
---

# Force visualizer scripts

## `software/force_visualizer/app.py`

- **Role:** interactive Matplotlib front end for gravity, light, charge, and
  magnetism.
- **Inputs:** the rendering-independent `force_visualizer.physics` core and
  default visualization parameters.
- **Behavior:** constructs one view per sector, supplies live controls, updates
  field/trajectory artists, and exposes each sector's characterized departure
  by default.
- **Output:** an interactive desktop window; in a headless environment it
  exits with guidance to use `render_all`.
- **Run:** from `software/`, `python -m force_visualizer.app`.

Source: `software/force_visualizer/app.py` — classes `_GravityView`,
`_LightView`, `_ChargeView`, `_MagnetismView`, and `ForceVisualizerApp`.

## `software/force_visualizer/render_all.py`

- **Role:** headless renderer using Matplotlib's Agg backend.
- **Inputs:** output directory, animation parameters, optional
  `--include-25pn-benchmark`, and optional counterfactual
  `--hide-departures`.
- **Outputs:** `gravity.gif`, `light.gif`, `charge.gif`, and `magnetism.gif`.
- **Run:** `python -m software.force_visualizer.render_all --output-dir DIR`.

Source: `software/force_visualizer/render_all.py` — functions `render_all` and
`main`.

## `software/force_visualizer/report.py`

- **Role:** deterministic headless numeric verifier; it does not invoke the
  rendering layer.
- **Checks:** orbit energy and 1PN precession, gravity field direction, light
  dispersion/polarizations/lensing, electric falloff and signs, and magnetic
  falloff/signs/scalar ratio/circulation.
- **Output:**
  `software/force_visualizer/output/verification_report.txt` and the same
  transcript on stdout.
- **Run:** `python -m software.force_visualizer.report`.
- **Current stored result:** 18/18 implementation checks pass. This verifies
  the configured effective laws, not their derivation from the parent model.

Source: `software/force_visualizer/report.py` — functions `build_report`,
`write_report`, and `main`; stored report final summary.

Related source: `memory/sources/software/force-visualizer.md`.
