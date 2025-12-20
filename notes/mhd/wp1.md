# WP1 — Topology Switching Gate: Implementation and Verification Report

## WP1 Goal

The purpose of WP1 was to establish a **minimal, controllable “topology switch”** in the toy superfluid model: a mechanism that can be turned **ON** and **OFF** in time (and optionally localized in space) to transition the system between:

* an **ideal / frozen-in** regime where vortex topology is largely conserved, and
* a **non-ideal / topology-change-accessible** regime where reconnection-like events and structural changes can occur.

WP1 was not intended to prove “fusion relevance” or optimize energy yield. Its purpose was to deliver a **reproducible experimental control knob** plus **diagnostics** showing that toggling this knob produces a measurable, time-localized change in system behavior.

---

## Conceptual Model: What the “Switch” Represents

In ideal MHD, topology is conserved by the frozen-in condition ( \mathbf{E} + \mathbf{v}\times\mathbf{B} = 0 ). Topology change requires a non-ideal term (e.g., resistivity) so that field lines (or vortex lines in the analogy) can reconnect.

In the toy superfluid setting, the WP1 switch is implemented as a **pulsed diffusion / damping channel** (an “effective resistivity”) that is activated for a finite time interval and can be spatially localized. Operationally, this creates a controlled interval during which the system can violate strict frozen-in behavior and undergo topology-changing dynamics.

---

## Implementation Summary

### A) Gate Definition (Pulse + Optional Localization)

The gate is defined as a time-dependent amplitude multiplied by a spatial window:

* **Temporal gate:** a finite-duration pulse beginning at (t_0) and lasting (\tau).
* **Spatial gate:** a Gaussian window (W_x(\mathbf{x})) centered at a chosen location with width (\sigma).

The effective diffusivity takes the form:
[
\eta(\mathbf{x},t) ;=; \eta_{\text{peak}},W_t(t;t_0,\tau),W_x(\mathbf{x};\sigma)
]
In the implemented script, (\eta_{\text{peak}}) is set by `D_ETA` and (W_t) is implemented as a smooth trapezoid (ramped on/off) to reduce numerical artifacts.

### B) Dynamics: Unitary Evolution + Non-Ideal Slot

The core evolution is a split-step (spectral) real-time Gross–Pitaevskii-style update (unitary part), augmented by a **non-unitary diffusion step** during the gate window:

[
\psi \leftarrow \psi + \Delta t,\eta(\mathbf{x},t),\nabla^2\psi
]

This diffusion term is the explicit “switch slot” that can open access to topology change (analogous to resistive MHD terms enabling reconnection).

### C) Stability Fix: Automatic Diffusion Subcycling

A key practical WP1 requirement was numerical stability. An explicit diffusion term can blow up if (\Delta t) is too large relative to (\eta) and the grid’s maximum wavenumber. Early tests produced a catastrophic energy spike (orders of magnitude), indicating an explicit diffusion instability.

WP1 resolved this by implementing **automatic diffusion subcycling** during the gate:

* When the gate is ON, the diffusion step is split into `n_sub` smaller substeps.
* `n_sub` is chosen so that a diffusion CFL-like stability bound is satisfied:
  [
  \Delta t_{\text{sub}},\eta_{\text{peak}},k^2_{\max} ;\le; \text{DIFF_LIMIT}
  ]

This eliminated blow-ups and enabled stable parameter exploration.

---

## Diagnostics Added in WP1

WP1 required at least one **measurable indicator** that changes during/after the gate window in a reproducible way.

### 1) Acoustic / Compressible Energy Proxy (“Sound”)

A compressible component of the mass current is computed via a Helmholtz-like projection in Fourier space. This yields a scalar compressible energy proxy (E_c(t)) (“sound”) used as a candidate reconnection/activity indicator.

### 2) Core Volume Fraction (Topology Proxy)

We compute the fraction of grid cells satisfying:
[
\rho(\mathbf{x},t) = |\psi|^2 < \rho_{\text{th}}
]
This metric (“core fraction”) measures how much of the domain is in low-density core-like state.

### 3) Connected Component Count of the Core Mask (`n_big`)

To detect discrete topological/structural changes, WP1 added a connected component analysis of the low-density mask. Components are counted using 26-connectivity in 3D. Only components above a minimum voxel count are considered “macroscopic”:

* threshold: `RHO_TH`
* minimum size: `MIN_VOX`
* output: `n_big(t)` = number of large core components

This provides a coarse but highly practical topology-proxy event detector: when structures split/merge, `n_big(t)` can change.

### 4) Visualization

WP1 outputs:

* a 3D scatter of low-density voxels (for qualitative inspection),
* a time series plot of compressible energy with the **gate overlay**, and
* a time series plot of `core_frac(t)` and `n_big(t)` with the **gate overlay**.

---

## Canonical Test Setup (WP1 Baseline)

WP1 used a consistent “canonical initial condition” so runs are comparable:

* **Canonical IC:** orthogonal vortex tubes
  (two vortex tubes initialized orthogonally in the 3D domain)

Representative successful gate configuration:

* `D_ETA = 0.08`
* `PULSE_T0 = 8.0`
* `PULSE_TAU = 8.0`
* spatial localization via Gaussian window (parameter `PULSE_SIGMA` used in the script)

This run produced a dimensionless and operational summary:

* **diffusion length during gate:**
  [
  \ell_d=\sqrt{\eta_{\text{peak}}\tau}\approx \sqrt{0.08\times 8}\approx 0.80
  ]
* **χ parameter (with δ=1 convention):**
  [
  \chi=\ell_d/\delta\approx 0.80
  ]
* **local gate RMS speed estimate:**
  [
  U_0 \approx 2.36
  ]
* **local Rm estimate (δ=1 convention):**
  [
  Rm\sim U_0\delta/\eta_{\text{peak}}\approx 286
  ]
  These values indicate the system remains largely “ideal” globally (large (Rm)) but the gate introduces a controlled non-ideal interval with meaningful diffusion length compared to the assumed thickness scale.

---

## Observed WP1 Results

### A) Stable Gate Operation

After introducing diffusion subcycling, the gate no longer produced catastrophic numerical blow-ups. The evolution remained stable for the WP1 parameter choices used in the test runs.

### B) Gate-Synchronized Behavioral Change in Diagnostics

Across the tested runs, the diagnostics show time-localized changes correlated with the gate window:

1. **Compressible energy (“sound”) response:**
   The compressible energy time series remains oscillatory, but shows a gate-correlated change (including post-gate suppression/settling in some runs). This demonstrates that the gate modifies compressible activity, although sound alone is not a uniquely diagnostic reconnection indicator.

2. **Topology proxy response (`n_big`):**
   The connected-component proxy exhibits excursions during the gate window—i.e., `n_big(t)` shows discrete changes while the gate is active compared to baseline behavior. This is consistent with structural splitting/merging of low-density regions during the non-ideal interval.

3. **Core fraction response (`core_frac`):**
   `core_frac(t)` stays within a bounded range but responds dynamically during the run, providing context for interpreting `n_big(t)`.

Collectively, these confirm that turning the gate ON produces a measurable, time-localized change in system behavior detectable by at least one topology-proxy diagnostic.

---

## WP1 Completion Criteria and Status

WP1 is considered complete because it delivered:

1. **A defined switch:** pulsed, optionally localized non-ideal diffusion channel.
2. **A stable implementation:** diffusion subcycling prevents blow-ups.
3. **A measurement suite:** sound proxy + core fraction + component count, all plotted with gate overlay.
4. **Empirical evidence:** gate-window-correlated changes in the topology proxy (`n_big`) and in compressible behavior.

Therefore, WP1 successfully establishes a controllable mechanism to open/close topology-change accessibility and an instrumented diagnostic framework to detect its effects.

---

## WP1 Artifacts Produced

* Updated WP1 diagnostics runner (Python):

  * gate implementation (time + space)
  * stability subcycling
  * plots for sound and topology proxies
* A verified canonical IC test harness (orthogonal vortex tubes)
* A standard per-run “scalar summary” printout (e.g., (\ell_d), (\chi), (U_0), (Rm))

---

## Known Limitations (Accepted in WP1)

WP1 intentionally uses **coarse topology proxies**, not full vortex-line tracking. Connected components of a thresholded density field can be sensitive to threshold choice and to density waves. For WP1 this is acceptable because the goal is **switch verification**, not definitive reconnection counting.

Later work can strengthen topology measurement (e.g., filament extraction, phase-winding line tracking, reconnection event counting), but WP1 already provides an actionable experimental control knob and a working detection scaffold.

---

## Transition to WP2

With the switch and diagnostics in place, WP2 can now proceed to parameter sweeps and phase diagrams:

* sweep (\eta_{\text{peak}}), (\tau), (\sigma), and gate timing (t_0)
* define outcome regimes using:

  * changes in `n_big(t)` (event counts or excursions),
  * changes in compressible energy statistics (pre/gate/post),
  * and stability constraints (subcycling requirement)

WP1 thus provides the foundation for systematic exploration of “controlled single-event vs multi-event cascade” regimes in the toy model.

---

