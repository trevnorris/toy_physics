# Toy Superfluid Gravity

Lightweight research sandbox for a toy superfluid model that reproduces Newtonian gravity and leading post-Newtonian effects.

> **Where is the work right now?** See [`STATUS.md`](STATUS.md) — the single front door to current state, the next step, and a map of which doc holds what.

## Repository Layout

### `research/`
Each published paper lives in its own folder with subdirectories for the paper source, Mathematica derivations, Python scripts, and notes (where applicable):

```
research/<paper_name>/
  paper/        .tex and .pdf
  mathematica/  .wl derivation scripts
  scripts/      Python verification / simulation scripts
  notes/        working notes (when present)
```

Papers currently in `research/`:
- `1pn_orbital_dynamics` — Paper 1: Newtonian and 1PN orbital dynamics
- `1pn_optics` — Paper 2: 1PN light propagation, lensing, Shapiro delay
- `1pn_spin_and_nbody` — Paper 3: Spin precession, Lense-Thirring, EIH N-body
- `em_fields` — Paper 4: Emergent electromagnetism from superfluid vacuum
- `brane_bulk_ontology` — Paper 5: Brane-bulk throat ontology
- `1pn_hybrid` — Paper 6: Unified scalar/vector/optical 1PN model
- `4d` — 4D bulk superfluid construction
- `4d_1pn_bridge` — Bridge from 4D bulk to 1PN observables
- `4d_1pn_full` — Full 4D 1PN derivation
- `4d_2pn` — 2PN extensions
- `4d_2_5pn` — Conditional 2.5PN dissipative sector
- `4d_3pn` — Full conservative 3PN sector
- `4d_4pn` — Conditional full conservative 4PN sector
- `4d_em_fields` — 4D electromagnetic sector
- `4d_plasma` — Plasma dynamics in the 4D model
- `pde` — Moving-throat PDE paper: geometry lift, reduced wall-support-gauge dynamics, quadrupole response (draft)
- `pde_ledger` — Moving-throat PDE derivation companion archive ledger
- `research_proposal` — Far-field-first research proposal for the one-medium brane--bulk model (draft)

### Work in progress (not yet published)
Files for ongoing work that doesn't have a paper yet remain in the original directories:
- `mathematica/` — `moving_throat/`, `inner_throat/`
- `scripts/` — `moving_throat/`, `inner_throat/`, `4pn/`, and loose audit scripts for 3pn, lepton, atom work
- `notes/` — `moving_throat/`, `inner_throat/`, `4pn/`, and various top-level working notes
- `notes/summaries/` — Per-paper summaries (kept together for easy batch upload)

### `superfluid_lib/`
WIP physics engine (scalar Poisson + wave solvers, particle dynamics, 1PN orbit integrator); uses CuPy if available, otherwise NumPy.
- `core.py` — Unit scaling helper and grid classes (3D FFT-ready; 4D stub)
- `solvers.py` — Spectral Poisson solver and scalar wave equation solver with CFL guard
- `dynamics.py` — Particle ensemble, force accumulator, leapfrog stepper
- `pn_orbit.py` — Fast two-body 1PN toy integrator with scalar 1/r^3 correction

### `experiments/`
Python experiments that exercise the engine:
- `verify_radial_law.py` — Fits inverse-square slope from a Gaussian source
- `verify_mercury_calibrated.py` — Mercury perihelion: Newtonian vs scalar-only vs full 1PN

Both use CuPy on GPU if installed; otherwise CPU.

## Python Environment
- Core deps: `numpy`, `scipy`, `matplotlib`
- Optional GPU accel: `cupy` (falls back to NumPy when unavailable)
- Python 3.9+

## Mathematica Notes
The `.wl` files in each `research/<paper>/mathematica/` directory are runnable scripts. Each file's output is stored as a final comment block within the file itself, so you can inspect results without re-running.

## Papers List

All papers are archived on Zenodo under author Norris, T. (2026). Within each section, papers are listed in reading order (foundations first, extensions later).

### Derivation archive

- *Moving-Throat PDE Derivation Companion Archive Ledger.* [10.5281/zenodo.19699523](https://doi.org/10.5281/zenodo.19699523)

### Unified 4D Toy Model series

- *4D Toy Model — Action, Projections, and Controlled Brane Limits.* [10.5281/zenodo.19449589](https://doi.org/10.5281/zenodo.19449589)
- *Deriving Key Post-Newtonian Coefficients from the Unified 4D Toy Model.* [10.5281/zenodo.19449653](https://doi.org/10.5281/zenodo.19449653)
- *Maxwell from the Unified 4D Toy Model: Localized 5D Gauge Dynamics, Brane Reduction, and KK Corrections.* [10.5281/zenodo.19449834](https://doi.org/10.5281/zenodo.19449834)
- *Plasma Dynamics from the Unified 4D Toy Model: A 4+1D Alternative to MHD with a Controlled MHD Limit and 4D Interaction Corrections.* [10.5281/zenodo.19449935](https://doi.org/10.5281/zenodo.19449935)
- *A Controlled Derivation of the Full First Post-Newtonian Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19450102](https://doi.org/10.5281/zenodo.19450102)
- *A Controlled Derivation of the Full Conservative Second Post-Newtonian Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19450284](https://doi.org/10.5281/zenodo.19450284)
- *A Conditional Derivation of the Point-Particle 2.5PN Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19492270](https://doi.org/10.5281/zenodo.19492270)
- *A Full Conservative Derivation of the Two-Body 3PN Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19501724](https://doi.org/10.5281/zenodo.19501724)
- *A Conditional Full Conservative Derivation of the Two-Body 4PN Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19561056](https://doi.org/10.5281/zenodo.19561056)

### Superfluid Defect Toy Model series

- *Newtonian and 1PN Orbital Dynamics from a Superfluid Defect Toy Model.* [10.5281/zenodo.19449058](https://doi.org/10.5281/zenodo.19449058)
- *Gravitational Optics and Soliton Geodesics in a Superfluid Defect Toy Model.* [10.5281/zenodo.19449182](https://doi.org/10.5281/zenodo.19449182)
- *Spin, Vorticity, and N-Body Dynamics in a Superfluid Defect Toy Model.* [10.5281/zenodo.19449261](https://doi.org/10.5281/zenodo.19449261)
- *Electromagnetic Fields and Charged Defects in a Superfluid Defect Toy Model.* [10.5281/zenodo.19449355](https://doi.org/10.5281/zenodo.19449355)
- *Brane-Bulk Throat Ontology for a Superfluid Defect Toy Universe.* [10.5281/zenodo.19449422](https://doi.org/10.5281/zenodo.19449422)
- *Hybrid 1PN Dynamics and the Acoustic Horizon in a Superfluid Defect Toy Model.* [10.5281/zenodo.19449512](https://doi.org/10.5281/zenodo.19449512)
