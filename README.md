# Toy Superfluid Gravity

Lightweight research sandbox for a toy superfluid model that reproduces Newtonian gravity and leading post-Newtonian effects.

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
- `4d_em_fields` — 4D electromagnetic sector
- `4d_plasma` — Plasma dynamics in the 4D model
- `atomic_p22_bridge` — Atomic P22 bridge construction

### Work in progress (not yet published)
Files for ongoing work that doesn't have a paper yet remain in the original directories:
- `mathematica/` — `moving_throat/`, `inner_throat/`
- `scripts/` — `moving_throat/`, `inner_throat/`, `4pn/`, and loose audit scripts for 2.5pn, 3pn, lepton, atom work
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
