# Toy Superfluid Gravity

Lightweight research sandbox for a toy superfluid model that reproduces Newtonian gravity and leading 1PN effects discussed in the accompanying papers.

## Repository Layout
- `papers/` – Source and PDFs for:
  - *Newtonian and 1PN Orbital Dynamics from a Superfluid Defect Toy Model* (`1pn_orbital_dynamics.*`)
  - *1PN Optics* (`1pn_optics.*`) – extends the model to light propagation, lensing, and Shapiro delay.
- `mathematica/` – Mathematica notebooks/scripts that carry out the derivations in the papers.
- `scripts/` – Small Python proofs-of-concept (Shapiro delay partitioning, ray-bending tests).
- `superfluid_lib/` – WIP physics engine (scalar Poisson + wave solvers, particle dynamics, 1PN orbit integrator); uses CuPy if available, otherwise NumPy.
- `experiments/` – Python experiments that exercise the engine (radial force law fit, Mercury perihelion calibration/precession).

## Papers Snapshot
**Paper 1 (Orbital Dynamics):** Constructs a superfluid defect toy model with two scalar potentials—an instantaneous Poisson sector and a finite-speed lag sector. The lag sector adds an attractive \(1/r^3\) correction that yields one half of the GR 1PN precession; introducing a mild position-dependent inertia with \(\beta = 3/2\) recovers the full GR 1PN perihelion shift. See `papers/1pn_orbital_dynamics.pdf`.

**Paper 2 (Optics):** Extends the model to light propagation, deriving the optical metric and recovering standard 1PN lensing and Shapiro delay results. See `papers/1pn_optics.pdf`.

## Python Environment
- Core deps: `numpy`, `scipy`, `matplotlib`.
- Optional GPU accel: `cupy` (falls back to NumPy when unavailable).
- Experiments assume Python 3.9+; install deps via `pip install -r requirements.txt` if you add one, or install the packages above manually.

## Running the Experiments
From the repo root:
- Radial law check (fits the inverse-square slope from a Gaussian source):  
  `python experiments/verify_radial_law.py`  
  Saves `experiments/radial_law_nearfield.png`.
- Mercury perihelion calibration (Newtonian vs scalar-only vs full 1PN toy):  
  `python experiments/verify_mercury_calibrated.py`  
  Prints orbit period/precession stats and saves mode-specific orbit plots.

Both experiments will use CuPy on GPU if installed; otherwise they run on CPU.

## Library Notes (`superfluid_lib/`)
- `core.py` – Unit scaling helper and grid classes (3D FFT-ready; 4D stub).
- `solvers.py` – Spectral Poisson solver and scalar wave equation solver with CFL guard.
- `dynamics.py` – Particle ensemble container plus force accumulator (Poisson + optional wave potential, β-based inertia correction, Magnus/halo forces) and a simple leapfrog stepper.
- `pn_orbit.py` – Fast two-body 1PN toy integrator with the scalar \(1/r^3\) correction.

## Mathematica Derivations
Files in `mathematica/` mirror the analytic steps in the paper (Lagrangian/Hamiltonian derivations, beta consistency checks, added-mass pieces). The `.wl` files are runnable notebooks; `.out` files store cached outputs so you can inspect results without re-running heavy steps.

## Misc Scripts (`scripts/`)
- `ray_tracing.py` – Hamiltonian ray-bending toy tests (refraction vs flow vs split model).
- `shapiro_delay.py` – Shapiro delay partitioning into refraction and flow contributions.

## Status / Next Steps
The PDE/ODE engine is still under active development (e.g., 4D grid support, tighter coupling between scalar sectors, better interpolation). Experiments currently focus on validation of the toy 1PN model; feel free to add new scenarios or extend `superfluid_lib` as the model evolves.
