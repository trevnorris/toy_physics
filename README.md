# Toy Superfluid Gravity

Lightweight research sandbox for a toy superfluid model that reproduces Newtonian gravity and leading 1PN effects discussed in the accompanying papers.

## Repository Layout
- `papers/` – Source and PDFs for:
  - *Paper 1: Newtonian and 1PN Orbital Dynamics from a Superfluid Defect Toy Model* (`1pn_orbital_dynamics.*`)
  - *Paper 2: 1PN Optics* (`1pn_optics.*`) – extends the model to light propagation, lensing, and Shapiro delay.
  - *Paper 3: Spin, Vorticity, and N-Body Dynamics* (`1pn_spin_and_nbody.*`) – spin precession, Lense–Thirring, and the EIH N-body Lagrangian.
  - *Paper 4: Electromagnetic Fields and Charged Defects* (`em_fields.*`) – emergent electromagnetism from the superfluid vacuum.
  - *Paper 5: Brane–Bulk Throat Ontology* (`brane_bulk_ontology.*`) – resolves the sphere–cylinder geometric tension via 4D bulk construction.
- `mathematica/` – Mathematica notebooks/scripts organized by paper (`1pn_orbital_dynamics/`, `1pn_optics/`, `1pn_spin_and_nbody/`, `em_fields/`, `brane_bulk_ontology/`, `2pn/`).
- `scripts/` – Small Python proofs-of-concept (Shapiro delay partitioning, ray-bending tests).
- `superfluid_lib/` – WIP physics engine (scalar Poisson + wave solvers, particle dynamics, 1PN orbit integrator); uses CuPy if available, otherwise NumPy.
- `experiments/` – Python experiments that exercise the engine (radial force law fit, Mercury perihelion calibration/precession).

## Papers Snapshot
**Paper 1 (Orbital Dynamics):** Constructs a superfluid defect toy model with two scalar potentials—an instantaneous Poisson sector and a finite-speed lag sector. The lag sector adds an attractive \(1/r^3\) correction that yields one half of the GR 1PN precession; introducing a mild position-dependent inertia with \(\beta = 3/2\) recovers the full GR 1PN perihelion shift. See `papers/1pn_orbital_dynamics.pdf`.

**Paper 2 (Optics):** Extends the model to light propagation, deriving the optical metric and recovering standard 1PN lensing and Shapiro delay results with PPN parameters \(\beta = \gamma = 1\). See `papers/1pn_optics.pdf`.

**Paper 3 (Spin & N-Body):** Promotes defects to composite "dyons" (flux-tube sink + vortex ring) whose far-field vorticity defines a gravitomagnetic vector potential. Reproduces Lense–Thirring frame dragging, spin precession, and the full Einstein–Infeld–Hoffmann N-body Lagrangian. Key finding: matching the EIH tensor requires an effective Lorentzian signature in the longitudinal sector, encoded by \(\alpha^2 = -2/5\). See `papers/1pn_spin_and_nbody.pdf`.

**Paper 4 (Electromagnetism):** Introduces a hydrodynamic dictionary mapping enthalpy and velocity fields to EM potentials—magnetic field becomes vorticity, electric field becomes minus the Euler acceleration. The breathing mode of a defect throat generates a Coulomb \(1/r\) potential, and the Lorentz force emerges from Magnus and pressure forces on vortices. Explains the EM/gravitational hierarchy as scaling with \(1/a^2\) (throat radius). See `papers/em_fields.pdf`.

**Paper 5 (Brane–Bulk Ontology):** Resolves the geometric tension between spherical (gravity) and cylindrical (EM) defect requirements by promoting defects to brane–bulk throats connecting the 3D brane to a 4D superfluid bulk. Dimensional reduction recovers the spherical monopolar far-field for 1PN gravity, while internal 4D acoustic modes retain cylindrical Bessel profiles for EM. Enthalpy minimization at fixed charge selects the preferred aspect ratio \(L/a \approx 1.85\). See `papers/brane_bulk_ontology.pdf`.

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
Files in `mathematica/` mirror the analytic steps in the papers, organized by subdirectory:
- `1pn_orbital_dynamics/` – Lagrangian/Hamiltonian derivations, beta consistency checks, added-mass pieces.
- `1pn_optics/` – Optical metric, polytropic index constraints, lensing/Shapiro derivations.
- `1pn_spin_and_nbody/` – Gravitomagnetic sector, EIH Lagrangian matching, alpha constraint derivation.
- `em_fields/` – Hydrodynamic-EM dictionary, breathing mode, Lorentz force from Magnus/pressure.
- `brane_bulk_ontology/` – Brane-bulk mode resonance, dimensional reduction, quadratic form diagonalization.
- `2pn/` – Work-in-progress 2PN extensions.

The `.wl` files are runnable scripts; each file's output is stored as a final comment block within the file itself, so you can inspect results without re-running.

## Misc Scripts (`scripts/`)
- `ray_tracing.py` – Hamiltonian ray-bending toy tests (refraction vs flow vs split model).
- `shapiro_delay.py` – Shapiro delay partitioning into refraction and flow contributions.

## Status / Next Steps
The paper series now covers 1PN orbital dynamics, optics, spin/N-body, electromagnetism, and brane–bulk ontology (Papers 1–5). Current work focuses on 2PN extensions (see `mathematica/2pn/`). The PDE/ODE engine in `superfluid_lib/` is still under active development (e.g., 4D grid support, tighter coupling between scalar sectors, better interpolation). Experiments currently focus on validation of the toy 1PN model.
