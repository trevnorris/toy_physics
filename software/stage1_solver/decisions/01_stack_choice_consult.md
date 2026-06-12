# Stage-1 solver — numerical-stack methodology consult (Claude → Codex)

**Type:** read-only methodology consult. Do NOT write or edit any code/files. Return a
recommendation with reasoning. This is a Claude+Codex math/methodology decision (the user
delegated the numerical-stack choice to us); escalate to the user only if something here
forces a change to what the model physically *claims* (it should not).

## What we are building

The **layer-3 Stage-1 branch-realization solver** for the moving-throat superfluid analog.
The pre-registration is **frozen** (`docs/stage1_preregistration.md`, target-blind). We must
now BUILD the solver that did not previously exist. Governing docs (already written, do not
re-derive): `docs/branch_realization_execution_plan.md` (engineering plan) and the frozen
pre-registration. Canonical equations: `notes/moving_throat_pde_program_compact.md`.

Physics shape (plan §2):
- **WP1** = stationary, isotropic (spherically symmetric) moving-throat branch. Stationary
  drops `t`; isotropy drops the angular `Omega`. Fields become functions of `(r, w)` only —
  a **2D mesh parameterizing 3D-with-axial-symmetry physics, full `r^2` measure factor**
  (NOT 2D physics). Fields: `psi_real, psi_imag` (Madelung ok), `A_0, A_r, A_w`, `R_0(w)`
  throat shape, one gauge multiplier → ~5–7 coupled fields/cell. Solved by **Newton**
  (no time-stepping). Gauge H=Z. Open finite throat: Dirichlet mouth, Robin/open-impedance
  exit, PML/sponge on the tangent.
- **WP3** = grouped real `P2` weak-axisymmetric tangent, LINEAR PDE linearized around the
  converged WP1 branch, same `(r, w)` mesh (2–5× cheaper than WP1).
- Natural units `a = c_s = c = G = 1`; **all targets are dimensionless O(1)–O(10^2) ratios**
  (so float64 dynamic range is NOT a problem here — unlike the prior 4d_1pn_sim attempt).
- Primary observable: `R_norm = 0` (externally GR-benchmarked). Extraction is **post-solve**.

## Build sequence (plan §7) — what each backend must support

1. **GPE soliton/vortex benchmark** on a 2D grid — KNOWN-ANALYTIC. Certifies discretization
   + boundary treatment. (standard cubic GPE, not the stiff EOS yet.) — *first deliverable.*
2. **Manufactured-solution per-operator tests** (matter / Maxwell / wall blocks) — discrete
   operator == continuum operator.
3. Stationary isotropic branch **Newton solve**, coarse grid — end-to-end smoke packet.
4. **Convergence study**, 3–5 grid levels — per-observable order of convergence.
5. **PML/sponge** characterization on the tangent — reflections below target signal.
6. **Conservation diagnostics** — mass/charge/energy drift.
7. **Error budget / noise floor** per observable.
8. **WP3 P2 tangent** on the converged branch.
9. Stage-1 verdict.

## The real failure mode (plan §6) — the deciding numerical question

**Stiff `P = K rho^5` polytrope.** Enthalpy `h(rho) = (5K/4) rho^4` becomes near-singular as
`rho -> 0` in the throat / near the wall → Newton Jacobian condition number ~1e6 (unfriendly,
not catastrophic). Mitigations on the table: geometric/algebraic **multigrid** on the
linearized operator; **change of variables** (`log rho` / regularized density); **continuation
in K** (mild EOS → stiff target). Plan §3.2 says the stack is chosen by *which path gives the
cleanest multigrid on the linearized operator*. **This stiffness only bites at step 3+, not
steps 1–2.**

## Hard environment facts (this machine, where steps 1–2 run for free)

- Python 3.10.12. **Installed:** numpy 1.26, scipy 1.15, **torch 2.11**, sympy 1.14, h5py,
  matplotlib. **NOT installed:** JAX, Dedalus, Firedrake/FEniCS/dolfinx, PETSc/petsc4py, mpi4py.
- **No local GPU** (nvidia driver absent). 16 CPU cores, 30 GB RAM.
- Production (steps 3+): runpods.io single A100/H100 (GPU spend is a USER gate; steps 1–2 are
  free local CPU).

## Two non-negotiable process constraints (from the freeze)

- **TARGET-BLINDNESS is the #1 audit axis.** The solver will be reviewed line-by-line for any
  target (the GR constant `54 G c_s^5/5 a^5 c^5`, `chi_Q=1`, the §H targets) leaking into
  initialization, branch selection, stopping criterion, or tuning. **Transparent custom code
  we can fully audit is worth more than a black-box framework here.**
- **VALIDATION-FIRST.** Known-answer certification (step 1–2) before ANY physical-branch
  extraction. The prior attempt failed by being noise-dominated/uncertifiable.

## The decision, split in two

**D1 (immediate — blocks step 1):** backend for steps 1–2 (local, free, known-analytic).
**D2 (deferred, informed by 1–2 conditioning data):** production GPU stack + the `n=5`
preconditioner strategy. Plan §3.2 explicitly says D2 is "not pre-committed."

Candidate stacks (plan §3.2): Dedalus (spectral; strong off-the-shelf nonlinear PDE + Newton +
continuation; install friction; spectral awkward for open-throat+PML?), Firedrake/FEniCS (FEM;
natural for r^2 measure + Robin BC + PML; heavy PETSc install, finicky), custom JAX (autodiff
JFNK, CPU-now/GPU-later, install needed, no off-the-shelf multigrid), custom torch (autodiff
JFNK, **already installed**, CPU-now/GPU-later same code, no off-the-shelf multigrid),
PETSc/Trilinos (serious multigrid, CPU-cluster).

## Claude's going-in recommendation (ratify or counter)

**D1: build steps 1–2 as a custom finite-difference GPE solver in torch** (with a thin
backend-abstraction so the array ops aren't torch-locked). Reasoning: (a) zero install — runs
here today on CPU; (b) **same code runs CPU-now and GPU-later** (device-agnostic tensors),
no rewrite; (c) autodiff gives the JFNK Jacobian for free; (d) fully transparent → satisfies
the target-blindness audit far better than a framework; (e) the manufactured-solution harness,
analytic references, observable-extraction, and convergence machinery are backend-agnostic and
carry forward regardless of the eventual D2 pick → bounded regret.

**D2: defer the production-stack lock-in until after step 2** (honors plan §3.2/§6). When the
step-1–2 conditioning data exists, decide between (i) staying in torch and adding a custom
geometric-multigrid V-cycle / continuation-in-K, or (ii) porting just the linearized operator
to PETSc / a spectral Dedalus solve if torch FD conditioning proves inadequate.

## Questions for Codex (please answer each)

1. **Rewrite-trap check:** does starting in custom torch FD create a throwaway trap, or does it
   genuinely carry forward to a viable production solve? If a trap, name the specific wall.
2. **Should we instead start directly in Dedalus or JAX** for steps 1–2, accepting install
   friction now to avoid a later port? Concrete pros/cons.
3. **float64 concern:** the stiff Jacobian wants float64; torch float64 on A100 GPUs is much
   slower than float32. Is that a strong enough argument to prefer a CPU-multigrid (PETSc) or
   spectral (Dedalus) production path — and does it change D1?
4. **Discretization family** for WP1's open-finite-throat `(r,w)` domain with full `r^2`
   measure, Robin exit impedance, and PML on the tangent: finite-difference vs finite-element
   vs spectral — which is least likely to fight the geometry/BCs?
5. **Preconditioner plan for `n=5`:** rank multigrid vs `log rho` change-of-variables vs
   continuation-in-K as the FIRST thing to try at step 3, given the chosen D1 backend.
6. Anything in the §7 sequence or the frozen pre-registration that the stack choice could
   accidentally compromise (esp. target-blindness)?

Return: a clear D1 recommendation (ratify/counter + why), a D2 direction (with the trigger that
would settle it), and concise answers to Q1–Q6. No code.
