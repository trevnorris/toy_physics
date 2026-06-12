# Stage-1 Branch-Realization Solver

A small, self-contained **numerical PDE solver and validation harness** for the Stage-1
*branch-realization* test of a toy superfluid-analog model. It is written in
[PyTorch](https://pytorch.org/) (float64, CPU-now / GPU-later), with a clean separation between
the physics-residual layer and the array backend.

> **Scope / claims.** This is a *toy analog* program: the goal is an *operationally demonstrable*
> mathematical bridge between EM-like and gravity-like derivations, **not** a claim that the model
> describes physical reality. Nothing in this directory computes or asserts a physical result yet —
> the code here only **validates the numerical machinery** on problems whose answers are already
> known. See `../../docs/branch_realization_execution_plan.md` for the full program.

---

## Why this exists

The model is being subjected to a three-layer falsification program:

1. **Layer 1** — line-by-line red-team of the symbolic PDE ledger (complete).
2. **Layer 2** — adversarial "fit-vs-derive" audit of every constant (complete).
3. **Layer 3** — *this* solver: actually realize the pre-registered branch numerically and compare
   the extracted observables to their targets, **within validated numerical error**.

A **pre-registration** (`../../docs/stage1_preregistration.md`) freezes — *before any solve* —
which branch is realized, the gauge/boundary conventions, the extraction formulas, and which
observable is primary. The freeze makes Stage 1 a **target-blind, falsifiable test**: a miss on the
pre-registered branch falsifies *that branch*; it may not be retuned to manufacture a pass.

Because the targets are frozen and the solve must be blind to them, the **first job is trust**: we
certify the solver on known answers before we point it at anything unknown. That is what the code in
this directory currently does.

---

## Build status

The build follows the sequence in `../../docs/branch_realization_execution_plan.md` §7:

| Step | What it certifies | Status |
|------|-------------------|--------|
| 1 | Solver skeleton + known-analytic GPE benchmark (linear HO eigenproblem; cubic GPE vs an independent SciPy BVP) | ✅ done |
| 2 | **Manufactured-solution (MMS) tests per operator block** — each discrete operator converges to its continuum operator at formal order | ✅ done |
| 3 | Coarse stationary isotropic branch Newton solve | ⏳ not started |
| 4 | 3–5 level convergence study | ⏳ |
| 5 | PML / sponge characterization | ⏳ |
| 6 | Conservation diagnostics | ⏳ |
| 7 | Error-budget statement | ⏳ |
| 8 | WP3 weak-axisymmetric P2 tangent | ⏳ |
| 9 | Stage-1 verdict | ⏳ |

No physical branch is solved yet, and **no physical observable has been extracted.**

---

## How to run

Requires Python 3 with `torch`, `numpy`, `scipy`, `sympy` (all already in the project env). All runs
are CPU-only, deterministic (fixed seed, single-threaded, deterministic kernels), and finish in
seconds.

From the repository root (`/var/projects/toy_physics`):

```bash
# Run the full validation gate (step-1 benchmarks + step-2 MMS suite).
# Writes a markdown report and exits 0 on pass, 1 on fail.
PYTHONPATH=software/stage1_solver/src python3 -m stage1_solver.harness

# Run the unit tests.
PYTHONPATH=software/stage1_solver/src python3 -m pytest -q software/stage1_solver/tests
```

The harness prints `validation gate passed` and writes
`reports/step2_manufactured_solutions_validation.md` (the step-1 report is
`reports/step1_gpe_benchmark_validation.md`). Per-run manifests and any generated arrays go under
`runs/` (git-ignored — see *Repo hygiene* below).

---

## What to expect

Every operator block is certified by the **method of manufactured solutions** (MMS): pick a smooth
analytic field, apply the *continuum* operator to it symbolically (SymPy) to get an exact forcing,
then confirm the *discrete* operator reproduces that forcing and that the error falls at the
scheme's formal order under grid refinement.

The current suite reports clean **2nd-order** convergence for all five blocks, e.g.:

| Operator block | Observed order | Finest-grid error |
|----------------|---------------:|------------------:|
| Matter — quintic gauged-GNLS (`P = K ρ⁵`, `h = 5K/4 ρ⁴`) | → 2.00 | ~3e-5 |
| 2D `(r,w)` tensor Laplacian | → 2.00 | ~2e-4 |
| Complex current `j = (ℏ/m) Im(ψ̄ ∂ψ)` (nonzero analytic current) | → 2.00 | ~1e-5 |
| Localized Maxwell (gauge `H=Z`) | → 2.00 | ~7e-4 |
| Wall effective-closure `S_η⁽²⁾` (ℓ=2 lane) | → 2.00 | ~7e-5 |

plus a Newton/Jacobian self-consistency check (matrix-free JVP vs a finite-difference probe agreeing
to ~1e-11) and machine-precision mass conservation on the GPE benchmarks.

### What the numbers mean — and what they do **not**

- A clean 2nd-order MMS result certifies that the **discretization is faithful to the continuum
  operator it encodes** — the building blocks (measures, boundary operators, Newton machinery) are
  correct to the scheme's formal order.
- It does **not** by itself assert any physics: MMS proves *discrete ≈ continuum*, so it is only as
  meaningful as the continuum operator being faithfully transliterated from the frozen source. Each
  operator here is checked against its source equation independently (compact-program §2.4–§2.5,
  stage-001 wall action, pre-registration §D gauge), and the MMS forcing is derived from the
  continuum operator — never from the discrete code under test — so the certification is not
  circular.
- The **wall** constitutive coefficients used in step 2 are explicitly-labelled **MMS-only
  placeholders**. The physical wall packet is `free_choice` (pre-registration §E) and gets frozen
  per-run at solve time; step 2 only certifies the wall operator's *discretization*, which is
  independent of the coefficient values.
- **Target-blindness is deliberate.** No physical target value appears anywhere in the code,
  config, manufactured fields, or tolerances. This is enforced and reviewed; it is the firewall the
  pre-registration freeze protects.

---

## Layout

```
software/stage1_solver/
├── src/stage1_solver/      # the solver + validation harness (tracked)
│   ├── backend.py          # float64 + determinism configuration
│   ├── config.py           # explicit, hashable run configuration
│   ├── grid.py             # radial / (r,w)-tensor / wall grids; conservative-FV measure & face areas
│   ├── boundaries.py       # Dirichlet / Neumann / Robin boundary operators (ghost values)
│   ├── operators.py        # conservative-FV operators: radial & tensor Laplacian, current,
│   │                       #   localized Maxwell, wall; integration & current diagnostics
│   ├── physics.py          # residuals: linear HO, cubic GPE (step-1 stand-in), quintic gauged-GNLS
│   ├── references.py       # INDEPENDENT references (closed-form linear; SciPy solve_bvp cubic)
│   ├── newton.py           # Newton–Krylov (GMRES on JVP) + Armijo + finite-difference JVP check
│   ├── benchmarks.py       # step-1 known-analytic benchmark suite
│   ├── mms.py              # MMS convergence harness (weighted-norm error → observed order)
│   ├── mms_benchmarks.py   # manufactured fields + SymPy continuum forcing + per-sector MMS runs
│   ├── report.py           # validation-report generator
│   ├── manifest.py         # per-run reproducibility manifest
│   └── harness.py          # CLI entrypoint (runs step-1 + MMS, writes report, sets the gate)
├── tests/                  # unit tests (tracked)
├── reports/                # small markdown validation reports (tracked)
├── directives/             # per-step build directives — requirements & acceptance criteria (tracked)
├── decisions/              # numerical-stack and methodology decisions (tracked)
├── runs/ data/ figures/    # generated output — GIT-IGNORED
└── _scratch/               # scratch / logs — GIT-IGNORED
```

---

## Numerical approach

A **custom conservative finite-difference / finite-volume** scheme on a logically-rectangular
`(r,w)` grid, in **float64**, solved by a **matrix-free Newton–Krylov** method (GMRES on
Jacobian-vector products from autodiff, with an Armijo line search). The conservative form carries
the full `4πr²` radial measure through face/cell weights, so `r=0` regularity is automatic (the
inner face has zero area). Rationale and the (deferred) production-stack escape hatch are in
`decisions/01_stack_choice.md`.

---

## Repo hygiene (firewall)

The solver can generate large mesh/field/checkpoint data. To keep the repository small, **all run
output is git-ignored**; only source, tests, the small markdown reports, directives, decisions, and
this README are tracked.

- **Tracked:** `src/`, `tests/`, `reports/*.md`, `directives/`, `decisions/`, `README.md`.
- **Ignored:** `runs/`, `data/`, `figures/`, `_scratch/`, and large data formats
  (`*.npz/*.h5/*.pt/*.vtk/...`) anywhere.
- Do **not** `git add` generated data, and never `git add -A`. Stage explicit paths.

---

## Governing documents

- `../../docs/branch_realization_execution_plan.md` — the Stage-1 build plan (§7 step sequence).
- `../../docs/stage1_preregistration.md` — the frozen, target-blind test definition.
- `decisions/01_stack_choice.md` — the numerical-stack decision.
- `directives/step1_gpe_benchmark.md`, `directives/step2_manufactured_solutions.md` — per-step build
  directives (requirements + acceptance criteria).
