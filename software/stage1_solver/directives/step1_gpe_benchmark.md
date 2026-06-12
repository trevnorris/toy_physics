# BUILD DIRECTIVE — Step 1: GPE known-analytic benchmark + solver skeleton

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 1; pre-registration
`docs/stage1_preregistration.md` §J.1 ("validate before trusting"). **Stack:** decision
`software/stage1_solver/decisions/01_stack_choice.md` (torch, conservative FD/FV, float64).

**Contract.** Codex DESIGNS and WRITES all code (discretization scheme, benchmark profiles,
Newton/linear solver, harness). Claude REVIEWS code + output and iterates with you until the
acceptance criteria are met. This directive states *requirements and acceptance criteria only*
— it does NOT prescribe the discretization route; that is your design call. Run everything
**locally on CPU** (free); no GPU, no runpods, no network. Iterate until the harness exits 0
and the criteria below are met.

---

## 1. Objective

Stand up the **reusable solver skeleton** and certify its **discretization + boundary +
Newton machinery on known-analytic GPE problems**, BEFORE any moving-throat / Maxwell / wall /
target work. This is the validation-first gate: nothing physical is extracted until the machine
is proven on answers we already know.

The machinery you build here is **production machinery, not a throwaway toy** (stack-decision
audit risk #2). The discrete differential operators, the `r^2`/radial measure, the `r=0`
regularity treatment, the boundary-operator framework, the Newton/Jacobian harness, the
manifest/determinism discipline, and the convergence/conservation tooling MUST all be the same
ones WP1 will use later. Only the *nonlinearity* differs: step 1 uses the **standard cubic GPE**
(`g|psi|^2 psi`) because it has known-analytic references; the production **quintic**
`P=K rho^5` (compact §2.4.1, `h=5K/4 rho^4`, `c_s^2=5K/m rho^4`) and the Maxwell + wall sectors
arrive in steps 2–3.

Canonical equation source for the operators that must transfer:
`notes/moving_throat_pde_program_compact.md` §2.4.1 (gauged GNLS), §2.5.4 (Madelung/Euler form).
Read these so the discrete Laplacian / measure / current operators are continuum-faithful.

---

## 2. Scope

**IN:** the solver skeleton (config, grid, conservative FD/FV operators on a `(r,w)`-style
logically-rectangular grid with full radial measure and `r=0` regularity, boundary-operator
framework incl. Dirichlet + Robin, a Newton–Krylov core with autodiff/JVP Jacobian); a
**known-analytic GPE benchmark suite**; the **validation harness** (convergence, conservation,
Jacobian cross-check); the run-manifest/determinism discipline.

**OUT (do not touch this step):** the moving-throat physical branch (WP1), the quintic EOS,
the Maxwell/mixed sector, the wall/BdG sector, PML on the tangent (step 5), any §H target,
GPU/runpods, any observable from §G of the pre-registration.

---

## 3. Requirements

**R1 — Backend & structure.** torch; **float64 everywhere**; device-agnostic tensors (CPU now,
GPU later, same code). Keep the **physics-residual layer cleanly separated from the backend
array layer**. The validation harness must be backend-neutral. Conservative FD/FV form for any
weighted-divergence/measure term.

**R2 — Known-analytic benchmark suite.** Provide at least these two, both exercising the
production-transferable operators:
  - **(a) Linear/exactly-solvable limit** that pins the discrete operator + measure + `r=0`
    regularity + boundary treatment to a *closed-form* answer (e.g. a linear-Schrödinger /
    `g→0` eigenproblem with exact eigenvalues, or another exact solution of your choice).
    Purpose: certify the linear machinery to near-machine precision with a measured convergence
    order.
  - **(b) Nonlinear stationary cubic-GPE known answer** solved by the **same Newton machinery
    WP1 will use** (e.g. a trapped ground state compared in the Thomas–Fermi + healing-layer
    sense, or a quantized vortex profile `psi=f(r)e^{i l theta}`), validated against an
    **independent high-accuracy reference you compute separately** (e.g. a 1-D radial BVP /
    shooting solve at much higher accuracy) — NOT against itself.
  You may add more; these two are the floor.

**R3 — Newton / Jacobian integrity.** The autodiff (or matrix-free JVP) Jacobian MUST be
cross-checked against an independent finite-difference / directional-derivative probe and the
agreement reported (audit risk #6). All solver controls — tolerances, line-search, max iters,
stopping criterion, dtype — are **explicit in config and logged**, never left to hidden backend
defaults (audit risk #3).

**R4 — Convergence study.** ≥3 grid levels per benchmark; report the **observed order of
convergence per benchmark observable** and show the error driving down toward a noise floor.
The observed order must be consistent with the scheme's formal order; if it is not, that is a
finding to resolve, not to paper over.

**R5 — Conservation / symmetry.** Report norm/mass (and any other relevant invariant or exact
symmetry) drift for each run; demonstrate it sits at the noise floor.

**R6 — Determinism & manifest.** Runs are deterministic (fixed seeds / deterministic kernels).
Every run emits a **manifest** recording: dtype, all solver tolerances, mesh spec, config hash,
git revision, and library versions — the frozen-manifest discipline of
`research/pde_audit/simulation/NONLINEAR_PROTOCOL_V2.md` (Freeze Boundary), **minus any target
fields** (there are none at this step).

**R7 — TARGET-BLINDNESS (the #1 review axis).** No §H target constant — `54 G c_s^5/5 a^5 c^5`,
`chi_Q=1`, `R_pole/R_norm/P_2/P_4`, the GR quadrupole, any benchmark target value — may appear
anywhere in code, config, initialization, stopping criteria, or tuning. This step is a STANDARD
GPE benchmark with no connection to the physical targets; the discipline starts here.

**R8 — Environment & deps.** Use the installed env (numpy/scipy/torch/sympy/h5py/matplotlib).
If you believe a new dependency is genuinely needed, **STOP and flag it in your response** with
the reason — do not silently install.

---

## 4. Acceptance criteria (the gate to step 2)

1. Benchmark (a) reaches its closed-form answer at the scheme's formal convergence order, error
   → floor under refinement.
2. Benchmark (b) Newton-converges to the independent high-accuracy reference within a stated
   tolerance, at the expected order.
3. Jacobian cross-check (R3) passes within a stated tolerance.
4. Conservation/symmetry drift (R5) at the noise floor.
5. Runs deterministic; manifest (R6) emitted for every run.
6. A concise **validation report** is written (markdown or YAML) stating, per benchmark: grid
   levels, observed order, final error vs reference, Jacobian-check residual, conservation
   drift, noise floor, and the explicit solver config used.
7. Target-blindness (R7) holds — verifiable by `grep`.

Claude reviews code + report; we iterate until all hold.

---

## 5. Repo hygiene (the firewall is live — do not blow up the repo)

- **Code** → `software/stage1_solver/src/` and tests → `software/stage1_solver/tests/`
  (both tracked). Choose a clean module layout.
- **All run output** (fields, checkpoints, large arrays, plots) → `software/stage1_solver/runs/`
  / `figures/` / `data/` — these are **gitignored**; nothing there gets committed.
- **The one tracked result artifact** is the small validation report → put it under
  `software/stage1_solver/reports/` (tracked; keep it KB-scale, markdown/YAML).
- **Do NOT commit** `*.npz/*.h5/*.pt/*.vtk/*.pkl/...` or any large data — the `.gitignore`
  firewall blocks these; do not `git add -f` around it. Do not run `git add -A`.
- You apply + run + iterate; the orchestrator (Claude) handles staging/commits with explicit
  paths after review.

---

## 6. Deliverable for review

When the acceptance criteria are met, report back with: the module layout, the discretization
scheme you chose and why, the two+ benchmarks and their references, the convergence table, the
Jacobian-check result, the conservation numbers, and the path to the validation report. Then
Claude reviews before anything is committed.
