# BUILD DIRECTIVE — Step 4 (CPU): grid-convergence study of the coupled-branch solve

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 4 ("Convergence study at 3–5 grid
levels … establishes noise floor and convergence order for each observable", brief §5.2). Pre-reg
`docs/stage1_preregistration.md` **§J.2 (non-negotiable):** "refine grid (≥3 levels); report
observed order of convergence per observable; an observable that drifts under refinement is not a
measurement" + §J.5 (state the numerical noise floor explicitly). **Stack:**
`decisions/01_stack_choice.md`. **Builds on:** step 3a (`dcd6a1b`) coupled stationary residual +
step 3b (`affc745`) preconditioner — REUSE both; the preconditioner is what makes the finer levels
of this ladder feasible on the laptop.

**Purpose (why now, on CPU).** The user chose to run this convergence study on CPU *before* any GPU
pivot, to MEASURE the discretization-error-vs-resolution behaviour of the coupled solve and thereby
size the resolution the physics will need (and hence the GPU requirement) from data rather than a
guess. This is the engineering-smoke, placeholder-parameter, **target-blind** precursor to the frozen
§J.2 study.

**Contract.** Codex DESIGNS and WRITES the study (refinement ladder, the self-convergence /
Richardson machinery, the diagnostics, the report). Claude REVIEWS code + output with clean agents +
an independent arbiter re-run and iterates with you until acceptance. This directive states
requirements + acceptance criteria only; the discretization/refinement/extrapolation design is your
call. Run everything LOCALLY ON CPU; no GPU, no network. Never wrap your own session in a shell
`timeout`.

---

## 1. Objective

On the coupled stationary isotropic branch (matter quintic gauged-GNLS ⊗ localized-Maxwell H=Z, the
step-3a residual) with the step-3b preconditioner, run a systematic grid-refinement ladder and
establish, for a set of **target-blind surrogate observables**:

1. the **observed order of convergence** per observable (Richardson p from level triples),
2. the **discretization error vs resolution** (per-level error estimate + a Richardson-extrapolated
   continuum estimate with an error bar),
3. a **preliminary numerical noise/error floor** (where successive-level differences stop shrinking —
   set by Newton tolerance, GMRES tolerance, round-off, extrapolation uncertainty), and
4. an **extrapolated resolution read**: from the error-vs-DOF curve, the DOF needed to reach a couple
   of stated *illustrative* error levels (e.g. 1e-3, 1e-4) — i.e. the resolution-sizing number that
   informs the later GPU decision.

This certifies that the SOLVE (operator + boundary treatment + Newton fixed point), not just the
isolated operators (step 2 MMS), converges under refinement — and quantifies how fast.

---

## 2. Scope

### IN — mandatory

**A. A systematic refinement ladder.** A fixed refinement ratio (e.g. 2× per dimension per level),
**≥3 levels (target 4–5)**, run as far as the laptop allows with the preconditioner. Successive-level
differences must be well-defined: either nest the grids (coarse nodes ⊂ fine nodes) or interpolate
each solution to a common fine reference grid with an interpolation order ≥ the scheme order so the
interpolation does not pollute the measured order. State which, and why it is consistent.

**B. Self-convergence / Richardson per observable.** For each surrogate observable, report: the value
at each level, the successive-level change, the **observed order p** (from level triples at the fixed
ratio), and a **Richardson-extrapolated continuum value + per-level error estimate**. Where a fine
reference solve is the practical truth proxy, report error against it AND the ratio-based order (so
the two agree or the discrepancy is diagnosed).

**C. Target-blind surrogate observables ONLY.** Choose a defensible set of functionals of the RAW
solved fields, each independent of the extraction map and of every §H target — e.g.: a solution-field
norm (L2/Linf self-difference between levels), total mass `∫ρ dV`, peak density, the cross-sector
gauge L2 and current L2 already reported in step 3, a field-energy-like integral, and the converged
Newton residual floor. **NONE** may be a §H observable (`chi_Q`, `R_norm`/`R_pole`, `P_2`/`P_4`),
an extraction coefficient (`D0/D2/D4`, `P0/P2/P4`, `N0/...`), or anything compared to a target value.

**D. Honesty gate (§J.2).** "An observable that drifts under refinement is not a measurement." For
EACH surrogate, state plainly whether it converges at the expected order, at a reduced order, or not
at all — and DIAGNOSE any shortfall (boundary-induced order reduction at the open exit / Robin
mouth, the r=0 axis, near-vacuum density positivity, throat-geometry under-resolution, continuation
stiffness). A reduced or non-converging observable, honestly diagnosed, is a legitimate and useful
result; a cherry-picked order-2 table is not.

**E. Resolution-sizing read.** From the error-vs-DOF curve, extrapolate the DOF/grid needed to reach
a couple of illustrative target error levels, and state where the laptop limit was hit (DOF /
wall-clock / peak memory per level). This is the data that will inform the later CPU-vs-GPU /
scalable-preconditioner decision — report it as an engineering read, not a target comparison.

### OUT — do not touch this step

The field→coefficient extraction map; ANY §H observable or target (firewalled); the frozen physical
values (gated); the frozen §J.2 per-§H-observable study (deferred to the gated physical run — it
needs the extraction map + frozen values); PML/sponge (step 5); conservation budget (step 6 — though
mass drift may appear as a surrogate here); the WP3 tangent (step 8); GPU / PETSc port; any packet
export.

---

## 3. Requirements

**R1 — Reuse, don't fork.** Build on `coupled_branch.py`, `newton.py` + the step-3b preconditioner,
`operators.py`, `grid.py`, `manifest.py`, and the existing resolution-ladder/report scaffolding. Add
the convergence-study machinery as new functions; do not duplicate the residual, the Newton core, or
the preconditioner.

**R2 — Same problem at every level.** Every level solves the SAME coupled problem (same placeholder
parameters, same BC class, same continuation schedule, same gauge) — ONLY the grid changes. Any
per-level change other than resolution invalidates the study; if a coarse level needs a different
continuation path to converge, report it explicitly as a caveat rather than silently differing.

**R3 — Defensible difference norm.** The norm used for self-convergence must be measure-consistent
(use the conservative cell volumes / r² measure, not raw point counts) and identical across levels.
State it. If interpolating to a common grid, the interpolation order must not cap the measured order.

**R4 — No fabricated convergence (the acceptance integrity axis).** Do NOT pick the subset of levels,
the norm, or the observable that happens to show order 2 and hide the rest. Report every chosen
surrogate's honest order. Do not loosen Newton/GMRES tolerances to flatter the curve; if the solver
floor limits the finest useful level, say so and treat it as the noise floor (§J.5), not as
convergence.

**R5 — Determinism & manifest.** Deterministic (fixed seed, single-thread, float64). Every solve
emits a manifest (grid, solver controls incl. preconditioner, surrogate values). The study emits a
machine-readable convergence table.

**R6 — TARGET-BLINDNESS (the #1 review axis).** No §H target value anywhere — not in code, config,
params, observables, norms, tolerances, the extrapolation, or the resolution-sizing read. The
illustrative error levels (1e-3, 1e-4) are generic numerical thresholds, not targets. Surrogate
observables are extraction-free and target-free. Verifiable by grep.

**R7 — NON-FROZEN smoke discipline / guard firewall intact.** Placeholder parameters only, labelled
engineering-smoke / target-blind. No physical packet; do NOT write under
`research/pde_audit/simulation/output/`; do NOT import/modify/satisfy/trip the export guard scripts;
no frozen-packet schema; no extraction.

**R8 — Environment & deps.** Self-contained on the existing stack (torch/numpy/scipy). If you believe
a new dependency is genuinely needed, STOP and FLAG it with the reason + tradeoff — do not silently
install.

---

## 4. Acceptance criteria

1. A refinement ladder with **≥3 levels (target 4–5)** at a fixed ratio, run as far as the laptop
   allows with the preconditioner; per-level DOF / wall-clock / peak memory / Newton iters / final
   residual / GMRES counts reported, and the laptop limit stated.
2. Per surrogate observable: value per level, successive change, **observed order p**, and a
   **Richardson-extrapolated continuum estimate + per-level error estimate** (R2/R3 honoured).
3. **Honest §J.2 verdict** per observable (converges at order / reduced order / drifts), with a
   diagnosis for any shortfall (R4).
4. A stated **numerical noise/error floor** (§J.5) and an **extrapolated resolution-sizing read**
   (DOF for illustrative target errors; laptop limit).
5. Determinism + manifest + machine-readable convergence table (R5).
6. Target-blind (R6, grep-verifiable); guard firewall intact; non-frozen smoke discipline (R7); no
   new dependency (or STOP-flagged).
7. A concise **report** (markdown/YAML, KB-scale) under `reports/` capturing all of the above.

Claude reviews code + report (clean agent — convergence integrity / no fabricated order /
target-blindness / firewall) + an independent arbiter re-run reproducing the convergence table; we
iterate with you until all hold.

---

## 5. Repo hygiene (firewall is live)

- **Code** → `software/stage1_solver/src/`, **tests** → `software/stage1_solver/tests/` (tracked).
- **All run output** → `software/stage1_solver/{runs,data,figures}/` — **gitignored**.
- **The one tracked artifact** → a small report under `software/stage1_solver/reports/`.
- Do **not** commit data; do **not** `git add` anything (orchestrator stages + commits after
  review); **never** `git add -A`; **never** write under `research/pde_audit/simulation/`.

---

## 6. Deliverable for review

Report back with: the refinement-ladder design (ratio, levels, nested-vs-interpolated + why) and the
difference norm; the per-observable convergence table (value / order / Richardson estimate / error)
and the honest §J.2 verdict per observable with diagnoses; the noise-floor read; the
extrapolated resolution-sizing read + the laptop limit; the report path; and any STOP-and-flag
items. Then Claude reviews before anything is committed.
