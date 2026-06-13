# BUILD DIRECTIVE — Step 3b (CPU): preconditioner for the coupled-branch JFNK solve

**Plan ref:** `docs/branch_realization_execution_plan.md` §6 (stiff n=5 strategy — "geometric or
algebraic multigrid on the linearized operator … this is why the software stack choice matters") +
§7 step 3 ("preconditioning is the main project", per `decisions/01_stack_choice.md`). Pre-reg
`docs/stage1_preregistration.md` §J.2 (convergence). **Stack:** `decisions/01_stack_choice.md`
(torch, conservative FD/FV, float64, JFNK; **D2 production stack deferred with a PETSc-bias escape
hatch** — this step is the CPU experiment that informs the D2 trigger). **Builds on:** step 3a
(committed `dcd6a1b`) — the coupled stationary residual (`coupled_branch.py`), the Newton–Krylov
core (`newton.py`, GMRES-on-JVP + Armijo), and the resolution ladder + diagnostics.

**Motivating finding (step 3a).** Unpreconditioned JFNK does **not** scale: GMRES iteration counts
grow steeply under refinement (≈50→200 over 241→1401 DOF), so the laptop stalls at ~1400 DOF in the
time budget. This step adds a **preconditioner** to flatten the GMRES counts and unlock materially
higher resolution **on this laptop, before any GPU spend**.

**Contract.** Codex DESIGNS and WRITES the preconditioner (the interface, the concrete method, the
re-run). Claude REVIEWS code + output with clean agents + an independent re-run and iterates with
you until acceptance. **This directive states requirements and acceptance criteria only — it does
NOT prescribe the preconditioner method (geometric multigrid vs block-Jacobi/ILU vs operator-split
vs assembled-Jacobian algebraic) or whether to stay matrix-free vs assemble the sparse Jacobian;
those are your design call.** Run everything **locally on CPU**; no GPU, no runpods, no network.
Never wrap your own session in a shell `timeout`.

---

## 1. Objective

Make the **linear (Krylov) solve inside each Newton step** scale: introduce a **preconditioner** for
the coupled-branch linearized operator so that **GMRES iteration counts per Newton step become
roughly grid-independent (flatten) under refinement**, and re-run the resolution ladder to reach
**materially higher DOF on the laptop** with bounded GMRES counts. This is the plan §6 "main
project" and the experiment that informs the deferred **D2 decision trigger** — *"does the linearized
solve show scalable, target-blind Newton–Krylov behaviour under refinement?"* (`decisions/01_stack_choice.md`).

This is still the **engineering smoke** regime: natural O(1) **placeholder** parameters, **no**
physical packet, export guard untouched, **target-blind**. The preconditioner changes only *how fast
the linear systems are solved* — it must not change *what* is solved (the Newton fixed point).

---

## 2. Scope

### IN — mandatory

**A. A pluggable preconditioner interface** in the Newton–Krylov path (`newton.py`): right- or
left-preconditioned GMRES on the JVP, with the preconditioner as a swappable component (so future
preconditioners — and eventually the D2 PETSc escape hatch — drop in without touching the residual,
Jacobian, or extraction layers). Backend-neutral; explicit in config.

**B. At least one concrete preconditioner** for the coupled linearized operator. Plan §6's target is
geometric/algebraic multigrid; a **cheaper first cut is acceptable** (e.g. block-Jacobi / block-ILU
on the per-cell field block, a geometric two-/V-cycle with a simple smoother, or a sector/operator-
split preconditioner) **provided it demonstrably flattens the GMRES-vs-DOF curve**. You choose the
method and whether to stay matrix-free (multigrid smoothers / operator application) or assemble the
sparse Jacobian (for ILU/AMG). The acceptance test is the *scaling behaviour*, not a specific method.

**C. Correctness preservation (non-negotiable).** The preconditioner accelerates the Krylov solve
**only**; the residual definition, the Jacobian (JVP) definition, and the Newton fixed point are
unchanged. Demonstrate that the **converged solution at a fixed grid is identical (to solver
tolerance) with vs without the preconditioner**, the **coupled MMS still converges at order 2**, and
the **coupled JVP-vs-FD check still passes**.

**D. Re-run the resolution ladder** at the **same grids as step 3a** (241→1401 DOF) to show the
**before/after GMRES-count improvement directly**, then **extend the ladder as far as the laptop
allows** within a comparable per-level time budget. Report the **GMRES-count-vs-DOF curve** (the
flattening evidence) and the **new maximum DOF reached**.

### OUT — do not touch this step

The GPU / PETSc port (that is D2, later — this step only *informs* the trigger); the frozen physical
run, the `free_choice` value freeze, and the export-guard flip (all gated/deferred); the
field→coefficient extraction map; the WP3 tangent; any §H target or packet export; full PML/sponge.

---

## 3. Requirements

**R1 — Reuse, don't fork.** Build on `newton.py`, `coupled_branch.py`, `operators.py`, `grid.py`,
`manifest.py`. Add the preconditioner as a new abstraction wired into the existing JFNK loop; do not
duplicate the Newton/GMRES core or the coupled residual.

**R2 — Correctness-preserving (the #1 review axis).** The preconditioner must not alter the Newton
fixed point. Provide explicit evidence: the converged field at a fixed grid matches the
unpreconditioned step-3a solve to tolerance (same solution), coupled MMS stays order 2, coupled JVP
check passes. If preconditioning is found to change the converged solution, that is a bug to fix, not
a result.

**R3 — Effectiveness / scaling (the acceptance signal).** Report GMRES iteration counts per Newton
step vs DOF, **before vs after** preconditioning, on the shared 241→1401 ladder. The counts must
become **substantially flatter** (trending toward grid-independent) — this is the D2-trigger
evidence. Quantify the improvement (e.g. counts and wall-clock per DOF).

**R4 — Reach higher resolution.** With the preconditioner, push the ladder to **materially higher
DOF than the unpreconditioned ~1400** on this laptop (as far as it allows in a comparable time
budget). Report DOF / wall-clock / peak memory / Newton iters / final residual / GMRES counts per
level, and where the new laptop limit is.

**R5 — Determinism & manifest.** Deterministic; every run emits a manifest including the
preconditioner type + its parameters and all solver controls (explicit in config, no hidden
defaults).

**R6 — TARGET-BLINDNESS.** No §H target value anywhere; placeholder params only; the preconditioner
choice and its parameters are driven by conditioning/convergence, never by any target.

**R7 — NON-FROZEN smoke discipline / guard firewall intact.** No physical packet; do not write under
`research/pde_audit/simulation/output/`; do not import/modify/satisfy/trip the export guard; no
frozen-packet schema. Output labelled engineering smoke / placeholder parameters / target-blind.

**R8 — Environment & deps.** Prefer a **self-contained implementation** (torch / numpy / scipy.sparse
already present) so the torch stack keeps owning the operators. If you believe a **new dependency**
is genuinely the right call (e.g. `pyamg` for algebraic multigrid), **STOP and flag it** with the
reason and the tradeoff — do not silently install; adding a dependency is a methodology decision for
review, not a default.

**R9 — No fabricated success.** If a preconditioner does not flatten the GMRES counts, report it
honestly and diagnose why (spectrum, sector coupling, smoother choice); do not cherry-pick grids or
loosen tolerances to fake a flat curve. A *partial, understood* improvement is a legitimate result
and a basis for the next iteration; an unexplained or cherry-picked one is not.

---

## 4. Acceptance criteria

1. A **pluggable preconditioner interface** in the Newton–Krylov path, with **≥1 concrete
   preconditioner** implemented and selectable via config.
2. **Correctness preserved (R2):** converged solution identical (to tolerance) with vs without
   preconditioning at a fixed grid; coupled MMS still order 2; coupled JVP check passes.
3. **GMRES counts flatten (R3):** before/after GMRES-vs-DOF curve on the shared 241→1401 ladder
   shows substantially reduced growth (trending grid-independent); improvement quantified.
4. **Higher resolution (R4):** the ladder reaches materially higher DOF than ~1400 on the laptop,
   with bounded GMRES counts; per-level cost reported and the new limit stated.
5. Determinism + manifest (preconditioner type/params recorded); target-blind; guard firewall
   intact; non-frozen smoke discipline.
6. A concise **report** (markdown/YAML, KB-scale): the preconditioner method + why, the before/after
   GMRES-vs-DOF curve, the new max DOF, the correctness-preservation evidence, and an explicit read
   on the **D2 trigger** (does the preconditioned solve look scalable, or is a PETSc/multigrid port
   indicated?).
7. Target-blindness (R6) holds — verifiable by `grep`.

Claude reviews code + report (clean agent — correctness-preservation, target-blindness, guard
firewall) + an independent arbiter re-run reproducing the GMRES-vs-DOF curve; we iterate with you
until all hold.

---

## 5. Repo hygiene (firewall is live)

- **Code** → `software/stage1_solver/src/`, **tests** → `software/stage1_solver/tests/` (tracked).
- **All run output** → `software/stage1_solver/{runs,data,figures}/` — **gitignored**.
- **The one tracked artifact** → a small report under `software/stage1_solver/reports/`.
- Do **not** commit data; do **not** `git add` anything (orchestrator stages + commits after
  review); **never** `git add -A`; **never** write under `research/pde_audit/simulation/`.

---

## 6. Deliverable for review

Report back with: the preconditioner interface + the concrete method chosen and why (and matrix-free
vs assembled, with any new-dependency flag), the correctness-preservation evidence (same converged
solution; MMS order 2; JVP pass), the **before/after GMRES-vs-DOF curve**, the new maximum DOF
reached on the laptop, the report path, and your read on the **D2 trigger**. Then Claude reviews
before anything is committed.
