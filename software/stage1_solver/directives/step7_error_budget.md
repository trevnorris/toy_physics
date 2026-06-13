# BUILD DIRECTIVE — Step 7 (CPU): error-budget statement for the coupled branch

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 7 ("Error budget statement. Quantitative
uncertainty per extracted observable. *Brief §5.5.*"); brief `docs/branch_realization_brief.md` §5.5
("**Error budget and noise floor.** State the numerical noise floor explicitly. Every extracted observable
must carry a quantitative uncertainty. A match or mismatch with a target is only meaningful relative to
that floor."). Pre-reg `docs/stage1_preregistration.md` **§J.5 (non-negotiable):** "Error budget & noise
floor — state the numerical noise floor explicitly; every extracted observable carries a quantitative
uncertainty; a match/mismatch is meaningful only relative to that floor." **Stack:**
`decisions/01_stack_choice.md`. **Builds on / CONSUMES:** step 4 (`c9b8b2c`) convergence study —
per-observable observed order + Richardson/finest error estimate + the preliminary numerical floor + the
null/floor vocabulary; step 5 (`63cd885`) boundary characterization — the interior boundary-error floor;
step 6 (`4a03797`) conservation diagnostics — the Gauss-closure / conservation floor. This step is a
**composition + reporting** step over those three; it introduces no new heavy physics.

---

## Methodology scope (read carefully — this is the load-bearing framing)

**Step 7 is the capstone of the NUMERICAL-VALIDITY chain (brief §5 / pre-reg §J), not the physical
measurement.** It does exactly two things and claims exactly two things:

1. **State the numerical noise floor explicitly** — as a per-component breakdown (solver floor;
   discretization / self-convergence floor; boundary floor; conservation-closure floor), each traced to
   the step that established it, and the governing (largest) floor per observable class.
2. **Attach a quantitative uncertainty to every observable** — composed from those identified numerical
   error components by an explicit, justified combination rule, with a per-observable statement of which
   component dominates.

**What this step is NOT, and must say it is not (the honesty core — this is the step-6 lesson applied):**

- It is **NOT** the physical observable uncertainty of a frozen run. We are still in the **engineering-
  smoke** regime: placeholder parameters, **target-blind**, **no extraction map**, **no physical
  packet**, export guard False. There are no §G observables (`D0/D2/D4, N0/N2/N4, P0/P2/P4, m̂0, S_port,
  chi_Q, N_Q, Xi_1, varrho`) extracted yet, and there will be none until the user flips GATE A (freeze
  free_choice values target-blind) and GATE B (export guard). The budget is therefore **demonstrated on
  the same target-blind surrogate observables step 4 already produced** (the raw-field functionals: density
  mass integral, peak/min density, raw-field L2, A0 L2, field-energy-like integral, chemical potential,
  the two identically-zero spatial channels, final residual L∞). The **methodology** is the deliverable —
  it is what will assign each real observable a defensible uncertainty once extraction exists.
- It covers **NUMERICAL** error only — discretization, solver/Newton-GMRES floor, boundary truncation,
  conservation closure. It does **NOT** cover **physical-parameter / free_choice-value uncertainty**
  (placeholder params; GATE A — *none frozen anywhere*; this will be the DOMINANT uncertainty of the
  eventual real observables) or **model-form / parent-action uncertainty**. These omissions must be
  stated loudly and explicitly in the report; do not let a small numerical budget read as "the observable
  is known to this precision."

**Honesty constraints on the budget itself (overclaim guards — every one is a review axis):**

- **Floor it; never beat the floor.** `u(O) ≥ the solver floor scale` for every observable. An observable
  whose discretization error estimate sits *below* the solver/Newton floor is **floor-limited**: its
  uncertainty is the floor, NOT the smaller Richardson number. Reuse step 4's existing
  `"solver-floor diagnostic"` verdict — do not invent a smaller uncertainty for a quantity step 4 already
  flagged as floor-limited (density mass integral, final residual L∞).
- **No cherry-picking.** The budget for each observable must include **every** identified component that
  applies to it, not the smallest. State the combination rule.
- **Justify the combination rule; do not RSS correlated components as if independent.** If two components
  are independent error axes (grid spacing vs boundary truncation vs solver tolerance), root-sum-square is
  defensible — state that assumption. If two components are *not* plausibly independent for a given
  observable, use the conservative `max` instead and say why. Default to RSS of the independent axes with a
  `max`-fallback sanity check, and report both if they differ materially.
- **Apply relative floors to the right observable class.** The step-5 boundary floor and step-6
  conservation/Gauss floor are *relative* quantities; convert them to a per-observable absolute
  contribution honestly (`u_component(O) = relative_floor × |O|`), and apply the conservation/Gauss floor
  only to the charge/energy/Maxwell-coupled observables it actually bounds — null/identity for the others;
  do not silently smear one floor across observables it doesn't bound.
- **Null sectors stay null.** The identically-zero spatial gauge / spatial current channels get the
  existing `"null diagnostic"` label and **no** manufactured uncertainty.

**Contract.** Codex DESIGNS and WRITES the error-budget module, harness, report, and tests. Claude REVIEWS
code + output with clean agents + an independent arbiter re-run, and the standing transliteration-fidelity
audit runs on any new math (here the audit target is the *error-composition arithmetic + the
floor/relative-scaling logic*, since there is no new parent-action physics). This directive states
requirements + acceptance criteria only; the mechanism is your design call. Run everything LOCALLY ON CPU;
no GPU, no network. Never wrap your own session in a shell `timeout` (the `timeout 600` cap is only for the
short scripts you invoke). Still the **engineering-smoke** regime: placeholder params, **NO physical
packet**, export guard untouched, **target-blind**.

---

## 1. Objective

Produce, for the converged coupled stationary isotropic branch, a single machine-readable + tracked
**error budget** that:

- states the **numerical noise floor** explicitly, broken into its components (solver, discretization,
  boundary, conservation) with each component's source step and value;
- assigns **every surrogate observable a quantitative uncertainty** `u(O)` composed from the components
  that apply to it, by an explicit justified rule, floored at the solver floor;
- gives each observable a **dominant-component verdict** (which error source governs its uncertainty) and
  reuses step 4's null/floor verdict vocabulary where applicable;
- states plainly the **scope limits**: numerical-only; surrogate (not §G) observables; parameter/model
  uncertainty deferred to the frozen run (GATE A/B).

---

## 2. The error components and their sources (the spec — take the numbers from the prior steps, not here)

Compose `u(O)` from these identified numerical error components. Pull each from the step that established
it (re-run the module as a library call where cheap and deterministic, or consume its committed
machine-readable table with an explicit provenance tag — see §4):

- **Solver / Newton-GMRES floor** `u_solver` — the residual / mass-constraint stopping floor from step 4
  (`convergence.run_step4` → `noise_floor.preliminary_numerical_floor`, `solver_residual_floor_linf`,
  `mass_constraint_floor`; step-4 report: residual L∞ floor ≈ `1.07e-8`). This is the absolute hard floor:
  no observable's uncertainty may be quoted below it.
- **Discretization / self-convergence** `u_disc` — the per-observable finite-grid error: step 4's
  `finest_error_estimate` (Richardson, from `richardson_estimate` + `observed_order_from_three`) per
  observable, and the raw-field self-convergence floor (step-4 report: finest relative self-difference
  ≈ `3.48e-4` at 30721 DOF). For floor-limited observables this estimate is *below* `u_solver` — then the
  observable is solver-floor-limited (see honesty constraints).
- **Boundary truncation** `u_bdy` — the interior boundary-induced error from step 5
  (`boundary_characterization.run_step5`; step-5 report: max interior relative L2 boundary error
  ≈ `6.85e-5`, the real-sponge result). A *relative* floor → apply as `u_bdy(O) = boundary_rel × |O|`.
- **Conservation closure** `u_cons` — the Gauss-law / conservation closure floor from step 6
  (`conservation_diagnostics.run_step6`; step-6 report: independent Gauss-closure relative residual at the
  finest step-6 level, order-converging). A *relative* floor bounding the charge / energy / Maxwell-coupled
  observables → apply `u_cons(O) = cons_rel × |O|` only to those; null/identity elsewhere.

These components are numerical-validity quantities derived from the solver's own convergence behaviour.
They touch **no** §H target and **no** §F extraction coefficient (see §5).

---

## 3. Scope — IN (mandatory)

**A. Component noise-floor table (REQUIRED).**
A single explicit statement of the numerical noise floor as the four components above: name, value,
units (absolute vs relative), and the source step + commit for each. Identify the governing (largest)
floor for each observable class (mass/density; gauge/charge; energy; solver-floor diagnostics).

**B. Per-observable error budget (REQUIRED — the deliverable).**
For each of step 4's surrogate observables, a row giving: the finest-grid value; each applicable component
contribution (`u_solver`, `u_disc`, `u_bdy`, `u_cons`) as an *absolute* number on that observable; the
composed total `u(O)`; the relative uncertainty `u(O)/|O|` where `|O|` is non-zero; the dominant component;
and the reused verdict label (`"expected-order convergence"` / `"solver-floor diagnostic"` /
`"null diagnostic"`). The composition rule must be the one declared in §2 / the honesty constraints, with
`u(O) ≥ u_solver` enforced and floor-limited observables flagged.

**C. Combination-rule statement + sensitivity (REQUIRED).**
State the combination rule (RSS of independent axes, floored at `u_solver`) and its independence
assumption explicitly in the report. Also report the conservative `max`-of-components alternative for each
observable; if RSS and max differ materially for any observable, surface it rather than hide it.

**D. Scope-limits statement (REQUIRED — the honesty core).**
An explicit paragraph: this is a *numerical* budget on *target-blind surrogate* observables; it excludes
free_choice/physical-parameter uncertainty (GATE A — none frozen) and model-form uncertainty; the real
§G-observable budget is produced only after extraction + GATES A/B. No precision is claimed below the
solver floor.

**E. Machine-readable table + tracked report** (see §6).

---

## OUT — explicitly deferred / forbidden

- **No physical-parameter / model-form uncertainty.** Numerical-only; state the deferral.
- **No §G-observable extraction, no extraction map, no real targets.** Surrogate observables only;
  methodology transfers to the frozen run later.
- **No new heavy solves required.** Consume step 4/5/6 outputs (re-run reduced/deterministic ladders or
  read their machine-readable tables with provenance — your call, §4). Do not run the full step-4 5-level
  ladder *and* the full step-5 sweep *and* step 6 inside one script if that breaches `timeout 600`.
- **No re-derivation of a parallel error estimator.** Reuse `convergence.py`'s
  `richardson_estimate` / `observed_order_from_three` / `finest_error_estimate` / verdict vocabulary.
- **No physical packet, no export.** Export guard untouched (`physical_export_permitted` stays False).
- **No new dependencies.** torch + numpy (+ sympy if needed) only — already in the stack.
- **Never write under `research/pde_audit/simulation/`.** All run output under the gitignored
  `runs/`/`data/`/`figures/`/`_scratch/`; only code + tests + the small report + this directive are tracked.

---

## 4. Reuse requirements (fidelity — compose, do NOT reimplement)

- **Error estimates / orders / verdicts:** use `convergence.py` (`run_step4`,
  `observed_order_from_three`, `richardson_estimate`, the `finest_error_estimate` field, the
  `"null diagnostic"` / `"solver-floor diagnostic"` / `"expected-order convergence"` labels). Do **not**
  invent a second error estimator or a parallel verdict string set.
- **Floors:** read the solver floor from step 4's `noise_floor` block; the boundary floor from
  `run_step5`'s result; the conservation/Gauss floor from `run_step6`'s result. If you consume a committed
  JSON table instead of re-running, tag it with a **provenance block** (source report path + commit hash)
  so the input is auditable and the budget is reproducible. Prefer the smallest deterministic re-run that
  stays under `timeout 600`; pin any recorded scalar as a documented constant with provenance, not a
  magic number.
- **Observable set + labels:** the same surrogate observables and `OBSERVABLE_LABELS` from
  `convergence.py`; do not add or rename observables.
- **Config / harness pattern:** mirror the step-4/5/6 harness + report shape (a `run_step7` library
  function + a thin `error_budget_harness.py` CLI writing a machine-readable table under `runs/…`).
- Reproducible to the digit on re-run; a diagnostics digest like step 6 is encouraged.

---

## 5. Firewall / target-blindness (the #1 review axis)

- The error budget uses only solver-convergence quantities (orders, Richardson errors, floors) and
  target-blind surrogate observables. It must touch **no** §H target (R_norm, R_pole, P2=P4=0, chi_Q, GR
  quadrupole, Λ₀) and **no** §F extraction coefficient (D0/D2/D4, N0/N2/N4, P0/P2/P4, chi_Q, N_Q, m̂0,
  S_port, Xi_1, varrho). Do not import, read, or compare against the extraction map or `benchmarks.py`
  target values. (`references.py`/`benchmarks.py` may hold targets — do not pull target numbers into the
  budget path.)
- All normalization / relative scales are **target-blind aggregate norms** (the observable's own value or
  an aggregate norm), never an extraction coefficient.
- Export guard stays False; emit no physical packet.

---

## 6. Deliverables

1. An error-budget module (your naming; e.g. `error_budget.py`) with a `run_step7(config)` library
   function + a thin harness (e.g. `error_budget_harness.py`) mirroring the step-4/5/6 pattern, writing a
   machine-readable table (JSON) under `runs/…` (gitignored).
2. `software/stage1_solver/reports/step7_error_budget.md` — the tracked result report (§7).
3. Tests in `tests/test_stage1_solver.py` (see §8).
4. Do **NOT** `git add` or commit — the orchestrator stages + commits after re-review.

---

## 7. Report requirements (`reports/step7_error_budget.md`)

- State the **scope framing** up front: numerical-only budget on target-blind surrogate observables;
  parameter/model uncertainty deferred to the frozen run (GATE A/B). No precision below the solver floor.
- The **component noise-floor table** (§3-A): the four components, values, abs/rel, source step + commit,
  governing floor per observable class.
- The **per-observable error-budget table** (§3-B): value, each component contribution, composed `u(O)`,
  relative `u(O)/|O|`, dominant component, reused verdict label, floor-limited flag.
- The **combination-rule statement + sensitivity** (§3-C): the rule, its independence assumption, the
  `max`-alternative, any material RSS-vs-max divergence surfaced.
- An explicit **honest-limitations** paragraph (§3-D): the numerical-only scope; surrogate-not-§G; the
  GATE-A/B parameter-uncertainty deferral; null-sector observables; no sub-floor precision.
- Reproduction command(s); confirm target-blind + export guard untouched + provenance of any recorded
  input.

---

## 8. Acceptance criteria

1. **Composition-faithful & non-reimplemented:** error estimates / orders / verdicts via `convergence.py`;
   floors read from step 4/5/6 outputs with provenance; no parallel estimator. (The
   transliteration-fidelity audit checks the composition + floor/relative-scaling arithmetic term-by-term.)
2. **Floor enforced:** every `u(O) ≥ u_solver`; floor-limited observables flagged and carry the floor, not
   a smaller Richardson number.
3. **No cherry-picking; rule justified:** every applicable component included per observable; the
   combination rule + independence assumption stated; the `max`-alternative reported.
4. **Relative floors applied to the correct observable class:** boundary/conservation relative floors
   converted per-observable and applied only where they bound; null elsewhere.
5. **Scope stated honestly:** numerical-only; surrogate-not-§G; parameter/model uncertainty deferred;
   null sectors labelled; no sub-floor precision.
6. **Target-blind:** no target / extraction-coefficient access; export guard untouched; no new deps.
7. **`pytest` passes** (expect the prior count + the new test(s)); the error-budget harness **exits 0**;
   results reproduce to the digit on a re-run (only wall-clock may differ).
8. No script runs > 10 min (`timeout 600` on each script you invoke); a timeout is a failure to
   reformulate, not a cap to raise.

---

## 9. Tests (§8.3 detail)

Add focused, **non-tautological** tests to `tests/test_stage1_solver.py`:
- **Floor enforcement (load-bearing):** feed an observable whose discretization estimate is below the
  solver floor and assert the reported `u(O)` equals the floor (not the smaller estimate) and the
  floor-limited flag is set.
- **Monotonic composition:** assert that adding a non-zero component never *decreases* `u(O)` (guards
  against an accidental `min`/cherry-pick); and that RSS ≥ each single component and `max` ≤ RSS ≤ sum.
- **Relative-floor scope mapping:** assert the conservation/Gauss relative floor contributes to a
  charge/energy-coupled observable and contributes **zero** to an observable it does not bound (guards a
  silently mis-wired smear).
- **Null-sector guard:** assert an identically-zero observable yields `"null diagnostic"` and no
  manufactured (non-zero) uncertainty.
- **Determinism:** the budget table reproduces to the digit across two runs (digest), like step 6.

**When done**, report: the new module + harness names; the component noise-floor table (the four values +
sources); the per-observable budget headline (each observable's `u(O)`, relative uncertainty, dominant
component, verdict); the combination rule chosen + the RSS-vs-max sensitivity; the new test names + what
each asserts; confirmation `pytest` passes (state the count) and the harness exits 0 and reproduces to the
digit; and confirmation target-blindness + export guard are intact and any recorded input is
provenance-tagged.
