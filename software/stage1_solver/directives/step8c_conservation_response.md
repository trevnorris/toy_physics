# BUILD DIRECTIVE — Step 8c (CPU): current-carrying conservation + low-frequency SURROGATE response-coefficient methodology

> **Scope confirmed by the Claude+Codex `S_η` read-only consult (2026-06-14).** Determination + the two folded-in
> guardrails are recorded in `software/stage1_solver/decisions/02_step8c_s_eta_scope.md`. Conceptual flag CLEAR;
> proceed surrogate-only.

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 6 ("Conservation diagnostics — mass/charge/energy
drift over each run") + step 8 ("WP3 — P2 tangent on the converged stationary branch") + §2.2 row "The grouped P2
response operator's low-frequency expansion (`R_pole`, `R_norm`)". Brief §5.4 (conservation diagnostics) + §5.5
(error budget). Prereg `docs/stage1_preregistration.md` §J.4 (conservation) + §G (observables) + §H (the WP3 tangent
targets `d ln R_tr = d ln R_target = d ln ε_η = 0`, `N_Q = 1`). **Canonical PDE source:**
`notes/moving_throat_pde_program_compact.md` §4.7 lines 1383-1455 (the linearized hierarchy + the wall-equation RHS
`S_η^(ψ)+S_η^(A)+f_ext`), §4.6 lines 1377-1381 (the "still open" list), §13.4 lines 6846-6856 (the 5PN actual-branch
WP3 card).

**Builds on:** step 6 (`4a03797`) — the conservation machinery (`conservation_diagnostics.py`/`conservation_harness.py`:
`_matter_number_current`, the conservative FV face-flux/flux-divergence helpers in `operators.py`, the INDEPENDENT
center-gradient Gauss-law reconstruction, the Noether energy flux `S=μ·j`, the sponge accounting, the
`identity_checks` vs physics-gate separation); step 8a (`f67c9bc`) — the static ℓ=2 tangent (`p2_tangent.py`:
`p2_surrogate_forcing`, `p2_response_observables`); step 8b (`be39ba6`) — the DRIVEN ω≠0 tangent
(`p2_driven_absorber.py`: `apply_p2_driven_tangent`, the complex128 solve path, the raw response-vs-ω diagnostic
table, the CAP, `run_step8b`/`write_step8b_report`); step 4 (`c9b8b2c`) — the convergence/Richardson surrogate
machinery; step 7 (`470ea52`) — the error-budget composition (`error_budget.py`). **REUSE all of it.** Do NOT fork the
residual, the Newton/linear-solve core, the preconditioner, the conservation operators, or the driven tangent.

**Methodology scope (load-bearing).** Step 6 reported conservation as **null-by-symmetry** on the static isotropic
branch (no spatial transport: real amplitude ⇒ no phase winding ⇒ current/spatial-gauge vanish) and EXPLICITLY
deferred the genuinely current-carrying test to step 8 (`reports/step6_conservation_diagnostics.md` L36, L119). 8b's
driven ω≠0 tangent is **complex** — it carries a genuine non-null δJ (phase structure + a spatial-gauge perturbation).
8c is therefore the first place the conservation diagnostics have a real transport to measure. Separately, §2.2 / the
WP3 card frame the response operator's **low-frequency expansion** as the WP3 deliverable; 8c builds + validates the
**extraction methodology** on target-blind surrogate functionals (8b deferred any fit; its
`raw_response_vs_omega_not_fitted_not_a_physics_gate` asserted check is the seam 8c picks up). Still the
**engineering-smoke** regime: placeholder params, target-blind surrogate forcing, NO physical packet, export guard
untouched, no hard caps.

**The hard scope boundary (the predicted `S_η` STOP — RESOLVED to "defer", user-confirmed 2026-06-14).** The
**physical** WP3 tangent targets `d ln R_tr = d ln R_target = d ln ε_η = 0` are the induced-P₂₂ response — the
SELF-CONSISTENT matter/gauge ⇄ wall coupling. 8a/8b implemented the wall→matter direction (`δV_conf`, the `C_η[η]`
term). The return direction — the matter/gauge→wall source `S_η^(ψ)+S_η^(A)` (§4.7 L1449) — is **OPEN**: `S_η^(A)` has
no formula anywhere; `S_η^(ψ)` is schematic only; and §4.6 lists "the full coupled matter/gauge renormalization of
these reduced lanes" + "the outgoing odd response and quadrupole-normalization branch" as still open. Completing it is
fresh derivation work = **Path-A material** (deferred per the frozen `parent_action_status = effective_closure`
decision, `docs/branch_realization_parent_status_decision.md`). **Therefore the physical `d ln R_tr / R_target / ε_η`
are BLOCKED in the effective-closure scope and 8c does NOT attempt them** — it does the SURROGATE methodology only and
DOCUMENTS the deferral.

**Reachability split (consult-corrected — do NOT over-credit the secondaries):**
- **Truly S_η-independent (ω⁰ / source-map, WP1-extractable) — unaffected, carry the falsification load:** the PRIMARY
  `R_norm`, plus `chi_Q` and `N_Q = chi_Q⁻¹`.
- **Low-frequency-response-gated (NOT safe-WP1):** `R_pole`, `P_2`, `P_4` depend on the ω²/ω⁴ response-bundle
  coefficients `D2/D4/N2/N4` (compact §7A.1 L2677-2710); their physical extraction needs the field→coefficient
  extraction map (TODO) and is at minimum extraction-map-gated, potentially S_η-gated. **8c does NOT extract them.**
- **Induced self-consistent (explicitly S_η-blocked):** `d ln R_tr / R_target / ε_η`. **Deferred to Path A.**

**Do NOT invent `S_η`.** **Biggest trap (consult): LABEL DRIFT** — surrogate coefficients, closed-conservation
language, or partial-`S_η` structure must NEVER be presented as physical WP3 target extraction.

**Contract.** Codex DESIGNS and WRITES the current-carrying conservation diagnostics on the driven solution, the
low-frequency surrogate response-coefficient fit + its validation, the report, and the tests. Claude REVIEWS code +
output with clean agents (a mandated transliteration-fidelity audit on any new math + an adversarial
overclaim/target-blindness/can't-fail-gate review — with special attention to the step-6 "is this independent or is it
x−x?" diagnostic-design discipline and to fit non-circularity) + an independent arbiter re-run, and iterates with you
until acceptance. This directive states requirements + acceptance only; the mechanism is your design call. Run
everything LOCALLY ON CPU; no GPU, no network. Every script you RUN stays wrapped in `timeout 600`. Never wrap your own
Codex session in a shell `timeout`.

---

## 1. Objective

On the 8b **driven (ω≠0) current-carrying** P2 tangent, (A) measure genuine **current-carrying conservation**
diagnostics — the transport that step 6 reported null-by-symmetry — as an HONEST source-balanced statement (the driven
tangent is forced + absorbed, so the conserved-current statement is `∂_M J^M = injected − absorbed`, measured and
shown to close, NOT a closed conservation dressed up), reusing the step-6 machinery; and (B) build + validate a
**low-frequency response-coefficient EXTRACTION METHODOLOGY** on **target-blind surrogate functionals** of the driven
response — the fit, its convergence under refinement, and a per-coefficient uncertainty (consuming step 7) — as the
methodology the eventual physical run would apply. **Explicitly DEFER** the physical WP3 tangent targets
(`d ln R_tr / R_target / ε_η`) to Path A (the open `S_η^(ψ,A)`), and record that deferral for the step-9 verdict. NO
physical packet, NO target, export guard untouched.

---

## 2. Scope

### IN — mandatory

**A. Current-carrying conservation on the driven solution.** Run the step-6 conservation diagnostics on the 8b driven
ω≠0 tangent solution (reuse `conservation_diagnostics.py`):
- the **matter/gauge number + charge current** — now GENUINELY NON-NULL (verify the driven δJ is non-zero; if it is
  null, that is itself a finding to report, not to hide). **⚠️ CONSULT GUARDRAIL: do NOT reuse `_matter_number_current`
  as-is** — its formula is for the REAL stationary branch fields (`coupled_branch.py:210`) and is WRONG on the complex
  phasor lanes of the driven tangent. **Build a properly LINEARIZED PHASOR current around the WP1 background** (the
  correct time-averaged / phasor matter+gauge current for the driven complex perturbation `δψ, δA` about `ψ0, A0`);
  charge `= q·number`. A clean transliteration-fidelity audit on this new current construction is mandatory (this is
  exactly the faithful-but-wrong-operator class the audit exists to catch).
- the **Noether energy flux `S = μ·j`** on the driven solution;
- an **INDEPENDENT** Gauss-law closure (center-gradient `E=−∇A0` reconstruction — a DIFFERENT stencil from the
  solver's own flux operator; this is the step-6 independence discipline, do not regress to the integrated-Newton-
  residual tautology that step-6's adversarial review caught);
- the **honest source-balance framing.** Because the tangent is DRIVEN (the surrogate `f_ext` injects) and ABSORBED
  (the CAP/sponge removes), the statement is a BALANCE `transport divergence = injected source − absorbed sink`, not a
  closed `∂_M J^M = 0`. Measure each term (injection from `f_ext`, absorption from the CAP/sponge — mirror step-6's
  sponge accounting), show the balance closes to the discretization floor, and label it a *balance with measured
  source/sink*, NOT a closed conservation. Keep the step-6 `identity_checks` vs physics-gate `_not_a_physics_gate`
  separation (a budget that closes by FV construction is an identity, not a physics gate).
- a **convergence read** (≥3 levels) of the non-null transport residuals + the Gauss closure (step-6 reported
  order ~2–4 for the independent Gauss closure; current-carrying residuals should converge or be honestly labelled).

**B. Low-frequency surrogate response-coefficient methodology.** Build + validate an extraction of a low-frequency
expansion from the driven response:
- consume the 8b raw response-vs-ω diagnostic (or recompute via `run_step8b`'s machinery) over a **low-ω band chosen
  on numerical grounds** (well-resolved — avoid the ω=6 resolution-limited regime; the low-freq expansion lives where
  resolution is clean) with enough ω samples to fit;
- fit a low-frequency expansion (e.g. `O(ω⁰)+O(ω²)+...` per the relevant parity, or a low-order
  Taylor/Padé as you justify) of **target-blind SURROGATE functionals** of the solved response `δu(ω)` — NOT
  `R_pole`/`R_norm`/`D2`/`D4`/`N2`/`N4`, NOT `d ln R_tr/R_target/ε_η`, NOT any §H target;
- **validate the methodology**: the extracted surrogate coefficients must CONVERGE under grid refinement (a coefficient
  that drifts under refinement is not a measurement — §J.2 discipline) AND be stable under reasonable changes to the
  ω-sampling set; attach a **per-coefficient uncertainty** by composing the step-7 error-budget axes (solver /
  discretization-Richardson / boundary / conservation) — extend `error_budget.py`'s methodology, do not re-pin floors
  by hand if reusable;
- state precisely what the surrogate coefficients ARE (functionals of the driven solution), that they are NOT the
  physical response coefficients, and how the SAME machinery would extract the physical ones once GATE A/B + the
  extraction map + (for the WP3 sector) the open `S_η` exist.

**C. Honest scope + the `S_η` deferral record.** The report MUST contain a clearly-labelled **"WP3 physical-target
reachability"** section stating: the physical `d ln R_tr / R_target / ε_η` require the OPEN matter/gauge→wall source
`S_η^(ψ,A)` (§4.7 RHS / §4.6 "full coupled matter/gauge renormalization"), which is Path-A material per the frozen
`effective_closure` decision; therefore they are DEFERRED (the step-9 verdict reports them blocked/deferred, NOT
pass/fail). State the consult-corrected reachability split: the PRIMARY `R_norm` + `chi_Q` + `N_Q = chi_Q⁻¹` are
S_η-independent (ω⁰/source-map, WP1-extractable) and carry the intact falsification load; `R_pole`/`P_2`/`P_4` are
low-frequency-response-gated (depend on `D2/D4/N2/N4`), NOT safe-WP1, and are reported extraction-map-gated /
potentially S_η-gated — not over-credited. Cite `docs/branch_realization_parent_status_decision.md` + the Claude+Codex
`S_η` consult record (`software/stage1_solver/decisions/02_step8c_s_eta_scope.md`).

**D. Certification + determinism.** Keep the **pass_checks / asserted_checks `_not_a_physics_gate` separation**
(`passed = all(pass_checks.values())`); the deterministic **diagnostics digest** (reproduces byte-for-byte on re-run);
a manifest per solve; a machine-readable 8c table. Any NEW operator/diagnostic term that enters a convergence claim is
MMS-certified or convergence-tested (the conservation diagnostics reuse already-certified operators; the fit is a
post-processing step validated by refinement-convergence, not MMS).

### OUT — do not touch this step

- **Inventing / positing `S_η^(ψ)` or `S_η^(A)`** — the matter/gauge→wall source stays OPEN; drive with the
  target-blind surrogate `f_ext` only (as 8a/8b). The physical WP3 tangent targets stay DEFERRED to Path A.
- **Any §H target or the physical response coefficients** — `R_pole`, `R_norm`, `D2/D4/N2/N4`, `P_2/P_4`,
  `d ln R_tr/R_target/ε_η`, `N_Q`, `chi_Q`, the induced-`P_{22}` numbers. Do not name, compute, fit-to, or compare
  against any of them. (The SURROGATE fit extracts functionals that are NOT these.)
- The **field→coefficient extraction map** (still TODO, gated); the **frozen physical run** / **free_choice freeze
  (GATE A)** / **export-guard flip (GATE B)** (all gated, USER calls); GPU / PETSc port.
- The **exact ℓ=2 vector-harmonic Maxwell reduction** (8a F2 / 8b temporal-truncation deferral — derivable, deferred);
  the conservation diagnostics inherit the 8b engineering-smoke scalarization + temporal-truncation caveats — keep them
  labelled, do not relabel as exact.
- **Hard caps are FORBIDDEN** (prereg §D); do NOT write under `research/pde_audit/simulation/`; do NOT import / modify /
  satisfy / trip the export-guard scripts (`physical_export_permitted` stays asserted-by-construction, NOT made "real"
  by importing the firewalled model).

---

## 3. Requirements

**R1 — Reuse, don't fork.** Conservation: extend/reuse `conservation_diagnostics.py` + the `operators.py` FV helpers +
the 8b driven solve path. Response fit: reuse the 8b response machinery + the step-4 convergence/Richardson helpers +
the step-7 error-budget composition. Mirror the `run_step8b`/`write_step8b_report` structure for
`run_step8c`/`write_step8c_report`. Do not duplicate the residual, the Newton/linear-solve core, the preconditioner,
the conservation operators, or the driven tangent.

**R2 — Genuinely current-carrying (teeth).** Demonstrate the driven solution carries a NON-NULL current/transport
(contrast against step-6's null-by-symmetry static result — the analog of the 8b reflecting-control teeth: the test is
only credible if the static branch gives null and the driven branch gives non-null). State the contrast.

**R3 — Independent, non-tautological diagnostics (the step-6 discipline).** The Gauss closure must use a stencil
INDEPENDENT of the solver's own fluxes (center-gradient reconstruction). Any budget that closes by FV construction is
an `identity_check`, NOT a physics gate. The source-balance must measure injection (`f_ext`) + absorption (CAP/sponge)
as DISTINCT measured quantities, not algebraically cancel them into `x−x`. No hardcoded nulls; the energy flux must be
a genuine `μ·j` measurement, not a discarded-`μ` zero (the exact sins step-6's adversarial review caught — do not
repeat them).

**R4 — Non-circular surrogate fit (consult guardrail).** The fitted surrogate functionals must be **predeclared /
fixed before the runs**, **target-blind**, and **independent of the forcing overlap** — labelled CAP/operator-
methodology coefficients only. The fit maps forcing→solve→functional and extracts the ω-expansion of the SOLVED
response; it may NOT fit a functional that is trivially proportional to / an overlap with the forcing (that would make
the "coefficient" a tautology of the input). State what each functional is, why it is a genuine response (not an echo
of `f_ext`), that it was declared before the runs, and how the fit was validated (refinement-convergence + ω-sampling
stability + uncertainty). Report the honest result even if a coefficient drifts/under-resolves — a drifting coefficient
is "not a measurement", reported, not hidden.

**R5 — No fabricated success (integrity axis).** Do not tune the window/frequency band/metric/sampling to manufacture
clean convergence or a closed balance. Do not loosen Newton/GMRES/linear-solve tolerances. Report any
non-convergence/non-closure. Do not certify any term not covered by MMS/convergence. The convergence/closure gates
must be genuinely computed (the recurring can't-fail-gate sin: no hardcoded-`True` checks counted in `passed`; no
cherry-picked observable subset; no `min_order=0`).

**R6 — TARGET-BLINDNESS (the #1 review axis).** No §H target value anywhere — not in code, config, params, ω band,
fit functionals, thresholds, or the conservation references. ω band + fit choices are numerically motivated, never
target-driven. Verifiable by grep + the AST/regex import+literal scan in `test_p2_solve_path_is_target_blind` — extend
it to cover the new 8c modules.

**R7 — Non-frozen smoke discipline / guard firewall intact.** Placeholder parameters only, labelled engineering-smoke /
target-blind. No physical packet; no response-coefficient/DtN/`chi_Q` extraction; do NOT write under
`research/pde_audit/simulation/`; do NOT import / modify / satisfy / trip the export-guard scripts. No hard caps.

**R8 — Environment & deps.** Self-contained on the existing stack (torch [incl. complex128] / numpy / scipy). If you
believe a genuinely new dependency is needed, STOP and FLAG it. Every audit/benchmark script you RUN stays wrapped in
`timeout 600`; a timeout (exit 124) is a failure → reformulate, never raise the cap.

---

## 4. Acceptance criteria

1. **Current-carrying conservation** on the 8b driven solution: non-null number/charge current + energy flux `S=μ·j`
   + an INDEPENDENT Gauss closure, framed as an HONEST source-balance (`divergence = injected − absorbed`, each term
   measured) that closes to the discretization floor, with the step-6 `identity_checks` vs physics-gate separation, a
   convergence read (≥3 levels), and the non-null-vs-static teeth (R2).
2. **Low-frequency surrogate response-coefficient methodology**: a validated fit of target-blind surrogate functionals
   over a well-resolved low-ω band — converges under refinement, stable under ω-sampling, with a per-coefficient
   uncertainty composed from the step-7 axes; explicitly NOT the physical response coefficients; non-circular (R4).
3. **The `S_η` deferral recorded**: a "WP3 physical-target reachability" report section stating the physical
   `d ln R_tr/R_target/ε_η` are blocked by the open `S_η^(ψ,A)` (Path-A material) and deferred for the step-9 verdict,
   with the consult-corrected split (`R_norm`/`chi_Q`/`N_Q` S_η-independent/intact; `R_pole`/`P_2`/`P_4`
   low-freq-response-gated, not safe-WP1, not over-credited), citing the parent-status decision + the consult record.
4. **Target-blind** (R6, grep- + test-verifiable); guard firewall intact; non-frozen smoke discipline; no hard caps;
   no response-coefficient extraction of any §H target; the **OPEN `S_η^(ψ,A)` NOT invented**.
5. **Honest scope/labels**: the conservation diagnostics inherit + keep the 8b scalarization + temporal-truncation
   caveats (not relabeled exact); the surrogate fit is labelled methodology-only.
6. **Determinism + manifest + machine-readable table**; the diagnostics digest reproduces byte-for-byte; the
   **pass_checks / asserted_checks `_not_a_physics_gate` separation** with `passed = all(pass_checks)`.
7. A concise **report** (markdown/YAML, KB-scale) under `reports/` capturing all of the above, with a forward note to
   step 9 (the Stage-1 verdict).

Claude reviews code + report with **two clean agents** — (i) the mandated **transliteration-fidelity audit** (any new
math: the current-construction / energy-flux / source-balance terms vs §4.7 + the step-6 source; the fit vs its stated
expansion) and (ii) an **adversarial overclaim / target-blindness / can't-fail-gate / diagnostic-independence review**
(the step-6 "independent or x−x?" lens + the fit-non-circularity lens + the non-null-current teeth) — plus an
**independent arbiter re-run** reproducing the 8c table + digest. We iterate with you until all hold. (MMS/arbiter
cannot catch a faithful-but-wrong operator, a tautological diagnostic, or a circular fit — the fidelity + adversarial
audits are the only checks that do; they caught 2 MAJOR in 8a, 5 in 8b, and 3 overclaims in step 6.)

---

## 5. Repo hygiene (firewall is live)

- **Code** → `software/stage1_solver/src/stage1_solver/`, **tests** → `software/stage1_solver/tests/` (tracked).
- **All run output** → `software/stage1_solver/{runs,data,figures}/` — **gitignored**.
- **The one tracked artifact** → a small report under `software/stage1_solver/reports/`.
- Do **not** commit data; do **not** `git add` anything (the orchestrator stages explicit paths + commits after
  review); **never** `git add -A`; **never** write under `research/pde_audit/simulation/`. If you update `README.md`,
  keep it purely additive (a step-8c row + run command) and flag it for review.

---

## 6. Deliverable for review

Report back with: the current-carrying conservation design (how the non-null driven current is constructed + verified
non-null vs the static null; the source-balance framing + each measured term; the independent Gauss closure; the
convergence read); the low-freq surrogate methodology (the ω band + why target-blind/well-resolved; the fit form; the
surrogate functionals + why genuine-response-not-forcing-echo; the refinement-convergence + ω-stability + per-coeff
uncertainty); the "WP3 physical-target reachability" / `S_η` deferral section; the pass/asserted check listing; the new
pytest count; harness exit 0 + the diagnostics digest reproduced; confirmation of target-blindness + firewall + the
asserted/counted separation + that `S_η` was NOT invented; the report path; and any STOP-and-flag items. Then Claude
reviews before anything is committed.
