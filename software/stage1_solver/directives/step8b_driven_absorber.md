# BUILD DIRECTIVE — Step 8b (CPU): driven (ω≠0) P2 tangent + outgoing-wave absorber + genuine reflection coefficient

**Plan ref:** `docs/branch_realization_execution_plan.md` §7 step 8 ("WP3 — P2 tangent on the converged
stationary branch"); §2.2 row "Outgoing DtN data via Robin impedance / PML on the linear tangent" and
the matching "does not capture: Long-time dynamic radiation propagation". Brief §5.3 (boundary control —
demonstrate reflections sit below the target signal). Boundary lit synthesis
`docs/boundary_and_noise_methods.md` **§0 decision of record** ("Defer the CAP/PML to step 8 … the
genuine outgoing-wave reflection coefficient is a property of the wave-like WP3 linearized tangent
(step 8); the proper absorber (CAP/PML) is designed and characterized there") and **§3** ("Step-8 plan:
generic smooth CAP inspired by Manolopoulos paired with the Robin exit; escalate to a self-adaptive higher-order
ABC or a (stationary-adapted) matter-wave PML if residual low-momentum reflection persists").
Step-5 report `reports/step5_boundary_characterization.md` **Deferral Note** ("The outgoing-wave
reflection coefficient is deferred to the step-8 linearized tangent because the stationary isotropic
branch has no propagating waves to reflect"). **Canonical PDE source:**
`notes/moving_throat_pde_program_compact.md` §4.7 lines 1383-1455 (the time-derivative structure:
matter `iℏ∂_t[δψ,δψ*]=L_BdG+C_A+C_η`, boxed Maxwell `∂_M(Z δF^{MN})+(1/ξ)∂^N(∂·δA)=μ₀δJ^N`, boxed wall
`μ_η ∂_t² η − ∂_w(T_w ∂_w η) − T_Ω Δ_{S²} η + K_η η = S_η^(ψ)+S_η^(A)+f_ext`).

**Builds on:** step 8a (`f67c9bc`) — the static (ω=0) grouped-real P2 tangent
(`p2_tangent.py`: `p2_mode_residual`, `apply_p2_tangent` [= JVP of the step-3 residual + the linear ℓ=2
angular terms], `solve_static_p2_tangent`, `p2_surrogate_forcing`, `p2_response_observables`,
`run_step8a`), the ℓ=2 operators (`operators.py`: `tensor_laplacian_with_spherical_l`,
`localized_maxwell_operator_with_spherical_l`), the P2 MMS pieces (`mms_benchmarks.py`), the colored
sparse Jacobian (`preconditioners.py`), and `P2TangentConfig`. **REUSE all of it.** Step 5 (`63cd885`)
provides the real soft sponge + the interior-window/difference-norm machinery
(`boundary_characterization.py`); step 4 (`c9b8b2c`) provides the convergence/Richardson surrogate
machinery. REUSE those too.

**Methodology scope (load-bearing).** 8a built the *static* ℓ=2 tangent — there are no waves at ω=0,
so it could only show interior insensitivity (step 5 did the stationary contamination study). 8b is the
step where **waves actually exist**: driving the linearized tangent at a nonzero frequency ω turns the
elliptic 8a operator into a wave-supporting (complex, indefinite) operator with **genuinely outgoing
solutions at the exit**. This is the step the plan and the boundary doc reserved the real absorber and
the real reflection coefficient for. Three deferred pieces converge here and are all in scope:
(1) the **driven ω≠0 tangent operator** (the temporal terms 8a set to zero); (2) the **real
outgoing-wave absorber** (CAP/PML — deferred from step 5); (3) the **genuine wave-reflection
coefficient** (deferred from step 5). Still the **engineering-smoke** regime: placeholder params,
target-blind surrogate forcing, NO physical packet, export guard untouched.

**The sponge×MMS gap goes LIVE here (read `docs/boundary_and_noise_methods.md` §3 + the step-8a report's
forward note).** 8a's sponge was a static additive term that never entered a certified differential
stencil. In 8b the absorber AND the new ω-frequency terms enter the operator that gets convergence-tested.
**Therefore every new operator term — the temporal/frequency terms and the absorber term — MUST be added
to the MMS continuum forcing and certified at the expected order BEFORE it is used in any
convergence/reflection claim.** A smooth CAP `−iσ(x)` is an additive zeroth-order complex
term (composes cleanly, easy to MMS). A PML modifies the differential operator (coordinate stretching) —
if you choose PML you MUST MMS-force the modified stencil and you should STOP-and-flag if it balloons
beyond a focused build.

**Contract.** Codex DESIGNS and WRITES the driven operator, the absorber, the reflection-coefficient
measurement, the MMS for the new terms, the report, and the tests. Claude REVIEWS code + output with
clean agents (a mandated transliteration-fidelity audit on the new math + an adversarial
overclaim/target-blindness review) + an independent arbiter re-run, and iterates with you until
acceptance. This directive states requirements + acceptance only; the mechanism (representation,
discretization, measurement) is your design call. Run everything LOCALLY ON CPU; no GPU, no network.
Every script you RUN stays wrapped in `timeout 600`. Never wrap your own Codex session in a shell
`timeout`.

---

## 1. Objective

On the converged stationary isotropic background, build the **driven (ω≠0) grouped-real P2 linearized
tangent** by adding the §4.7 time-derivative terms (matter `iℏ∂_t` BdG-doublet, Maxwell temporal part of
`∂_M(Z δF^{MN})`, wall `μ_η ∂_t² η`) that step 8a set to zero, making the operator wave-supporting at a
chosen drive frequency. Add a **genuine outgoing-wave absorber** at the exit (default: generic smooth
CAP inspired by Manolopoulos 2002, without claiming the analytic pole-profile guarantee) upgrading the step-5 Robin/sponge treatment, and **measure a genuine, non-
tautological wave-reflection coefficient** that demonstrably distinguishes the absorber from a reflecting
control and sits below a stated target-blind floor. Certify the new operator terms by MMS (the
sponge×MMS gap is live), gate-A internal consistency, well-posedness, and convergence — all on
**target-blind surrogate forcing**. This closes the boundary doc's step-8 absorber plan and the step-5
deferral. NO physical response-coefficient extraction, NO target, export guard untouched.

---

## 2. Scope

### IN — mandatory

**A. Driven (ω≠0) tangent operator.** Extend the static 8a tangent to a time-harmonic drive at frequency
ω by transliterating the §4.7 (lines 1406-1451) temporal terms the static build dropped:
- matter: the `iℏ∂_t[δψ,δψ*]` BdG-doublet term (in the frequency domain a `±ℏω`-type doublet term — get
  the sign/conjugate structure from the source, not from intuition);
- Maxwell: the temporal part of the boxed `∂_M(Z δF^{MN})` (and of `(1/ξ)∂^N(∂·δA)`) that contributes
  ω-dependent terms the static H=Z reduction omitted;
- wall: `μ_η ∂_t² η → −μ_η ω² η`.

The result is in general **complex and indefinite** (wave-supporting). The representation — complex128
fields vs a real BdG `(u,v)`-doublet / real 2-block — is your design call, but whatever you choose:
- the background-dependent (8a) Jacobian part must be **unchanged** (still the JVP of the step-3
  residual — re-verify gate-A on it);
- the new ω-terms must be **MMS-certified** (§D) and **fidelity-audited** (term-by-term vs §4.7);
- include **≥3 drive frequencies ω** chosen on **numerical grounds** (target-blind), with **at least one
  in the genuinely propagating regime** (an oscillatory interior solution whose wavelength is well
  resolved yet small enough to reach the exit as an outgoing wave — i.e. a real wave to absorb), and at
  least one near-static for a sanity limit (ω→0 should recover the 8a static operator; assert this).

**B. Genuine outgoing-wave absorber.** Implement a real absorber at the exit in the propagating
sector(s), upgrading step 5's Robin/sponge. **Default = generic smooth CAP** (in the spirit of
Manolopoulos 2002, but not claiming the analytic pole-profile guarantee): additive `−iσ(x)` with a single width parameter (`docs/boundary_and_noise_methods.md` §3, §5). A PML /
self-adaptive higher-order ABC is permitted **only** if you justify it and MMS-force the modified stencil
(§D); STOP-and-flag if a PML balloons beyond this focused build. The absorber profile parameters are
**conditioning / numerically-motivated free_choice**, never target-driven.

**C. Genuine wave-reflection coefficient (with teeth).** Measure a quantitative, **non-tautological**
reflection (or boundary-contamination-of-an-outgoing-wave) coefficient for the driven solution, and show:
- it sits **below a stated target-blind floor** (e.g. the step-4 discretization floor / the interior
  signal magnitude — NOT any §H target), AND
- a **reflecting control** (absorber OFF, or a deliberately reflecting truncation) **demonstrably shows
  materially higher reflection** — i.e. the metric has teeth and is not a quantity that passes
  regardless. This control is the analog of the 8a F1 discipline ("the test must FAIL on the wrong
  construction"): the metric is only credible if the reflecting case fails it.

The measurement method is your design call (candidates in `docs/boundary_and_noise_methods.md` §3-§4:
domain-doubling reference comparison / interior invariance as the truncation moves out; an analytic
Manolopoulos reflection bound only if the true wavelength-sized profile is implemented; an incident/reflected amplitude decomposition
where a clean wave region exists). Pick one as primary, state precisely what it measures, and justify why
it is not circular.

**D. MMS for every new operator term (the live sponge×MMS gap).** Add MMS coverage for the new
continuum terms — the ω-frequency terms (matter/Maxwell/wall) and the absorber term — manufacturing the
forcing from the **same continuum operator including those terms**, and certify each at the expected
order (~2), non-circularly, on a refinement ladder (mirror the existing `run_p2_centrifugal_mms` /
`run_p2_localized_maxwell_angular_mms` pattern). No new term may enter a convergence/reflection claim
before it is in the MMS forcing.

**E. Certification + honest scope.** Retain/extend: the **gate-A FD check** (state explicitly it
verifies the assembled operator IS the JVP of the driven mode-residual — *internal consistency only*;
the ω/absorber *physics* is certified by MMS §D + the fidelity audit, NOT by gate-A); a **well-posedness
check** generalized to the complex/indefinite driven operator (e.g. smallest singular value / a
defensible conditioning measure for the complex system); **convergence** of the driven solve + the
surrogate response observables; the 8a **asserted_checks / pass_checks `_not_a_physics_gate` separation**
(`passed = all(pass_checks.values())`); and the deterministic **diagnostics digest**. Add a raw
surrogate **response-magnitude-vs-ω diagnostic table** (target-blind functionals of `δu` at each ω) —
as a diagnostic only.

### OUT — do not touch this step

- The **OPEN matter/gauge→wall source `S_η^(ψ,A)`** (compact lines 1377-1381 / §4.7). 8b drives with the
  **target-blind surrogate `f_ext`** (reuse/extend `p2_surrogate_forcing`), exactly as 8a did. Do **not**
  invent the physical source — it stays deferred to 8c, where it will likely force a Claude+Codex
  methodology decision or a user-scope call.
- **Response-coefficient extraction / low-frequency-expansion fitting** — `R_pole`, `R_norm`, `D2`, `D4`,
  `N2`, `N4`, the induced-`P_{22}` response numbers (§F extraction / §H targets, firewalled). A raw
  response-vs-ω diagnostic (E) is allowed; **fitting it to extract expansion coefficients is OUT** (→ 8c
  surrogate / the gated extraction map). Do not name or compute any §H target.
- The **exact ℓ=2 vector-harmonic Maxwell reduction** (8a F2 deferral — derivable but deferred). The new
  ω² Maxwell term **inherits the same engineering-smoke scalarization caveat** (exact for the scalar
  lanes A0/Aw; a componentwise scalarization for the radial vector lane Ar). Keep that label honest;
  do not relabel the scalarization as "Pinned"/exact.
- **Current-carrying conservation** (→ 8c); any §H target; the **field→coefficient extraction map**; the
  **frozen physical run / free_choice freeze (GATE A) / export-guard flip (GATE B)** (all gated); GPU /
  PETSc port. **Hard caps are FORBIDDEN** (§D of the pre-reg) — never clamp a field to a hard value to
  fake absorption.

---

## 3. Requirements

**R1 — Reuse, don't fork.** Build the driven operator as an extension of `p2_mode_residual` /
`apply_p2_tangent` (add the ω-terms + absorber as new terms / a frequency parameter), reuse
`solve_static_p2_tangent`'s solve path (generalized to the complex/driven system), `p2_surrogate_forcing`,
`p2_response_observables`, the step-5 sponge + interior-window/difference-norm helpers
(`boundary_characterization.py`), the convergence/Richardson machinery, the colored sparse Jacobian
preconditioner (extend to complex if needed), `manifest.py`, `grid.py`. Mirror the `run_step8a` /
`write_step8a_report` structure for `run_step8b` / `write_step8b_report`. Do **not** duplicate the
residual, the Newton/linear-solve core, or the preconditioner.

**R2 — Same operator except the piece under test.** In the reflection study, change ONLY the absorber
(on/off, profile, width) — drive frequency, grid, placeholder params, gauge, and all other BCs stay
fixed, so the measured reflection change is attributable to the absorber. In the convergence study,
change only the grid. State this.

**R3 — Defensible, non-tautological reflection metric (R-with-teeth).** The reflection metric must be
measure-consistent (conservative r² cell volumes), identical across the on/off settings, and accompanied
by the reflecting control (C) that demonstrably fails / shows materially higher reflection. State the
metric, what it measures, the target-blind floor, the achieved ratio, and the control contrast.

**R4 — No fabricated success (integrity axis).** Report the honest reflection even if the absorber
underperforms — an absorber that does not beat the control is a real finding to diagnose, not to hide. Do
not tune the measurement window / frequency / metric to manufacture a low number. Do not loosen
Newton/GMRES/linear-solve tolerances. Report any non-convergence; do not drop a failing setting silently.
Do not certify any term not covered by MMS (D).

**R5 — Determinism & manifest.** Deterministic (fixed seed, single-thread, float64 / complex128). Every
solve emits a manifest (grid, solver controls incl. preconditioner, drive frequency, absorber settings).
Emit a machine-readable 8b table. The diagnostics digest must reproduce byte-for-byte on a re-run.

**R6 — TARGET-BLINDNESS (the #1 review axis).** No §H target value anywhere — not in code, config,
params, drive frequencies, absorber widths, reflection thresholds, the reflection floor, or the surrogate
observables. Drive frequencies and absorber widths are numerically/dispersion-motivated, never
target-driven. The reflection floor reference is target-blind (step-4 floor / interior signal magnitude).
Verifiable by grep, and by the AST/regex import+literal scan already in
`test_p2_solve_path_is_target_blind` — extend that test to cover the new driven/absorber modules.

**R7 — Non-frozen smoke discipline / guard firewall intact.** Placeholder parameters only, labelled
engineering-smoke / target-blind. No physical packet; no response-coefficient/DtN/`chi_Q` extraction; do
NOT write under `research/pde_audit/simulation/`; do NOT import / modify / satisfy / trip the export-guard
scripts (`physical_export_permitted` stays asserted-by-construction, NOT made "real" by importing the
firewalled model); no frozen-packet schema. No hard caps.

**R8 — Environment & deps.** Self-contained on the existing stack (torch [incl. complex128] / numpy /
scipy [scipy sparse LU handles complex]). If you believe a genuinely new dependency is needed, STOP and
FLAG it (reason + tradeoff); do not silently install. Every audit/benchmark script you RUN stays wrapped
in `timeout 600`; a timeout (exit 124) is a failure → reformulate, never raise the cap.

---

## 4. Acceptance criteria

1. A **driven ω≠0 tangent operator** transliterated from §4.7 (matter `iℏ∂_t` doublet, temporal Maxwell,
   wall `μ_η∂_t²`), with the background Jacobian unchanged from 8a (gate-A re-verified), ≥3 drive
   frequencies incl. ≥1 genuinely propagating + a near-static ω→0 sanity limit that recovers the 8a
   static operator.
2. A **genuine outgoing-wave absorber** (default generic CAP inspired by Manolopoulos) upgrading the step-5 exit treatment,
   with stated conditioning-motivated parameters.
3. A **genuine, non-tautological wave-reflection coefficient** below a stated **target-blind floor**,
   **AND** a reflecting control that demonstrably shows materially higher reflection (teeth, R3/C).
4. **MMS coverage for every new operator term** (ω-frequency terms + absorber term), certified at the
   expected order non-circularly (the live sponge×MMS gap closed, §D).
5. **gate-A** (internal-consistency, with the explicit caveat it does not certify the ω/absorber physics),
   **well-posedness** (complex operator), **convergence** of the driven solve + surrogate observables,
   and the **pass_checks / asserted_checks `_not_a_physics_gate` separation** with `passed=all(pass_checks)`.
6. **Target-blind** (R6, grep- + test-verifiable); guard firewall intact; non-frozen smoke discipline; no
   hard caps; no response-coefficient extraction; the **OPEN `S_η^(ψ,A)` source NOT invented** (surrogate
   forcing only).
7. Honest **scope/labels**: the new ω² Maxwell term carries the 8a F2 engineering-smoke scalarization
   caveat (not relabeled "Pinned"); an explicit note that the physical response-coefficient extraction +
   the OPEN wall source remain deferred (8c / gated).
8. **Determinism + manifest + machine-readable table**; the diagnostics digest reproduces byte-for-byte.
9. A concise **report** (markdown/YAML, KB-scale) under `reports/` capturing all of the above, with a
   forward note to 8c.

Claude reviews code + report with **two clean agents** — (i) the mandated **transliteration-fidelity
audit** (term-by-term: the driven ω-terms vs §4.7; the absorber profile vs the stated CAP provenance; the MMS forcing
vs the continuum operator) and (ii) an **adversarial overclaim / target-blindness / can't-fail-gate
review** (the reflection metric must have teeth vs the control) — plus an **independent arbiter re-run**
reproducing the 8b table + digest. We iterate with you until all hold. (MMS/arbiter/gate-A cannot catch a
faithful-but-wrong operator — the fidelity audit is the only check that does; it caught both 8a MAJORs.)

---

## 5. Repo hygiene (firewall is live)

- **Code** → `software/stage1_solver/src/stage1_solver/`, **tests** → `software/stage1_solver/tests/`
  (tracked).
- **All run output** → `software/stage1_solver/{runs,data,figures}/` — **gitignored**.
- **The one tracked artifact** → a small report under `software/stage1_solver/reports/`.
- Do **not** commit data; do **not** `git add` anything (the orchestrator stages explicit paths +
  commits after review); **never** `git add -A`; **never** write under `research/pde_audit/simulation/`.
  If you update `README.md`, keep it purely additive (a step-8b row + run command) and flag it for review.

---

## 6. Deliverable for review

Report back with: the driven-operator design (representation choice; the exact §4.7 terms added + how the
sign/conjugate structure was fixed; the ω→0 sanity-limit result); the drive-frequency set + why each is
target-blind and which is propagating; the absorber choice + parameters; the reflection-coefficient
method + what it measures + why non-circular + the target-blind floor + the achieved ratio + the
reflecting-control contrast (the teeth); the new MMS pieces + observed orders; the gate-A result (+ its
internal-consistency caveat); the well-posedness + convergence results; the pass/asserted check listing;
the new pytest count; harness exit 0 + the diagnostics digest reproduced; confirmation of
target-blindness + firewall + the asserted/counted separation + that the OPEN `S_η` source was NOT
invented; the report path; and any STOP-and-flag items. Then Claude reviews before anything is committed.
