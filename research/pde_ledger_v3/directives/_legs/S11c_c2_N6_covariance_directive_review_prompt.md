# Decision review — the S11c-c2 N6 SUFFICIENT-TEST (R_cov) build directive (physics-bearing pre-builder directive)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_covariance_directive.md` — orchestrator-written
directive for an astra-authored companion `scripts/S11c_c2_N6_covariance_sympy.py` that computes the
source-naturality residual `R_cov = ms − source_terms(μ_E.subs(Φ), V_E)` to close the OPEN per-engine N6 question.
Working dir `/var/projects/toy_physics`; paths under `research/pde_ledger_v3/`.

## Role + what this review is (and is NOT)
Pre-builder directive review of physics-bearing content: check the requested DECISIONS + physics. ⛔ No fictional-
script ablation (the deliverable does not exist yet), ⛔ no CAS run, ⛔ no build. Executable script-control tests are
the build legs' job. The test itself was already framed by an E1 vet sent to BOTH engines (they converged —
`_measurements/S11c_c2_N6_sufficient_test_vet_adjudication.md`); this review confirms the DIRECTIVE faithfully carries
it and is not leak-prone or under-specified.

## Method — DOCUMENT branch (read the sources first)
- The converged framing: `_measurements/S11c_c2_N6_sufficient_test_vet_adjudication.md` and the two vet reports
  `scratchpad/{codex,grok}_N6_suff_vet.log`.
- The open finding: `_measurements/S11c_c2_N6_reconcile_adjudication.md` (the over-clear correction + why localization
  ≠ sufficiency).
- The objects: `scripts/S11c_c2_N6_reconcile_sympy.py` (`source_terms` use, `μ_E`, `C_M`, `Φ`/frozen relations),
  `scripts/S11c_c2_N6_diagnostic_sympy.py` (`constitutive` :318-360 for `μ_E`/`μ_M`, `source_terms` :378-389); route-2
  spec §2 (`Φ`: `a_ρ`, `h_α`, prolongations) `_measurements/S11c_c2_N6_route2_spec_astra.md:39-83`.

## What to check (substantiate independently)
1. **Non-circularity — the crux.** Does the directive REQUIRE `ms_pred` be built from the IMPORTED Eulerian `μ_E`
   (`inputs.mu/ε`) with the SUPPLIED `Φ` substituted (`.subs`), and ⛔ forbid `b.material_pullback` / subtracting
   material from Eulerian / reusing `ms`/`ΔS`? Is the "vary-then-substitute-Φ vs pull-back-then-vary" commuting square
   the right non-circular test, or is there a hidden path by which `ms_pred` re-imports the builder under test (making
   it vacuous)?
2. **Is `Φ` correctly specified?** `θ↦θ+a_ρ`, `e_W↦e_W+h_α`, prolongations `D_i(·)`, with `a_ρ=u_i D_iρ₄/ρ₄`,
   `h_α=u_i D_iW_bg/W_bg` (LAB_HELD) / 0 (MATERIAL_ADVECTED), RHO4⇒`g_i=0`. Faithful to route-2 §2? Is using `Φ` as a
   `.subs` (the declared premise) legitimate, and is θ-independence of `Φ` (so vary/substitute commute) correctly the
   condition being tested rather than assumed away?
3. **Object correctness.** Is `R_cov = ms − ms_pred` per face/wave/retained-grade, emitted BEFORE Z/resolvent/weak
   extraction (carrier-class), the right object? Is the end-to-end guard `R_COV_INCREMENT = B(C_M, R_cov)` correctly
   the source channel propagated through the reconciled carrier? Do the velocities (`V_E≡V_M` SHA-equal per the
   reconcile) correctly isolate the μ channel?
4. **Able-to-fail.** Is the one-sided `Φ`-coefficient knife (actual path uses `2·a_ρ`/drops `h_α`, prediction keeps
   declared `Φ`) genuinely something the `a_ρ+h_α` truth table CANNOT see but `R_cov` MUST move? Is the θ-independent-
   junk knife right? Is the instrument required to be structured so both bite (separately parameterizable actual-Φ vs
   prediction-Φ)?
5. **Leakage / emit contract / tractability.** Any supplied expected `R_cov` value, any residual-zero exit/assert, any
   answer-bearing tag name, any A−A, any place emission is value-conditioned? Is the shared-sample in-process PIT +
   one-sided-certificate semantics correct? Is it genuinely carrier-class tractable (⛔ not an F/G-style full-symbolic
   wall)? Is the strict-vs-covariance interpretation correctly WITHHELD from the builder?
6. **Under-specification.** Anything astra must guess that changes what is computed — how `μ_E.subs(Φ)` threads the
   prolongations, how `ms_pred`'s `source_terms` call mirrors the real one, the `R_cov` keying, the knife structure.

## Physics filter
Report a finding only if it catches a way the sufficient test could be circular, vacuous, leak-prone, intractable, or
mislabel a covariance defect as covariant (or vice versa) — not "wrong on a different input."

## Output
Findings each with `file:line` + why it matters + minimal fix. End with **DIRECTIVE SOUND — CLEAR TO BUILD** or the
exact fold list. Brief, evidence-first.
