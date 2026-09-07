# S11c-c2 N6 SUFFICIENT TEST (R_cov, source-naturality commuting square) — build directive (a SCRIPT; astra)

**Role.** Close the OPEN per-engine N6 question the reconcile left: is the nonzero constitutive residual
`R_N6 = B(C_M, ΔS)` the **covariant transformation** of the source under the declared field map `Φ`, or a real
non-covariance / `material_pullback` defect? The reconcile proved the geometric carrier reconciles (`C_E=C_M`, live)
and localized the residual to the constitutive source — but that is localization only (tautological given `ΔC=0`), and
the `a_ρ+h_α` truth table is necessary-only (cannot tell `1·Φ` from `2·Φ`). This instrument computes the
**source-naturality residual** `R_cov` that IS decisive. ⛔ NOT pre-judged; the script PRINTS and decides nothing.

**Author = Codex `gpt-6-astra` (high).** ⛔ Orchestrator does not author (E1). This directive was framed by an E1
vet sent to BOTH engines (`_measurements/S11c_c2_N6_sufficient_test_vet_adjudication.md` — they CONVERGED on this
exact test) and gets its 2 decision legs before you build. Your job ends at **build → verify → report**. ⛔ No review
legs, ⛔ no `.claude/skills/`, ⛔ no scratch files beyond the deliverable + its `.out`.

## The object (implement THIS; derive from it, ⛔ not from any output)
Deliverable: a NEW companion `scripts/S11c_c2_N6_covariance_sympy.py` importing the reconcile/diagnostic machinery
(`S11c_c2_N6_reconcile_sympy` and/or `S11c_c2_N6_diagnostic_sympy`). Per fixed `(α,ρ)`:
```
ms      = source_terms(inputs, α, ρ, face, μ_M, V_M)     # ACTUAL material source: pull-back-then-vary
                                                          #   (μ_M = material_pullback(E) then vary; the diagnostic ms)
ms_pred = source_terms(inputs, α, ρ, face, μ_E.subs(Φ), V_E)   # PREDICTION: Eulerian amplitude μ_E with the SUPPLIED
                                                          #   field map Φ substituted in — vary-then-substitute-Φ
R_cov   = ms − ms_pred                                    # source-naturality residual, per face / wave / retained grade
```
- **`μ_E`** = the imported Eulerian constitutive amplitude `es` already uses (`inputs.mu/ε`, reconcile `:199,:207`) —
  ⛔ NOT a parallel reconstruction of μ.
- **`Φ` (the crux — prolong to the jet order PRESENT IN `μ_E`, ⛔ NOT the energy-level 0+1 dict).** `μ_E = EL(E) =
  ∂E/∂θ − Σ_i D_i(∂E/∂(∂_iθ))` is already varied, so `D_i` of the first-jet flux writes **second** jets
  (`theta_second`/`theta_didi`, `e_W_didi`) into `μ_E` (confirmed present in the imported composite,
  `S11c_b_exports.py:5629`; `DERIVATIVE_MAP` `S11c_b_brane_operator_sympy_audit.py:733`). The energy-level frozen dict
  (`frozen_relations` reconcile `:53-57` / `material_pullback` `b:1946-1968`) is fields + first jets ONLY — complete
  for the scalar density `E`, **incomplete for `μ`**. ⛔⛔ Do NOT reuse it as the amplitude substitution. Instead
  `μ_E.subs(Φ)` = the **differential prolongation of the declared field map through EVERY θ/`e_W` jet symbol that
  occurs in imported `μ_E`** (currently rank-2): `θ_I ↦ D_I(θ+a_ρ)`, `e_{W,I} ↦ D_I(e_W+h_α)` for each multi-index `I`
  present, with `a_ρ = u_i D_iρ₄/ρ₄`, `h_α = u_i D_iW_bg/W_bg` at LAB_HELD else 0 (RHO4 ⇒ `g_i=0⇒a_ρ=0`). Generate the
  prolonged map with the LIVE `b.total_derivative(..., background_depth=3)` + `DERIVATIVE_MAP` chain (⛔ not a hand-typed
  extension), apply with `simultaneous=True`, and emit a **pre-substitution domain-coverage census** proving every
  imported θ/`e_W` jet atom in `μ_E` has a map entry (an uncovered atom = a silent covariance mislabel). `Φ` is the N4
  map the reconcile predeclared as SUPPLIED/unfalsifiable; the prolongation adds no new premise, only completes it.
- **`V_E`, `V_M`** are the Eulerian and material face velocities; the reconcile established `V_E`≡`V_M` (SHA-equal), so
  `R_cov` isolates the μ channel — carry both to keep the source circuits honest.
- ⛔⛔ **NON-CIRCULARITY IS THE WHOLE POINT.** `ms_pred` must be built by substituting the SUPPLIED `Φ` into the
  IMPORTED Eulerian `μ_E` — ⛔ NOT by calling `b.material_pullback` (that is the builder under test; matching a second
  pullback is vacuous), ⛔ NOT by subtracting the material result from the Eulerian result, ⛔ NOT reusing `ms`/`ΔS`.
  The test is the commuting square **pull-back-then-vary (`ms`) vs vary-then-substitute-Φ (`ms_pred`)**; they commute
  iff `Φ` is θ-independent (it is: `a_ρ=u·g`) and the Jacobian/quadratic projection leaves no linear-wave residue.
  That commutation is the physics being tested.
- **Emit** (tag prefix `S11CC2_N6COV_`), each columns + `numerator_denominator` + `nonzero_modular_numerator` under one
  in-process joint PIT per case (⛔ never a cross-run join): `SOURCE_ACTUAL` (`ms`), `SOURCE_PREDICTED` (`ms_pred`),
  `R_COV` (`ms−ms_pred`) keyed face/wave/grade **before** `Z`/resolvent/weak extraction (carrier-class); and the
  end-to-end guard `R_COV_INCREMENT` = the reconcile **`closed_response`** on `(m_coeff, R_cov)` (signatures 6/9/12
  ONLY, same column keying as `SOURCE_CHANNEL`). ⛔⛔ NOT `build_increment(C_M, R_cov)` / `I(C_M, R_cov)` — that re-adds
  the bare `−C_M·p` signature-0 term and would emit a certified-nonzero increment when `R_cov=0` (the affine trap the
  reconcile already forbade). Also emit the frozen `Φ` used (the substitution map), a
  provenance census that `ms_pred` is built from imported `μ_E` + supplied `Φ` (⛔ not `material_pullback`), and the
  PIT provenance (primes/draws/δ, `family·max(per_prime)`).

## Supplied vs withheld
- **Supplied (unfalsifiable in this build):** the reconcile/diagnostic machinery; `μ_E`; the field map `Φ` (the N4
  premise); `source_terms`; `C_M`. Reuse by import — ⛔ do not re-derive Z/resolvent/μ.
- **Withheld:** ⛔ NO expected value for `R_cov` (or its increment). ⛔ No residual-zero exit/assert/loop. A certified
  nonzero numerator is a one-sided certificate of a covariance discrepancy; all-zero = "no nonzero found" at the
  conditional δ, ⛔ never "certified zero." The strict-vs-covariance interpretation of the result is the orchestrator's
  (⛔ do NOT label `R_cov=0` "N6 passes" or `R_cov≠0` "N6 fails" — PRINT the residual).

## Control — able-to-fail, one-sided, ⛔ never A−A (build-leg ablation; ship uncorrupted)
- **Φ-coefficient knife:** in a /tmp COPY, make the ACTUAL material path use `2·a_ρ` (or drop `h_α`) while `ms_pred`
  keeps the DECLARED `Φ` (`κ=1`) ⇒ `R_cov` MUST move, even though the `a_ρ+h_α` 4-case truth table is unchanged. Use
  `RHOBR_CONSTANT` for the `a_ρ` knife, `LAB_HELD` for the `h_α` knife. This is what the truth table cannot see.
- **θ-independent-junk knife (make it concretely executable):** ⚠ `wave_terms` (diagnostic `:363`) requires every
  nonzero term to carry exactly one recognized wave, so a literal θ-independent constant ERRORS or is projected away
  (vacuous control). Instead, an ACTUAL-only, wave-LINEAR perturbation inserted after `μ_M` and before `source_terms`:
  add `κ_j · J_μ · e_W` to `μ_M` with a fresh nonzero-sampled `J_μ`, `DIMENSION_SCHEMA[J_μ]=(-1,-2,1)`, baseline
  `κ_j=0`, corrupted `κ_j=1`; keep `μ_E∘Φ` (the prediction) UNCHANGED ⇒ `R_cov` MUST move (nonzero at
  `MATERIAL_ADVECTED.RHO4` where `R_N6=0` and the truth table is dark). Emit the source-level control delta.
- ⛔ The instrument must be structured so these corruptions bite (the actual-material Φ and the prediction Φ are
  SEPARATELY parameterizable); ⛔ do not bake a corruption into the shipped emit. If a knife's derivative is
  annihilated by retained projection at source level, that component is untestable at this order — emit the computed
  fact, ⛔ do not force it.

## ⭐⭐⭐ THE THREE SCRIPT CLAUSES + ⭐⭐ FOUR COROLLARIES (verbatim — same as the reconcile directive)
1. PRINT computed objects, ⛔ never state conclusions. 2. PRINT the residual, ⛔ never assert it (⛔ no residual-zero
exit). 3. Interpretation belongs to the record. Corollaries: (1) every object reached by COMPUTATION from the imported
`μ_E`+supplied `Φ` / material builders, ⛔ never hand-typed; (2) the tag NAME names the object, ⛔ never its value;
(3) `R_cov` operands come from INDEPENDENT routes (imported-Eulerian+Φ vs material-pullback), verified by one-sided
corruption — ⛔ never differencing an object against its own substitution; (4) emission conditional only on
package/quantity, ⛔ never on a payload's value.

## Deliverable + verification
- `scripts/S11c_c2_N6_covariance_sympy.py`; run per case with `timeout` guards, `.out` to an ABSOLUTE path OUTSIDE the
  repo, ONE case at a time (carrier-class, so light — ⛔ still never a full-symbolic zero-test).
- Verify: exists, non-empty, imports cleanly, `.out` shows `SOURCE_ACTUAL`/`SOURCE_PREDICTED`/`R_COV`/
  `R_COV_INCREMENT` with PIT tables + the frozen `Φ` census. Report token usage. ⛔ Do NOT interpret the `R_cov`
  disposition — that is the orchestrator's, after the build legs.
