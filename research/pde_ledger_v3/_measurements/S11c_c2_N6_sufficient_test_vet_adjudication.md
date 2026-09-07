# S11c-c2 N6 sufficient-test framing vet — adjudication (2026-09-06)

Contested framing sent to BOTH engines (`directives/_legs/S11c_c2_N6_sufficient_test_vet.md`; reports
`scratchpad/{codex,grok}_N6_suff_vet.log`). Both EXIT=0. **They CONVERGED on the concrete test** (a strong result);
they differ only on a downstream naming/interpretation I defer.

## Converged (both engines)
- **A sufficient test IS needed.** Grok RETRACTS its earlier adjudication-review "SOUND": its (i)/(ii)/(iii)
  non-vacuity arguments miss a source-proportional coefficient error; the `a_ρ+h_α` truth table is **necessary-only**
  (cannot distinguish `1·Φ` from `2·Φ`). Carrier-reconciliation + source-localization is the complete statement only
  of the c2 **fold** ("Claim G"), ⛔ NOT of the constitutive covariance ("Claim I").
- **Option (a) quotient reduction: ill-posed / wrong identity — ⛔ DO NOT DO.** `material_pullback` already substitutes
  `θ_M=θ_E+a_ρ`, `e_M=e_E+h_α` (+prolongations) into the material operand, so `R_N6`'s normal form is `R_N6` itself;
  a lifted independent-field redesign either reproduces the current `R_N6` (does not vanish when sanctioned) or
  double-counts `Φ`. Not a missing Gröbner step; the wrong identity, and it risks the F/G wall.
- **Option (b) naive predicted-source: circular** (re-deriving `E(θ+a_ρ)` then varying IS `material_pullback`;
  `ΔS = b(Φ)−b` by construction; `V_E` SHA = `V_M` SHA so the velocities already agree).
- ⭐ **The sufficient test = the source-naturality commuting-square `R_cov`:**
  ```
  ms_pred = source_terms(μ_E.subs(Φ), V_E)   # Eulerian amplitude μ_E (= inputs.mu/ε, the one es uses) + supplied Φ,
                                              # vary-then-evaluate (NOT a second material_pullback)
  R_cov   = ms − ms_pred                      # actual material source (pull-back-then-vary) vs that prediction
  ```
  `Φ` = the already-emitted frozen field+jet maps (`θ↦θ+a_ρ`, `e↦e+h_α`, `D_i` of both) used as a **substitution**,
  not metadata. This tests `I_M ?= I_E∘Φ` on the ONLY remaining channel (carrier live-reconciled, velocities SHA-equal
  ⇒ the difference is μ). It is **non-circular** (shares the declared premise `Φ`, ⛔ not the builder under test) and
  **carrier-class tractable** (the 8-column `source_terms` circuit + the existing finite-field PIT; ⛔ not a symbolic
  slab). Optionally emit `B(C_M, R_cov)` as an end-to-end guard.
- **Able-to-fail discriminator:** one-sided `Φ`-coefficient corruption — `material_pullback` uses `2·a_ρ` (or drops
  `h_α`) while `ms_pred` uses the declared `Φ`: the truth table still lights the same 3/4 cases, but `R_cov` MUST move.
  Also θ-independent junk in `μ_M` ⇒ `R_cov` nonzero at `MATERIAL_ADVECTED.RHO4` (where `R_N6=0`). Codex's `κ_a,κ_h`
  knives are the same idea. The discriminator lives at retained order (`a_ρ,h_α` are `O(σ_W)`; the 18 live columns are
  exactly `(η,σ_W)∈{(0,1),(1,1)}`).

## Deferred (the naming/interpretation dispute — resolve AFTER R_cov)
- **Codex:** per-engine N6 = **strict** `NF(I_E−I_{M→E})=0`; `R_N6≠0` ⇒ not met; `R_cov` is a PARENT/source-naturality
  audit (validates `material_pullback`'s implementation of `Φ`), ⛔ not a nonzero offset to subtract from N6.
- **Grok:** per-engine N6 as specified = **covariance** `I_M = I_E∘Φ` (S11c_decisions N4 "must agree after that field
  redefinition"; directive "= 0 in the quotient"); `R_cov=0` closes it (⇒ `R_N6 = B(C_M,ΔS)` is the sanctioned Φ
  content).
- **My adjudication (G4):** build `R_cov` — it is decisive under BOTH readings. `R_cov ≠ 0` ⇒ a genuine defect
  (`material_pullback`/c2-binding does NOT implement the declared `Φ`) — a real finding both engines accept. `R_cov =
  0` ⇒ the material construction correctly implements `Φ`; THEN whether "covariant-via-Φ but not strictly invariant"
  counts as N6 passing is a clean **spec-interpretation decision** (strict R_N6=0 vs covariance I_M=I_E∘Φ), taken with
  the result in hand (and likely the user, since it defines what N6 asserts for the physics). ⛔ Do NOT pre-decide it.

⇒ NEXT: `R_cov` directive (tight — both engines specified it) → 2 decision legs → astra extends the reconcile
instrument → 2 build legs → adjudicate `R_cov` → resolve the strict-vs-covariance framing.
