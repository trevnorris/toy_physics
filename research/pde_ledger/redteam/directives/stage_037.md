---
unit_id: 037
batch: III.1
created_at: 2026-05-22T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-22T18:16:38Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 037

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl:95-188`

**Issue:** The Mathematica audit script is a section-by-section, variable-by-variable transliteration of the SymPy audit. In particular, Section 3 (Schur factorization) and Section 4 (continuum closed forms) mirror the SymPy choreography: the same intermediate variables (`deltaUW`, `xiTerm`, `alphaTerm`) are introduced in the same order with the same hand-supplied forms (`xiTerm = gU^2/aU`, `alphaTerm = gB^2/aPhi + (aU*gW + gR*gU)^2/(aU*deltaUW)`), and the same `*Expected` closed forms are hand-supplied as targets (e.g. `aExpected = (kU*kEtaEff - cEtaU^2)/(muEta*kU)`, `chiExpected = (kU*cEtaW + cUW*cEtaU)/(muU*Sqrt[muEta*muW])`). This means a common algebra error would be reproduced identically in both engines and the cross-engine PASS provides no independent corroboration. The Mathematica engine must derive the headline closed forms by a route that does not consist of hand-supplying the same path-2 expression that SymPy uses.

**Required change:**

1. In Section 3 of the `.wl` file (lines 95-123), replace the current Schur choreography with an independent reconstruction of `Xi` and `alpha`. Concretely:

   - Remove the hand-supplied closed forms for `xiTerm`, `alphaTerm`, and `sigmaExpected` (current lines 115-121).
   - Replace them with code that:
     a. Computes `sigmaWall = FullSimplify[cMat . LinearSolve[bMat, Transpose[cMat]], Assumptions -> $Assumptions]` (keep current line 114 as-is).
     b. Posits the ansatz `Sigma_wall = xi * I_2 + alpha * v . Transpose[v]` with `v = {{kappa0}, {kappa1}}` and solves for `xi` and `alpha` using two entries of `sigmaWall`. Use `Solve` on the `(1,1)` and `(1,2)` entries: from `sigmaWall[[1,2]] == alpha * kappa0 * kappa1`, recover `alphaSolved = FullSimplify[sigmaWall[[1,2]] / (kappa0 * kappa1), Assumptions -> $Assumptions]`; from `sigmaWall[[1,1]] == xi + alpha * kappa0^2`, recover `xiSolved = FullSimplify[sigmaWall[[1,1]] - alphaSolved * kappa0^2, Assumptions -> $Assumptions]`.
     c. Assert `expectZero["Sigma_wall (2,2) consistency with ansatz", sigmaWall[[2,2]] - xiSolved - alphaSolved * kappa1^2]`. This is the substantive cross-check: if the ansatz `xi I + alpha v v^T` is the correct factorization, the `(2,2)` entry must agree with the `xi`, `alpha` solved from the `(1,1)` and `(1,2)` entries.
     d. Print `xiSolved` and `alphaSolved` via `Print["Xi (recovered) = ", fmt[xiSolved]]` and `Print["alpha (recovered) = ", fmt[alphaSolved]]`. The verifier will compare these printed forms to the SymPy `Xi`, `alpha` printout for engine cross-check.

2. In Section 4 of the `.wl` file (lines 125-177), derive the closed forms for `A` and `delta` by an alternative path that does not pre-define the `*Expected` literals:

   - Remove the hand-supplied `aExpected` and `deltaExpected` assignments (current lines 146, 161).
   - Replace each with a derivation:
     - `aDerived = FullSimplify[Together[k0 - gUCont^2/omegaU2], Assumptions -> $Assumptions]`. Then assert `expectZero["A reduces to closed form", a - aDerived]` (already true by definition of `a`; this is sanity — keep it as the first step). Then assert `expectZero["A numerator matches Schur form", FullSimplify[Numerator[Together[aDerived]] * (muEta*kU) - Denominator[Together[aDerived]] * (kU*kEtaEff - cEtaU^2), Assumptions -> $Assumptions]]`. This is the substantive identity: it confirms `a` equals `(kU*kEtaEff - cEtaU^2)/(muEta*kU)` *without* hand-supplying that form as `aExpected`; the closed form emerges as a consequence of the `Together` reduction.
     - Similarly for `delta`: `deltaDerived = FullSimplify[Together[deltaKAx/aDerived], Assumptions -> $Assumptions]`. Assert `expectZero["delta numerator matches closed form", FullSimplify[Numerator[Together[deltaDerived]] * (ell^2 * (kU*kEtaEff - cEtaU^2)) - Denominator[Together[deltaDerived]] * (kU*Pi^2*tw), Assumptions -> $Assumptions]]`.

3. Keep Sections 1, 2, 5 of the `.wl` file unchanged (lines 1-94 and 179-198 currently). Those sections do not contain transliterated algebraic choreography that violates the second-engine policy — Section 1 evaluates explicit integrals (independent verification suffices), Section 2 only prints definitions (no assertions), and Section 5 restates A and delta in their (now-independently-derived) closed-form variant.

4. Keep the existing assertions for `Chi`, `beta0`, `alphaMix`, `mMix` (current lines 165-168) since they are second-order derivatives of `A`, `delta0`, `Chi` and re-deriving them adds complexity without changing the transliteration verdict — once `A` and `delta` are derived independently, the rest follow.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 037` and confirm:
1. The new assertion `Sigma_wall (2,2) consistency with ansatz` appears in the saved output and prints `= 0`.
2. The new assertions `A numerator matches Schur form` and `delta numerator matches closed form` appear in the saved output and both print `= 0`.
3. The script exits 0.
4. The lines `aExpected = (kU*kEtaEff - cEtaU^2)/(muEta*kU)` and `deltaExpected = ...` no longer appear in the `.wl` file (grep returns empty).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl`
- summary: Replaced the Schur and continuum closed-form checks with recovered `Xi`/`alpha` and numerator identities derived from `Together` reductions.
- deviation: none
