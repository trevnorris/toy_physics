---
unit_id: 024
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-21T22:02:20Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 024

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:1-211`

**Issue:** The Mathematica audit is a line-by-line port of the SymPy audit. Both scripts share the same recursive `pairings` helper, the same hand-decomposed `i4`/`i6` index-pair sums, the same Maxwell/mixed coupled-pair shorthand symbols (`Delta, S, Q, H, P` in `.py` mirrored as `deltaPair, sPair, qPair, hPair, pPair` in `.wl`), the same Section IV.2 pre-substituted lane forms, and the same `lane_ratio` series construction. The second-engine policy requires the `.wl` to derive the angular algebra by a structurally independent route. The required change preserves all existing assertion content (every `expectZero` line stays — same name, same mathematical residual) and only restructures the setup layer.

**Required change:**

Edit the file in place. The following five concrete edits are required.

1. **Replace the angular-moment helpers with direct surface integrals over the 2-sphere.** Currently lines 31-55 define `pairings[list_List]`, `i4[i,j,k,l]`, and `i6[i,j,k,l,m,n]` via a Wick-pair recursion. Replace this block with the following native-Mathematica formulation:

   ```mathematica
   n[1] = Sin[theta] Cos[phi];
   n[2] = Sin[theta] Sin[phi];
   n[3] = Cos[theta];

   i4[i_, j_, k_, l_] := Integrate[
     n[i] n[j] n[k] n[l] Sin[theta],
     {theta, 0, Pi}, {phi, 0, 2 Pi}
   ];

   i6[i_, j_, k_, l_, m_, nn_] := Integrate[
     n[i] n[j] n[k] n[l] n[m] n[nn] Sin[theta],
     {theta, 0, Pi}, {phi, 0, 2 Pi}
   ];
   ```

   Rename the last argument of `i6` from `n` to `nn` to avoid shadowing the new `n[_]` Cartesian-component function. Update the call site at line 65 accordingly: `i6[i, j, k, l, m, n]` becomes `i6[i, j, k, l, m, nn]` inside `tripleOverlap`, with the outer summation loop variable renamed from `n` to `nn` as well (the `{n, 1, 3}` becomes `{nn, 1, 3}` in the `Sum[]`). The `pairings` helper is no longer needed and must be removed.

2. **Remove the Maxwell/mixed shorthand symbols and expand them inline.** Currently lines 132-136 define `gU, gW, rPair, deltaPair, sPair, qPair, hPair, pPair`. Keep `gU = lambdaU*iEtaU; gW = lambdaW*iEtaW; rPair = lambdaR*iUW;` (these are physical group amplitudes, not algebraic shorthand) but delete the four lines defining `deltaPair, sPair, qPair, hPair, pPair`. In each subsequent target on lines 151-156, expand those symbols in place. The six target residuals must read mathematically:

   - `z0 - (gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2) / (omegaU^2*omegaW^2 - rPair^2)`
   - `z2 - ((gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2)*(omegaU^2 + omegaW^2) - (gU^2 + gW^2)*(omegaU^2*omegaW^2 - rPair^2)) / (omegaU^2*omegaW^2 - rPair^2)^2`
   - `z4 - ((gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2)*((omegaU^2 + omegaW^2)^2 - (omegaU^2*omegaW^2 - rPair^2)) - (omegaU^2 + omegaW^2)*(gU^2 + gW^2)*(omegaU^2*omegaW^2 - rPair^2)) / (omegaU^2*omegaW^2 - rPair^2)^3`
   - `n0 - (omegaU^2*gW + rPair*gU)^2 / (omegaU^2*omegaW^2 - rPair^2)^2`
   - `n2 - 2*(omegaU^2*gW + rPair*gU)*((omegaU^2*gW + rPair*gU)*(omegaU^2 + omegaW^2) - (omegaU^2*omegaW^2 - rPair^2)*gW) / (omegaU^2*omegaW^2 - rPair^2)^3`
   - `n4 - ((omegaU^2*omegaW^2 - rPair^2)^2*gW^2 - 2*(omegaU^2*omegaW^2 - rPair^2)*(omegaU^2*gW + rPair*gU)^2 - 4*(omegaU^2*omegaW^2 - rPair^2)*(omegaU^2*gW + rPair*gU)*(omegaU^2 + omegaW^2)*gW + 3*(omegaU^2*gW + rPair*gU)^2*(omegaU^2 + omegaW^2)^2) / (omegaU^2*omegaW^2 - rPair^2)^4`

   The `zResp` and `nResp` numerator/denominator expressions at lines 138 and 142 must also be rewritten to use the inline forms: e.g.,

   ```mathematica
   zResp = Expand[Normal[Series[
     ((gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2) - (gU^2 + gW^2)*omega^2) /
     ((omegaU^2*omegaW^2 - rPair^2) - (omegaU^2 + omegaW^2)*omega^2 + omega^4),
     {omega, 0, 4}
   ]]];
   ```

   and the parallel rewrite for `nResp`. The `dResp` line 146 already uses `k - m*omega^2 - bResp - zResp` (no shorthand) and can stay as-is. Use `Together` to canonicalize each residual before the `FullSimplify` if needed for the equality test to land.

3. **Replace the Section IV.2 pre-substituted lane forms with a lambda-parameterized constructor.** Currently lines 171-173 read:

   ```mathematica
   x20ax = x0 + eps*x1;
   x21ax = x0 + eps*x1/2;
   x22ax = x0 - eps*x1;
   ```

   Replace with:

   ```mathematica
   xLane[lam_] := x0 + eps*lam*x1;
   x20ax = xLane[1];
   x21ax = xLane[1/2];
   x22ax = xLane[-1];
   ```

   The subsequent `xbarAx`, `axAx`, `bxAx` definitions and the four `expectZero` assertions on lines 174-180 remain unchanged.

4. **Preserve all existing `expectZero` assertion names and residual mathematical content.** Do not add, remove, rename, or reorder any of the 28 PASS lines visible in the saved output. The verifier matches the output line-by-line. The only changes are the helper definitions and the inline expansion of shorthand symbols; the asserted residuals themselves stay mathematically identical.

5. **Preserve the file's overall layout.** Keep the banner strings, the section ordering (I → II → III → IV → V), the `$Assumptions` blocks, the `Clear[...]` calls, the `FINAL STAGE-007 LEDGER` print block, and the `Exit[0]` line. Do not introduce new features or scope extensions.

**Claim manifest:** N/A (this is a transliteration finding, not a missing-script finding; the assertion list is preserved.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 024` and confirm:
- The script exits 0.
- The saved output contains all 28 PASS lines from the current run, in the same order, with the same names.
- Inspection shows: (a) `Integrate[..., {theta, 0, Pi}, {phi, 0, 2 Pi}]` calls present in `i4` and `i6` definitions, (b) no `deltaPair`, `qPair`, `hPair`, `pPair`, `sPair` bindings remain in the file, (c) the `pairings[list_List]` helper has been removed, (d) `xLane[lam_]` is defined once and called with the three λ values.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`
- summary: Replaced Wick-pair angular helpers with sphere integrals, expanded Maxwell/mixed shorthand inline, and introduced the lambda-parameterized lane constructor.
- deviation: none
