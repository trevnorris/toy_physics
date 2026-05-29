---
unit_id: 154
batch: IV.6
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 154

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl:26-51`

**Issue:** The Mathematica audit script is a line-by-line transliteration of the SymPy script. It (a) replicates the SymPy variable choreography (same names, same order, same intermediate expressions), (b) types in the same hand-written closed-form `rShiftExpected = 1/4 - dg/Sqrt[1+r^2] + dg^2/(1+r^2)` rather than obtaining it via `Series[]`, (c) drops cross-terms with the same substitution dictionary `{dSigma0*dR -> 0, dSigma0*dS -> 0, dR*dS -> 0, dSigma0*dR*dS -> 0}` (line 47) that the .py uses verbatim (.py line 49–53), and (d) copies the SymPy script's incorrect banner text `"STAGE 137 — EXACT CO-EVOLVING CORE-MOUTH MAP"` (line 26). The .wl therefore provides no independent algebraic derivation — any algebra error in the .py would be silently mirrored.

**Required change:**

Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl` as follows (paths and line numbers refer to the current file state shown in the audit report):

1. Fix the banner string at line 26.

   Before:
   ```
   banner["STAGE 137 — EXACT CO-EVOLVING CORE-MOUTH MAP"];
   ```
   After:
   ```
   banner["STAGE 154 — EXACT CO-EVOLVING CORE-MOUTH MAP"];
   ```

2. Replace the shifted-R block (lines 37–39) so the expansion is derived by `Series` rather than typed in by hand.

   Before:
   ```
   rShift = Expand[rFun /. g -> gStar + dg];
   rShiftExpected = Expand[1/4 - dg/Sqrt[1 + r^2] + dg^2/(1 + r^2)];
   expectZero["exact shifted R formula", rShift - rShiftExpected];
   ```
   After:
   ```
   rShiftSeries = Normal[Series[rFun /. g -> gStar + dg, {dg, 0, 2}]];
   rShiftExpected = 1/4 - dg/Sqrt[1 + r^2] + dg^2/(1 + r^2);
   expectZero["exact shifted R formula", rShiftSeries - rShiftExpected];
   ```
   The `rShiftExpected` here is the paper-side claim from notes §3; the `rShiftSeries` is now obtained by Mathematica's own series expansion of `R(g_* + dg)` to second order, not by a substitution that mirrors the .py.

3. Replace the dPi linearization (lines 46–51) so the linear part is extracted by `Series` rather than by a hand-built substitution dictionary.

   Before:
   ```
   piExpr = (sigma0 + dSigma0) * (1 - (rStar + dR) * (sStar + dS));
   piLin = Expand[piExpr] /. {dSigma0*dR -> 0, dSigma0*dS -> 0, dR*dS -> 0, dSigma0*dR*dS -> 0};
   pi0 = sigma0*(1 - rStar*sStar);
   dPi = Expand[piLin - pi0];
   dPiExpected = Expand[(1 - rStar*sStar)*dSigma0 - sigma0*(rStar*dS + sStar*dR)];
   expectZero["dPi identity", dPi - dPiExpected];
   ```
   After:
   ```
   piExpr = (sigma0 + dSigma0) * (1 - (rStar + dR) * (sStar + dS));
   piLin = Normal[Series[piExpr, {dSigma0, 0, 1}, {dR, 0, 1}, {dS, 0, 1}]];
   pi0 = sigma0*(1 - rStar*sStar);
   dPi = Expand[piLin - pi0];
   dPiExpected = Expand[(1 - rStar*sStar)*dSigma0 - sigma0*(rStar*dS + sStar*dR)];
   expectZero["dPi identity", dPi - dPiExpected];
   ```
   `piLin` is now the genuine linear part of `piExpr` in `(dSigma0, dR, dS)` as computed by Mathematica's multivariate `Series` machinery, independent of any choice made in the .py. The `dPiExpected` line still mirrors the paper claim (notes §4) — keep it; it is the *target identity*, not an algebraic step copied from the .py.

Do not touch any other lines in the .wl. Do not touch the .py.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 154` and confirm:
- the script still Exit 0,
- the saved output's banner now reads `STAGE 154 — EXACT CO-EVOLVING CORE-MOUTH MAP`,
- both `expectZero` checks still PASS (residual zero),
- the `.wl` no longer contains the substitution rules `dSigma0*dR -> 0`, `dSigma0*dS -> 0`, `dR*dS -> 0`, `dSigma0*dR*dS -> 0`,
- the `.wl` contains at least two `Series[...]` calls (one for `rShift`, one for `piLin`).
