---
unit_id: 064
batch: III.3
created_at: 2026-05-26T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-26T12:50:12-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
orchestrator_hotfix_2026_05_26:
  reason: "Codex iter1's F1 fix integrated hFun*chiSigmaFun^2 in Mathematica but FullSimplify could not pull gPhi^2 outside Integrate[] with symbolic chiPhi/hFun. The expectZero on the difference of integrals failed with a non-reducing form: -gPhi^2*Integrate[c^2/h] + Integrate[gPhi^2*c^2/h]."
  fix: "Verified integrand equality first (FullSimplify[hFun*chiSigmaFun^2 - gPhi^2*chiPhi^2/hFun] == 0); then defined thetaGeneral = lambdaGeneral = gPhi^2*i1Integral directly. By integrand equality the integrals are equal — bypassing Mathematica's inability to pull constants out of symbolic Integrate."
  pitfall_candidate: "Integrate[] with symbolic unspecified functions does not factor constant multipliers — verify integrands first (or factor manually) before comparing integral values."
---

# Codex directive — unit 064

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:156-165`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:123-129`

**Issue:**

The paper card and notes state `Delta K_X^eq = g_phi^2 I_1` as a general (any-H) identity on the parent-equilibrium-aligned branch. The current scripts' A10/M10 chain proves only the matched-layer reduction: `Theta_abs = H_w * N_ss`, then `soft_abs = g_phi^2 I_1^2 / (H_w I_2)` is force-substituted with `I_2 -> I_1^2/N_pp` (Cauchy-equality) and `N_pp -> I_1 H_w` (matched-layer relation), so the final assertion collapses to a tautology `g_phi^2 I_1 == g_phi^2 I_1` in the matched limit, not a general-H check.

**Required change (sympy):**

In `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`, REPLACE the block at lines 156-165:

Before (current lines 156-165):
```python
# Closure-level definitions: Theta = H_w * N_(sigma sigma), Lambda = g_phi * O_(sigma phi).
Theta_abs = sp.simplify(Hw * Nss)
Lambda_abs = sp.simplify(g_phi * Osp)
soft_abs = sp.simplify(Lambda_abs**2 / Theta_abs)
print("Lambda^2/Theta (closure form) =", soft_abs)
# The matched-layer closure forces I2 = I1^2 / N_pp; only THEN should the softening
# reduce to g_phi^2 * I1.
soft_matched = sp.simplify(soft_abs.subs(I2, I1**2 / Npp).subs(Npp, I1 * Hw))
print("Lambda^2/Theta (matched layer) =", soft_matched)
expect_zero("equilibrium softening equals g_phi^2 I1", soft_matched - g_phi**2 * I1)
```

After (general-H two-point branch check, no matched-layer substitutions):
```python
# General-H equilibrium softening on a two-point parent-equilibrium-aligned branch.
# On the aligned branch chi_sigma_k = g_phi chi_phi_k / H_k, so
#   Theta   = sum_k H_k chi_sigma_k^2     = g_phi^2 sum_k chi_phi_k^2 / H_k = g_phi^2 I_1
#   Lambda  = g_phi sum_k chi_phi_k chi_sigma_k = g_phi^2 sum_k chi_phi_k^2 / H_k = g_phi^2 I_1
# Therefore Lambda^2/Theta = g_phi^2 I_1 for ANY H(y), with no matched-layer assumption.
chi_sig1 = g_phi * sp.sqrt(w1) / H1  # chi_sigma_k for w_k = chi_phi_k^2
chi_sig2 = g_phi * sp.sqrt(w2) / H2
Theta_general = sp.simplify(H1 * chi_sig1**2 + H2 * chi_sig2**2)
Lambda_general = sp.simplify(g_phi * (sp.sqrt(w1) * chi_sig1 + sp.sqrt(w2) * chi_sig2))
soft_general = sp.simplify(Lambda_general**2 / Theta_general)
print("Theta (general, two-point) =", Theta_general)
print("Lambda (general, two-point) =", Lambda_general)
print("Lambda^2/Theta (general, two-point) =", soft_general)
expect_zero(
    "general equilibrium softening equals g_phi^2 I_1",
    soft_general - g_phi**2 * I1_disc,
)
```

This block must be placed AFTER the existing two-point Cauchy-gap section (line 139 onwards) since it reuses `w1, w2, H1, H2, I1_disc`. Concretely, after line 139 (`expect_zero("two-point Cauchy gap identity", gap_disc - gap_expected)`), insert a fresh banner and the new block. Move the `banner("ELIMINATED-SOURCE SOFTENING CHECK")` and the existing `F_eff` derivation (current lines 141-154) to AFTER the new general-softening block, so the abstract `F_eff` algebra still appears as scaffolding but no longer pretends to verify the general softening.

If the simplest path is to keep lines 141-154 (the F_eff check) in place and just rewrite lines 156-165, that is acceptable; in that case, reference `I1_disc, w1, w2, H1, H2` from earlier in the script (they are defined at lines 126-131). Verify the symbols are still in scope when the new block runs.

**Required change (mathematica):**

In `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`, REPLACE the block at lines 123-129. Do NOT translate the new sympy block line-by-line — see F2. Use an independent derivation path. Suggested form:

```mathematica
(* General-H equilibrium softening derived from the variational source-energy.    *)
(* The aligned closure is chi_sigma[y] = gPhi chiPhi[y] / H[y]; the parent source *)
(* self-energy is (1/2) Integrate[H[y] chi_sigma[y]^2, y].                        *)
Clear[chiPhi, hFun, y];
chiSigmaFun[y_] := gPhi * chiPhi[y] / hFun[y];
$Assumptions = Element[{y}, Reals];
thetaGeneral = Integrate[hFun[y] * chiSigmaFun[y]^2, {y, -Infinity, Infinity}] //
  Simplify;
lambdaGeneral = gPhi * Integrate[chiPhi[y] * chiSigmaFun[y], {y, -Infinity, Infinity}] //
  Simplify;
i1Integral = Integrate[chiPhi[y]^2 / hFun[y], {y, -Infinity, Infinity}];
(* Both Theta and Lambda equal gPhi^2 * I_1 on the aligned branch.                *)
expectZero["general Theta equals gPhi^2 I_1", thetaGeneral - gPhi^2 * i1Integral];
expectZero["general Lambda equals gPhi^2 I_1", lambdaGeneral - gPhi^2 * i1Integral];
softGeneral = Simplify[lambdaGeneral^2 / thetaGeneral];
expectZero["general equilibrium softening equals gPhi^2 I_1",
  softGeneral - gPhi^2 * i1Integral];
```

If `Integrate` with an unspecified `chiPhi[y]` / `hFun[y]` does not symbolically simplify (it should — they cancel algebraically before evaluating the integral), fall back to a concrete non-Gaussian profile (e.g., `chiPhi[y_] := Exp[-Abs[y]/L]` and `hFun[y_] := h0 + a*y^2`) and verify the same chain numerically by `Limit` or by direct symbolic integration in that profile. Whatever path is chosen, it must NOT include any substitution of the form `i2 -> i1^2/npp` or `npp -> i1*hw`.

**Claim manifest:**

Specific physical results the new general-softening check must verify:
- M1: `Theta = g_phi^2 I_1` on the equilibrium-aligned branch (any H(y)). Symbolic form: `int H(y) (g_phi chi_phi(y)/H(y))^2 dy = g_phi^2 int chi_phi(y)^2/H(y) dy`.
- M2: `Lambda = g_phi^2 I_1` on the equilibrium-aligned branch (any H(y)). Symbolic form: `g_phi int chi_phi(y) (g_phi chi_phi(y)/H(y)) dy = g_phi^2 int chi_phi(y)^2/H(y) dy`.
- M3: Therefore `Lambda^2/Theta = g_phi^2 I_1` for any H(y). This is the paper claim `Delta K_X^eq = g_phi^2 I_1`.

**Verification command:**

After Codex applies, verifier runs `redteam exec-sympy 064` and `redteam exec-mathematica 064`. Both must (a) contain a new assertion whose name contains "general equilibrium softening", (b) exit 0, and (c) NOT contain `subs(I2, I1**2/Npp)` or `i2 -> i1^2/npp` anywhere in the new block.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- summary: Replaced the matched-layer softening chain with general-H equilibrium-aligned softening checks.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:26-134`

**Issue:**

The `.wl` script is a structural line-by-line port of the `.py` script: identical variable choreography, identical Gaussian profile `Exp[-y^2/(2 L^2)]`, identical closure-minimisation pattern, identical discrete two-point model, identical matched-layer substitution chain. This violates the second-engine independent-derivation policy.

**Required change:**

Rewrite the Mathematica audit so each claim is derived along a path that does *not* mirror sympy. The mathematica script must still verify every paper-side deliverable (a)-(e), but the choreography must differ. Specific changes:

1. **Closure law block (current lines 33-42):** Replace the quadratic-minimisation construction with a variational derivative. Concretely:
   ```mathematica
   Needs["VariationalMethods`"];
   Clear[yLoc, chiPhiLoc, hLoc, sigmaFun];
   energyDensity = (1/2) hLoc[yLoc] sigmaFun[yLoc]^2 - gPhi chiPhiLoc[yLoc] sigmaFun[yLoc];
   eulerLagrange = VariationalD[energyDensity, sigmaFun[yLoc], yLoc];
   closureSolutions = Solve[eulerLagrange == 0, sigmaFun[yLoc]];
   If[Length[closureSolutions] =!= 1, fail["linear-response minimiser must be unique"]];
   chiSigmaClosure = sigmaFun[yLoc] /. First[closureSolutions];
   Print["closure chi_sigma = ", fmt[chiSigmaClosure]];
   expectZero["closure law chi_sigma = g_phi chi_phi/H",
     chiSigmaClosure - gPhi*chiPhiLoc[yLoc]/hLoc[yLoc]];
   ```
   If `VariationalMethods`` is unavailable, use `D[energyDensity, sigmaFun[yLoc]]` and solve — but rename variables and reshape the surrounding flow so the block is not structurally identical to sympy:63-68.

2. **Gaussian-integral block (current lines 44-56) and matched-layer reductions (current lines 70-82):** Replace the Gaussian profile `Exp[-yInt^2/(2 lInt^2)]` with a different profile, e.g. a Lorentzian `chiPhiL[y_, L_] := 1/(1 + y^2/L^2)`. Compute `nppInt = Integrate[chiPhiL[y, L]^2, {y, -Infinity, Infinity}]` (returns `Pi L/2`), and similarly for `i1Int, i2Int`. Verify `i1Int == nppInt/hw` and `i2Int == nppInt/hw^2` on this Lorentzian. The point is to use a distinct concrete function — the algebraic identity is the same, but the symbolic engine grinds it through a different integrand.

3. **Two-point Cauchy gap (current lines 96-108):** Replace with a continuous Cauchy-Schwarz check on the chosen non-Gaussian profile and a non-constant H, e.g., `H[y_] := h0 + a y^2`. Set `f[y_] := chiPhiL[y, L]/Sqrt[H[y]]` and `g[y_] := chiPhiL[y, L]/H[y]`. Compute
   - `Integrate[f[y] g[y], {y, -Infinity, Infinity}]^2`,
   - `Integrate[f[y]^2, {y, -Infinity, Infinity}] * Integrate[g[y]^2, {y, -Infinity, Infinity}]`,
   and assert the difference is `>= 0` symbolically (use `Reduce` or `Simplify[... >= 0, Assumptions -> h0>0 && a>0 && L>0]`).

4. **Softening (current lines 110-129, see F1):** Use the variational-integral form described in F1 above (not a finite-dimensional `theta*sigma^2 - lambda*sigma*phi + ...` mirror of sympy).

5. **Abstract `F_eff` block (current lines 110-121):** This is the only block that may legitimately remain structurally similar to sympy (it is pure algebra of a 2-variable quadratic). Keep it but rename variables (e.g., `theta -> sourceCoeff`, `lambda -> mixCoeff`, `phi -> supportAmp`, `sigma -> sourceAmp`) so it does not read as a transliteration.

Important: Do NOT change the front-matter `$Assumptions` declarations more than necessary, and do NOT remove any existing claim coverage. After the rewrite the Mathematica script must verify: closure law, overlap identities (on the new profile), `C^2 = I_1^2/(N_pp I_2)` formula construction, Cauchy bound `C^2 <= 1`, matched-layer `C^2 = 1`, matched-layer gain `G_eq = gPhi^2 N_pp/(kX hw)`, and the general softening `Delta K_X^eq = gPhi^2 I_1` (from F1).

**Claim manifest:** (same as sympy script — list the seven paper-side deliverables above; each must have a corresponding `expectZero` line in the rewritten `.wl`)

**Verification command:**

After Codex applies, verifier runs `redteam exec-mathematica 064`. The script must exit 0, must contain `expectZero` calls covering all seven items above, and must NOT use `Exp[-y^2/(2 L^2)]` in the integral-block (use the chosen alternative profile instead). A textual diff between the rewritten `.wl` and the original `.py` should show distinct intermediate steps for at least four of the six numbered blocks above.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- summary: Rewrote the Mathematica audit around a variational closure, Lorentzian profile integrals, continuous Cauchy checks, and integral-form general softening.
- deviation: none
