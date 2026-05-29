---
unit_id: 083
batch: III.4
created_at: 2026-05-22T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-25T00:26:45-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 083

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:**
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:51-95`

**Issue:** The Mathematica script copies the SymPy script's closed forms for `Delta_0(F1)`, `Delta_inf(F1)`, `A_F1`, and `Omega(Pe)` verbatim. Both engines derive nothing independently. An independent re-derivation route on the Mathematica side is needed so a typo in either closed form would surface as engine disagreement.

**Required change:**

1. **Independent derivation of `Delta_0(F1)` and `Delta_inf(F1)` via the defining linear BVP.**
   Insert after the line `deltaInfF1 = FullSimplify[...]` (i.e., after current line 60) a new block:

   ```mathematica
   (* Independent BVP derivation: Delta_0 and Delta_inf are the two boundary
      values of the linear two-point BVP
        u''(s) - alphaF1^2 * u(s) == 0,  0 < s < 1,
        u(0) - (1/etaF1) * u'(0) == 0,   (Robin at the wall)
        u(1) == 1.                       (unit value at outer edge)
      Then Delta_0 = u(0) / (alphaF1^2 + etaF1 * something)  -- BUT the
      cleanest invariant is: Delta_0 and Delta_inf each satisfy a single
      algebraic equation that does NOT use the closed-form Cosh/Sinh
      shortcut.  Concretely we verify they each satisfy the defining
      identity below. *)

   (* Defining identity for Delta_0(F1):
        (alphaF1^2 * (alphaF1 * Sinh[alphaF1] + etaF1 * Cosh[alphaF1])) * Delta_0
          == etaF1 * (Cosh[alphaF1] - 1)
      Verify the closed-form delta0F1 above satisfies it. *)
   delta0Residual = FullSimplify[
     (alphaF1^2*(alphaF1*Sinh[alphaF1] + etaF1*Cosh[alphaF1])) * delta0F1
       - etaF1*(Cosh[alphaF1] - 1)
   ];
   expectZero["Delta_0(F1) defining-equation residual", delta0Residual];

   (* Defining identity for Delta_inf(F1):
        (alphaF1 * Sinh[alphaF1] + etaF1 * Cosh[alphaF1]) * Delta_inf
          == Cosh[alphaF1] + (etaF1/alphaF1)*Sinh[alphaF1] - 1
      Verify deltaInfF1 above satisfies it. *)
   deltaInfResidual = FullSimplify[
     (alphaF1*Sinh[alphaF1] + etaF1*Cosh[alphaF1]) * deltaInfF1
       - (Cosh[alphaF1] + (etaF1/alphaF1)*Sinh[alphaF1] - 1)
   ];
   expectZero["Delta_inf(F1) defining-equation residual", deltaInfResidual];
   ```

2. **Independent re-statement of `Omega(Pe)` via its rational decomposition.**
   Insert after the line `omega[pp_] := Pi*pp*(2*pp*Exp[pp] + Pi)/((4*pp^2 + Pi^2)*(Exp[pp] - 1));` (i.e., after current line 88) a check that the closed form satisfies the partial-fraction identity that originally generated it:

   ```mathematica
   (* Independent identity for Omega(Pe):
        Omega(Pe) * (4 Pe^2 + Pi^2) (Exp[Pe] - 1)
          == Pi * Pe * (2 Pe Exp[Pe] + Pi).
      This is a structural identity the closed form must satisfy.  It is
      mathematically the definition we typed in, but verifying it as an
      identity (rather than as the function body) catches paste errors of
      the form "exp(Pe) -> exp(2 Pe)" or "4 Pe^2 -> 4 Pe" that would
      otherwise pass silently. *)
   With[{pp = pTest},
     omegaResidual = FullSimplify[
       omega[pp]*(4*pp^2 + Pi^2)*(Exp[pp] - 1)
         - Pi*pp*(2*pp*Exp[pp] + Pi),
       Assumptions -> pTest > 0
     ];
   ];
   expectZero["Omega(Pe) identity residual", omegaResidual];
   ```

   Before this block, declare `pTest` as a symbolic positive: insert immediately after the `omega[pp_] := ...` line:
   ```mathematica
   ClearAll[pTest];
   ```

   (Note: this check is structurally equivalent to verifying the function body, but it forces FullSimplify to combine the closed form back into the product side — a transcription error in the body would not vanish under FullSimplify.)

3. **Independent computation of `A_F1` via `y*tan(y) = eta` substitution.**
   Insert after the line `aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yRoot^2), 40];` (i.e., after current line 63):

   ```mathematica
   (* Independent check on A_F1: y*Tan[y] == eta must hold at yRoot, and
      A_F1 must be derivable directly from kappa and that eigenvalue.
      Re-derive A_F1 via a separate symbolic substitution and verify it
      matches the closed form numerically. *)
   yRootResidual = N[yRoot*Tan[yRoot] - etaF1, 40];
   expectApprox["y_F1 satisfies y Tan[y] = eta", yRootResidual, 0, 10^-20];
   aF1Indep = N[(kappaF1 + (Pi/2)^2)/(kappaF1 + yRoot^2), 40];
   expectApprox["A_F1 independent vs closed-form", aF1 - aF1Indep, 0, 10^-30];
   ```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 083` (single-seat Mathematica) and confirm the new `expectZero["Delta_0(F1) defining-equation residual", ...]`, `expectZero["Delta_inf(F1) defining-equation residual", ...]`, `expectZero["Omega(Pe) identity residual", ...]`, and `expectApprox["y_F1 satisfies y Tan[y] = eta", ...]` lines appear and PASS, and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- summary: Added independent Delta defining-equation residuals, y-root/A_F1 checks, and the Omega identity check.
- deviation: none

## F2 — hardcoded_result

**Target:**
- `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:75-78`
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:70-73`

**Issue:** `Theta_chi_coeff = 4.06863235008162`, `Theta_J_coeff = 0.927552032539308`, and the integer factor `136900` used in `Xi_chi = 136900 * Theta_chi_coeff` are typed-in literals with no derivation, no anchoring equation, and no upstream citation. Downstream windows depend on them.

**Required change:**

1. In `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`, immediately above line 75 (`Theta_chi_coeff = sp.Float("4.06863235008162")`), insert a comment block:

   ```python
   # SOURCE-ANCHOR (operator selectors):
   #   Theta_chi_coeff (= 4.06863235008162) and Theta_J_coeff (= 0.927552032539308)
   #   are the natural Family-1 operator selectors carried forward from earlier
   #   stages.  The integer prefactor 136900 = 370^2 = eta_F1 * (kappa_F1 + ...?)
   #   originates upstream as well.  If an upstream verifying script anchors
   #   these (e.g., a root of an operator equation), reference it here.  As of
   #   this audit no upstream sympy/mathematica script in this unit verifies
   #   them, so any transcription typo here will propagate silently.
   #   TODO(upstream-anchor): cite the producing stage's verifying script.
   ```

2. In `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`, immediately above line 70 (`thetaChi = N[4.06863235008162, 30];`), insert the analogous `(* SOURCE-ANCHOR ... *)` block (same content, Mathematica comment syntax `(* ... *)`).

3. Do NOT invent an equation for Theta_chi / Theta_J. The point of this finding is to make the unverified status explicit in-script, not to fabricate a derivation.

**Verification command:**
After Codex applies, the verifier will diff the two scripts and confirm both `SOURCE-ANCHOR` comment blocks are present immediately above the constant definitions. No new assertion is required for this finding; it is a documentation-of-trust-boundary fix.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- summary: Added SOURCE-ANCHOR comments above the unverified operator selector constants in both scripts.
- deviation: none

## F3 — insufficient_verification

**Target:**
- `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:64,110` (insertion points)
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:113` (insertion point)

**Issue:** The SymPy script has no assertion on items 2 and 3 of the docstring (wall constants, Pe windows, zeta windows, Pi/C_mix windows). The monotonicity claim in the theorem ledger (sympy 120, mathematica 124) is never exercised.

**Required change:**

1. **Defining-equation residuals in SymPy** — immediately after line 64 (the close of the `DeltaInf_F1 = sp.simplify(...)` block), insert:

   ```python
   # Defining-equation residuals for Delta_0 and Delta_inf (Family-1):
   #   (alpha^2 * (alpha sinh(alpha) + eta cosh(alpha))) * Delta_0 == eta (cosh(alpha) - 1)
   #   (alpha sinh(alpha) + eta cosh(alpha))           * Delta_inf == cosh(alpha) + (eta/alpha) sinh(alpha) - 1
   # These are non-tautological because they express Delta_* as the unique
   # solution of a linear equation, not as a chosen closed form.
   delta0_residual = sp.simplify(
       (alpha_F1**2 * (alpha_F1 * sp.sinh(alpha_F1) + eta_F1 * sp.cosh(alpha_F1))) * Delta0_F1
       - eta_F1 * (sp.cosh(alpha_F1) - 1)
   )
   expect_zero("Delta_0(F1) defining-equation residual", delta0_residual)

   deltaInf_residual = sp.simplify(
       (alpha_F1 * sp.sinh(alpha_F1) + eta_F1 * sp.cosh(alpha_F1)) * DeltaInf_F1
       - (sp.cosh(alpha_F1) + (eta_F1 / alpha_F1) * sp.sinh(alpha_F1) - 1)
   )
   expect_zero("Delta_inf(F1) defining-equation residual", deltaInf_residual)
   ```

2. **`y*tan(y) = eta` residual in SymPy** — immediately after line 67 (the `y_F1 = sp.nsolve(...)` line), insert:

   ```python
   y_residual = sp.N(y_F1 * sp.tan(y_F1) - eta_F1, 40)
   # nsolve solved this to tol=1e-30; verify residual is below 1e-25.
   if abs(y_residual) > sp.Float("1e-25"):
       raise AssertionError(f"y_F1 fails defining equation: residual = {y_residual}")
   print(f"y_F1 defining-equation residual = {y_residual}")
   ```

3. **Monotonicity check on `zeta_F1` in SymPy** — immediately after line 110 (the `zeta_max_F1 = ...` line), insert:

   ```python
   # Monotonicity check: zeta_F1(Pe) = A_F1 * Omega(Pe)^2 should be monotone
   # increasing in Pe over the relevant range, since Omega is positive and
   # increasing on (0, infinity).  The ledger claim "monotone in Xi on the
   # stable branch" follows because Pe = Xi * Delta is monotone in Xi on the
   # stable branch.  Check d zeta / d Pe > 0 at a sample of representative Pe
   # values spanning the chi and J windows.
   dzeta_dPe = sp.diff(zeta_F1, Pe_sym)
   for Pe_test_val in [sp.Float("10"), sp.Float("100"), sp.Float("1000"), sp.Float("10000")]:
       deriv_num = sp.N(dzeta_dPe.subs(Pe_sym, Pe_test_val), 30)
       print(f"  d zeta / d Pe at Pe = {Pe_test_val} : {deriv_num}")
       if deriv_num <= 0:
           raise AssertionError(
               f"zeta_F1 not monotone increasing at Pe = {Pe_test_val}: deriv = {deriv_num}"
           )
   ```

4. **Monotonicity check on `zeta_F1` in Mathematica** — immediately after line 113 (the `expectApprox["zeta_max^(F1) numeric check", ...]` line), insert:

   ```mathematica
   (* Monotonicity check: d zeta / d Pe must be positive at sample Pe values
      spanning the chi and J windows.  This exercises the "monotone in Xi on
      the stable branch" ledger claim (Pe = Xi*Delta is monotone in Xi). *)
   dzetaDpe[pp_] := D[zetaF1[ppSym], ppSym] /. ppSym -> pp;
   Module[{vals, signs},
     vals = N[dzetaDpe[#], 30] & /@ {10, 100, 1000, 10000};
     Print["  d zeta / d Pe samples = ", vals];
     signs = Sign[Re[#]] & /@ vals;
     If[!AllTrue[signs, # === 1 &],
       fail["zeta_F1 monotone increasing", vals],
       pass["zeta_F1 monotone increasing"]
     ];
   ];
   ```

   Note: the symbol `ppSym` is used as a fresh dummy here. Add `ClearAll[ppSym, dzetaDpe];` immediately above the new block to avoid contamination.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 083` and `redteam exec-mathematica 083` (sequential, single-seat Mathematica) and confirm:
- SymPy output now contains "Delta_0(F1) defining-equation residual = 0", "Delta_inf(F1) defining-equation residual = 0", "y_F1 defining-equation residual = ..." (printed), and four "d zeta / d Pe at Pe = ..." lines with positive values, plus exit 0.
- Mathematica output now contains "d zeta / d Pe samples = ..." and "PASS: zeta_F1 monotone increasing", plus exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- summary: Added SymPy Delta/y-root residual assertions and sampled zeta monotonicity checks in both SymPy and Mathematica.
- deviation: Raised the SymPy y-root nsolve precision so the required 1e-25 defining-equation residual check can pass.

## F4 — tautological_check

**Target:**
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:103-113`

**Issue:** The eleven `expectApprox` numeric self-checks compare the script's own computed value against a literal decimal that was generated from the same formula. They cannot fail unless someone edits the formula and forgets to also edit the literal — they verify file-edit-consistency, not physics.

**Required change:**

Do NOT delete the lines. Instead, after F1 lands, modify each numeric target so it references an **independent** value produced by F1's BVP-residual / identity-residual route. Concretely, for each `expectApprox` on lines 103-113, replace its numeric literal target with a freshly-evaluated symbolic expression that does not share the same algebraic body as the value being checked.

Mechanical edit, line by line:

- **Line 103** — currently:
  ```mathematica
  expectApprox["Delta_0(F1) numeric check", delta0F1, 1.73302079021525*10^-4, 10^-16];
  ```
  Replace literal target `1.73302079021525*10^-4` with the rearranged independent value:
  ```mathematica
  expectApprox["Delta_0(F1) numeric check",
    delta0F1,
    N[etaF1*(Cosh[alphaF1] - 1) /
      (alphaF1^2*(alphaF1*Sinh[alphaF1] + etaF1*Cosh[alphaF1])), 40],
    10^-16];
  ```
  (This is still derived from the same closed form, so the gain is small. Leave as-is if F1's `delta0Residual` `expectZero` lands; this `expectApprox` is then redundant but harmless. The acceptable alternative is to delete line 103. Codex picks ONE of: (a) delete, (b) keep as numeric sanity. Either is acceptable. Document the choice in the Applied block.)

- **Lines 104-113**: same pattern. Either delete each line, OR leave it; do not "fix" by changing the literal target to a different decimal — the issue is structural, not a digit mismatch.

If Codex is unsure, prefer keeping the lines as low-cost numerical sanity checks and let F1's `expectZero` lines carry the substantive verification.

**Verification command:**
Verifier will accept either: (a) lines 103-113 deleted, with F1's `expectZero` lines present and PASS, or (b) lines 103-113 unchanged but F1's `expectZero` lines added and PASS. The unit's overall information content is what matters; F4 is a clean-up directive that depends on F1 having landed.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- summary: Kept the existing numeric sanity checks unchanged while F1's new residual identities carry the substantive verification, per accepted option (b).
- deviation: none

## Self-test notes

- For F1, the new `expectZero` lines verify identities like `(alpha^2 (alpha sinh + eta cosh)) Delta_0 - eta (cosh - 1) == 0`. Substituting the closed form for `Delta_0` (which is exactly the ratio above) makes the residual zero by construction — but this is the correct check: it verifies that the closed form *is* the solution of the defining linear equation, which is the substantive physical claim. A transcription error in `delta0F1` (e.g., a missing factor of alpha) would make this residual nonzero.
- For F3, the monotonicity samples are at `Pe ∈ {10, 100, 1000, 10000}`, all positive. `Omega(Pe) = pi*Pe*(2 Pe e^Pe + pi)/((4 Pe^2 + pi^2)(e^Pe - 1))` is positive on `(0, infinity)` and approaches `pi/2` as `Pe -> infinity` from below — so `d Omega / d Pe > 0`, hence `d zeta / d Pe = 2 A Omega Omega' > 0`. The assertion is non-trivial (would catch a sign error in Omega) and non-tautological (the script does not currently know Omega is monotone).
- For F2, no executable test is added — only `SOURCE-ANCHOR` comment blocks. The verifier checks comment presence, not assertion residuals.
- Path specifications: all four findings edit files that already exist at the named paths. No new `.py` or `.wl` files are created.
