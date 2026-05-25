---
unit_id: 076
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
verification_status: pending
applied_at: 2026-05-22T23:21:39-06:00
findings_applied: 3
findings_blocked: 0
---

# Codex directive — unit 076

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:30-63`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl:26-61`

**Issue:** All four `expect_zero` calls verify identities forced by definitions made one to three lines earlier (see report F1 for the substitution chains). To make the assertions substantive, anchor the load-bearing pieces (`U`, `Theta_w` definition, the healing-length substitution) to an independent derivation that does not coincide with the assertion target.

**Required change (SymPy `.py`):**

Edit the block starting at line 30 (after the `banner(...)` call) through line 63. The new structure should be:

1. **After line 32 (symbol declaration)**, add a symbolic polytropic-index parameter `n_poly` so the EOS is index-parameterized:

```python
n_poly = sp.symbols("n_poly", positive=True, real=True)
```

2. **Replace line 35 (`P = K * rho**5`)** with an index-parameterized form:

```python
P = K * rho**(1 + 1/n_poly)
```

3. **Replace line 36 (`U = K * rho**5 / 4`)** with the polytrope enthalpy-density derivation. For a polytrope `P = K rho^(1 + 1/n)`, the per-mass enthalpy is `h_mass = ((n+1)/n) K rho^(1/n)`, so the per-particle enthalpy is `h = m * h_mass = m (n+1) K rho^(1/n) / n`. Encode this and then check that for `n_poly = 5` it equals `5 K rho^4 / 4`:

```python
# Polytrope enthalpy per particle: h = m * (n+1)/n * K * rho^(1/n)
h_polytrope = sp.simplify(m * (n_poly + 1) / n_poly * K * rho**(sp.Integer(1)/n_poly))
# Sanity: at n_poly = 5, h_polytrope = 6 m K rho^(1/5) / 5  (per-mass enthalpy form)
# But the script's claimed identity uses the energy-density convention:
#   U(rho) = integral of P/rho^2 dV over rho  (per-mass)
#   times rho to get per-volume; then differentiate in rho to recover h(rho).
# Derive U from the integral and then h = dU/drho:
U_per_mass = sp.simplify(sp.integrate(P / rho**2, rho))
U = sp.simplify(rho * U_per_mass)
```

4. **Replace line 38 (`h = sp.diff(U, rho)`)** with the derivation kept, but specialize to `n_poly = 5` BEFORE the assertion so that the printed `cs2` and `h` retain their current numeric form:

```python
cs2 = sp.diff(P, rho) / m
h = sp.diff(U, rho)
# Specialize to the n=5 case for the printed numbers and the assertion
cs2_n5 = sp.simplify(cs2.subs(n_poly, 5))
h_n5 = sp.simplify(h.subs(n_poly, 5))
print("c_s^2(rho)  [n=5] =", cs2_n5)
print("h(rho)      [n=5] =", h_n5)
expect_zero("n=5 enthalpy identity", h_n5 - m * cs2_n5 / 4)
# Non-tautology check: the same identity must FAIL for n=3 (gives factor 1/3, not 1/4)
residual_n3 = sp.simplify(h.subs(n_poly, 3) - m * cs2.subs(n_poly, 3) / 4)
print("n=3 residual (should be NONZERO) =", residual_n3)
if residual_n3 == 0:
    raise AssertionError("n=3 residual is zero — identity does not actually distinguish n=5")
```

Delete the original line 40 (`print("c_s^2(rho) =", sp.simplify(cs2))`) and line 41 (`print("h(rho)     =", sp.simplify(h))`); they are replaced by the n=5-specialized prints above.

5. **Lines 44-48 (Theta_w block)**: keep `mu_star = lambda_mu * m * csw**2 / 4` and `Theta_w = 4 * rho_w^2 * mu_star^2 / (hbar^2 * csw^2)` definitions, but replace the hand-typed `Theta_expected` on line 46 with an independent route: solve `mu_star - lambda_mu * m * csw**2 / 4 = 0` for `mu_star`, substitute, and check the result. Concretely, replace lines 44-48 with:

```python
# Define mu_star symbolically (independent symbol), then impose the enthalpy-lock condition
mu_star_sym = sp.symbols("mu_star_sym", positive=True, real=True)
enthalpy_lock = mu_star_sym - lambda_mu * m * csw**2 / 4  # set this to zero
mu_star_solved = sp.solve(enthalpy_lock, mu_star_sym)[0]
Theta_w = sp.simplify(4 * rho_w**2 * mu_star_solved**2 / (hbar**2 * csw**2))
# Independent route: Theta_w as (2 rho_w mu_star / (hbar c_sw))^2
Theta_w_alt = sp.simplify((2 * rho_w * mu_star_solved / (hbar * csw))**2)
print("Theta_w (enthalpy lock) =", Theta_w)
expect_zero("Theta_w vs alternative-form derivation", Theta_w - Theta_w_alt)
```

This routes `Theta_w` through a `solve` call and cross-checks against an alternative algebraic form `(2 rho_w mu_star / (hbar c_sw))^2`, neither of which appears in the original block.

6. **Lines 50-56 (healing-lock block)**: replace the hand-typed `healing_sub = {ell: hbar / (2 * m * csw)}` with a `solve` call:

```python
healing_condition = csw - hbar / (2 * m * ell)  # the healing-length defining relation
ell_solved = sp.solve(healing_condition, ell)[0]
print("ell from healing condition =", ell_solved)
# Substitute c_sw in Theta_w using the inverse relation (solve for c_sw):
csw_from_ell = sp.solve(healing_condition, csw)[0]
Theta_w_in_ell = sp.simplify(Theta_w.subs(csw, csw_from_ell))
Theta_heal_target = sp.simplify(lambda_mu**2 * rho_w**2 / (16 * ell**2))
expect_zero("healing-lock reduction", Theta_w_in_ell - Theta_heal_target)
print("Theta_w (healing lock) =", Theta_heal_target)
```

Delete the old `healing_sub`, `Theta_heal`, `Theta_heal_expected`, `cs_sub` lines (50-56).

7. **Lines 58-63 (reference-branch block)**: keep the structure but make the `1/20` factor explicit by tying `25` to `(20)^2/16` algebraically so a future change to one updates the other:

```python
ref_factor = sp.Rational(1, 20)  # reference-branch convention: ell = a * ref_factor  (see F2 below for provenance)
ref_sub = {ell: a * ref_factor}
Theta_ref = sp.simplify(Theta_heal_target.subs(ref_sub))
print("Theta_w (reference branch, general a) =", Theta_ref)
Theta_ref_norm = sp.simplify(Theta_ref.subs(a, 1))
print("Theta_w (reference branch, normalized wall units) =", Theta_ref_norm)
# The "25" target is derived from ref_factor as (1/ref_factor)^2 / 16
ref_target = sp.simplify((1 / ref_factor)**2 / 16) * lambda_mu**2 * rho_w**2
expect_zero("normalized reference factor", Theta_ref_norm - ref_target)
```

This makes the `25` flow from `ref_factor` rather than being typed independently.

**Required change (Mathematica `.wl`):**

Mirror the SymPy changes, but use Mathematica's native machinery to ensure independent derivation:

1. **After line 28 (Clear[...])**, add `nPoly` to the cleared and declared symbol list, and to `$Assumptions` add `nPoly > 0`.

2. **Replace lines 33-36** with:

```mathematica
press = kConst*rho^(1 + 1/nPoly);
uPerMass = FullSimplify[Integrate[press/rho^2, rho], Assumptions -> $Assumptions];
uRho = FullSimplify[rho*uPerMass, Assumptions -> $Assumptions];
csSq = FullSimplify[D[press, rho]/mpsi, Assumptions -> $Assumptions];
hRho = FullSimplify[D[uRho, rho], Assumptions -> $Assumptions];
csSqN5 = FullSimplify[csSq /. nPoly -> 5, Assumptions -> $Assumptions];
hRhoN5 = FullSimplify[hRho /. nPoly -> 5, Assumptions -> $Assumptions];

Print["c_s^2(rho)  [n=5] = ", fmt[csSqN5]];
Print["h(rho)      [n=5] = ", fmt[hRhoN5]];
expectZero["n=5 enthalpy identity", hRhoN5 - mpsi*csSqN5/4];

residualN3 = FullSimplify[hRho /. nPoly -> 3, Assumptions -> $Assumptions] - mpsi*FullSimplify[csSq /. nPoly -> 3, Assumptions -> $Assumptions]/4;
residualN3 = FullSimplify[residualN3, Assumptions -> $Assumptions];
Print["n=3 residual (should be NONZERO) = ", fmt[residualN3]];
If[TrueQ[residualN3 === 0], fail["n=3 residual is zero — identity does not actually distinguish n=5"]];
```

3. **Replace lines 42-47** (mu_star and theta block) with a `Solve`-routed derivation, mirroring SymPy step 5. **Rename `thetaExpected` to `thetaCanonical`** to break the line-by-line correspondence:

```mathematica
Clear[muStarSym];
$Assumptions = $Assumptions && muStarSym > 0;
enthalpyLock = muStarSym - lambdaMu*mpsi*cSw^2/4;
muStarSolved = First[muStarSym /. Solve[enthalpyLock == 0, muStarSym]];
thetaW = FullSimplify[4*rhoW^2*muStarSolved^2/(hbar^2*cSw^2), Assumptions -> $Assumptions];
thetaCanonical = FullSimplify[(2*rhoW*muStarSolved/(hbar*cSw))^2, Assumptions -> $Assumptions];
Print["Theta_w (enthalpy lock) = ", fmt[thetaW]];
expectZero["Theta_w vs alternative-form derivation", thetaW - thetaCanonical];
```

4. **Replace lines 49-55** (healing-lock block) with a `Solve` route:

```mathematica
healingCondition = cSw - hbar/(2*mpsi*ell);
cSwFromEll = First[cSw /. Solve[healingCondition == 0, cSw]];
ellSolved = First[ell /. Solve[healingCondition == 0, ell]];
Print["ell from healing condition = ", fmt[ellSolved]];
thetaWInEll = FullSimplify[thetaW /. cSw -> cSwFromEll, Assumptions -> $Assumptions];
thetaHealTarget = FullSimplify[lambdaMu^2*rhoW^2/(16*ell^2), Assumptions -> $Assumptions];
expectZero["healing-lock reduction", thetaWInEll - thetaHealTarget];
Print["Theta_w (healing lock) = ", fmt[thetaHealTarget]];
```

5. **Replace lines 57-61** (reference-branch block) with the `ref_factor`-derived form:

```mathematica
refFactor = 1/20;  (* reference-branch convention: ell = a * refFactor  (see F2 below for provenance) *)
thetaRef = FullSimplify[thetaHealTarget /. ell -> a*refFactor, Assumptions -> $Assumptions];
thetaRefNorm = FullSimplify[thetaRef /. a -> 1, Assumptions -> $Assumptions];
Print["Theta_w (reference branch, general a) = ", fmt[thetaRef]];
Print["Theta_w (reference branch, normalized wall units) = ", fmt[thetaRefNorm]];
refTarget = FullSimplify[(1/refFactor)^2/16*lambdaMu^2*rhoW^2, Assumptions -> $Assumptions];
expectZero["normalized reference factor", thetaRefNorm - refTarget];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 076` and `redteam exec-mathematica 076` and confirm:
(a) both scripts exit 0;
(b) new printed lines `n=3 residual (should be NONZERO) = ...`, `ell from healing condition = ...` appear in both `.txt` outputs;
(c) the existing four `expect_zero`/`expectZero` PASS lines still appear with residual 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl`
- summary: Replaced the tautological enthalpy, Theta_w, healing-lock, and reference-factor checks with independent integral/Solve-derived checks and nonzero n=3 residuals.
- deviation: Used an exponent-indexed EOS `P = K*rho^n_poly`/`press = kConst*rho^nPoly` because the directive's literal `1 + 1/n` EOS contradicts the required n=5 identity and fails validation.

## F2 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:9` (docstring line 4) and the `ref_factor` introduction added by F1
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl` (the `refFactor` introduction added by F1)

**Issue:** Two items: (a) the docstring's line 4 says `Theta_w = (25/4) lambda_mu^2 rho_w^2` but the algebra in the script (and the paper throughout — `paper/appendices/stage_ledger.tex:221`, `paper/appendices/stage_appendix_part03.tex:130`, `paper/stages/stage_077.tex:24`, `paper/stages/stage_082.tex:46`) consistently says `25 lambda_mu^2 rho_w^2`. Human-resolved: the docstring is a stale typo; the math and the paper are correct. Fix the docstring to read `25`. (b) The `1/20` factor in `ell = a/20` is hardcoded with no provenance — add a comment so a future reader knows it is a convention rather than a derived number.

**Required change (SymPy `.py`):**

1. **Edit the docstring line 9**. Replace:

   Before:
   ```
       4. Reference-branch form Theta_w = (25/4) lambda_mu^2 rho_w^2 in normalized Family-1 wall units.
   ```

   After:
   ```
       4. Reference-branch form Theta_w = 25 lambda_mu^2 rho_w^2 in normalized Family-1 wall units.
   ```

   (Drop the parenthesized `/4`; the assertion and the downstream paper both use `25`.)

2. **Above the `ref_factor = sp.Rational(1, 20)` line** added by F1, insert a comment:

   ```python
   # Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20.
   # TODO(provenance): cite the upstream stage that fixes ref_factor. This factor is
   # the load-bearing piece of the "25" in the normalized reference identity.
   ```

**Required change (Mathematica `.wl`):**

Mirror the same comment addition above the `refFactor = 1/20` line added by F1:

```mathematica
(* Reference-branch convention: ell = a * refFactor with refFactor = 1/20.
   TODO(provenance): cite the upstream stage that fixes refFactor. This factor is
   the load-bearing piece of the "25" in the normalized reference identity. *)
```

**Verification command:**
After Codex applies, the verifier will run `git diff` on the two files and confirm:
(a) the docstring on line 9 of the `.py` reads `Theta_w = 25 lambda_mu^2 rho_w^2` (no `(25/4)`);
(b) both files have the `TODO(provenance)` comment block above the `ref_factor`/`refFactor` line;
(c) the saved `.txt` outputs still show the normalized reference factor as `25*lambda_mu**2*rho_w**2` (no change in computed value); both scripts still exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl`
- summary: Corrected the stale Python docstring reference factor and added provenance TODO comments above the reference-factor conventions in both scripts.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl:26-61`

**Issue:** The `.wl` script reproduces the SymPy script's variable choreography, definitions, and assertion ordering line-by-line. F1's required `.wl` edits already introduce `Integrate[]`, `Solve[]`, and the renamed intermediate `thetaCanonical`, which break the transliteration. F3 codifies one additional rename so the variable mapping is not 1-to-1 even after F1.

**Required change (Mathematica `.wl`):**

After F1's edits land, rename **`thetaHealTarget` to `thetaHealReduced`** throughout the script (it appears in the healing-lock block and the reference-branch block — about 4 occurrences). Do NOT rename anything in the SymPy script. This ensures the variable mapping `Theta_heal_expected (sympy) <-> thetaHealTarget (math)` is broken; the Mathematica intermediate now has a distinct name driven by its post-Solve provenance.

Concretely, in the F1-edited `.wl`:

- Healing-lock block (the `thetaHealTarget = FullSimplify[...]` line introduced by F1 item 4): rename to `thetaHealReduced`.
- The `expectZero["healing-lock reduction", thetaWInEll - thetaHealTarget]` line: rename the second operand to `thetaHealReduced`.
- The `Print["Theta_w (healing lock) = ", fmt[thetaHealTarget]]` line: rename to `thetaHealReduced`.
- The reference-branch block's `thetaRef = FullSimplify[thetaHealTarget /. ell -> a*refFactor, ...]` line: rename to `thetaHealReduced`.

**Verification command:**
After Codex applies, the verifier will run `grep -c thetaHealTarget` on the `.wl` (must return 0) and `grep -c thetaHealReduced` (must return >= 4). Both scripts must still exit 0 and the four `expectZero` residuals must still be 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl`
- summary: Renamed the Mathematica healing-lock target intermediate to `thetaHealReduced` throughout the healing and reference blocks.
- deviation: none
