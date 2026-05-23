---
unit_id: 064
batch: III.3
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-22T19:47:25-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 064

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:67-78` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:43-50`

**Issue:** The matched-layer assertions `C2_const - 1 == 0` and `Geq_const - g_phi^2 Npp/(KX Hw) == 0` are tautologies forced by the substitution dictionary `const_subs = {I1: Npp/Hw, I2: Npp/Hw**2}`. The substitution itself is exactly the result that needs proving; the script never instantiates a `chi_phi(y)` or `H(y)` to derive the integral reductions.

**Required change:**

In the sympy file, after line 66 ("print(\"G_eq =\", Geq)") and before the `banner("CONSTANT-COMPRESSIBILITY / MATCHED-LAYER LIMIT")` call on line 67, insert a concrete-profile reduction block:

```python
banner("MATCHED-LAYER INTEGRAL REDUCTION (concrete Gaussian profile)")

y, L = sp.symbols("y L", positive=True, real=True)
chi_phi_y = sp.exp(-y**2 / (2 * L**2))
H_const = Hw  # constant compressibility on the active layer

Npp_int = sp.integrate(chi_phi_y**2, (y, -sp.oo, sp.oo))
I1_int = sp.integrate(chi_phi_y**2 / H_const, (y, -sp.oo, sp.oo))
I2_int = sp.integrate(chi_phi_y**2 / H_const**2, (y, -sp.oo, sp.oo))

print("Npp_int =", Npp_int)
print("I1_int  =", I1_int)
print("I2_int  =", I2_int)
expect_zero("matched-layer I1 reduction", I1_int - Npp_int / Hw)
expect_zero("matched-layer I2 reduction", I2_int - Npp_int / Hw**2)
```

Keep the existing block at lines 67-78 unchanged; this new block precedes it and verifies that the `const_subs` relations actually follow from the integral definitions.

In the Mathematica file, after line 41 (`Print["G_eq = ", fmt[gEq]];`) and before line 43, insert the same construction:

```
banner["MATCHED-LAYER INTEGRAL REDUCTION (concrete Gaussian profile)"];

Clear[y, L];
$Assumptions = Element[{y, L, hw}, Reals] && L > 0 && hw > 0;

chiPhiY = Exp[-y^2/(2 L^2)];
nppInt = Integrate[chiPhiY^2, {y, -Infinity, Infinity}, Assumptions -> L > 0];
i1Int  = Integrate[chiPhiY^2/hw, {y, -Infinity, Infinity}, Assumptions -> L > 0 && hw > 0];
i2Int  = Integrate[chiPhiY^2/hw^2, {y, -Infinity, Infinity}, Assumptions -> L > 0 && hw > 0];

Print["Npp_int = ", fmt[nppInt]];
Print["I1_int  = ", fmt[i1Int]];
Print["I2_int  = ", fmt[i2Int]];
expectZero["matched-layer I1 reduction", i1Int - nppInt/hw];
expectZero["matched-layer I2 reduction", i2Int - nppInt/hw^2];
```

Re-establish the prior $Assumptions for gPhi, kX, npp, i1, i2, hw before the existing matched-layer block resumes at the next line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 064` and `redteam exec-mathematica 064` and confirm the new `matched-layer I1 reduction` and `matched-layer I2 reduction` lines appear in both outputs and both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- summary: Added concrete Gaussian matched-layer integral reductions before the existing constant-compressibility checks.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:112-118` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:79-84`

**Issue:** The "equilibrium softening" assertion `soft_eq - g_phi^2 * (Npp/Hw) == 0` is tautological because `Theta_eq` and `Lambda_eq` are both built from `const_subs` and the assertion compares the resulting expression against itself. A genuine check should derive `Lambda^2/Theta` symbolically from the closure (`Theta = Hw*Nss`, `Lambda = g_phi*Osp`) and only then apply the matched-layer closure `I2 = I1^2/Npp` to land on `g_phi^2 * I1`.

**Required change:**

In the sympy file, replace lines 112-118 (the block beginning `# Under the equilibrium-matched branch with constant H:` and ending with the second `expect_zero` call) with:

```python
# Closure-level definitions: Theta = H_w * N_(sigma sigma), Lambda = g_phi * O_(sigma phi).
Theta_abs = sp.simplify(Hw * Nss)
Lambda_abs = sp.simplify(g_phi * Osp)
soft_abs = sp.simplify(Lambda_abs**2 / Theta_abs)
print("Lambda^2/Theta (closure form) =", soft_abs)
# The matched-layer closure forces I2 = I1^2 / N_pp; only THEN should the softening
# reduce to g_phi^2 * I1.
soft_matched = sp.simplify(soft_abs.subs(I2, I1**2 / Npp))
print("Lambda^2/Theta (matched layer) =", soft_matched)
expect_zero("equilibrium softening equals g_phi^2 I1", soft_matched - g_phi**2 * I1)
```

In the Mathematica file, replace lines 79-84 with:

```
thetaAbs = FullSimplify[hw*nss, Assumptions -> gPhi > 0 && kX > 0 && npp > 0 && hw > 0 && i1 > 0 && i2 > 0];
lambdaAbs = FullSimplify[gPhi*osp, Assumptions -> gPhi > 0 && i1 > 0];
softAbs = FullSimplify[lambdaAbs^2/thetaAbs, Assumptions -> gPhi > 0 && kX > 0 && npp > 0 && hw > 0 && i1 > 0 && i2 > 0];
Print["Lambda^2/Theta (closure form) = ", fmt[softAbs]];
softMatched = FullSimplify[softAbs /. i2 -> i1^2/npp, Assumptions -> gPhi > 0 && npp > 0 && hw > 0 && i1 > 0];
Print["Lambda^2/Theta (matched layer) = ", fmt[softMatched]];
expectZero["equilibrium softening equals g_phi^2 I1", softMatched - gPhi^2*i1];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 064` and confirm both `Lambda^2/Theta (closure form)` (which must print a non-trivial expression involving `I1, I2, N_pp, H_w`, NOT a final answer) and `Lambda^2/Theta (matched layer)` (which must print `g_phi^2 * I1`) appear, and the new assertion passes. The same applies to the Mathematica output.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- summary: Replaced the tautological equilibrium softening check with a closure-level derivation followed by matched-layer reduction.
- deviation: Added the matched-layer relation `N_pp = I1 H_w` after the requested `I2 = I1^2/N_pp` substitution, because the requested substitution alone leaves a nonzero `N_pp - H_w I1` residual.

## F3 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:30-65` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:26-41`

**Issue:** The docstring's main checks 1 ("`chi_sigma = g_phi chi_phi/H`") and 2 (overlap invariants from the integral definitions of `I1`, `I2`) are never executed. `I1`, `I2`, `Npp`, etc. are abstract positive symbols; `Osp = sp.simplify(g_phi*I1)` is a definition not a derivation; no local-linear-response argument is realised in code.

**Required change:**

In the sympy file, immediately after line 54 (`Hw = sp.symbols("H_w", positive=True, real=True)`) and before line 56 (`# Derived overlap invariants from chi_sigma = g_phi * chi_phi / H.`), insert:

```python
# --- Local static linear-response closure ---
# At fixed y, the local sigma-energy is (1/2) H(y) sigma^2 - g_phi chi_phi(y) sigma.
# Minimising over sigma yields the closure law sigma(y) = g_phi chi_phi(y) / H(y).
y_loc = sp.symbols("y_loc", real=True)
chi_phi_loc = sp.Function("chi_phi")(y_loc)
H_loc = sp.Function("H")(y_loc)
sigma_loc = sp.symbols("sigma_loc", real=True)
F_loc = sp.Rational(1, 2) * H_loc * sigma_loc**2 - g_phi * chi_phi_loc * sigma_loc
closure_solutions = sp.solve(sp.diff(F_loc, sigma_loc), sigma_loc)
assert len(closure_solutions) == 1, "linear-response minimiser must be unique"
chi_sigma_closure = closure_solutions[0]
print("closure chi_sigma =", chi_sigma_closure)
expect_zero("closure law chi_sigma = g_phi chi_phi/H", chi_sigma_closure - g_phi * chi_phi_loc / H_loc)

# --- Integral-level overlap invariants (concrete Gaussian profile) ---
y_int, L_int = sp.symbols("y_int L_int", positive=True, real=True)
chi_phi_g = sp.exp(-y_int**2 / (2 * L_int**2))
H_g = Hw
Npp_int_check = sp.integrate(chi_phi_g**2, (y_int, -sp.oo, sp.oo))
I1_int_check = sp.integrate(chi_phi_g**2 / H_g, (y_int, -sp.oo, sp.oo))
I2_int_check = sp.integrate(chi_phi_g**2 / H_g**2, (y_int, -sp.oo, sp.oo))
chi_sigma_g = g_phi * chi_phi_g / H_g
Osp_int_check = sp.integrate(chi_sigma_g * chi_phi_g, (y_int, -sp.oo, sp.oo))
Nss_int_check = sp.integrate(chi_sigma_g**2, (y_int, -sp.oo, sp.oo))
print("Npp_int =", Npp_int_check, "  I1_int =", I1_int_check, "  I2_int =", I2_int_check)
expect_zero("overlap O = g_phi * I1 (integral form)", Osp_int_check - g_phi * I1_int_check)
expect_zero("overlap N_ss = g_phi^2 * I2 (integral form)", Nss_int_check - g_phi**2 * I2_int_check)
```

(If F1's block is already inserted, this F3 block may share the `chi_phi`-`H_g`-`L` symbols; in that case keep the F3 block's checks but use distinct names like `Osp_int_check`, `Nss_int_check` so neither block overwrites the other's verification.)

In the Mathematica file, insert the corresponding block immediately after line 31 (the `$Assumptions = ... && hw > 0;` line) and before line 33:

```
(* --- Local static linear-response closure --- *)
Clear[yLoc, chiPhiLoc, hLoc, sigmaLoc];
Block[{},
  fLoc = (1/2) hLoc[yLoc] sigmaLoc^2 - gPhi chiPhiLoc[yLoc] sigmaLoc;
  closureSolutions = Solve[D[fLoc, sigmaLoc] == 0, sigmaLoc];
  If[Length[closureSolutions] =!= 1, fail["linear-response minimiser must be unique"]];
  chiSigmaClosure = sigmaLoc /. First[closureSolutions];
  Print["closure chi_sigma = ", fmt[chiSigmaClosure]];
  expectZero["closure law chi_sigma = g_phi chi_phi/H", chiSigmaClosure - gPhi*chiPhiLoc[yLoc]/hLoc[yLoc]];
];

(* --- Integral-level overlap invariants (concrete Gaussian profile) --- *)
Clear[yInt, lInt];
chiPhiG = Exp[-yInt^2/(2 lInt^2)];
hG = hw;
nppIntCheck = Integrate[chiPhiG^2, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0];
i1IntCheck  = Integrate[chiPhiG^2/hG, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
i2IntCheck  = Integrate[chiPhiG^2/hG^2, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
chiSigmaG = gPhi chiPhiG/hG;
ospIntCheck = Integrate[chiSigmaG*chiPhiG, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
nssIntCheck = Integrate[chiSigmaG^2, {yInt, -Infinity, Infinity}, Assumptions -> lInt > 0 && hw > 0];
Print["Npp_int = ", fmt[nppIntCheck], "  I1_int = ", fmt[i1IntCheck], "  I2_int = ", fmt[i2IntCheck]];
expectZero["overlap O = g_phi * I1 (integral form)", ospIntCheck - gPhi*i1IntCheck];
expectZero["overlap N_ss = g_phi^2 * I2 (integral form)", nssIntCheck - gPhi^2*i2IntCheck];
```

**Verification command:**
After Codex applies, the verifier will run both `redteam exec-sympy 064` and `redteam exec-mathematica 064` and confirm `closure chi_sigma`, `overlap O = g_phi * I1 (integral form)`, and `overlap N_ss = g_phi^2 * I2 (integral form)` lines all appear and all `expect_zero` calls pass.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- summary: Added local static linear-response and concrete Gaussian overlap-invariant checks before the abstract invariant definitions.
- deviation: none

## F4 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:43-50` and `:79-84`

**Issue:** The `.wl` mirrors the `.py` line-for-line, including the `constSubs` substitution dictionary and the `Solve[D[f, sigma] == 0, sigma]` choreography. This means the second engine reaches `PASS` via the same algebraic path as the first, providing no adversarial cross-check.

**Required change:**

After F1 lands (which adds concrete-profile integrals to the `.wl`), refactor the matched-layer block at lines 43-50 so that `c2Const` and `gEqConst` consume the integral-derived `nppInt, i1Int, i2Int` symbols rather than the abstract `i1, i2`:

Replace lines 43-50 with:

```
c2Const = FullSimplify[(gPhi*i1Int)^2/(nppInt*gPhi^2*i2Int), Assumptions -> lInt > 0 && hw > 0 && gPhi > 0 && npp > 0];
gEqConst = FullSimplify[gPhi^2*i1Int/kX, Assumptions -> lInt > 0 && hw > 0 && gPhi > 0 && kX > 0];

Print["C^2 | H=const = ", fmt[c2Const]];
Print["G_eq | H=const = ", fmt[gEqConst]];
expectZero["matched-layer coherence", c2Const - 1];
expectZero["matched-layer gain vs Stage-45 best-alignment formula", gEqConst - gPhi^2*nppInt/(kX*hw)];
```

(`nppInt`, `i1Int`, `i2Int` are the integrals from the F1 block.)

For the softening block (lines 79-84 after F2 is applied), the structure already differs from the `.py` if F2 was applied (using `softAbs`/`softMatched` rather than `constSubs`). No further refactor required there.

**Verification command:**
After Codex applies, the verifier confirms `redteam exec-mathematica 064` exits 0 and the `.wl` output's `C^2 | H=const`, `G_eq | H=const`, and `matched-layer coherence`/`matched-layer gain ...` lines reference the integral-derived quantities (not abstract `i1*i2` ratios). Visually compare to the `.py` flow: the `.wl` matched-layer block must use `nppInt, i1Int, i2Int`, while the `.py` matched-layer block continues to use `const_subs` over abstract symbols — the two engines should now reach the same numerical conclusion via different intermediate computations.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- summary: Refactored the Mathematica matched-layer check to use the concrete integral-derived `nppInt`, `i1Int`, and `i2Int` values.
- deviation: none
