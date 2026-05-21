---
unit_id: 017
batch: I.2
created_at: 2026-05-21T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-21T13:25:18-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 017

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py:19-26`

**Issue:** The `real_y20_square_ratio` helper short-circuits the m=0 case by returning `sp.Integer(1)` without invoking SymPy's `gaunt`. The downstream assertion `assert_zero('Y20 overlap lane 20', lam20 - 1)` at line 44 therefore reduces to a literal `1 - 1 == 0` check that cannot fail. The m=0 ratio should be computed by the same machinery as m=1, 2.

**Required change:**

Replace the current function body at lines 19-26:

```
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
    if same_sign != 0:
        raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```

with:

```
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m != 0:
        same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
        if same_sign != 0:
            raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```

This removes the `if m == 0: return sp.Integer(1)` shortcut. For m=0 the function now returns `(-1)^0 * gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0)` which simplifies to `1` via SymPy, exercising the same Gaunt machinery as m=1, 2.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 017` and confirm:
- The script still exits 0.
- The output transcript still contains `lambda_(20) = 1, lambda_(21) = 1/2, lambda_(22) = -1.`.
- The function source no longer contains the `if m == 0: return sp.Integer(1)` lines.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py`
- summary: Removed the m=0 shortcut so the lane-20 ratio is computed through the same Gaunt overlap path as the other lanes.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py:79-80, 88`

**Issue:** Lines 79-80 assert `K1_wall - (-dMsym + dKsym/9) == 0` and `Hev_wall - (2/3 dMsym - dKsym/27) == 0`, but `K1_wall` and `Hev_wall` were constructed at lines 74-78 by literal `expand(...).subs(wall_only)` substitution into expressions whose only non-zero remaining symbols *are* `dMsym` and `dKsym`. The asserts thus merely re-verify the construction. The real claim the script makes (lines 140-141 of the source) is that the wall-only K1 / H_even formulas govern the three lane-specific gates `K1_gate_2m` / `Hev_2m`; that's what should be checked.

**Required change:**

Step 1: Delete lines 79 and 80 entirely:

```
    assert_zero("wall-only K1 specialization", K1_wall - (-dMsym + dKsym / 9))
    assert_zero("wall-only H_even specialization", Hev_wall - (sp.Rational(2, 3) * dMsym - dKsym / 27))
```

The `K1_wall` and `Hev_wall` symbolic forms are still used by `wall_matrix` (lines 91-94) and `sol_even` (lines 95-99), so do not delete the variables themselves — only the two asserts.

Step 2: After line 88 (which defines `Hev_22`), and before line 90 (`# Solve the wall-only even gates.`), insert the following six new assertions:

```
    # Cross-check: the generic lane formulas K1_wall = -M1 + K1w/9 and
    # Hev_wall = 2 M1 / 3 - K1w / 27 must reproduce the three explicit lane
    # gates K1_gate_2m and Hev_2m via X_A = eps * lambda_A * X1.
    generic_K1 = -M1 + K1w / sp.Integer(9)
    generic_Hev = sp.Rational(2, 3) * M1 - K1w / sp.Integer(27)
    assert_zero("generic K1 vs lane 20", K1_gate_20 - eps * lam20 * generic_K1)
    assert_zero("generic K1 vs lane 21", K1_gate_21 - eps * lam21 * generic_K1)
    assert_zero("generic K1 vs lane 22", K1_gate_22 - eps * lam22 * generic_K1)
    assert_zero("generic Hev vs lane 20", Hev_20 - eps * lam20 * generic_Hev)
    assert_zero("generic Hev vs lane 21", Hev_21 - eps * lam21 * generic_Hev)
    assert_zero("generic Hev vs lane 22", Hev_22 - eps * lam22 * generic_Hev)
```

Indentation must match the surrounding function body (4 spaces inside `main()`).

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 017` and confirm:
- The script still exits 0.
- The output transcript still contains the K1_(20,21,22) and H_even_(20,21,22) lines unchanged.
- The two old `wall-only K1/H_even specialization` `assert_zero` lines are gone.
- The six new `generic K1 vs lane *` / `generic Hev vs lane *` `assert_zero` lines appear in source.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py`
- summary: Replaced the wall-only specialization asserts with lane-gate cross-checks against the generic K1 and H_even formulas.
- deviation: none

## F3 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl` (new file)

**Issue:** No Mathematica partner script exists for stage 017. Per manifest, this unit requires both engines. The SymPy claims must be independently re-derived in Mathematica using `ThreeJSymbol[]` / `ClebschGordan[]` / direct integration of `SphericalHarmonicY[]` — NOT a line-by-line port of the SymPy code.

**Required change:**

Create the file at the path above. The script must use Mathematica's native conventions (`ThreeJSymbol[{l1,m1},{l2,m2},{l3,m3}]` or a direct integration of `SphericalHarmonicY[2,m1,th,ph] SphericalHarmonicY[2,m2,th,ph] SphericalHarmonicY[2,m3,th,ph]` over the sphere with `Sin[th]` measure) — not a transliteration of `sympy.physics.wigner.gaunt`.

Recommended structure (this is structural guidance, not literal code to be copied):

1. Implement an independent computation of `gauntMath[l1_,l2_,l3_,m1_,m2_,m3_]` either via `ThreeJSymbol` and the standard Gaunt-coefficient formula, or via direct integration of three `SphericalHarmonicY` factors over `{th, 0, Pi}, {ph, 0, 2 Pi}` with the `Sin[th]` measure. Pick the *direct integration* route to maximize independence from SymPy's symbolic Wigner code.
2. From `gauntMath`, derive `lambda[m_] := (-1)^m gauntMath[2,2,2,0,m,-m] / gauntMath[2,2,2,0,0,0]`. Compute `lambda[0]`, `lambda[1]`, `lambda[2]`.
3. Verify the same-sign cross term `gauntMath[2,2,2,0,m,m]` vanishes for m=1 and m=2 (selection rule m1+m2+m3=0).
4. Build the symbolic gate expressions independently: declare scalar symbols `dMsym`, `dKsym`, `M1`, `K1w`, `eps`, `D0`, `N0`. Define `K1wall = -dMsym + dKsym/9` and `Hevwall = 2 dMsym / 3 - dKsym/27` as the wall-only formulas. Compute the 2x2 matrix `{{D[K1wall,dKsym], D[K1wall,dMsym]}, {D[Hevwall,dKsym], D[Hevwall,dMsym]}}` and verify `Det[]` equals `1/27`.
5. Lane gates: define `dM[m_] := eps * lambda[m] * M1`, `dK[m_] := eps * lambda[m] * K1w`, `D01[m_] := dK[m]`, `D21[m_] := -dM[m]`, `D41[m_] := 0`. Then `K1gate[m_] := D21[m] + D01[m]/9` and `Hevgate[m_] := D41[m] - (2/3) D21[m] - D01[m]/27`. Verify `FullSimplify[K1gate[m] - eps * lambda[m] * (-M1 + K1w/9)] == 0` for each m in {0,1,2} and similarly for Hev.
6. Grouped trace/anomaly: define `Xbar[x20_,x21_,x22_] := (x20 + 2 x21 + 2 x22)/5`, `aX[x20_,x21_,x22_] := (2 x20 - x21 - x22)/10`, `bX[x20_,x21_,x22_] := (x21 - x22)/2`. Apply to {dM[0], dM[1], dM[2]}, {dK[0], dK[1], dK[2]}, {Xi[0], Xi[1], Xi[2]} := {-D01[0]/D0, -D01[1]/D0, -D01[2]/D0}, and {dP[0], dP[1], dP[2]} := {-N0 D01[0]/D0^2, ...}. Verify trace = 0 and b - 3 a = 0 in each case via `FullSimplify[...] == 0`.
7. For every check, if `FullSimplify[...]` is not zero, call `Print["FAIL: ", label]; Exit[1]`. Otherwise `Print["PASS: ", label]`. After all checks, `Print["STATUS: PASS"]` and exit 0.

Do NOT mimic the SymPy script's `grouped_trace_anomaly` Python function name or its variable order; use Mathematica-native naming. Do NOT call `sympy` or shell out. Do NOT copy the SymPy line structure.

**Claim manifest:**

The new Mathematica script must independently verify each of the following:

- M1: `gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0) = 1` (m=0 lane ratio, computed via `ThreeJSymbol` or direct `SphericalHarmonicY` integration, NOT returned as a literal `1`).
- M2: `(-1)^1 * gaunt(2,2,2,0,1,-1) / gaunt(2,2,2,0,0,0) = 1/2`.
- M3: `(-1)^2 * gaunt(2,2,2,0,2,-2) / gaunt(2,2,2,0,0,0) = -1`.
- M4: `gaunt(2,2,2,0,1,1) = 0` and `gaunt(2,2,2,0,2,2) = 0` (selection rule).
- M5: Grouped trace of {dM_20, dM_21, dM_22} = 0 and `b_M - 3 a_M = 0` where dM_2m = eps * lambda[m] * M1.
- M6: Grouped trace of {dK_20, dK_21, dK_22} = 0 and `b_K - 3 a_K = 0` where dK_2m = eps * lambda[m] * K1w.
- M7: For each m in {0,1,2}, `K1_gate_2m = D21_2m + D01_2m / 9` agrees with `eps * lambda[m] * (-M1 + K1w / 9)`.
- M8: For each m in {0,1,2}, `Hev_2m = D41_2m - (2/3) D21_2m - D01_2m / 27` agrees with `eps * lambda[m] * (2 M1 / 3 - K1w / 27)`.
- M9: The 2x2 Jacobian matrix `[[d K1_wall / d dKsym, d K1_wall / d dMsym], [d Hev_wall / d dKsym, d Hev_wall / d dMsym]]` has determinant `1/27`.
- M10: The homogeneous linear system `{K1_wall == 0, Hev_wall == 0}` has only the trivial solution `{dKsym -> 0, dMsym -> 0}`.
- M11: `Xi_load_2m = -D01_2m / D0` has grouped trace = 0 and `b_Xi - 3 a_Xi = 0`.
- M12: `dP0_2m = -N0 * D01_2m / D0^2` has grouped trace = 0 and `b_P - 3 a_P = 0`.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 017` and confirm:
- The new `.wl` file exists at the target path.
- The script exits 0.
- The transcript prints `PASS:` for each labeled check above, and ends with `STATUS: PASS`.
- The script does NOT contain the strings `sympy`, `gaunt(` (Python-style), or any direct port of the Python function names `real_y20_square_ratio`, `grouped_trace_anomaly`, `assert_zero`, `assert_nonzero`.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl`
- summary: Added a Mathematica audit that independently verifies the stage-017 overlap ratios, grouped trace/anomaly identities, gate formulas, wall Jacobian, and trivial wall-only solve.
- deviation: none
