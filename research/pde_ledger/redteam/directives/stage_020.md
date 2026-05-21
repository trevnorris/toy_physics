---
unit_id: 020
batch: I.2
created_at: 2026-05-21T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-21T19:39:47Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 020

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl` (create new file)

**Issue:** Unit 020 has no Mathematica script, but the unit's manifest entry has `is_status_only_candidate: False` and the unit is not a checkpoint, so both engines are required by the second-engine policy. The SymPy script alone proves nothing against a second engine; the determinant 1/27, the closed-form wall slopes, and the compensated-branch identities are currently single-engine.

**Required change:**

Create the file at the Target path. The script must use Mathematica's own symbolic engine (`Solve`, `Det`, `FullSimplify`, `ThreeJSymbol`/`ClebschGordan` as needed for Gaunt — do not import any value from SymPy or hard-code SymPy output). It must independently derive each claim from the symbolic definitions of D0, D01, D21, D41, K1, H_even, Xi1, then assert each via `If[!TrueQ[FullSimplify[lhs - rhs] === 0], Print["FAIL: <label>"]; Exit[1]]`. Final `Print["STATUS: PASS"]; Exit[0]`.

Structurally the script MUST NOT mimic the SymPy file's function/variable choreography (no helper named `assert_zero`, no Python-style ordering of symbol declarations). Use Mathematica idioms: top-level `Module`, `SetAttributes`, `ClearAll`. The independence test is whether the algebra reads as a fresh Mathematica derivation, not a port.

**Claim manifest:**

M1. Real-Y20 self-overlap ratios via Gaunt coefficients. Let `g[m1, m2, m3] = ThreeJSymbol[{2,0},{2,m2},{2,m3}]` (or equivalent `ClebschGordan` / explicit closed form for the Gaunt integral). Verify:
  - `lambda_0 := gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0) == 1` (compute uniformly — do not short-circuit on m=0).
  - `lambda_1 := -1 * gaunt(2,2,2,0,1,-1) / gaunt(2,2,2,0,0,0) == 1/2`.
  - `lambda_2 := +1 * gaunt(2,2,2,0,2,-2) / gaunt(2,2,2,0,0,0) == -1`.
  - Same-sign cross terms vanish: `gaunt(2,2,2,0,1,1) == 0` and `gaunt(2,2,2,0,2,2) == 0`.

M2. Determinant of the even-gate Jacobian.
  Define D0 = KSigma - B0 - Z0, D01 = dKSigma - B01 - Z01, D21 = -(dMSigma + B21 + Z21), D41 = -(B41 + Z41), K1 = D21 + D01/9, H_even = D41 - (2/3) D21 - D01/27. Then
  `Det[{{D[K1,dKSigma], D[K1,dMSigma]}, {D[H_even,dKSigma], D[H_even,dMSigma]}}] == 1/27`.

M3. Closed-form solve of the even gates.
  `Solve[{K1 == 0, H_even == 0}, {dKSigma, dMSigma}]` yields a unique solution with
  `dKSigma -> B01 + Z01 + 27 (B41 + Z41)` and
  `dMSigma -> -(B21 + Z21) + 3 (B41 + Z41)`.

M4. Compensated-branch deficit identities. Substituting the solution back:
  `D01 |_{sol} == 27 (B41 + Z41)`
  `D21 |_{sol} == -3 (B41 + Z41)`
  `D41 |_{sol} == -(B41 + Z41)`.

M5. Compensated normalization defect.
  `Xi1 := N01/N0 - D01/D0` evaluated on the solved branch (with `D0 = KSigma - B0 - Z0`) satisfies
  `Xi1 |_{sol} == N01/N0 - 27 (B41 + Z41) / (KSigma - B0 - Z0)`.

The `.wl` script must independently encode and assert each of M1-M5.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 020` and confirm the new file exists, runs to completion, prints each of M1-M5's pass labels, prints `STATUS: PASS`, and exits 0. The verifier will additionally inspect the `.wl` source for transliteration patterns (presence of `assert_zero`-shaped helper, identical sequential ordering of symbol declarations, identical variable names with Python-style snake_case where Mathematica would use CamelCase) and treat any such pattern as a regression.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl`
- summary: Created the missing Mathematica audit script deriving and checking M1-M5 with Mathematica symbolic operations.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:20-27`

**Issue:** `real_y20_square_ratio(0)` returns the literal `sp.Integer(1)` before any Gaunt evaluation, so the subsequent `assert_zero('Y20 overlap lane 20', lam20 - 1)` at line 49 evaluates `1 - 1 == 0` regardless of what `gaunt(2,2,2,0,0,0)` actually equals. The m=0 lane is a tautology while m=1, m=2 are non-tautological.

**Required change:**

Replace lines 20-28 (the entire `real_y20_square_ratio` function) with a version that computes the ratio uniformly for all `m in {0, 1, 2}`, calling `gaunt` for both numerator and denominator. Keep the same-sign cross-term guard for `m != 0` only (the m=0 "same-sign" case IS the base itself and must not be flagged as a cross term).

Before:

```python
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
    if same_sign != 0:
        raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```

After:

```python
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m != 0:
        same_sign = sp.simplify(gaunt(2, 2, 2, 0, m, m))
        if same_sign != 0:
            raise AssertionError(f"Real-harmonic same-sign cross term should vanish for m={m}: {same_sign}")
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```

Do not change any other line of the file. The function signature, name, and call sites at lines 46-48 must remain identical.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 020` and confirm:
  - the script exits 0 with `STATUS: PASS`,
  - the source no longer contains `if m == 0: return sp.Integer(1)`,
  - the source contains the new `if m != 0:` guard wrapping the same-sign check,
  - the printed output is unchanged (all three lane assertions still pass).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py`
- summary: Replaced the m=0 shortcut with a uniform Gaunt numerator-over-denominator ratio and guarded same-sign checks only for nonzero m.
- deviation: none
