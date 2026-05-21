---
unit_id: 017
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 017 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

The script asserts that the grouped-lane ratios (lambda_(20), lambda_(21), lambda_(22)) = (1, 1/2, -1) arise from real spherical-harmonic Y20-triple overlap integrals (computed via SymPy `gaunt(2,2,2,0,m,-m)`); that under those ratios the wall-only inertia and stiffness anisotropies obey trace=0 and b=3a; that the wall-only contributions to the Stage-5 weak-axisymmetric gates K1 and H_even reduce to the formulas K1_wall = -delta_M + delta_K/9 and H_even,wall = 2 delta_M / 3 - delta_K / 27, whose linear system has determinant 1/27 (so only the trivial branch closes both gates); and that the induced Xi_load and prefactor shifts likewise satisfy trace=0 / b=3a. Mathematica is absent, so no cross-engine verification is performed.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 44 | `assert_zero('Y20 overlap lane 20', lam20 - 1)` | no (tautological — `real_y20_square_ratio(0)` returns `sp.Integer(1)` directly) |
| A2 | sympy | 45 | `assert_zero('Y20 overlap lane 21', lam21 - 1/2)` | yes (real Gaunt-ratio numerical check) |
| A3 | sympy | 46 | `assert_zero('Y20 overlap lane 22', lam22 + 1)` | yes (real Gaunt-ratio numerical check) |
| A4 | sympy | 58 | `assert_zero('wall inertia trace', Mbar)` | yes (depends on lam21, lam22 from gaunt) |
| A5 | sympy | 59 | `assert_zero('wall inertia b=3a', bM - 3*aM)` | yes |
| A6 | sympy | 60 | `assert_zero('wall stiffness trace', Kbar)` | yes |
| A7 | sympy | 61 | `assert_zero('wall stiffness b=3a', bK - 3*aK)` | yes |
| A8 | sympy | 79 | `assert_zero("wall-only K1 specialization", K1_wall - (-dMsym + dKsym/9))` | no (tautological — K1_wall constructed via direct substitution of wall_only=0 into a literal expression) |
| A9 | sympy | 80 | `assert_zero("wall-only H_even specialization", Hev_wall - (2/3 dMsym - dKsym/27))` | no (tautological by the same construction) |
| A10 | sympy | 100 | `assert_zero("wall-only even-gate determinant", wall_matrix.det() - 1/27)` | partial (the matrix entries are derivatives of tautologically-constructed expressions; det=1/27 is an arithmetic consequence) |
| A11 | sympy | 101 | `assert_nonzero("mutated wall-only determinant should fail", wall_matrix.det() + 1/27)` | yes (sign sanity check) |
| A12 | sympy | 114 | `if sol_even != [{dKsym:0, dMsym:0}]: raise` | yes (linear-system uniqueness given det != 0) |
| A13 | sympy | 116 | `assert_zero('Xi trace', Xibar)` | yes (depends on lane ratios) |
| A14 | sympy | 117 | `assert_zero('Xi b=3a', bXi - 3*aXi)` | yes |
| A15 | sympy | 118 | `assert_zero('prefactor trace', Pbar)` | yes |
| A16 | sympy | 119 | `assert_zero('prefactor b=3a', bP - 3*aP)` | yes |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py:19-26, 44`

**What's wrong:**

`real_y20_square_ratio` short-circuits the m=0 branch:

```
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    ...
    return sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)
```

It then immediately asserts `lam20 - 1 == 0` (line 44). Because `lam20` is hard-coded to `sp.Integer(1)` by the function (the m=0 branch never touches SymPy's `gaunt`), the assertion reduces to `1 - 1 == 0` and cannot fail no matter what the Gaunt-coefficient library produces. The lane-20 ratio is exactly what the script claims to verify; that verification needs to come from the same machinery (Gaunt symbols) used for m=1 and m=2.

**Why this matters:**

The script's grouped-lane signature (1, 1/2, -1) is the bedrock for everything downstream (Mbar, Kbar, Xibar, Pbar, the b=3a relations). If the m=0 normalization were ever to drift in SymPy's `gaunt` convention or get re-defined, the script would silently report PASS even though the other lane ratios would be off by the same factor. A real check should compute the m=0 ratio identically to the m!=0 branch.

**Required change:**

Drop the `if m == 0: return sp.Integer(1)` shortcut at lines 21-22 and let the function fall through to the same `(-1)^m * gaunt(2,2,2,0,m,-m) / gaunt(2,2,2,0,0,0)` formula it uses for m!=0. For m=0 this evaluates to `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)` and the `same_sign` check (line 24) is vacuous (gaunt(2,2,2,0,0,0) is the base itself, not a same-sign cross term). Wrap the `same_sign` check inside `if m != 0:` so it only fires for the m!=0 branches it was meant for.

**Verification:**

After the fix, the function body for m=0 runs `sp.simplify(gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0))` and returns `1`. The existing `assert_zero('Y20 overlap lane 20', lam20 - 1)` (line 44) still passes, but now exercises SymPy's `gaunt` rather than a literal `Integer(1)`. The output `lambda_(20) = 1` in the txt should still appear identically.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py:71-80`

**What's wrong:**

Lines 71-78 construct symbolic expressions and then immediately substitute to obtain the "wall-only" specializations:

```
D01_full = dKsym - B01 - Z01
D21_full = -(dMsym + B21 + Z21)
D41_full = -(B41 + Z41)
K1_full = sp.expand(D21_full + D01_full / 9)
Hev_full = sp.expand(D41_full - sp.Rational(2, 3) * D21_full - D01_full / 27)
wall_only = {B01: 0, B21: 0, B41: 0, Z01: 0, Z21: 0, Z41: 0}
K1_wall = sp.expand(K1_full.subs(wall_only))
Hev_wall = sp.expand(Hev_full.subs(wall_only))
```

Lines 79-80 then assert:

```
assert_zero("wall-only K1 specialization", K1_wall - (-dMsym + dKsym / 9))
assert_zero("wall-only H_even specialization", Hev_wall - (sp.Rational(2, 3) * dMsym - dKsym / 27))
```

But by direct construction, `K1_full.subs({B01:0, B21:0, Z01:0, Z21:0}) = -dMsym + dKsym/9` algebraically — sympy substitution and `expand` cannot give any other result. The reference RHS on the assert line was hand-written by reading the very LHS construction. The assertion is guaranteed by sympy semantics and cannot fail no matter what the physics is.

**Why this matters:**

The script's docstring/comments claim the wall-only gates reduce to "K1_wall = -delta_M + delta_K/9" and "H_even,wall = 2 delta_M / 3 - delta_K / 27" (lines 140-141). The current assertion does not test that claim — it only tests that sympy's `subs` and `expand` are working. A meaningful check would tie these formulas back to the lane-resolved quantities `K1_gate_2m` and `Hev_2m` (already computed at lines 82-88) by checking that they equal `eps * lambda_(2m) * (-M1 + K1w/9)` and `eps * lambda_(2m) * (2 M1 / 3 - K1w / 27)` respectively. That ties the "generic lane formula" claim to the actual three lanes via the verified gaunt ratios.

**Required change:**

Replace the two tautological asserts at lines 79-80 with cross-checks against the lane-resolved gates. Concretely, append (after line 88 where `K1_gate_*` and `Hev_*` are defined):

```
# Cross-check the generic lane formulas against the three explicit lanes.
generic_K1 = -M1 + K1w / sp.Integer(9)
generic_Hev = sp.Rational(2, 3) * M1 - K1w / sp.Integer(27)
assert_zero("generic K1 vs lane 20", K1_gate_20 - eps * lam20 * generic_K1)
assert_zero("generic K1 vs lane 21", K1_gate_21 - eps * lam21 * generic_K1)
assert_zero("generic K1 vs lane 22", K1_gate_22 - eps * lam22 * generic_K1)
assert_zero("generic Hev vs lane 20", Hev_20 - eps * lam20 * generic_Hev)
assert_zero("generic Hev vs lane 21", Hev_21 - eps * lam21 * generic_Hev)
assert_zero("generic Hev vs lane 22", Hev_22 - eps * lam22 * generic_Hev)
```

Then delete the two tautological assertions at lines 79 and 80. The `K1_wall`, `Hev_wall`, `wall_matrix`, `sol_even`, and the lines 100-101 determinant checks may all remain — those are downstream uses of the same construction and are fine because they test the linear-system uniqueness, not the construction itself.

**Verification:**

After the fix, the script should still exit 0. The six new `assert_zero` lines must appear in the source and the output transcript should now also confirm `generic K1 vs lane *` and `generic Hev vs lane *` indirectly (by virtue of the script exiting 0 and the existing prose lines for K1_gate_* and Hev_* remaining unchanged).

### F3 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no stage_017 script)

**What's wrong:**

Per the audit prompt's manifest entry, this unit is not a status-only candidate and not a checkpoint pass-through. Both engines are required. No `.wl` script exists for stage 017. The SymPy check is therefore single-engine, with no independent re-derivation of (a) the real-Y20 triple-overlap ratios (1, 1/2, -1), (b) the wall-only K1/H_even formulas, (c) the determinant 1/27 of the wall-only 2x2 even-gate system, or (d) the Xi_load and prefactor b=3a relations.

**Why this matters:**

Single-engine verification has been the source of several silent-pass failures in this project. The Y20 lane ratios in particular are computed entirely inside SymPy's `gaunt`; if SymPy's convention ever shifted (e.g., a Condon-Shortley sign change or a √(4π) normalization tweak), the SymPy script would still pass but no second source would catch it. Mathematica's `ThreeJSymbol[]` or `SphericalHarmonicY[]` provides an independent path.

**Required change:**

Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage017_parent_throat_action_weak_axisym_mathematica_audit.wl`. See the directive for the full claim manifest.

**Verification:**

The verifier runs `redteam exec-mathematica 017`; the new file must exit 0 and produce a transcript that lists each named claim with a PASS marker.

## Independent-derivation check (Mathematica)

Not applicable — no Mathematica script exists. See F3.

## Engine cross-check

Not applicable — no Mathematica script. See F3.

## Verdict justification

The script's m=1 and m=2 Gaunt ratios are real numerical checks against SymPy's `gaunt` and survive adversarial reading; so do the trace/b=3a relations for inertia, stiffness, Xi_load, and prefactor — those follow from the lane ratios via concrete arithmetic on the `grouped_trace_anomaly` formulas, and the determinant check at line 100 (det = 1/27) plus the sign-mutation guard at line 101 form a reasonable test of the wall-only 2x2 linear system. What does not hold up: the m=0 lane ratio is hard-coded to `Integer(1)` and never exercises gaunt, and the two "wall-only specialization" asserts at lines 79-80 just re-verify their own construction. Combined with the absence of a Mathematica engine, this earns three findings; none of them is `stop_cold` because the underlying math is correct on the substantive lanes (m=1, 2), and the fixes are local edits inside the script (plus creating a `.wl` partner).

## Self-test notes

I checked the three self-test traps and one additional one. (1) Variable independence: `wall_matrix` takes derivatives with respect to `dKsym` and `dMsym`; both K1_wall and Hev_wall genuinely depend on those symbols (1/9, -1, -1/27, 2/3 are the four entries), so no zero-derivative trap. The two trivial cases substituted mentally: `dKsym=9, dMsym=1` gives K1_wall = -1 + 1 = 0 and Hev_wall = 2/3 - 9/27 = 6/9 - 3/9 = 3/9 != 0, so K1_wall and Hev_wall are independent — the determinant 1/27 is consistent. (2) Symmetry/parity: not applicable — no integrals over an unbounded domain; all sums are finite Gaunt symbols. (3) Trivial-case pre-check for the proposed F2 cross-checks: with `lam20=1, lam21=1/2, lam22=-1`, `K1_gate_20 = eps*(K1w - 9*M1)/9 = eps*(-M1 + K1w/9)*1` matches `eps*lam20*(-M1+K1w/9)`; similarly `K1_gate_21 = eps*(K1w-9*M1)/18 = eps*(-M1+K1w/9)*(1/2)` matches `eps*lam21*generic_K1`; and `K1_gate_22 = eps*(-K1w+9*M1)/9 = -eps*(-M1+K1w/9) = eps*lam22*generic_K1`. The Hev rows match by the same arithmetic. (4) Path specifications for F3 directive: `.wl` lives in `/var/projects/toy_physics/research/pde_ledger/mathematica/`, confirmed by `ls`.
