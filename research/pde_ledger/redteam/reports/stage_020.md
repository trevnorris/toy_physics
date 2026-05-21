---
unit_id: 020
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 020 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

The script asserts that, with the parent throat action combined with the support/mixed-slope packets, the two "even gates" K1 = D21 + D01/9 and H_even = D41 - (2/3) D21 - D01/27 uniquely fix the parent wall slopes (dKSigma, dMSigma). Specifically, it claims the 2x2 coefficient matrix has determinant 1/27, the solved slopes equal dKSigma = B01+Z01+27(B41+Z41) and dMSigma = -(B21+Z21)+3(B41+Z41), and on this compensated branch the residual deficits collapse to D01 = 27(B41+Z41), D21 = -3(B41+Z41), D41 = -(B41+Z41), and Xi1 = N01/N0 - 27(B41+Z41)/(KSigma-B0-Z0). In a separate isolated block it claims the real-Y20 self-overlap ratios are (lambda_0, lambda_1, lambda_2) = (1, 1/2, -1) using sympy's Gaunt coefficients, and that the same-sign m≠0 cross terms vanish.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 49 | `assert_zero('Y20 overlap lane 20', lam20 - 1)` | no (tautological — see F2) |
| A2 | sympy | 50 | `assert_zero('Y20 overlap lane 21', lam21 - sp.Rational(1, 2))` | yes |
| A3 | sympy | 51 | `assert_zero('Y20 overlap lane 22', lam22 + 1)` | yes |
| A4 | sympy | 57 | `assert_zero('even-gate solve determinant', coeff_matrix.det() - sp.Rational(1, 27))` | yes |
| A5 | sympy | 66 | `assert_zero('dKSigma closed form', sol[dKSigma] - (B01 + Z01 + 27*(B41 + Z41)))` | yes |
| A6 | sympy | 67 | `assert_zero('dMSigma closed form', sol[dMSigma] - (-(B21 + Z21) + 3*(B41 + Z41)))` | yes |
| A7 | sympy | 68 | `assert_zero('compensated D01', D01_comp - 27*(B41 + Z41))` | yes |
| A8 | sympy | 69 | `assert_zero('compensated D21', D21_comp + 3*(B41 + Z41))` | yes |
| A9 | sympy | 70 | `assert_zero('compensated D41', D41_comp + B41 + Z41)` | yes |
| A10 | sympy | 71 | `assert_zero('compensated Xi1', Xi1_comp - (N01/N0 - 27*(B41 + Z41)/(KSigma - B0 - Z0)))` | yes |

The same-sign vanishing checks at line 26 (raise AssertionError) for m=1,2 are additional non-tautological Gaunt-coefficient assertions and were considered when scoring A2/A3.

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl` (does not exist)

**What's wrong:**
The unit's manifest entry has `is_status_only_candidate: False` and the unit is not a checkpoint, so both engines are required by the second-engine policy. Only the SymPy script exists. The `mathematica/` directory contains no file matching the stage 020 pattern (the listing skips from `stage012` to `stage021`). No upstream unit's `.wl` script independently re-derives the K1/H_even solve, the determinant 1/27, the compensated branch identities, or the Xi1 closed form that this SymPy script asserts; consequently, every result the SymPy script claims is currently single-engine.

Subtype: `missing_mathematica`.

**Why this matters:**
Without an independent Mathematica derivation, a SymPy-only algebraic identity (e.g., dKSigma = B01 + Z01 + 27(B41+Z41), determinant = 1/27) is unverified against a second symbolic engine. Sign/factor errors in the K1, H_even definitions (e.g., a swapped sign on -D01/27, an off-by-3 factor in 2/3 D21) would be invisible because the SymPy script self-consistently propagates them. The Gaunt-coefficient ratios (1, 1/2, -1) for the real-Y20 overlap also rely on a single library's convention.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl` that independently derives, from the symbolic definitions of D0, D01, D21, D41 and the even-gate combinations K1, H_even, the same identities the SymPy script claims (see Claim manifest in the directive). The Mathematica script must use `Solve` (not import or hard-code the SymPy output), perform its own `Det` on the Jacobian, and verify each compensated identity by `FullSimplify[... == 0]` followed by an explicit `If[!TrueQ[...], Print["FAIL ..."]; Exit[1]]`. It must NOT be a line-by-line transliteration of the SymPy script.

**Verification:**
After the new `.wl` lands and `redteam exec-mathematica 020` runs, the saved output `mathematica/output/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.txt` must show: determinant = 1/27, dKSigma = B01+Z01+27(B41+Z41), dMSigma = -(B21+Z21)+3(B41+Z41), D01_comp = 27(B41+Z41), D21_comp = -3(B41+Z41), D41_comp = -(B41+Z41), Xi1_comp = N01/N0 - 27(B41+Z41)/(KSigma-B0-Z0), and the Gaunt-ratio claims (1, 1/2, -1). Exit code 0.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:20-28, 46, 49`

**What's wrong:**
`real_y20_square_ratio(0)` short-circuits before touching `gaunt`:

```python
def real_y20_square_ratio(m: int) -> sp.Expr:
    base = sp.simplify(gaunt(2, 2, 2, 0, 0, 0))
    if m == 0:
        return sp.Integer(1)
    ...
```

So `lam20 = real_y20_square_ratio(0)` returns the literal `1` without any Gaunt-coefficient evaluation in the m=0 branch. The subsequent check `assert_zero('Y20 overlap lane 20', lam20 - 1)` on line 49 then evaluates `1 - 1`, which is algebraically guaranteed to be zero regardless of whether the real-Y20 m=0 self-overlap actually equals `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)`. The m=1 and m=2 branches DO call `gaunt(...)` for both numerator and denominator and are non-tautological; only the m=0 lane is a self-check of a constant.

**Why this matters:**
The script's `Y20 overlap lane 20` assertion appears in the inventory as evidence that the m=0 real-Y20 self-overlap matches the (1, 1/2, -1) pattern. As written, that line cannot fail even if Gaunt's m=0 self-overlap were broken or sign-flipped in the underlying library. The bottom-line claim "the real-Y20 self-overlap ratios at m=0,1,2 are 1, 1/2, -1" is thus only verified at m=1, 2.

**Required change:**
Make the m=0 branch also compute the ratio explicitly so the assertion exercises `gaunt(2,2,2,0,0,0)` non-tautologically. Replace the body of `real_y20_square_ratio` at lines 20-27 with a uniform formula that does not branch on `m==0`. Concretely, remove the `if m == 0: return sp.Integer(1)` short-circuit and use the same `(-1)^m * gaunt(2,2,2,0,m,-m) / base` expression for all m in {0, 1, 2} (when m=0 this becomes `gaunt(2,2,2,0,0,0)/gaunt(2,2,2,0,0,0)` which equals 1 only if the Gaunt evaluation is consistent). The `same_sign` cross-term guard should be retained for m != 0 only.

**Verification:**
After Codex applies, `redteam exec-sympy 020` should still exit 0 with the printed STATUS: PASS, and the body of `real_y20_square_ratio` at lines 20-27 should no longer contain the `if m == 0: return sp.Integer(1)` short-circuit. The function should call `gaunt(2,2,2,0,0,0)` once and compute the ratio uniformly.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists for this unit. This is captured by F1 (`missing_mathematica`).

## Engine cross-check

Not applicable — only the SymPy engine is present. See F1.

## Verdict justification

The SymPy script's core algebra (K1, H_even, determinant 1/27, the closed-form solve, and the compensated-branch identities) is internally consistent and the assertions are non-tautological — I checked the determinant by hand from the partial derivatives ((1/9)(2/3) - (-1)(-1/27) = 2/27 - 1/27 = 1/27) and confirmed the compensated branch substitutions match the asserted closed forms by direct substitution. The Gaunt-coefficient ratios (1, 1/2, -1) for m=0,1,2 are physically correct for the real-Y20 self-overlap. However: (1) no Mathematica counterpart exists, making this a single-engine result for a non-status-only, non-checkpoint unit; (2) the m=0 Gaunt lane short-circuits to a literal `1`, so its `assert_zero(lam20 - 1)` is tautological. Verdict: `findings`. Not `stop_cold` — both findings are correctable without invalidating any downstream result (a new `.wl` re-derivation is additive, and the F2 fix only tightens an existing check).

## Self-test notes

- Variable independence: A4's Jacobian uses `sp.diff(K1, dKSigma)`, `sp.diff(K1, dMSigma)`, etc.; both K1 and H_even genuinely depend on dKSigma and dMSigma (via D01 and D21 respectively), so no zero-derivative trap. Confirmed by re-reading lines 42-43.
- Symmetry/parity: no integrals over unbounded domains in this unit; trap not applicable.
- Trivial-case pre-check: substituting B01=Z01=B21=Z21=B41=Z41=0 gives dKSigma=dMSigma=0, D01=D21=D41=0, Xi1 = N01/N0 — all assertions reduce to 0 = 0, as expected. Substituting B41=1, others=0 gives dKSigma=27, dMSigma=3, D01=27, D21=-3, D41=-1, matching the asserted closed forms in A7-A9.
- Path: the new `.wl` target path uses the `mathematica/` directory and the `_mathematica_audit.wl` suffix per the convention visible in the directory listing (stages 001-012, 021+). Confirmed against `ls mathematica/`.
- F2 fix self-test: after removing the `if m==0` short-circuit, `real_y20_square_ratio(0)` evaluates `(-1)^0 * gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0) = 1`, which still passes `assert_zero('... lane 20', lam20 - 1)`. The `same_sign` check guard must be skipped when m=0 (since gaunt(2,2,2,0,0,0) is the base itself, not a cross term); the directive notes this explicitly.
