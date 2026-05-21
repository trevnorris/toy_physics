---
unit_id: 018
batch: I.2
auditor_model: claude-opus-4-7
audit_date: 2026-05-21T12:29:31-06:00
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 018 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

Per the docstring ("Master-note audit for step_16_parent_throat_action_bundle_master_notes.md") and the inline assertions, the script checks three things on the parent-throat action bundle. First, an "isotropic compatibility" block: it derives `u4 - 4*u2^2 = (D0*(B4+Z4) - 3*(MSigma+B2+Z2)^2)/D0^2` and confirms two ways of solving for `KSigma` (one-pole closure and normalization to `N0/Ptarget`) are mutually consistent. Second, an "exact weak-axisymmetric wall-slope solve": linearizing in `dKSigma, dMSigma` it computes a 2x2 coefficient matrix whose determinant is `1/27`, then solves the gate equations `K1 = 0`, `H_even = 0` for the first-order slopes and reports `dKSigma = B01+Z01+27(B41+Z41)`, `dMSigma = -(B21+Z21)+3(B41+Z41)`, along with the residual amplitude `Xi1`. Third, a concrete Gaussian-profile check that the integrals defining `MSigma` and `KSigma` reduce to `sqrt(pi)` and `3*sqrt(pi)/2` respectively. Mutation rows (rows that should *not* simplify to zero) are paired with several of the substantive checks.

## Assertion inventory

| #  | Script | Line  | Form | Anchored to claim? |
|----|--------|-------|------|--------------------|
| A1 | sympy  | 30-33 | `assert_zero (u4 - 4*u2^2) - (D0*(B4+Z4) - 3*(M+B2+Z2)^2)/D0^2` | yes |
| A2 | sympy  | 37    | `assert_zero (u4 - 4*u2^2).subs(KSigma, K_from_one_pole)` | yes |
| A3 | sympy  | 38    | `assert_zero (N0/D0).subs(KSigma, K_from_norm) - Ptarget` | yes (definitional) |
| A4 | sympy  | 39    | `assert_nonzero (u4 - 3*u2^2).subs(KSigma, K_from_one_pole)` | yes (mutation) |
| A5 | sympy  | 40    | `assert_nonzero (N0/D0).subs(K_from_norm) - 2*Ptarget` | yes (mutation) |
| A6 | sympy  | 42    | `assert_zero simplify(K_from_norm - K_from_one_pole) - compatibility` | **no — tautology** |
| A7 | sympy  | 59    | `assert_zero coeff_matrix.det() - 1/27` | yes |
| A8 | sympy  | 64    | `assert_zero sol[dKSigma] - (B01+Z01+27(B41+Z41))` | yes |
| A9 | sympy  | 65    | `assert_zero sol[dMSigma] - (-(B21+Z21)+3(B41+Z41))` | yes |
| A10| sympy  | 66    | `assert_nonzero coeff_matrix.det() + 1/27` | yes (mutation) |
| A11| sympy  | 69    | `assert_zero D01_comp - 27(B41+Z41)` | yes (corollary of A8) |
| A12| sympy  | 70    | `assert_zero K1.subs(sol)` | **no — tautology** |
| A13| sympy  | 71    | `assert_zero H_even.subs(sol)` | **no — tautology** |
| A14| sympy  | 72-75 | `assert_zero Xi1.subs(sol) - (N01/N0 - 27(B41+Z41)/(K-B0-Z0))` | partial (corollary of A11) |
| A15| sympy  | 76-79 | `assert_nonzero Xi1.subs(sol) - (N01/N0 + 27(B41+Z41)/(K-B0-Z0))` | yes (mutation) |
| A16| sympy  | 85    | `assert_zero MSigma_example - sqrt(pi)` | yes |
| A17| sympy  | 86    | `assert_zero KSigma_example - 3*sqrt(pi)/2` | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:41-42`

**What's wrong:**
The script defines

```
compatibility = N0 / Ptarget - 3 * (MSigma + B2 + Z2) ** 2 / (B4 + Z4)
assert_zero("compatibility equality", sp.simplify(K_from_norm - K_from_one_pole) - compatibility)
```

But three lines above, `K_from_norm = B0 + Z0 + N0/Ptarget` and `K_from_one_pole = B0 + Z0 + 3*(MSigma+B2+Z2)**2/(B4+Z4)`. Their difference, by elementary cancellation of the common `B0 + Z0`, is exactly `N0/Ptarget - 3*(MSigma+B2+Z2)**2/(B4+Z4)`, which is literally how `compatibility` is defined on the preceding line. The `assert_zero` reduces to `expr - expr == 0`, which is algebraically guaranteed by construction and cannot fail no matter what the physics is. It does not verify "isotropic compatibility" — that compatibility content was already baked into the choice of `compatibility`'s definition.

**Why this matters:**
The docstring/print line advertises an "isotropic compatibility" check, but A6 contributes nothing to that verification. The substantive compatibility work is done by A2 and A3 (each closure gives the same `Ptarget` / one-pole residual). A6 must either (a) be replaced with a substantive check that the two `KSigma` closures imply a non-trivial constraint, or (b) be removed.

**Required change:**
Replace the tautology with a genuine cross-closure check. The natural one: when both closures hold simultaneously, the *equality* `K_from_norm = K_from_one_pole` is equivalent to `compatibility = 0`. Assert that — but treat `compatibility` as the *content* to verify, not as a relabel of the difference. Concretely, solve `K_from_norm - K_from_one_pole = 0` for `N0` (or for any one of `B4+Z4`, `MSigma+B2+Z2`, `Ptarget`) and assert the solution equals what `compatibility = 0` predicts.

Suggested replacement at line 41-42:

```python
compatibility = N0 / Ptarget - 3 * (MSigma + B2 + Z2) ** 2 / (B4 + Z4)
# When both K closures hold, compatibility=0 must follow.
N0_from_compat = sp.solve(compatibility, N0)[0]
N0_from_equality = sp.solve(K_from_norm - K_from_one_pole, N0)[0]
assert_zero("compatibility equality", N0_from_compat - N0_from_equality)
assert_nonzero("mutated compatibility equality should fail", N0_from_compat - 2 * N0_from_equality)
```

**Verification:**
The line at 42 should no longer expand to `compatibility - compatibility == 0`. The verifier reruns `redteam exec-sympy 018`; the new check should appear in the script and the exit code stays 0.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:61, 70-71`

**What's wrong:**
At line 61 the script runs `sp.solve([sp.Eq(K1, 0), sp.Eq(H_even, 0)], [dKSigma, dMSigma], dict=True)[0]`, then at lines 70-71 it asserts

```
assert_zero("compensated K1", K1.subs(sol))
assert_zero("compensated H_even", H_even.subs(sol))
```

By construction, `sol` is the solution that makes `K1` and `H_even` both zero — that is `sp.solve`'s contract. Substituting `sol` back into `K1` and `H_even` is algebraically guaranteed to give zero (modulo simplification), so neither assertion can fail. Both rows are tautological with the `sp.solve` call that produced `sol`.

**Why this matters:**
These two assertions add no verification value to the wall-slope-solve claim. The substantive checks for that block are A8 and A9 (slopes equal the expected closed forms), the determinant in A7, the mutation rows A10 and A15, and the compensated `D01` check A11. The two tautological rows obscure where the real coverage is.

**Required change:**
Replace each tautological row with a substantive consistency or mutation check. Two natural replacements:

1. Cross-check that the closed-form expected slopes (`expected_dK`, `expected_dM`) actually satisfy `K1 = 0` and `H_even = 0` when substituted directly (no `sp.solve` round-trip):

```python
assert_zero("expected_dK satisfies K1", K1.subs({dKSigma: expected_dK, dMSigma: expected_dM}))
assert_zero("expected_dK satisfies H_even", H_even.subs({dKSigma: expected_dK, dMSigma: expected_dM}))
```

2. Pair them with mutation checks (one is enough but two is preferred): `assert_nonzero` for `K1.subs({dKSigma: expected_dK + 1, dMSigma: expected_dM})` and similar for `H_even`.

Concretely, replace lines 70-71 with:

```python
assert_zero("expected_dK satisfies K1", K1.subs({dKSigma: expected_dK, dMSigma: expected_dM}))
assert_zero("expected_dK satisfies H_even", H_even.subs({dKSigma: expected_dK, dMSigma: expected_dM}))
assert_nonzero("mutated expected_dK fails K1", K1.subs({dKSigma: expected_dK + 1, dMSigma: expected_dM}))
```

**Verification:**
The verifier reruns `redteam exec-sympy 018`; the new substantive checks should appear in place of the tautological ones and the script should exit 0.

### F3 — missing_verification_script

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl` (must be created)

**What's wrong:**
No Mathematica audit script exists for stage 018. The manifest lists this unit with `is_checkpoint: false` and `is_status_only_candidate: false`, so the two-engine policy applies: an independent Mathematica re-derivation of the same physical claims is required. The unit's SymPy script encodes seven substantive physical results (one-pole numerator identity, one-pole `KSigma` closure, normalization closure consistency, even-gate determinant = `1/27`, the two wall-slope formulas, and the Gaussian-profile integrals for `MSigma` and `KSigma`) — none of these are independently checked.

**Why this matters:**
A single-engine verification of a non-checkpoint, non-status-only unit cannot certify the unit. The whole point of the second engine is to catch SymPy-specific bugs (incorrect simplification under aggressive assumptions, branch-cut surprises, accidental tautologies of the kind already flagged in F1 and F2) by re-deriving from the physical premises in a different CAS.

**Required change:**
Create the file at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl`. The script must independently re-derive (not transliterate) and verify the following claims; see the Claim manifest in the directive for symbolic forms.

**Verification:**
The verifier runs `redteam exec-mathematica 018`. The new `.wl` file must be present, must contain assertions (`If[FullSimplify[...] =!= 0, Exit[1]]` or equivalent), and must exit with code 0. The verifier will also check that the assertions are independent of the SymPy script's intermediate choreography (see `mathematica_transliteration` category).

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:72-75`

**What's wrong:**
The assertion at lines 72-75 is

```
assert_zero(
    "residual Xi1",
    Xi1.subs(sol) - (N01 / N0 - 27 * (B41 + Z41) / (KSigma - B0 - Z0)),
)
```

`Xi1` is defined as `N01/N0 - D01/D0`. Given A11 has already established `D01.subs(sol) = 27*(B41+Z41)` and `D0 = KSigma - B0 - Z0`, this assertion follows by elementary substitution from A11. It is not literally tautological with the `sp.solve` call (the right-hand side `27*(B41+Z41)/(KSigma-B0-Z0)` is a non-trivial form), but it does not exercise any claim beyond A11. The paired mutation A15 partially compensates by checking the sign.

**Why this matters:**
This is the only check that anchors `Xi1` to the wall-slope solution — it is supposed to certify that the *amplitude residual* of the perturbation is `N01/N0` minus the closed-form expression for the inhomogeneous correction. If the only verification reuses `D01_comp` from A11, then a SymPy bug in the `D01.subs(sol)` simplification path propagates silently through both A11 and A14.

**Required change:**
Add (do not replace; A14 is fine as is) an independent recomputation of `Xi1.subs(sol)` from a different path. Concretely, recompute `Xi1` after substituting the *expected closed-form* slopes (`expected_dK`, `expected_dM`) rather than the `sp.solve` output `sol`, and assert the result equals the same RHS:

Insert after line 79:

```python
Xi1_from_expected = Xi1.subs({dKSigma: expected_dK, dMSigma: expected_dM})
assert_zero(
    "residual Xi1 from expected slopes",
    Xi1_from_expected - (N01 / N0 - 27 * (B41 + Z41) / (KSigma - B0 - Z0)),
)
```

This routes through the closed-form expected slopes rather than `sp.solve`'s output, giving an independent path for the same claim.

**Verification:**
The verifier reruns `redteam exec-sympy 018`; the new check should appear after the existing mutation row and the script should exit 0.

## Independent-derivation check (Mathematica)

No `.wl` script exists. See F3.

## Engine cross-check

Not applicable. Only the SymPy engine is present. See F3.

## Verdict justification

The SymPy script is mostly sound and the substantive claims (isotropic one-pole numerator identity, both `KSigma` closures, even-gate determinant = `1/27`, the wall-slope closed forms, the Gaussian integrals) hold up under attack: I verified the algebra in each case (e.g., `u4 - 4*u2^2` reduction, the 2x2 determinant computed by hand, and the Gaussian integrals `int exp(-w^2) = sqrt(pi)`, `int w^2*exp(-w^2) = sqrt(pi)/2`). The mutation rows (A4, A5, A10, A15) correctly catch sign and coefficient errors. However, three assertions are tautological by construction (A6, A12, A13) — they cannot fail regardless of the physics — and A14 is a corollary of A11 along a single path. More importantly, the unit is non-checkpoint and non-status-only, so the missing Mathematica script is a high-severity gap that the two-engine policy requires closed. Verdict is `findings`, not `clean` and not `stop_cold`.

## Self-test notes

I checked the four trap categories. (1) Variable independence: F2's proposed `K1.subs({dKSigma: expected_dK, dMSigma: expected_dM})` is well-formed because `K1` depends on both `dKSigma` and `dMSigma`, so the substitution is non-trivial; the assertion vanishes because that is exactly the algebraic identity we want and I verified it by hand (`K1 = -(dM+B21+Z21) + (dK-B01-Z01)/9` with `dK = B01+Z01+27(B41+Z41)` and `dM = -(B21+Z21)+3(B41+Z41)` gives `-3(B41+Z41) + 27(B41+Z41)/9 = -3(B41+Z41) + 3(B41+Z41) = 0`). (2) Parity: F1's proposed `sp.solve(compatibility, N0)` is straightforward and not an integral, so parity is not at issue; the Gaussian-integral assertions (A16, A17) involve `beta^2` (even) and `(d beta / dw)^2` (also even since the square of an odd function), so both integrate to nonzero finite values — already correctly handled. (3) Trivial-case pre-check: for F1, `compatibility = 0` implies `N0/Ptarget = 3(M+B2+Z2)^2/(B4+Z4)`, and `K_from_norm - K_from_one_pole = 0` implies the same identity for `N0`, so the two `N0` solutions match; mutation by factor 2 obviously fails. (4) Path specification: F3's target path is `mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl`, matching the naming pattern of the other `.wl` files I listed under `mathematica/`.
