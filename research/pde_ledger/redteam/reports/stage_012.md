---
unit_id: 012
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 012 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.txt`
- mathematica output: (missing)

## What the script claims to verify

The script perturbs the Stage-4/5 primitive one-port Maxwell/mixed data `(Q, S2, H, Delta, P, Gw)` by small linear slippages `(q1, s1, h1, d1, p1, g1)` weighted by `ell`, and derives the first-order corrections `z0, z2, z4, n0, n2, n4` to the bundle quantities `Z0, Z2, Z4, N0, N2, N4`. The script then claims (i) those linearizations equal the Frechet derivatives of the closed-form primitives, (ii) the static `Xi1 = n0/N0 + z0/D0` reduces to a stated explicit form, (iii) the isotropic compatibility surface obtained by eliminating `K` between the normalization and one-pole K surfaces is `(N0+ell n0)/Ptarget - 3*(S+ell z2)^2/(T+ell z4)`, with linear shift `n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2`, and the same for the transported-target variant, (iv) `z0` cancels from the compatibility surfaces but not from the normalization surface itself.

## Assertion inventory

| #   | Script | Line    | Form                                                                                          | Anchored to claim? |
|-----|--------|---------|-----------------------------------------------------------------------------------------------|--------------------|
| A1  | sympy  | 88      | `z0 - frechet(Z0) == 0`                                                                       | partial            |
| A2  | sympy  | 89      | `z2 - frechet(Z2) == 0`                                                                       | partial            |
| A3  | sympy  | 90      | `z4 - frechet(Z4) == 0`                                                                       | partial            |
| A4  | sympy  | 91      | `n0 - frechet(N0) == 0`                                                                       | partial            |
| A5  | sympy  | 92      | `n2 - frechet(N2) == 0`                                                                       | partial            |
| A6  | sympy  | 93      | `n4 - frechet(N4) == 0`                                                                       | partial            |
| A7  | sympy  | 97-100  | `Xi1_static - (2p1/P - 2d1/Delta + (Delta q1 - Q d1)/(D0 Delta^2)) == 0`                      | yes                |
| A8  | sympy  | 124     | `(K_norm_p - K_one_p) - compat_direct_p == 0`                                                 | no (tautological)  |
| A9  | sympy  | 125     | `dCompat - dCompat_direct == 0`                                                               | no (tautological)  |
| A10 | sympy  | 126-129 | `dCompat_direct - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2) == 0`                                | yes                |
| A11 | sympy  | 130     | `K_norm_transport_p - (B0 + Z0slot + ell z0 + D0target) == 0`                                 | partial            |
| A12 | sympy  | 132-134 | `compat_transport_p - (D0target - 3 (S+ell z2)^2/(T+ell z4)) == 0`                            | partial            |
| A13 | sympy  | 136-138 | `dCompat_transport - (-6 S z2/T + 3 S^2 z4/T^2) == 0`                                         | yes                |
| A14 | sympy  | 143     | `diff(K_norm_probe - K_one_probe, z0_probe) == 0`                                             | no (tautological)  |
| A15 | sympy  | 144-147 | `diff(K_norm_transport_probe - K_one_probe, z0_probe) == 0`                                   | no (tautological)  |
| A16 | sympy  | 148     | `assert_nonzero diff(K_norm_probe, z0_probe)` (= ell)                                         | no (tautological)  |
| A17 | sympy  | 149-152 | `assert_nonzero` mutated dCompat_direct with wrong sign on z4 term                            | yes (cheap)        |
| A18 | sympy  | 153-156 | `assert_nonzero` mutated dCompat_transport with wrong sign on z4 term                         | yes (cheap)        |

## Findings

### F1 — missing_verification_script

**Severity:** high
**Files:**
- `(missing)` — no Mathematica script for unit 012

**What's wrong:**
The manifest entry for unit 012 has `is_status_only_candidate: False`, so it is a real verification unit and both engines are required. Only a SymPy script exists at `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py`. There is no `.wl` companion. Consequently `engines_agree = n/a`, and every claim in this stage is verified by a single algebra system.

**Why this matters:**
The unit derives closed-form first-order corrections `z0..z4, n0..n4` and explicit compatibility-shift formulas that downstream units consume. With only one engine, an error in sympy's `simplify`/`series`/`solve` (or in the typed-in closed forms) would not be caught.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl` that independently re-derives — not transliterates — the claims in the manifest below, using Mathematica primitives (`Series`, `Coefficient`, `D`, `Solve`, `FullSimplify`). It must terminate with an explicit pass/fail (e.g. `If[FullSimplify[...] =!= 0, Print["FAIL"]; Exit[1]]`) for each claim.

**Verification:**
After Codex creates the `.wl`, the verifier runs `redteam exec-mathematica 012` and expects (i) the file to exist, (ii) the script to exit 0, (iii) all eight claim checks listed in the directive's claim manifest to print PASS lines.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py:124`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py:125`

**What's wrong:**
At line 103-106 the script solves `(K - B0 - (Z0slot + ell*z0))*(T + ell*z4) = 3*(S + ell*z2)^2` for `K` to obtain `K_one_p = B0 + Z0slot + ell*z0 + 3*(S+ell*z2)^2/(T+ell*z4)`. At line 107-110 it solves `(N0base+ell*n0)/(K - B0 - (Z0slot+ell*z0)) = Ptarget` for `K` to obtain `K_norm_p = B0 + Z0slot + ell*z0 + (N0base+ell*n0)/Ptarget`. At line 111 it defines `compat_direct_p = (N0base+ell*n0)/Ptarget - 3*(S+ell*z2)^2/(T+ell*z4)`. By construction `K_norm_p - K_one_p == compat_direct_p` algebraically — there is no physics in the equality. The assertion at line 124, `(K_norm_p - K_one_p) - compat_direct_p == 0`, cannot fail regardless of any error elsewhere. Line 125's `dCompat - dCompat_direct == 0` is a direct consequence of line 124 (series expansion is linear), so it is equally tautological.

**Why this matters:**
These checks consume audit budget while testing nothing about z0,z2,z4,n0 correctness; a reader inspecting the script's PASS line is misled into thinking the compatibility-surface derivation is independently verified when it is only a rename of a difference.

**Required change:**
Replace the line 124 assertion with one that exercises the compatibility surface against an independent derivation — concretely, expand `K_norm_p - K_one_p` directly and compare its closed form, term by term, to the right-hand side `(N0base + ell*n0)/Ptarget - 3*(S + ell*z2)**2/(T + ell*z4)` only after substituting the solved K's back into the *original* defining equations (e.g. verify `(K_norm_p - B0 - (Z0slot + ell*z0))*Ptarget - (N0base + ell*n0) == 0` and `(K_one_p - B0 - (Z0slot + ell*z0))*(T + ell*z4) - 3*(S + ell*z2)**2 == 0`). Then delete the line 125 assertion (it duplicates A10 at line 126-129 once A10 is the substantive check).

**Verification:**
After patching, lines 124-125 either (i) are removed entirely in favor of two new "solve round-trip" checks, or (ii) read as direct round-trip residual checks against the original solve equations. The output file should no longer print the labels `primitive compatibility surface after eliminating K` and `primitive compatibility shift from eliminated surface` unless those labels now front substantive checks.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py:139-148`

**What's wrong:**
At lines 140-142 the script types `K_norm_probe = B0 + Z0slot + ell*z0_probe + (N0base + ell*n0)/Ptarget`, `K_one_probe = B0 + Z0slot + ell*z0_probe + 3*(S + ell*z2)^2/(T + ell*z4)`, `K_norm_transport_probe = B0 + Z0slot + ell*z0_probe + D0target`. The `ell*z0_probe` term appears identically in all three. The assertions at lines 143 and 144-147, `diff(K_norm_probe - K_one_probe, z0_probe) == 0` and `diff(K_norm_transport_probe - K_one_probe, z0_probe) == 0`, are forced by construction — `z0_probe` cancels because the same literal `ell*z0_probe` summand was typed into both terms of each difference. The assertion at line 148, `assert_nonzero diff(K_norm_probe, z0_probe)`, evaluates to `ell != 0`, which is also guaranteed by the typed-in form. None of these touches `z0` (the derived bundle correction) at all.

**Why this matters:**
The script's stage-5 prose claims "z0 cancels from both the fixed-target and transported-target compatibility shifts, even though the normalization K surface still carries z0 before compatibility elimination." That claim is about the derived `z0 = (Delta*q1 - Q*d1)/Delta^2`, not about a free symbol `z0_probe` that the author typed identically into both sides. The current test cannot detect any of the errors it appears to guard against (e.g. a sign error in z0 derivation, or a Ptarget-vs-Ptarget_transport substitution slip).

**Required change:**
Replace lines 139-148 with a check that uses the actual derived `z0` (from line 76). Concretely, substitute the explicit `z0 = (Delta*q1 - Q*d1)/Delta**2` into `K_norm_p` and `K_one_p` from lines 103-110, then verify `diff(K_norm_p - K_one_p, q1) == diff(compat_direct_p, q1)` and similarly for `d1` — i.e. the q1 and d1 partials of the compatibility shift must come only from `n0` and `z2/z4`, not from `z0`. Or, equivalently, verify that after substituting the derived `z0`, the symbolic dependence of `K_norm_p - K_one_p` on `q1, d1` arises solely through `n0` and the (z2,z4) of the compatibility surface, with the `z0` channel contributing zero. Then assert that `K_norm_p` (alone, no subtraction) does depend on `q1` and `d1` through the explicit `z0` channel — which is the actual non-tautological version of A16.

**Verification:**
After patching, lines 139-148 reference the bound symbol `z0` (and through it `q1, d1`) rather than an unbound `z0_probe`. The output transcript prints a check label such as `primitive compatibility surface has no q1/d1 channel from z0` and a complementary `primitive normalization K surface retains q1/d1 channel from z0`, both with PASS.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.py:73-93`

**What's wrong:**
The Frechet block (A1-A6) compares `dlin(expr)` (via `sp.series(expr.subs(subs), ell, 0, 2).removeO() - expr) / ell`) against `frechet(expr)` (via `sum(diff(expr, var) * slip)`). Both are sympy-internal first-order Taylor coefficients of the same polynomial-rational `expr`; they must agree for any well-formed `expr` regardless of whether the typed-in closed forms `Z0..N4` are physically correct. The script's docstring says the unit "derives the induced bundle corrections z0, z2, z4, n0, n2, n4 exactly at first order" — but A1-A6 verify only that two equivalent sympy methods give the same answer, not that the resulting `z0..n4` match any externally-anchored closed form.

**Why this matters:**
If a typo or sign error were introduced into `Z2`, `Z4`, `N2`, or `N4` at lines 57-62, A1-A6 would still pass — both `dlin` and `frechet` would silently linearize the wrong primitive. Section 2 of the output transcript (the closed forms for `z0..n4`) is therefore the actual claim, but nothing checks those closed forms against an independent derivation.

**Required change:**
Add explicit closed-form checks for `z0, z2, z4, n0, n2, n4` against typed-in expected expressions, mirroring the style used at line 99 for `Xi1_static`. For example, after line 81, add:
```
assert_zero("z0 closed form",  z0 - (Delta*q1 - Q*d1)/Delta**2)
assert_zero("z2 closed form",  z2 - (-Delta**2*h1 + Delta*(H*d1 + Q*s1 + S2*q1) - 2*Q*S2*d1)/Delta**3)
assert_zero("z4 closed form",  z4 - (-Delta**2*H*s1 - Delta**2*S2*h1 - Delta**2*q1 + 2*Delta*H*S2*d1 + 2*Delta*Q*S2*s1 + 2*Delta*Q*d1 + Delta*S2**2*q1 - 3*Q*S2**2*d1)/Delta**4)
assert_zero("n0 closed form",  n0 - 2*P*(Delta*p1 - P*d1)/Delta**3)
assert_zero("n2 closed form",  n2 - (-(2*Delta**2*(Gw*p1 + P*g1) - 2*Delta*P*(2*Gw*d1 + P*s1 + 2*S2*p1) + 6*P**2*S2*d1)/Delta**4))
assert_zero("n4 closed form",  n4 - 2*(Delta**3*Gw*g1 - Delta**2*Gw**2*d1 - 2*Delta**2*Gw*P*s1 - 2*Delta**2*Gw*S2*p1 - 2*Delta**2*P*S2*g1 - 2*Delta**2*P*p1 + 6*Delta*Gw*P*S2*d1 + 3*Delta*P**2*S2*s1 + 3*Delta*P**2*d1 + 3*Delta*P*S2**2*p1 - 6*P**2*S2**2*d1)/Delta**5)
```
These RHSs are exactly the closed forms the script prints in section 2 of the output file (transcript lines 36-47). Replacing or supplementing A1-A6 with these elevates them from method-consistency to closed-form anchors. Keep A1-A6 if desired, but add the six closed-form anchors.

**Verification:**
After patching, the output transcript prints six new lines labelled `z0 closed form` ... `n4 closed form`, each passing.

## Independent-derivation check (Mathematica)

No `.wl` exists. See finding F1. There is no second engine to compare against.

## Engine cross-check

`engines_agree = n/a` — only sympy is present.

## Verdict justification

The SymPy script is internally consistent and many of its assertions are substantive (A7, A10, A13, the closed-form Xi1 check, and the mutated-sign nonzero sanity checks at lines 149-156). However, four substantive defects exist: (1) the Mathematica engine is entirely absent for a non-status-only unit; (2) the "compatibility surface = K_norm - K_one" check at line 124 and its derivative at line 125 are tautological renamings; (3) the entire z0_probe block at lines 139-148 uses a free symbol identically typed on both sides of each difference and exercises no physics; (4) the Frechet block A1-A6 verifies only sympy-internal method consistency, not that the closed forms `z0..n4` printed in the output are correct. Attacks I tried that failed: I checked whether `Z0..N4` and the Frechet vs series outputs would actually disagree under a typo (they would, but only A1-A6 catches structural mismatches, not closed-form errors); I checked whether the compatibility-shift linear coefficients `n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2` and `-6 S z2/T + 3 S^2 z4/T^2` are computed correctly (they are — series expansion of `3*(S+ell*z2)^2/(T+ell*z4)` at ell=0 has linear coefficient `6 S z2/T - 3 S^2 z4/T^2`, with the sign discipline carried correctly into A10 and A13); I checked symbol assumptions (`Q, S2, H, Delta, P, Gw, K, B0, Z0slot, N0base, D0target, Ptarget, S, T, D0sym` are all declared `nonzero=True`, no further assumptions made, no divide-by-zero hazard given the rational structure). The unit is repairable; verdict is `findings`, not `stop_cold`.
