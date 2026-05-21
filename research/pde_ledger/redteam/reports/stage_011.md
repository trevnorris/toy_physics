---
unit_id: 011
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 011 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.txt`
- mathematica output: `(missing)`

## What the script claims to verify

The script connects projected-Maxwell near-throat corrections to the grouped real-P2 / full-bundle bookkeeping used by the moving-throat PDE program. It perturbs the isotropic bundle moments by `D_n -> D_n - eps*z_n`, `N_n -> N_n + eps*n_n`, derives the induced first-order shifts in `u2, u4, P0, P2, P4`, and asserts closed-form expressions for each. It also verifies (a) the first variation of the one-pole defect, (b) that the isotropic compatibility surface between the one-pole and normalization K-surfaces is independent of `z0` after K-elimination, both for fixed and transported targets, (c) a constant-prefactor branch factorization, (d) real-Y20 weak-axisymmetric grouped lane ratios `(1, 1/2, -1)` with `xbar = x0`, `b = 3 a`, and (e) the projected-Maxwell static contributions `Xi1^(proj,static) = na/N0 + za/D0` and `du2/d(ea*lam) = (D0*z2a - D2*za)/D0^2`. Two negative-control assertions (`assert_nonzero` with deliberately sign-flipped target expressions) are also exercised.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 102 | `assert_zero("delta u2", du2 - (D0*z2 - D2*z0)/D0**2)` | yes |
| A2 | sympy | 103 | `assert_zero("delta P0", dP0 - (D0*n0 + N0*z0)/D0**2)` | yes |
| A3 | sympy | 113 | `assert_zero("one-pole first variation", d_pole - (D0*z4 - z0*(B4+Z4) - 6*z2*(M+B2+Z2)))` | yes |
| A4 | sympy | 141 | `assert_zero("compatibility surface after eliminating K", (K_norm_p - K_one_pole_p) - compat_direct_p)` | yes |
| A5 | sympy | 142 | `assert_zero("one-pole K shift", dK_one_pole - (z0 + 6*S*z2/T - 3*S**2*z4/T**2))` | yes |
| A6 | sympy | 143 | `assert_zero("normalization K shift", dK_norm - (z0 + n0/Ptarget))` | yes |
| A7 | sympy | 144 | `assert_zero("compatibility shift from competing K surfaces", d_compat - (dK_norm - dK_one_pole))` | yes |
| A8 | sympy | 145 | `assert_zero("compatibility first variation from eliminated surface", d_compat - d_compat_direct)` | yes |
| A9 | sympy | 146 | `assert_zero("compatibility first variation", d_compat_direct - (n0/Ptarget - 6*S*z2/T + 3*S**2*z4/T**2))` | yes |
| A10 | sympy | 147 | `assert_zero("transported-target normalization K surface", K_norm_transport_p - (B0+Z0slot+eps*z0+D0target))` | yes |
| A11 | sympy | 148-151 | `assert_zero("transported-target compatibility surface", compat_transport_p - (D0target - 3*(S+eps*z2)**2/(T+eps*z4)))` | yes |
| A12 | sympy | 152-155 | `assert_zero("transported-target compatibility first variation", d_compat_transport - (-6*S*z2/T + 3*S**2*z4/T**2))` | yes |
| A13 | sympy | 156 | `assert_zero("fixed-target compatibility z0 cancellation", sp.diff(d_compat_direct, z0))` | partial (see below) |
| A14 | sympy | 157 | `assert_zero("transported-target compatibility z0 cancellation", sp.diff(d_compat_transport, z0))` | partial (see below) |
| A15 | sympy | 158 | `assert_nonzero("normalization K surface still carries z0", sp.diff(dK_norm, z0))` | yes |
| A16 | sympy | 159-162 | `assert_nonzero("mutated compatibility transport should fail", ...)` | yes |
| A17 | sympy | 163-166 | `assert_nonzero("mutated transported-target compatibility should fail", ...)` | yes |
| A18 | sympy | 175 | `assert_zero("constant-prefactor P2 factorization", P2_base - (N2 - N2_const)/D0)` | yes |
| A19 | sympy | 176 | `assert_zero("constant-prefactor P4 factorization", P4_base.subs(N2, N2_const) - (N4 - N4_const)/D0)` | yes |
| A20 | sympy | 183 | `assert_zero("Y20 overlap lane 20", lam20 - 1)` | yes |
| A21 | sympy | 184 | `assert_zero("Y20 overlap lane 21", lam21 - 1/2)` | yes |
| A22 | sympy | 185 | `assert_zero("Y20 overlap lane 22", lam22 + 1)` | yes |
| A23 | sympy | 195 | `assert_zero("weak-axisymmetric grouped trace", xbar - x0)` | yes |
| A24 | sympy | 196 | `assert_zero("weak-axisymmetric b=3a", bx - 3*ax)` | yes |
| A25 | sympy | 204 | `assert_zero("static Xi1 slope", Xi1_proj - (na/N0 + za/D0))` | yes |
| A26 | sympy | 211 | `assert_zero("u2 projected-Maxwell slope", u2_slope - (D0*z2a - D2*za)/D0**2)` | yes |

Note on A13/A14 marked "partial": both `d_compat_direct` and `d_compat_transport` are constructed by directly differentiating closed-form expressions in `S, T, n0, z2, z4` (`compat_direct_p`, `compat_transport_p`) that do not contain `z0` to begin with. The non-trivial K-elimination step is verified separately by A4, A7, A8 (which chain together to show that the K-eliminated surface equals `compat_direct_p`). Given that chain, A13/A14 are structurally guaranteed once A4/A7/A8 hold — but they are not standalone tautologies because they document the result the chain produces. Acceptable as kept.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- mathematica script: `(missing)` — no `.wl` companion to `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py`

**What's wrong:**
The manifest entry for unit 011 has `is_status_only_candidate: False` and is not a checkpoint, so the second-engine policy applies: both SymPy and Mathematica must independently derive the unit's claims. Only the SymPy script is present. The audit prompt explicitly names the mathematica path as `(missing)`. No `.wl` script in the scripts directory verifies any of the 26 SymPy assertions above. Sub-type: `missing_mathematica`.

**Why this matters:**
The whole point of the dual-engine policy is to catch engine-specific simplification bugs or hidden assumptions in a single algebra system (e.g. SymPy's `series` / `simplify` quietly collapsing a factor it shouldn't, or `gaunt` returning a value in a convention different from what the script expects). With only SymPy present, every assertion in A1-A26 — including the non-trivial Y20 lane ratios `(1, 1/2, -1)` from `sympy.physics.wigner.gaunt`, the K-elimination chain (A4, A7, A8), and the linearized constant-prefactor factorizations (A18, A19) — has zero independent verification.

**Required change:**
Create a new file `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl` that independently re-derives, from the same physical premises (the isotropic full-bundle definitions of `D0, D2, D4, P0, P2, P4`, the projected-Maxwell perturbations `D_n -> D_n - eps*z_n`, `N_n -> N_n + eps*n_n`, the one-pole defect, the two K-surfaces, the real-Y20 grouped pattern `(1, 1/2, -1)`, and the static lane ansatz), each of the following claims (M1-M11 in the directive). The Mathematica script must not be a line-by-line port of the SymPy choreography — it must derive the variations using Mathematica-native tooling (`Series`, `Simplify`, `Solve`, `ThreeJSymbol` or `ClebschGordan`) and write its results to `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.txt`. Each claim must be enforced by `If[FullSimplify[lhs - rhs] =!= 0, Print["FAIL ..."]; Exit[1]]` or equivalent so that a wrong result aborts with non-zero exit code.

**Verification:**
After Codex creates the file, the verifier runs `redteam exec-mathematica 011` and confirms: (i) the script exists, (ii) it contains executable `If[... Exit[1]]` checks for each of M1-M11 in the directive, (iii) it exits 0, (iv) the saved output transcript contains a final `STATUS: PASS` line.

## Independent-derivation check (Mathematica)

Not applicable — no Mathematica script exists. See F1.

## Engine cross-check

Not applicable — only one engine is present. See F1.

## Verdict justification

The SymPy script's 26 assertions hold up under attack. I verified by hand: the linear coefficients in `eps` of the perturbed bundle moments match the asserted closed forms (A1, A2, A3); the K-elimination chain A4 + A7 + A8 + A9 is internally consistent (both `K_one_pole_p` and `K_norm_p` share the `B0 + Z0slot + eps*z0` constant, so subtracting them removes `z0` and reproduces `compat_direct_p`); the constant-prefactor factorizations A18, A19 follow by direct substitution and cancellation; the weak-Y20 lane algebra A23, A24 with `(1, 1/2, -1)` yields `xbar = x0`, `ax = ea*x1/4`, `bx = 3*ea*x1/4`; the static Xi1 and u2 lane slopes A25, A26 are direct first derivatives. Symbol assumptions (`nonzero=True` on the bundle backbones, no assumptions on the perturbations) are appropriate; no `positive=True` is invoked aggressively. The negative-control assertions A15-A17 protect against trivial passes. Output transcript is newer than the script (`May 11 vs May 4`), so freshness is fine.

The one substantive finding is structural: the second engine is absent and this unit is not status-only, so the dual-engine requirement is unmet. Verdict: `findings`. No stop-cold flag — the missing engine does not invalidate downstream results; it just means independent confirmation is owed.
