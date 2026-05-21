---
unit_id: 014
batch: I.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-21T15:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 014

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

`scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:97-104` — deleted the two trivial `z0 derivative map` / `n0 derivative map` assertions and inserted a comment block plus a double-loop that iterates over the six primed slips `(Qx, Sx, Hx, Dx, Px, Gx)` and the four non-trivial bundle derivatives `(z2d, z4d, n2d, n4d)`, asserting `sp.diff(_expr, _slip, 2) == 0` for each pair. This is 24 individual structural assertions.

**Assessment:**

The change matches the directive exactly. The new assertions are non-tautological: the bundle expressions z2/z4/n2/n4 are *literal* linear combinations of single primitive slips (e.g. `Delta*q1*S2`, `Q*s1`, `Gw*p1`, `P*g1`), each containing exactly one slip variable per term, so after `subs_der` and `/mu1` they become linear in the primed slips. If a future edit introduced a quadratic-in-slips term anywhere in the definitions on lines 50, 51, 53, 54 (e.g. accidentally writing `q1*s1` instead of `q1+s1`), the second-derivative assertion would fire. The exec log shows the script exited 0 after these checks were added. No collateral edits beyond what was asked.

### F2 — tautological_check

**Classification:** resolved

**What changed:**

`scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:106-125` — inserted the directive's verbatim comment block re-framing the next five `sp.diff` independence checks as construction-level. After the loop, added three bundle pull-back consistency assertions:

```
Xi_bundle = sp.simplify(2*Px/P - 2*Dx/Delta + Qx/(D0*Delta) - Q*Dx/(D0*Delta**2))
assert_zero("Xi bundle pull-back consistency", Xi - Xi_bundle)
K1_bundle = sp.simplify(-(z2d + z0d/sp.Integer(9)))
assert_zero("K1 bundle pull-back consistency", K1 - K1_bundle)
He_bundle = sp.simplify(-(z4d) + sp.Rational(2,3)*z2d - z0d/sp.Integer(27))
assert_zero("H_even bundle pull-back consistency", He - He_bundle)
```

**Assessment:**

Matches the directive verbatim. The Xi_bundle check is the weakest of the three (it just re-states the algebraic form Xi reduces to after substitution), but K1_bundle and He_bundle are real bundle-vs-primitive consistency checks: K1 was constructed by substituting subs_der into `-(z2+z0/9)` then dividing by mu1, while K1_bundle is constructed by combining the already-built z2d and z0d. These two paths differ in when sympy distributes simplify across addition, so the assertion catches any inconsistency between the bundle and primitive constructions of K1 (e.g. a typo that creeps into either the K1 line or the z2d/z0d lines). The construction-level caveat comment correctly flags that the preceding five `sp.diff(EXPR, SYM) == 0` checks are structural, not physical. Exec log: passes.

### F3 — tautological_check

**Classification:** resolved

**What changed:**

`scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:139-146` — deleted the `A + (-A) = 0` line that pretended to track the Z2 slot. In its place: (a) added the directive's verbatim comment block above the two solve-roundtrip lines marking them as tautological-by-sp.solve-contract, (b) defined `Z2_slot = (Q*S2 - Hport*Delta)/Delta**2` afresh from primitive symbols, and (c) added two new assertions linking the compensation denominators to `-9*Delta**4*Z2_slot` and `-9*Delta**3*Z2_slot`.

**Assessment:**

Matches the directive. The new check is not strictly stronger than the existing lines 131-132 (which already factored Hx_den and Sx_den as `9*Delta**2*(Delta*Hport - Q*S2)` and `9*Delta*(Delta*Hport - Q*S2)`), but it does provide an independent link via a fresh `Z2_slot` definition pinned to the upstream Z2 slot semantics. As the directive notes, the assertion would fire if a future edit changed the sign convention of Z2 (e.g. wrote `(Hport*Delta - Q*S2)/Delta**2` instead). This is meaningfully non-tautological for sign-convention drift, even if it is algebraically derivable from the existing lines 131-132. The deletion of the pure `A + (-A)` line removes the padding the auditor flagged. Exec log: passes.

### F4 — missing_verification_script

**Classification:** resolved

**What changed:**

Created `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl` (182 lines). The script defines the upstream primitive identities Z0, Z2, Z4, N0, N2, N4 from physical sources, then performs the Taylor lift via formal differentiation: `D[Z[Q + mu1*Qx*ell, ...], ell] /. ell -> 0 / mu1`. From the lifted bundle slips it derives `XiLoad = n0d/N0[P,Delta] + z0d/D0`, `K1 = -(z2d + z0d/9)`, `He = -z4d + (2/3)*z2d - z0d/27`, `deltaP2`, `deltaP4`. The 10 claim manifest items M1-M10 are then guarded by `expectZero` / `expectNonzero` / `expectEqual` helpers that `Exit[1]` on failure. A saved output transcript was produced at `mathematica/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.txt` (mtime 2026-05-21 15:01:14, newer than the .wl script mtime 2026-05-21 12:54:26).

**Assessment:**

The Mathematica script is a genuinely independent derivation rather than a transliteration: it defines the primitive identities (Z0, Z2, etc.) and uses `D[..., ell]` formal differentiation, structurally different from SymPy's `subs_der` substitution. Cross-engine agreement on the non-trivial results is confirmed by the saved output:

- M2 (`D[XiLoad, Px] == 2/P`): residual 0, PASS — matches SymPy's `coef_Xi_Px - 2/P == 0` (A5).
- M3 (`D[deltaP2, Gx] == -2*P/(D0*Delta**2)`): residual 0, PASS — matches SymPy's A6.
- M4 (`D[deltaP4, Gx] != 0`): residual `(2*D0*Delta*Gw + 4*D2*Delta*P - 4*D0*P*S2)/(D0^2*Delta^3)`, expand SymPy A7 transcript `2*(D0*Delta*Gw - 2*D0*P*S2 + 2*D2*Delta*P)/(D0**2*Delta**3)` and the two agree exactly.
- M5 (`Solve[{K1==0, He==0} /. {Sx->0, Hx->0}, {Qx, Deltax}] == {{Qx -> 0, Deltax -> 0}}`): PASS — matches SymPy `qd_only == [{Deltax: 0, Qx: 0}]`.
- M6 (analogous spectral solve): PASS — matches SymPy `sh_only`.
- M7, M8 (Jacobian determinants non-zero): both engines independently agree the determinants do not vanish.
- M9 (compensation denominators factor as `9*Delta^2*(Delta*Hport - Q*S2)` and `9*Delta*(Delta*Hport - Q*S2)`): residual 0 in both Hx and Sx denominators — matches SymPy A10, A11.
- M10 (sign-flip mutation residual `-18*Delta^2*Q*S2` is non-zero): PASS — matches SymPy A12.

The script is well-formed: it uses `Exit[1]` guards on every assertion via the `expectZero` / `expectNonzero` / `expectEqual` helpers, prints `STATUS: PASS` only on full success, and exits 0. No transliteration: the bundle slips are derived by formal `D[..., ell]` rather than re-substituting SymPy's pre-built z2/z4/n2/n4 expressions.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

- Line 6: `1) Mouth-local Taylor ansatz` — script reaches readout section.
- Lines 60-65: useful direct coefficients agree with assertion targets (`2/P`, `-2*P/(D0*Delta**2)`, non-zero deltaP4 form).
- Lines 71-75: mechanism sieve results `[{Deltax: 0, Qx: 0}]` and `[{Hx: 0, S2x: 0}]` confirm A13, A14 negative-existence.
- Line 106: `STATUS: PASS`, `# exit_code: 0` on line 107.

The new assertions (24 linearity checks + 3 bundle pull-back checks + 2 Z2-linked denominator checks) all run to completion without raising `AssertionError`, since the script reaches the final `STATUS: PASS`.

**Mathematica:** exit=0 (inferred from saved output transcript; orchestrator-captured `stage_014_mathematica.log` is NOT present in `redteam/exec_logs/`, but the saved transcript at `mathematica/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.txt` shows all 17 named guards as `PASS` and ends with `STATUS: PASS`). Notable lines from the saved transcript:

- Line 2-15: M1 sector decomposition — all seven sub-checks PASS with residual 0.
- Line 16-19: M2, M3 direct coefficients — residual 0, PASS.
- Line 20-21: M4 deltaP4 Gx dependence — non-zero residual exactly matching the SymPy A7 transcript form.
- Line 22-25: M5, M6 sieve solves — residual True (i.e. equality holds), PASS.
- Line 26-29: M7, M8 Jacobian determinants — non-zero residuals, PASS.
- Line 30-33: M9 denominator factorizations — residual 0, PASS.
- Line 34-35: M10 sign-flip mutation — residual `-18*Delta^2*Q*S2` (non-zero), PASS.
- Line 36: `STATUS: PASS`.

**Output freshness:**
- SymPy script: mtime 2026-05-21 12:53:05. SymPy output: mtime 2026-05-21 15:00:39. Output is newer than script. Fresh.
- Mathematica script: mtime 2026-05-21 12:54:26. Mathematica output: mtime 2026-05-21 15:01:14. Output is newer than script. Fresh.

## Material-change assessment

`material_change`: false.

The script's *derived results* (XiLoad, K1, He, deltaP2, deltaP4, comp_surface) are unchanged by the Codex edits. F1 replaced trivial assertions with structural assertions; F2 added bundle-pull-back consistency assertions; F3 replaced a tautology with a Z2-slot linkage assertion; F4 added a new independent Mathematica script that confirms (rather than alters) the SymPy results. No bottleneck definition, no coefficient, no comp_surface element changed. Downstream units depending on stage 014's derived results are not affected.

## Side observations (non-blocking)

1. The orchestrator did not capture `stage_014_mathematica.log` in `redteam/exec_logs/`. The saved output transcript at `mathematica/output/.../mathematica_audit.txt` IS fresh and shows `STATUS: PASS`, so we can confirm the script ran cleanly post-fix, but the missing exec log is a small gap relative to other stages (e.g. stage 001, 004-011 all have mathematica logs in `exec_logs/`). Not a blocker for verification since the saved output transcript serves the same role here.

2. The new linearity check loop (F1) uses Python loop-variable leakage: `_slip` and `_name` / `_expr` remain in scope after the loop. This is cosmetic and harmless; mentioned only for completeness.

3. The Mathematica's `XiLoad` derivation is reconstructed from `n0d/N0[P, Delta] + z0d/D0` rather than from a primitive-slip pull-back. Algebraically this matches SymPy's Xi (confirmed by M2). The construction-level M1 checks for XiLoad are structural in Mathematica too (XiLoad as built contains no Sx, Hx, Gx), as the directive itself acknowledged in its self-test rationale. Acceptable per the directive's design.

## Verdict justification

All four findings are resolved with the changes the directive asked for, applied verbatim with no deviations and no collateral edits. The SymPy script's new assertions are structurally non-tautological (24 second-derivative-zero linearity checks, 3 bundle-pull-back consistency checks, 2 Z2-slot denominator checks) — they would fire on plausible future typos. The new Mathematica script is genuinely independent (formal `D[..., ell]` Taylor lift starting from upstream primitive identities, not a transliteration) and confirms M2/M3/M4/M5/M6/M7/M8/M9/M10 in agreement with SymPy A5-A14. Both engines exited 0 with fresh saved-output transcripts. No regressions appear in the diff. Verdict: `verified`.
