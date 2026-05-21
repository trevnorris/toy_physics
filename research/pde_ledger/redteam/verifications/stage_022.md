---
unit_id: 022
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 022

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**

In `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`:

- Section I (lines 33-48): the `Series[d0/dCons, ...]` route was replaced with an inverse-relation `Solve[coeffEqs, {u2Sym, u4Sym}]` on `Expand[yRespCand*dCons - d0]` truncated at omega^4. `u2, u4` are now extracted from the `Solve` solution.
- Section II (lines 50-69): the `Series[d0*nFac/dCons^2, ...]` route was replaced with `Solve[coeffEqsPref, {p0Sym, p2Sym, p4Sym}]` on `Expand[prefCand*dCons^2 - d0*nFac]`. The branch step now uses a direct `Expand[prefTrunc*y2Out]` (no `Series`) on the truncated polynomial `prefTrunc = p0 + p2 omega^2 + p4 omega^4`, exactly as the directive required.
- Section III (lines 89-92): the requested parallel-route comment was inserted before the existing six `expectZero` calls; no algebra changed (3x3 inverse-map has no independent route).
- Section IV (lines 105-116): the `Series[nProto, ...]` route was replaced with `Solve[coeffEqsN, {n0Cand, n2Cand, n4Cand}]` on `Expand[nCand*dProto^2 - pProto^2]`.
- Section V (lines 140-143): the hand-typed `j2 = ((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2`, `y2 = ...`, `h2 = FullSimplify[j2 + I y2, ...]` block was removed and replaced by `h2 = SphericalHankelH1[2, z]`. The subsequent `lambda2`, `lambda2Series`, `y2Resp`, `aStage4`, `bStage4`, `g5Stage4` derivation is unchanged but now consumes the built-in special function rather than a hand-typed polynomial-rational mirror of SymPy.

**Assessment:**

The change matches the directive exactly. The route by which each LHS reaches the `expectZero` check is now genuinely independent of SymPy's algebraic path: SymPy expands `D0/Dcons` as a Taylor series and reads coefficients; Mathematica posits a symbolic ansatz and solves a coefficient-by-coefficient linear system on the product. Likewise for Section II's `Pref` and Section IV's `Nproto`. Section V's `h2` now comes from Mathematica's built-in `SphericalHankelH1[2, z]` (which internally returns the closed-form `-((3 - 3 I z - z^2) E^(I z))/z^3`) rather than the explicit `j2 + i y2` mirror.

The new checks are non-tautological. `Solve[Expand[yRespCand*dCons - d0] coefficients == 0, {u2Sym, u4Sym}]` only returns `u2Sym -> -d2/d0, u4Sym -> (d2^2 - d0 d4)/d0^2` if the algebra is correct; the subsequent `expectZero["u2 formula", u2 + d2/d0]` is then a genuine consistency check between the inverse-relation solve and the claim target. If the target literal were wrong (e.g. wrong sign of `d2/d0`), the assertion would fail. The same non-tautology argument carries through Sections II and IV.

All target literals on the RHS (`d2/d0`, `(d2^2 - d0*d4)/d0^2`, `n0/d0`, `(d0*n2 - 2*d2*n0)/d0^2`, `(d0^2*n4 - 2*d0*(d2*n2 + d4*n0) + 3*d2^2*n0)/d0^3`, `p0proto^2/delta0^2`, `2*p0proto*(p0proto*s2 - delta0*gW)/delta0^3`, the `N4` long form, `a^2/(9*cS^2)`, `4*a^4/(81*cS^4)`, `a^5/(27*cS^5)`, `54*G*cS^5/(5*a^5*c^5)`, `6*G*cS^3/(5*a^3*c^5)`, `8*G*cS/(15*a*c^5)`) are preserved verbatim — matching directive verification point (d).

The Mathematica output (lines 9-10) shows the solved values `u2 = -(d2/d0)`, `u4 = (d2^2 - d0*d4)/d0^2` matching the targets, and every subsequent `expectZero` reports `0` and `PASS:`.

No collateral edits in this file beyond the requested ones. The `Clear[z, j2, y2, h2, ...]` defensive symbol-clear at line 135 still mentions `j2, y2` as cleared names, which is harmless (the symbols are no longer assigned, just cleared).

### F2 — tautological_check

**Classification:** resolved

**What changed:**

- SymPy (`scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`, was lines 245-261): the `Nproto_stage4 = ...`, `Nseries_stage4 = ...`, and the three `expect_zero("N0/N2/N4 round-trip into Stage-4 symbols", ...)` calls have been deleted. The dictionary printouts at lines 231-243 (`subbanner`, `dict_back = {...}`, `print/sp.pprint` lines) remain in place. The function now ends at the `return {"N0": N0, "N2": N2, "N4": N4}` line.
- Mathematica (was lines 104-109): the `nStage4 = Expand[Normal[Series[...]]]` construction and the three `expectZero["N0/N2/N4 round-trip", ...]` calls have been deleted. The `deltaBack`, `s2Back`, `p0Back` dictionary definitions at lines 126-128 remain (the directive explicitly allowed keeping them).

**Assessment:**

The deletions match the directive exactly. The remaining IV.1 prototype assertions (SymPy lines 217-222, Mathematica lines 119-124) still fire and still print residual `0`, confirming the non-tautological prototype-formula claims are preserved. The saved outputs contain `N0 prototype = 0`, `N2 prototype = 0`, `N4 prototype = 0` in both engines, and `grep round-trip` against both transcripts returns nothing — the redundant assertions are gone from the ledger.

Section V in both engines does not depend on `Nseries_stage4` / `nStage4` (it uses only the spherical-Hankel derivation), so the deletion is safe and introduces no broken references. I verified this by grepping the post-fix files: `Nseries_stage4`, `Nproto_stage4`, and `nStage4` are no longer referenced anywhere.

No collateral edits beyond the requested deletions.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

- Line 20-21: `u2 formula = 0`, `u4 formula = 0` — Section I assertions pass.
- Lines 40-42, 62-65: `P0/P2/P4/K0/K2/K4/Gamma5 formula = 0` — Section II assertions pass.
- Lines 96-98: `x20/x21/x22 recovered = 0` — Section III inverse-map passes.
- Lines 116-118: `N0/N2/N4 prototype = 0` — Section IV.1 assertions pass; the round-trip lines that previously followed are correctly absent.
- Lines 157-159: `Stage-4 A/B/G5 coefficient = 0` — Section V.1 fingerprint passes.
- Lines 199, 222-223: `mhat=1 target = 0`, `mhat=1 K2 target = 0`, `mhat=1 K4 target = 0` — Section V.2/V.3 invariant-product targets pass.
- Line 248: `# exit_code: 0`.

**Mathematica:** exit log file is not present at `redteam/exec_logs/stage_022_mathematica.log`. However, the saved output transcript at `mathematica/output/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.txt` records every `PASS:` for the 22 assertions (mtime 2026-05-21 15:02, newer than the script's 14:25 mtime), and the prompt confirms the script exited 0 on its most recent run. Notable lines from the output transcript:

- Lines 11-14: `u2 formula = 0` / `PASS: u2 formula`, `u4 formula = 0` / `PASS: u4 formula`.
- Lines 19-32: all seven Section II assertions `PASS:`.
- Lines 37-48: all six Section III assertions `PASS:`.
- Lines 53-58: all three Section IV.1 prototype assertions `PASS:` (no `round-trip` lines follow, as required).
- Lines 63-68: Section V Stage-4 `A/B/G5` coefficient assertions `PASS:`.
- Lines 75-80: Section VI `mhat=1 K0/K2/K4 target` assertions `PASS:`.

The print label at line 82 reads `FINAL STAGE-005 LEDGER:` — flagged as a stale copy-paste-style label from a prior stage naming. It is a print-statement string only and does not affect any assertion. Per the verifier instructions, this is noted but does not block verification.

**Output freshness:**

- `scripts/output/moving_throat_pde_stage022_grouped_p2_sympy_audit.txt` mtime 1779397235 vs script mtime 1779392714 — output is newer than script. Fresh.
- `mathematica/output/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.txt` mtime 1779397366 vs script mtime 1779392705 — output is newer than script. Fresh.

## Material-change assessment

`material_change`: false.

No derived result changed. The closed-form targets (`u2 = -d2/d0`, `u4 = (d2^2 - d0 d4)/d0^2`, `P0/P2/P4` rational expressions, `K0..K4, Gamma5` branch coefficients, `N0/N2/N4` prototype formulas, `A_stage4 = a^2/(9 c_s^2)`, `B_stage4 = 4 a^4/(81 c_s^4)`, `G5_stage4 = a^5/(27 c_s^5)`, `mhat^2 P0 = 54 G c_s^5/(5 a^5 c^5)`, `K2_target = 6 G c_s^3/(5 a^3 c^5)`, `K4_target = 8 G c_s/(15 a c^5)`) are all preserved verbatim on the RHS of the assertions. F1 changed only the route by which the LHS is computed (still produces the same simplified value, since `SphericalHankelH1[2, z]` is mathematically identical to `j_2 + i y_2`, and the inverse-relation solves recover the same coefficient expressions). F2 deleted a redundant check; it did not modify any active claim. Downstream units that consume stage 022's results (the invariant normalization product and the `K0/K2/K4` targets) see no change.

## Side observations (non-blocking)

- The `Clear[z, j2, y2, h2, lambda2, lambda2Series, y2, y2Static, y2Hat, y2HatOmega, ...]` symbol clear on Mathematica line 135 still lists `j2, y2` and lists `y2` twice. This is harmless defensive cleanup of names that are no longer assigned, but a future stylistic pass could prune the dead names. Not a finding.
- The Section VI `expectZero["mhat=1 K0 target", ...]` label differs from the SymPy script's label `mhat=1 target`. The RHS literal (`54*G*cS^5/(5*a^5*c^5)`) is unchanged and matches the SymPy claim, so the substance is identical; the label naming is a cosmetic deviation introduced in an earlier stage and is not within the scope of the F1/F2 directive.
- The print label at line 172 of the Mathematica file says `FINAL STAGE-005 LEDGER:` — the user-supplied note explains this is a codex copy-paste artifact in a print statement, not a script failure. Logged here as instructed; non-blocking.

## Verdict justification

Both findings are `resolved`. Codex applied the directive verbatim: Sections I, II, IV of the Mathematica script now use inverse-relation `Solve` on coefficient-equation lists rather than `Series` extraction; Section V now uses `SphericalHankelH1[2, z]` instead of the hand-typed `j_2 + i y_2`; and the tautological IV.2 round-trip block was removed from both engines. The exec log (SymPy exit 0) and the saved Mathematica transcript (script exit 0 confirmed by the prompt, 22 `PASS:` lines, no `round-trip` lines remaining) confirm every active assertion still reduces to residual `0`. Target literals on the RHS are preserved. No regressions appear in the diff. Verdict: `verified`.
