---
unit_id: 020
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

# Verification — unit 020

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
Codex created the new file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl` (155 lines). It defines its own `GauntIntegral[la,lb,lc,ma,mb,mc]` from `ThreeJSymbol` and a `gauntWeight` prefactor (lines 21-30), then independently derives and asserts M1-M5 inside a top-level `Module`. Each claim has an explicit `If[!TrueQ[FullSimplify[...] === 0], Print["FAIL: ..."]; Exit[1]]` guard (lines 49-50, 53-55, 57-59, 61-63, 65-67, 91-93, 98-100, 110-112, 114-116, 127-129, 131-133, 135-137, 146-149) and ends with `Print["STATUS: PASS"]; Exit[0]`.

**Assessment:**
The Mathematica file is not a transliteration of the SymPy script. It uses Mathematica idioms (`ClearAll`, `Module`, `SetAttributes`, `Listable`, `$Assumptions`, `ThreeJSymbol`, `Solve`, `Det`, `FullSimplify`), defines its own Gaunt integral from first principles (`Sqrt[(2L1+1)(2L2+1)(2L3+1)/(4 Pi)]` times two `ThreeJSymbol` factors) rather than calling `sympy.physics.wigner.gaunt`, uses CamelCase (`HEven`, `gateJacobian`, `branchSolve`, `branchRules`) where the Python uses snake_case (`H_even`), and reorders operations so the angular block is computed before any wall-slope symbols are declared. M1 covers the five Gaunt-ratio claims including the non-tautological m=0 lane (line 35: `lambda0 = FullSimplify[(-1)^0 GauntIntegral[2,2,2,0,0,0]/overlapBase]` — it does NOT short-circuit; it computes the numerator and denominator separately and lets `FullSimplify` reduce them). M2 verifies the 1/27 determinant via `Det` on the symbolic Jacobian. M3 checks both the uniqueness (`Length[branchSolve] - 1 === 0`) and the explicit closed forms. M4 verifies the three deficit identities by substitution. M5 verifies the Xi1 closed form. The saved output (`mathematica/output/.../stage020...txt`, mtime newer than the script) shows every residual = 0 and `STATUS: PASS`.

Spot-checked the Gaunt ratios against known values: gaunt(2,2,2,0,1,-1)/gaunt(2,2,2,0,0,0) = -1/2 (so `(-1)^1 * (-1/2) = 1/2` ✓), gaunt(2,2,2,0,2,-2)/gaunt(2,2,2,0,0,0) = -1 (so `(-1)^2 * (-1) = -1` ✓). The same-sign m=1,2 cross terms vanish identically (selection rule on m1+m2+m3). All consistent.

No collateral edits — only the new `.wl` file was added.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:20-26`: the `if m == 0: return sp.Integer(1)` short-circuit was removed and replaced with an `if m != 0:` guard wrapping the same-sign cross-term check. The return statement is now reached for all m in {0, 1, 2} and uniformly evaluates `(sp.Integer(-1) ** m) * gaunt(2,2,2,0,m,-m) / base`. The diff (`stage_020_diff.patch`) shows exactly the swap of the m=0 short-circuit for the m!=0 guard with no other edits.

**Assessment:**
The change matches the directive's "After" code block exactly. For m=0 the function now computes `gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0)` and `assert_zero('Y20 overlap lane 20', lam20 - 1)` will fail if SymPy's `gaunt(2,2,2,0,0,0)` ever returns 0 or an inconsistent value — the assertion is no longer a literal `1 - 1`. Same-sign cross-term guard is correctly gated on `m != 0` so it isn't triggered by the m=0 base self-overlap. Function signature, name, and call sites unchanged. The SymPy exec log shows `STATUS: PASS` and `exit_code: 0`, confirming the new m=0 lane still passes — i.e., the Gaunt evaluation is self-consistent.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `dKSigma_from_even_gates = B01 + 27*B41 + Z01 + 27*Z41`
- `D21_on_compensated_branch = -3*B41 - 3*Z41`
- `Check_D21_plus_D01_over_9 = 0` and `Check_D41_minus_2D21_over_3_minus_D01_over_27 = 0`
- `STATUS: PASS`

**Mathematica:** exit=0 (inferred from `STATUS: PASS` line in saved output; no `stage_020_mathematica.log` in `redteam/exec_logs/`, but the output txt is post-fix). Notable lines:
- `M1 lambda_0 residual = 0`, `M1 lambda_1 residual = 0`, `M1 lambda_2 residual = 0`
- `M1 same-sign m=1 residual = 0`, `M1 same-sign m=2 residual = 0`
- `M2 even-gate Jacobian determinant residual = 0`
- `M3 solution count residual = 0`, `M3 dKSigma residual = 0`, `M3 dMSigma residual = 0`
- `M4 D01 residual = 0`, `M4 D21 residual = 0`, `M4 D41 residual = 0`
- `M5 Xi1 residual = 0`
- `STATUS: PASS`

**Output freshness:** Confirmed.
- sympy script mtime: 1779392353, sympy output mtime: 1779397233 (output newer ✓)
- mathematica script mtime: 1779392466, mathematica output mtime: 1779397357 (output newer ✓)

## Material-change assessment

`material_change`: false.

The F1 edit is purely additive (new independent Mathematica derivation that reproduces existing SymPy results — no new claim, no changed numeric value). The F2 edit tightens an existing assertion's non-tautology guarantee but does not change any computed value: `lam20` still equals `1` post-fix (it just now arrives through a real Gaunt division rather than a literal short-circuit), and all downstream assertions on the compensated branch reference unchanged closed forms. No downstream unit depends on a different result.

## Side observations (non-blocking)

- The orchestrator did not save a `redteam/exec_logs/stage_020_mathematica.log` file, only the saved output under `mathematica/output/`. The mathematica exit code is therefore inferred from the `STATUS: PASS` terminator in the .txt rather than a top-of-file `# exit_code: 0` line. This does not block verification (the freshness and content of the output are direct evidence of a successful run), but the missing log is worth flagging for the orchestrator's bookkeeping.
- The `expectedRules` list and `solveResiduals` vector in the .wl are computed but only used for printing — the actual pass/fail guards re-run `FullSimplify[(dKSigma /. branchRules) - (B01 + Z01 + 27 (B41 + Z41))]` directly, which is correct and non-tautological. No change requested.

## Verdict justification

Both findings are fully resolved. F1's new `.wl` is an independent Mathematica derivation (not a SymPy port) that exercises `ThreeJSymbol`, `Det`, `Solve`, and `FullSimplify` to confirm the same five claim blocks M1-M5; the saved output shows every residual evaluates to 0 with `STATUS: PASS`. F2's `if m == 0: return sp.Integer(1)` short-circuit is gone, replaced by a uniform Gaunt-ratio computation that exercises `gaunt(2,2,2,0,0,0)` non-tautologically; the SymPy script still exits 0. No regressions in the diff, no collateral edits, and the material content of the unit (wall-slope closed forms, determinant 1/27, compensated identities, Xi1) is unchanged.
