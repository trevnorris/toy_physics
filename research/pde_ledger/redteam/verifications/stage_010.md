---
unit_id: 010
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 010

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl` (223 lines). The directive nominally pointed at `scripts/`, but `.redteam-config.yaml` (lines 9-10) places Mathematica audit files under `mathematica/`. Codex documented this in the `## Applied: F1` deviation block, and the runner's path glob is `mathematica/moving_throat_pde_stage{N}_*_mathematica_audit.wl`. The deviation is required by the runner layout, not a workaround.

The script verifies each of M1-M17 with `If[FullSimplify[expr] =!= 0, Print["FAIL: M<n>"]; Exit[1]]` (resp. `=== 0` for the M16/M17 negative-control mutations). Mathematica-native constructs used as required: `Series[..., {eps, 0, 1}]` + `Coefficient[..., eps, 1]` (lines 39-41, 80, 99, 112, 139, 202) for first-order expansion; `Solve` (lines 68, 88, 122) for K surfaces with uniqueness check via `Length[poleSolutions] - 1 === 0` (lines 69, 89); `ThreeJSymbol`-composed Gaunt at lines 145-156 rather than calling a Gaunt routine directly; bundle slot definitions written first as closed-form functions `slot0[e_], slot2[e_], slot4[e_]` (lines 31-37) and then differentiated, rather than building `P0p/P2p/P4p` intermediates as in the SymPy script. Variable naming differs throughout (`den0/num0/slot0/onePoleSurface/normSurface/transportSurface/gauntByThreeJ/lambda2/laneN/meanTrace/axisTrace/branchTrace/primZ/primN`).

**Assessment:**
Per-claim substance check, looking for shortcuts vs. genuine independent derivation:

- **M1 (dP0):** `slot0Linear` independently derived via `Coefficient[Normal[Series[slot0[eps], {eps, 0, 1}]], eps, 1]`, compared to closed form `n0/D0 + N0 z0/D0^2`. Independent of SymPy's `sp.diff`. Non-tautological.
- **M2 (dP2):** Same independent `Series` path; RHS matches my hand-derivation `dP2 = n2/D0 + N2 z0/D0^2 + 2 N0 z2/D0^2 - 2 D2 n0/D0^2 - 4 D2 N0 z0/D0^3` from the perturbation definition (the directive's illustrative form had `2 N2 z0/D0^2`; codex correctly used `N2 z0/D0^2` per the actual derivative — confirmed by both engines independently).
- **M3 (dP4):** Same independent `Series` derivation; closed form has 10 terms in `eps^0` shifts, all cross-validated against the SymPy `sp.diff` result via both engines hitting 0.
- **M4 (K_one_pole closed form):** Independent `Solve` of `(K - B0 - Z0slot - eps z0)(T + eps z4) == 3 (S + eps z2)^2`; uniqueness check `Length[poleSolutions] - 1 === 0`; compared to closed form.
- **M5 (dK_one_pole):** Series expansion of the M4 surface, compared to closed form. Non-tautological.
- **M6 (K_norm closed form):** Independent `Solve` of `(N0 + eps n0)/(K - B0 - Z0slot - eps z0) == Ptarget`; uniqueness check; compared to closed form.
- **M7 (dK_norm):** Series expansion of M6, compared to closed form.
- **M8 (compatibility surface identity):** `compatSurface = FullSimplify[normSurface - onePoleSurface]` compared to `compatDirect = (N0 + eps n0)/Ptarget - 3 (S + eps z2)^2/(T + eps z4)`. Genuine: M4 and M6 were independently derived via `Solve`, so the M8 identity is computed by subtracting two `Solve` results, not by assuming.
- **M9 (dcompat first variation):** Coefficient of Series of M8 expression, compared to closed form `n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2` (with the load-bearing `+` sign on the z4 term — opposite the M16 mutation).
- **M10 (K_norm_transport):** Independent `Solve` against transported target `(N0 + eps n0)/D0target`.
- **M11 (compat_transport z0 cancellation):** `transportCompat = transportSurface - onePoleSurface` compared to `D0target - 3(S+eps z2)^2/(T+eps z4)`. The RHS contains no z0, so the assertion implicitly verifies the cancellation (non-tautological because `transportSurface` and `onePoleSurface` each contain `eps z0` from M10/M4).
- **M12 (dcompat_transport):** Series coefficient of M11, compared to closed form `-6 S z2/T + 3 S^2 z4/T^2` (sign-matched to M9 modulo the n0/Ptarget term).
- **M13 (Y20 overlap lanes):** `lambda2[m]` defined via `ThreeJSymbol`-composed Gaunt; checks `lambda2[0] - 1, lambda2[1] - 1/2, lambda2[2] + 1`, plus the two same-sign vanishings `gauntByThreeJ[2,2,2,0,1,1] == 0` and `[2,2,2,0,2,2] == 0`. This is the requested independent re-implementation, not calling sympy's `gaunt`.
- **M14 (weak-axisymmetric trace decomposition):** Uses M13's `lambda2[m]` outputs to build lanes, computes `meanTrace, axisTrace, branchTrace` from the grouped recipe `(lane0 + 2 lane1 + 2 lane2)/5`, `(2 lane0 - lane1 - lane2)/10`, `(lane1 - lane2)/2`, and checks against the four target identities including `branchTrace - 3 axisTrace`.
- **M15 (primitive static Xi):** Independent substitution `N0sym -> P^2/Delta^2` and direct comparison to the claimed RHS.
- **M16/M17 (negative-control mutations):** Encoded with the correct polarity (`=== 0` triggers FAIL). Both produce nonzero residual `(6*S^2*z4)/T^2` per the output transcript line 17-18, confirming the unmutated assertions in M9/M12 are not vacuous.

No M is verified via a shortcut that assumes its own claim. The script genuinely re-derives each identity from the underlying defining equations using Mathematica-native primitives.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:58-83`, the two `if not {...}.issubset(...)` symbol-presence blocks are gone and replaced with two `assert_zero` calls. The dP2 closed form is `n2/D0 + N2 z0/D0^2 + 2 N0 z2/D0^2 - 2 D2 n0/D0^2 - 4 D2 N0 z0/D0^3`. The dP4 closed form is the 10-term expansion derived from `P4p = (D0p^2 N4p - 2 D0p (D2p N2p + D4p N0p) + 3 D2p^2 N0p)/D0p^3`.

**Assessment:**
I re-derived dP2 by quotient rule from `P2p = (D0p N2p - 2 D2p N0p)/D0p^2`: at eps=0 the derivative is `n2/D0 + z0 N2/D0^2 + 2 N0 z2/D0^2 - 2 D2 n0/D0^2 - 4 D2 N0 z0/D0^3`, matching codex's RHS exactly. The directive's illustrative form had `2 N2 z0/D0^2` (coefficient 2 on the first cross term), which is wrong; codex correctly trusted its own derivation per the directive's "if the hand-derivation disagrees, trust the hand-derivation" instruction. The dP4 closed form is cross-checked by the independent Mathematica `Series` expansion (M3) hitting 0 against the same RHS. The new asserts are non-tautological — `dP2` and `dP4` are computed by `sp.diff(P{2,4}p, eps).subs(eps, 0)` from the perturbation definition, so a sign-flipped or wrong-coefficient RHS would leave a nonzero residual.

Deviation noted in `## Applied: F2`: "Used expansion-derived dP2/dP4 targets because the illustrative forms in the directive disagreed with the perturbation definitions." This is the correct response to the directive's caveat.

Diff is contained to lines 58-82 of the .py file; no collateral edits.

## Exec log assessment

**SymPy:** stage_010_sympy.log is not present in `redteam/exec_logs/`; only `stage_010_diff.patch` is. However, the saved transcript at `scripts/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.txt` (mtime 2026-05-21 11:51:54 MDT, 23 minutes after Codex's apply timestamp 11:29:33) ends with `STATUS: PASS`, indicating the script ran to completion with exit 0. Treating user-reported "Both engines exit 0; orchestrator confirmed" as authoritative for the missing per-stage log file.

Notable lines from the output transcript:
- `STEP 08 PROJECTED MAXWELL PUSH MASTER AUDIT`
- `Checked bundle perturbation slots, z0 cancellation from compatibility, grouped signature, and primitive static dependencies.`
- `Target-transport z0 cancellation guards = PASS`
- `STATUS: PASS`

**Mathematica:** stage_010_mathematica.log is similarly not in `redteam/exec_logs/`. The saved transcript at `mathematica/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.txt` (mtime 2026-05-21 11:51:36 MDT) ends with `STATUS: PASS`. M1-M15 residuals all 0; M16/M17 mutation residuals both `(6*S^2*z4)/T^2` (correctly nonzero, the encoded `=== 0` polarity would have failed if these had cancelled).

Notable lines:
- `M1 residual = 0` ... `M15 residual = 0`
- `M13 residuals = {0, 0, 0, 0, 0}`
- `M14 residuals = {0, 0, 0, 0}`
- `M16 mutation residual = (6*S^2*z4)/T^2`
- `M17 mutation residual = (6*S^2*z4)/T^2`
- `STATUS: PASS`

**Output freshness:** Confirmed.
- Script mtimes: sympy .py = 11:28:21, mathematica .wl = 11:29:13.
- Output mtimes: sympy .txt = 11:51:54, mathematica .txt = 11:51:36.
Both outputs are ~22 minutes newer than the corresponding scripts.

## Material-change assessment

`material_change`: false.

The F2 fix tightens an under-verified assertion but does not change any derived closed form that downstream units depend on — the dP2/dP4 RHS values were already implicit in the `P2p, P4p` definitions; the prior symbol-presence stubs would have caught only catastrophic omissions, not the actual algebra. F1 adds a second-engine cross-check but does not alter any value. No symbolic identity changed; no constant was renumbered. Downstream units that consume dP2/dP4, the compatibility surfaces, the Y20 overlaps, or the Xi formula see the same values they would have seen before this iteration.

## Side observations (non-blocking)

- The per-stage exec logs `stage_010_sympy.log` and `stage_010_mathematica.log` are absent from `redteam/exec_logs/` even though sibling stages (e.g., 005, 006, 011) have theirs. The saved `output/*.txt` transcripts compensate, but the orchestrator's log-collection step seems to have skipped this stage. Not a blocker since user message confirms orchestrator already saw exit 0 from both engines.
- The Mathematica script uses `Clear[...]` followed by `$Assumptions = Element[{...}, Reals] && D0 != 0 && T != 0 && ...`. This is broader than the SymPy script's `nonzero=True` flags but is appropriate for `FullSimplify` and does not narrow the assertion's scope.
- The directive said M2's illustrative form was `2 N2 z0/D0^2`. Codex's expansion-derived `N2 z0/D0^2` (coefficient 1) is the correct one, confirmed by hand. The directive's parenthetical "Wait — let me redo this" was the prompt author already catching the slip; codex resolved it correctly.

## Verdict justification

Both findings are resolved with non-tautological, independently-derived checks. The .wl genuinely re-implements each of M1-M17 using Mathematica-native primitives (`Series`/`Coefficient`/`Solve`/`ThreeJSymbol`), variable naming and derivation order differ substantially from the SymPy script, and the M16/M17 negative-control mutations correctly leave a nonzero residual. The F2 fix replaces symbol-presence stubs with closed-form `assert_zero` checks whose RHS I confirmed by hand for dP2 and which Mathematica independently confirms via Series expansion. Both engines exit 0; outputs are fresh post-edit. No regressions visible in the diff (which is tightly scoped to the two specified locations). Verdict: `verified`.
