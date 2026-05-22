---
unit_id: 043
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 043

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The Mathematica audit (`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`) was reworked along all five sub-axes the directive prescribed:
- Line 50: `rPhi = FullSimplify[(y[[2]]/kappa1) / (y[[1]]/kappa0), ...]` (residue-ratio form, no longer `(y1/y0)/(kappa1/kappa0)`).
- Line 52: `dPhi = FullSimplify[Det[{{y0, y1}, {kappa0, kappa1}}], ...]` (no literal `kappa0 y1 - kappa1 y0`); the matching `dPhiExpected` on line 53 was flipped from negative to positive to honour `Det`'s row orientation (declared deviation in the `## Applied: F1` block).
- Lines 79-86: new endpoint checks `sUEndpointZero` (`deltaU -> 0`) and `sUEndpointInf` (`Limit[..., deltaU -> Infinity]`) with closed-form values `sigma/kU` and `(9/11) sigma/kU` respectively; both assert as `expectZero`.
- Line 135: `dPhiZ = FullSimplify[Det[{{y0, y1}, {z0, z1}}], ...]` (no literal `y0 z1 - y1 z0`).
- Lines 154-157: new `Series[mismatch, {deltaU, 0, 1}] // Normal` leading-coefficient check against `deltaU (rho0 - sigma0)/((1+rho0)(1+sigma0))`.

The Mathematica F3 block (lines 109-129) also replaces the SymPy-mirrored quotient comparison with a free-baseline `bBaseline` structural form plus mass-derivative independence checks (no longer a transliteration).

**Assessment:**
All five required structural replacements landed exactly where the directive specified. The literal expressions `kappa0 y1 - kappa1 y0` and `y0 z1 - y1 z0` are absent from the file (grep-confirmed by reading; only `Det[...]` forms appear). The Mathematica path is now genuinely distinct from the SymPy path at the algebra-construction level: SymPy expands `kappa0 y1 - kappa1 y0` directly, Mathematica computes a `Det`; SymPy substitutes `kappa0^2 -> sigma - kappa1^2` then `kappa1^2 -> (2/11) sigma`, Mathematica directly substitutes `(9/11) sigma`; SymPy uses `.subs(delta_U, 0)`, Mathematica uses `Limit[..., deltaU -> 0]`; mismatch is exercised by `Series` expansion in Mathematica with no SymPy counterpart. The directive's declared sign-flip deviation on `dPhiExpected` is internally consistent with the chosen `Det` row order (the output `D_phi = (deltaU*gS*gU*kappa0*kappa1)/(kU + deltaU*kU)` is `+`-sign, matching the positive expected, and the residual is `0`). Three new diagnostic lines (`overlap endpoint deltaU=0`, `overlap endpoint deltaU->Infinity`, `mismatch leading-in-deltaU coefficient`) are all printed with `PASS` in the Mathematica output.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy (`scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:103-115`): inserted the F2 block after the original `expect_zero("A_phi^(eff) - expected", ...)`. `Aphi_eff_min = Aphi_eff.subs(delta_U, 0)` is checked against `Kphi_eff - cUphi**2 * sigma/KU`; the ratio `(Kphi_eff - Aphi_eff)/(Kphi_eff - Aphi_eff_min)` is checked against `1 - (2/11) delta_U/(1+delta_U)`.
- Mathematica (`mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:90-98`): `aPhiEffMin = Limit[aPhiEff, deltaU -> 0]` is asserted against `kPhiEff - cUphi^2 sigma/kU`, and `overlapRatio = (kPhiEff - aPhiEff)/(kPhiEff - aPhiEffMin)` against the same `(2/11) deltaU/(1+deltaU)` factor.

**Assessment:**
Both new assertions are non-tautological. `Aphi_eff_min` is `Aphi_eff` with `delta_U -> 0` (or `Limit[..., deltaU -> 0]` in Mathematica, structurally distinct from SymPy's `subs`); a sign or factor-of-2 error in the `(2/11)` overlap weight would propagate into `Aphi_eff` but cancel out of `Aphi_eff_min` (since `Aphi_eff_min` collapses the `delta_U` dependence entirely), so the ratio is sensitive to `(2/11)` specifically. The SymPy output shows the ratio as `(9 delta_U + 11)/(11 (delta_U + 1))` which equals `1 - (2/11) delta_U/(1+delta_U)` — confirmed `= 0` residual on line 52 of the transcript. The Mathematica output line 39 shows the same ratio `(9 + 2/(1 + deltaU))/11` and residual `= 0` (line 40). The expected baseline `Aphi_eff_min = K_phi_eff - c_Uphi^2 sigma/KU` is set up from the physical `eps_phi = c_Uphi^2 sigma/(KU Kphi_eff)` identification independently of the `(2/11)` overlap structure, so it genuinely anchors the pole-shift identification. Original `A_phi^(eff) - expected` line is preserved (per directive); the new lines are additive.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy lines 121-148 replace the original `Msupp_cont_eval` vs `Msupp_expected` quotient comparison with four checks: `M_supp independent of mu_eta` (via `sp.diff`), `M_supp independent of mu_phi` (via `sp.diff`), `M_supp structural form (free baseline)` (with free symbol `B`), and `M_supp at baseline B = 8/pi^2` (the isolated baseline-value identification).
- Mathematica lines 102-129 mirror this with `D[mSuppCont, muEta]`, `D[mSuppCont, muPhi]`, `mSuppContInB` via `kappa0^2 -> bBaseline`, and the final `bBaseline -> 8/Pi^2` substitution. The structural check uses the new free baseline symbol.
- The original `expect_zero("M_supp - expected", ...)` / `expectZero["M_supp - expected", ...]` assertions are no longer present, as specified by the directive.

**Assessment:**
The mu-independence checks are genuine: `Msupp_cont` explicitly contains `mu_eta`, `mu_phi` in both numerator and denominator (visible at script lines 122-123 / wl 103-105), so `sp.diff(..., mu_eta) == 0` is a non-trivial cancellation test (a mistake in factoring `mu_eta` would produce a non-zero derivative). The free-baseline structural check is non-tautological because `B` (resp. `bBaseline`) is a free symbol — both sides must reduce to the same product structure for every value of `B`, which would catch a sign error in `(1 + ceU cUphi/(KU cB))^2` or in the `(1 - eps_eta)(1 - eps_phi_split)` denominator independently of any choice of baseline. The isolated `B = 8/pi^2` step finally records the hard-coded baseline value as a separate assertion. SymPy output lines 57-58, 71, 84 show all four checks passing; Mathematica output lines 46-55 show the same four `PASS` lines. Note: the residual SymPy output for `M_supp structural form` and `M_supp at baseline B = 8/pi^2` still prints a non-zero-looking factored expression (lines 60-70 and 73-83) — that's the pre-simplification *display* of the LHS, while the assertion compares LHS-RHS which reduces to `0` (printed on lines 71 and 84). This matches the script's intent (`sp.pprint(sp.factor(Msupp_cont_in_B))` then `expect_zero(..., Msupp_cont_in_B - Msupp_struct_expected)`).

## Exec log assessment

**SymPy:** exit=0 (transcript ends cleanly with the theorem ledger; no `FAIL`/`Traceback`/`AssertionError`). Notable lines:
- Line 47: `A_phi^(eff) at delta_U=0 (minimal) = 0`
- Line 52: `split-vs-minimal overlap ratio = 0`
- Lines 57-58: `M_supp independent of mu_eta = 0`, `M_supp independent of mu_phi = 0`
- Line 71: `M_supp structural form (free baseline) = 0`
- Line 84: `M_supp at baseline B = 8/pi^2 = 0`

**Mathematica:** exit=0 (transcript ends with `Stage 043 Mathematica audit passed.`; no `$Failed`/`FAIL`). Notable lines:
- Lines 27, 29: `PASS: overlap endpoint deltaU=0`, `PASS: overlap endpoint deltaU->Infinity`
- Line 38: `PASS: A_phi^(eff) at deltaU=0 (minimal)`
- Line 41: `PASS: split-vs-minimal overlap ratio`
- Lines 47, 49, 52, 55: all four `M_supp` checks `PASS`
- Line 71: `PASS: mismatch leading-in-deltaU coefficient`
- Two `Limit::alimv` warnings (lines 23, 35) are benign — Mathematica notes that assumptions on `deltaU` are ignored inside `Limit`, but the limits evaluate correctly and the residuals are `0`.

**Output freshness:** confirmed via `stat`:
- sympy script mtime 1779474987 < sympy output mtime 1779475219 (output newer by ~4 min)
- mathematica script mtime 1779475148 < mathematica output mtime 1779475229 (output newer by ~1.3 min)

Both outputs were re-generated after the Codex edits.

## Material-change assessment

`material_change`: false.

The closed-form results that downstream units could consume — `R_phi`, `D_phi`, `v.D_U.v`, `A_phi^(eff)`, `D_(phi z)`, the mismatch formula, and the `M_supp = 8/(pi^2) cB^2 (1+sigma0)^2/[Keta_eff Kphi_eff (1-eps_eta)(1-eps_phi^split)]` baseline — are all unchanged. F1's edits replaced one algebraic path (literal expansion) with another (`Det`) while keeping the *values* of `D_phi`, `dPhiZ`, etc. identical (only the orientation/sign convention on the Mathematica side was reflected in `dPhiExpected`, with no propagation outside this unit). F2 added new assertions without modifying any earlier expression. F3 split one assertion into four without changing the underlying `Msupp_cont` quotient or the final `8/pi^2` baseline. No downstream unit can see a different value than before.

## Side observations (non-blocking)

- The `Limit::alimv` warnings in the Mathematica transcript (lines 23, 35) are emitted whenever a `Limit[..., deltaU -> ...]` is invoked inside the global `$Assumptions` context (which constrains `deltaU > 0`). They're cosmetic; the limit values evaluated by Mathematica are correct. No action required, but a future cleanup could wrap the `Limit` calls in `Quiet[..., Limit::alimv]` if the warnings get annoying.
- The new `$Assumptions = $Assumptions && bBaseline > 0;` line in the Mathematica F3 block (line 115) mutates the global assumption context. This is local to stage 043 (the script ends before any cross-unit reuse) so it is benign here, but worth noting as a pattern.
- SymPy still uses `subs({kappa0**2: sigma - kappa1**2, kappa1**2: ...})` (line 89). This was not part of any finding and is correctly left alone.

## Verdict justification

All three findings are resolved with the directive's required edits applied verbatim (modulo the explicitly declared `dPhiExpected` sign flip needed to honour the new `Det` row convention, which is internally consistent). Both engines exit cleanly; the SymPy transcript ends with the stage-26 theorem ledger and no `AssertionError`, the Mathematica transcript ends with `Stage 043 Mathematica audit passed.` and no `$Failed`. Every new diagnostic line called out in the directive's verification section (`overlap endpoint deltaU=0`, `overlap endpoint deltaU->Infinity`, `A_phi^(eff) at deltaU=0 (minimal)`, `split-vs-minimal overlap ratio`, `M_supp independent of mu_eta`/`mu_phi`, `M_supp structural form (free baseline)`, `M_supp at baseline B = 8/Pi^2`, `mismatch leading-in-deltaU coefficient`) appears with a `= 0` residual in SymPy and a `PASS:` line in Mathematica. The Mathematica script no longer transliterates the SymPy algebra — it now uses `Det`, `Limit`, endpoint values, free-baseline structural reduction, and `Series` expansion as independent checks. No regressions visible in the diff or transcripts.
