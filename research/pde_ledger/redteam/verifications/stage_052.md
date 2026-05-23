---
unit_id: 052
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 052

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py:85-102`: the two old `expect_zero(zeta0_phys.subs(...) - zeta_req)` calls were replaced with a `sp.solve(KW*Omega0_sq_sym/Kphi0 - zeta_req, Omega0_sq_sym)` route (uniqueness asserted) compared against the pre-existing `Omega0_req_sq` definition, plus the analogous `sp.solve((KW*Omega0**2/Kphi0) - zeta_req, Kphi0)` route compared against `Kphi0_req`. The pre-existing `Omega0_req_sq` and `Kphi0_req` definitions (lines 77-78) are retained as the reference forms, exactly as the directive specifies.
- Mathematica `mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl:77-83`: the old `expectZero["threshold equality at fixed stiffness", (zetaPhys /. omega0^2 -> omega0ReqSq) - zetaReq]` and `... fixed overlap` lines were replaced with `omegaSqSol = First[omega0Sq /. Solve[(kW omega0Sq/kPhi0) - zetaReq == 0, omega0Sq]]` and `kPhi0Sol = First[kPhi0 /. Solve[(kW omega0^2/kPhi0) - zetaReq == 0, kPhi0]]` compared to `omega0ReqSq` and `kPhi0Req` respectively, with the required `solve(zeta_phys = zeta_req) for ...` labels.

**Assessment:**
The replacement is non-tautological: each `solve(...)` call derives the rescue variable independently from `zeta_phys - zeta_req == 0` and then is compared to the engine's pre-existing rearrangement `Omega0_req_sq = zeta_req*Kphi0/KW` / `Kphi0_req = KW*Omega0**2/zeta_req`. A wrong rearrangement on either side would now produce a nonzero residual. The `assert len(...) == 1` uniqueness guards are appropriate (SymPy returns a list, so this is a meaningful sanity check, not a tautology). The output transcripts confirm the new lines (`solve(zeta_phys = zeta_req) for Omega0^2 - expected = 0` and `... for Kphi0 - expected = 0`) and the absence of the prior `threshold equality at fixed stiffness/overlap` lines.

The Mathematica orchestrator-applied idiom patch (the `ConditionalExpression[e_, _] :> e` strip inside `expectZero` on lines 22-27) is mechanical and load-bearing only because Mathematica's `Solve` can wrap the dummy `omega0Sq` solution in a `ConditionalExpression` under aggressive `$Assumptions`; the strip preserves the spirit of `expectZero` (a residual that is identically zero under the declared domain). It does not weaken the substantive check — both residuals must still simplify to literal `0` after the strip-then-simplify pass, which the transcript confirms.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage052_nontwin_asymmetry_threshold_mathematica_audit.wl:39-44`: `zetaReq` is now produced by `Solve[(sReq - 1) - zetaSym (1 + eps (sReq - 2)) == 0, zetaSym]` with uniqueness check, instead of being typed in.
- Lines 51-59: `dZdPiAlt` is built via logarithmic differentiation `zetaReq (D[numZ, piTr]/numZ - D[denZ, piTr]/denZ)` of `Together[zetaReq]`'s numerator/denominator, and the new `expectZero["dZdPi vs dZdPiAlt (independent path)", dZdPi - dZdPiAlt]` cross-check runs *before* the comparison against the hand-typed `dZExpected`.
- Lines 64-68: `deltaZetaDerived = FullSimplify[Together[zetaReq - 1], ...]` replaces the direct `zetaReq - 1` subtraction; the subsequent `expectZero` uses `deltaZetaDerived - deltaExpected`.
- Lines 89-92: `softFracDerived = FullSimplify[Together[1 - 1/zetaReq], ...]` is added with a new `expectZero["softFrac vs Together[1 - 1/zetaReq] (independent path)", softFrac - softFracDerived]` before the comparison against the hand-typed `softExpected`.

**Assessment:**
The four insertions match the directive verbatim. The transcript shows all four new and preserved PASS lines:
- `PASS: dZdPi vs dZdPiAlt (independent path)` (output line 12)
- `PASS: dzeta_req/dPi - expected` (output line 15)
- `PASS: Delta_zeta - expected` (output line 18)
- `PASS: softFrac vs Together[1 - 1/zetaReq] (independent path)` (output line 27)
- `PASS: softening fraction - expected` (output line 34)

Mathematica now has independent algebraic paths for `zetaReq` (via `Solve`), `dZdPi` (via logarithmic differentiation of the rational form), `Delta_zeta` (via `Together`), and the softening fraction (via `Together[1 - 1/zetaReq]`). The hand-typed `dZExpected`, `deltaExpected`, and `softExpected` remain in place as the audit targets, but each is now compared to an engine-derived form. A typo in any of the three hand-typed forms would now fail the corresponding `expectZero` (because the engine-derived form would no longer match it). No collateral edits beyond what the directive requested.

## Exec log assessment

**SymPy:** exit=0 (per fix_batch_III.2.log; transcript ends cleanly with the STAGE 35 THEOREM LEDGER block). Notable lines:
- `solve(zeta_phys = zeta_req) for Omega0^2 - expected = 0` (output line 52)
- `solve(zeta_phys = zeta_req) for Kphi0 - expected = 0` (output line 53)
- `dzeta_req/dPi - expected = 0` (output line 22)
- `Delta_zeta - expected = 0` (output line 31)
- `softening fraction - expected = 0` (output line 74)
- The old `threshold equality at fixed stiffness/overlap` lines are absent.

**Mathematica:** exit=0 (transcript ends with `Stage 052 Mathematica audit passed.` on line 36, and the script's final `Exit[0]` is reached). Notable lines:
- `PASS: dZdPi vs dZdPiAlt (independent path)` (line 12)
- `PASS: dzeta_req/dPi - expected` (line 15)
- `PASS: Delta_zeta - expected` (line 18)
- `PASS: solve(zeta_phys = zeta_req) for Omega0^2 - expected` (line 23)
- `PASS: solve(zeta_phys = zeta_req) for Kphi0 - expected` (line 25)
- `PASS: softFrac vs Together[1 - 1/zetaReq] (independent path)` (line 27)
- `PASS: zeta_0^(twin) - 1` (line 29)
- `PASS: softening fraction - expected` (line 34)
- The old `threshold equality at fixed stiffness/overlap` lines are absent.

**Output freshness:** confirmed. Script mtimes: `*.py` = 17:19, `*.wl` = 17:22. Output mtimes: SymPy `.txt` = 17:35, Mathematica `.txt` = 17:24. Both outputs are newer than their corresponding scripts.

## Material-change assessment

`material_change`: false.

No derived result downstream of unit 052 changed value. The closed forms `zeta_req`, `dZ/dPi`, `Delta_zeta`, `Omega0_req_sq`, `Kphi0_req`, and `softFrac` are identical to their pre-edit values; only the verification path used to confirm them was strengthened (tautological substitutions replaced with `solve(...)` routes, hand-typed Mathematica expressions replaced with engine-derived equivalents). No downstream stage that depends on the numerical or symbolic value of any unit-052 output will see a different value.

## Side observations (non-blocking)

- The fresh symbol declaration `zetaSym = zetaSym;` on line 41 of the `.wl` is an idiomatic no-op (it assigns the symbol to itself, leaving it as a `Clear`-able free symbol). It is not strictly necessary because `Clear[piTr, cMix, eps, kW, kPhi0, omega0]` on line 34 does not touch `zetaSym`; but it is harmless and matches the directive's wording verbatim.
- The orchestrator-applied `ConditionalExpression[e_, _] :> e` strip in `expectZero` (lines 22-27) is a generic idiom patch and is not specific to the F1 change. It is a defensible defensive measure: `Solve` and `Reduce` commonly emit `ConditionalExpression[0, ...]` under nontrivial `$Assumptions`, and on the declared domain such a wrapper is identically zero. The strip is applied before the final `FullSimplify`, so a nonzero conditional value would not be silently passed (the post-strip `FullSimplify` and the `TrueQ[res === 0]` test would still fail).
- The Mathematica `softFrac` and `softFracDerived` are simplified independently, and the transcript shows `softFrac - softFracDerived = 0` — i.e., the engine confirms both paths agree before the hand-typed `softExpected` is compared. Good.
- The `kPhi0Sol` `Solve` in the `.wl` shadows the global `kPhi0` symbol inside the `Solve` call (since `kPhi0` is also the variable being solved for). Mathematica handles this correctly because `Solve[expr == 0, kPhi0]` localizes `kPhi0` automatically inside the solver; the result `First[kPhi0 /. Solve[...]]` is then a function of the still-global `omega0`, `kW`, `cMix`, `piTr`, `eps`. The transcript confirms the residual simplifies to zero. Non-blocking observation only.

## Verdict justification

Both findings are resolved with non-tautological, engine-independent checks; the diffs match the directive's "Required change" blocks verbatim with no collateral edits, both exec logs exit 0, all required new PASS lines are present, all pre-existing PASS lines are preserved, and the outputs are fresh relative to the scripts. The orchestrator's mechanical `ConditionalExpression` strip in `expectZero` is a benign, scoped idiom patch that does not weaken the substantive checks. No regressions, no downstream value changes — `material_change` is false.
