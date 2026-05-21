---
unit_id: 013
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T15:10:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 013

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
A new file `mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl` (169 lines) was created. It implements `assertZero`/`assertNonzero` helpers (lines 8-22), executes the six-claim manifest in order:
- M1 one-sided Taylor first moment via `Integrate[Exp[-u] X, {u, 0, Infinity}]` and `Series` (lines 39-44)
- M2 W2 second moment and first-moment recovery (lines 46-54)
- M3 z-tower derived from `Zsource[x_] := (qFun[t] - hFun[t] x^2)/(dFun[t] - sFun[t] x^2 + x^4)` via `Series` and `timeDerivativeAtZero` (lines 66-83)
- M4 n-tower derived from `Nsource[x_] := (pFun[t] - gFun[t] x^2)^2/(dFun[t] - sFun[t] x^2 + x^4)^2` (lines 85-106)
- M5 `D[Xi, Px] - 2/P`, `D[deltaP2Der, Gx] + 2 P/(D0 Delta^2)`, and nonzero check on `D[deltaP4Der, Gx]` (lines 141-143)
- M6 2x2 source/denominator and spectral determinants nonzero, `Solve` returns trivial solutions only (lines 145-166)

**Assessment:**
The edit addresses the finding correctly. The chain-rule polynomials are genuinely re-derived from a primitive via Mathematica's `Series` rather than transliterated. Codex's documented deviation (using `(Q - Hport ell^2)/(Delta - S2 ell^2 + ell^4)` instead of the directive's denominator-only `Q/(Delta - S2 ell^2 - Hport ell^4)`) is legitimate: the directive's primitive places Hport at order ell^4 in the denominator, which cannot produce the `-h1` term that appears at order ell^2 in the expected z2. Codex's primitive places Hport at order ell^2 in the numerator, which correctly contributes `-h1*Delta^2/Delta^3` to the time-derivative at t=0. By hand-check, expanding `(Q - Hport ell^2)/(Delta - S2 ell^2 + ell^4)` to order ell^2 gives `Q/Delta + (Q*S2/Delta^2 - Hport/Delta)*ell^2 + ...`; differentiating the ell^2 coefficient w.r.t. t at t=0 yields exactly the expected z2 numerator `-Delta^2 h1 + Delta(Hport d1 + Q s1 + S2 q1) - 2 Q S2 d1`. The Mathematica output `OK M3 z*` with `residual = 0` confirms this match — `FullSimplify[derived - expected]` actually evaluated to 0, not a tautology. The same logic applies to the n-tower (whose primitive squares the bracket, producing the doubled-Hport-style contributions through `2*p1` and `2*S2*p1` terms in n2). M5's `delta P4 Gprime dependence` residual is `(2 D0 Delta Gw + 4 D2 Delta P - 4 D0 P S2)/(D0^2 Delta^3)` — a non-trivial closed form that a `FullSimplify`-aware tautology check would not produce. M6 determinants similarly evaluate to non-trivial closed forms `(-(Delta Hport) + Q S2)/Delta^4` and a higher-degree analog. No regression.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Lines 86-90 of the SymPy script — the two `for sym in (...)` loops with five `assert_zero` calls testing `sp.diff(EXPR, VAR)` against variables not in `EXPR` — were deleted. The diff (`stage_013_diff.patch`, lines 41-45) shows exactly the five lines removed, no other change. The substantive `assert_zero("dXi/dPprime", sp.diff(Xi, Px) - 2 / P)` remains at the new line 104 of the script.

**Assessment:**
Minimal, exact application of the directive. No collateral edits. The five tautologies (which could not fail) are gone. The remaining `dXi/dPprime`, `d(delta P2)/dGprime`, `delta P4 should depend on Gprime`, `qd_matrix.det()`, `sh_matrix.det()`, and `qd_only`/`sh_only` linear-solve assertions still exercise the K1/H_even/Xi forms through routes that can fail. The SymPy exec log confirms `exit_code: 0` and the final `STATUS: PASS` line still appears, so removing the tautologies did not destabilize the script.

### F3 — hardcoded_result

**Classification:** resolved

**What changed:**
Two comment blocks were inserted: (a) lines 49-57 immediately above `z0 = ...`, naming the master primitive `(Q(t) - Hport(t) ell^2)/(Delta(t) - S2(t) ell^2 + ell^4)` and citing the `.wl` file's claims M3 and M4; (b) lines 71-79 immediately above `n0 = ...`, naming the analogous primitive `(P(t) - Gw(t) ell^2)^2/(Delta(t) - S2(t) ell^2 + ell^4)^2` and citing the same `.wl` file's claim M4.

**Assessment:**
The anchor is concrete and points to a real file that actually does derive the polynomials (per F1 verification). The deviation Codex noted — that the cited primitive form differs from the directive's exemplar — is internally consistent: the SymPy comment names the same primitive the WL file uses, and the WL file reproduces the literal z/n polynomials with residual=0. The comments are not silent or generic; they name the variables, the order in `ell`, the leading-value-and-derivative split, and the destination of the verification. A future reader hitting an error in `z4` would land at the `.wl` file's M3 block, which is precisely the canonical place. The script's executable behavior is unchanged.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STEP 11 PROJECTED MAXWELL MOUTH-TAYLOR MASTER AUDIT`
- `Checked one-sided Taylor projection, bottleneck dependencies, G_W transport entry, and mechanism sieve.`
- `STATUS: PASS`

**Mathematica:** exit=0 (inferred — the `mathematica/output/.txt` file ends with `STAGE 013 MATHEMATICA AUDIT: PASS` produced only by the `Exit[0]` branch at the end of the script; no stage_013_mathematica.log was placed in `redteam/exec_logs/` but the saved transcript at `mathematica/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.txt` is the canonical output). Notable lines:
- `OK M3 z2 chain-rule coefficient residual = 0` (derived from Series, not transliterated)
- `OK M5 d(delta P2)/dGprime residual = 0`
- `OK M5 delta P4 Gprime dependence residual = (2*D0*Delta*Gw + 4*D2*Delta*P - 4*D0*P*S2)/(D0^2*Delta^3)` (non-trivial closed-form residual, exercises a real algebraic dependence)
- `OK M6 spectral determinant residual = (-(Delta*Hport) + Q*S2)/Delta^4`
- `STAGE 013 MATHEMATICA AUDIT: PASS`

**Output freshness:** Confirmed. The script mtimes are 2026-05-21 12:48 (`.py`) and 12:50 (`.wl`); the output transcripts are 15:00 (`.txt` for sympy) and 15:01 (`.txt` for mathematica). Both output files are newer than their corresponding scripts, so the saved transcripts reflect the post-fix scripts.

## Material-change assessment

`material_change`: false.

Justification: F2 removed five tautological assertions that could not change any derived quantity (they returned 0 by symbol-presence rules in SymPy). F3 added comments only. F1 created a new Mathematica file that is verification-only — it derives the same z/n polynomials the SymPy script uses and confirms they match, but does not export new values that downstream units would consume. No downstream-visible numerical result changed.

## Side observations (non-blocking)

1. The orchestrator did not write `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_013_mathematica.log`. The canonical Mathematica transcript at `mathematica/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.txt` is present, fresh, and clean, so the verification is unaffected, but the orchestrator's log-routing convention has a gap here that the user may want to check independently.
2. The SymPy script now has an unused `Ptarget` symbol declaration on line 42 (`D0, D2, D4, N0, Ptarget = sp.symbols(...)`) — `Ptarget` is never used downstream. This is pre-existing and not introduced by the F2 edit, but it is dead weight worth a future cleanup.
3. The directive's "M5 substantive coefficient checks" item lists `D[deltaP2_der, Gx] + 2*P/(D0*Delta^2) == 0` (matches SymPy line 109). The current WL implementation also includes that check (line 142) but does NOT add a counterpart to the SymPy script's existing `assert_nonzero("Xi should depend on Pprime", sp.diff(Xi, Px))` at line 144. The WL covers M5 fully via M5 dXi/dPprime / delta P2 / delta P4 G dependence; the missing duplication is intentional given the scope of F1 — no rework needed.

## Verdict justification

All three findings were applied with edits that match the directive's required changes. F1's new Mathematica audit genuinely derives the z- and n-towers via `Series` of a master primitive and produces non-trivial closed-form residuals on M5/M6 (so the checks could fail and didn't); the documented primitive deviation is correct and internally consistent with the polynomials the SymPy script uses. F2's five tautological assertion lines are cleanly excised, leaving the substantive `dXi/dPprime` line intact and the script still passing. F3's two comment blocks correctly anchor the hardcoded literals to the new Mathematica file's M3/M4 claims. Both exec outputs are PASS, both saved transcripts are fresher than their scripts, and no material-change cascade is triggered. Verdict: `verified`.
