---
unit_id: 249
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T13:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 249

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- §1 (sympy_audit.py:61-71): added distinct source symbols `S_full, S_res`, defined the covariance source `S_cov = S_full - S_res`, formed the two ledger RHSs `rhs_proj = -2*S_full`, `rhs_res = -2*S_res`, the subtracted RHS `rhs_sub = -2*S_full + 2*S_res`, and asserted `simplify(rhs_sub - (-2*S_cov)) == 0`. The original placeholder-linearity assert (line 59) is retained as documentation.
- §2 (sympy_audit.py:107-120): added the branch-gap identity `(Gplus2 - Gminus2) - 2*Gamma0*alpha_h == 0`, two concrete sign-flip substitutions (gap = +1 at alpha_h=1/2, = -1 at alpha_h=-1/2), and the both-positive / anti-aligned-negative checks (`Gminus2 < 0` at alpha_h=3/2). The pre-existing factoring asserts (incl. the `Gamma1->alpha_h*Gamma0` self-substitution at line 91) remain as documentation.

**Assessment:**
Non-tautological. In §1 the two ledger RHSs are genuinely distinct (`-2*S_full` vs `-2*S_res`); the load-bearing claim is that subtracting them yields exactly `-2*(S_full − S_res) = -2*S_cov`, which represents the covariance-source content of eq:hsub-eq. No X−X, no content-free placeholders, no `diff` w.r.t. an absent variable. The transcript confirms `subtracted RHS = -2*S_full + 2*S_res` matched against `-2*S_cov`. In §2 the substantive checks now exercise sign/magnitude/ordering (the gap is `2*Gamma0*alpha_h`, sign tracks `alpha_h`; anti-aligned branch goes negative once `|alpha_h|>1`) and can genuinely fail. The `Gamma1->alpha_h*Gamma0` self-substitution is no longer the load-bearing step — it survives only on the retained documentation asserts.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New independent `.wl` created at the named path (`mathematica/..._mathematica_audit.wl`), untracked/new (not a diff of the `.py`). Uses the project `expectZero`/`expectTrue`/`expectApprox` harness with `Exit[1]` on failure.

**Assessment:**
Genuinely independent decomposition, not a transliteration:
- M3 derives the peak Möbius inverse via `Solve[Rpksym == RpkExpr, ah]` with `RpkExpr` built from a `Hdot[s_] := Gamma0 + s Gamma1` function and `/. Gamma1 -> ah Gamma0` — derived, not restated; differs from the `.py`'s direct `Gplus/Gminus` form.
- M4 `FreeQ[Rint, etah]` is the legitimate scale-cancellation independence test: `etah` is genuinely present in `Hplus = etah(I0+I1)` and `Hminus = etah(I0-I1)` and divides out of the constructed ratio. This is NOT the differentiate-w.r.t.-absent-variable trap. The integrated inverse is also derived via `Solve`, not restated.
- M1 anchors the transfer-law source to the covariance `Scov = Sfull - Sres` directly.
- M5 recomputes the benchmark packet with `expectApprox`, including the independent final-load-ratio cross-check `RfinalBench` vs `RintReport`, plus the `0 < abar_h < alpha_pk < 1` ordering.
Benchmark literals match the `.py` and the auditor's packet exactly. All six M-checks print `PASS`; exit 0.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- sympy_audit.py:179 re-anchored `alpha_int_num` from `(ratio_integrated_report - 1)/(... + 1)` to `(ratio_final - 1)/(ratio_final + 1)`, deriving it from the independent final-load ratio `hint_aligned/hint_antialigned`.
- sympy_audit.py:214 tolerance relaxed `5e-13 -> 5e-9` to match the reported-data precision.

**Assessment:**
The near-self-confirming loop is broken. `alpha_int_num` is now a function of `ratio_final` (the independent final-load measurement, `20.58070146/5.00843357 = 4.10920923126...`), checked against the published `0.6085499908172678`; the transcript shows the derived `0.6085499909138585`, differing at ~1e-10, within the relaxed `5e-9`. The consistency tie (line 213, `ratio_final ≈ ratio_integrated_report`) is retained, and the `alpha_peak_num > alpha_int_num` ordering (line 216) still holds against the independently-anchored value.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`subtracted RHS = -2*S_full + 2*S_res` / `matched -2*S_cov = -2*S_full + 2*S_res` (F1 §1 covariance anchor), `aligned-minus-anti gap = 2*Gamma0*alpha_h` (F1 §2), `alpha_int(report) = 0.6085499909138585` (F3 final-load-derived), `All symbolic and numerical checks passed.`

**Mathematica:** exit=0. Notable lines:
`PASS: M1 projected-minus-resolved source`, `M4 eta_h cancels = True` / `PASS: M4 eta_h cancels`, `ah solved from Rpk = (-1 + Rpksym)/(1 + Rpksym)` / `PASS: M3 peak Mobius inverse`, `M5 final-load ratio diff = 1.26...*^-9` / `PASS`, `PASS: M5 asymmetry ordering`. All six M-checks PASS.

**Output freshness:** confirmed. Saved `.txt` outputs (both 2026-06-03 13:05:18) are newer than the `.py` (12:53:43) and `.wl` (12:54:41) scripts — regenerated post-fix.

## Material-change assessment

`material_change`: false. No published/forwarded value changed. F1 and F2 only add verification surface (new asserts, new engine); F3 re-anchors a sub-check to an independent quantity and relaxes a tolerance to the precision the reported data supports without altering any hardcoded value. The six-number Session-II packet and all carried-forward run outputs are byte-identical to the prior version and to the auditor's documented packet. Downstream units do not depend on any newly-introduced quantity.

## Side observations (non-blocking)

- `alpha_final_num` (line 180) and `alpha_int_num` (line 179) are now syntactically identical expressions, so the line-215 check `abs(alpha_final_num - alpha_int_num) < 1e-10` is trivially 0. This pre-existed and is not part of any finding; the directive explicitly only asked to confirm it still holds (it does). Harmless documentation, not a regression.

## Verdict justification

All three findings are resolved with no regressions. The diff contains exactly the F1/F3 `.py` edits and no collateral changes; the F2 `.wl` is a genuinely independent new engine. F1's covariance-source identity and §2 sign/ordering checks are non-tautological and can fail; the prior placeholder-linearity and factoring asserts survive only as documentation, and the `Gamma1->alpha_h*Gamma0` self-substitution is no longer load-bearing. F3 breaks the self-confirming loop by anchoring the integrated asymmetry to the independent final-load ratio. The `.wl` does not re-introduce the variable-independence trap — its `FreeQ[Rint, etah]` is a true scale-cancellation test (etah present then divides out). Both engines exit 0, outputs are fresh, and `material_change` is false.
