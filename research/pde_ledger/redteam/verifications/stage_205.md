---
unit_id: 205
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 5
findings_total: 5
material_change: false
---

# Verification — unit 205

The original report raised 4 findings (F1 missing_mathematica, F2 missing_branch, F3 tautological_check, F4 insufficient_verification); the directive carries a 5th mechanical item (F5 banner mislabel). All five have `## Applied:` blocks. Verified each below.

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl` (345 lines) covering M1–M11. Exits 0 with 22 PASS, 0 FAIL (`stage_205_mathematica.log`).

**Assessment — independent-derivation (load-bearing, not rubber-stamped):** The `.wl` is a genuinely independent derivation, not a transliteration of the `.py`. Three corresponding sections confirm a different decomposition with native primitives:

- *Log-Hessian (M1/M2).* `.py:66-67` assembles the closed form `Hlog = Hchi/Phi0 - (gvec*gvec.T)/Phi0**2` and contracts it. `.wl:143-154` instead builds `chiLocal = phiBase + gradVec.vars + 1/2 vars.hessMat.vars` and computes the log-Hessian by *actual second differentiation*: `logHessianByD = Table[D[Log[chiLocal], vars[[i]], vars[[j]]] /. (vars->0), ...]`, with slope/curvature from `D[chiAlong,u]`/`D[chiAlong,{u,2}]`. Different route (differentiation vs. algebraic assembly) — exactly what the directive's anti-transliteration guard demanded.
- *Predictors (M3–M6).* `.py:88,104,121,137` hard-substitute the closed forms `2(1-Φ₀)/(Φ₁±√Δ)`. `.wl:110-126` defines `ordinaryBranch`/`logBranch` that `Solve[model==0, tau]` and use `selectByLimit` to pick the branch whose `Φ₂→0`/`L₁→0` limit matches the affine/log target. Solve-driven, not substituted.
- *Gap coefficient (M11).* `.py:224` divides by `eps**3`. `.wl:334-341` uses `Series[..., {eps,0,3}]` then `Coefficient[gapSeries, eps, k]` for k=0..3, additionally pinning eps⁰/¹/² to zero. The `Coefficient`-extraction route required by the directive.

Native primitives throughout (`Solve`, `Reduce`, `Series`, `Coefficient`, `D`, `Limit`, `FullSimplify`+`Assumptions`, `Table`, matrix `.`), distinct variable naming (`phiBase`, `p0ap`, `l0lp`, …). `expectZero` strips `ConditionalExpression[0,_]` per project idiom. Independent derivation confirmed.

### F2 — missing_branch

**Classification:** resolved

**What changed:**
SymPy: affine negative branch added at `:100-112` (`Phi1a_n = -Phi1a_abs`, predictor denominator `Phi1a_n - sqrt(Delta_aff_n)` = the `sgn(Φ₁)=−1` form), residual + `Φ₂→0` limit both asserted zero; log negative branch added at `:133-145` (`L0l_n = -L0l_abs`, denominator `L0l_n - sqrt(Delta_log_n)`), residual + `L₁→0` limit asserted. Positive-branch checks at `:97-98,130-131` retained. Mathematica `.wl` independently exercises both branches: `tauAffNeg` (`:185-202`, `-p1anAbs`) and `tauLogNeg` (`:224-241`, `-l0lnAbs`), Solve-selected.

**Assessment:** Both slope signs are genuinely exercised — the negative-branch symbols are strictly negative (negation of a `positive=True` symbol), and the `sgn`-correct minus-sign denominator is used, matching appendix lines 710/714. Output (`stage_205_sympy.log:97-98,118-119`) shows the negative-slope residual and limit checks passing. Both branches exercised: yes.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
SymPy `:151-179` replaced the two by-construction root residuals (old `res_tp_plus`/`res_tp_minus`, both identically 0) with: (a) a radicand/product sign-bridge `expect_zero` (`:154-157`); (b) `positive_product_cases` (sign-aligned `(apos,bpos)`,`(-apos,-bpos)`) asserted `Q.positive` and `negative_product_cases` (opposite-sign) asserted `Q.negative` of the radicand (`:159-175`); (c) one residual retained as secondary, non-load-bearing confirmation (`:177-179`). Mathematica M7 (`.wl:243-285`) uses `Reduce` to establish `Equivalent[radicandRegion, criterionRegion]`, `Equivalent[rootExistRegion, criterionRegion]`, and `negativeRegion === False`, with τ± residuals only as secondary checks inside the criterion region.

**Assessment — is the new criterion genuinely falsifiable?** Yes. The SymPy check would fail if the sign relation were wrong: with a wrong-signed radicand (e.g. `2(Φ₀−1)/Φ₂`), the positive-product cases would evaluate negative and `sp.ask(Q.positive(...))` would not return True, raising the AssertionError. The output confirms the substitution resolved to concrete sign-bearing values `[2*apos/bpos, 2*apos/bpos]` (positive product → real radicand) and `[-2*apos/bpos, -2*apos/bpos]` (negative product → no real root). The Mathematica side is even stronger — a genuine logical `Equivalent` of two `Reduce` regions plus an existential `Exists[tauRoot, ...]` region check, none of which is by-construction zero. The load-bearing check is now the reality equivalence `(1−Φ₀)Φ₂>0 ⟺ real root`, exactly the headline theorem; the surviving residual is explicitly downgraded to confirmation. Tautology now falsifiable: yes.

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy `:181-189`: introduced `quad_model_g = (Φ₀g−1) + Φ₁g·τ + ½Φ₂g·τ²` and an `expect_zero("tangency model from quadratic closure", quad_model_g.subs({Φ₀g:1, Φ₁g:0}) − ½Φ₂g·τ²)`. The tangency form is now derived from the quadratic model under `Φ₀→1,Φ₁→0`, not just printed. Mirrored in `.wl:287-291` (M8) with `phi0tan→1, phi1tan→0`. Symbol `Phi2g` unified across the print and the assertion.

**Assessment:** The form is now genuinely derived (substitution into the general model and subtracted against the claimed `½Φ₂τ²`), non-tautological. Passes in both engines.

### F5 — banner mislabel (mechanical)

**Classification:** resolved

**What changed:**
SymPy banners changed `STAGE 188` → `STAGE 205` at `:35,228` (diff lines 9-10, 110-111). New `.wl` banner labeled `STAGE 205` (`:71,343`). No computation altered. Transcript headers now read STAGE 205.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `:97 quadratic affine residual = 0` and `:99 quadratic affine residual (negative slope) = 0` (F2 both branches).
- `:124 turning-point radicand/product sign bridge = 0`, `:125-128` positive-product radicand `[2*apos/bpos, 2*apos/bpos]`, `:129-132` negative-product `[-2*apos/bpos, -2*apos/bpos]`, `:133 turning-point root residual on real-root branch = 0` (F3 falsifiable criterion).
- `:139 tangency model from quadratic closure = 0` (F4).
- `:159 STAGE 205 SYMPY AUDIT PASSED`, `# exit_code: 0`.

**Mathematica:** exit=0, 22 PASS / 0 FAIL. Notable lines:
- `M7 radicand criterion via Reduce = True`, `M7 no real root on negative product = True`, `M7 real-root existence region = True` (the load-bearing reality equivalences, not by-construction residuals).
- `M3/M4` and `M5/M6` each show positive- AND negative-slope PASS lines (both branches).
- `M11 gap coefficient eps^0..eps^3` all PASS (Coefficient route).
- `STAGE 205 MATHEMATICA AUDIT PASSED`, `# exit_code: 0`.

PASS count is 22 (both the exec log and the saved `.txt` agree). This is full coverage of manifest M1–M11: M1,M2 (2) + M3,M4 pos/neg (4) + M5,M6 pos/neg (4) + M7 five checks (5) + M8 (1) + M9,M10 (2) + M11 four coefficients (4) = 22. The orchestrator's stated expectation of 23 is an off-by-one in the expected count, not a missing check — every manifest item has its own pass/fail check. Not a defect.

**Output freshness:** confirmed. `.py` mtime 11:05:22, its output `.txt` 11:12:51; `.wl` mtime 11:09:41, its output `.txt` 11:12:51. Both outputs are newer than their scripts — regenerated post-fix.

## Material-change assessment

`material_change`: false.

All five edits either add new falsifiable assertions (F3 reality criterion, F4 tangency) around already-correct quantities, add coverage of the previously-untested negative-slope branch (F2), add a new independent engine (F1), or relabel a banner (F5). No derived/computed value that a downstream unit consumes was changed — the predictor formulas `τ_quad`, `τ_log,2`, the discriminants, the gap coefficient `(L₀²+3L₁)/(6L₀³)`, and the bridge identities are all unchanged. The F3/F4 changes added checks; they did not alter a computed value. No downstream unit is affected.

## Side observations (non-blocking)

- `.py:159` declares `apos, bpos = positive`; the F3 sign-bridge relies on SymPy assumption propagation (`Q.positive`/`Q.negative` on `±2*apos/bpos`), which the output confirms resolves correctly. Sound as written.
- The `.wl` `selectByLimit` branch-selection helper (`:48-69`) is a clean, genuinely Solve-driven mechanism and a notable strength of the independent derivation; no concern.

## Verdict justification

All four report findings (plus the mechanical F5) are resolved. The new Mathematica `.wl` is a genuinely independent derivation — log-Hessian by symbolic second differentiation, predictors by `Solve`+limit-based branch selection, gap coefficient by `Series`+`Coefficient` — not a transliteration of the SymPy choreography. The F3 tautology is replaced by a falsifiable reality criterion (`Reduce`-based `Equivalent` regions in Mathematica; sign-bearing `Q.positive`/`Q.negative` radicand checks in SymPy) that would fail under a wrong sign relation, with the old by-construction residuals downgraded to secondary confirmation. F2 now exercises both slope branches with the `sgn`-correct denominators. Both engines exit 0, outputs are fresh, 22 PASS lines fully cover manifest M1–M11. No material change to any downstream-consumed result. Verdict: verified.
