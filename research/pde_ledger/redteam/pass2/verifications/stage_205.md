---
unit_id: 205
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 205

The original report carries exactly one finding (F1 = `paper_misalignment`, card Verification line "Mathematica audit: none yet" though a passing `.wl` exists). It is a non-script, card-text-lag-class discrepancy USER-DEFERRED (directive `## Resolve before fix_loop`, `needs_user_resolution: true`); the directive explicitly directs Codex to apply nothing this batch. No source edit was directed and none was needed. Verification therefore confirms (A) the refreshed outputs are clean/fresh, (B) the audit disposition (independent `.wl`, non-tautological checks) holds on the refreshed artifacts, (C) the card-text-lag deferral to paper-cleanup (P4-51) is correctly classed, and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (card Verification line "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to paper-cleanup P4-51)

**What changed:**
Nothing in scripts. No diff patch exists (`redteam/pass2/exec_logs/stage_205_diff.patch` is absent), consistent with the directive's "Codex applies nothing on this unit." The finding routes `paper/stages/stage_205.tex:11` to the user/paper-cleanup tracker to either cite the present passing `.wl` (expected direction (a), with a MATHEMATICA_MIRROR_POLICY sync) or quarantine it (direction (b)). paper.tex is off-limits to the red-team (scripts-only), so this is correctly out-of-scope for any Codex/verifier source change.

**Assessment:**
Correctly classified and correctly deferred. This is the card-text-lag class: a pass-1 dual-engine retrofit added the `.wl` but the card's Verification pointer was not updated, so the card UNDER-reports coverage — it contradicts no math result. The `.wl` exists (9367 bytes) and passes M1-M11 cleanly (Mathematica output every line `PASS`, no `FAIL`, "STAGE 205 MATHEMATICA AUDIT PASSED", exit 0). Non-blocking; deferring the card-text correction to P4-51 leaves no math defect. No Codex edit was the right call.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed artifacts. Every load-bearing object is extracted by a method the `.py` does not contain:
  - log-Hessian / bridge — `.wl` DERIVES the log curvature by symbolic second differentiation of an explicit local Taylor model: `logHessianByD = Table[D[Log[chiLocal], vars[[i]], vars[[j]]] /. Thread[vars -> ConstantArray[0,5]], ...]` (wl:149-150), whereas the `.py` POSITS the closed form `Hlog = sp.simplify(Hchi/Phi0 - (gvec*gvec.T)/Phi0**2)` (py:66). Derive-by-`D[Log[...]]` vs posit-the-formula → independent.
  - quadratic predictors — `.wl` obtains both roots via `Solve` and selects the physical branch by zero-curvature limit: `tau /. Solve[ordinaryModel[...]==0, tau]` (wl:113) / `tau /. Solve[logModel[...]==0, tau]` (wl:122) feeding `selectByLimit[...]` (wl:116,125), whereas the `.py` POSITS the quadratic-formula radical `tau_quad = sp.simplify(2*(1-Phi0a)/(Phi1a+sp.sqrt(Delta_aff)))` (py:88) and `tau_log2 = sp.simplify(-2*sp.log(Phi0l)/(L0l+sp.sqrt(Delta_log)))` (py:121). Solve+branch-select vs posited radical → independent.
  - turning-point reality — `.wl` proves semialgebraic region equivalence via `Reduce`/`Exists`: `radicandRegion = Reduce[...]` (wl:246), `rootExistRegion = Reduce[Exists[tauRoot, ordinaryModel[...]==0]&&...]` (wl:261-262), whereas the `.py` substitutes case symbols and queries `sp.ask(sp.Q.positive(...))` / `sp.ask(sp.Q.negative(...))` (py:170,174). `ask`/case-substitution vs `Reduce`/`Exists` → independent (the `.wl` is the stronger check). The shared Section V `Series` step expands objects extracted differently (posited radical in py vs Solve-derived limit-selected root in wl), so it is not the "both Series-expand the same closed form" port pattern.
- **Non-tautological checks:** confirmed. Every predictor residual substitutes a candidate predictor (posited in py, Solve-derived in wl) back into a closure model defined separately, so a wrong predictor yields a nonzero residual. The trivial-case pre-checks are real can-fail statements: `Phi2->0` collapses `tau_quad` to `tau_aff` (py A4 / wl M4) and `L1->0` collapses `tau_log2` to `tau_log` (A6/M6), both confirmed zero, so the predictor checks are not vacuously satisfied. Both reality branches `(1-Φ0)Φ2>0` and `<0` are exercised in both engines (py L160-175 / wl M7).
- **0 result-value misalignments:** confirmed. The report reconciles 15/15 deliverable values MATCH (the F1 finding is a Verification-pointer staleness, not a result-value mismatch). The agreement-theorem cubic coefficient `(L0²+3L1)/(6L0³)` matches between engines and the notes — SymPy output L150-155 prints `eps³(L0e²+3L1e)/(6L0e³)`, Mathematica output L70-77 prints `M11 gap coefficient eps^0..eps^3 = 0` after subtracting it.

## Exec log assessment

**SymPy:** exit=0. Notable lines: "L1 - (Phi2/Phi0 - Phi1^2/Phi0^2) = 0" (L76); "Phi2 - Phi0*(L1 + L0^2) = 0" (L77); "quadratic affine residual = 0" / "...(negative slope) = 0" (L95,97); "turning-point root residual on real-root branch = 0" (L133); "quadratic predictors agree through O(eps^2) = 0" (L156); "STAGE 205 SYMPY AUDIT PASSED" (L159). All equalities hold to zero.

**Mathematica:** exit=0. Notable lines: "PASS: M1 log-curvature identity" (L15); "PASS: M3 affine residual positive slope" / "...negative slope" (L23,27); "M7 radicand criterion via Reduce = True" / "PASS" (L46-47); "M7 real-root existence region = True" / "PASS" (L50-51); "PASS: M11 gap coefficient eps^3" (L77); "STAGE 205 MATHEMATICA AUDIT PASSED" (L80). Every check prints PASS / `= 0` / `True`; no FAIL.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime 2026-06-09 16:51, newer than the `.py` (2026-06-02 11:05) and the `.wl` (2026-06-02 11:09). The SymPy output was already fresh (byte-equal expected) and re-run; the Mathematica output likewise refreshed. Both banners read "STAGE 205".

## Material-change assessment

`material_change`: false. No source code changed (no diff patch exists; the directive directed no Codex edit). The only artifacts touched were the regenerated committed `.txt` outputs from the orchestrator's independent re-run, which carry identical results. No derived result changed, so no downstream unit (>205) is affected.

## Side observations (non-blocking)

None beyond the single reported finding. The Section V / M9-M11 shared `Series` step was re-examined and correctly not flagged as a transliteration port: the object being expanded is extracted differently in each engine (posited radical vs Solve-derived, limit-selected root). I concur and add nothing.

## Verdict justification

`verified`. The sole finding (F1) is a non-script, card-text-lag `paper_misalignment` — the card's Verification line says "Mathematica audit: none yet" though a passing `.wl` exists — correctly classed and USER-DEFERRED to paper-cleanup (P4-51), with the directive directing no Codex edit; no source change was needed and none was made (no diff patch). The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (`D[Log[chiLocal],...]` log-Hessian, `Solve`+`selectByLimit` predictors, `Reduce`/`Exists` turning-point reality — all absent from the `.py`, which posits the corresponding closed forms), the checks are non-tautological (predictor substituted into an independently-defined closure model; trivial-limit collapses confirmed), and reconciliation is 15/15 MATCH with 0 result-value misalignments. Both engines pass at exit 0 with no FAIL; outputs are fresh (mtime 16:51 > scripts). No regressions; `material_change: false`.
