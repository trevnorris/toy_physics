---
unit_id: 187
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T02:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 187

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py:55-86`: a new block declares eight positive primitive ratios `rL,rC,rG,rU,rEta,rW,rM,rT`, builds the three invariant *ratios* directly from the §1 monomial forms — `Ctr_ratio = (rG*rC/rU)**(1+deltaU)*(rT/rU)**(1+chi)`, `Cnt_ratio = (rL**2*rM/(rEta*rW**2))*(rG**2*rL**2/(rU*rW))**E*(rT/rU)**(-F)`, `eps_ratio = rC**2/(rU*rEta)` — and asserts `expand_log(log(ratio)) - row.subs(log_subs) == 0` for all three.
- Mathematica `mathematica/...mathematica_audit.wl:38-59`: the tautological `ctrRatio = FullSimplify[Exp[row...]]` definitions were deleted and replaced by native primitive-ratio monomials (`ctrRatio = (rG*rC/rU)^(1+deltaStar)*(rT/rU)^(1+chiStar)`, etc.), with `expectZero["log C_tr ratio - row_tr", PowerExpand[Log[ctrRatio]] - (rowTr /. logSubs)]` for all three.

**Assessment:**
Correct and non-tautological. The decisive change: the ratio is now built from the monomial primitive factors, while the comparison `row` is the independently-posited linear expression. `log(ratio)` and `row.subs(log_subs)` therefore come from two separate code sources, so a mistyped exponent in any monomial factor (e.g. swapping `(1+chi)` for `(1+deltaU)`) would now make the residual nonzero and the assertion fail. I hand-checked the three log-expansions against the posited rows: `log(Cnt_ratio)` yields DL:`2(1+E)`, DM:`+1`, DEta:`−1`, DW:`−(2+E)`, DG:`2E`, DU:`(F−E)`, DT:`−F` — matching `row_nt` term-for-term; `log(Ctr_ratio)` and `log(eps_ratio)` likewise match `row_tr`/`row_eta`. This is the exact opposite of the prior `Log[Exp[row]] − row` (X−X), which had no fail mode. Monomial forms match the directive's M1/M2/M3 and the report's reconstruction of notes §1; the `pi^2/L^2` and `sigma` constants correctly cancel in the ratios and do not appear. Exec logs show the three new PASS lines in both engines.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The Mathematica ratio block (wl 38-59) is now written natively against the primitive ratios `rL,rC,rG,rU,rEta,rW,rM,rT` using `PowerExpand[Log[...]]` and a `logSubs` rule, rather than positing `Exp[row]`. The monomial expressions are coded directly in Wolfram syntax, not transcribed from the SymPy text.

**Assessment:**
Resolved as a consequence of an *independent* F1 implementation, exactly as the directive specified (F2 required no edit beyond F1 done independently per engine). Each engine now derives `log(ratio)` from the monomial structure by its own algebra (`PowerExpand` in WL, `expand_log(force=True)` in SymPy) and compares to its own posited row. A coefficient error copied into both row definitions would still be caught by each engine's monomial check independently, so the cross-check is now genuine for the load-bearing monomial→row step. The `.wl` no longer contains any `Exp[(1+deltaStar)*...]`-style row-echoing definition (confirmed in the diff: those lines are deleted). The `$Assumptions` positivity was expressed as `Element[{...},Reals] && rL>0 && ...` rather than the directive's `(rL|rC|...)>0` shorthand — a benign, equivalent transcription, not a defect.

**On the noted deviation (dropped ratio-after-solve checks):** Acceptable, directive-sanctioned, and removes no load-bearing coverage. The dropped checks were `(ctrRatio/.sol)−1`, `(cntRatio/.sol)−1`, `(epsEtaRatio/.sol)−1`. Under the *old* definitions (`ratio = Exp[row]`), `ratio/.sol == 1` was algebraically identical to `row/.sol == 0`, so they were already redundant with the retained `row_tr/nt/eta after solve` checks (wl 102-104, all PASS in the log). Under the *new* definitions the ratios live in `r`-variables while `sol` is in `d`-variables, so substituting `sol` into them would not even apply — the directive explicitly permitted dropping them rather than awkwardly rebinding. Solve consistency is fully retained via the three `row after solve = 0` PASS lines.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `log C_tr ratio - row_tr = 0`
- `log C_nt ratio - row_nt = 0`
- `log epsilon_eta ratio - row_eta = 0`
- `matrix row 1/2/3 - exact row_* = 0`, `Delta_eta/T/mu finite law = 0`, `row_* after solve = 0`

**Mathematica:** exit=0. Notable lines:
- `PASS: log C_tr ratio - row_tr`
- `PASS: log C_nt ratio - row_nt`
- `PASS: log epsilon_eta ratio - row_eta`
- `PASS: selected minor determinant`, `PASS: Delta_*/finite law`, `PASS: row_* after solve`

All checks pass; both engines exit 0. The three new monomial-built ratio checks are present and reference primitive ratios, not `Exp[row]`.

**Output freshness:** confirmed. `scripts/output/...sympy_audit.txt` (mtime 01:41:57) and `mathematica/output/...mathematica_audit.txt` (01:42:05) are both newer than their scripts (01:34:00 and 01:34:18). Exec logs (01:41:55 / 01:42:03) are post-fix. Outputs were regenerated after Codex's edits.

## Material-change assessment

`material_change`: false.

The edits add *new verification coverage* (the monomial→row derivation) and remove three *redundant* solve-consistency checks. No derived result changed: the rows, matrix, minor determinant (`1+chi`), and the three boxed orbit laws (Delta_eta/T/mu) are identical pre- and post-fix. Downstream units depending on Stage 187's carry-forward formulas see no change. No upstream-stale flag is warranted on numerical grounds.

## Side observations (non-blocking)

- The `.wl` `logSubs`/positivity declaration uses `Element[{...},Reals] && rL>0 && ...` instead of the directive's `(rL|rC|...)>0` form — functionally equivalent; noted only for completeness.
- Both engines' banner still reads "STAGE 170 — EXACT ORBIT-QUOTIENT CLOSURE" (a pre-existing mislabel for unit 187, present before this fix and unchanged by it). Not introduced by Codex; out of red-team scripts-math scope. Flagging only as an observation.

## Verdict justification

Both findings are resolved. The fix replaces the genuinely tautological `Log[Exp[row]] − row` (X−X) construction with a monomial-derived `log(ratio) − row` check that has a real fail mode — the ratio is built from the §1 monomial primitive factors independently of the posited rows, so an exponent transcription error would now be caught. The Mathematica side is coded natively rather than transliterated, restoring a genuine two-engine cross-check of the stage's distinctive finite-linearity claim (D1). The single noted deviation (dropping the three ratio-after-solve checks) was directive-sanctioned and removes only coverage already provided by the retained `row after solve` checks. Both scripts exit 0, all assertions pass, and outputs are freshly regenerated. Verdict: verified.
