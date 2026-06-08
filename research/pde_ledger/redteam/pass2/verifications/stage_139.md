---
unit_id: 139
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 139

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED FULL RE-AUTHOR)

**Classification:** resolved

**What changed:**
The `.wl` (`moving_throat_pde_stage139_..._mathematica_audit.wl`) was fully re-authored so its
load-bearing inputs are derived in-engine instead of imported as Stage-134 literals:

- The two hardcoded imports `piStar = SetPrecision[1.50882951349316, 30]` and
  `sQStar = SetPrecision[0.658075937605429, 30]` (old wl:31-32) are GONE.
- `Pi_*` is now derived by a cleared-denominator bracketed `FindRoot` (wl:46-50):
  `gThresholdResidual[z_] := 2*z*(2*z*Exp[z] + Pi) - gMinusRoot*(4*z^2 + Pi^2)*(Exp[z] - 1)`,
  solved on `{p, 1.4, 1.6}` at `WorkingPrecision -> 80`. The residual numerator carries
  `(Exp[p] - 1)` — the correct sign, NOT the §131-F3/144 `(1 - Exp[p])` trap.
- `S_q(Pi_*)` is now derived by direct source-moment quadrature (wl:53-63):
  `NIntegrate[Sigma*K_q, {x,0,1}]` with `Sigma = piStar*Exp[-piStar*x]/(1-Exp[-piStar])`,
  `K_q = Cosh[(Pi/2)(1-x)]/Cosh[Pi/2]` at `WorkingPrecision -> 70`. A second independent
  symbolic route (`Integrate` of the parametric moment, then substitute `p -> piStar`, wl:65-69)
  is cross-checked against the quadrature.
- The gain algebra (`mSNat`, `mQNat`, `mSComp`, `mQComp`) is unchanged in form and is now fed by
  the derived `piStar`/`sQStar` (wl:71-77).
- The old hardcoded Stage-134 kernel closed form that used to FEED the S_q check (old wl:90-94,
  `sQRecon = ...`) was DELETED. The kernel literals `1.50882951349316` and `0.658075937605429`
  now survive ONLY as expectApprox cross-check TARGETS (wl:98-99), not as sources.
- The "sanctioned mirror of the SymPy route" comment (old wl:42) is REMOVED; the new comment
  block reads as an independent route.

**Assessment:**
The edit fully matches all four required-change items and all acceptance criteria.

- Independence: confirmed. `Pi_*` comes from a root of `g(Pi)=g_-` (FindRoot), and `S_q` from a
  quadrature — neither reuses the Stage-134 hardcoded literal as its value. The two literals appear
  only as RHS targets of `expectApprox["Pi_* value", ...]` and `expectApprox["S_q(Pi_*) value", ...]`,
  i.e. corroboration, not source. `grep`-level: no `SetPrecision[1.50882951349316` and no
  `SetPrecision[0.658075937605429` survive as assignments.
- New assertions are non-tautological and load-bearing:
  - `Pi_* cleared-denominator residual` (residual `~1.6e-58`, PASS) is a genuine root self-check.
  - `Pi_* value` vs literal (residual `-4.44e-15`, PASS) — the FindRoot result corroborates Stage 134.
  - `S_q(Pi_*) value` vs literal (residual `2.22e-16`, PASS) — quadrature corroborates Stage 134.
  - `S_q quadrature vs direct Integrate` (residual `0`, PASS) — two distinct S_q routes agree.
- Preserved content: `g_-^F1 = 0.758035078944662826919680890414` branch-discrimination anchor
  (residual `0`, PASS) and `R_q^nat = 0.1454…` (≠ 1/4) are both retained and pass. All four gains
  (M_s/M_q nat/comp), r_F1, and `R_q^comp - 1/4` (1e-25) pass.
- No collateral edits beyond F1 in the stage-139 files; the SymPy `.py` was not touched.

## Exec log assessment

**SymPy:** exit=0. The `.py` is byte-identical to HEAD (`git diff HEAD` empty; it does not appear in
`git status`). Re-run still passes:
```
S_q(Pi_*) = 0.658075937605429000000000000000
g_-^F1 = 0.758035078944662826919680890414
all assertions passed
```

**Mathematica:** exit=0. All expectApprox lines PASS, including the new ones:
```
PASS: Pi_* cleared-denominator residual   (residual ~1.6e-58)
PASS: Pi_* value                          (residual -4.44e-15)
PASS: S_q(Pi_*) value                     (residual  2.22e-16)
PASS: S_q quadrature vs direct Integrate  (residual  0)
PASS: g_-^F1 value                        (residual  0)
PASS: M_s^nat,*  M_q^nat,*  M_s^comp,*  M_q^comp,*  R_q^nat  R_q^comp - 1/4
Stage 139 Mathematica audit passed.
```

**Output freshness:** confirmed. `.wl` mtime 2026-06-08 10:34:25; its `.txt` output 10:53:35 (newer
than script). The `.py` mtime is 2026-05-29 16:51 (HEAD, untouched) while its `.txt` was regenerated
2026-06-08 10:53:35. Both exec logs dated 2026-06-08T10:46, exit 0.

## Material-change assessment

`material_change`: false.

This is a method-only re-author of the Mathematica engine. Every downstream-relevant deliverable
value is preserved at the checked 1e-12 precision (all expectApprox literal targets unchanged and
passing). The only numerical movement is sub-15th-digit: the quadrature now reports
`S_q(Pi_*) = 0.658075937605429274…` vs the old hardcoded `…6580759376054290000…`, and the FindRoot
`Pi_*`/the propagated gains differ from the old printout only beyond the 15th significant figure
(M_s/M_q at WorkingPrecision 30 differ ~1e-15). These are exactly the acceptable high-precision-route
differences the directive/prompt sanctioned; no deliverable moved at 1e-12. No downstream unit that
consumes the boxed Stage-139 literals is affected.

## Side observations (non-blocking)

- The git working tree also contains modifications to stage 146 (`.wl`/`.py`/outputs) and pass-2
  scaffolding files. These are outside unit 139 and were not assessed; the stage-139 `.py` is
  confirmed clean (no diff vs HEAD).
- `tolAlg = 10^-25` is defined (wl:94) and used by the `R_q^comp - 1/4` check; `tolLit` is used by the
  rest. No unused-symbol concern.

## Verdict justification

`verified`. The single finding F1 is fully resolved: the `.wl` no longer mirrors the `.py` — it
derives `Pi_*` by a cleared-denominator bracketed FindRoot (correct `(Exp[p]-1)` sign, with a genuine
root residual self-check at `~1.6e-58`) and derives `S_q(Pi_*)` by independent source-moment
quadrature (cross-validated against a symbolic `Integrate` route at residual 0). The Stage-134
literals survive only as expectApprox corroboration targets, the "sanctioned mirror" comment is gone,
and the gain algebra is now fed by the in-engine-derived values. The SymPy `.py` is byte-identical to
HEAD and still passes. Both engines exit 0 with all assertions passing; outputs are fresh. No
deliverable value moved at the checked 1e-12 precision, so `material_change: false`.
