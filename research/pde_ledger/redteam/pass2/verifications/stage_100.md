---
unit_id: 100
batch: IV.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T23:50:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 100

## Per-finding outcomes

### F1 — mathematica_transliteration (independent-route required)

**Classification:** resolved

**What changed:**
`mathematica/...stage100...audit.wl:40-55`. Codex removed the transliterated
choreography (`yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5)`
followed by `ySeries = Expand[Normal[Series[yRet, {omega, 0, 5}]]]` and
`Coefficient[ySeries, omega, n]`) and replaced it with a hand-collected geometric
expansion. New body:
- `u2Coeff = 1/omegaQ^2`, `u5Coeff = I*chiQ*sigmaCan`, `uSquared4Coeff = u2Coeff^2` (wl:41-43)
- order coefficients formed by hand: `y2Coeff = u2Coeff/4`, `y4Coeff = uSquared4Coeff/4`,
  `y5Coeff = u5Coeff/4` (wl:48-50), with a comment that `1 + u + u^2` is the expansion and
  higher powers start at omega^6/omega^7.
- `k2 = k0*y2Coeff`, `k4 = k0*y4Coeff`, `gamma5 = (y5Coeff/I)*k0` (wl:53-55).
The diff (`exec_logs/stage_100_diff.patch`) touches ONLY this block; targets,
ratio checks, and closure construction (wl:57-81) are unchanged. The `.py`
reference engine was correctly left untouched (no diff for the `.py`). No
collateral edits.

**Assessment:**
The edit is correct and addresses the finding.

INDEPENDENCE — genuinely independent. The decisive structural difference: the
`.wl` no longer constructs the rational `yRet` and never calls `Series[]` on it.
The coefficients are obtained by the analytic geometric-series identity
`1/(1-u) = 1 + u + u^2 + ...` with the order-bookkeeping done by hand (only
`1+u+u^2` contributes through omega^5; `u` starts at omega^2). Compare:

- `.py:33-34`:
  `Y = 3/4 + (1/4)/(1 - omega**2/Omega**2 - I*chiQ*sigma_can*omega**5)`
  `Yser = sp.expand(sp.series(Y, omega, 0, 6).removeO())`  ← Series on full rational.
- `.wl:48-51`:
  `y2Coeff = u2Coeff/4; y4Coeff = uSquared4Coeff/4; y5Coeff = u5Coeff/4;`
  `ySeries = Expand[1 + y2Coeff*omega^2 + y4Coeff*omega^4 + y5Coeff*omega^5];`
  ← hand-built from the geometric coefficients; no `Series`, no `yRet`.
- coefficient read: `.py:36-38` does `K0 * Yser.coeff(omega, 2/4)` and
  `sp.im(Yser.coeff(omega,5))*K0`; `.wl:53-55` does `k0*y2Coeff`, `k0*y4Coeff`,
  `(y5Coeff/I)*k0` — read from the hand coefficients, not from a series-expanded
  rational. The `omega^4` term in particular is sourced from `u2Coeff^2`
  (the `u^2` term of the geometric series), which is a manifestly different
  derivation than letting `Series` collect it.

This is the geometric route explicitly sanctioned by the directive (one of the
named acceptable methods). `Series[yRet, ...]` no longer appears anywhere in the
file (acceptance criterion 1 met).

Deliverables (criterion 2): all four checks remain, can-fail and non-tautological
— each LHS K2/K4/Gamma5 is built from the geometric coefficients, so a wrong
factor would break the `==0`. Exec log shows M1 `K2/K2_target - NQ = 0` PASS,
M2 `K4/K4_target - NQ = 0` PASS, M3 `Gamma5/Gamma5_target - chiQ*NQ = 0` PASS,
M4 `closure_ratio - (mhat0^2 chi_Q N_Q - 1) = 0` PASS. Same values as the `.py`
reference. The closure is non-tautological for the same reason noted in the
audit: gamma5 feeds both M3 and M4.

chiQ (criterion 3): remains a free real symbol — `Element[chiQ, Reals]`
(wl:37), NOT in the positivity list (wl:38), never substituted to 1; comment at
wl:34 preserves the rationale. Met.

No value change (criterion 4): confirmed — see freshness/material-change below.

Exits 0 (criterion 5): both engines exit 0.

## Exec log assessment

**SymPy:** exit=0 (reference engine, untouched). Notable lines:
- `K2 = K0/(4*Omega**2)` / `K4 = K0/(4*Omega**4)` / `Gamma5 = 9*K0*chiQ/(32*Omega**5)`
- `K2/K2_target - NQ = 0`, `K4/K4_target - NQ = 0`, `Gamma5/Gamma5_target - chiQ*NQ = 0`
- `closure_ratio - (mhat0^2 chi_Q N_Q - 1) = 0`; `STAGE 100 AUDIT PASSED`.

**Mathematica:** exit=0. Notable lines:
- `Yhat_Q^ret series = 1 + (((9*I)/32)*chiQ*omega^5)/omegaQ^5 + omega^4/(4*omegaQ^4) + omega^2/(4*omegaQ^2)`
- `K2 = k0/(4*omegaQ^2)`, `K4 = k0/(4*omegaQ^4)`, `Gamma5 = (9*chiQ*k0)/(32*omegaQ^5)`
- `PASS: K2/K2_target - NQ`, `PASS: K4/K4_target - NQ`, `PASS: Gamma5/Gamma5_target - chiQ*NQ`
- `PASS: closure_ratio - (mhat0^2 chi_Q N_Q - 1)`; `Stage 100 Mathematica audit passed.`

These values are byte-identical (modulo symbol names) to the prior committed
output, including the hand-built `ySeries` reproducing the exact same series form.

**Output freshness:** confirmed. Both `.txt` outputs (mtime 2026-06-05 23:42:03)
are newer than the edited `.wl` (mtime 23:38:49). `git diff HEAD` shows ZERO diff
for both output `.txt` files (they match the committed HEAD byte-for-byte) and a
14-insert/5-delete diff confined to the `.wl` derivation block — i.e. the
derivation method changed but no emitted value did.

## Material-change assessment

`material_change`: false. Only the Mathematica derivation *method* changed; every
emitted/derived value (series form, K2, K4, Gamma5, NQ, all four PASS lines, the
closure residual form) is identical to HEAD. The `.py` reference was not touched.
No downstream unit can observe a different value, so no `upstream_stale`
propagation is warranted on numeric grounds.

## Side observations (non-blocking)

- The new route encodes the order-bookkeeping by human reasoning (the comment at
  wl:45-47 asserts that powers `u^3+` start at omega^6 and cannot affect K2/K4/Gamma5).
  This is correct for this specific `u` (lowest term omega^2, so `u^3` ~ omega^6),
  and it is the intended trade-off of the geometric route: independence is bought by
  hand-collection rather than machine truncation. Not a defect — flagging only that
  the truncation correctness now rests on the comment's reasoning, which I verified
  holds. Non-blocking.
- `ySeries` is now printed purely for display parity with the `.py`; the actual
  coefficients used in the checks are the hand-built `y{2,4,5}Coeff`, so the printed
  series is not on the assertion path. Harmless.

## Verdict justification

All five acceptance criteria are met: (1) no `Series[yRet, ...]` remains as the
coefficient source — replaced by a hand-collected geometric `1+u+u^2` expansion;
(2) all four deliverables (M1–M4) verified by genuine can-fail checks at the same
values; (3) chiQ stays a free real symbol; (4) no emitted value changed (output
`.txt` byte-identical to HEAD, only the `.wl` derivation block differs); (5) both
engines exit 0. The route is now auditably independent of the `.py`'s
Series-on-the-full-rational choreography. Verdict: verified.
