---
unit_id: 185
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-30T15:50:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
iteration: 2
checkpoint: true
---

# Verification — unit 185 (iteration 2)

CHECKPOINT stage: higher bar applied (both engines required; assertions must be
substantive and load-bearing; no rubber-stamping; exact paper alignment). F3 is a
non-blocking hygiene item (banner label) and is not counted in `findings_total`;
it is confirmed applied below. Iteration 1 closed F2 + F3 and returned F1
`needs_rework` because the observable coefficients `C_tr_star`/`A_tr_star` were not
load-bearing (the law residuals multiplied the coefficient onto a quantity already
proven zero, and the "coefficient form" anchors were byte-identical X−X). The
iteration-2 fix (route (b): reconstruct `Theta1`/`Xi1` independently from the
microscopic slippage drifts) is verified resolved here.

## Per-finding outcomes

### F1 — insufficient_verification (observable coefficients not exercised)

**Classification:** resolved

**What changed (iteration 2):**
- SymPy `scripts/...sympy_audit.py:199-213` — DELETED the two tautological
  `C_tr,* coefficient form` / `A_tr,* coefficient form` anchors (the iter-1
  X−X re-transcriptions). `C_tr_star`/`A_tr_star` remain typed literals
  (L197-198). `Theta1` (L201-204) and `Xi1` (L205-208) are now reconstructed
  INDEPENDENTLY from the slippage drifts:
  `chi1_indep = chi0s*Sigma_chi`, `deltaU1_indep = deltaUs*Sigma_delta`, then
  `Theta1 = -(chi0s(1+chi0s)deltaU1_indep + deltaUs(1+deltaUs)chi1_indep)/D`
  and `Xi1 = Sigma_Z + (2chi0s/(1+chi0s))Sigma_chi + E_star*Sigma_eps
  - [4 epsWs deltaUs/(11(1-epss)(1+deltaUs)^2)]Sigma_delta`. Neither references
  `C_tr_star`/`A_tr_star`. Four checks added/retained (L210-213):
  `Theta_1 independent slippage law` (`Theta1 - (-C_tr_star*Sigma_tr)`),
  `Xi_1 independent slippage law` (`Xi1 - (A_tr_star*Sigma_tr + Sigma_nt)`), plus
  the two `monomial law` checks against `Sigma_tr_compiled`/`Sigma_nt_compiled`.
- Mathematica `mathematica/...mathematica_audit.wl:188-210` — identical mirror:
  `chi1Indep = chi0s*sigmaChi`, `deltaU1Indep = deltaUs*sigmaDelta`, same closed
  forms for `theta1`/`xi1` (no `cTrStar`/`aTrStar` inside), same four `expectZero`
  checks. `cTrStar`/`aTrStar` typed; no coefficient-form anchors remain.
- `grep` confirms no `coefficient form` or `normalization` string survives in
  either script — the iter-1 tautological anchors are gone.

**Assessment (checkpoint bar — load-bearing reasoning, not executed):**

(a) INDEPENDENCE — confirmed. The new `Theta1`/`Xi1` are functions of the
microscopic slippage drifts (`Sigma_chi`, `Sigma_delta`, `Sigma_Z`, `Sigma_eps`)
and the parameters `chi0s`, `deltaUs`, `epsWs`, `epss`, `E_star` only. They do NOT
reference `C_tr_star`/`A_tr_star`. So the comparison is between an independently
built quantity and the coefficient-typed RHS — not the circular X−X of iter 1.

(b) TEETH — confirmed `(C_typed − C_true)*Sigma_tr` form. Because `Theta1` does not
contain `C_tr_star`, the residual `Theta1 - (-C_tr_star*Sigma_tr)
= Theta1_indep + C_tr_star*Sigma_tr`. The identity holds because, factoring the
numerator, `Theta1_indep = -chi0s*deltaUs*[(1+chi0s)Sigma_delta + (1+deltaUs)Sigma_chi]/D
= -chi0s*deltaUs*Sigma_tr/D = -C_true*Sigma_tr` (using `Sigma_tr` def at sympy:61).
Hence with a WRONG typed coefficient `C_bad`, the residual is
`(-C_true + C_bad)*Sigma_tr = (C_bad − C_true)*Sigma_tr`, which is NONZERO (see
(c)) and FAILS. Same for `Xi_1`: `Xi1_indep = A_true*Sigma_tr + Sigma_nt` (verified
by expanding `A_tr_star*Sigma_tr` against the `Sigma_chi`/`Sigma_delta` coefficients
and `Sigma_nt = Sigma_Z + E_star*Sigma_eps - F_star*Sigma_delta`, with the
`Sigma_delta` coefficient `2chi0s/(1+deltaUs) - F_star = -4epsWs deltaUs/(11(1-epss)(1+deltaUs)^2)`
matching `Xi1_indep`), so a wrong `A_tr_star` leaves `(A_bad − A_true)*Sigma_tr ≠ 0`.
The two `monomial law` checks additionally exercise the primitive-exponent compiler
(`Sigma_tr_compiled`/`Sigma_nt_compiled`, sympy:164/181 — built by differentiating
`exp` of raw primitives) against the now-independent `Theta1`/`Xi1`, so a wrong
compiler exponent OR a wrong coefficient fails there too.

(c) `Sigma_tr` NOT identically zero — confirmed. The log prints
`Sigma_tr = -(chi0s+1)*(kU - tau1) + (deltaUs+1)*(c1+gam1-kU)`, a nonzero symbolic
expression in the free drifts `tau1, kU, c1, gam1` and parameters. So the
`(C_bad − C_true)*Sigma_tr` fail-residual cannot collapse to 0 — the check has
genuine teeth.

(d) NO REGRESSION — confirmed. F2 (det minor) and F3 (banner) are byte-unchanged
from iter 1 per the diff (the `## Iteration 2` edits touched only the observable
block). Both engines still PASS the full prior check set; the printed
`tau1`/`kappa_eta`/`mu1` carry-forward forms are unchanged and engine-agree.

This fully closes F1's headline concern: a transcription slip in
`C_tr_star`/`A_tr_star` now fails an assertion rather than passing silently.
Resolved on the checkpoint bar.

### F2 — insufficient_verification (det M_*^(τ,κη,μ) = 1+χ0* unverified)

**Classification:** resolved (carried clean from iteration 1; unchanged in iter 2)

`Mstar_minor` (sympy:222-229; wl:217-222) is BUILT by differentiating the
independently-defined drifts `Sigma_tr/Sigma_nt/Sigma_eta` w.r.t. the dependent
triple `(tau1,keta,mu1)` — entries not hardcoded. Rows compute to
`[1+chi0s,0,0]`, `[-F_star,-1,1]`, `[0,-1,0]`; cofactor along row 1 gives
`(1+chi0s)·det([[-1,1],[-1,0]]) = (1+chi0s)·1 = 1+chi0s`. Residual
`det - (1+chi0s)` is identically 0 (passes) while the value `1+chi0s` is the
load-bearing nonzero positivity literal (`chi0s>0 ⇒ >1`). A wrong dependent-block
drift coefficient changes the determinant and fails. Both engines print
`det M_*^(tau,keta,mu) - (1+chi0s) = 0` + PASS. Anchors the appendix
`det M_*^{(τ,κη,μ)} = 1+χ0*` in-stage. Resolved.

### F3 — stale banner label (non-blocking hygiene)

**Classification:** resolved (cosmetic). Both banners read
`STAGE 185 — DIRECT MICROSCOPIC MONOMIALS` (sympy:41, wl:26); saved `.txt`
outputs echo "STAGE 185" (confirmed). No assertion depends on it.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Theta_1 independent slippage law = 0`, `Xi_1 independent slippage law = 0`
  (the new load-bearing checks)
- `Theta_1 monomial law = 0`, `Xi_1 monomial law = 0` (compiler-exercising checks)
- `det M_*^(tau,keta,mu) - (1+chi0s) = 0` (F2)
- All prior checks still `= 0` (`d ln C_tr,* - Sigma_tr`, primitive-compiler,
  substitution checks); banner `STAGE 185`.

**Mathematica:** exit=0. 36 PASS / 0 FAIL. Notable lines:
- `PASS: Theta_1 independent slippage law`, `PASS: Xi_1 independent slippage law`
- `PASS: Theta_1 monomial law`, `PASS: Xi_1 monomial law`
- `PASS: det M_*^(tau,keta,mu) - (1+chi0s)`
- All prior PASS lines retained; printed `tau1`/`kappa_eta`/`mu1` engine-agree
  with SymPy; banner `STAGE 185`.

**Output freshness:** confirmed. sympy script mtime 09:03:03 < output 15:33:52;
mathematica script 09:03:17 < output 15:34:07. Saved `.txt` outputs carry the new
`independent slippage law` lines, the monomial-law and det checks, and the
"STAGE 185" banner (grep-verified). Outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false. The iteration-2 edits are assertion-only. `Theta1`/`Xi1`
were re-routed through the slippage drifts but their simplified values are
identical to the paper's `-C_tr,*Sigma_tr` / `A_tr,*Sigma_tr + Sigma_nt` (the log
prints the same closed forms in both engines), and the printed
`tau1`/`kappa_eta`/`mu1`/`mu1-full` carry-forward forms are unchanged. No derived
result that a downstream unit depends on changed. Stage 186 (which consumes the
rank-3 / `det>0` fact, now asserted by F2) is re-confirmed, not altered. No
narrow re-audit of downstream units needed on math grounds.

## Side observations (non-blocking)

- The two `monomial law` checks (`Sigma_tr_compiled`/`Sigma_nt_compiled`) are now
  genuinely additive coverage on top of the independent-slippage-law checks: they
  exercise the primitive-exponent compiler against the independent `Theta1`/`Xi1`,
  so they fail on either a wrong compiler exponent or a wrong coefficient. Not a
  defect — extra teeth.
- `expect_zero` (sympy:30) runs `sp.simplify(sp.expand(expr))`; the new rational
  constructions simplify cleanly, no `ConditionalExpression` artifacts in either
  log; checks closed well under the 600s cap.

## Verdict justification

`verified`. On the checkpoint higher bar, F1 is now genuinely resolved: `Theta1`
and `Xi1` are reconstructed independently from the microscopic slippage drifts
(`Sigma_chi`, `Sigma_delta`, `Sigma_Z`, `Sigma_eps`) with no reference to
`C_tr_star`/`A_tr_star`, the iter-1 tautological coefficient-form anchors are
deleted, and the `Theta_1/Xi_1 independent slippage law` residuals take the form
`(C_typed − C_true)*Sigma_tr` with `Sigma_tr` symbolically nonzero — so a wrong
coefficient literal now FAILS rather than passing silently (the exact concern F1
raised). F2 (differentiated det minor = 1+chi0s) and F3 (banner) carry clean and
unchanged. Both engines exit 0, agree, and the saved outputs are fresh.
`material_change: false` — no carry-forward derived value changed; downstream units
are re-confirmed, not affected.
