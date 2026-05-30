---
unit_id: 147
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T22:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 5
findings_total: 5
material_change: false
---

# Verification — unit 147

Source of findings: Codex review `redteam/codex_reviews/stage_147.md` (5 findings R1-R5),
refined by the batch-6 Claude+Codex consult (`_consult_batch6.md` Q4/Q5/Q6 — added the
source-centering assertion and the ratio-label fix). Verified against the live `.py`/`.wl`,
the exec logs, the diff patch, and (per consult Q4 authorization) the owning prose
`paper/appendices/stage_appendix_part04.tex:798-848`.

## Per-finding outcomes

### R1 — tautological_check (SymPy A_T chain route rebuilt the same algebra)

**Classification:** resolved

**What changed:** `scripts/...sympy_audit.py:66-80`. The old hand-typed
`dSigma_dPi_at_star = 1/(1-S_star/4) + Pi_star*Sp_star/(4*(1-S_star/4)**2)` chain factor is
gone. Live code builds `Tm_of_Pi = sp.sqrt(sp.Rational(9,20)*(Pi/(1-Sformula/4)))` and lets
SymPy autodifferentiate: `dTm_dPi = sp.diff(Tm_of_Pi, Pi)`, `dg_dPi = sp.diff(gPi, Pi)`,
then `AT_autodiff = -(dTm_dPi.subs(Pi,Pi_star))/(dg_dPi.subs(Pi,Pi_star))`, asserted equal to
the closed-form `AT` to 1e-20.

**Assessment:** Genuinely independent of the target's primitive. The closed-form `AT`
(py:33-38) is a hand-expanded chain rule with literal factors `1/(1-S/4)`,
`Pi S'/(4(1-S/4)^2)`, `9/(40 T_*)`; `sp.diff` regenerates those factors symbolically from
`Sformula` and the sqrt structure — it never types them. A sign/power error localized to the
hand-written `AT` would NOT be reproduced by autodiff, so the assertion can fail. The only
shared inputs (`gPi`, `Sformula`, `Pi_star`) are upstream primitives audited elsewhere.
Matches the appendix retuning identity `A_T=-(dT_m/dPi)/(dg/dPi)` (appendix:827-831). Output:
`PASS: A_T closed form agrees with autodiff of T_m(Pi) (residual < 1e-20)`. Non-tautological.

### R2 — tautological_check (SymPy centered kernel checked against its own typed form) + CONSULT-ADDED source-centering

**Classification:** resolved

**What changed:** `scripts/...sympy_audit.py:93-128`. The old `Wcenter_const` extract-and-
compare-to-typed-constant tautology is replaced by TWO checks:
(a) Projection identity — for non-canonical normalized deformation `smallsigma=2x`:
`lhs_proj = ∫ W_*·(smallsigma - Sigma_*) dx` (full-kernel quadrature) compared to
`rhs_moment = A_T(ḡ_s - g_*) + B_T(S̄_s - S_*)` (algebraic moment formula), with a
normalization sanity guard, asserted to 1e-22.
(b) CONSULT Q5 source-centering — `center_resid = ∫ Sigma_*·W_* dx` asserted `< 1e-22`.

**Assessment:** Both routes use a different primitive than the target. The projection LHS
integrates the full typed kernel `W_*` numerically; the RHS assembles A_T/B_T against
moment integrals of the test profile only (it never integrates `W_*`), so a wrong sign,
dropped term, or c↔Kq swap breaks the equality. The consult correctly identified that the
projection alone is BLIND to the centering constants `-A_T g_*, -B_T S_*` (they cancel against
`(smallsigma - Sigma_*)` since both integrate to 1). The source-centering assertion
`∫ Sigma_* W_* dx = 0` genuinely tests those constants: dropping them leaves
`∫ Sigma_*(A_T c + B_T Kq) = A_T g_* + B_T S_* ≠ 0`, which would fail. This is the check that
actually exercises the `-A_T g_* - B_T S_*` constants. Both present and passing (see log).
Non-tautological. Appendix defs (Σ_*, ḡ/S̄, c/Kq, normalization) faithfully implemented.

### R3 — insufficient_verification (Mathematica kernel checked only at x=1/2) + CONSULT-ADDED source-centering

**Classification:** resolved

**What changed:** `mathematica/...mathematica_audit.wl:91-129`. The old single-sample
`(wCenter - (aT*c+bT*kq)) /. x -> 1/2` is replaced by THREE checks:
(a) Full-symbolic x-independence: `expectZero[..., Chop[FullSimplify[D[wStar-(aT*c+bT*kq), x]], 10^-25]]`
— a symbolic derivative over all x, not a sample.
(b) NIntegrate projection identity mirroring R2 (LHS `NIntegrate[wStar*(smallSigmaX-sigmaStarX)]`
vs algebraic `rhsMoment`).
(c) CONSULT Q5 source-centering `NIntegrate[sigmaStarX*wStar, ...] == 0`.

**Assessment:** The single-point blindness is closed — `D[..., x]` proves the offset is
constant for ALL x, so a residue that merely happens to vanish at x=1/2 now fails. The
projection and centering checks use Wolfram's NIntegrate, an implementation independent of
SymPy's integrator and of the typed kernel. The centering assertion is present and passing.
SANCTIONED deviation verified: `Chop[..., 10^-25]` is applied to the symbolic zero-derivative
residual (cleaning FullSimplify near-zero noise, NOT loosening — the comparison is still
`=== 0` via expectZero); the centering NIntegrate uses higher WorkingPrecision -> 60 with
AccuracyGoal/PrecisionGoal -> 30 and MaxRecursion -> 30; comparison tolerances are unchanged
(projection 10^-20, centering 10^-20, x-independence exact). No threshold was loosened — the
higher precision tightens, not loosens. The `NIntegrate::precw` warnings in the log are
precision-tracking advisories; both integrals still return exact-zero residuals (`= 0`) and PASS.

### R4 — tautological_check (SymPy g_*/S_* resubstitution repeated own definitions)

**Classification:** resolved

**What changed:** `scripts/...sympy_audit.py:130-143`. The old
`g_star_resub = gPi.subs(Pi,Pi_star)` / `S_star_resub = Sformula.subs(Pi,Pi_star)`
(byte-for-byte repeats of py:25-26) is replaced by quadrature of the source-moment integrals:
`g_star_moment = ∫ Sigma_*·c dx`, `S_star_moment = ∫ Sigma_*·Kq dx`, compared to the closed
forms `g_star`/`S_star` to 1e-25.

**Assessment:** Independent integrator route. `Sigma_star_x` is built directly from `Pi_star`
and the appendix kernels `c`/`Kq`; the integral never references `gPi`/`Sformula`, so a
transcription error in either closed form (e.g. wrong `4Pi^2+pi^2` denominator) would be
caught. The old resub could only catch mutation between blocks. Matches appendix
`eq:app-part04-gbar-Sbar`. Output: `PASS: g_*, S_* equal their source-moment integrals
(residual < 1e-25)`. Non-tautological.

### R5 — transliteration (Mathematica block a verbatim port of SymPy)

**Classification:** resolved

**What changed:** `mathematica/...mathematica_audit.wl:74-83` (A_T) and `:131-141`
(g_*/S_*). The verbatim chain-rule port (`dTmDSigma`/`dSigmaDPi`) is replaced by Wolfram
`D[]` autodiff: `tmOfP = Sqrt[(9/20)*(p/(1-sFormula/4))]`, `dTmDp = D[tmOfP, p]`,
`dgDp = D[gFormula, p]`. The resub port is replaced by `NIntegrate` moment checks
(`gStarMoment = NIntegrate[sigmaStarXm*c, ...]`, `sStarMoment = NIntegrate[sigmaStarXm*kq, ...]`).

**Assessment:** The Mathematica engine now uses independent IMPLEMENTATIONS (Wolfram
`D[]`/`NIntegrate` vs SymPy `sp.diff`/`sp.integrate`) and both are independent of the typed
closed forms `aT`/`gFormula`/`sFormula`. The shared cross-engine STRATEGY (CAS-diff +
quadrature) is by design (consult LOW item 3) — the anti-tautology requirement is independence
from the TARGET's primitive, which is met in both engines. Outputs:
`PASS: A_T closed form vs autodiff of T_m(p)`, `PASS: g_* equals source-moment integral`,
`PASS: S_* equals source-moment integral`. Non-transliteration.

### R6 / ratio label (consult Q6)

**Classification:** resolved

**What changed:** SymPy py:55,62-64 renamed `ratio_paper -> ratio_crosscheck`; print/assert
text now reads `|A_T|/B_T computed ratio cross-check (not a paper literal) matches 31.6785 to
1e-3`. Mathematica wl:66,71-72 renamed `ratioPaper -> ratioCrosscheck` with the matching
relabel. Confirmed in the diff (py:119-132, wl:9-25).

**Assessment:** The misleading paper-quoted claim is removed; it now describes the script's
own computed ratio cross-check. The check still can-fail (a corrupted A_T or B_T moves the
ratio out of 1e-3). The genuine paper literals A_T `-4.27263956256927` and B_T
`0.134875005736706` are KEPT (py:53-54 / wl:64-65) and confirmed verbatim against
`stage_appendix_part04.tex:846,848`. I independently confirmed `31.6785` does NOT appear in
that appendix, corroborating the honesty flag.

## Exec log assessment

**SymPy:** exit=0. PASS lines:
- `PASS: A_T matches paper-quoted -4.27263956256927 to 1e-12`
- `PASS: B_T matches paper-quoted 0.134875005736706 to 1e-12`
- `PASS: |A_T|/B_T computed ratio cross-check (not a paper literal) matches 31.6785 to 1e-3`
- `PASS: A_T closed form agrees with autodiff of T_m(Pi) (residual < 1e-20)`  [R1]
- `PASS: kernel projection of W_* reproduces two-moment traction shift (residual < 1e-22)`  [R2 projection]
- `PASS: rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0, residual < 1e-22)`  [R2 CONSULT centering]
- `PASS: g_*, S_* equal their source-moment integrals (residual < 1e-25)`  [R4]
The autodiff, projection, source-centering, and moment PASS lines are all present. (The
implicit normalization assert is inside the projection block, no separate print.)

**Mathematica:** exit=0. PASS lines:
- `PASS: A_T vs paper -4.27263956256927`
- `PASS: B_T vs paper 0.134875005736706`
- `PASS: |A_T|/B_T computed ratio cross-check (not a paper literal) vs 31.6785`
- `PASS: A_T closed form vs autodiff of T_m(p)`  [R5a autodiff]
- `PASS: W_* centering offset is x-independent`  [R3 x-independence — replaces x=1/2 sample]
- `PASS: deformation normalized`
- `PASS: kernel projection reproduces two-moment traction shift`  [R3 projection]
- `PASS: rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0)`  [R3 CONSULT centering]
- `PASS: g_* equals source-moment integral`  [R5b]
- `PASS: S_* equals source-moment integral`  [R5b]
Autodiff, x-independence, projection, source-centering, and both moment PASS lines present.
Two `NIntegrate::precw` precision-tracking warnings appear; both affected integrals still
return exact-zero residuals (`= 0`) and PASS — advisory only, not failures.

## Material-change assessment

`material_change`: false. No derived numeric result changed: A_T = -4.27263956256927...,
B_T = 0.134875005736706..., |A_T|/B_T = 31.6785..., Pi_*, Sigma_*, T_*, g_*, S_* are all
identical to pre-fix values (this is verification-surface strengthening — the tautological
checks were replaced by independent ones, and one centering assertion + one label fix were
added). No downstream unit's inputs are affected.

## Side observations (non-blocking)

- The two `NIntegrate::precw` warnings (SymPy log clean of warnings) arise because the
  high-WorkingPrecision integrands lose interior precision; harmless here since the residuals
  evaluate to exact 0. If a future pass wants a cleaner transcript, wrapping the integrand in
  `SetPrecision` (already done for the centering integral) on the projection integral would
  suppress them. Non-blocking.
- The print-only `Wcenter`/`wCenter` `simplify`/`FullSimplify` dumps (py:89-91, wl:88-89)
  were correctly left intact per the directive — they are output, not assertions.

## Verdict justification

All five Codex findings (R1-R5) are resolved with genuinely independent primitives: A_T via
CAS autodiff of T_m(Pi)/T_m(p) (not the hand-typed chain rule); the kernel via numerical
projection against a non-canonical 2x deformation plus full-symbolic x-independence (not a
single x=1/2 sample); g_*/S_* via quadrature of the source moments (not resubstitution). The
consult-added source-centering assertion `∫ Σ_* W_* dx == 0` is present AND passing in BOTH
engines, and it is correctly reasoned to catch dropped centering constants. The ratio
mislabel (R6/Q6) is corrected to a non-paper-sourced computed cross-check while the genuine
A_T/B_T paper literals are kept and confirmed against appendix:846-848. The sanctioned R3
deviation (Chop + higher-precision NIntegrate) loosens no threshold. Both engines exit 0 with
the full expected PASS-line sets and freshly regenerated outputs. No regressions in the diff.
