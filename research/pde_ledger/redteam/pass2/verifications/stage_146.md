---
unit_id: 146
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 146

The original report raised two findings (F1 over-claimed g-side label, F2 S-side
tautology). The directive consolidated both into a single user-authorized
de-tautologization F1 applied symmetrically to g and S in BOTH engines. I verify
the consolidated F1, which subsumes the report's F1 and F2.

## Per-finding outcomes

### F1 — tautological_check / insufficient_verification (consolidated; was report F1+F2)

**Classification:** resolved

**What changed:**
Both affine-law blocks were rewritten (sympy py:96-156, mathematica wl:88-123).

- **Closed-form intercepts, no gminus.** sympy py:105-106 sets
  `g_star_closed = sp.N(gPi.subs(Pi, Pi_star), 50)` and
  `S_star_closed = sp.N(Sformula.subs(Pi, Pi_star), 50)`; mathematica wl:96-97 sets
  `gStarClosed = N[gFormula /. p -> pStar, 50]`, `sStarClosed = N[sFormula /. p -> pStar, 50]`.
  The affine RHS at py:123 / wl:116 is `g_*+eps*(gbar_v-g_*)` (resp. S). `gminus` no
  longer appears as the affine intercept anywhere in the affine block; it survives
  only in the upstream Π_* root-solve (py:72-73 / wl:65-66) and the δΠ/δS symbolic
  laws, which are unchanged.
- **Direct assembled-profile quadrature.** The deformed moment is a single quadrature
  of the assembled `Sigma_eps = (1-eps)*Sigma(Π_*)+eps*varsigma_test`:
  py:120-122 `gbar_eps = quad_moment(Sigma_eps_sample*cos(...))`,
  `Sbar_eps = quad_moment(Sigma_eps_sample*Kq)`; wl:113-115 the same via `momentQuad`
  (NIntegrate WP80/AG35). This replaces the old symbolic-linearity decomposition that
  collapsed the residual algebraically.
- **Own-quadrature gbar_v / Sbar_v.** py:107-108 / wl:98-99 compute the test-profile
  moments by their own quadrature.
- **Non-vacuity slope guards.** py:115-118 raise AssertionError if
  `|g_slope| <= 1e-3` or `|S_slope| <= 1e-3`; wl:104-109 `fail[...]` if
  `!TrueQ[Abs[..] > 10^-3]`. Both engines print the slope magnitudes.
- **Honest relabel.** All PASS/label strings now read "convex affine moment law via
  direct quadrature with closed-form g_*/S_* intercept and nonzero slope"; no string
  says "affine law (integral form)" any longer.

**Assessment:**
The edit matches the directive precisely. The residual is now genuinely discriminating:
with `gbar_eps`/`Sbar_eps` computed by an independent quadrature of `Sigma_eps` and the
RHS built from the closed-form `gFormula`/`sFormula` evaluated at Π_*, the residual reduces
to `(1-eps)*(quadrature(Sigma(Π_*)) − closed_form(Π_*))`. A wrong `gFormula`, `sFormula`,
or `Kq` would no longer cancel — the kernel identity is now the falsifiable content,
evaluated AT Π_*. This is no longer x−x (the old S-side defect) nor pure root-closure
`(1-eps)(g(Π_*)−gminus)` (the old g-side mislabel), since the intercept is the closed-form
moment, not the independent algebraic target `gminus`. The two slope guards prove the
eps-term is actually present (slopes ~0.0936 and ~0.0969, far above 1e-3), so neither
check can pass via a vanishing slope.

The moment-formula symbolic verifications at py:33-53 (sympy `integrate(Sigma*cos)` /
`integrate(Sigma*Kq)` vs `gPi`/`Sformula`, with the 4-sample numeric fallback) and
wl:44-51 (`expectZero[gDirect-gFormula]`, `expectZero[sDirect-sFormula]`) are byte-for-byte
intact — the diff begins at py:93 / wl:85, well below those lines. They remain the
cross-engine-independent core.

No collateral edit: the diff touches only the affine block and its leading comment. The
Π_* root-solve, the g_*/S_*/g'_*/S'_* reporting, and the δΠ/δS symbolic laws are untouched.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `g_eps convex affine moment law nonzero slope |gbar_v - g_*| = 0.09359618...` (l.204)
- `S_eps convex affine moment law nonzero slope |Sbar_v - S_*| = 0.09688093...` (l.205)
- `g_eps ... via direct quadrature at eps=1/10: residual = 0` (l.206); `eps=1/2: residual = 1.34E-51` (l.208)
- `S_eps ... via direct quadrature at eps=1/10: residual = -5.35E-51` (l.207); `eps=1/2: -2.67E-51` (l.209)
- Both PASS lines (l.210-211). Π_*, g_*, S_*, g'_*, S'_* (l.32-36) match the report values exactly.

**Mathematica:** exit=0. Notable lines:
- `g_eps ... nonzero slope |gbar_v - g_*| = 0.09359618...` (l.38)
- `S_eps ... nonzero slope |Sbar_v - S_*| = 0.09688093...` (l.39)
- residuals at eps=1/10 and eps=1/2 all print `0` (round-off zero) (l.40-43)
- Both PASS lines (l.44-45). g(Pi)/S_q(Pi) direct-formula PASS (l.13,15), Π_* compensation point PASS (l.31). Π_*, g_*, S_*, g'_*, S'_* (l.25-29) match the report.

Slope magnitudes agree across engines to ~50 digits; all residuals ≤ ~5e-51, well inside
the 1e-25 budget.

**Output freshness:** confirmed. sympy .txt mtime 10:53:35 > script 10:37:11; mathematica
.txt mtime 10:53:35 > script 10:33:44. Both outputs regenerated post-fix.

## Material-change assessment

`material_change`: false.

This is a verification-strength de-tautologization only. No deliverable VALUE moved: Π_*,
g_*, S_*, g'_*, S'_* are bit-identical to the pre-fix audit-recorded values, and the δΠ/δS
symbolic retuning laws were not edited. The change only hardens HOW the affine moment law is
exercised; nothing a downstream unit could consume changed. No paper/notes file was touched
(confirmed: `git status` shows no modified paper/*.tex or notes/*.md).

## Side observations (non-blocking)

- The sympy `quad_moment` uses `sp.Integral(expr,(x,0,1)).evalf(50)` (numerical evaluation),
  i.e. a numeric quadrature distinct from the symbolic-linearity route — consistent with the
  directive's intent that the deformed moment be an independent quadrature. The g-side residual
  printing exactly `0` at eps=1/10 reflects mpmath's numeric quadrature landing on the closed
  form at full precision; the eps=1/2 case shows a genuine ~1e-51 residual, so the check has
  live headroom rather than being a forced zero. Non-blocking.
- The Mathematica `momentQuad` (NIntegrate WP80) is a distinct evaluator from the symbolic
  `Integrate` used in the intact moment-formula core, preserving the second engine's
  independence for the affine-law check as well.

## Verdict justification

Both engines now use closed-form g_*/S_* intercepts (no gminus as the affine intercept),
compute the deformed moment by one direct quadrature of the assembled profile, take gbar_v/Sbar_v
by their own quadrature, and carry the two non-vacuity slope guards. The residual genuinely tests
the quadrature-vs-closed-form kernel identity at Π_* and would fail on a wrong gPi/Sformula/K_q;
it is no longer x−x or pure root-closure. The intact moment-formula symbolic verifications
(py:33-53 / wl:44-51) are untouched. Both exec logs exit 0 with the new PASS lines, nonzero
slopes (~0.094, ~0.097), and residuals ≤ ~5e-51. No deliverable value moved and no prose file
was edited. Verdict: verified; material_change false.
