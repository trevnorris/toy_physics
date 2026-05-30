---
unit_id: 146
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T22:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 146

Remediation pass (overwrites the stale 2026-05-27 verification, which
rubber-stamped the loose 1e-15 / 10^-6 tolerances that the Codex review then
flagged). Source of findings is `redteam/codex_reviews/stage_146.md` (R1 SymPy
affine-law tolerance, R2 Mathematica Chop-mask, R3 banner), as encoded in
`directives/stage_146.md` with `## Applied:` blocks, and refined by
`redteam/codex_reviews/_consult_batch6.md` Q3. Scripts read post-fix; not
executed (used captured exec logs + committed outputs).

## Per-finding outcomes

### R1 — insufficient_verification (SymPy affine-law tolerance)

**Classification:** resolved

**What changed:**
- `scripts/...stage146..._sympy_audit.py:73` — root is now
  `Pi_star = sp.nsolve(gPi - gminus, 1.5, prec=50)`. The old
  `sp.N(sp.nsolve(gPi - gminus, 1.5), 30)` zero-pad-only form is GONE
  (diff hunk line 73). `prec=50` solves the root at ~50 digits, not the
  ~15-digit `nsolve` default.
- `:114-127` — residual now evaluated at `sp.N(..., 50)`, compared with sympy
  `Float` arithmetic against `TOL = sp.Float("1e-25", 50)` via
  `if sp.Abs(g_res) >= TOL:` / `if sp.Abs(S_res) >= TOL:`. The old
  `abs(float(g_res)) > 1e-15` / `> 1e-15` (which collapsed to a ~15-16 digit
  Python float before comparing) is GONE.

**Assessment:** Edit matches directive Edits 1a + 1b exactly; no collateral
change in the touched hunk. Non-tautological: the residual is the
directly-integrated `gbar_phys - (gminus + eps*(gbar_v - gminus))`, an
independent quadrature against the physical kernels, not an algebraic
restatement of the affine form. Printed raw residuals are genuinely `< 1e-25`:
- g_eps eps=1/10: `-1.43...E-51`; eps=1/2: `-5.92...E-52`
- S_eps eps=1/10: `2.18...E-50`; eps=1/2: `5.87...E-51`
All ~1e-50/1e-51, deep inside the 1e-25 budget — confirming `prec=50` took
(old residuals sat at ~3.6e-18 at the 15-digit floor). The CONSULT Q3 risk
(`prec=50` failing to tighten the root) did NOT materialize.

### R2 — insufficient_verification (Mathematica Chop hides residual)

**Classification:** resolved

**What changed:**
- `mathematica/...stage146..._mathematica_audit.wl:103-121` — both
  `Chop[..., 10^-6]` wrappers REMOVED; raw residuals printed via
  `fmt[gEpsSample1]` etc.; both `If[... < 10^-6 ...]` guards replaced by
  `Abs[...] < 10^-25`. The two stale precision-9 / 10^-6 justification
  comments are gone. A grep for `Chop` in the live `.wl` returns only one hit,
  inside the new comment ("no Chop") — no live `Chop[...]` call survives on
  these residuals.

**Assessment:** Matches directive R2. SANCTIONED DEVIATION confirmed and
honest: the `## Applied: R2` block declares Codex used the directive's
symbolic endpoint-integral fallback —
`gEpsRes = ((1 - eps)*(gDirect /. p -> pStar) + eps*gBarV) - (gMinus + eps*(gBarV - gMinus))`
(and the `sDirect` analogue) — because numeric-`pStar` `Integrate` of the
combined `sigmaEps` profile produced only low-precision zeros. This is exactly
the fallback the directive authorized (lines 316-319): substitute `pStar` into
the F2-anchored closed-form endpoint integrals (`gDirect`/`sDirect`) and
combine. It did NOT re-loosen the tolerance (assert is `< 10^-25`) and did NOT
re-introduce a Chop mask. The fallback genuinely achieves `< 1e-25`:
- g_eps eps=1/10: `8.03...*^-58`; eps=1/2: `4.46...*^-58`
- S_eps eps=1/10 and 1/2: literal `0` (printed `0``40.`)
All `< 10^-25` with large headroom.

### R3 — paper_misalignment (banner)

**Classification:** resolved

**What changed:** `...stage146..._mathematica_audit.wl:32` now reads
`banner["STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];`
(was `FINITE-CORRECTION`). Em-dash and `STAGE 146` preserved (diff confirms
only `FINITE-CORRECTION` → `FIRST-ORDER`).

**Assessment:** Matches directive R3 byte-for-byte. The committed Mathematica
output (line 3) and exec log (line 8) print the corrected banner
`STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS`.

## CONSULT Q3 caveat (scope of the affine check)

Confirmed correctly stated in the directive's "Anti-tautology guard" /
"CONSULT Q3 caveat" section: the affine residual collapses to
`(1-eps)*(gPi(Pi_*) - gminus)`, and `gminus` ALSO feeds
`Pi_star = nsolve(gPi - gminus, ...)`, so a wrong `gminus` shifts `Pi_star` to
compensate and the residual stays small. This check is therefore NOT an
independent guard against `gminus`; it is credited only with testing
intercept-vs-direct-integral plus the kernels (`cos(pi x/2)` / `Kq`) and source
(`Sigma_*`), which do NOT feed the root. `gminus` is anchored at its owning
branch stage; stage 146 does not re-derive it. Per the verifier overrides this
is a correctly-documented scope limit and does NOT fail verification — R1/R2
only required closing the tolerance/Chop finding, which the fixes do.

## Reconcile anchors (untouched)

Confirmed the direct-integral PASS anchors (review rows 1-3) were NOT
disturbed. The diffs touch only `.py:73` and `.py:103-127`, plus `.wl:32` and
`.wl:92-121`; SymPy lines 33-53 (`g(Pi)`/`S_q(Pi)` direct-formula) and
Mathematica lines 44-51 (`expectZero` direct-formula anchors) are unchanged.
Outputs still show `g(Pi) direct-formula = 0`, the four `S_q(Pi)` numeric
samples `diff=0E-165`, and both Mathematica `direct-formula = 0` PASS lines.

## Exec log assessment

**SymPy:** exit=0. No warnings/errors/tracebacks. Notable lines:
- `g(Pi) direct-formula = 0`
- `g_eps affine law (integral form) at eps=1/2: residual = 5.87...E-51`
- `PASS: g_eps affine law (integral form) at eps=1/10 and eps=1/2`
- `PASS: S_eps affine law (integral form) at eps=1/10 and eps=1/2`

**Mathematica:** exit=0. No messages/errors. Notable lines:
- `STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS`
- `g_eps affine law (integral form) at eps=1/10: 8.02...*^-58`
- `PASS: g_eps affine law (integral form)`
- `S_eps affine law (integral form) at eps=1/10: 0``40...` then
  `PASS: S_eps affine law (integral form)`

**PASS-line counts:** SymPy committed output has 3 `^PASS:` lines (S_q
direct-formula fallback + the two affine-law PASS lines; the `g(Pi)` anchor
prints `= 0` via `expect_zero`, not a `PASS:` prefix). Mathematica committed
output has 8 `^PASS:` lines (g/S direct-formula, 3 kernel samples, Pi_*
compensation, g_eps, S_eps). Counts are consistent with script structure;
passing logs are necessary-but-not-sufficient and the substance checks above
confirm the PASSes are earned (raw residuals < tolerance, no masking).

**Output freshness:** committed outputs (both 21:46) are newer than the live
scripts (`.py` 21:37, `.wl` 21:39) — regenerated post-fix.

## Material-change assessment

`material_change`: false. Only tolerances and root-finder precision tightened
(`nsolve` 15→50 digits, asserts 1e-15/10^-6 → 1e-25/10^-25, residual eval
30→50 digits, Chop removed). No carried/derived constant moved: `Pi_*`, `g_*`,
`S_*`, `g_*'`, `S_*'` and the retuning-law coefficients match prior printed
values (higher root precision only adds trailing digits below the 30-digit
print). No downstream unit input changed.

## Side observations (non-blocking)

- The two engines use equivalent-but-different residual constructions: SymPy
  uses the directly-integrated `gbar_phys` (combined-profile integral
  evaluates symbolically there), while Mathematica uses the declared
  sanctioned fallback `(1-eps)*(gDirect/.pStar)+eps*gBarV` (combined-profile
  numeric `Integrate` lost precision). Both reduce to the same root-error
  quantity and both clear `< 1e-25`; the divergence is documented in
  `## Applied: R2 deviation`.
- The SymPy `kernel check` retains its pre-existing `1e-12` tolerance
  (`.py:67`), explicitly listed under the directive's "Reconcile … untouched"
  section as an independent cross-check; not in scope for R1/R2, correctly
  left alone.

## Verdict justification

All three Codex-review findings are resolved with edits matching the directive
(R1 root precision + 1e-25 Float assert, R2 Chop removal + 10^-25 raw-residual
assert, R3 banner), no collateral edits, both engines exit 0 with clean logs,
and the printed raw residuals are genuinely far below 1e-25 (SymPy ~1e-50/51,
Mathematica ~1e-58 / literal 0). No loose tolerance, Chop mask, or missed
`prec=50` survives. The reconciled direct-integral anchors are undisturbed and
the CONSULT Q3 coupled-`gminus` scope limit is correctly documented. No derived
numeric result moved → material_change false.
