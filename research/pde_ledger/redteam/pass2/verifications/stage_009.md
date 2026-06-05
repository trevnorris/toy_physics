---
unit_id: 009
batch: I.1
verifier_model: claude-opus-4-8
verify_date: 2026-06-04T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 009

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
Codex added a generic-kernel first-moment check in BOTH engines (direction (a)),
using a normalized Gamma one-sided kernel `w(u)=u·e^{-u}` (μ₁≠1) alongside the
intact exponential representative:

- SymPy `scripts/...stage009_..._sympy_audit.py:119-128`:
  - `gamma_norm = sp.integrate(u·e^{-u}, (u,0,∞))` → asserted `=1` (line 125)
  - `gamma_mu1 = sp.integrate(u·(u·e^{-u}), (u,0,∞))` → asserted `=2` (line 126)
  - `avg_Q_gamma = ∫₀^∞ (w/ℓ)e^{-w/ℓ}/ℓ · Qw dw`, leading term via `series(...,ell,0,2)`,
    asserted `avg_Q_gamma_leading - (q0 + ell*gamma_mu1*q1) == 0` (lines 122-128).
- Mathematica `mathematica/...stage009_..._mathematica_audit.wl:31-42` (M1d/M1e/M1f):
  the same kernel, via independent `Integrate`/`Series`, asserting `gammaNorm-1`,
  `gammaMu1-2`, and `avgQGammaLeading - (q0 + ell·gammaMu1·q1)` all `==0`.
  (`ClearAll` at line 18 also adds `u`; benign hygiene.)

**Assessment:**
Correct and matches the directive's direction (a) requirement. The new check is
NON-TAUTOLOGICAL and genuinely μ₁-DEPENDENT:

- `gamma_mu1` is COMPUTED by independent symbolic integration of `∫₀^∞ u·w(u) du`
  (yielding 2), not hardcoded. The Gamma mouth average is COMPUTED by an independent
  integral of the kernel against Q. The two sides reach the assertion via separate routes.
- The independently-integrated leading mouth coefficient is `2·q1` (SymPy output L47:
  `<Q>_ell = ... + 2*ell*q1 + q0`). The RHS uses the computed `gamma_mu1=2`. Had the
  coefficient been the kernel-independent exponential value `q1` (μ₁=1, RHS `q0+ell*q1`),
  the residual would be `2·ell·q1 − ell·q1 = ell·q1 ≠ 0` → the assertion would FAIL.
  So the coefficient is genuinely tied to `∫₀^∞ u w(u) du`, not the collapsed 1.
- The pre-existing exponential-kernel checks (SymPy 108-117 / `.wl` M1a-M1c) are intact
  and unchanged; no previously-asserted value changed. The added print lines (134-135)
  are informational only.

No collateral edits beyond the additive checks and explanatory prints.

## Exec log assessment

The named exec logs (`stage_009_sympy.log`, `stage_009_mathematica.log`) are absent
from `redteam/pass2/exec_logs/` (only `stage_009_diff.patch` present). Per the
orchestrator note, both engines were re-run and exit 0; I assessed the committed,
refreshed `.txt` outputs instead (allowed substitute — outputs are the run record).

**SymPy:** exit=0 (inferred from refreshed output ending `STATUS: PASS`, no AssertionError).
Notable lines: `μ1 = 2`; `<Q>_ell = 5*ell**4*q4 + 4*ell**3*q3 + 3*ell**2*q2 + 2*ell*q1 + q0`
(Gamma kernel, leading `2*ell*q1`); `so its O(ell) term is ell*μ1*Q'(0), not the
exponential value ell*Q'(0).`; `STATUS: PASS`.

**Mathematica:** exit=0. Notable lines:
`OK: M1d Gamma one-sided kernel normalization residual = 0`,
`OK: M1e Gamma one-sided first moment residual = 0`,
`OK: M1f Gamma mouth first-moment leading term residual = 0`, `STATUS: PASS`
(all 14 M-checks residual 0).

**Output freshness:** confirmed. Both output `.txt` mtimes are 2026-06-04 20:26;
both scripts are 2026-06-04 19:15 — outputs are newer than scripts (re-generated post-fix).

## Material-change assessment

`material_change`: false.

Purely additive coverage. The new Gamma-kernel μ₁ check introduces no change to any
previously-asserted value (exponential μ₁=1 results, even-kernel moments, μ_eff/ξ_eff
series, and Gaussian fingerprints are byte-for-byte the same in the refreshed output).
No downstream unit can depend on a value that did not change. No `upstream_stale`
concern beyond the orchestrator's default bookkeeping.

## Side observations (non-blocking)

- The `.wl` `ClearAll[w, u, ...]` (line 18) now also clears `u`; harmless, just keeps
  the Gamma kernel's `u` symbol clean. No effect on prior M-blocks.
- The SymPy Gamma leading term is taken via `series(...,ell,0,2).removeO()` (O(ℓ²)
  truncation) which correctly isolates the `q0 + ell·μ₁·q1` leading structure; the
  full `<Q>_ell` printed (`5ℓ⁴q4+4ℓ³q3+3ℓ²q2+2ℓq1+q0`) confirms μ₁=2 acts on the
  O(ℓ) term specifically.

## Verdict justification

The single low-severity finding F1 is resolved exactly as the user-approved directive
(direction (a)) required: a generic-kernel first-moment check was added in both engines
using a Gamma one-sided kernel with μ₁=2, and the assertion ties the O(ℓ) mouth
coefficient to the independently-computed `∫₀^∞ u w(u) du` rather than the collapsed
exponential value 1 — it would fail if the coefficient were the kernel-independent
`q1`. Both engines exit 0 with all checks (including the new M1d/M1e/M1f and the SymPy
Gamma asserts) at residual 0, outputs are freshly regenerated, the pre-existing
exponential-kernel checks are intact, and no prior value changed. Verdict: verified;
material_change: false.
