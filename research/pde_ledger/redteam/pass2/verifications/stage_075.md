---
unit_id: 075
batch: III.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T22:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 075

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy (`scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:117-136`): the tautological round-trip block (old lines 117-125 — the misleading "not the trivial identity 100*Theta == 100*Theta" comment, `Upsilon_*_from_Theta = alpha_r**2 * Theta_*`, and both `assert sp.simplify(Upsilon_* - Upsilon_*_from_Theta) == 0`) was deleted. In its place a `expect_close(name, value, target, rel_tol)` helper was added that compares `sp.N(value, 30)` against `sp.Float(target, 30)` with limit `rel_tol * max(1, |target|)`, plus eight calls anchoring `Delta_0`, `Delta_inf`, `Upsilon_fail/suff`, `Xi_fail/suff`, `Theta_fail`, `Theta_suff` against the fixed paper/notes literals (`1.73302079021525e-4`, `2.01447565540522e-2`, `0.0362605617972939`, `4.21495341569977`, `49.6407091004953`, `5770.27122609299`, `3.62605617972939e-4`, `4.21495341569977e-2`).
- Mathematica (`mathematica/...audit.wl:113-116` region): the two `expectZero["Upsilon_* - alphaR^2 * Theta_*", upsilon* - alphaR^2*theta*]` round-trip lines were deleted; nothing else in the file changed. The eight `expectApprox` numeric anchors (wl:116-123) and the `expectZero["alpha_r^2 - 100 ...", alphaR^2 - 100]` lock (wl:41) remain intact.

**Assessment:**
Correct and complete on both required axes.
(a) The tautology is genuinely gone from both engines — `grep` of both the live scripts and the committed `.txt` outputs shows no remaining `Upsilon - alpha_r^2*Theta` round-trip. The diff touches only the two intended regions; no collateral edits.
(b) The reduction's real content stays locked: SymPy `assert alpha_r**2 == 100` (py:28) and Mathematica `expectZero[alphaR^2 - 100]` (wl:41) are both still present and pass.
(c) The new anchors are non-tautological. `target` is a hardcoded string literal parsed by `sp.Float(target, 30)` / a Mathematica numeric literal — fully independent of the computed closed form `value`. I checked adversarially that the comparison is value-vs-external-literal, not value-vs-itself: e.g. `Theta_fail` is computed as `Upsilon_fail/alpha_r**2` from the closed form, but is checked against the fixed literal `3.62605617972939e-4`; a wrong `alpha_r^2` factor (e.g. dividing by 10 instead of 100) would shift it ~10x, a diff ~3.6e-5 that blows past the `1e-14*|target|` ≈ 3.6e-18 limit. The literals match the paper-stated deliverables per the directive's source map (tex:17/22/28/34, md:45/47/51/53).
(d) The committed outputs show all eight new anchors PASS with tiny residuals (largest `Xi_suff` diff 9.38e-13 vs its `1e-14*5770≈5.8e-11` limit; smallest `Delta_0`/`Theta_fail` diffs ~1.3-1.5e-19), and the printed deliverable numerics are byte-identical to the pre-fix values — no deliverable value moved.
(e) The `.wl` retains its eight independent `expectApprox` anchors (wl:116-123), all reporting `diff = 0` to displayed precision.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `alpha_r^2 (paper Inputs line lock) = 100   PASS` (reduction lock intact)
- `Theta_fail / Pe_req diff = 1.30...E-19` → `Theta_fail / Pe_req PASS` (boxed FAIL endpoint now asserted, previously only printed)
- `Theta_suff / Pe_req diff = 2.87...E-17` → `Theta_suff / Pe_req PASS` (boxed SUCCEED endpoint now asserted)
- No `Upsilon_fail - alpha_r^2 * Theta_fail` round-trip line remains.

**Mathematica:** exit=0. Notable lines:
- `PASS: alpha_r^2 - 100 (paper Inputs line lock)` (lock intact)
- `Theta_fail / Pe_req numeric check diff = 0``23.` → `PASS`; `Theta_suff ... diff = 0``21.` → `PASS`
- `Stage 075 Mathematica audit passed.` and no round-trip lines.

**Output freshness:** confirmed. SymPy `.txt` mtime 1780695899 > script 1780695497; Mathematica `.txt` mtime 1780695923 > script 1780695502. Both committed outputs were regenerated post-fix and their content matches the exec logs (new anchor PASS lines present, round-trip lines absent).

## Material-change assessment

`material_change`: false. The edit only swapped a tautological assertion for genuine literal anchors and deleted two no-op Mathematica checks. Every printed deliverable (`Delta_0`, `Delta_inf`, `Upsilon/Xi/Theta_fail/suff`, `alpha_r^2`) is unchanged digit-for-digit from the pre-fix outputs; no derived result moved, so no downstream unit can be affected.

## Side observations (non-blocking)

- The SymPy `expect_close` tolerance is `rel_tol * max(1, |target|)`, i.e. relative for |target|≥1 and absolute (`1e-14`) for the sub-unity targets like `Delta_0`. Since the sub-unity literals are only carried to ~15 significant figures, the effective check on those is generous but still factor-sensitive (a wrong factor shifts the value by orders of magnitude, far beyond 1e-14). Acceptable; flagged only for completeness.
- `Theta` symbol (`sp.symbols("Theta_w")`, py:98) is now unused after the round-trip removal, but it was already unused before (it never entered the round-trip); harmless leftover, not introduced by this fix.

## Verdict justification

The single finding (F1) is fully resolved. The tautological `Upsilon - alpha_r^2*Theta` round-trip and its false "not the trivial identity" comment are removed from both engines; the reduction's genuine content (`alpha_r^2 == 100`) remains locked in both; and the SymPy script gains eight non-tautological anchors comparing each deliverable against a fixed external paper/notes literal — verified adversarially to be value-vs-literal, not value-vs-self, and tight enough to catch a real factor error. Both engines exit 0, all new anchors pass with tiny residuals, no deliverable value moved, and the committed outputs are fresh. No regressions in the diff. Verified.
