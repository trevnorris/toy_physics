---
unit_id: 231
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T22:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 231

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved

**What changed:**
`notes/stages/moving_throat_pde_stage231_..._sympy_audit.md:98` — the `\partial_\xi F`
numerator polynomial was corrected from `81\delta^3+240\delta^2\xi+72\delta^2+297\delta\xi^2+189\xi^3`
to `81\delta^3+189\delta^2\xi+72\delta^2+297\delta\xi^2+121\xi^3`. Git diff confirms exactly the
two authorized coefficient swaps (`240\delta^2\xi → 189\delta^2\xi`, `189\xi^3 → 121\xi^3`),
one line changed, nothing else in `notes/` or `paper/` touched for unit 231.

**Assessment:**
Correct and exactly matches the directive's RESOLVED-F1 instruction. The script polynomial
(`.py:81-84`) was NOT altered — it already held the correct 189/121 coefficients. The corrected
notes line now agrees with: (a) SymPy's own `sp.factor(sp.diff(F,xi))` (sympy output line 11
shows `189*delta**2*xi ... 121*xi**3`); and (b) the new `.wl`'s independent `Factor[D[F,x]]`
(mathematica output line 8: `... + 189*d^2*x + ... + 121*x^3`). The notes were the erroneous
side; bringing prose into line with the already-verified script is a non-tautological,
correctly-targeted fix.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage231_..._mathematica_audit.wl` (197 lines)
establishing claim manifest M1–M7 with hard `Exit[1]`-on-mismatch failure (`fail[]`, lines 8-12;
`expectZero`/`expectTrue`/`expectClose` all route to it).

**Assessment — genuinely independent route, not a transliteration:**
- Monotonicity signs (M1 dF>0/dG>0/dRND<0, M4 dPc>0) are proved by symbolic
  `Resolve[ForAll[{x,d}, Implies[...]], Reals]` — a strictly stronger universal-quantifier proof,
  whereas the `.py` proves the same signs by *numerical grid sampling* of lambdified handles.
  Different decomposition.
- Threshold roots (M5 R_flip/R_den) are found by native `NSolve[poly==0, x, Reals]` with a
  unique-root `Select`, whereas the `.py` uses a hand-coded `bisect_increasing` bisection.
- Soft limit (M2) uses `Limit[..., Direction->"FromBelow"]` + the reciprocal==0 pole idiom;
  the `.py` uses `sp.limit(..., dir="-") == oo`.
- Variable names differ (x/d/c vs xi/delta/c); grid points are exact rationals
  ({1/100,1/5},…) vs the `.py`'s floats. Same physical points, different representation.
- M1 derives `Factor[D[F,x]]` itself (line 76 `dF = Factor[D[f,x]]`; lines 79-81 strip the
  `(9d+11x)^3` factor and compare to `dFPublishedPoly = 81 d^3 + 189 d^2 x + 72 d^2 + 297 d x^2
  + 121 x^3`). It REPORTS Mathematica's own factor (output line 8) carrying coefficients 189 and
  121 — independently adjudicating the F1 notes discrepancy as required by the directive.

No transliteration markers (var names, assertion order, and choreography all diverge from the
`.py`). FAIL path is non-vacuous: a wrong `dFPublishedPoly` coefficient or a diverged threshold
would trip `Exit[1]`/`expectClose`. All M1–M7 pass; exit 0.

### F3 — hardcoded_result

**Classification:** resolved

**What changed:**
`.py:180` — provenance comment added naming the Stage 230 source script
(`scripts/moving_throat_pde_stage230_..._sympy_audit.py`) from which `s_minus_den`/`s_minus_num`
are carried. Values unchanged (`0.411024574532864`, `-0.334368725711457`).

**Assessment:**
Annotation-only, exactly as directed; the numeric output is unchanged (sympy output line 30
still `R_* = 1.22925543846333...`). The R_* derivation from the s_minus pair remains
independently pinned by `assert_close(..., 1.229255438463336)` at line 187. Non-blocking and correct.

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
`.py:271` — provenance comment naming the upstream budget source; `.py:282-286` — a new loop over
`probe_grid` asserting `0 <= R_num(xi_val, delta_val) <= R_num(0.0, delta_val)`
(i.e. sampled `R_phys ∈ [0, R_ND(0,delta)]`, the physical subset the survival argument relies on).
The four budget literals (lines 272-275) are unchanged.

**Assessment:**
The subset-bound check is non-tautological: `R_num` is a genuine non-constant rational function
of both xi and delta; sampling at six distinct grid points exercises the classifier (a wrong
classifier or a point outside the physical subset would raise AssertionError). It connects the
otherwise free-floating literal-vs-literal inequality to the pullback, as the directive required.
The `.wl` independently mirrors this with its own M7 subset samples (output lines 91-108, all
`sample <= onset` True). Codex's documented `deviation` on F4 (citing "Stage 247 old numbering /
canonical Stage 230" because current Stage 247 script does not define these budgets) is an honest
provenance annotation, not a math change — acceptable.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L11: `dF/dxi = ...(81*delta**3 + 189*delta**2*xi + 72*delta**2 + 297*delta*xi**2 + 121*xi**3)...`
  — confirms the canonical 189/121 polynomial the F1 notes edit was aligned to.
- L30-31: `R_* = 1.22925543846333...`, `delta_dyn_* = 0.72311161787501...` (F3 carry intact).
- L53: `All Stage 231 symbolic and numerical audits passed.`

**Mathematica:** exit=0. Notable lines:
- L8: `M1 ... D[F,x] numerator polynomial factor = 72*d^2 + 81*d^3 + 189*d^2*x + 297*d*x^2 + 121*x^3`
  — Mathematica's own `Factor`, independently corroborating 189/121.
- L14: `PASS: M1 D[F,x] numerator polynomial factor`.
- L57-77: M5 R_flip/R_den match SymPy sample rows to ~2e-10 (`diff ... *^-10`), `R_flip <= R_den` True.
- L109: `All Stage 231 Mathematica audits passed.`
- PASS count: 39 PASS, 0 FAIL. Per claim: M1=5, M2=4, M3=6, M4=4, M5=11, M6=1, M7=8 — every
  manifest item present, no silent parser-skip.

**Output freshness:** confirmed. Output mtimes (1780459295 for both `.txt`) are newer than the
`.py` (1780458714) and `.wl` (1780458700) source mtimes — both outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false.

The only value-bearing edit is the F1 notes coefficient correction, which brought erroneous prose
into agreement with the already-correct (and already-passing) script polynomial. No derived or
carried *numeric* result changed: the `.py` was unchanged on the polynomial, R_* / delta_dyn_* /
sample thresholds / budgets are all identical to the pre-fix output. F3 and F4 are
annotation + an additive non-vacuous check; the new `.wl` is additive cross-check coverage. No
downstream unit depends on a changed value, so no `upstream_stale` propagation is warranted on
numeric grounds (any blanket `>231` staleness the orchestrator applies is procedural, not
value-driven here).

## Side observations (non-blocking)

- The working tree also contains modified/untracked files for unit 232 (notes, `.wl`, outputs,
  MANIFEST, reports/directives 232-242). These are outside unit 231's scope and were not assessed;
  the 231-scoped changes are clean and self-contained (one notes line, the `.py`, its `.txt`, and
  the new `.wl` + its `.txt`).
- M7 is the one place both engines compare the same four hardcoded budget literals, but this is
  exactly what the directive's M7 specified (carried-forward upstream budgets → assert the
  inequalities), and both engines back it with the independent subset-inclusion corroboration.
  Not a defect.

## Verdict justification

All four findings are `resolved`. The F1 notes edit is exactly the two authorized coefficient
swaps and nothing else in prose was touched; the script polynomial was correctly left unchanged
and is independently confirmed by both engines. The new Mathematica `.wl` is a genuinely
independent route (symbolic `Resolve[ForAll]` vs numeric grid; `NSolve` vs bisection; native
`Limit`; different var names and exact-rational grid), it derives `Factor[D[F,x]]` itself and
reports the corrected 189/121 coefficients, and its FAIL path is a non-vacuous `Exit[1]`. F3/F4
are correct annotation plus an additive, non-tautological subset-bound check. Both engines exit 0
with all manifest claims present (39 PASS / 0 FAIL on the Mathematica side; 20 asserts on the
SymPy side) and outputs are freshly regenerated. Verdict: verified.
