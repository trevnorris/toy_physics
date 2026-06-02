---
unit_id: 210
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 210

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the second-engine script
`mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl`
(registered/canonical path after the orchestrator's filename-convention rename; the
rename is purely a glob-convention fix and is NOT treated as a defect). The `.wl`
verifies the directive's full claim manifest M1–M9 with symbolic residual checks
(`expectZero`) plus one `expectTrue`, aborting via `Exit[1]` on any nonzero residual.
The script carries the required symbol domains in `$Assumptions`
(`ki,kj,kk>0`, `ai,aj,ak≥0`, `H0>0`, `r,s,nu≥0`, six `u_••` real). The independent
re-run exits 0 with 36 PASS / 0 FAIL.

**Assessment:**
The edit fully discharges the dual-engine requirement and is correct. Every manifest
item is represented by an explicit, non-tautological check:

- M1 — three edge normalizations `a·a−1` (l.92–94).
- M2 — boundary slopes obtained by `Limit[…, eps→0, FromAbove]` on a perturbed 3D
  patch, checked against the pairwise cone forms (l.98–113).
- M3 — gradient-optimal point DERIVED via an explicit Lagrangian + `D[]` +
  `Solve[…,Reals]` + `SelectFirst` for the positive branch, then matched to `k/|k|`;
  also the Cauchy–Schwarz identity (l.118–157).
- M4 — max-slope value and the two interior ratios (l.161–163).
- M5 — three Pythagorean decompositions encoding the strict triple-over-pairwise gain
  (l.167–169).
- M6 — `w_Σ` sum-square identity, Cauchy slack, slack-bound, equal-mix (=2) and
  pairwise-edge (=1) screen values, and the equal-mix UNIQUENESS via `Reduce`
  (l.191–203).
- M7 — discriminant coefficients A–F extracted by `Series`+`CoefficientList`/
  `Coefficient` and matched term-by-term (l.227–243). Load-bearing item handled by
  coefficient extraction, exactly as the anti-transliteration guard prescribed.
- M8 — three curvature edge reductions (again via `Limit`) plus the diagonal-neutral
  case (l.248–268).
- M9 — ratio-coordinate τ map and its three boundary reductions (s=0 → ij, r=0 → ik,
  and the direct jk-edge form) (l.272–313).

No check is of the form `x=expr; assert x==expr`: each object is defined one way
(e.g. κ as `vec.H.vec`, τ as the closed root map, A–F as series coefficients) and
checked against an independently written target form. No collateral edits; the
`.wl` and its `.txt` are the only new files. `deviation: none` in the directive is
accurate.

## Independent-derivation check (CRITICAL — load-bearing)

**Result: confirmed independent.** The `.wl` is a genuine re-derivation using native
Mathematica primitives, not a line-by-line port of the SymPy `.py`. Three corresponding
sections demonstrate a different decomposition:

1. **Gradient-optimal point (M3).** The `.wl` constructs an explicit Lagrangian and
   solves the stationarity system:
   `lagrangian = x ki + y kj + z kk - lag (x^2+y^2+z^2-1)`;
   `Solve[{D[lagrangian,x]==0, …, x^2+y^2+z^2==1}, {x,y,z,lag}, Reals]`, then
   `SelectFirst` for the strictly-positive branch.
   The SymPy `.py` instead merely POSITS the closed form
   `avec_grad = Matrix([ki/Kgrad, kj/Kgrad, kk/Kgrad])` and checks normalization/slope.
   Different route (constrained optimization vs. asserted closed form).

2. **Edge reductions (M2 / M8).** The `.wl` uses limiting patches:
   `Limit[slopeOf[normalizeRaw[{1, r, eps}]], eps -> 0, Direction -> "FromAbove"]`.
   The SymPy `.py` uses direct substitution of the 2D edge vector
   `k_simplex.subs({ai: avec_ij[0], aj: avec_ij[1], ak: avec_ij[2]})`.
   Different decomposition (eps→0 limit of a 3D interior patch vs. direct 2D subs).

3. **Discriminant numerator (M7).** The `.wl` extracts coefficients:
   `deltaSeries = Series[…,{r,0,2},{s,0,2}]`; `CoefficientList[deltaSeries,{r,s}]` and
   per-monomial `Coefficient[Coefficient[deltaSeries,r,i],s,j]`, matching A–F one at a
   time. The SymPy `.py` does a single monolithic clearing identity
   `sp.simplify((1+r²+s²)(k_rs²−2H₀κ_rs) − Delta_sharp)`.
   Different decomposition (coefficient extraction vs. cleared-and-simplified residual).

The `.wl` also adds checks the `.py` does not contain (M3 explicit Lagrange branch,
M3 Cauchy–Schwarz identity, M6 equal-mix uniqueness via `Reduce`), further confirming
it was written from the math rather than transcribed. No shared intermediate-variable
choreography or substitution order. This is independent, not transliteration.

## Exec log assessment

**SymPy:** exit=n/a. The SymPy side is unchanged for this unit (F1 only added a
second engine); no new SymPy run was required.

**Mathematica:** exit=0. Notable lines:
- `M3 positive Lagrange branch = {ki/Sqrt[ki^2+kj^2+kk^2], …}` and
  `M3 Lagrange stationarity residuals = {0, 0, 0, 0}` — the optimizer was genuinely
  solved, not asserted.
- `M7 CoefficientList in {r,s} = {{ki^2-2*H0*uii, …}, …}` followed by
  `M7 coefficient A/B/C/D/E/F residual = 0` — per-monomial match.
- `M6 equal-mix uniqueness condition = True` (the `Reduce` output equals the
  barycenter).
- Terminating banner `STAGE 210 MATHEMATICA AUDIT PASSED`, `# exit_code: 0`.
36 PASS lines, 0 FAIL. (Prompt's "~37" is approximate; the exact count is 36 = 35
`expectZero` calls + 1 `expectTrue` M6-uniqueness call. Every M1–M9 item is covered.)

**Output freshness:** confirmed.
- `.wl` mtime 2026-06-02 10:15:42
- `.txt` mtime 2026-06-02 10:41:41 (newer than the script — fresh)
The saved transcript reflects the current script.

## Material-change assessment

`material_change`: false. F1 only ADDED a second-engine verification script; it changed
no derived result. The SymPy `.py` and all upstream/downstream numeric content are
untouched. No downstream unit (>210) depends on any new value introduced here.

## Side observations (non-blocking)

- The `.wl` and its `.txt` are git-untracked (`??`), which is why
  `stage_210_diff.patch` is empty (git diff vs. HEAD shows nothing for new files).
  This is the expected state for a brand-new second engine and is not a defect; the
  orchestrator will commit them with the batch.
- The SymPy banner reads "STAGE 193 … (canonical triple-screen)" while the unit is
  paper-stage 210 — the known internal-renumbering cosmetic label, noted in the
  original report, carries no math, and is out of red-team scope. The new `.wl`
  correctly banners "STAGE 210", so the second engine does not propagate the label.

## Verdict justification

The single finding (missing Mathematica engine) is resolved: a correct, independent
`.wl` exists at the registered path, exits 0 with 36/36 PASS and no FAIL, covers the
entire M1–M9 manifest with non-tautological assertions, and is a genuine re-derivation
(explicit Lagrange optimization, limiting-patch edge reductions, and `Series`/
`CoefficientList` coefficient extraction) rather than a port of the SymPy choreography.
Output is fresh, no regressions, no material change. Verdict: verified.
