---
unit_id: 210
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.md]
  paper_appendix: present
---

# Audit unit 210 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_210.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (row 51, narrative rows 236/1091, input row 1343)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.txt`

## What the paper claims

Stage 210 promotes the mixed-ray search to genuine three-coordinate positive simplices `Δ⁺_{ijk} = {a∈ℝ³_{≥0}: a_i²+a_j²+a_k²=1}` and proves seven exact deliverables (notes §intro list): (1) the positive spherical simplex itself; (2) the **boundary-reduction theorem** that each of the three simplex edges (a_k=0, a_j=0, a_i=0) is exactly the corresponding Stage 209 pairwise cone; (3) the **gradient-synergy theorem** giving the unique interior gradient-optimal ray `a_grad = (k_i,k_j,k_k)/√(k_i²+k_j²+k_k²)`, max slope `√(k_i²+k_j²+k_k²)`, and interior ratios `r=k_j/k_i, s=k_k/k_i`, strictly beating every pairwise edge; (4) the **curvature law** `H_1 = aᵀHa` with total off-diagonal weight `w_Σ = 2(a_i a_j+a_i a_k+a_j a_k) = (Σa)²−1 ≤ 2`, maximized uniquely at the equal-mix barycenter `(1,1,1)/√3` (giving w_Σ=2 vs w_Σ=1 on a pairwise equal edge); (5) the **fixed-simplex certified bracket** `τ_{ijk,★}(a)` and its interior ratio-coordinate form with discriminant numerator `Δ# = A+Br+Cs+Dr²+Ers+Fs²` (six exact coefficients A..F) reducing to the Stage 209 pairwise formulas on each edge; (6) the canonical five-row triple-screen packet; (7) the interior-screen dominance criterion. The card's `\stagefield{Output}` reads verbatim: "Boundary reduction for triples, three-coordinate gradient and curvature synergy, canonical triple-screen audit, interior-screen dominance criterion, and non-improvement filter." Status is `Mixed: ExactClosure, Numerical`. Deliverables (6) and (7) (the screen packet bundling and the dominance/non-improvement comparison rules) are organizational/inequality statements, not closed-form algebraic identities, so they are not expected to carry a CAS assertion.

## What the script claims to verify

The SymPy script verifies the algebraic content of deliverables (1)-(5): edge normalizations and the three slope reductions (I), the gradient-optimal vector's normalization/slope/ratios and the three Pythagorean inequalities-of-content `Kgrad²−Kij²−k_k²=0` etc. (II), the cross-leverage identity `w_Σ=(Σa)²−‖a‖²`, the Cauchy slack identity, and the equal-mix vs pairwise screen values 2 and 1 (III), the curvature law's three edge reductions and the diagonal-neutral reduction (IV), the τ root map, interior ratio normalization, the discriminant-numerator reduction, the ratio-coordinate τ form, and three edge reductions of τ (V), and the two canonical interior slope values (VI). The Mathematica script verifies the same five deliverables but, crucially, re-derives several load-bearing objects by independent methods (Lagrange `Solve` for a_grad; `Reduce` for barycenter uniqueness; `Series`+`CoefficientList` for the A..F coefficients; `Limit`-from-interior for the edge reductions). Every check is an `expect_zero`/`expectZero` residual-to-zero or (M6) an `expectTrue` equivalence; none is a bare print.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) positive spherical simplex Δ⁺ | edge normalizations I / M1; ratio-patch normalization V / (implicit) | match |
| (2) boundary reduction = Stage 209 cones | edge slope reductions I / M2 (limit), edge τ reductions V / M9 | match |
| (3) gradient-synergy: a_grad, max slope, ratios, strict gain | II / M3-M5 (Lagrange Solve in .wl) | match |
| (4) curvature law + w_Σ=(Σa)²−1, ≤2, unique barycenter, 2 vs 1 | III+IV / M6-M8 (Reduce uniqueness in .wl) | match |
| (5) fixed-simplex τ, Δ# coeffs A..F, edge reductions | V / M7+M9 (Series/CoefficientList in .wl) | match |
| (6) canonical 5-row screen packet S^tri | (organizational bundle — no algebraic identity) | n/a (not a CAS claim) |
| (7) interior-screen dominance / non-improvement filter | (inequality comparison rule — no closed-form identity) | n/a (not a CAS claim) |
| card `\stagefield{Verification}`: "Mathematica audit: none yet" | a `.wl` exists and passes | mismatch (F1) |

Set `paper_alignment: partial` — the algebraic deliverables all match, but the card's Verification field denies the existence of the present, passing Mathematica audit.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 63-65 | `simplify(edge·edge − 1)==0` | (1) | yes |
| A2 | sympy | 67-78 | `simplify(k_simplex(edge) − closed)==0` | (2) | yes |
| A3 | sympy | 93-96 | grad normalization/slope/ratios `==0` | (3) | yes |
| A4 | sympy | 101-103 | Pythagorean content `Kgrad²−Kij²−k_k²==0` | (3) | yes |
| A5 | sympy | 114-123 | `w_Σ−((Σa)²−‖a‖²)==0`, Cauchy slack `==0` | (4) | yes |
| A6 | sympy | 127-130 | barycenter/edge normalization + w_Σ=2,1 | (4) | yes |
| A7 | sympy | 151-169 | curvature edge + diagonal-neutral reductions `==0` | (4) | yes |
| A8 | sympy | 196-212 | ratio normalization, Δ# reduction, τ form `==0` | (5) | yes |
| A9 | sympy | 227-239 | τ edge reductions (s=0, r=0, jk) `==0` | (5) | yes |
| A10 | sympy | 252-254 | canonical slope values + k_grad²−‖k‖² | (3)/(4) | yes |
| M1 | math | 92-94 | `expectZero` edge normalizations | (1) | yes |
| M2 | math | 111-113 | `expectZero` slope via `Limit` from interior | (2) | yes (independent route) |
| M3 | math | 154-157 | Lagrange `Solve` stationarity residuals + Cauchy identity | (3) | yes (derive-not-posit) |
| M4 | math | 161-163 | max slope + ratios from solved branch | (3) | yes |
| M5 | math | 167-169 | Pythagorean content | (3) | yes |
| M6 | math | 191-203 | w_Σ identities + screen values + `Reduce` uniqueness `expectTrue` | (4) | yes (stronger than .py) |
| M7 | math | 238-243 | `Series`+`CoefficientList` coefficient extraction A..F | (5) | yes (independent route) |
| M8 | math | 261-268 | curvature edge via `Limit` + diagonal-neutral | (4) | yes |
| M9 | math | 310-313 | τ ratio residual + 3 boundary reductions | (5) | yes |

All rows non-tautological: each defines an object by one route and checks it against an independently-stated closed form, or solves/reduces and compares.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_210.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl:1-317`

**What's wrong:**
The card's Verification field states (verbatim, line 11):
"SymPy audit: \StageFile{...stage210...sympy_audit.py}.  Mathematica audit: none yet."
But a Mathematica audit `.wl` exists, runs, and passes (the committed `.txt` ends `STAGE 210 MATHEMATICA AUDIT PASSED`). The `.wl` was added/strengthened in the pass-1 dual-engine retrofit; the card was not updated to cite it. The notes file does not mention scripts at all, so the `.tex` Verification line is the sole prose carrier and it is stale.

**Why this matters:**
The card under-reports the verification coverage of the stage; a reader auditing trust would conclude the Mathematica engine is absent when it is in fact present and independent. This is a paper-side staleness, not a math error, so direction routes to the user.

**Required change:**
None by Codex. Routes to user via `## Resolve before fix_loop` (the user owns paper/ edits). Expected resolution: update the Verification line to cite `mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl`.

**Verification:**
After user resolution, line 11 of the card names the `.wl` path.

### F2 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.txt` (mtime 2026-05-11 12:49)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py` (mtime 2026-06-03 15:59)

**What's wrong:**
The committed SymPy output predates the SymPy script by ~3 weeks AND its content disagrees with the current script: the saved `.txt` prints `STAGE 193 — THREE-COORDINATE MIXED-SIMPLEX...` and `STAGE 193 SYMPY AUDIT COMPLETED SUCCESSFULLY` (lines 11, 173), whereas the current `.py` banners read `STAGE 210` (lines 35, 256). This is the known P4-52 stale-banner set: the committed output predates the banner fix. The numeric/symbolic check results in the `.txt` (all `= 0`) are still consistent with the current script's checks; only the banner label is stale, plus the freshness gap. The Mathematica output (`.txt` mtime 2026-06-02 10:41) is NEWER than its `.wl` (2026-06-02 10:15) and shows the correct `STAGE 210` banner → fresh.

**Why this matters:**
Informational; the verifier will re-run SymPy and regenerate the output with the correct `STAGE 210` banner. The check substance is unaffected.

**Required change:**
Re-run `python3 scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py` and recommit the refreshed `.txt`. No script edit required.

**Verification:**
Refreshed `.txt` shows `STAGE 210` banner at top and `STAGE 210 SYMPY AUDIT COMPLETED SUCCESSFULLY` at bottom, mtime newer than the `.py`.

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.** The `.wl` is not a transliteration of the `.py`; it derives the load-bearing objects by methods the `.py` does not use. Three discriminating examples:

1. **Gradient-optimal point (deliverable 3).** The `.py` *posits* the answer and verifies it (lines 85-96):
   ```python
   avec_grad = sp.Matrix([ki / Kgrad, kj / Kgrad, kk / Kgrad])
   expect_zero("gradient-optimal ratio r", sp.simplify(avec_grad[1] / avec_grad[0] - kj / ki))
   ```
   The `.wl` *derives* it from a Lagrangian and solves (lines 118-140, 154-155):
   ```
   lagrangian = (x ki + y kj + z kk - lag (x^2 + y^2 + z^2 - 1));
   stationaritySolutions = Solve[stationarityEquations, {x, y, z, lag}, Reals];
   ...
   expectZero["M3 stationary branch minus gradient vector", stationaryVec - gradVec];
   ```
   Posit-and-verify vs solve-the-stationarity-system = independent derivation of the load-bearing vector.

2. **Discriminant coefficients A..F (deliverable 5).** The `.py` defines `Delta_sharp` as a hand-written polynomial and checks the whole expression simplifies to it (lines 192, 200-203). The `.wl` instead extracts each coefficient from a Taylor expansion (lines 233-243):
   ```
   deltaSeries = Expand[Normal[Series[deltaCleared, {r, 0, 2}, {s, 0, 2}]]];
   ...
   expectZero["M7 coefficient A residual", coeffRS[0, 0] - coefA];
   ```
   Whole-expression `simplify`-to-zero vs `Series`/`CoefficientList` coefficient extraction = different extraction operation on the load-bearing object.

3. **Barycenter uniqueness (deliverable 4).** The `.py` only checks the *value* w_Σ=2 at the barycenter (line 129). The `.wl` proves the maximizer is *unique* via `Reduce` (lines 181-203):
   ```
   uniqueEqualMixCondition = FullSimplify[Reduce[... normA2 == 1 && wSigma == 2, {ai, aj, ak}, Reals], ...];
   expectTrue["M6 equal-mix uniqueness condition", Equivalent[uniqueEqualMixCondition, ai == 1/Sqrt[3] && ...]];
   ```
   This is a strictly stronger, independently-derived claim than anything in the `.py`.

Additionally, the boundary reductions M2/M8 use `Limit[..., eps -> 0, Direction -> "FromAbove"]` (approaching the edge from the simplex interior with an `eps` perturbation), whereas the `.py` uses direct substitution `ak=0` via `.subs(...)`. Different operation, same physical premise. The shared monomial definitions (kVec, hMat, the unit-sphere constraint) are legitimately common; the method that extracts each load-bearing object differs. Not a port.

## Engine cross-check

Both engines reach the same conclusions and all residuals are `0` / all `expectTrue` are `True`. The SymPy output reports every Section I-VI residual as `= 0`; the Mathematica output reports every M1-M9 residual `= 0` with `PASS`, plus `M6 equal-mix uniqueness condition = True`. No sign, factor, or domain disagreement. `engines_agree: true`. (The stale SymPy banner does not affect the residual content; the checks themselves agree.)

## Verdict justification

`findings`. The algebraic substance of all five closed-form deliverables holds up against the paper under adversarial reading: I tried to break the gradient-optimal claim (the .wl independently solves the Lagrange system and lands on the same vector), the uniqueness of the barycenter (the .wl's `Reduce` confirms it is the sole maximizer), the discriminant coefficients (.wl extracts them by a different `Series`/`CoefficientList` route), the Lagrange/Cauchy identity (a genuine 3-vector identity, non-tautological), and the edge reductions (.wl approaches via interior limit, .py via direct substitution — both agree). No tautological, hardcoded, symbol-assumption, missing-branch, engine-disagreement, transliteration, or insufficient-verification finding survives. The two findings are both low-severity and non-blocking: F1 is a paper-side staleness (the card says "Mathematica audit: none yet" while a passing `.wl` exists — routes to user), and F2 is the known P4-52 stale SymPy output banner (`STAGE 193` in the committed `.txt` vs `STAGE 210` in the current `.py`). I confirm I read the card, the notes, and the appendix rows, and the scripts' algebraic claims match the paper.

## Value Reconciliation (pass-2 augmentation)

All Stage 210 deliverables are **symbolic** (no pinned floating-point constants); the only literals are the bound `w_Σ ≤ 2`, the edge value `w_Σ = 1`, and the barycenter `1/√3`. Each is a stated deliverable carried in the notes (boxed) and summarized in the card Output line.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Δ⁺_{ijk}` (unit-sphere simplex) | py 42-49 / wl 79-81; sympy out 59-61, math out 9-13 | card l.9 (Inputs), notes §1 boxed l.73-81 | MATCH |
| edge param `a_ij=(1,r,0)/√(1+r²)` (+ik,jk) | py 52-54 / wl 86-88; sympy out 17-58 | notes §2.1-2.3 boxed l.113-153 | MATCH |
| edge slope `(k_i+r k_j)/√(1+r²)` (+ik,jk) | py 67-78 / wl 111-113; math out 19-23 | notes §2.1-2.3 boxed | MATCH |
| `a_grad=(k_i,k_j,k_k)/√(Σk²)` | py 86 / wl 141,153; math out 29 | notes §3 boxed l.182-186 | MATCH |
| max slope `√(k_i²+k_j²+k_k²)` | py 85,247 / wl 117,161; sympy out 87-90 | notes §3 boxed l.189-194 | MATCH |
| interior ratios `r=k_j/k_i, s=k_k/k_i` | py 95-96 / wl 162-163; math out 44-46 | notes §3 boxed l.198-202 | MATCH |
| `w_Σ=2(a_i a_j+a_i a_k+a_j a_k)` | py 110 / wl 173; sympy out 102-103 | notes §4 boxed l.258-261 | MATCH |
| `w_Σ=(Σa)²−1` (on simplex) | py 116 / wl 191; sympy out 104 | notes §4 boxed l.265-269 | MATCH |
| `w_Σ ≤ 2` bound | py 129 / wl 194; sympy out 108 | notes §4 l.283 | MATCH |
| barycenter `a_eq=(1,1,1)/√3` | py 125 / wl 179; math out 72 | notes §4 boxed l.291-294 | MATCH |
| `w_Σ(a_eq)=2`, `w_Σ(pair)=1` | py 129-130 / wl 194-195; sympy out 108-109 | notes §4 boxed l.298-305 | MATCH |
| barycenter uniqueness `a_i=a_j=a_k=1/√3` | wl 196-203 (math out 72-74); not in py | notes §4 l.285-288 | MATCH |
| `H_1=aᵀHa` curvature law | py 148 / wl 212,247; sympy out 114-117 | notes §4 boxed l.241-253 | MATCH |
| τ root map `2H₀/(k+√(k²−2H₀κ))` | py 176 / wl 272; sympy out 126-143 | notes §5 boxed l.349-355 | MATCH |
| six coeffs `A..F` (Δ# numerator) | py 185-190 / wl 219-224; math out 79 | notes §5.1 boxed l.371-385 | MATCH |
| `Δ#=A+Br+Cs+Dr²+Ers+Fs²` | py 192 / wl 225; sympy out 144-149 | notes §5.1 boxed l.390-393 | MATCH |
| τ ratio-coord form `2H₀√(1+r²+s²)/(...)` | py 205-208 / wl 276-279; sympy out 152 | notes §5.1 boxed l.398-404 | MATCH |
| edge-τ reductions → Stage 209 pairwise | py 227-239 / wl 311-313; math out 110-115 | notes §5.1 l.406 | MATCH |

Internal scaffolding (no finding): pass/fail flags, `expect_zero`/`expectZero`/`expectTrue` residuals (all 0/True), the `eps` limit perturbation symbol, the Lagrange multiplier `lag`, the intermediate `slopeRS`/`kappaRS`/`deltaCleared`/`deltaSeries`/`deltaCoefficientList`, and `cauchyResidual`/`cauchySlack` identities.

reconciliation: complete; 18 deliverable values checked, 0 misaligned. The only paper-side discrepancy is non-numeric: the Verification field's "Mathematica audit: none yet" (F1), which is a coverage-statement staleness, not a value mismatch.

## Self-test notes

Checked: (1) no `sp.diff`/`D[expr,var]` derivative-vs-independent-variable trap — the only differentiation is `D[lagrangian, x/y/z]` in the .wl, and the Lagrangian genuinely depends on x,y,z, so the stationarity equations are non-trivial. (2) No unbounded-domain integrals (no parity concern). (3) Trivial-case pre-check: substituting `k_i=k_j=k_k` gives `a_grad=(1,1,1)/√3` (the barycenter coincides with the gradient point) and `w_Σ=2`, consistent; the Cauchy/Lagrange identities reduce to `0` as required. (4) No missing-script directive (both engines present), so no path-spec trap. (5) Both findings route to user (F1 paper_misalignment) or to a re-run (F2 stale_output); neither prescribes a Codex script edit, so there is no paper round-trip risk.
