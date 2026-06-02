---
unit_id: 213
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md]
  paper_appendix: present
---

# Audit unit 213 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_213.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows at lines 57, 236, 1029, 1307)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The paper card's `\stagefield{Output}` states verbatim: "Four-coordinate simplex gate, face-reduction theorem, four-way gradient/equal-mix screens, support-cardinality-four improvement/non-improvement filters." The card's `\stagefield{Derivation ledger}` enumerates: "proves exact face reduction, derives the four-coordinate gradient-synergy law, computes six-way cross leverage, constructs fixed-simplex brackets, and states the support-cardinality-four improvement gate against the imported boundary ledger." The notes file is far more detailed and is authoritative on intent; it lists seven exact deliverables (Section "Purpose"): (1) the combinatorial ledger for the five primitive quadruples, (2) the positive spherical four-simplex and its face reduction to the imported triple packets, (3) the gradient-synergy theorem with unique interior gradient-optimal ray `a_grad = (k_i,k_j,k_k,k_l)/||k||` and `max k = ||k||`, (4) the curvature law and the theorem that the equal-mix barycenter `(1/2,1/2,1/2,1/2)` maximizes total six-way off-diagonal leverage at `w_Σ = 3` (vs 2 on a triple face, 1 on a pair edge), (5) the fixed-simplex certified bracket `τ = 2H₀/(k + √(k²−2H₀κ))` with the ten discriminant coefficients A…J and the ratio-patch form, (6) the canonical quadruple-screen audit (four imported faces + grad row + eq row), and (7) the support-cardinality-4 theorem gate (interior-screen dominance, non-improvement filter, and the global gate against the imported support-≤3 ledger). The appendix row (part06:57) summarizes the stage as "Reduces four-simplex faces to triple packets and introduces four-way gradient/equal-mix canonical screens."

## What the script claims to verify

The single SymPy script verifies, in six sections, exactly the seven note deliverables: (I) the combinatorial ledger — #quadruples = C(5,4) = 5, each triple in exactly 2 quadruples, each axis in exactly 4; (II) the four face parametrizations, their unit-norm, and the slope reductions `k_simplex|_face = (k_i + r k_j + s k_k)/√(1+r²+s²)` etc.; (III) `a_grad` normalization, its slope value `= ||k||`, the three optimal ratios `k_j/k_i, k_k/k_i, k_l/k_i`, and the four synergy gaps `||k||² − ||k_face||² = (dropped k)²`; (IV) `w_Σ = (Σa)²−||a||²`, the Cauchy slack identity, and `w_Σ = 3/2/1` at the four-way/triple/pair equal-mix points; (V) the Hessian quadratic form `κ = aᵀHa`, the diagonal-neutral reduction, the ratio-patch parametrization, the discriminant numerator `Δ♯` built from the ten coefficients A…J, the `τ` bracket form, and four face-reduction collapses (t=0→ijk, s=0→ijl, r=0→ikl, and the jkl face); (VI) `k_eq = (Σk)/2`, the gradient slope, and four finite integer-sampling checks of the interval-ordering gate theorems (boundary splice, screen dominance, non-improvement filter, support-4 improvement/non-improvement).

## Paper ↔ script cross-check

| Paper deliverable (notes §) | Script-side check | Verdict |
|---|---|---|
| (1) 5 primitive quadruples, triple↔2, axis↔4 | I (lines 59–71) | match |
| (2) face reductions to triple packets | II (lines 101–125) | match |
| (3) gradient-synergy: a_grad, max=‖k‖, ratios, synergy gaps | III (lines 140–153) | match |
| (4) w_Σ identity, Cauchy slack, barycenter maximizer (3/2/1) | IV (lines 166–192) | match |
| (5) curvature law κ=aᵀHa, Δ♯ coeffs A…J, τ bracket, face collapses | V (lines 211–302) | match |
| (6) canonical screen packet (grad + eq rows; k_eq=(Σk)/2) | VI (lines 309–311) | match |
| (7) interval-gate theorems (dominance, filters, global support-4 gate) | VI (lines 318–400) | match (weak; see note) |
| second-engine (Mathematica) independent re-derivation | (none) | missing |

`paper_alignment: aligned` — every paper-side deliverable has a faithful SymPy counterpart with the correct constants and forms. The only gap is the absent second engine.

Note on (7): the four nested-loop integer-sampling checks (lines 327–400) exercise interval-ordering implications that are essentially true by the loop's own construction (e.g. line 354 gates on `can_hi < blo` and then asserts `c_star < b_star` for `b_star ∈ [blo,bhi]`, `c_star ∈ [can_lo,can_hi]` — guaranteed since `c_star ≤ can_hi < blo ≤ b_star`). These are weak (they re-confirm transitivity of `<` rather than any physics), but they are not false and they do mirror the paper's stated gate theorems; the same weakness would be present in any engine. I do not raise this as a separate finding because the load-bearing physics of the stage lives in Sections II–V (the algebraic identities), which the script exercises non-tautologically, and the gate theorems are genuinely just interval comparisons in the paper too.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `#quadruples − C(5,4) == 0` | (1) | yes |
| A2 | sympy | 64–65 | each triple in exactly 2 quads | (1) | yes |
| A3 | sympy | 69–71 | each axis in exactly 4 quads | (1) | yes |
| A4 | sympy | 101–104 | face unit-norm (4 faces) | (2) | yes |
| A5 | sympy | 106–125 | face slope reductions (4 faces) | (2) | yes |
| A6 | sympy | 140 | a_grad normalization | (3) | yes |
| A7 | sympy | 141 | a_gradᵀk − ‖k‖ == 0 | (3) | yes |
| A8 | sympy | 142–144 | grad ratios k_j/k_i, k_k/k_i, k_l/k_i | (3) | yes |
| A9 | sympy | 150–153 | synergy gaps ‖k‖²−‖k_face‖² = dropped² | (3) | yes |
| A10 | sympy | 166–169 | w_Σ = (Σa)²−‖a‖² | (4) | yes |
| A11 | sympy | 170–182 | Cauchy slack identity | (4) | yes |
| A12 | sympy | 187–192 | w_Σ = 3 / 2 / 1 at eq points | (4) | yes |
| A13 | sympy | 214–218 | diagonal-neutral κ reduction | (5) | yes |
| A14 | sympy | 248–252 | Δ♯ discriminant numerator reduction | (5) | yes |
| A15 | sympy | 259–262 | τ ratio-patch bracket form | (5) | yes |
| A16 | sympy | 265–302 | τ face collapses (ijk/ijl/ikl/jkl) | (5) | yes |
| A17 | sympy | 310 | k_eq = (Σk)/2 | (6) | yes |
| A18 | sympy | 311–312 | grad slope, ‖k‖² | (3)/(6) | yes |
| A19 | sympy | 327–400 | interval-gate integer sampling (×4) | (7) | partial (weak; tautological ordering) |

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 213 is non-status-only (`is_status_only_candidate: False`) and non-checkpoint (`is_checkpoint: False`), yet it has only a SymPy engine. The paper card's `\stagefield{Verification}` line itself flags this: "Mathematica audit: none yet." Per the project dual-engine contract and the rendered prompt (line ~118), both engines are required wherever Mathematica *can* independently verify the math. Stage 213's deliverables are entirely closed-form symbolic identities — symmetric `4×4` quadratic forms `aᵀHa`, surd-normalized face parametrizations, a degree-2 discriminant polynomial in `(r,s,t)`, the `w_Σ = (Σa)²−‖a‖²` rewrite, the Cauchy slack sum-of-squares identity, and a per-point `τ` bracket — all of which Mathematica can verify natively (`FullSimplify`, `Solve`, `Series`/`Coefficient`, matrix products, `Maximize`/`Reduce`). There is no impossibility carve-out: Mathematica can re-derive every load-bearing identity here. The absence is therefore a finding, not a tolerated single-engine case.

**Why this matters:**
The second engine is the cross-check that catches SymPy-specific simplification artifacts and transcription errors in the ten discriminant coefficients and the τ bracket. Without it, Section V's `Δ♯` reduction and the τ face-collapse identities rest on a single CAS. The dual-engine policy exists precisely so an algebra-level slip in one engine is caught by an independent decomposition in the other.

**Required change:**
Create the missing Mathematica audit at the target path with an independent re-derivation (claim manifest M1–M9 in the directive). Must use native Mathematica primitives via a different decomposition than the `.py`; a line-by-line port is rejected as transliteration.

**Verification:**
After Codex applies, `mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl` exists, runs under `math -script`, exits 0, and each claim M1–M9 prints a zero residual / confirmed-True line.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed yet. The directive's claim manifest specifies an anti-transliteration guard requiring a different decomposition (e.g., `Maximize`/Lagrange for the gradient optimum and the `w_Σ ≤ 3` bound rather than re-asserting the closed forms; `Coefficient` extraction on the expanded `(1+r²+s²+t²)(k²−2H₀κ)` for the discriminant rather than re-typing A…J).

## Engine cross-check

Only one engine present; no cross-check possible. The SymPy output (`...sympy_audit.txt`) is fresh: output mtime `2026-05-11 12:49:27` is later than script mtime `2026-05-11 11:58:52`, and every printed residual is `0` / `True` with `EXIT_CODE: 0`. No `stale_output` finding.

## Verdict justification

The SymPy script is a faithful, non-tautological verification of all seven note deliverables with the correct constants (C(5,4)=5; incidences 2 and 4; `a_grad=(k_i..k_l)/‖k‖`; `max k = ‖k‖`; synergy gaps equal the squared dropped slope; `w_Σ` bound 3 at the barycenter, 2/1 on lower faces; the ten A…J coefficients matching the boxed paper forms; the τ bracket and its four face collapses). I attacked the algebra (face-reduction surds, the discriminant numerator `(1+r²+s²+t²)(k²−2H₀κ)−Δ♯`, the Cauchy slack identity, the diagonal-neutral reduction) and the constants against the paper — all hold. The integer-sampling gate checks in Section VI are weak (they re-confirm interval transitivity) but are not false and faithfully mirror the paper's interval theorems, so they are noted, not flagged. The only defect is the absent second engine, which the dual-engine contract requires here because Mathematica demonstrably *can* verify this stage. Verdict: `findings` (one `missing_verification_script` / `missing_mathematica`), no stop-cold, paper alignment exact.

## Self-test notes

I checked the self-test traps against the claim manifest I am prescribing. (1) Variable independence: the manifest contains no `D[expr,var]` derivative checks (the gradient optimum is verified via constrained `Maximize`/Lagrange, not by asserting a derivative vanishes), so the identically-zero-derivative trap does not arise; the `w_Σ ≤ 3` claim is a bounded `Maximize` over the simplex, whose maximizer `(1/2,1/2,1/2,1/2)` I confirmed gives `w_Σ=3`. (2) Symmetry/parity: no unbounded-domain integrals appear; all integration-free. (3) Trivial-case pre-check: substituting `k_i=k_j=k_k=k_l=1` gives `a_grad=(1/2)(1,1,1,1)` with slope `‖k‖=2` and `w_Σ=3`, and the discriminant identity reduces to `0` (both sides expand identically), so each `expectZero` target is genuinely 0 and each confirmed value is the stated literal — none are trivially-true for the wrong reason. (4) Path: the target `.wl` lives under `mathematica/` with the `_mathematica_audit.wl` suffix matching the stage218 sibling convention. (5) Paper round-trip: the manifest uses only constants the paper states (C(5,4), 2/4 incidences, ‖k‖, 3/2/1, the boxed A…J coefficients, the τ bracket), introducing no new constant — no new paper_misalignment.
