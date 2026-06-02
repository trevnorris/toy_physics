---
unit_id: 210
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.md]
  paper_appendix: present
---

# Audit unit 210 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_210.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (read intro rows + table row 210 at line 51 + the three-coordinate simplex sections lines 920–1091)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The stage card (`\stagefield{Output}`) states verbatim: "Boundary reduction for triples, three-coordinate gradient and curvature synergy, canonical triple-screen audit, interior-screen dominance criterion, and non-improvement filter." The purpose is to "promote the search to genuine three-coordinate positive simplices and prove that their boundaries are already the Stage 209 pairwise ledgers." The notes and the Part VI appendix (lines 923–986) make the algebra explicit: (1) the positive spherical simplex `Δ⁺_ijk = {a∈ℝ³_{≥0} : a_i²+a_j²+a_k²=1}` with oriented ray `ŝ = a_i ê_i + a_j ê_j + a_k ê_k`; (2) each of the three boundary edges reduces exactly to a lower-support pairwise cone, with slope `k_simplex = a·k → (k_i + r k_j)/√(1+r²)` etc.; (3) the gradient-optimal interior point `a^grad = (k_i,k_j,k_k)/√(k_i²+k_j²+k_k²)` with max slope `√(Σk²)`, strictly exceeding every pairwise edge max; (4) the curvature law `H₁ = aᵀHa` with off-diagonal weight `w_Σ = 2(a_ia_j+a_ia_k+a_ja_k) = (Σa)²−1 ≤ 2`, maximized (=2) uniquely at the equal-mix barycenter `(1,1,1)/√3`, vs `w_Σ=1` on a pairwise equal edge; (5) the fixed-simplex certified root `τ = 2H₀/(k_simplex + √(k_simplex²−2H₀κ))` and its interior ratio-coordinate form with quadratic discriminant numerator `Δ^♯(r,s) = A + Br + Cs + Dr² + Ers + Fs²`, `A=k_i²−2H₀u_ii`, `B=2k_ik_j−4H₀u_ij`, …; (6) the five-row canonical triple-screen packet; (7) the interior-screen dominance criterion and non-improvement filter (these last two are inequality-comparison logic, not closed-form identities).

## What the script claims to verify

The SymPy script verifies the exact symbolic algebra underpinning deliverables (1)–(5). It checks: edge parametrization normalization and the three edge slope reductions to the pairwise forms (Sec I); gradient-optimal normalization, slope value `=√(Σk²)`, ratios `r=k_j/k_i, s=k_k/k_i`, and the Pythagorean relations `Kgrad²−Kij²−k_k²=0` etc. (Sec II); the cross-leverage identity `w_Σ−((Σa)²−||a||²)=0`, the Cauchy slack identity, and the `w_Σ=2` (barycenter) vs `w_Σ=1` (edge) values (Sec III); the quadratic-form curvature scalar `κ=aᵀHa`, its three edge reductions, and the diagonal-neutral reduction (Sec IV); the fixed-simplex root map `τ`, the ratio-coordinate discriminant reduction `(1+r²+s²)(k_rs²−2H₀κ_rs) − Δ^♯ = 0`, the ratio-coordinate `τ` form, and the three boundary reductions of `τ` to the pairwise forms (Sec V); and the equal-mix/gradient slope values (Sec VI). The script does not (and need not — these are comparison inequalities, not identities) numerically encode deliverables (6)/(7).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Positive spherical simplex + oriented ray | edge parametrizations + normalizations (lines 52–65) | match |
| Boundary reduction: each edge = pairwise cone (slope) | edge ij/ik/jk slope reductions (lines 67–78) | match |
| Gradient-optimal point + max slope `√(Σk²)` | a_grad normalization, slope value, ratios (lines 93–96, 252–254) | match |
| Interior strictly beats pairwise edges (first order) | Pythagorean `Kgrad²−Kij²−k_k²=0` family (lines 101–103) | match (encodes the strict-gain decomposition) |
| Curvature law `H₁=aᵀHa`, weight `w_Σ` identity, ≤2 | `w_Σ` identity, Cauchy slack, edge curvature reductions (lines 114–169) | match |
| Equal-mix barycenter maximizes `w_Σ`(=2); edge=1 | `w_Σ(eq)−2`, `w_Σ(edge)−1`, diagonal-neutral (lines 127–130, 166–169) | match |
| Fixed-simplex root `τ` + ratio-coord discriminant `Δ^♯` | `τ` map, discriminant reduction, ratio τ form, boundary reductions (lines 176–239) | match |
| Canonical 5-row triple-screen packet | (none) | not a closed-form identity; structural/definitional, no script check needed |
| Interior-screen dominance + non-improvement filter | (none) | inequality comparison logic; not symbolically asserted, acceptable |
| Mathematica independent verification | (none) | missing — see F1 |

`paper_alignment: aligned` — every closed-form identity the paper states is faithfully exercised by a non-tautological SymPy assertion; the only items without a script check (deliverables 6/7) are definitional/inequality content, not algebraic identities, so their absence is not a `script_missing_paper_claim` defect. No numeric constants disagree; no extra (unmotivated) script content. The only gap is the absent second engine.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 63–65 | `expect_zero(aᵀa − 1)` for 3 edges | simplex normalization | yes |
| A2 | sympy | 67–78 | `expect_zero(a·k|_edge − (k_i+rk_j)/√(1+r²))` ×3 | boundary slope reduction | yes |
| A3 | sympy | 93 | `expect_zero(a_grad·a_grad − 1)` | grad point normalization | yes |
| A4 | sympy | 94 | `expect_zero(a_grad·k − √(Σk²))` | gradient-optimal max slope | yes |
| A5 | sympy | 95–96 | `expect_zero(ratio − k_j/k_i)`, `(− k_k/k_i)` | interior optimal ratios | yes |
| A6 | sympy | 101–103 | `expect_zero(Kgrad² − Kij² − k_k²)` ×3 | strict interior-over-edge gain | yes |
| A7 | sympy | 114–117 | `expect_zero(w_Σ − ((Σa)²−||a||²))` | cross-leverage identity | yes |
| A8 | sympy | 119–123 | `expect_zero(3||a||² − (Σa)² − Σ(a_i−a_j)²)` | Cauchy slack (≤2 proof) | yes |
| A9 | sympy | 127–130 | `expect_zero(w_Σ(eq)−2)`, `(edge)−1`, norms | barycenter maximizes leverage | yes |
| A10 | sympy | 151–169 | `expect_zero(κ|_edge − edge form)` ×3 + diagonal-neutral | curvature law + edge reduction | yes |
| A11 | sympy | 196–203 | `expect_zero((1+r²+s²)(k_rs²−2H₀κ_rs) − Δ^♯)` | discriminant numerator (load-bearing) | yes |
| A12 | sympy | 209–212 | `expect_zero(τ|_rs − τ_rs_expected)` | ratio-coordinate root form | yes |
| A13 | sympy | 227–239 | `expect_zero(τ_rs|_{s=0} − pairwise ij)` ×3 | τ boundary reduction to pairwise | yes |
| A14 | sympy | 252–254 | `expect_zero(eq/grad slope values)`, `k_grad²−||k||²` | screen slope values | yes |
| (none) | mathematica | — | — | all of the above | MISSING |

All SymPy rows are non-tautological: each defines an object one way and checks it against an independently-written closed form (e.g. `w_Σ` defined as `2(a_ia_j+…)` and checked against `(Σa)²−||a||²`; `Δ^♯` written term-by-term and checked against the substituted-and-cleared `(1+r²+s²)(k²−2H₀κ)`). None is a `x=expr; assert x==expr` self-check.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `.wl` for unit 210)

**What's wrong:**
Stage 210 is non-status-only (`is_status_only_candidate: False`) and non-checkpoint (`is_checkpoint: False`). The paper card's `\stagefield{Verification}` line itself records the gap: "Mathematica audit: none yet." The stage performs substantial, fully verifiable symbolic algebra — matrix quadratic forms `aᵀHa`, a square-root root map `τ = 2H₀/(k + √(k²−2H₀κ))`, a polynomial discriminant clearing of square roots, Lagrange/Cauchy-optimal points, and edge-limit reductions. Mathematica can independently verify every one of these from the physical premises. The project dual-engine contract therefore requires a `.wl`; the existing passing SymPy script does not discharge that requirement.

**Why this matters:**
The whole point of the second engine is to catch an algebra error or a SymPy-specific simplification artifact that a single engine would silently pass. Several of this stage's checks rely on `sp.simplify(sp.expand(...))` over radical expressions (the `τ` map and the discriminant clearing in Sec V); a second engine using a different decomposition is exactly the cross-check designed to surface a hidden branch/sign issue in those radicals. Without it the stage is single-engine.

**Required change:**
Add an independent Mathematica script `mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.wl` verifying the claim manifest in the directive. It must derive each result independently (native Mathematica primitives, a different decomposition than the `.py`), not transliterate the SymPy algebra.

**Verification:**
The verifier runs the new `.wl` (via `redteam exec-mathematica 210` / `math -script`), confirms it exits 0 with all in-file checks passing, and confirms (per the independent-derivation check) it is not a line-by-line port of the `.py`.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. The directive's anti-transliteration guard requires Codex to use a different decomposition (Series/Coefficient extraction of `Δ^♯`, eigen/matrix primitives for the quadratic form, FindRoot or Reduce for the τ map limits) rather than echoing the SymPy `subs` choreography.

## Engine cross-check

Only one engine present; not applicable. SymPy output (mtime 2026-05-11 12:49) is newer than the script (mtime 2026-04-14 08:42), so the saved transcript is fresh and its 30+ `= 0` lines plus `EXIT_CODE: 0` reflect the current script. Note: the SymPy banner and output header read "STAGE 193" while the unit is paper-stage 210 — this is the known internal-renumbering cosmetic label (notes call it "Stage 244"); it carries no math and is out of red-team scope (the red-team does not edit prose/labels for paper-alignment, and the banner string is not load-bearing). Recorded here for visibility only.

## Verdict justification

The SymPy script is solid: every closed-form deliverable the paper states (boundary reduction, gradient synergy, curvature/cross-leverage law, fixed-simplex root map and its ratio-coordinate discriminant, and the screen slope values) is exercised by a non-tautological, well-anchored assertion, and all constants and forms match the appendix verbatim. I tried to break it on (a) tautology — each object is checked against an independently written closed form, not against itself; (b) the radical simplifications in Sec V — the discriminant-clearing identity `(1+r²+s²)(k²−2H₀κ) = Δ^♯` is a genuine polynomial identity that holds for all `(r,s)`, not a limiting-case artifact; (c) symbol domains — `k_i,k_j,k_k>0`, `a_i,a_j,a_k≥0`, `H₀>0`, `u_••` real, all consistent with the physical setup. No paper misalignment. The single defect is the absent second engine, which the project contract requires. Verdict: `findings` (one `missing_mathematica`), not stop_cold.

## Self-test notes

Variable-independence trap: N/A — the script contains no `sp.diff`/`D[]`; all checks are static algebraic identities, so no derivative can collapse to identically zero. Parity/integral trap: N/A — no integrals over unbounded domains. Trivial-case pre-check: I mentally substituted `a_grad=k/|k|` into A4 (`a·k=|k|`✓), `a=(1,1,1)/√3` into A9 (`w_Σ=2·3·(1/3)=2`✓ and `2(1/3+1/3+1/3)=2`✓), and `s=0` into A13 (τ reduces to the pairwise-ij form✓). Path spec for F1: target is `mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit.wl` (the `.wl` directory), stated in the directive. The fix prescribes no new constant, so no new paper_misalignment is introduced.
