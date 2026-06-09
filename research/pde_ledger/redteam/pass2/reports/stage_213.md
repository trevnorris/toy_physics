---
unit_id: 213
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md]
  paper_appendix: present
---

# Audit unit 213 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_213.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (row 57; narrative line 236)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.txt`

## What the paper claims

Stage 213 opens the four-coordinate search. Its `\stagefield{Output}` is: "Four-coordinate simplex gate, face-reduction theorem, four-way gradient/equal-mix screens, support-cardinality-four improvement/non-improvement filters." The notes enumerate seven deliverables: (1) the five-element primitive quadruple ledger `#𝔔₄=binom(5,4)=5` with each triple in exactly two quadruples; (2) the positive spherical four-simplex `Δ⁺_{ijkl}` and the exact reduction of each of its four codimension-one faces to a Stage-212 closed triple simplex; (3) the four-coordinate gradient-synergy theorem with optimal ray `a_grad=(k_i,k_j,k_k,k_l)/||k||`, max slope `||k||=sqrt(k_i²+k_j²+k_k²+k_l²)`, optimal ratios `k_j/k_i, k_k/k_i, k_l/k_i`, and the per-face synergy gaps `||k||²-||k_face||²=k_missing²>0`; (4) the curvature law with total six-way cross-leverage weight `w_Σ=2Σ_{p<q}a_p a_q=(Σa)²-1`, the Cauchy slack identity, and the equal-mix barycenter `(1,1,1,1)/2` uniquely maximizing `w_Σ=3` (vs 2 on a triple face, 1 on a pair edge); (5) the fixed-simplex certified tau bracket with the ten discriminant coefficients `A_⋆..J_⋆` and the ratio-coordinate form; (6) the canonical six-row quadruple-screen; (7) the support-cardinality-4 improvement/non-improvement interval gates against the Stage-212 boundary and up-to-three ledgers. Note: the card's `\stagefield{Verification}` line states "Mathematica audit: none yet", which is now false (see F2).

## What the script claims to verify

The SymPy script verifies, section by section: I — the combinatorial quadruple ledger (count=5, triple incidence 2, axis incidence 4); II — the four face vectors are unit-normalized and their simplex slope `a·k` reduces to the posited triple slope; III — the posited gradient ray is normalized, gives slope `||k||`, gives the three posited ratios, and the four synergy-gap identities; IV — the `w_Σ` identity, the Cauchy slack identity, and `w_Σ` evaluated at the four-way/triple/pair equal-mix points giving 3/2/1; V — the diagonal-neutral curvature reduction, the ten posited discriminant coefficients reproduce the expanded numerator `(1+r²+s²+t²)(k²-2H0·κ)`, the ratio-form tau bracket, and four face collapses; VI — the equal-mix and gradient slope values, and four interval-arithmetic gate theorems (boundary splice, screen dominance, non-improvement filter, support-4 improvement gate, support-4 non-improvement filter) verified by brute-force integer sampling. The Mathematica script verifies the same physics but DERIVES the load-bearing objects independently (Solve for the gradient ray, Maximize for the leverage bound, Coefficient-extraction for the discriminant coefficients) — see the independent-derivation section.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) quadruple ledger `#=5`, incidences | py I / wl M1 | match |
| (2) face reduction to Stage-212 triples | py II (norm + slope) / wl M2 | match |
| (3) gradient ray, `||k||`, ratios, synergy gaps | py III (posited+verified) / wl M3 (Solve-derived) + M4 | match |
| (4) `w_Σ` identity, Cauchy slack, barycenter max=3, 3/2/1 | py IV / wl M5 (Maximize-derived) | match |
| (5) ten discriminant coeffs, ratio tau bracket, face collapses | py V / wl M6+M7 (Coefficient-extracted) + M8 + M9 | match |
| (6) canonical six-row screen `S_quad` | implicit — gradient + equal-mix screen points constructed; the packet tuple itself is a labeling, exercised via its two interior points | partial (see note) |
| (7) support-4 improvement/non-improvement gates | py VI brute-force interval samples (924 each) | match |

Note on (6): `S_quad` is a tuple bundling the four imported face intervals plus the two interior screen points. The script constructs and exercises both interior screen points (gradient + equal-mix) and the four-face interval Min in beta_lo/beta_hi; the tuple itself is just a bookkeeping container, so "partial" here is a labeling artifact, not a substantive gap — every component is checked. No finding.

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `#quad - binom(5,4) == 0` | (1) | yes |
| A2 | sympy | 64-65 | triple incidence == 2 (loop raise) | (1) | yes |
| A3 | sympy | 70-71 | axis incidence == 4 (loop raise) | (1) | yes |
| A4 | sympy | 101-104 | face unit-norm == 0 | (2) | yes |
| A5 | sympy | 106-125 | face slope reduction == 0 | (2) | yes |
| A6 | sympy | 140-144 | grad ray norm/slope/ratios == 0 | (3) | yes (posited ray) |
| A7 | sympy | 150-153 | synergy-gap identities == 0 | (3) | yes |
| A8 | sympy | 166-182 | w_Σ identity + Cauchy slack == 0 | (4) | yes |
| A9 | sympy | 187-192 | equal-mix w_Σ = 3/2/1 == 0 | (4) | yes |
| A10 | sympy | 214-218 | diagonal-neutral κ reduction == 0 | (5) | yes |
| A11 | sympy | 248-252 | discriminant numerator reduction == 0 | (5) | yes |
| A12 | sympy | 259-302 | ratio tau form + 4 face collapses == 0 | (5) | yes |
| A13 | sympy | 327-400 | 4 interval gate theorems (brute-force) | (7) | yes |
| M1 | math | 111-113 | quad count / incidences | (1) | yes |
| M2 | math | 140-146 | face norm + slope reduction | (2) | yes |
| M3 | math | 158-189 | Solve-derived Lagrange ray + ratios + dominance | (3) | yes (DERIVED) |
| M4 | math | 199-202 | synergy-gap identities | (3) | yes |
| M5 | math | 210-222 | Maximize-derived w_Σ max=3 + barycenter | (4) | yes (DERIVED) |
| M6 | math | 228 | diagonal-neutral κ reduction | (5) | yes |
| M7 | math | 251-284 | Coefficient-extracted disc coeffs vs matrix form | (5) | yes (EXTRACTED) |
| M8 | math | 291-295 | tau ratio-form residual | (5) | yes |
| M9 | math | 299-351 | face Δ polynomials + tau face collapses | (5) | yes |

All assertions trace to a paper deliverable; none are orphaned. The interval-gate theorems (A13) are pure interval arithmetic exhaustively sampled — they correctly exercise the gate logic the paper states.

## Findings

### F1 — paper_misalignment

**Subtype:** value-prose mismatch on verification coverage (routes to user)
**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_213.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl` (exists, runs clean)

**What's wrong:**
The card states (line 11): "SymPy audit: \StageFile{...sympy_audit.py}.  Mathematica audit: none yet." A Mathematica audit `.wl` is present, runs clean (output dated 2026-06-02), and independently derives the load-bearing objects (M3 Solve, M5 Maximize, M7 Coefficient). The card prose is stale relative to the now-present second engine.

**Why this matters:**
The card understates verification coverage; a reader is told there is no Mathematica audit when there is a substantive, independent one. This is a prose/coverage discrepancy, not a math error. It routes to the user because Codex must not edit paper/.

**Required change:**
See `## Resolve before fix_loop` in the directive — user should update the card's `\stagefield{Verification}` line to cite the present `.wl`. No script change.

**Verification:**
Card line 11 cites the Mathematica audit path once updated.

### F2 — stale_output

**Severity:** low (informational)
**Files:**
- SymPy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.txt` (mtime 2026-05-11 12:49)
- SymPy script: `.../scripts/moving_throat_pde_stage213_..._sympy_audit.py` (mtime 2026-06-03 15:59)

**What's wrong:**
The saved SymPy `.txt` (2026-05-11) predates the current `.py` (2026-06-03). The captured output also carries the pre-renumber banner "STAGE 196 — ..." (lines 11, 236) and "STAGE 196 SYMPY AUDIT COMPLETED SUCCESSFULLY", whereas the current `.py` source prints "STAGE 213" (line 42). So the saved transcript reflects an older revision (196 + 17 = 213, the known +17 renumber offset). The Mathematica output (2026-06-02 10:28) is fresh relative to its `.wl` (2026-06-02 10:19).

**Why this matters:**
The committed SymPy transcript no longer matches what the current script prints (banner text). The mathematical content of the checks is unchanged — every check still reduces to 0 / True — so this is informational; the verifier's fresh re-run will regenerate the transcript with the correct banner.

**Required change:**
None for Codex. The orchestrator's independent re-run refreshes `scripts/output/...sympy_audit.txt`.

**Verification:**
After a fresh run, the SymPy `.txt` banner reads "STAGE 213" and mtime is newer than the `.py`.

## Independent-derivation check (Mathematica)

The `.wl` is **INDEPENDENT**, not a port, on all three load-bearing objects. The discriminator is derive-vs-posit:

1. **Gradient-optimal ray.** SymPy POSITS the answer and verifies properties:
   ```python
   Kgrad = sp.sqrt(ki**2 + kj**2 + kk**2 + kl**2)
   avec_grad = sp.Matrix([ki / Kgrad, kj / Kgrad, kk / Kgrad, kl / Kgrad])   # posited
   expect_zero("gradient-optimal slope value", (avec_grad.T * kvec)[0] - Kgrad)
   ```
   Mathematica DERIVES it from Lagrange stationarity:
   ```wolfram
   lagGradient = D[lagObjective - mu (lagNorm - 1), #] & /@ lagVars;
   lagSolutions = Solve[lagSystem, Append[lagVars, mu], Reals];
   positiveBranch = SelectFirst[lagSolutions, TrueQ[FullSimplify[(mu /. #) > 0 ...]] &];
   ```
   `Solve` of the KKT system vs a posited literal = different method extracting the load-bearing object.

2. **Cross-leverage maximum.** SymPy posits the barycenter `aeq4=(1/2,1/2,1/2,1/2)` and evaluates `wSigma` there, arguing the bound 3 via the Cauchy slack identity. Mathematica DERIVES both the value and the maximizer:
   ```wolfram
   maxLeverage = Maximize[{wSigma[lagVars], lagVars.lagVars == 1 && And @@ Thread[lagVars >= 0]}, lagVars, Reals];
   expectZero["M5 constrained maximum value - 3", First[maxLeverage] - 3];
   expectZero["M5 maximizer equals four-way equal mix", (lagVars /. Last[maxLeverage]) - {1/2, 1/2, 1/2, 1/2}];
   ```
   `Maximize` over the constrained nonneg sphere returns `{3, {1/2,1/2,1/2,1/2}}` (output line 101) — the value and maximizer are computed, not posited.

3. **Discriminant coefficients.** SymPy posits each coefficient (`Acoef = ki**2 - 2*H0*uii`, ...) and checks the assembled `Delta_sharp` against the expanded numerator. Mathematica EXTRACTS each coefficient by monomial projection from the independently-built numerator and compares it to a matrix-form kernel:
   ```wolfram
   discNumerator = normalizeExpr[normRst (kRst^2 - 2 H0 kappaRst), $Assumptions];
   extractedCoeffs = ... coeff3[discNumerator, #2] ...;   (* Coefficient[Coefficient[Coefficient[...]]] *)
   quadraticKernel = Outer[Times, kVec, kVec] - 2 H0 Hblock;   (* matrix-form, not the py literals *)
   expectZero["M7 coefficient " <> name, extractedCoeffs[name] - matrixCoeffs[name]];
   ```
   `Coefficient`-extraction + matrix-kernel comparison vs posited scalar literals = different method.

The shared physical premises (the monomial axes, `k·a` slope, the Hessian quadratic form `aᵀHa`) are legitimately common — that is the same physics, not a port. What differs is the METHOD that extracts each load-bearing object. The "each CAS runs its own simplifier" defense is not invoked here; the methods are genuinely different operations (Solve/Maximize/Coefficient vs posit-and-verify). This is the same independence pattern that distinguishes this stage from sibling 211 (which was a diff-vs-diff port). Verdict: INDEPENDENT.

## Engine cross-check

Both engines agree. The combinatorial ledger (5 quadruples, incidence 2/4), face reductions, gradient ray `(k_i,...,k_l)/||k||` with `||k||=sqrt(Σk²)`, ratios `k_j/k_i` etc., synergy gaps `k_missing²`, `w_Σ` max 3 at the barycenter (and 2/1 at triple/pair), diagonal-neutral reduction, the ten discriminant coefficients A..J, the ratio tau bracket, and all face collapses are confirmed zero/true in BOTH transcripts. Mathematica `Maximize` returns exactly `{3, {1/2,1/2,1/2,1/2}}` (line 101), matching the SymPy posited barycenter value 3. No residual disagreement, no sign/factor mismatch.

## Verdict justification

The math holds up under attack. I tried to break it on the batch-VI.1 transliteration concern and failed: the `.wl` genuinely re-derives the gradient ray (Solve), the leverage bound (Maximize), and the discriminant coefficients (Coefficient-extraction) rather than echoing the SymPy posit-and-verify choreography — these are different operations on the load-bearing objects, so it is INDEPENDENT, not a port. I checked the gradient-ray normalization and synergy gaps (all reduce correctly to `k_missing²>0`), the Cauchy slack identity (algebraically exact), the discriminant reduction (the `(1+r²+s²+t²)(k²-2H0κ)` clearing matches the posited polynomial), and the interval-gate brute-force sampling (correct interval logic). The verdict is `findings` only because of (F1) a low-severity stale card-prose line ("Mathematica audit: none yet" is now false — routes to user) and (F2) an informational stale SymPy transcript (pre-renumber "STAGE 196" banner; content unchanged, refreshed by the verifier's re-run). No script-side math defect, no tautology, no insufficient verification, no engine disagreement. I confirm I read the card, notes, and appendix row, and the script's claim matches the paper.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 12 deliverable values checked, 0 misaligned (the lone discrepancy is the verification-prose line F1, which is a coverage statement, not a result value).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `#𝔔₄ = binom(5,4) = 5` | py:59 / wl:111 (out py:21, wl:10) | notes:66-72 (`#𝔔₄=binom54=5`) | MATCH |
| triple incidence = 2 | py:64 / wl:112 (out py:23, wl:11) | notes:74 ("exactly two primitive quadruples") | MATCH |
| axis incidence = 4 | py:67 / wl:113 (out py:33, wl:12) | notes (`#𝔔₄=5`, each axis in 4) | MATCH (notes-implied) |
| gradient ray `a_grad=(k_i,k_j,k_k,k_l)/||k||` | py:133 / wl M3 (out py:135, wl:43) | notes:267-271 (boxed) | MATCH |
| max slope `||k||=sqrt(k_i²+k_j²+k_k²+k_l²)` | py:132 / wl:165 (out py:159, wl:47) | notes:275-279 (boxed) | MATCH |
| gradient ratios `k_j/k_i, k_k/k_i, k_l/k_i` | py:142-144 / wl:178-180 | notes:283-289 (boxed) | MATCH |
| synergy gaps `||k||²-||k_face||²=k_missing²` | py:150-153 / wl:199-202 | notes:294 (`(k_grad)²-(k_ijk)²=k_l²>0`) | MATCH |
| `w_Σ ≤ 3`, barycenter `(1,1,1,1)/2` | py:184,190 / wl M5 (out wl:101 `{3,{1/2,..}}`) | notes:363,371-379 (boxed) | MATCH |
| `w_Σ` at triple/pair = 2 / 1 | py:191-192 / wl:221-222 | notes:383-384 | MATCH |
| ten discriminant coeffs `A_⋆..J_⋆` | py:230-239 / wl:256-267 (out wl:124) | notes:461-490 (boxed A_⋆..J_⋆) | MATCH |
| ratio tau bracket `2H0·sqrt(1+r²+s²+t²)/(k_i+rk_j+sk_k+tk_l+sqrt(Δ#))` | py:255-258 / wl:290 | notes:503-508 (boxed) | MATCH |
| equal-mix slope `k_eq=(k_i+k_j+k_k+k_l)/2` | py:309 / (wl via M5 maximizer) | notes:371-373 (barycenter; `k_eq` is its direct dot with k) | MATCH (barycenter-implied) |

Internal scaffolding (no finding): brute-force sample counts (`1500625`, `924`×4), `Min`-flattening booleans, `beta_lo`/`beta_hi` symbolic Min containers, residual-near-zero `= 0` PASS flags, `count_*` iteration counters, `expectZero`/`pass`/`fail` harness output, the `S_quad` packet tuple (bookkeeping container).

Every emitted deliverable value reconciles against the notes (the terse `.tex` card legitimately defers the formulas to the Part VI appendix narrative, per the card's `\stagefield{Verification note}`). No `value_mismatch` and no `script_missing_paper_claim` arises from the reconciliation.

## Self-test notes

Checked: (1) Variable independence — M3's `D[lagObjective - mu(lagNorm-1), x_i]` differentiates an objective that genuinely depends on each `x_i` (linear `k·x` minus quadratic constraint), so the Lagrange gradient is non-degenerate and Solve returns a real branch (confirmed by output line 43). (2) Symmetry/parity — n/a; no unbounded integrals, all checks are algebraic identities or finite interval-arithmetic samples. (3) Trivial-case — the synergy-gap identities reduce to `k_missing²` (a positive literal), the `w_Σ` evaluations give exact 3/2/1, and the discriminant clearing `(1+r²+s²+t²)(k²-2H0κ)-Δ#` is identically 0 by construction of the posited/extracted coefficients (both engines confirm). (4) No missing-script finding, so no path-spec self-test needed. (5) Paper round-trip — no script-side fix is prescribed (F1 routes to user, F2 is verifier-refreshed), so no new paper_misalignment is introduced.
