---
unit_id: 209
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem.md]
  paper_appendix: present
---

# Audit unit 209 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_209.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (subsection `app-part06-pairwise-ratio-optimizer`, lines 878–918; row line 49)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: *"Finite pairwise candidate set, optimized pairwise bracket, special diagonal-neutral and pair-symmetric reductions, pairwise promotion theorem, and mixed-pair winner theorem."* The card + notes + appendix establish a chain of exact algebraic identities for the one-parameter pairwise ratio search on a compact window `[0,R_ij]`: (1) the certified root function `tau_{ij,*}(r) = 2H0√(1+r²)/(k_i + r k_j + √(A+Br+Cr²))` with `A=k_i²−2H0 u`, `B=2k_i k_j−4H0 v`, `C=k_j²−2H0 w`; (2) the discriminant-numerator reduction `(1+r²)(k_ij²−2H0 κ_ij) = A+Br+Cr²`; (3) the denominator functional `Phi` with `tau = 2H0/Phi`; (4) the exact stationary numerator `N = 2(k_j−k_i r)√(A+Br+Cr²) + B + 2(C−A)r − Br²` and the derivative law `dPhi/dr = N/(2(1+r²)^{3/2}√disc)`; (5) the quartic elimination `Q = [B+2(C−A)r−Br²]² − 4(k_j−k_i r)²(A+Br+Cr²)`, degree 4, with `N` as one factor (extraneous roots = the other factor); (6) two special reductions — diagonal-neutral (`u=w=κ_*, v=0` → gradient-optimal ray `r=k_j/k_i`) and pair-symmetry (`k_i=k_j, u=w` → `r→1/r` invariance, equal-mix ray `r=1` critical). The finiteness/bracket/promotion/winner theorems are combinatorial wrappers on top of these identities; the script verifies the underlying exact algebra (items 1–6), which is the StatusExactClosure content the card marks. Items 5/7 (finite-set counting, promotion/winner inequalities) are not symbolic identities and are not separately scriptable; they follow from "Q is degree 4" (verified) plus the bracket monotonicity (a stated lemma, not an algebra check).

## What the script claims to verify

Both scripts verify the six exact-closure identities above as symbolic identities in `(k_i,k_j,H0,u,v,w,r)`: the explicit `tau` form (I/M1), the discriminant-numerator reduction and A,B,C extraction (I,VI/M2), the `tau=2H0/Phi` functional and `Phi` derivative law (II/M3,M4), the quartic degree-4 + factorization + `N`-is-a-factor identity (III/M5), and the two special reductions' stationarity (IV,V/M6,M7). No numeric benchmark is emitted; every emitted result is a symbolic closed form. There are no tautological "x defined as expr then assert x==expr" blocks: each posited closed form is checked against an independently assembled object (e.g. `tau` built from `kij`/`kappa` is checked against the A,B,C-coefficient form; `Phi'` from `sp.diff` checked against the manifest `N`).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Root function `tau_{ij,*}(r)` (eq app-part06-pairwise-root-function) | sympy I `explicit algebraic tau form`; wl M1 | match |
| A,B,C + discriminant numerator (eq …-ABC, …-discriminant-numerator) | sympy I `discriminant numerator reduction`, VI; wl M2 | match |
| Denominator functional `Phi`, `tau=2H0/Phi` | sympy II `tau=2H0/Phi`; wl M3 | match |
| Stationary numerator `N` + derivative law (eq …-stationary-numerator) | sympy II `Phi derivative law`; wl M4 | match |
| Quartic `Q`, degree 4, `N` a factor (eq …-pairwise-quartic) | sympy III; wl M5 (resultant) | match |
| Diagonal-neutral → gradient ray `k_j/k_i` (notes 8.1) | sympy IV; wl M6 | match |
| Pair-symmetry → equal-mix ray `r=1` critical (notes 8.2) | sympy V; wl M7 | match |
| Finite candidate set / 12-eval count (notes 6) | implied by `degree(Q)=4`; no direct count check | partial (acceptable — counting, not an identity) |
| Optimized bracket / promotion / winner theorems (notes 7,9) | none (combinatorial wrappers) | missing-but-acceptable (not symbolic identities) |

The card's StatusExactClosure content is the six exact identities; all six are faithfully exercised by both engines. The bracket/promotion/winner deliverables are order/comparison logic on top of those identities (StatusNumerical insertions per the card), not algebra a CAS can certify — their absence is appropriate, not a finding.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 68 | `tau - tau_expected == 0` | root function (1) | yes |
| A2 | sympy | 69–72 | `(1+r²)(kij²−2H0κ) − (A+Br+Cr²) == 0` | disc numerator (2) | yes |
| A3 | sympy | 88 | `tau − 2H0/Phi == 0` | functional (3) | yes |
| A4 | sympy | 89 | `diff(Phi,r) − N/(2(1+r²)^{3/2}S) == 0` | derivative law / N (4) | yes |
| A5 | sympy | 105 | `degree(Q) − 4 == 0` | quartic degree (5) | yes |
| A6 | sympy | 106 | `Q − (J−2(kj−kir)S)(J+2(kj−kir)S) == 0` | quartic factorization (5) | yes |
| A7 | sympy | 107 | `N − (J + 2(kj−kir)S) == 0` | N is the +radical factor (5) | yes |
| A8 | sympy | 121–122 | `κ_diag − κ_*==0`, `diff(tau_diag,r)|_{kj/ki}==0` | diagonal-neutral (6) | yes |
| A9 | sympy | 132–133 | `tau_sym − tau_sym(1/r)==0`, `diff(tau_sym,r)|_1==0` | pair-symmetry (6) | yes |
| M1 | wl | 102–105 | `closureRoot − coefficientRoot == 0` | root function (1) | yes |
| M2 | wl | 108–115 | `abcByExtraction − abcClosedForm == 0`, disc residual | disc numerator (2) | yes |
| M3 | wl | 118–121 | `closureRoot − 2H0/Phi == 0` | functional (3) | yes |
| M4 | wl | 140–154 | `D[Phi]≠0`; `scaledDerivNum − manifestN == 0`; deriv-law residual | derivative law / N (4) | yes |
| M5 | wl | 166–183 | resultant degree 4; factorization; `+radical factor == derivNum` | quartic (5) | yes |
| M6 | wl | 188–206 | gradient ratio from `Solve[D[raySlope]==0]`; stationarity at `kj/ki` | diagonal-neutral (6) | yes |
| M7 | wl | 211–222 | reciprocal invariance; stationarity at `r=1` | pair-symmetry (6) | yes |

Every script-side check traces to a specific paper deliverable. No orphaned assertions, no `paper_missing_script_claim`.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.txt` (mtime 2026-06-02 10:54:33)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage209_pairwise_ratio_optimizer_and_mixed_ray_winner_theorem_sympy_audit.py` (mtime 2026-06-03 15:59:11)

**What's wrong:**
The committed SymPy output `.txt` is older than the `.py` (output 06-02 10:54 < script 06-03 15:59), so it does not reflect the current script. Concretely, the captured output banners read `STAGE 192 — EXACT PAIRWISE RATIO OPTIMIZER…` / `STAGE 192 SYMPY AUDIT COMPLETED SUCCESSFULLY` (txt lines 3, 149), whereas the current `.py` prints `STAGE 209 …` (py line 35) / `STAGE 209 SYMPY AUDIT COMPLETED SUCCESSFULLY` (py line 148). This is the known P4-52 stale-banner artifact in the script/output numbering band. The Mathematica output is fresh (wl 06-02 10:47 < txt 06-02 10:54) and already prints `STAGE 209`.

**Why this matters:**
Only the banner text differs; all mathematical content in the stale SymPy `.txt` (the symbolic forms, the `= 0` residuals) matches what the current `.py` would emit — the algebra was unchanged between captures. So this is informational, not blocking: the captured residuals are still all zero and agree with Mathematica. The stale banner label belongs to the deferred dedicated script/output-band numbering pass, not to this audit's math scope.

**Required change:**
No Codex script edit. The orchestrator's independent re-run of `python3 …stage209…_sympy_audit.py` refreshes the `.txt` (the current `.py` already emits the correct STAGE 209 banner). No directive is written; the numbering-band banner concern is owned by the separate deferred plan, not this stage's fix loop.

**Verification:**
After the orchestrator re-run, the SymPy `.txt` mtime > `.py` mtime and its first banner reads `STAGE 192` → `STAGE 209`, with all `… = 0` residual lines unchanged.

## Independent-derivation check (Mathematica)

**Verdict: INDEPENDENT.** The `.wl` is not a transliteration; it derives the load-bearing objects by methods that differ from the `.py`'s posit-and-verify route at every key step.

1. **A,B,C / discriminant numerator (load-bearing object 2).** SymPy *posits* the literals `A = ki**2 - 2*H0*u`, `B = 2*ki*kj - 4*H0*v`, `C = kj**2 - 2*H0*w` (py 50–52) and verifies them by an identity (py 71). Mathematica *extracts* them from the assembled polynomial: `discPoly = Expand[Numerator[Together[den[r](raySlope[r]^2 - 2 H0 rayCurvature[r])]]]; discCoefficients = CoefficientList[discPoly, r]` (wl 68–78), then checks the extracted coefficients equal the closed form (M2, wl 108–111). Derive (coefficient extraction) vs posit (literal) — independent.

2. **Quartic Q (load-bearing object 5).** SymPy *posits* the difference-of-squares form directly: `Q = sp.expand(J**2 - 4*(kj-ki*r)**2*(A+B*r+C*r**2))` (py 97). Mathematica *eliminates the radical via a resultant*: `quarticByResultant = Factor[Resultant[polyPartFromDifferentiatedDiscriminant + 2 (kj - ki r) z, z^2 - discFromCoefficients[r], z]]` (wl 157–163), introducing `z=√disc` and computing `Resultant[…, z]` to clear it. The resultant elimination is a genuinely different algebraic route to the quartic than writing the difference of squares by hand. This is the single strongest discriminator.

3. **Stationary numerator N / gradient ray (objects 4, 6).** SymPy posits `N` (py 80) and posits the gradient ratio `r_grad = kj/ki` (py 115). Mathematica *solves out* N by scaling the differentiated functional — `stationaryNumeratorFromDerivative = FullSimplify[2 den[r]^(3/2) Sqrt[disc] * D[denominatorFunctional[r], r]]` (wl 124–132) — and *solves* the slope-stationarity equation for the gradient ratio — `gradientRatioFromSlopeSolve = r /. Solve[D[raySlope[r], r] == 0, r, Reals]` (wl 187), checking `kj/ki` is a member (wl 192–198). Derive vs posit again.

The shared physical premises (the ray slope `(ki + r kj)/√(1+r²)`, the ray curvature `(u+2vr+wr²)/(1+r²)`, the closure-root form) are the legitimately-shared monomial definitions; the methods that EXTRACT the load-bearing objects (A,B,C by `CoefficientList`; Q by `Resultant`; N by scaled differentiation; gradient ray by `Solve`) differ from the SymPy posit-and-check route. The "each CAS runs its own simplifier" defense is not invoked — the independence is structural, not simplifier-based.

## Engine cross-check

Both engines pass all checks with zero residuals. Cross-checking the load-bearing quartic: the Mathematica M5 resultant (wl/`.txt` line 39) yields, at order r⁴, `8 H0 (-2 ki kj v + 2 H0 v^2 + ki^2 w)`, and at order r⁰, `8 H0 (kj^2 u - 2 ki kj v + 2 H0 v^2)`. The SymPy section VI emits exactly the same leading coefficient `8⋅H₀⋅(2⋅H₀⋅v² + kᵢ²⋅w − 2⋅kᵢ⋅k_j⋅v)` and constant coefficient `8⋅H₀⋅(2⋅H₀⋅v² − 2⋅kᵢ⋅k_j⋅v + k_j²⋅u)` (sympy `.txt` lines 134–146). The two independent routes (difference-of-squares vs resultant) produce the identical quartic. Engines agree.

## Verdict justification

`findings`, count 1 — a single low-severity informational `stale_output` (SymPy `.txt` predates `.py`; only the STAGE-192-vs-209 banner differs, all math residuals match, refreshed by the orchestrator re-run). I read the paper card, the notes, and the Part VI appendix subsection first, and the script's verified claim matches the paper exactly: all six exact-closure identities are present and faithfully exercised; the absent bracket/promotion/winner theorems are combinatorial/comparison wrappers (StatusNumerical), not CAS-certifiable identities, so their absence is correct. Attacks tried and failed: (a) tautology hunt — every posited closed form is checked against an independently assembled object, none is asserted against itself; (b) `assume(positive=True)` abuse — `ki,kj,H0,r>0` is justified by the physical setup (positive slopes, positive throat parameter, ratio `r≥0`), and `u,v,w` are correctly left unsigned (py 53–55, wl 60–63) so the radicand sign is not smuggled in; (c) zero-derivative trap — both engines explicitly assert `D[Phi]≠0` / non-trivial derivative before the reductions (wl 140–143, 188–191), so the stationarity checks are not vacuously satisfied; (d) hand-checked the N derivative law and the Q=N·(other factor) factorization algebraically — both are correct; (e) transliteration attack — the `.wl` derives A,B,C by `CoefficientList`, Q by `Resultant`, N by scaled differentiation, and the gradient ray by `Solve`, none of which mirrors the `.py`'s posit route → INDEPENDENT.

## Value Reconciliation (pass-2 augmentation)

All emitted deliverables are symbolic closed forms (no numeric benchmarks at this stage). Each reconciles against the appendix equations and/or notes.

| value | source (py/wl + output) | .tex/.md location | status |
|---|---|---|---|
| `tau_{ij,*}(r) = 2H0√(1+r²)/(ki + r kj + √(A+Br+Cr²))` | py 58 / `.txt` 24–31; wl 94–96 | appendix eq `…-pairwise-root-function` (part06 l.896–901); notes sec.2 (l.122–129) | MATCH |
| `A = ki²−2H0 u` | py 50; wl 80 | appendix eq `…-pairwise-ABC` (l.884); notes sec.2 (l.99–106) | MATCH |
| `B = 2 ki kj − 4H0 v` | py 51; wl 81 | appendix l.886; notes l.99–106 | MATCH |
| `C = kj²−2H0 w` | py 52; wl 82 | appendix l.888; notes l.99–106 | MATCH |
| disc numerator `Δ♯ = A+Br+Cr²`, with `(1+r²)(kij²−2H0κ)=Δ♯` | py 56,71 / `.txt` 33; wl 113–114 | appendix eq `…-discriminant-numerator` (l.893); notes sec.2 (l.116–120) | MATCH |
| `Phi = (ki + r kj + √disc)/√(1+r²)`, `tau = 2H0/Phi` | py 79,83 / `.txt` 38–45; wl 97–99,118 | notes sec.4 (l.161–172) | MATCH |
| `N = 2(kj−ki r)√disc + B + 2(C−A)r − B r²` | py 80 / `.txt` 46–57; wl 133–138 | appendix eq `…-stationary-numerator` (l.905–907); notes sec.4 (l.186–192) | MATCH |
| derivative law `dPhi/dr = N/(2(1+r²)^{3/2}√disc)` | py 81,89; wl 149–154 | notes sec.4 (l.177–182) | MATCH |
| `Q = [B+2(C−A)r−Br²]² − 4(kj−ki r)²(A+Br+Cr²)`, degree 4 | py 97,104 / `.txt` 67–82; wl 157,166–168 / `.txt` 39 | appendix eq `…-pairwise-quartic` (l.912–914); notes sec.5 (l.218–225) | MATCH |
| `N = J + 2(kj−ki r)√disc` (N is the +radical factor of Q) | py 107; wl 178–183 | notes sec.5 (l.240, extraneous-root caveat) | MATCH |
| diagonal-neutral optimizer `r = kj/ki` (gradient ray) | py 115,122; wl 187,199–206 | notes sec.8.1 (l.343–346) | MATCH |
| pair-symmetry: `r→1/r` invariance, equal-mix `r=1` critical | py 132,133; wl 211–222 | notes sec.8.2 (l.360–376) | MATCH |
| explicit quartic coefficients (highest→constant, factored) | py 145–146 / `.txt` 133–146 | not boxed in card; matches notes sec.5 quartic implicitly | MATCH (intermediate elaboration of the boxed Q) |

INTERNAL scaffolding (no finding): `kij`, `kappa`, `S`/`Sqrt[disc]`, `J`/`polyPartFromDifferentiatedDiscriminant`, `tau_expected`/`coefficientRoot`, `Phi_prime_expected`, `Qpoly`/`quarticByResultant`, `r_grad`/`gradientRatioFromSlopeSolve`, `tau_diag`/`tau_sym`, all `expect_zero`/`expectZero`/`expectTrue` pass flags.

reconciliation: complete; 13 deliverable values checked, 0 misaligned.
