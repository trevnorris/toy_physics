---
unit_id: 035
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: misaligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md
  paper_appendix: present
---

# Audit unit 035 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_035.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 60 summary; `\input{stages/stage_035}` at line 108)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.txt`

## What the paper claims

Stage 035 substitutes the exact D/N constants `kappa_0^2 = 8/pi^2` and `kappa_1^2 = 16/(9 pi^2)` into the Stage-034 softening-depth normal form, defines the dimensionless variables `xi = x/A`, `delta = DeltaK_ax/A` (with `0 <= xi < 1`, `delta > 0`), and factorises `N_-(x) = N_-(0) F(xi, delta)` with `N_-(0) = beta_0 kappa_0^2 / A`. `\stagefield{Output}` states verbatim: "Stage~035 outputs the shape function \eqref{eq:app-stage035-F}, the target ratio \eqref{eq:app-stage035-Rtarget}, the monotonicity derivative \eqref{eq:app-stage035-F-derivative}, the existence/uniqueness theorem \eqref{eq:app-stage035-existence}, and the required total loading \eqref{eq:app-stage035-alpha-req}." The card boxes (i) `F = (9 delta + 11 xi)^4 / [81 (1-xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2]`, (ii) `R_target = N_Q^target A / (beta_0 kappa_0^2)`, (iii) the explicit `partial_xi F` formula with bracket polynomial `81 delta^3 + 72 delta^2 + 206 delta^2 xi + 297 delta xi^2 + 138 xi^3`, (iv) the trichotomy `R_target <,=,> 1`, and (v) `alpha_req = 9 pi^2 A xi (xi+delta) / [8 (9 delta + 11 xi)]`. The notes additionally specify `F(0,delta) = 1`, `lim_{xi -> 1^-} F = +infinity`, the asymptotic constant `C(delta)`, the near-onset Taylor expansion through `O(xi^2)`, and the near-softening tail `1 - xi_req ~ C(delta)/R_target`. The notes' Section 3 polynomial matches the paper card's (206, 138).

## What the script claims to verify

Both engines independently (a) construct `N_-(x)` from `kappa0^2, kappa1^2, A, x, DeltaK_ax` and form `F = N_-(x) / [beta_0 kappa_0^2 / A]`, asserting it equals the closed `F_target`; (b) differentiate that `F_target` and assert the result equals a hand-encoded "manifestly positive" form with bracket polynomial `81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3` — note the literal coefficients are `189` and `121`, not the paper's `206` and `138`; (c) verify `F(0,delta) = 1` and the softening leading coefficient `C(delta) = (9 delta + 11)^4 / [81 (9 delta^2 + 18 delta + 11)^2]` via a direct limit and via `eps_soft * F(1 - eps_soft) -> C(delta)`; (d) check `alpha_req` and its `xi -> 1^-` limit `alpha_crit`; (e) print (without assertion) the support-coupling combination `g_B^2/varpi^2 = alpha_req - alpha_mix`; (f) verify the near-onset Taylor expansions of `F` and `alpha_req` through `O(xi^2)`. The Mathematica script additionally asserts `R_target == Pi^2 A NQ / (8 beta0)`, which is arithmetic on the substituted constant `kappa0^2 = 8/Pi^2`. All assertions pass in both saved transcripts.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Closed-form `F(xi,delta)` (eq. app-stage035-F) | A1 (sympy) / B1 (mma) | match |
| `R_target = N_Q^target A / (beta_0 kappa_0^2)` (eq. Rtarget) | sympy line 60 print; mma B2 asserts equality to the rationalised form | match |
| Monotonicity derivative `partial_xi F` (eq. F-derivative) | A2 / B3 | **mismatch** — script asserts polynomial `... 189 delta^2 xi + 297 delta xi^2 + 121 xi^3`; paper card boxes `... 206 delta^2 xi + 297 delta xi^2 + 138 xi^3`; notes match the paper's wrong values |
| Existence/uniqueness theorem (eq. existence) | Strict positivity of `dF/dxi` (manifest from coefficients) + endpoints | partial — no explicit trichotomy assertion, but every input the theorem needs is verified |
| `alpha_req` (eq. alpha-req) | A6 / B7 | match |
| `alpha_crit` (notes Sec. 4; paper eq. alpha-crit) | A7 / B8 | match |
| Endpoint `F(0,delta) = 1` (notes Sec. 3; paper eq. F-endpoints) | A3 / B4 | match |
| `C(delta)` (notes Sec. 3; paper eq. Cdelta) | A4-A5 / B5-B6 | match |
| Near-onset Taylor `F = 1 + (1 + 8/(9 delta)) xi + (1 + 8/(9 delta) - 28/(27 delta^2)) xi^2 + O(xi^3)` (notes 6.1; paper eq. onset-expansion) | A8 / B9 | match |
| `alpha_req` near-onset leading term (notes 6.1) | A9 / B10 | extra — verifies through `O(xi^2)`; notes only state leading `~ pi^2 A xi / 8`. Not a defect; just over-coverage. |
| Support coupling `g_B,req^2/varpi^2 = alpha_req - alpha_mix` (notes Sec. 5) | Computed and printed in both engines; no assertion | partial — carry-forward only; no algebraic identity to assert |

Dominant pattern: aligned except for a single concrete polynomial-coefficient discrepancy in `partial_xi F` between the paper card / notes (which agree with each other) and both engines' scripts. Setting `paper_alignment: misaligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `expect_zero("F - closed D/N form", F - F_target)` | F closed form | yes |
| A2 | sympy | 73 | `expect_zero("dF/dxi - manifestly positive form", dF - dF_target)` | monotonicity derivative | yes, but `dF_target` literal polynomial differs from paper/notes coefficients (189 vs 206, 121 vs 138) |
| A3 | sympy | 74 | `expect_zero("F(0,delta) - 1", ...)` | endpoint `F(0,delta) = 1` | yes |
| A4 | sympy | 78 | `expect_zero("softening constant - closed form", soft_const - soft_const_target)` | `C(delta)` | yes |
| A5 | sympy | 82 | `expect_zero("near-softening asymptotic coefficient", ...)` | softening leading coefficient | yes |
| A6 | sympy | 89 | `expect_zero("alpha_req - closed D/N form", ...)` | `alpha_req` closed form | yes |
| A7 | sympy | 94 | `expect_zero("alpha_crit - closed form", ...)` | `alpha_crit` | yes |
| A8 | sympy | 105 | `expect_zero("near-onset series through O(xi^2)", ...)` | near-onset Taylor for F | yes |
| A9 | sympy | 110 | `expect_zero("alpha_req near-onset series through O(xi^2)", ...)` | extends notes Sec. 6.1 leading term | partial (over-coverage, not a defect) |
| B1 | mma | 61 | `expectZero["F - closed D/N form", f - fTarget]` | F closed form | yes |
| B2 | mma | 62 | `expectZero["R_target - Pi^2 A NQ/(8 beta0)", rTarget - rTargetClosed]` | `R_target` rationalisation | partial — `rTarget = NQ*A/(beta0*(8/Pi^2))` and `rTargetClosed = Pi^2*A*NQ/(8*beta0)` are the same rational expression. The check verifies the constant substitution arithmetically, not a new claim. |
| B3 | mma | 71 | `expectZero["dF/dxi - manifestly positive form", dF - dFTarget]` | monotonicity derivative | same status as A2 |
| B4 | mma | 72 | `expectZero["F(0,delta) - 1", ...]` | endpoint | yes |
| B5 | mma | 83 | `expectZero["softening constant - closed form", ...]` | `C(delta)` | yes |
| B6 | mma | 90 | `expectZero["near-softening asymptotic coefficient", ...]` | softening leading coefficient | yes |
| B7 | mma | 105 | `expectZero["alpha_req - closed D/N form", ...]` | `alpha_req` | yes |
| B8 | mma | 106 | `expectZero["alpha_crit - closed form", ...]` | `alpha_crit` | yes |
| B9 | mma | 122 | `expectZero["near-onset F series through O(xi^2)", ...]` | near-onset Taylor for F | yes |
| B10 | mma | 123 | `expectZero["alpha_req near-onset series through O(xi^2)", ...]` | over-coverage | partial |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** `target_mismatch` (and equivalently `notes_contradicts_script`, since the notes carry the same wrong polynomial as the paper card)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_035.tex:65-73` (boxed equation `eq:app-stage035-F-derivative`)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md:85-87`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py:65-70`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:64-69`

**What's wrong:**

The paper card boxes (lines 65-73):

```
\frac{\partial F}{\partial\xi}
=
\frac{(9\delta+11\xi)^3
\left(81\delta^3+72\delta^2+206\delta^2\xi+297\delta\xi^2+138\xi^3\right)}
{81(1-\xi)^2(9\delta^2+18\delta\xi+11\xi^2)^3}>0.
```

The notes (lines 85-87) say:

```
dF/dxi
= (9 delta + 11 xi)^3
  [ 81 delta^3 + 72 delta^2 + 206 delta^2 xi + 297 delta xi^2 + 138 xi^3 ]
  / [ 81 (1 - xi)^2 (9 delta^2 + 18 delta xi + 11 xi^2)^3 ]
> 0.
```

Both engines, however, assert the bracket polynomial is `81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3`. SymPy `dF_target` (lines 67-69):

```
(9 * delta + 11 * xi) ** 3
* (81 * delta**3 + 72 * delta**2 + 189 * delta**2 * xi + 297 * delta * xi**2 + 121 * xi**3)
/ (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)
```

Mathematica `dFTarget` (lines 66-68) is identical in coefficients. The Mathematica saved transcript (line 12) prints the expanded numerator as

```
(9*delta + 11*xi)^3*(81*delta^3 + 297*delta*xi^2 + 121*xi^3 + 9*delta^2*(8 + 21*xi))
```

i.e., `9*delta^2*(8 + 21*xi) = 72 delta^2 + 189 delta^2 xi`, confirming the literal `189` and `121`. The script's `expect_zero("dF/dxi - manifestly positive form", dF - dF_target)` passes (sympy output line 34), so the engine-derived `D[F_target, xi]` reduces to exactly the script's hand-encoded form, not the paper's.

Direct hand derivation confirms the script. Let `M = 9 delta^2 + 18 delta xi + 11 xi^2`, `M' = 18 delta + 22 xi = 2 (9 delta + 11 xi)`. From `F = (9 delta + 11 xi)^4 / [81 (1-xi) M^2]`, the quotient rule and pulling out the common factor `(9 delta + 11 xi)^3 / [81 (1-xi)^2 M^3]` leaves the bracket

```
44 (1-xi) M + (9 delta + 11 xi) M - 4 (1-xi) (9 delta + 11 xi)^2.
```

Expanding (full term-by-term arithmetic):

- `M (44 + 9 delta - 33 xi) = 81 delta^3 + 396 delta^2 - 135 delta^2 xi + 792 delta xi - 495 delta xi^2 + 484 xi^2 - 363 xi^3`.
- `- 4 (1-xi) (9 delta + 11 xi)^2 = -324 delta^2 - 792 delta xi - 484 xi^2 + 324 delta^2 xi + 792 delta xi^2 + 484 xi^3`.

Adding: `81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3`. Exactly the script's polynomial.

Quick numerical sanity at `delta = 0`: `F = (11 xi)^4 / [81 (1-xi)(11 xi^2)^2] = 121 / [81 (1-xi)]`, so `dF/dxi = 121 / [81 (1-xi)^2]`. The script's form at `delta = 0` reduces to `(11 xi)^3 * 121 xi^3 / [81 (1-xi)^2 (11 xi^2)^3] = 121 / [81 (1-xi)^2]` ✓. The paper's form at `delta = 0` gives `138 / [81 (1-xi)^2]` ≠ true derivative.

Therefore the scripts verify the correct derivative; the paper card and notes have arithmetic errors in two coefficients (`206 -> 189` and `138 -> 121`). Positivity itself is unaffected — both polynomials have strictly positive coefficients on `delta, xi > 0`, so the monotonicity claim still holds either way — but the boxed identity in the paper is numerically wrong.

**Why this matters:**

`eq:app-stage035-F-derivative` is referenced from `\stagefield{Output}` and is a load-bearing boxed identity. A boxed identity with wrong coefficients undermines the exact-closure claim and would propagate to any downstream consumer that quotes the polynomial. Independent readers reproducing the derivative will see a mismatch.

**Required change:**

This is `paper_misalignment` — Codex does NOT auto-resolve. The orchestrator halts for user resolution. The expected direction (based on independent hand derivation plus agreement of both engines) is paper-side: update both the .tex (line 71) and the notes (line 86) so the bracket polynomial reads `81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3`. The scripts should be left as-is. The directive routes the question to the user; see `## Resolve before fix_loop` in the directive.

**Verification:**

After the user approves the paper-side fix and a follow-up directive applies it, `paper/stages/stage_035.tex` line 71 must read `\left(81\delta^3+72\delta^2+189\delta^2\xi+297\delta\xi^2+121\xi^3\right)`, and the notes' Section 3 polynomial must match. No re-run of either engine is required; the script assertions already pass against the correct polynomial.

## Independent-derivation check (Mathematica)

The Mathematica script mirrors the SymPy choreography closely — same intermediate variables (`kappa0Sq`, `kappa1Sq`, `x = A xi`, `deltaKSub`, `alphaX`, `nX`, `f`), same hand-encoded literal target polynomials, same assertion order. This is consistent with the project's documented Mathematica-mirror policy (per the `MATHEMATICA_MIRROR_POLICY` tracker noted in user memory), and the substantive checks still call Mathematica's own `D[f, xi]`, `Limit[...]`, `Series[...]`, and `FullSimplify` — i.e., the algebraic engines are exercised independently even though the targets are typed in the same form.

Specifically, for the load-bearing `dF/dxi` check: SymPy computes `sp.diff(F_target, xi)` and `sp.simplify` reduces it against the literal `dF_target`; Mathematica computes `D[f, xi]` (where `f` is derived from `nX`, the constructed `N_-(x)`, not from the literal closed form) and `FullSimplify` reduces it against the literal `dFTarget`. The literal target polynomials are typed independently in each script, but they happen to be the same (and correct) polynomial. If a typo had been introduced into the SymPy literal alone, the Mathematica engine would still catch it only via the Mathematica literal — i.e., the two literals are independent typings.

I do not flag this as `mathematica_transliteration`. The two engines independently expand the derivative and series and obtain the same simplified residual against independently typed literals. The literal polynomial is the claim under test, not derived from the SymPy script's output. The cross-engine agreement on `D[f, xi]` simplifying to the same polynomial is a genuine independent check.

## Engine cross-check

Both engines produce equal closed forms and equal zero residuals:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `F(xi,delta)` | `-(9δ + 11ξ)^4 / [81(ξ-1)(2ξ² + 9(δ+ξ)²)²]` | `-1/81*(9*delta + 11*xi)^4/((-1 + xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2)` (algebraically equal) |
| `dF/dxi - hand-encoded target` | 0 | 0 |
| `C(delta)` | `(9δ+11)^4 / [81(9δ²+18δ+11)^2]` | `(11+9δ)^4 / [81(11+9δ(2+δ))^2]` (algebraically equal) |
| `alpha_req` | `9π²A ξ(δ+ξ)/[8(9δ+11ξ)]` | `9 A Pi^2 xi(δ+xi)/(72δ+88xi)` (algebraically equal) |
| `alpha_crit` | `9π²A(δ+1)/[8(9δ+11)]` | `9 A(1+δ)Pi²/(88+72δ)` (algebraically equal) |
| F series O(xi^2) residual | 0 | 0 |
| alpha series O(xi^2) residual | 0 | 0 |

`engines_agree: true`.

## Verdict justification

Every script-side assertion correctly verifies a non-trivial identity, and the residuals are correctly zero. No tautological checks (B2 is borderline — comparing two equivalent rationalisations of `R_target` — but it does exercise the constant substitution and is not circular). No symbol-domain errors (`delta > 0`, `0 <= xi < 1` are explicitly declared and match the paper). No missing branches. No `hardcoded_result` (every literal closed form is independently expanded by each engine from its constructed `f`/`alphaX`). Outputs are fresh (both transcripts younger than their scripts).

The one defect is a paper-side polynomial-coefficient mismatch in the boxed `dF/dxi` formula: paper says `206 delta^2 xi` and `138 xi^3`; both engines (and independent hand derivation, verified numerically at `delta = 0`) give `189` and `121`. The script side is correct; the paper card and notes need updating. Verdict: `findings`, one `paper_misalignment` of subtype `target_mismatch`. Not `stop_cold`: the math is reconcilable and the resolution direction is clear, but user sign-off is required because Codex must not silently edit paper or notes.

## Self-test notes

- Variable independence: every `sp.diff(EXPR, xi)` and `D[expr, xi]` is taken w.r.t. a variable that appears in `EXPR` (verified by inspection of `F_target`, `alphaX`, etc.); no zero-derivative trap.
- Symmetry/parity: not applicable (no unbounded-domain integrals; all checks are algebraic identities, limits, or Taylor series).
- Trivial-case pre-check: substituted `xi = 0` mentally — `F(0, delta) = (9 delta)^4 / [81 * 1 * (9 delta^2)^2] = 6561 delta^4 / (6561 delta^4) = 1` ✓; substituted `delta = 0, xi != 0` into both polynomial forms to confirm the script's `121` reproduces the true `dF/dxi` while the paper's `138` does not (see F1).
- Path specifications: no missing-script findings, so no path-only edits.
- Paper round-trip: the only finding routes to user resolution; no Codex-applied script-side change is contemplated, so there is no risk of introducing a new `paper_misalignment` by mis-citing a constant.
