---
unit_id: 109
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage109_linearized_branch_selection.md]
  paper_appendix: present
---

# Audit unit 109 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_109.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage109_linearized_branch_selection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (present; only references this unit via `\input{stages/stage_109}` at line 1252 — no separate prose row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.txt`

## What the paper claims

Stage 109 linearizes the exact Stage-90 outgoing-quadrupole normalization
`chi_Q = 3(S*beta^5 + 9*Sigma_5)/(3S - Sigma_0)` around the canonical outgoing
branch by writing `S=1+eps*s`, `beta=1+eps*b`, `Sigma_0=eps*a_0`,
`Sigma_5=eps*a_5`. The card's derivation-ledger box (stage_109.tex:16) states the
bottom-line deliverable verbatim: *"First-order outgoing defect is
\(\Delta_Q=5b+a_0/3+9a_5\)."* The notes (lines 27–34) give the same as a boxed
series `chi_Q = 1 + eps*(5b + a_0/3 + 9a_5) + O(eps^2)`. Three distinct
deliverables follow: (i) the overall scale `s` cancels to first order (Class-A
invariance, notes §1); (ii) the first-order defect coefficient is
`5b + a_0/3 + 9a_5` with minimal branch-selection triple `(b, a_0, a_5)`; and
(iii) the linearized preservation condition `5b + a_0/3 + 9a_5 = 0`, equivalently
the boxed closed form `a_5 = -5b/9 - a_0/27` (notes lines 60–73). The card's
secondary `Checks` (Robin/mixed-pole limits, compensated even-coefficient
preservation) are explicitly delegated to downstream stages 110/111/112.

## What the script claims to verify

Both scripts build `chi_Q` from the four linearized inputs and verify four
things: (A1) the eps-series of `chi_Q` equals the literal card claim
`1 + eps*(5b + a0/3 + 9a5)`; (A2) the first-order defect coefficient is
`s`-independent (`d/ds = 0`), i.e. overall scale cancels; (A3) solving the
preservation equation `coeff == 0` for `a5` reproduces the independent closed
form `-5b/9 - a0/27`; and (A4) substituting that independent closed form back
into the defect coefficient annihilates it. The SymPy docstring and the `.wl`
header comment both state that the secondary Robin/mixed-pole/compensation checks
are exercised downstream (110/111/112), matching the card's delegation.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| `chi_Q = 1 + eps*(5b + a0/3 + 9a5)` defect (card:16, notes:27–34) | A1 py:44 / wl:51 — series vs literal `expected` | match |
| overall scale `s` cancels (notes §1) | A2 py:48 / wl:57 — `diff(coeff, s) == 0` | match |
| defect coefficient `= 5b + a0/3 + 9a5`, triple `(b,a0,a5)` (notes §2/§3) | printed py:47 / wl:56 + A1 anchors it | match |
| preservation `a5 = -5b/9 - a0/27` (notes:70–73) | A3 py:54 / wl:62 + A4 py:56 / wl:64 | match |
| Robin/mixed-pole/even-coeff Checks | delegated to 110/111/112 (docstring py:4–11, wl:28–30) | match (out of scope by design) |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful,
non-tautological script-side counterpart; the delegated secondary checks are
correctly out of this stage's scope.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `expect_zero(chi_series - (1 + eps*(5b+a0/3+9a5)))` | defect deliverable (card:16) | yes |
| A2 | sympy | 48 | `expect_zero(diff(coeff, s))` | scale cancellation (notes §1) | yes |
| A3 | sympy | 54 | `expect_zero(a5_sol - (-5b/9 - a0/27))` | preservation closed form | yes (literal anchor) |
| A4 | sympy | 56 | `expect_zero(coeff.subs(a5, -5b/9 - a0/27))` | preservation condition | yes |
| B1 | mathematica | 51 | `expectZero[chiSeries - (1 + eps*(5b+a0/3+9a5))]` | defect deliverable | yes |
| B2 | mathematica | 57 | `expectZero[D[defectCoeff, s]]` | scale cancellation | yes |
| B3 | mathematica | 62 | `expectZero[a5Pres + 5b/9 + a0/27]` | preservation closed form | yes (literal anchor) |
| B4 | mathematica | 64 | `expectZero[(chiSeries - 1) /. a5 -> -5b/9 - a0/27]` | preservation condition | yes |

All eight assertions trace to a specific paper deliverable. None are tautological:
A1/B1 compare an independently-derived series against a hand-typed literal claim;
A2/B2 are real `s`-cancellation tests (the derivative is nonzero precisely if the
cancellation failed); A4/B4 substitute the *independent* literal closed form (not
the self-solved value) into the defect coefficient. A3/B3 are mildly
self-referential (`a5_sol` is solved from `coeff==0`) but anchored to an
independent literal and backstopped by A4/B4.

## Findings

None.

### Adversarial attacks attempted (all failed)

1. **Tautology on A1/B1.** `expected` is a hand-typed literal, while `chi_series`
   is derived by an actual eps-expansion of the rational `chi_Q`. The residual is
   nonzero if the algebra is wrong, so the check can fail. Confirmed by hand: the
   expansion gives `1 + eps*(5b + 9a5 + a0/3)` (the `s` terms cancel between
   numerator and denominator at first order). Not tautological.

2. **Self-solve trap on A3/B3.** `a5_sol`/`a5Pres` is solved from the same
   `coeff==0`/`chiSeries-1==0`, so comparing it to the closed form partly tests
   SymPy/MMA's solver rather than the physics. But this is explicitly
   de-tautologized by A4/B4, which substitute the INDEPENDENT literal
   `-5b/9 - a0/27` into the defect coefficient — a check that fails unless the
   closed form genuinely annihilates `5b + a0/3 + 9a5`. Hand-check:
   `5b + a0/3 + 9*(-5b/9 - a0/27) = 5b + a0/3 - 5b - a0/3 = 0`. Real check.

3. **Scale-cancellation degeneracy on A2/B2.** `diff(coeff, s)`: `coeff` is the
   genuine first-order coefficient of the full ratio, and `s` is one of its
   declared symbols, so the derivative is not identically zero by construction —
   it is zero only because the `s` terms cancelled in the ratio. (Self-test trap
   #1: `coeff` does depend on `s` BEFORE simplification of the ratio; the
   vanishing of `d/ds` after the ratio is the substantive content.) Real check.

4. **Symbol-assumption trap.** All symbols declared `real` (py:30–31, wl:33) with
   no positivity assumptions. No `assume(positive=True)` that could mask a branch.
   The denominator `3S - Sigma_0 = 3 + eps*(3s - a0)` is nonzero at `eps=0`
   (equals 3), so the series/`1/denLin` inversion is valid at the expansion point.
   No domain error.

5. **Transliteration.** See Independent-derivation check — the `.wl` uses a
   distinct denominator-inversion (geometric-series) route, not a port of the
   SymPy `series(whole-ratio)` route. Not transliteration.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent** route, not a transliteration of the `.py`. The
load-bearing chi-series step differs materially:

- **SymPy (py:39–40):** forms the *whole rational* `chi = simplify(3*(S*beta^5 +
  9*Sigma5)/(3*S - Sigma0))` and Taylor-expands the entire quotient in one shot:
  `chi_series = sp.series(chi, eps, 0, 2).removeO()`.
- **Mathematica (wl:42–47):** never series-expands the whole ratio. It expands
  `num` and `den` separately, series-truncates each to order 1 (`numLin`,
  `denLin`), then forms the denominator inverse via its own series
  `invDenLin = Normal[Series[1/denLin, {eps,0,1}]]`, and multiplies
  `numLin*invDenLin`. This is the geometric-series-of-the-denominator route —
  a genuinely different choreography (the same independence pattern accepted for
  stage 100).

The downstream three steps (coefficient read-off, `D[..,s]` scale check, `Solve`
for `a5`, substitution) are structurally parallel between engines, but those are
canonical single-line operations any independent author would write identically;
they are not evidence of transliteration. The distinguishing derivation — how
`chi_Q`'s first-order term is obtained — is independent. **Verdict: independent.**

## Engine cross-check

Both engines emit the same deliverables (modulo print ordering):

| quantity | SymPy output | Mathematica output |
|---|---|---|
| `chi_Q series` | `a0*eps/3 + 9*a5*eps + 5*b*eps + 1` (txt:5) | `1 + (a0*eps)/3 + 9*a5*eps + 5*b*eps` (txt:5) |
| `linearized chi law` residual | `0` (txt:6) | `0` + PASS (txt:6–7) |
| `defect coefficient` | `a0/3 + 9*a5 + 5*b` (txt:7) | `a0/3 + 9*a5 + 5*b` (txt:8) |
| `overall scale cancels` | `0` (txt:8) | `0` + PASS (txt:9–10) |
| `a5 preservation condition` | `-a0/27 - 5*b/9` (txt:9) | `-1/27*a0 - (5*b)/9` (txt:11) |
| `preservation substitution` | `0` (txt:11) | `0` + PASS (txt:14) |

All residuals are identically zero and all symbolic forms agree. Engines agree.

## Verdict justification

`clean`. The card and notes claim a single load-bearing deliverable — the
first-order defect `Delta_Q = 5b + a0/3 + 9a5` from linearizing the Stage-90
`chi_Q` formula — plus the `s`-cancellation and the `a5` preservation closed form.
Both engines build `chi_Q` from the four linearized inputs and verify all three
non-tautologically; the previously-suspect self-solve check is de-tautologized by
the independent-literal substitution A4/B4. I attacked the four assertion classes
(tautology, self-solve, scale-cancellation degeneracy, symbol-assumption) and
each held up, and verified the algebra by hand. The Mathematica route is
genuinely independent (separate num/den series + denominator inversion). Outputs
are fresh (txt mtimes 11:23 > script mtimes 10:57). Paper alignment is exact.

## Self-test notes

Checked: (1) variable independence — `diff(coeff, s)`/`D[defectCoeff, s]` operate
on a coefficient that genuinely contains `s` prior to the ratio's cancellation, so
the derivative is not vacuously zero; the vanishing IS the content. (3)
trivial-case — substituting `a5 = -5b/9 - a0/27` into `5b + a0/3 + 9a5` reduces by
hand to `0`, and a wrong closed form (e.g. `-5b/9`) would leave a nonzero `-a0/3`,
so A4/B4 can fail. (Symmetry/parity #2 n/a — no integrals.) No directive written;
zero findings.

## Value Reconciliation (pass-2 augmentation)

Every result value the two scripts emit as a labeled deliverable maps to the card
and/or notes. The scripts pin no numeric constants; all deliverables are symbolic.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `chi_Q = 1 + eps*(5b + a0/3 + 9a5)` | py:41 / wl:50, sympy txt:5, mma txt:5 | notes:27–34 (boxed); card:16 (`\Delta_Q=5b+a_0/3+9a_5`) | MATCH |
| first-order defect coeff `= 5b + a0/3 + 9a5` | py:47 / wl:56, sympy txt:7, mma txt:8 | card:16; notes:31, §2/§3 (45–58) | MATCH |
| `a5` preservation `= -5b/9 - a0/27` | py:52 / wl:61, sympy txt:9, mma txt:11 | notes:70–73 (boxed `a_5=-\frac{5b}{9}-\frac{a_0}{27}`) | MATCH |
| preservation condition `5b + a0/3 + 9a5 = 0` | A4 py:56 / B4 wl:64 (residual 0) | notes:64 (boxed) | MATCH |

INTERNAL (verification scaffolding, no prose expected): `linearized chi law`
residual, `overall scale cancels` residual, `a5 preservation closed-form`
residual, `preservation substitution` residual, PASS flags — all are zero-residual
check values, not deliverables.

reconciliation: complete; 4 values checked, 0 misaligned.
