---
unit_id: 046
batch: III.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md
  paper_appendix: present
---

# Audit unit 046 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_046.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 70)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt`

## What the paper claims

Stage 046 establishes the monotonicity of the tracking-branch functions in R and the residual-comparison theorem. The card's `\stagefield{Output}` reads verbatim: "Monotonicity \eqref{eq:app-stage046-monotonicity} and residual theorem \eqref{eq:app-stage046-residual}." The two boxed equations are (i) `partial G_tr/partial R < 0` and `partial F_tr/partial R > 0` and (ii) `E_tr - E_flat = F_flat - F_tr > 0`. The card also states the implied inequalities `G_tr > G_flat`, `F_tr < F_flat` at fixed `(xi, delta)` (eq:app-stage046-inequalities). The companion notes give explicit closed forms for `G_tr`, `F_tr`, their R-derivatives, the branch differences `G_tr - G_flat` and `F_flat - F_tr`, the strong-split endpoints `G_tr(R=0) = xi`, `F_tr(R=0) = 1/(1-xi)`, and the corresponding positive auxiliary polynomials `P_R`, `P_1`, `P_2`. The part-03 appendix row 70 calls the stage "Tracking branch bounds — Monotonicity in R and proof that first split-U tracking worsens the target." `\stagefield{Inputs}` names `F_tr`, `G_tr`, and the constructive inequality `R_tr < 1`.

## What the script claims to verify

Both scripts define `G_tr = 9 xi (xi + delta) / (9 delta + (9 + 2 R^2) xi)` and `F_tr = (9 delta + (9 + 2 R^2) xi)^2 (9 delta + (9 + 2 R) xi)^2 / [81 (1-xi) (9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2)^2]`. They verify (1) the endpoint identities `G_tr(R=0) = xi`, `F_tr(R=0) = 1/(1-xi)`; (2) the closed-form R-derivatives with explicit factored expressions; (3) the closed-form branch differences `G_tr - G_flat` and `F_flat - F_tr`; (4) positivity of the auxiliary polynomials `P_R`, `P_1`, `P_2` (SymPy via coefficient inspection, Mathematica via `Reduce[ForAll ...]`); (5) three numerical sample points showing both differences are strictly positive; (6) boundary-value checks at R=1. The Mathematica script additionally proves universal strict-sign claims (`dG/dR < 0`, `dF/dR > 0`, `G_tr - G_flat > 0`, `F_flat - F_tr > 0`) over the entire bound domain `0 < R, xi < 1`, `delta > 0` via `Reduce[ForAll, ..., Reals]`.

## Paper -> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `partial G_tr / partial R < 0` (monotonicity, eq:app-stage046-monotonicity) | sympy: `dG_tr/dR formula = 0` (closed form); mathematica: `Reduce[dG/dR < 0 on (0,1)^3] = True` | match |
| `partial F_tr / partial R > 0` (monotonicity, eq:app-stage046-monotonicity) | sympy: `dF_tr/dR formula = 0` + `P_R` positivity; mathematica: `Reduce[dF/dR > 0 on (0,1)^3] = True` | match |
| `G_tr > G_flat`, `F_tr < F_flat` (eq:app-stage046-inequalities) | sympy: `G_tr - G_flat formula = 0` + `P_1`, `P_2` positivity + 3 samples; mathematica: `Reduce[G_tr - G_flat > 0]`, `Reduce[F_flat - F_tr > 0]` | match |
| Residual theorem `E_tr - E_flat = F_flat - F_tr > 0` (eq:app-stage046-residual) | covered as direct algebraic consequence of `F_flat - F_tr > 0`: `E_tr - E_flat = (R_target - F_tr) - (R_target - F_flat) = F_flat - F_tr` is identity by definition; script verifies the substantive `F_flat - F_tr > 0` half | match (the residual identity itself is tautological from the definitions in notes section 5; the load-bearing content is `F_flat - F_tr > 0`) |
| Strong-split endpoints `G_tr(R=0)=xi`, `F_tr(R=0)=1/(1-xi)` (notes section 4) | both scripts: `expect_zero` at R=0 | match (paper card omits, notes assert) |
| Boundary `G_tr(R=1) = G_flat`, `F_tr(R=1) = F_flat` | both scripts: `expect_zero` at R=1 | match |
| Explicit coefficient values of `P_R`, `P_1`, `P_2` (notes section 2.2, 3.2) | script polynomials use different coefficients than the notes (see F1 below) | mismatch (notes <-> script) |

Set `paper_alignment: partial` because every paper-card-level deliverable matches, but the notes contain typo-level coefficient mismatches with what the scripts (and both engines, independently) actually compute.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `expect_zero("strong-split endpoint for G", G_tr(R=0) - xi)` | notes endpoint formula | yes |
| A2 | sympy | 59 | `expect_zero("strong-split endpoint for F", F_tr(R=0) - 1/(1-xi))` | notes endpoint formula | yes |
| A3 | sympy | 87 | `expect_zero("dG_tr/dR formula", dG_dR - dG_expected)` | monotonicity (paper eq:app-stage046-monotonicity) | yes (non-trivial: machine-computed Factor against typed closed form) |
| A4 | sympy | 88 | `expect_zero("dF_tr/dR formula", dF_dR - dF_expected)` | monotonicity | yes |
| A5 | sympy | 89 | `expect_positive_coefficients("P_R", ...)` | sign of dF/dR | yes |
| A6 | sympy | 138 | `expect_zero("G_tr - G_flat formula", delta_G - G_diff_expected)` | branch ordering (paper eq:app-stage046-inequalities) | yes |
| A7 | sympy | 139 | `expect_zero("F_flat - F_tr formula", delta_F - F_diff_expected)` | residual theorem (paper eq:app-stage046-residual) | yes |
| A8 | sympy | 140-141 | `expect_positive_coefficients("P1", ...)`, `("P2", ...)` | branch ordering / residual sign | yes |
| A9 | sympy | 147-162 | boundary `expect_zero` at R=0 and R=1 | endpoint sanity | yes |
| A10 | sympy | 179-186 | sample-point assertions `g_sample > 0`, `f_sample > 0` | branch ordering | yes (three interior rational points; not symbolic but real numerical tests) |
| B1 | math | 52-53 | `expectZero` strong-split endpoints | notes endpoint formula | yes |
| B2 | math | 64-67 | `Reduce[ForAll, dG/dR < 0, ...] = True` | monotonicity | yes (universal claim over the bound domain) |
| B3 | math | 69-73 | `Reduce[ForAll, dF/dR > 0, ...] = True` | monotonicity | yes |
| B4 | math | 82-86 | `PolynomialQuotientRemainder` proves `(1 - r^2)` divides numerator of `G_tr - G_flat` | branch ordering factorization | yes |
| B5 | math | 89-93 | `PolynomialQuotientRemainder` proves `(1 - r)` divides numerator of `F_flat - F_tr` | branch ordering / residual factorization | yes |
| B6 | math | 96-100 | `Reduce[ForAll, G_tr - G_flat > 0, ...] = True` | branch ordering (universal) | yes |
| B7 | math | 102-106 | `Reduce[ForAll, F_flat - F_tr > 0, ...] = True` | residual theorem (universal) | yes |
| B8 | math | 109-112 | boundary `expectZero` at r=0 and r=1 | endpoint sanity | yes |
| B9 | math | 119-129 | sample-point assertions | branch ordering | yes |

All listed assertions are anchored. None are tautological by construction: in each `expect_zero` the LHS is computed by SymPy/Mathematica from the defining `G_tr`, `F_tr` expressions and the RHS is a hand-typed closed form; a typo in either factor would surface as a nonzero residual. The Mathematica `Reduce[ForAll ...]` results are not formally certifiable from the saved output alone (Mathematica may return `True` cheaply or do real CAD), but the universal claim is at least the right shape, and the SymPy side independently nails down the same factorizations symbolically.

## Findings

### F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:87-91` (P_R)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:131-133` (P_1)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:136-139` (P_2)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:67-79` (P_R), `:98-109` (P_1), `:111-128` (P_2)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt:14, 20` (independent Mathematica derivation)

**What's wrong:**
The notes give explicit closed forms for the positive auxiliary polynomials `P_R` (in `dF_tr/dR`), `P_1` and `P_2` (in `F_flat - F_tr`). Five coefficients in the notes do not match what the scripts assert AND what both engines independently compute from `D[F_tr, R]` and `F_flat - F_tr`:

`P_R` notes lines 87-91 quote:
```
P_R
= 4 R^4 xi^3
  + 54 R^2 delta^2 xi + 90 R^2 delta xi^2 + 36 R^2 xi^3
  + 230 R delta^3 + 324 R delta^2 xi + 230 R delta xi^2
  + 81 delta^3 + 243 delta^2 xi + 243 delta xi^2 + 81 xi^3.
```
SymPy script (`scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:67-79`) uses `162 R delta^3` and `162 R delta xi^2` (not 230 and 230). Mathematica's machine-computed `D[fTr, r]` (mathematica output line 14) confirms `162*delta^3*r` and `162*delta*r*xi^2`.

`P_1` notes lines 131-133 quote:
```
P_1
= 18 R^2 delta^2 xi + 36 R^2 delta xi^2 + 22 R^2 xi^3
  + 81 R delta^3 + 248 R delta^2 xi + 99 R delta xi^2
  + 230 delta^3 + 423 delta^2 xi + 360 delta xi^2 + 99 xi^3,
```
SymPy script (`scripts/...:98-109`) uses `180 R delta^2 xi` (not 248) and `162 delta^3` (not 230). Mathematica's machine-computed `Factor[fFlat - fTr]` (mathematica output line 20) confirms `180*delta^2*r*xi` and `162*delta^3`.

`P_2` notes lines 136-139 quote (relevant fragment):
```
P_2 = ... + 237 R^2 xi^4 + ...
```
SymPy script (`scripts/...:111-128`) uses `220 R^2 xi^4` (not 237). Mathematica (line 20) confirms `220*r^2*xi^4`.

The paper stage card itself (`stage_046.tex`) does not enumerate these coefficients, so the card and the script are consistent; the disagreement is between the **notes** and the scripts. Because the notes are the authoritative derivation source per the audit protocol (item 4 under "What the paper claims"), a notes/script disagreement is a `paper_misalignment` finding.

**Why this matters:**
The script-side coefficients are double-verified: SymPy expands `dF_tr/dR` to exactly the form the script asserts (output line 20 shows the raw, un-touched SymPy output with `162` not `230`), and Mathematica independently computes `D[fTr, r]` from the same definition of `F_tr` and lands on the same `162` (output line 14). So the math is correct in the scripts and incorrect in the notes — the notes have arithmetic typos in three places that, if propagated to a careful reader trying to redo the algebra by hand, would not reproduce the script-asserted identities. The notes also feed the eventual `\stagefield`-formal write-up: the paper card today is terse, but if a future revision lifts these polynomials from the notes into the .tex, the error becomes a paper bug.

The positivity argument still works on either side (every coefficient is positive in both versions), so the qualitative sign claim is unaffected. The mismatch is purely a coefficient-of-record disagreement, not a sign or branch error.

**Required change:**
This is a `paper_misalignment` (notes_contradicts_script) finding. Codex must not silently rewrite either side. See the directive's `## Resolve before fix_loop` block.

**Verification:**
After the user picks a direction:
- If notes are corrected to match the scripts (the math-correct direction — both engines agree), no script change; verifier confirms the notes file now contains `162 R delta^3`, `162 R delta xi^2` in P_R; `180 R delta^2 xi`, `162 delta^3` in P_1; `220 R^2 xi^4` in P_2.
- If the user instead believes the notes' 230/248/237 are derived from a different (e.g., grouped or rescaled) form of `F_tr` than the scripts implement, the script's `F_tr` definition (`scripts/...:44-48`) needs revisiting and both engines must be re-run.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** a transliteration. The SymPy script verifies the polynomial factorizations by typing the candidate closed forms (`dG_expected`, `dF_expected`, `G_diff_expected`, `F_diff_expected`) and asserting their residual against the engine-computed expressions vanishes. The Mathematica script does the opposite: it never types `P_R`, `P_1`, or `P_2`; instead it uses `Reduce[ForAll, ..., < 0]`/`> 0]` over the entire bound domain `0 < r, xi < 1`, `delta > 0` to prove the inequality universally, and uses `PolynomialQuotientRemainder` to verify the structural factor `(1-r^2)` in `G_tr - G_flat` and `(1-r)` in `F_flat - F_tr`. Compare:

SymPy (lines 81-89):
```python
dG_expected = -36 * R * xi ** 2 * (delta + xi) / (2 * R ** 2 * xi + 9 * delta + 9 * xi) ** 2
dF_expected = (
    4 * xi * (2 * R * xi + 9 * delta + 9 * xi) * (2 * R ** 2 * xi + 9 * delta + 9 * xi) * P_R
    / (81 * (1 - xi) * (2 * R ** 2 * xi ** 2 + 9 * delta ** 2 + 18 * delta * xi + 9 * xi ** 2) ** 3)
)
expect_zero("dG_tr/dR formula", sp.simplify(dG_dR - dG_expected))
expect_zero("dF_tr/dR formula", sp.simplify(dF_dR - dF_expected))
expect_positive_coefficients("P_R", P_R, R, delta, xi)
```

Mathematica (lines 56-73):
```mathematica
dGdR = Together[D[gTr, r]];
dFdR = Together[D[fTr, r]];
...
dGSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, dGdR < 0], Reals];
...
dFSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, dFdR > 0], Reals];
```

These are genuinely different verification strategies on the same physical claim. No transliteration finding.

## Engine cross-check

Both engines compute the same `G_tr`, `F_tr`, `G_flat`, `F_flat`, the same factored derivatives, and the same factored branch differences. Numerical sample points produce identical rationals:

| Sample | SymPy `G_tr - G_flat` | Mma `G_tr - G_flat` | SymPy `F_flat - F_tr` | Mma `F_flat - F_tr` |
|---|---|---|---|---|
| 1 (R=1/4, xi=1/3, delta=1/2) | `225/8869` | `225/8869` | `38617837960/99381932001` | `38617837960/99381932001` |
| 2 (R=1/2, xi=1/2, delta=1/4) | `81/1736` | `81/1736` | `759648230/1473329763` | `759648230/1473329763` |
| 3 (R=3/4, xi=1/5, delta=2/3) | `91/21935` | `91/21935` | `5842146019415/70196178995856` | `5842146019415/70196178995856` |

`engines_agree: true`. Both engines also report the same internal coefficients (162, 180, 220) — disagreeing only with the notes, not with each other.

## Verdict justification

The script-side math is sound: both engines independently verify monotonicity, endpoint values, branch-difference factorizations, and universal positivity of the difference functions on the bound domain. The paper card's load-bearing claims (`\stagefield{Output}`: monotonicity + residual theorem) are fully exercised, with no tautological checks (every `expect_zero` is a hand-typed-closed-form-vs-engine-expansion comparison that could fail), no symbol-assumption errors (the script declares `xi, delta, R > 0, real` consistently with the physical setup `0 < xi < 1`, `delta > 0`, `R > 0`), no missing branch (the strong-split endpoint R=0 and the flat-branch endpoint R=1 are both tested, and Mathematica's `Reduce[ForAll, 0 < r < 1]` covers the open interior). Engines agree on numerical samples and on factored forms.

The only finding is that the **notes** contain three polynomial-coefficient typos (`P_R`: 230 vs script's 162 in two places; `P_1`: 248 and 230 vs script's 180 and 162; `P_2`: 237 vs script's 220) that disagree with what both engines independently compute. This is a `paper_misalignment` (notes_contradicts_script) that requires user resolution before any side is edited.

Attacks tried that failed:
- Trying to find a `simplify`-under-aggressive-assumptions that could hide a branch error: SymPy uses `positive=True, real=True` which matches the physical domain stated in the notes (`0 < xi < 1`, `delta > 0`, `R > 0`).
- Looking for sign-of-(1-xi) confusion: SymPy's factored form has `(xi-1)` in the denominator, but the typed expected form has `(1-xi)`; their ratio is `-1`, but the numerator carries the matching sign (`(R-1)(xi-1)` collapses to `(1-R)(1-xi)` after simplification). Both engines and the saved outputs are consistent.
- Checking whether the residual theorem `E_tr - E_flat = F_flat - F_tr` is independently asserted in the script: it is not, but it is identically true by the definitions in notes section 5 (`E_tr := R_target - F_tr`, `E_flat := R_target - F_flat`), so the substantive content reduces to `F_flat - F_tr > 0`, which is verified.
- Checking the strong-split endpoint by hand: G_tr(R=0) = 9 xi (delta+xi)/(9 delta + 9 xi) = xi (delta+xi)/(delta+xi) = xi (with delta+xi != 0). F_tr(R=0) = (9(delta+xi))^4 / [81 (1-xi) (9(delta+xi)^2)^2] = 9^4/(81 * 81 (1-xi)) = 1/(1-xi). Both correct.
- Checking whether the script's `expect_positive_coefficients` could rubber-stamp because of an accidentally-true expansion: I confirmed by reading the SymPy output (line 23) that the raw coefficient list is `[4, 54, 90, 36, 162, 324, 162, 81, 243, 243, 81]`; if the script's typed `P_R` were itself wrong, the prior `expect_zero` against `dF_expected` would have failed. So script's `P_R` is the SymPy-side `P_R`, not a hand transcription that escapes verification.

`verdict: findings`, `stop_cold: null`, `paper_alignment: partial`.

## Self-test notes

Variable-independence trap: every `sp.diff(EXPR, R)` and `D[expr, r]` is taken of an expression that explicitly contains `r` (both `G_tr` and `F_tr` depend on `R`), so derivatives are not identically zero. Symmetry/parity trap: no integrals appear, so this category does not apply. Trivial-case pre-check: substituted `R=0`, `R=1`, and three rational interior points by hand; all match the asserted closed forms. Path specification: not applicable — both scripts already exist. Paper round-trip: the proposed direction (correcting notes to match scripts) does not introduce a new `paper_misalignment` because the paper stage card never quotes the disputed coefficients; the notes correction is a notes-only edit and the user must approve it.
