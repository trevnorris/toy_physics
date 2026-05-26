---
unit_id: 042
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage042_rank2_selected_mode_normalization.md
  paper_appendix: present
---

# Audit unit 042 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_042.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage042_rank2_selected_mode_normalization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 62 + `\input` at line 202)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.txt`

## What the paper claims

The paper stage card (`stage_042.tex`) carries `\stagefield{Output}` (line 41): "The rank--2 eigenvector \eqref{eq:app-stage042-evec-rank2} and overlap formulas \eqref{eq:app-stage042-z-overlap}--\eqref{eq:app-stage042-s-overlap}." Boxed in the card body are three claims: (1) the selected eigenvector ratio `e1/e0 = [m(q-r)+r*xi] / [delta+xi - m*q*(q-r)]`, (2) the outgoing/mixed overlap `(z.e_-)^2/z0^2 = [delta+(1+qr)*xi]^2 / D_(q,r)`, and (3) the source overlap `(s.e_-)^2/s0^2 = [delta+(1+rt)*xi - m*(q-r)*(q-t)]^2 / D_(q,r)`, with `D_(q,r) = [delta+xi-m*q*(q-r)]^2 + [m*(q-r)+r*xi]^2`. The notes file expands these to six deliverables matching the six numbered sections of the script: (1) eigenvector ratio, (2) overlap formulas and the generalized normalization function `F_(q,r,t) = (z-overlap)(s-overlap)/(1-xi)`, (3) tracking-support collapse `r=q` back to Stage-23, (4) source-tied split-U specialization `F_src(xi,delta;m,R_U)` with `r=t=sqrt(lam0)`, `q=sqrt(lam0)*R_U`, (5) flat-U recovery `F_src(R_U=1) = F_flat`, and (6) first-order deformation coefficients `H_n^(src)` and `H_F^(src)` about `R_U=1`. The notes are authoritative on intent; the stage card .tex is terse but does not contradict them. The appendix row (line 62) summarizes the unit as "Selected eigenvector and source/loading overlaps with distinct outgoing, support, and source directions."

## What the script claims to verify

The sympy script's docstring (lines 4-16) enumerates six checks that map 1-to-1 onto the notes' six sections: (1) inserting Stage-24 `n_req` makes both rows of the eigenvalue equation collapse to the closed-form `e1/e0`; (2) the constructed Z and S overlaps reduce to the named expected closed forms and combine into `F_(q,r,t)`; (3) `F_(q,r,t)|_{r=q}` collapses to the Stage-23 normalization function; (4) the source-tied substitution `q->sqrt(lam0)*R_U, r->sqrt(lam0), t->sqrt(lam0)` reproduces an independently constructed `F_src`; (5) `F_src(R_U=1) = F_flat`; (6) the R_U-derivatives at `R_U=1` of `n_src` and `F_src/F_flat` match named closed forms `H_n_expected`, `H_F_expected`, and the truncated Taylor series of `n_src`, `F_src/F_flat` to first order in `eps = R_U - 1` match the linearized forms. The Mathematica script mirrors all six sections one-for-one.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Selected eigenvector ratio `e1/e0` (paper eq, notes Sec. 1) | A1, A2 (row1, row2 vs ratio_expected) | match |
| Z-overlap closed form (paper eq, notes Sec. 2) | A3 (Z_overlap vs Z_expected) | match |
| S-overlap closed form (paper eq, notes Sec. 2) | A4 (S_overlap vs S_expected) | match |
| Generalized `F_(q,r,t)` (notes Sec. 3 only; paper card omits) | A5 (F_general vs F_expected) | match (but tautological — see F1) |
| Tracking-support collapse to Stage 23 (notes Sec. 4) | A6 (F_track vs F_stage23) | match |
| Source-tied `F_src` (notes Sec. 5) | A7 (F_src_direct vs F_src constructed) | match |
| Flat-U recovery (notes Sec. 6) | A8 (F_src(R_U=1) vs F_flat) | match |
| First-order `H_n^(src)`, `H_F^(src)`, linear expansions (notes Sec. 6) | A9-A12 | match |

The paper card's `\stagefield{Output}` enumerates only the first three deliverables explicitly; the remaining four deliverables (F_(q,r,t), tracking collapse, source-tied `F_src`, flat-U recovery, first-order deformation) appear in the notes (authoritative on intent per audit protocol). Per audit-prompt order of authority (notes file is item 4), the additional script checks do not constitute `paper_missing_script_claim`. Paper alignment is **aligned**.

## Assertion inventory

| #   | Script       | Line    | Form                                                                | Exercises which paper claim?            | Anchored to claim? |
|-----|--------------|---------|---------------------------------------------------------------------|------------------------------------------|--------------------|
| A1  | sympy        | 72      | `expect_zero(ratio_row1 - ratio_expected)`                          | Selected eigenvector ratio              | yes                |
| A2  | sympy        | 73      | `expect_zero(ratio_row2 - ratio_expected)`                          | Selected eigenvector ratio (other row)  | yes                |
| A3  | sympy        | 89      | `expect_zero(Z_overlap - Z_expected)`                               | Z-overlap closed form                   | yes                |
| A4  | sympy        | 90      | `expect_zero(S_overlap - S_expected)`                               | S-overlap closed form                   | yes                |
| A5  | sympy        | 101     | `expect_zero(F_general - F_expected)`                               | F_(q,r,t) (notes Sec. 3)                | no — TAUTOLOGICAL  |
| A6  | sympy        | 112     | `expect_zero(F_track - F_stage23)`                                  | Tracking collapse (notes Sec. 4)        | yes                |
| A7  | sympy        | 130     | `expect_zero(F_src_direct - F_src)`                                 | Source-tied F_src (notes Sec. 5)        | yes                |
| A8  | sympy        | 138     | `expect_zero(F_src(R_U=1) - F_flat)`                                | Flat-U recovery (notes Sec. 6)          | yes                |
| A9  | sympy        | 152     | `expect_zero(H_n_src - H_n_expected)`                               | H_n^(src) (notes Sec. 6)                | yes                |
| A10 | sympy        | 165     | `expect_zero(H_F_src - H_F_expected)`                               | H_F^(src) (notes Sec. 6)                | yes                |
| A11 | sympy        | 171-174 | `expect_zero(series(n_src, eps, ord 2) - n_linear)`                 | First-order linear expansion of n_src   | yes                |
| A12 | sympy        | 175-178 | `expect_zero(series(F_ratio, eps, ord 2) - F_linear)`               | First-order linear expansion of F_src   | yes                |
| B1  | mathematica  | 56-57   | `expectZero` row1, row2 vs ratioExpected                            | Selected eigenvector ratio              | yes                |
| B2a | mathematica  | 80      | `expectZero` zOverlap - zExpected                                   | Z-overlap                               | yes                |
| B2b | mathematica  | 81      | `expectZero` sOverlap - sExpected                                   | S-overlap                               | yes                |
| B2c | mathematica  | 82      | `expectZero` fGeneral - fExpected                                   | F_(q,r,t)                               | no — TAUTOLOGICAL  |
| B3  | mathematica  | 94      | `expectZero` tracking collapse                                       | Tracking collapse                       | yes                |
| B4  | mathematica  | 112     | `expectZero` source-tied specialization                              | Source-tied F_src                       | yes                |
| B5  | mathematica  | 120     | `expectZero` F_src(rU=1) - F_flat                                    | Flat-U recovery                         | yes                |
| B6  | mathematica  | 146-147 | `expectZero` H_n, H_F                                                | First-order coefficients                | yes                |
| B7  | mathematica  | 148-155 | `expectZero` linear-series expansions                                | First-order linear expansions           | yes                |

## Findings

### F1 — tautological_check

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py:92` (definition of `F_general`) and `:101` (assertion)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl:70` (definition of `fGeneral`) and `:82` (assertion)

**What's wrong:**

The script computes
```
F_general = simplify(Z_expected * S_expected / (1 - xi))           # sympy line 92
```
and then asserts `F_general - F_expected == 0` (sympy line 101), where
```
F_expected = simplify(
    (delta + (1 + q*r)*xi)**2
    * (delta + (1 + r*t)*xi - m*(q-r)*(q-t))**2
    / ((1 - xi) * D_qr**2))
```
Now `Z_expected = (delta + (1+q*r)*xi)**2 / D_qr` and `S_expected = (delta + (1+r*t)*xi - m*(q-r)*(q-t))**2 / D_qr` (lines 81-82). Substituting these into `F_general` gives the numerator of `F_expected` divided by `(1-xi) * D_qr**2`, i.e. **literally** `F_expected`. The residual `F_general - F_expected` is identically zero by algebraic substitution — no simplification or rearrangement of physical quantities is involved. The assertion **cannot fail** regardless of whether `Z_expected`, `S_expected`, or `D_qr` are physically correct.

The Mathematica script repeats the same construction at lines 70-75 / 82.

The substantive checks for the Z and S overlaps (A3 and A4) are non-tautological and remain valid. Only the F_general consequence is tautological.

**Why this matters:**

The notes section 3 lists `F_(q,r,t)` as a distinct deliverable (the product of the two overlaps, divided by `(1-xi)`). The current check does not exercise that product structure; it only verifies an algebraic rearrangement of `a/b * c/b = a*c/b^2`. To exercise the deliverable, `F_general` must be built from the eigenvector-constructed forms (`Z_overlap`, `S_overlap`, not `Z_expected`, `S_expected`), so the assertion forces sympy/Mathematica to chain the simplifications established in A3/A4 (B2a/B2b) through the product.

**Required change:**

In the SymPy script, change line 92 from
```python
F_general = sp.simplify(Z_expected * S_expected / (1 - xi))
```
to
```python
F_general = sp.simplify(Z_overlap * S_overlap / (1 - xi))
```
This routes the F_general check through the constructed-from-eigenvector overlap forms, making the residual `F_general - F_expected` a real simplification test (it passes only after `Z_overlap` reduces to `Z_expected` and `S_overlap` reduces to `S_expected`).

In the Mathematica script, change line 70 from
```
fGeneral = FullSimplify[zExpected sExpected/(1 - xi), Assumptions -> $Assumptions];
```
to
```
fGeneral = FullSimplify[zOverlap sOverlap/(1 - xi), Assumptions -> $Assumptions];
```

**Verification:**

After Codex applies the change, re-running both scripts must still print `F_general - expected = 0` / `fGeneral - fExpected = 0` and exit 0. The new construction is provably equal to the old one (via A3, A4 / B2a, B2b), so no math identity changes. The saved output transcripts should be identical at the assertion-result level after fresh runs.

## Independent-derivation check (Mathematica)

The Mathematica script's structure is parallel to the SymPy script: same six-section partition, same intermediate variable names (`nReq <-> n_req`, `ratioRow1 <-> ratio_row1`, `dQR <-> D_qr`, `qTiedQR <-> q_tied_qr`, `a1`, `b1`, `dSrc`, `fSrc`, `fRatio`, `hNSrc`, `hFSrc`, `hNExpected`, `hFExpected`, `nLinear`, `fLinear`), same assertion targets, same evaluation order. Three corresponding examples:

- SymPy `n_req = sp.simplify((xi*(delta + xi) - m*(delta + (1 + q**2)*xi)) / (delta + (1 + r**2)*xi - m*(q - r)**2))` (lines 55-58) vs Mathematica `nReq = FullSimplify[(xi (delta + xi) - m (delta + (1 + q^2) xi))/(delta + (1 + r^2) xi - m (q - r)^2), Assumptions -> $Assumptions]` (lines 44-47) — same formula, same intermediate name.
- SymPy `Z_overlap = sp.simplify((1 + q*ratio_expected)**2 / (1 + ratio_expected**2))` (line 78) vs Mathematica `zOverlap = FullSimplify[(1 + q ratioExpected)^2/(1 + ratioExpected^2), Assumptions -> $Assumptions]` (line 62) — same construction from ratioExpected, same intermediate name.
- SymPy source-tied block (lines 117-121: `q_tied_qr = lam0 * R_U`, `a1`, `b1`, `D_src`, `F_src`) vs Mathematica (lines 98-105: `qTiedQR = lambda0 rU`, `a1`, `b1`, `dSrc`, `fSrc`) — identical structure.

This is parallel structure rather than independent first-principles re-derivation. However, the script's *claims* are themselves a set of named algebraic identities between rational functions ("`Z_overlap` equals `Z_expected`", etc.). Verifying such identities is inherently engine-bound: there is one input, two normalization paths to test agreement of. Both engines exercise their own simplification kernels (`sp.simplify(sp.expand(...))` at SymPy line 38 vs `FullSimplify[Together[Expand[...]], Assumptions -> $Assumptions]` at Mathematica line 28) — these are independent algorithms with independent normal-form choices. Evidence: the printed forms of `H_F^(src)` differ between engines (SymPy line 132-136 prints combined-fraction form `4*(18*d^2*m + ...)/((9*d+11*x)*(9*d^2+...))`; Mathematica line 62 of the wl output prints partial-fraction form `(4*(xi/(delta + 11*xi/9) + (2*delta*m)/((2*xi^2)/9 + (delta + xi)^2)))/9`), yet both simplify the *residual* to identically 0. The cross-check therefore has real content. **Not flagged as `mathematica_transliteration`**, but noted that the per-section structural parallelism is at the upper bound of acceptable for a second-engine audit; a cleaner future revision would solve the eigenvalue equation in Mathematica via `Solve[]` instead of asserting `n_req`'s closed form.

## Engine cross-check

Both engines report PASS on every assertion (sympy exit code 0 at output line 155; Mathematica exit code 0 at output line 74). Residuals all print as 0. Selected final forms:

- `e1/e0`: SymPy output line 27-31 prints `-(m*q - m*r + r*xi) / (-delta + m*q^2 - m*q*r - xi)`; Mathematica output line 19 prints `(m*(q - r) + r*xi) / (delta + m*q*(-q + r) + xi)` — these are identical after a sign flip of both numerator and denominator.
- `H_n^(src)`: both engines print `-4*m*xi / (9*delta + 11*xi)` (SymPy lines 127-129, Mathematica line 61).
- `F_src(R_U=1) - F_flat = 0` in both (SymPy line 121, Mathematica line 55).
- `F_general - expected = 0` in both — but this assertion is tautological by construction (see F1).

Engines agree: true.

## Verdict justification

The audit holds up on substance. The eigenvector / overlap / collapse / source-tied / flat-U / first-order-deformation claims are all anchored to non-tautological checks in both engines. The paper card .tex is terse (Output line lists only three of the six deliverables) but the notes file is authoritative and matches the script's six sections precisely; alignment is therefore `aligned`. Engine cross-check is genuine despite the parallel structure of the Mathematica script, because the two simplification kernels are independent and produce visibly distinct normal forms before reducing the asserted residuals to zero. The one defect is the F_general / fGeneral assertion (A5 / B2c), which is algebraically guaranteed by construction and adds no verification power beyond A3 + A4 (B2a + B2b). Fixing this is a one-line change in each script. No downstream propagation risk; no UNFIXABLE math. Verdict: `findings` with one low-severity tautology.

## Self-test notes

- **Variable independence**: For `sp.diff(n_src, R_U)` (line 148) and `sp.diff(F_ratio, R_U)` (line 155), confirmed that `n_src` (line 142-145) and `F_src` (line 121) explicitly depend on `R_U`; derivatives are non-trivial and reduce at `R_U=1` to the printed nonzero closed forms. The trap "derivative of an R_U-independent expression giving 0 trivially" does not apply.
- **Symmetry/parity**: No unbounded-domain integrals; n/a. R_U=1 substitution does not pinch any denominator (`dSrc|_{rU=1} = (delta+xi)^2 + lam0*xi^2 > 0` under `delta,xi > 0`).
- **Trivial-case substitution**: At R_U=1, `n_src` evaluates to `xi*(delta+xi)/(delta+(1+lam0)*xi) - m = G_flat - m`, the zeroth-order term of `n_linear` — independently asserted at A11 alongside the first-order coefficient. At r=q, `D_qr = (delta+xi)^2 + q^2*xi^2`, matching Stage-23 denominator — used in A6.
- **Paper round-trip (v2)**: The proposed F1 fix only changes which intermediate is multiplied into `F_general` (Z_overlap and S_overlap instead of Z_expected and S_expected). The resulting check still verifies the same paper/notes claim `F_(q,r,t) = (z-overlap)(s-overlap)/(1-xi)` and uses the same paper-side constants (`lam0 = 2/9` is unchanged, no new literals). No new `paper_misalignment` is introduced.
- **Soft note (not a formal finding)**: `eps` is declared `positive=True` in SymPy (line 49) and `eps > 0` in Mathematica `$Assumptions` (line 38), but the notes section 6 explicitly state that the constructive split-U branch has `R_U < 1`, i.e. `eps < 0`. Since the assertions A11/A12/B7 are polynomial in `eps` of degree at most 1 (after series truncation), the positivity assumption does not affect the residual being tested. No finding filed.
