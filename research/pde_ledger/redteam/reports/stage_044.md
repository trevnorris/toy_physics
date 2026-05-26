---
unit_id: 044
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T06:54:07Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage044_continuum_selected_rank2_closure.md
  paper_appendix: present
---

# Audit unit 044 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_044.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage044_continuum_selected_rank2_closure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 66, 206, 312)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.txt`

Freshness: sympy output mtime > sympy script mtime; mathematica output mtime > mathematica script mtime. Both outputs are fresh.

## What the paper claims

The paper card (stage_044.tex) states two boxed deliverables and two diagnostic special surfaces. (1) The physical softening depth satisfies the exact quadratic
`xi^2 + B_cont xi + C_cont = 0`
with
`B_cont = delta - M_mix*(1 + t^2 R_U^2) - M_supp*(1 + t^2 R_phi^2)`
and
`C_cont = -delta*(M_mix + M_supp) + t^2 M_mix M_supp (R_U - R_phi)^2`
(equation `eq:app-stage044-Ccont`, lines 25-30 of the .tex). (2) Once the physical root `xi_phys in (0,1)` is selected, the normalization gate is `R_target = F_cont(xi_phys)` (lines 31-36). The card also names two diagnostic surfaces (minimal-kernel source-tied at `R_phi=1`, and interference-matched tracking at `R_phi=R_U`), distinguished by whether the direction-splitting term `(R_U - R_phi)^2` survives (line 37). The notes file supplies the underlying Stage-24 `n_req` formula (notes §2, lines 67-90), the explicit `F_cont` and `D_cont` expressions (notes §4, lines 146-158), the `+` root-selection rule (notes lines 116-123, "the `+` branch is selected because it reduces continuously to xi_phys = 0 when M_mix = M_supp = 0"), the exact mismatch penalty `lambda_0 M_mix M_supp (R_U - R_phi)^2` (notes §6, lines 220-232), and the two surface limits `R_phi=1` (source-tied) and `R_phi=R_U` (tracking, where the branch collapses to `G_q(xi,delta) = xi(delta+xi) / [delta + (1 + q^2) xi]` with total loading `M_mix + M_supp`; notes §5, lines 170-214). The notes use `lambda_0 = t^2 = 2/9` and the paper card uses `t^2` directly; both refer to the same constant. Part appendix line 66 records the stage as exact closure with summary "Physical softening depth fixed by a quadratic and selected normalization gate."

## What the script claims to verify

Both engines build `n_req^(cont)` from the Stage-24 support theorem with `q = sqrt(lambda0) R_U`, `r = sqrt(lambda0) R_phi`, `lambda0 = 2/9` (sympy lines 56-61; mathematica lines 44-48), then assert: (1) the numerator of `n_req^(cont) - M_supp` factors as `9*(xi^2 + B_cont xi + C_cont)` with hand-typed `B_cont`, `C_cont` matching the paper (sympy line 72; mathematica line 56); (2) the closed-form `xi_phys = (-B_cont + sqrt(B_cont^2 - 4 C_cont))/2` vanishes at zero load (sympy line 78; mathematica line 61); the Mathematica engine adds an independent `Solve[branchEq == 0, xi]`-based derivation of `xi_phys` and verifies it matches the closed form (mathematica lines 65-77); (3) `F_cont` and `D_cont` are constructed and verified at a literal third slice `R_phi = 2` against an independently-typed template (sympy lines 99-109; mathematica lines 95-103); (4) the source-tied limit `R_phi = 1` yields the source-tied `n_source`, `F_source` matching independently-typed templates (sympy lines 113-132; mathematica lines 107-124); (5) the tracking limit `R_phi = R_U` collapses `n_req` to `G_q - M_mix` and `F_cont` to a reduced form with the `(R_U-R_phi)`/`(R_U-1)` factors vanishing (sympy lines 136-146; mathematica lines 128-138); (6) the xi-constant coefficient of `branch_eq / 9` equals `-delta*(M_mix + M_supp) + lambda0 M_mix M_supp (R_U - R_phi)^2`, i.e. the mismatch penalty lives in `C_cont` (sympy line 154; mathematica line 148).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Quadratic `xi^2 + B_cont xi + C_cont = 0` with stated `B_cont`, `C_cont` | sympy L65-72 / math L49-56 (`quadratic branch equation` PASS) | match |
| Physical root `xi_phys` as `+` branch reducing to 0 at zero load | sympy L74-78 (`zero-load root`); math L58-77 (closed form + `Solve` cross-check) | match |
| Normalization gate `R_target = F_cont(xi_phys)` (build `F_cont`, `D_cont`) | sympy L83-109 (build + third-slice); math L81-103 | match |
| Minimal-kernel source-tied surface (`R_phi = 1`) | sympy L113-132; math L107-124 | match |
| Interference-matched tracking surface (`R_phi = R_U`) collapse to one-direction `G_q` with total loading | sympy L136-146; math L128-138 (`tracking collapse of n_req`, `tracking F collapse` PASS) | match |
| Mismatch penalty `+t^2 M_mix M_supp (R_U-R_phi)^2` in `C_cont` | sympy L150-154; math L143-148 | match |
| Domain restriction `xi_phys in (0,1)` | not algebraically checked (physical condition, not a closure obligation) | n/a |

All six paper deliverables are covered. The unchecked `xi_phys in (0,1)` domain restriction is a downstream physical condition, not an algebraic claim of the stage closure. Paper alignment: `aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 72 | `expect_zero(branch_eq - 9*branch_expected)` | quadratic decomposition (B_cont, C_cont) | yes (numerator of n_req - M_supp built by an algebraically distinct path from B_cont/C_cont) |
| A2 | sympy | 78 | `expect_zero(xi_phys.subs(Mmix=0, Msupp=0))` | root selection rule | yes (closed-form sqrt evaluated at zero load) |
| A3 | sympy | 109 | `expect_zero(F_lit - F_lit_expected)` at Rphi=2 | F_cont construction | yes (drift-protective: F_cont definition and template are independently typed; a typo in F_cont would not propagate to F_lit_expected) |
| A4 | sympy | 131 | `expect_zero(n_source - n_source_expected)` at Rphi=1 | source-tied surface | yes (independent template at Rphi=1) |
| A5 | sympy | 132 | `expect_zero(F_source - F_source_expected)` at Rphi=1 | source-tied surface | yes (independent template at Rphi=1; (R_U-1)^2 form written explicitly) |
| A6 | sympy | 138 | `expect_zero(n_track - (G_q - Mmix))` | tracking surface collapse | yes (genuine collapse: 2-direction form reduces to 1-direction G_q minus Mmix) |
| A7 | sympy | 146 | `expect_zero(F_track - F_track_expected)` | tracking surface F-collapse | yes (target form drops factors that vanish at Rphi=R_U) |
| A8 | sympy | 154 | `expect_zero(C_from_branch_eq - C_expected)` | mismatch penalty in C_cont | yes (extract xi-constant coeff from numerator(n_req-Msupp), compare to hand-typed C_cont; distinct paths) |
| B1 | math | 56 | `expectZero[branchEq - 9 branchExpected]` | quadratic decomposition | yes |
| B2 | math | 61 | `expectZero[xiPhys /. {mMix->0, mSupp->0}]` | zero-load reduction | yes |
| B3 | math | 77 | `expectZero[xiPhysSolve - xiPhys]` (Solve-based independent root) | root selection, independent algebraic route | yes (genuinely independent — Solve derives the root from branchEq, then compared to manually-typed `(-B + Sqrt[disc])/2`) |
| B4 | math | 103 | `expectZero[fLit - fLitExpected]` at rPhi=2 | F_cont construction | yes (drift-protective) |
| B5 | math | 123 | `expectZero[nSource - nSourceExpected]` at rPhi=1 | source-tied surface | yes |
| B6 | math | 124 | `expectZero[fSource - fSourceExpected]` at rPhi=1 | source-tied surface | yes |
| B7 | math | 137 | `expectZero[nTrack - (gQ - mMix)]` | tracking collapse | yes |
| B8 | math | 138 | `expectZero[fTrack - fTrackExpected]` | tracking F-collapse | yes |
| B9 | math | 148 | `expectZero[cFromBranchEq - cExpected]` | mismatch penalty | yes |

Every assertion is anchored to a stated paper or notes claim, and each is non-tautological in the sense that a coefficient typo in any one of `n_req^(cont)`, `B_cont`, `C_cont`, `F_cont`, `D_cont`, `G_q`, or the source-tied / tracking templates would produce a nonzero residual on at least one of these checks. The slice checks (A3/B4, A4/B5, A5/B6) are drift-protective: they pair the constructed expression evaluated at a literal `R_phi` value against an independently-typed template at the same literal. The strong audit comes from A1/B1 (extracting `B_cont`, `C_cont` from the numerator of `n_req - M_supp`), B3 (Solve-based independent route), A6/B7 (genuine collapse), and A8/B9 (extracting `C_cont` from the polynomial coefficient of `branch_eq`).

## Independent-derivation check (Mathematica)

Both scripts share the same template definitions (the formulas have to be the same — they are the paper-stated `B_cont`, `C_cont`, `F_cont`, `D_cont`). However, the Mathematica script provides a genuine independent algebraic route for the root selection at lines 65-77:

```
xiRoots = xi /. Solve[branchEq == 0, xi];
xiPhysSolve = SelectFirst[xiRoots, TrueQ[FullSimplify[(# /. {mMix -> 0, mSupp -> 0}) === 0, Assumptions -> $Assumptions]] &];
...
expectZero["xiPhys solve match", xiPhysSolve - xiPhys];
```

The SymPy script forms `xi_phys = (-B_cont + sqrt(Delta_disc))/2` by hand. The Mathematica script does the same but ALSO solves the quadratic via `Solve` and picks the zero-load-vanishing root, then asserts the two forms agree (mathematica output line 16-17 shows the `Solve`-derived form printed verbatim; line 17-18 shows `xiPhys solve match = 0`, `PASS: xiPhys solve match`). A typo in either the manual coefficients of `B_cont`/`C_cont` or in the manual `(-B + sqrt(disc))/2` would be detected by this independent route.

The corresponding SymPy assertion A1 (`branch_eq - 9*branch_expected == 0`) is itself algebraically non-trivial: `branch_eq` is the numerator of `n_req - M_supp` (computed via `together` and `as_numer_denom`), while `branch_expected` is built from the explicit hand-typed `B_cont`, `C_cont`. The two are constructed by algebraically distinct paths.

Conclusion: the `.wl` is not a transliteration. It has at least one truly independent verification route (Solve) on top of the shared algebraic framework. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit identical assertion-level residuals (all 0). Cosmetic differences:

- SymPy's printed `n_req^(cont)` (output lines 9-14) factors out a `-1` from both numerator and denominator relative to Mathematica's form (output line 9): SymPy shows `(2 M_mix R_U^2 xi + 9 M_mix delta + 9 M_mix xi - 9 delta xi - 9 xi^2) / (2 M_mix R_U^2 - 4 M_mix R_U R_phi + 2 M_mix R_phi^2 - 2 R_phi^2 xi - 9 delta - 9 xi)`, Mathematica shows `(xi*(delta + xi) - mMix*(delta + xi + (2*rU^2*xi)/9))/(delta - (2*mMix*(rPhi - rU)^2)/9 + xi + (2*rPhi^2*xi)/9)`. After multiplying SymPy's numerator and denominator by `-1` and dividing both by `9`, the expressions are identical.
- `xi_phys` printed forms differ in presentation (Mathematica also prints a `Solve`-derived form on output line 16); the assertion `xiPhys solve match = 0` confirms the two forms are algebraically equal.
- `F_cont`, `n_source`, `F_source`, `n_track`, `F_track` printed forms are identical modulo internal simplifier conventions.

Final residual values:
- SymPy: `quadratic branch equation = 0`, `zero-load root = 0`, `third-slice F at Rphi=2 = 0`, `source-tied n = 0`, `source-tied F = 0`, `tracking collapse of n_req = 0`, `tracking F collapse = 0`, `mismatch penalty in C coefficient = 0`.
- Mathematica: same residuals all PASS, plus `xiPhys solve match = 0` and `PASS: xiPhys solve match`.

No engine disagreement.

## Verdict justification

The unit's two boxed paper deliverables (the quadratic `xi^2 + B_cont xi + C_cont = 0` and the normalization gate `R_target = F_cont(xi_phys)`) and the two diagnostic surfaces (source-tied at `R_phi=1`, tracking at `R_phi=R_U`) are all faithfully exercised by substantive, non-tautological assertions in both engines. The Mathematica script provides a genuinely independent algebraic route via `Solve[branchEq == 0, xi]` (lines 65-77). The mismatch penalty `+t^2 M_mix M_supp (R_U-R_phi)^2` is anchored to the actual `C_cont` coefficient extracted from `Poly(branch_eq, xi).nth(0) / 9` (rather than a variable-renaming tautology). All assertions pass with zero residual; both engines exit 0; outputs are fresh.

Attacks tried that the script survived:
- (i) checked the `t^2` vs `lambda_0` convention discrepancy between paper card and notes — both equal `2/9`, no inconsistency, and both scripts use `lambda0 = 2/9` consistently.
- (ii) hand-expanded `n_req - M_supp` and matched coefficient-by-coefficient against `9*(xi^2 + B_cont xi + C_cont)` — agrees with script.
- (iii) verified `xi_phys` at zero load reduces to 0 by hand (`B_cont = delta`, `C_cont = 0`, `Delta_disc = delta^2`, root = 0).
- (iv) verified tracking-surface limit `R_phi = R_U` does collapse `D_cont` to `(delta+xi)^2 + lambda0 R_U^2 xi^2` and the `F_cont` factor `(R_U - R_phi)(R_U - 1)` vanishes; matches `F_track_expected`.
- (v) verified mismatch penalty `+lambda0 M_mix M_supp (R_U - R_phi)^2` appears in `C_cont` by direct expansion of `n_req - M_supp` and matching the xi^0 coefficient.
- (vi) cross-checked the `Solve`-derived `xiPhysSolve` formula in the Mathematica output (line 16) against the manually-typed closed form: the `Solve` form is `xiPhys = (-9 delta + 9 mMix + 9 mSupp + 2 mSupp rPhi^2 + 2 mMix rU^2)/18 + Sqrt[...]/18`, which equals `(-B_cont)/2 + sqrt(B_cont^2 - 4 C_cont)/2` after multiplying both numerator and denominator by 9; matches the manual form.
- (vii) checked that the source-tied `F_source_expected` uses `(R_U - 1)^2` explicitly (script line 128), which is algebraically equal to `(R_U - Rphi)*(R_U - 1)|_{Rphi=1}` but written in a distinct simplified form, so a typo in F_cont's `(R_U - Rphi)*(R_U - 1)` factor would not silently propagate.
- (viii) checked the part appendix row 66 — summary text matches the paper card and notes; no contradiction.

Verdict: `clean`. The five prior-pass findings (mathematica transliteration, mismatch-penalty tautology, redundant tracking total-loading, insufficient F_cont verification, unused `sigma0`) have all been applied and the resulting script is paper-aligned and audit-tight.

## Self-test notes

- Variable independence: no `sp.diff` or `D[...]` calls in this unit; no derivative-trap risk.
- Symmetry/parity: no integrals over unbounded domains; no parity-based vanishing claims; n/a.
- Trivial-case pre-check: re-confirmed `xi_phys` at zero load reduces to 0 (B_cont=delta>0, C_cont=0, Delta_disc=delta^2, (-delta+delta)/2=0); re-confirmed tracking-surface F collapse at Rphi=RU eliminates the (R_U-Rphi)(R_U-1) cross term.
- Path specifications: no missing-script findings, no new paths proposed.
- Paper round-trip: paper card, notes, and part appendix row all agree on the two boxed deliverables, the two diagnostic surfaces, and the mismatch-penalty form; the scripts encode the same content. The `t^2` (paper) vs `lambda_0` (notes) notation difference is harmless since both denote the same numerical constant `2/9` and both scripts hard-code `lambda0 = 2/9`.
