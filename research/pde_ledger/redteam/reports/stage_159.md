---
unit_id: 159
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage159_hybrid_outlet_projection.md]
  paper_appendix: present
---

# Audit unit 159 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_159.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage159_hybrid_outlet_projection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_159}` row at line 1352; no row-level summary beyond that)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage159_hybrid_outlet_projection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage159_hybrid_outlet_projection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage159_hybrid_outlet_projection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage159_hybrid_outlet_projection_mathematica_audit.txt`

## What the paper claims

The paper card (stage_159.tex) is the linear projection of the co-evolving Family-1 defect onto the compensated Robin-mixed outlet. Its load-bearing claim is the quoted result:

> "Even preservation forces \(\delta\mathfrak g=0\), \(\delta\kappa_W=0\), leaving \(\Delta_Q=-9\sigma_*\delta\gamma_W/(1-\sigma_*)\)."

The accompanying notes (section 3-6) make this concrete with four deliverables on the compensated canonical branch (rho_R,* = 4 sigma_*, kappa_W,* = 1/3, gamma_W,* = 1/9):
(a) the linearized outlet defects delta_E2 = (delta_C - 9 sigma_* delta_kappa_W)/(27(1-sigma_*)), delta_E4 = (5 delta_C - 72 sigma_* delta_kappa_W)/(243(1-sigma_*)), Delta_Q = (delta_C - 27 sigma_* delta_gamma_W)/(3(1-sigma_*));
(b) the mouth-gain to outlet-loading transport with the cancellation delta_C := delta_rho_R - 4 delta_sigma_W = -4 Xi Sigma0_can delta_R (i.e., delta_Sigma0 drops out);
(c) the canonical-even gate (delta_E2 = delta_E4 = 0) forcing delta_C = delta_kappa_W = 0 on a nontrivial sigma_* != 0 branch;
(d) the collapse Delta_Q = -9 sigma_*/(1 - sigma_*) delta_gamma_W.
The checklist in the card also requires showing (i) deviations are taken about the renormalized canonical point, (ii) even-preservation is imposed before reading the odd defect, (iii) tangent motion on the parent compensation family is what cancels delta_Sigma0.

## What the script claims to verify

Both scripts (SymPy and Mathematica) construct the nonlinear hybrid-outlet primitives L0, L2, L4, chi, E2, E4; substitute rho = 4 sigma0 + drho, sigma = sigma0 + dsigma, kappa = 1/3 + dkappa, gamma = 1/9 + dgamma (the compensated canonical branch); linearize via Taylor expansion; then assert the linearized chi-1, E2, E4 equal the closed-form expressions (drho - 4 dsigma - 27 sigma0 dgamma)/(3(1-sigma0)), (drho - 4 dsigma - 9 sigma0 dkappa)/(27(1-sigma0)), (5 drho - 20 dsigma - 72 sigma0 dkappa)/(243(1-sigma0)). They then verify mouth-gain transport delta_C = Xi delta_Sigma0 + (-1)(... ) reduces to -4 Xi Sigma0_can delta_R (the delta_Sigma0 cancellation), substitute Xi = 4 sigma_*/Sigma0_can to get delta_C = -16 sigma_* delta_R, solve the canonical-even linear system to confirm the unique solution delta_C = dk = 0 (and print determinant = -27 sigma0), and finally check that with delta_C = 0 the residual Delta_Q + 9 sigma0 dgamma/(1 - sigma0) vanishes. Both scripts exit 0 with all assertions passing.

## Paper to script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| delta_E2 = (delta_C - 9 sigma_* delta_kappa_W)/(27(1-sigma_*)) | sympy L79 / wl L57 `expect_zero("delta E2 formula", E2_lin - E2_expected)` | match |
| delta_E4 = (5 delta_C - 72 sigma_* delta_kappa_W)/(243(1-sigma_*)) | sympy L80 / wl L58 `expect_zero("delta E4 formula", E4_lin - E4_expected)` | match |
| Delta_Q = (delta_C - 27 sigma_* delta_gamma_W)/(3(1-sigma_*)) | sympy L78 / wl L56 `expect_zero("delta chi formula", delta_chi - delta_chi_expected)` (delta_chi = chi_lin - 1 = Delta_Q) | match |
| delta_C = delta_rho_R - 4 delta_sigma_W cancels delta_Sigma0; = -4 Xi Sigma0_can delta_R | sympy L92 / wl L71 `expect_zero("deltaC mouth transport", ...)` | match |
| sigma_*-substitution: delta_C = -16 sigma_* delta_R | sympy L95-98 / wl L75-78 | match |
| Canonical-even gate forces delta_C = delta_kappa_W = 0 | sympy L103-108 / wl L84-88 (Solve + assertion + printed determinant -27 sigma0) | match |
| Collapse: Delta_Q = -9 sigma_*/(1 - sigma_*) delta_gamma_W | sympy L114-115 / wl L93-94 | match |
| Deviations about renormalized canonical point | Linearization is around (rho=4 sigma0, sigma=sigma0, kappa=1/3, gamma=1/9) | match |
| Even-preservation imposed before reading odd defect | "FINAL REDUCED DEFECT" section substitutes deltaC -> 0 then asserts the residual | match |

Paper alignment: aligned. No paper deliverable is missing a script-side check, and no script check is orphaned relative to the paper. The script's banner string says "STAGE 142" in both engines, but this is a cosmetic display string, not a math claim, and the docstring/filename correctly say stage 159.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 78 | `expect_zero("delta chi formula", delta_chi - delta_chi_expected)` | Delta_Q closed form | yes |
| A2 | sympy | 79 | `expect_zero("delta E2 formula", E2_lin - E2_expected)` | delta_E2 closed form | yes |
| A3 | sympy | 80 | `expect_zero("delta E4 formula", E4_lin - E4_expected)` | delta_E4 closed form | yes |
| A4 | sympy | 92 | `expect_zero("deltaC mouth transport", deltaC_expr - deltaC_expected)` | delta_Sigma0 cancellation, delta_C = -4 Xi Sigma0_can delta_R | yes |
| A5 | sympy | 95-98 | `expect_zero("sigma_star substitution", ...)` | sigma_* parameter substitution | yes |
| A6 | sympy | 107-108 | `if sol != [{deltaC: 0, dk: 0}]: raise` | canonical-even gate forces delta_C = delta_kappa_W = 0 | yes (sympy `solve` returns unique solution; determinant -27 sigma0 is also printed) |
| A7 | sympy | 115 | `expect_zero("final Delta_Q + 9 sigma* dgamma /(1-sigma*)", ...)` | final collapse Delta_Q = -9 sigma_* dgamma/(1-sigma_*) | yes |
| B1 | mathematica | 56 | `expectZero["delta chi formula", deltaChi - deltaChiExpected]` | Delta_Q closed form | yes |
| B2 | mathematica | 57 | `expectZero["delta E2 formula", e2Lin - e2Expected]` | delta_E2 closed form | yes |
| B3 | mathematica | 58 | `expectZero["delta E4 formula", e4Lin - e4Expected]` | delta_E4 closed form | yes |
| B4 | mathematica | 71 | `expectZero["deltaC mouth transport", ...]` | delta_Sigma0 cancellation | yes |
| B5 | mathematica | 75-78 | `expectZero["sigma_star substitution", ...]` | sigma_* substitution | yes |
| B6 | mathematica | 86-88 | `If[sol =!= {{deltaC -> 0, dk -> 0}}, fail[...]]` | canonical-even gate | yes |
| B7 | mathematica | 94 | `expectZero["final Delta_Q + 9 sigma* dgamma /(1-sigma*)", ...]` | final collapse | yes |

All assertions are non-tautological: each compares a linearization derived from the underlying nonlinear primitives (L0, L2, L4, chi, E2, E4) against a closed-form expression stated in the notes. The closed-form expressions are computed and substituted independently from the nonlinear ones, so an algebra error in either side would surface as a nonzero residual.

## Findings

None. The audit produced zero findings.

Trivial-case sanity checks I performed (mentally, no script execution):
- chi at the background point: chi(rho=4 sigma_0, sigma=sigma_0, kappa=1/3, gamma=1/9) = 3(1 - 9 sigma_0 * 1/9)/(3 - 4 sigma_0 + sigma_0) = 3(1 - sigma_0)/(3(1 - sigma_0)) = 1. So Delta_Q = chi - 1 vanishes at zeroth order, as expected.
- partial chi / partial gamma at the point: -27 sigma_0/(3(1-sigma_0))) = -9 sigma_0/(1-sigma_0), matching the final-collapse coefficient.
- partial chi / partial rho at the point = 1/(3(1-sigma_0)); partial chi / partial sigma at the point = -4/(3(1-sigma_0)). These give the (drho - 4 dsigma) combination in the numerator of the closed form. Consistent.
- Mouth transport algebraically: deltaC_expr = Xi delta_Sigma0 - 4 * (-Xi) * (-(1/4) delta_Sigma0 - Sigma0_can delta_R) = Xi delta_Sigma0 - Xi delta_Sigma0 - 4 Xi Sigma0_can delta_R = -4 Xi Sigma0_can delta_R. Confirms the delta_Sigma0 cancellation that the paper card check (iii) asks for.
- Canonical-even linear system determinant: 1*(-72 sigma_0) - (-9 sigma_0)*5 = -72 sigma_0 + 45 sigma_0 = -27 sigma_0. Matches the printed determinant and the notes section 5.

## Independent-derivation check (Mathematica)

The Mathematica script mirrors the SymPy script section-for-section and uses analogous variable names. Section banners, assertion labels, the linearize helper's structure, and the order of checks are essentially identical. However, the algebraic engines invoked are genuinely independent:
- Linearization: SymPy uses `sp.series(..., var, 0, 2).removeO()` and Poly filtering; Mathematica uses `D[f, var]` and a Taylor sum at zero. Different machinery.
- Simplification: SymPy uses `sp.simplify(sp.expand(...))`; Mathematica uses `FullSimplify[Together[Expand[...]], Assumptions -> $Assumptions]`. Different simplifiers.
- Solve: `sp.solve([...], [...], dict=True)` vs `Solve[{...}, {...}, Reals]`. Different solvers.
- Determinant: `sp.Matrix(...).det()` vs `Det[{...}]`. Different routines.

The setup of the nonlinear expressions (chi, E2, E4) and the closed-form expected expressions are stated identically because they are the paper's stated quantities — both engines must work with the same physical objects. What matters is that the simplification machinery that confirms residual = 0 is independent, and it is. I am NOT raising a `mathematica_transliteration` finding: the script choreography is similar by necessity (same paper-stated identities), but the engines independently reduce the residual.

## Engine cross-check

Both engines produce identical zero residuals on every named assertion (sympy output L13-15, L20-21, L32; mathematica output L13-18, L23-26, L37). Determinant printed by both: `-27*sigma0`. Solve / sp.solve both return the unique solution {deltaC -> 0, dk -> 0}. Engines agree.

## Verdict justification

The script's assertions exercise the paper's stated deliverables one-for-one. The linearized E2, E4, Delta_Q expressions are independently constructed by Taylor-expanding the full nonlinear chi, E2, E4 around the compensated canonical branch and compared to the closed forms from the notes — non-tautological. The mouth-transport check explicitly demonstrates the delta_Sigma0 cancellation that the paper card's check (iii) requires. The canonical-even Solve returns the unique solution {deltaC=0, dk=0} and the determinant -27 sigma_0 is printed, confirming the non-degeneracy on a loaded branch. The final-collapse residual Delta_Q + 9 sigma_0 dgamma/(1-sigma_0) is verified zero. Both engines pass with identical residuals. The only cosmetic blemish is the section banner string "STAGE 142" appearing in both scripts — this does not affect any math, the docstring and filename correctly identify stage 159, and per the audit policy cosmetic display strings are not findings. Attacks attempted that failed: (a) checking whether the linearized formulas might be tautological — they are not, because the closed forms are stated independently of the construction; (b) checking the canonical-even Solve for σ_0 = 0 degeneracy — the script prints the determinant -27 sigma_0 making this transparent, and the paper explicitly assumes a nontrivial loaded branch sigma_* != 0; (c) checking whether the Mathematica is just a port — engines independently reduce residuals via different simplification stacks. Paper-aligned, internally consistent, both engines verify, outputs are fresh (sympy output 2026-05-11T12:47; script 2026-05-11 11:56; mathematica output 2026-05-11T13:17; script 2026-05-11 11:56). Verdict: clean.

## Self-test notes

Checked: (1) variable independence in `linearize` — each derivative is taken w.r.t. one of {drho, dsigma, dkappa, dgamma} which the substituted expressions genuinely depend on, so no zero-derivative trap. (2) Symmetry/parity — not applicable (no integrals). (3) Trivial-case pre-check — at zeroth-order substitution chi = 1 exactly (verified manually above), so Delta_Q = chi - 1 vanishes; the linear residual against the closed form must therefore have no constant piece, which is consistent with the assertion. (4) Path specifications — both scripts present at canonical paths. (5) Paper round-trip — closed-form expressions in the script's `*_expected` definitions match the boxed formulas in notes sections 3 and 5 character-for-character (modulo notation: drho/dsigma/dkappa/dgamma in the script = delta_C-decomposition / delta_kappa_W / delta_gamma_W in the notes after the identification delta_C = drho - 4 dsigma). No new paper_misalignment introduced.
