---
unit_id: 020
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 020 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_020.tex`
- notes: `(none)` (no `notes/stages/moving_throat_pde_stage020_*.md` present)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 62)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.txt`

## What the paper claims

The card's `\paragraph{Output.}` states: "Stage~020 exports the weak-axisymmetric parent packet \eqref{eq:stage020-wall-slopes}--\eqref{eq:stage020-residual-xi}." The deliverables are: (i) the slope-packet definitions `D_{01}=δK_Σ-B_{01}-Z_{01}`, `D_{21}=-(δM_Σ+B_{21}+Z_{21})`, `D_{41}=-(B_{41}+Z_{41})` (eq:stage020-denominator-slopes); (ii) the live even gates `K_1=D_{21}+D_{01}/9` and `H_even=D_{41}-(2/3)D_{21}-D_{01}/27` (eq:stage020-gates) plus the obstruction `Ξ_1=N_{01}/N_0-D_{01}/D_0` (eq:stage020-xi); (iii) the exact wall-slope solve `δK_Σ=B_{01}+Z_{01}+27(B_{41}+Z_{41})`, `δM_Σ=-(B_{21}+Z_{21})+3(B_{41}+Z_{41})` (eq:stage020-wall-slope-solve); and (iv) the residual obstruction after that solve `Ξ_1=N_{01}/N_0-27(B_{41}+Z_{41})/(K_Σ-B_0-Z_0)` (eq:stage020-residual-xi). The claim is that the two even gates fix `(δK_Σ, δM_Σ)` uniquely and that the first-order normalization defect collapses to the single scalar `Ξ_1`. The appendix row marks the stage `\StatusExactClosure{}` — "Weak-axisymmetric packet exported from the parent throat-action bundle."

## What the script claims to verify

The SymPy script defines the gates exactly as the card does and verifies: (1) the even-gate Jacobian determinant `∂(K1,H_even)/∂(dKSigma,dMSigma) = 1/27` (line 39); (2) the linear `solve` of `K1=0, H_even=0` returns the card's closed forms for `dKSigma` and `dMSigma` (lines 41,48-49); (3) on that compensated branch `D01=27(B41+Z41)`, `D21=-3(B41+Z41)`, `D41=-(B41+Z41)` (lines 50-52); and (4) the compensated obstruction equals `N01/N0-27(B41+Z41)/(KSigma-B0-Z0)` (line 53). The Mathematica script (M2-M5) reproduces the same four claims via an independent `Det`, an independent `Solve` (with a uniqueness/solution-count check), and residual comparisons against the card's closed forms.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Slope-packet/denominator defs D01,D21,D41 (eq:stage020-denominator-slopes) | py L27-29 / wl L27-30 (definitions, exercised through all checks) | match |
| Even gates K1, H_even (eq:stage020-gates) | py L31-32 / wl L31-32; det=1/27 at py L39 / wl L39-43 | match |
| Ξ_1 = N01/N0 - D01/D0 (eq:stage020-xi) | py L33 / wl L33 | match |
| Wall-slope solve δK_Σ, δM_Σ (eq:stage020-wall-slope-solve) | py L41,48-49 / wl L46-66 | match |
| Compensated D01,D21,D41 | py L50-52 / wl L69-88 | match |
| Residual Ξ_1 (eq:stage020-residual-xi) | py L53 / wl L90-100 | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `assert_zero(coeff_matrix.det() - 1/27)` | gates invertibility (eq:gates) | yes |
| A2 | sympy | 48 | `assert_zero(sol[dKSigma] - (B01+Z01+27(B41+Z41)))` | wall-slope solve δK_Σ | yes |
| A3 | sympy | 49 | `assert_zero(sol[dMSigma] - (-(B21+Z21)+3(B41+Z41)))` | wall-slope solve δM_Σ | yes |
| A4 | sympy | 50 | `assert_zero(D01_comp - 27(B41+Z41))` | compensated D01 | yes |
| A5 | sympy | 51 | `assert_zero(D21_comp + 3(B41+Z41))` | compensated D21 | yes |
| A6 | sympy | 52 | `assert_zero(D41_comp + B41 + Z41)` | compensated D41 | yes |
| A7 | sympy | 53 | `assert_zero(Xi1_comp - (N01/N0 - 27(B41+Z41)/(KSigma-B0-Z0)))` | residual Ξ_1 (eq:residual-xi) | yes |
| A8 | math | 41-43 | `Det[gateJacobian] - 1/27 === 0 else Exit[1]` | gates invertibility | yes |
| A9 | math | 48-50 | `Length[branchSolve] - 1 === 0 else Exit[1]` | uniqueness of solve | yes |
| A10 | math | 60-66 | dKSigma/dMSigma residual vs closed forms `=== 0 else Exit[1]` | wall-slope solve | yes |
| A11 | math | 77-87 | D01/D21/D41 residual `=== 0 else Exit[1]` | compensated deficits | yes |
| A12 | math | 96-99 | Xi1 residual vs closed form `=== 0 else Exit[1]` | residual Ξ_1 | yes |

All rows trace to a specific card deliverable; none are orphaned.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is not a transliteration. It defines the gates from the same physical premises (which it must, since the gate forms are the paper's deliverable, not free choices) but then performs a genuinely independent computation: `Det[gateJacobian]` (wl L39) versus SymPy's `coeff_matrix.det()` (py L39) are independent linear-algebra calls; `Solve[{K1==0,HEven==0},{dKSigma,dMSigma}]` (wl L46) is an independent solve, and the `.wl` adds a uniqueness gate `Length[branchSolve]-1===0` (wl L47-50) that the `.py` does not have. Both engines then compare their *own* solve output against the paper's stated closed forms (wl L52-66 vs py L48-49) — they verify against the paper, not against each other's intermediate algebra. The residual checks (wl L69-99) substitute the engine's own `branchRules` into the deficit/obstruction expressions. This is two independent derivations converging on the card's closed forms, not an echo.

## Engine cross-check

Both engines report all residuals exactly 0 and `STATUS: PASS`. SymPy output (lines 10-18) prints the symbolic closed forms `dKSigma_from_even_gates = B01 + 27*B41 + Z01 + 27*Z41`, `dMSigma_from_even_gates = -B21 + 3*B41 - Z21 + 3*Z41`, compensated `D01=27*B41+27*Z41`, `D21=-3*B41-3*Z41`, `D41=-B41-Z41`, and `Xi1_on_compensated_branch = (27*N0*(B41+Z41)+N01*(B0-KSigma+Z0))/(N0*(B0-KSigma+Z0))`, which rationalizes to `N01/N0 - 27(B41+Z41)/(KSigma-B0-Z0)`. Mathematica output (M2-M5, all residuals = 0) confirms the same identities. The two engines agree.

## Verdict justification

`clean`. I read the card, the appendix row (no notes file exists for this stage), and both scripts. Attacks tried and failed: (1) the determinant `1/27` was recomputed by hand from the gate partials `∂K1/∂dKSigma=1/9, ∂K1/∂dMSigma=-1, ∂H_even/∂dKSigma=-1/27, ∂H_even/∂dMSigma=2/3` → `(1/9)(2/3)-(-1)(-1/27)=2/27-1/27=1/27`, matching A1/A8 — not tautological since the gates and the literal `1/27` come from independent definitions. (2) The `solve` assertions (A2/A3, A10) compare an in-script `sp.solve`/`Solve` result against the paper's closed forms; neither closed form is fed into the solve, so the check can genuinely fail. (3) The compensated-branch deficits and `Ξ_1` (A4-A7, A11-A12) were verified by hand on the branch: `K1=−3(B41+Z41)+27(B41+Z41)/9=0` and `H_even=−(B41+Z41)+2(B41+Z41)−(B41+Z41)=0`, consistent with output lines 16-17. (4) Symbol domains: `N0`, `KSigma-B0-Z0` are declared `nonzero` (py L23-24) / `!= 0` (wl L25), which is exactly what the `Ξ_1` and `D_0` denominators require — no over- or under-assumption. Every card deliverable maps to a substantive, non-tautological, dual-engine check. No findings.

## Self-test notes

Checked: (1) variable-independence of the Jacobian partials — `D41` correctly has no `dKSigma`/`dMSigma` dependence (its rows are zero in the Jacobian), and the nonzero partials reproduce `det=1/27` by hand. (2) No unbounded integrals here, so the parity trap does not apply. (3) Trivial-case substitution: setting `B41=Z41=0` collapses the compensated branch to `D01=D21=D41=0` and `Ξ_1=N01/N0`, consistent with all assertions. (4) Denominator-nonzero assumptions are present and justified. Conclusion: no script-side defect, full paper alignment.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 11 values checked, 0 misaligned.

Enumerated every RESULT/deliverable value the two scripts emit (from script source + saved outputs) and located each in the card. There are no notes file for this stage; the `.tex` card is the sole prose carrier and it states all of them.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `D01 = δK_Σ - B01 - Z01` | py L27 / out L3; wl L28 | stage_020.tex:23 | MATCH |
| `D21 = -(δM_Σ + B21 + Z21)` | py L28 / out L4; wl L29 | stage_020.tex:24 | MATCH |
| `D41 = -(B41 + Z41)` | py L29 / out L5; wl L30 | stage_020.tex:25 | MATCH |
| even-gate Jacobian det `= 1/27` | py L39; wl L39 / mout L2 | stage_020.tex:30-31 (gate forms; det is the unique-solvability of these gates) | MATCH |
| `δK_Σ = B01 + Z01 + 27(B41+Z41)` | py L48 / out L10; wl L53 | stage_020.tex:42 | MATCH |
| `δM_Σ = -(B21+Z21) + 3(B41+Z41)` | py L49 / out L11; wl L54 | stage_020.tex:44 | MATCH |
| compensated `D01 = 27(B41+Z41)` | py L50 / out L13; wl L71 | stage_020.tex:42-43 (= δK_Σ - B01 - Z01) | MATCH |
| compensated `D21 = -3(B41+Z41)` | py L51 / out L14; wl L72 | stage_020.tex:44 (= -(δM_Σ+B21+Z21)) | MATCH |
| compensated `D41 = -(B41+Z41)` | py L52 / out L15; wl L73 | stage_020.tex:25 | MATCH |
| `Ξ_1 = N01/N0 - D01/D0` (live form) | py L33 / out L8; wl L33 | stage_020.tex:35 | MATCH |
| residual `Ξ_1 = N01/N0 - 27(B41+Z41)/(K_Σ-B0-Z0)` | py L53 / out L18; wl L93 | stage_020.tex:49-50 | MATCH |

INTERNAL scaffolding (accounted for, no finding): `K1`/`H_even` printed intermediate expanded forms (out L6-7), `Check_D21_plus_D01_over_9 = 0` (out L16), `Check_D41_minus_2D21_over_3_minus_D01_over_27 = 0` (out L17), per-check residual prints (`M2..M5 ... residual = 0`, mout L2-13), and `STATUS: PASS` flags. These are verification residuals/consistency confirmations, not deliverable values.

All 11 emitted deliverable values reconcile against the card. Note: the compensated-branch `Xi1` was printed by SymPy in an unrationalized form `(27*N0*(B41+Z41)+N01*(B0-KSigma+Z0))/(N0*(B0-KSigma+Z0))` (out L18); this rationalizes exactly to the card's `N01/N0 - 27(B41+Z41)/(K_Σ-B0-Z0)` (the assertion at py L53 enforces this), so it is a MATCH, not a mismatch.
