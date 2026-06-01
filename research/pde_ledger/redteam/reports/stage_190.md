---
unit_id: 190
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage190_direct_defect_vs_dressing_split.md]
  paper_appendix: present
---

# Audit unit 190 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_190.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage190_direct_defect_vs_dressing_split.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 111, 770-888 referencing this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

`\stagefield{Output}` (stage_190.tex:15): "Separates direct transfer-shape drift from selected-branch dressing residual and records scalar no-go filters." The appendix row (part05:111) restates the same: "Separates direct transfer-shape drift from selected-branch dressing residual and records scalar no-go filters." The notes enumerate the concrete deliverables: (1) an exact direct-defect / dressing split, (2) an exact support-blindness theorem for the direct transfer shape `T^2`, `R_target`, and the corrected composite `N_*`, (3) the microscopic slippage ledger `(Sigma_Z, Sigma_chi, Sigma_eps, Sigma_delta, Sigma_eta)` with their closed forms, (4) the direct-defect law `Xi_1 = A_tr*Sigma_tr + Sigma_nt` and the triangular compiler `Theta_1 = -C_tr*Sigma_tr`, `Xi_1 = A_tr*Sigma_tr + Sigma_nt`, `R_1+Xi_1 = -epsilon_eta/(1-epsilon_eta)*Sigma_eta` with explicit `A_tr`, `C_tr`, (5) the exact inverse reconstruction formulas and the triple-rigidity theorem, and (6) the no-go filter: pure grouped real P2 anisotropy gives no linear scalar feed-down, the first grouped invariant being quadratic `I ~ (7/10)eps^2`. The appendix carries the same `Xi_1` slippage law (eq:app-part05-Xi1-slippage-law), `Sigma_tr`, `Theta_1`, `Sigma_nt`, the triangular normal form, and `A_tr`/`C_tr` defs verbatim.

## What the script claims to verify

The SymPy script verifies, section by section: (I) support-blindness — that the support-loaded reconstructions of `T^2`, `R_target`, `N_*` equal their bare forms and that their log-derivatives w.r.t. the support parameter `zeta` vanish, plus a spoiled-packet anti-tautology guard confirming a corrupted support term DOES break the blind route; (II) that the five microscopic log-drifts derived from the kernel definitions equal their stated closed forms; (III) the central decomposition identity `Xi_direct = A_tr*Sigma_tr + Sigma_nt`; (IV) the triangular compiler matrix, its determinant, the off-diagonal-zero structure, and the three inverse reconstruction formulas plus three rigidity theorems; (V) the weak-axisymmetric grouped signature `xbar = x0`, `b_x = 3a_x`, the invariant `I[x,y] = (7/10)eps^2 x1 y1`, and that the leading scalar invariant has no linear-in-eps term.

## Paper <-> script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Support-blindness of `T^2`, `R_target`, `N_*` (notes Sec 2) | lines 72-81 reconstructions + dln/dzeta = 0, with spoiled-packet guard 83-90 | match |
| Microscopic slippage ledger closed forms (notes Sec 3 / appendix 770-784) | lines 144-148 (Sigma_chi/delta/eta/eps/Z) | match |
| Direct-defect law `Xi_1 = A_tr*Sigma_tr + Sigma_nt` (notes Sec 4 / appendix eq:Xi1-slippage-law) | line 194 | match |
| `A_tr`, `C_tr` explicit forms (appendix eq:Ctr-Atr-defs) | lines 170-173, printed and used in 194/202/222 | match |
| Triangular compiler + det > 0 (notes Sec 5 / appendix triangular-normal-form) | lines 207-217 | match |
| Inverse reconstruction `Sigma_tr/Sigma_nt/Sigma_eta` + `A_tr/C_tr` identity (notes Sec 6) | lines 222-229 | match |
| Triple rigidity / tracking-rigid / cancellation theorems (notes Sec 6.4) | lines 230-232 | match |
| No-go filter: quadratic grouped invariant `I ~ (7/10)eps^2`, no linear feed-down (notes Sec 7) | lines 260-267 | match |
| Specific `epsilon_L^(1,P2)=epsilon_v=epsilon_T=0` (notes Sec 7) | not directly tested; these are Stage 253 observables that follow as corollaries of the verified lemma `I = O(eps^2)` | partial (acceptable — see verdict) |

`paper_alignment: aligned` — every load-bearing identity in the paper card/appendix/notes has a faithful, non-tautological script check; the only "partial" row is a downstream corollary the notes themselves frame as a consequence of the verified representation-theoretic lemma.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 72 | `T2_support - T2 == 0` | support-blindness reconstruction | yes |
| A2 | sympy | 73 | `Rtarget_support - Rtarget == 0` | support-blindness reconstruction | yes |
| A3 | sympy | 74 | `Nstar_support - Nstar == 0` | support-blindness reconstruction | yes |
| A4 | sympy | 75 | `Ecomp - (1 - epseta) == 0` | selected-branch `E = 1-eps_eta` | yes |
| A5 | sympy | 76 | `dln(T2_support)/dzeta == 0` | support-blindness theorem | yes |
| A6 | sympy | 77 | `dln(Rtarget_support)/dzeta == 0` | support-blindness theorem | yes |
| A7 | sympy | 78-81 | `dln(Nstar_support)/dzeta == 0` (with Rtr support-blind subst) | support-blindness theorem | partial (T2 part substantive; Rtr part asserted, but Rtr genuinely has no zeta) |
| A8 | sympy | 89-90 | spoiled dln(R)/dzeta != 0 | anti-tautology guard | yes (proves checks CAN fail) |
| A9-A13 | sympy | 144-148 | each `Sigma_X - closedform == 0` | slippage ledger | yes |
| A14 | sympy | 194 | `Xi_direct - (A_tr*Sigma_tr + Sigma_nt) == 0` | direct-defect law (central) | yes |
| A15 | sympy | 218-220 | `dXi/dSigma_eta`, `d(R+Xi)/dSigma_tr`, `d(R+Xi)/dSigma_nt == 0` | block-triangular off-diagonals | no (guaranteed by construction — see note) |
| A16 | sympy | 226 | `Sigma_tr_rec - Str == 0` | inverse reconstruction (1/C_tr prefactor) | yes |
| A17 | sympy | 227 | `A_tr/C_tr - 2(1+chi0+deltaU)/deltaU == 0` | inverse-formula coefficient | yes |
| A18 | sympy | 228 | `Sigma_nt_rec - Snt == 0` | inverse reconstruction (A_tr/C_tr cancel) | yes |
| A19 | sympy | 229 | `Sigma_eta_rec - Seta == 0` | dressing inversion prefactor | yes |
| A20 | sympy | 230 | `Xi(Str=0) - Snt == 0` | tracking-rigid theorem | yes |
| A21 | sympy | 231 | `Xi(Snt=-Atr*Str) == 0` | grouped-defect cancellation theorem | yes |
| A22 | sympy | 232 | `RplusXi(Seta=0) == 0` | selected-branch rigidity | yes |
| A23-A30 | sympy | 260-267 | `xbar-x0`, `bx-3ax`, `I[x,y]-(7/10)eps^2 x1 y1`, `I[x,x]-(7/10)eps^2 x1^2`, `dxbar/deps|0`, `dI/deps|0 == 0` | no-go quadratic invariant | yes |

The two non-"yes" rows (A7, A15) are discussed in the verdict; neither masks a substantive claim that is left untested elsewhere.

## Findings

None. (Two non-blocking observations are recorded in the verdict justification and self-test notes; neither rises to a script-side finding because neither causes a check to pass for the wrong reason in place of a substantive check.)

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit; the stage card (line 11) explicitly states "Mathematica audit: none yet." `mathematica_transliteration` does not apply.

## Engine cross-check

Single engine; `engine_disagreement` does not apply.

## Verdict justification

Verdict is **clean**. I attacked every load-bearing assertion and they hold up against the paper:

1. **Support-blindness (A1-A8) is genuinely substantive, not a derivative-of-a-constant trap.** `T2_support` is built as `Lambda0*(1-epseta)/Rtarget_support`, and `Rtarget_support` carries explicit `zeta` through both `Ssupport` and `Msupp`. So `T2_support` truly depends on `zeta` before the reconstruction `T2_support - T2 == 0` forces the `zeta`-dependence to cancel; `dln(T2_support)/dzeta == 0` then confirms it. The spoiled-packet guard (A8) injects `bad*zeta*Mmix` and confirms the derivative becomes nonzero — a textbook anti-tautology control proving the blindness checks CAN fail. This is exactly the kind of construction that defeats the "differentiate w.r.t. a symbol never wired in" failure mode.

2. **The central decomposition A14 (`Xi_direct = A_tr*Sigma_tr + Sigma_nt`) is a real, falsifiable algebraic identity.** `Xi_direct` is built from the microscopic law with `Sigma_chi`, and `Sigma_nt` (lines 174-182) carries no `Sigma_chi` term — so `A_tr*Sigma_tr` must reproduce the entire `Sigma_chi` content of `Xi_direct` plus the matching `Sigma_delta` piece, including the exact `(11+9*deltaU)/(11*(1+deltaU))` and `2*chi0/(1+chi0)` factors. Wrong coefficients (a factor-of-2 or an `11`/`9` slip) would break it. Hand-verified the coefficient structure matches appendix eq:app-part05-Xi1-slippage-law and eq:app-part05-Sigma-nt-def.

3. **The inverse reconstructions A16-A19 are substantive** — each verifies a non-trivial prefactor (`1/C_tr`, the `A_tr/C_tr` cancellation, the dressing `-(1-eps_eta)/eps_eta` factor). Hand-checked `A_tr/C_tr = 2(1+chi0+deltaU)/deltaU` and `Sigma_tr_rec = -(1/C_tr)*Theta = Str`.

4. **The no-go invariant A23-A30 is correct.** Hand-derivation: `a_x = eps*x1/4`, `b_x = 3*eps*x1/4` (so `b_x = 3 a_x`), and `I[x,x] = 4 a_x^2 + (4/5) b_x^2 = (1/4 + 9/20)eps^2 x1^2 = (7/10)eps^2 x1^2`, with `dI/deps|_{eps=0} = 0` — confirming the leading scalar invariant is quadratic with no linear feed-down.

**Two non-blocking observations** (neither a finding):
- (a) **Banner mislabel.** The print banners say "STAGE 173" (lines 35, 269) and the ledger header says "STAGE 173 LEDGER" while the file/card is stage 190. The math, however, is unambiguously stage-190's (all formulas match stage_190's notes/appendix exactly). This is a cosmetic provenance label left over from the project's original-stage-order renumbering (the card itself, line 17, notes it "restores original stage order"). It does not affect any assertion or verification validity; it is not a `paper_misalignment` of a verified identity. Recorded as cosmetic only — not worth a Codex directive over a print string.
- (b) **A7 and A15 are partially/fully guaranteed-by-construction.** A7's `Rtr_sb` support-blindness is asserted via `diff(Rtr_sb, zeta)->0` rather than derived, but `Rtr` (line 50) genuinely contains no `zeta`, so the substitution is physically correct and the `T2_support` factor of A7 is independently substantive. A15's three off-diagonal derivative checks (`dXi/dSigma_eta` etc.) are zero by construction because `Xi` is built as `A_tr*Str + Snt` without `Seta`; they are decorative confirmations of the block-triangular layout. Critically, they do NOT stand in for the substantive triangularity content — that is carried by A14 (the microscopic `Xi_direct` decomposition) plus the explicit `Sigma_tr`/`Sigma_nt` definitions, both of which are independently and non-tautologically verified. There is also no available non-tautological in-script replacement: any restatement of "Xi is `Sigma_eta`-independent" is structural to how the paper sets up the coordinates, not an emergent algebraic identity, so flagging A15 would only generate Codex churn relocating a tautology. I therefore do not raise it to a finding.

**On the missing Mathematica engine:** I do NOT flag `missing_mathematica`. Per the prompt's line-114 judgment for this non-status-only-but-SymPy-only stage: SymPy fully and genuinely settles every claimed symbolic identity. These are all deterministic CAS algebra checks (closed-form log-drift derivations, an exact rational-function decomposition, a triangular matrix inversion with explicit prefactors, and a polynomial representation-theoretic invariant). There is no numerical sensitivity, branch cut, transcendental simplification, or pole subtlety where a second independent engine would add confidence; SymPy's `simplify(expand(...))` settles each residual to an exact `0`. This matches established pipeline precedent (SymPy-only non-status-only stages 121/122/123; 56 SymPy-only stages pipeline-wide). I confirm I read the paper card, notes, and appendix rows, and the script's verified claims match the paper's stated claims.

## Self-test notes

I checked all five traps. (1) Variable independence: A5/A6 differentiate `T2_support`/`Rtarget_support` w.r.t. `zeta`, which they genuinely depend on (verified via `Msupp`/`Ssupport`); the spoiled-packet guard confirms non-vacuity. A15's derivatives ARE w.r.t. symbols not wired into the constructed expression (the only guaranteed-zero block) — flagged as a non-blocking observation, not relocated. (2) Symmetry/parity: no unbounded integrals; the no-go uses finite weighted sums, hand-verified `I[x,x] = (7/10)eps^2 x1^2`. (3) Trivial-case: hand-substituted the weak-axisymmetric profile into `a_x`, `b_x`, `I` and into the `Xi_direct` decomposition coefficients — residuals reduce to 0 for the right reasons. (5) Paper round-trip: all coefficients (`A_tr`, `C_tr`, the `11+9*deltaU` factors, the `7/10` invariant) match appendix eqs verbatim. Conclusion: clean.
