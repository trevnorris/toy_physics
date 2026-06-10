---
unit_id: 234
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 234 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_234.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row 80; narrative eqs lines 905-938)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.py`
- mathematica: `/var/projects/toy_projects` NOT used — `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage234_direct_branch_observable_static_gate_and_the_two_observable_kill_test_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Direct static gate: on the rigid-mouth branch, the static same-charge test is decided by the realized pair $(R_{\rm target},\epsilon_\eta)$, with tracking data factored out." The derivation ledger states the card "shows two cancellations, reduces the static gate to $R_{\rm target}$ and $\epsilon_\eta$, and states the two-observable kill test." The notes (the authoritative carrier) enumerate the deliverables: (1) the exact finite quotient chart $(q_{\rm tr},q_{\rm nt},q_\eta)$ in direct observables plus its exact inverse; (2) the first-order linearization; (3) the triangular compiler $\Theta_1=\delta\ln R_{\rm tr}$, $\Xi_1=-\delta\ln R_{\rm target}-c_\eta\delta\ln\epsilon_\eta$, $\mathcal R_1=\delta\ln R_{\rm target}$ with $c_\eta:=\epsilon_{\eta,*}/(1-\epsilon_{\eta,*})$; (4) the exact $R_{\rm tr}$ cancellation $\partial_{\delta\ln R_{\rm tr}}\Xi_1=0$; (5) the rigid-mouth reduction $q_{\rm nt}=\Xi_1$; (6) the robust ($0.367930328492646$) and nonempty ($0.737619063660757$) two-observable strip ceilings; (7) the canonical families (pure-target, pure-dressing, balanced minimal-norm). The appendix (eqs 901-938) reproduces the chart, compiler, and gate verbatim.

## What the script claims to verify

The SymPy script and the Mathematica script verify the same seven deliverables: (M1) the finite chart forward/inverse roundtrip; (M2) the first-order drifts via $\partial_\epsilon|_0$ (SymPy) / `Series` coefficient (Mathematica); (M3) the triangular compiler $\Theta_1,\Xi_1,\mathcal R_1$ plus the inverse $E_1$ map; (M4) $\partial_{r1}\Xi_1=0$ (SymPy `diff`; Mathematica `D` plus `FreeQ`); (M5) rigid-mouth reduction $\Xi_1|_{r1=0}=q_{\rm nt}^{(1)}|_{r1=0}$ and the finite $q_{\rm nt}$ form; (M6) the two-observable strip edges and ceiling ordering $0<{\rm robust}<{\rm nonempty}$; (M7) pure-target, pure-dressing, and balanced minimal-norm families. All assertions zero out in both saved outputs.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Finite chart + exact inverse | M1 roundtrip (both engines) | match |
| First-order linearization | M2 (both) | match |
| Triangular compiler + inverse | M3 (both) | match |
| $R_{\rm tr}$ cancellation $\partial\Xi_1=0$ | M4 (both) | match |
| Rigid-mouth $q_{\rm nt}=\Xi_1$ | M5 (both) | match |
| Two ceilings + strip form | M6 (both) | match |
| Canonical families (3) | M7 (both) | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 55-57 | `assert_zero(qX_roundtrip - q_X)` x3 | chart+inverse | yes |
| A2 | sympy | 88-90 | `assert_zero(qX1 - expected)` x3 | linearization | yes |
| A3 | sympy | 105-110 | compiler + inverse E1 | triangular map | yes |
| A4 | sympy | 122 | `assert_zero(sp.diff(Xi1, r1))` | $R_{\rm tr}$ cancel | yes (r1 genuinely enters qnt1/qtr1, cancels) |
| A5 | sympy | 126-135 | rigid finite form + reductions | rigid-mouth | yes |
| A6 | sympy | 148-163 | ceiling ordering/gap + strip widths | two ceilings | yes |
| A7 | sympy | 177-200 | three families + Lagrange solve | canonical families | yes |
| B1 | math | 98-100 | `expectLogZero` roundtrips | chart+inverse | yes |
| B2 | math | 120-122 | `expectZero` linearization | linearization | yes |
| B3 | math | 135-138 | compiler + inverse | triangular map | yes |
| B4 | math | 147-148 | `D[xi1Raw,r1]` + `FreeQ` | $R_{\rm tr}$ cancel | yes |
| B5 | math | 156-158 | rigid reductions | rigid-mouth | yes |
| B6 | math | 173-177 | strip edges + ordering | two ceilings | yes |
| B7 | math | 184-199 | families + `Minimize` | canonical families | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

Genuinely independent, not a transliteration. The `.wl` uses materially different machinery: (1) linearization via `linCoeff` = `Coefficient[Normal[Series[expr,{eps,0,1}]],eps,1]` (lines 66-69) vs SymPy's `sp.diff(qX_eps, eps).subs(eps,0)` (lines 82-84); (2) the balanced minimal-norm family via built-in `Minimize[{R1^2+E1^2, -R1-cEta*E1==xi},{R1,E1},Reals]` (line 187-190) vs SymPy's hand-built Lagrange-multiplier system `sp.solve(eqs,...)` (lines 184-193); (3) M4 adds an independent `FreeQ[FullSimplify[xi1Raw],r1]` structural check (line 148) that SymPy does not perform; (4) custom log-residual machinery `PowerExpand`/`cleanLogResidual` and `stripConditional` for the `Minimize` `Piecewise` output. Variable naming and overall stage structure are parallel (expected for the same physics), but the derivation routes differ — this passes the second-engine bar.

## Engine cross-check

Both outputs report all checks PASS / `= 0`. Final compiler forms agree: SymPy `Xi1 = -E1*eps_eta_star/(1 - eps_eta_star) - R1` (output line 16) vs Mathematica `Xi1 = (E1*epsEtaStar)/(-1 + epsEtaStar) - R1` (output line 44) — algebraically identical ($-1/(1-x)=1/(-1+x)$). Balanced family agrees (SymPy line 33; Mathematica `Piecewise` line 95 reduces to the same). No disagreement.

## Verdict justification

Clean. I read the card, notes, and appendix before the scripts, and the script-side bottom-line assertions match the paper's `Output` and all seven enumerated notes/appendix deliverables exactly. Attacks tried and failed: (1) the suspected variable-independence trap at M4 is NOT vacuous — r1 appears with nonzero coefficient `B_star` in `qnt1` and `-C_star` in `qtr1`, and the `(B_star/C_star)` weighting is exactly what makes the r1 terms cancel; output line 21 confirms `d Xi1/d(ln R_tr)=0` is a real cancellation, and Mathematica's added `FreeQ` corroborates structurally; (2) no hardcoded "answer" smuggling — the two ceilings are carried-forward corridor inputs, and the script only verifies their ordering/gap and the strip-width algebra, both of which reconcile with the notes; (3) no tautology — each roundtrip/compiler check substitutes an independently-defined inverse/perturbation and is free to fail; (4) symbol domains are physically correct (`positive` for radii/dressings, `nonzero` for `vareps`). Both engines present, independent, agreeing; outputs fresh (22:10 > scripts 22:05/22:07); no stale self-labels (both scripts say "Stage 234" throughout; the pass-1 "Stage 251" notes-title drift is no longer present).

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Compiler: $\Xi_1=-R_1-c_\eta E_1$ | py:112 / wl:166; out py:16, wl:44 | tex:13,15; appendix:920,932; notes:273 | MATCH |
| $c_\eta=\epsilon_{\eta,*}/(1-\epsilon_{\eta,*})$ | py:86 / wl:82 | appendix:927; notes:139 | MATCH |
| $\Theta_1=\delta\ln R_{\rm tr}$, $\mathcal R_1=\delta\ln R_{\rm target}$ | py:101,103; out:15,17 | appendix:918,922; notes:145,159 | MATCH |
| $R_{\rm tr}$ cancellation $\partial_{r1}\Xi_1=0$ | py:122 / wl:147; out:21 | notes:211 | MATCH |
| Rigid-mouth $q_{\rm nt}=\Xi_1$ | py:135 / wl:157 | notes:246 | MATCH |
| Robust ceiling $0.367930328492646$ | py:146 / wl:164 | notes:284 | MATCH |
| Nonempty ceiling $0.737619063660757$ | py:147 / wl:165 | notes:305 | MATCH |
| Ceiling gap $0.369688735168111$ | py:152 (assert) | not in docs (internal consistency literal) | INTERNAL |
| Balanced family $R_1=-\xi/(1+c_\eta^2)$, $E_1=-c_\eta\xi/(1+c_\eta^2)$ | py:202-203 / wl:194-195; out py:33 | notes:366-368 | MATCH |
| Pure-target $(R_1,E_1)=(-\xi,0)$; pure-dressing $(0,-\xi/c_\eta)$ | py:176,180; out:31-32 | notes:335-337,345-347 | MATCH |

INTERNAL (scaffolding, no finding): `assert_zero` residuals, the ceiling-gap consistency literal `0.369688735168111` (a derived cross-check of the two ceilings, not a stage deliverable), Lagrange multiplier `lam`, `Piecewise`/`ConditionalExpression` strippers, pass/fail flags.

reconciliation: complete; 10 deliverable values checked, 0 misaligned

## Self-test notes

Checked the variable-independence trap at M4: r1 is a real dependency of both `qnt1` (coeff `B_star`) and `qtr1` (coeff `-C_star`), so `diff(Xi1,r1)=0` reflects a genuine cancellation, not a derivative w.r.t. an absent symbol (output line 21 confirms). No unbounded-domain integrals exist in this stage. Trivial-case spot check: with `r1=0` the rigid-mouth reductions collapse correctly and `Theta1_rigid=0` as asserted. All `positive`/`nonzero` symbol assumptions are physically warranted by the chart (radii, dressings > 0; `vareps != 0` for the strip division). No directive written (zero findings).
