---
unit_id: 099
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage099_reduced_finish_line.md]
  paper_appendix: present
---

# Audit unit 099 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_099.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage099_reduced_finish_line.md` (only file matching the glob)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage099_reduced_finish_line_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage099_reduced_finish_line_mathematica_audit.txt`

## What the paper claims

Stage 099 is the "reduced finish line" geometry-lane firewall ledger step. The card states the only remaining reduced theorem gate is `N_Q = 1` (passive/outgoing quadrupole normalization), with everything else reduced away. The Derivation ledger isolates the forced conservative carrier `Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)` (Inputs carry the branch identity `K_0 K_4 = 4 K_2^2`). The card's Checks list three carry-in checks: (a) static limit `eps_2=eps_4=0` returns `c_pole=1/4`; (b) `l=0/l=2` orthogonality before the firewall; (c) any support/source success carries the minimal-module hypothesis. The notes add the canonical-invariant form `Kbar_Q^cons = Kbar_0[3/4 + (1/4)/(1 - omega^2/Omega_Q^2)]`, the target normalization `Kbar_0^target = 64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` (with `Omega_Q = 3 c_s/(2a)`), the definition `N_Q := Kbar_0/Kbar_0^target`, and the equivalent reductions `rho_alpha=4/3`, `zeta_req=1/3`. Appendix Part IV (Sec. retarded-2.5PN collapse) supplies the load-bearing equations eq:app-part04-kbar-cons, -kbar-even-moments, -kbar0-target, -NQ-def, -gamma5-chiN, and -factorized-defect-again (`Gammabar_5/(2G/(5c^5)) = chi_Q N_Q`).

## What the script claims to verify

The docstring states the stage exercises the static-slot value and pole residue of `Yhat_Q^cons` locally as a sanity anchor on the partial-fraction form, plus the structural even-moment relations forced by that module, and explicitly disclaims (i)/(ii)/(iii) as upstream carry-ins from stages 091/092/094/096 and 088/089/090. The six assertions verify: (1) `Yhat_Q^cons(0) = 1`; (2) the pole residue of `Yhat_Q^cons` at `omega=Omega_Q` is `-Omega_Q/8` (the `1/4` partial-fraction residue, card check (a)); (3) the geometric form of `K0_target`; (4) the minimal branch identity `K0 K4 = 4 K2^2` on the structural even moments; (5) the `Gamma_5` sqrt form equals the canonical odd-coeff form; (6) the normalization `Gamma5_struct/(2G/(5c^5)) = N_Q` on the `chi_Q=1` branch (with `K0 = N_Q*K0_target`). `K0_sym` is a free positive symbol, so the structural assertions are not forced by an `N_Q` recipe.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` (eq:app-part04-kbar-cons; card ledger) | `Yhat_Q_cons` def + static-slot=1 + residue=-Omega_Q/8 | match |
| card check (a): static eps=0 ⇒ `c_pole=1/4` (eq:app-part04-cpole-dynamic-geometry at eps=0) | pole residue `-Omega_Q/8` anchors the `1/4` partial-fraction weight | match (residue route to `c_pole=1/4`) |
| branch identity `K_0 K_4 = 4 K_2^2` (Inputs; eq:app-part04-minimal-branch-identity) | `K0_sym*K4_struct - 4*K2_struct^2 == 0` | match |
| even moments `Kbar_2=Kbar_0/(4Omega_Q^2)`, `Kbar_4=Kbar_0/(4Omega_Q^4)` (eq:app-part04-kbar-even-moments) | `K2_struct`, `K4_struct` defs (printed) | match |
| `Kbar_0^target = 64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` (eq:app-part04-kbar0-target; notes §2) | `K0_target geometric form == 0` | match |
| `Gammabar_5 = chi_Q 9 Kbar_0/(32 Omega_Q^5)` (eq:app-part04-gamma5-chiN) | `Gamma5_struct` def + sqrt-form identity | match |
| `Gammabar_5/(2G/(5c^5)) = chi_Q N_Q` (eq:app-part04-factorized-defect-again) | check (6) on chi_Q=1 branch ⇒ `= N_Q` | match (chi_Q=1 slice; chi_Q general is stages 100-106) |
| card check (b): l=0/l=2 orthogonality | (none — disclaimed carry-in to stage 094) | n/a (status carve-in, explicitly disclaimed) |
| card check (c): minimal-module hypothesis | (none — disclaimed carry-in to 088/089/090) | n/a (status carve-in, explicitly disclaimed) |
| `N_Q = 1` gate, `rho_alpha=4/3`, `zeta_req=1/3` | not computed (carry-forward statements, not stage-local) | n/a (notes-level carry-forward, not script deliverable) |

`paper_alignment: aligned`. Every algebraic identity the scripts emit corresponds verbatim to a card/notes/appendix equation; the two unmodelled card checks are explicitly and correctly flagged in the docstring as upstream carry-ins (the card itself labels (b)/(c) as carry-forward, not stage-099-local derivations), consistent with the card's own "geometry-lane firewall ledger step" framing.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50-51 | `Yhat_Q_cons.subs(omega,0) - 1 == 0` | eq:app-part04-kbar-cons static slot | yes |
| A2 | sympy | 56-57 | `residue(Yhat,omega,Omega_Q) - (-Omega_Q/8) == 0` | card check (a) / `c_pole=1/4` weight | yes |
| A3 | sympy | 59-61 | `K0_target_geom - 54 G c_s^5/(5 a^5 c^5) == 0` | eq:app-part04-kbar0-target | yes |
| A4 | sympy | 71 | `K0_sym*K4_struct - 4*K2_struct^2 == 0` | branch identity (Inputs) | yes |
| A5 | sympy | 74-75 | `9*K2^(5/2)/K0^(3/2) - Gamma5_struct == 0` | eq:app-part04-gamma5-chiN | yes |
| A6 | sympy | 79-80 | `Gamma5_struct[K0->N_Q K0_target]/(2G/(5c^5)) - N_Q == 0` | eq:app-part04-factorized-defect-again | yes |
| B1 | mathematica | 37 | `(yhatCons /. omega->0) - 1 == 0` | = A1 | yes |
| B2 | mathematica | 38-39 | `Residue[yhatCons,{omega,omegaQ}] - (-omegaQ/8) == 0` | = A2 | yes |
| B3 | mathematica | 41-43 | `k0TargetGeom - 54 G cS^5/(5 a^5 c^5) == 0` | = A3 | yes |
| B4 | mathematica | 52 | `k0Sym*k4Struct - 4*k2Struct^2 == 0` | = A4 | yes |
| B5 | mathematica | 53-54 | `9*Sqrt[k2Struct^5/k0Sym^3] - gamma5Struct == 0` | = A5 | yes |
| B6 | mathematica | 55-56 | `(gamma5Struct /. k0Sym->nQ*k0Target)/(2G/(5c^5)) - nQ == 0` | = A6 | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is structurally parallel to the `.py`: identical six checks in identical order, identical variable choreography (`yhatCons`↔`Yhat_Q_cons`, `k0Target`↔`K0_target`, `k2Struct`/`k4Struct`/`gamma5Struct`↔`K2_struct`/`K4_struct`/`Gamma5_struct`). I considered `mathematica_transliteration`. The stage's deliverables, however, are all closed-form algebraic identities lifted directly from the appendix (eq:app-part04-kbar-cons … -factorized-defect-again); there is no alternative "physical re-derivation route" to demand of a second engine — the only legitimate second-engine job is an independent CAS simplification of the same residuals, which is exactly what the `.wl` performs (its own `Residue`, `FullSimplify`, `/.` substitutions, distinct from SymPy's `residue`/`simplify`). The parallel structure is therefore inherent to an identity-verification stage, not an illegitimate echo of an algorithmic derivation, so it does not meet the bar for a transliteration finding. Note: residue/static-slot (A2/B2) and the sqrt-form identity (A5/B5) genuinely require each engine's own simplifier to evaluate, so the engines are not merely re-running shared scaffolding.

## Engine cross-check

Both outputs report `= 0` (sympy) / `= 0` + `PASS:` (mathematica) for all six checks and end with `STAGE 099 AUDIT PASSED`. The two engines render `Yhat_Q^cons` in different but equivalent forms — sympy: `(Omega_Q**2 - 3*omega**2/4)/(Omega_Q**2 - omega**2)`; mathematica: `(3 + (1 - omega^2/omegaQ^2)^(-1))/4` — confirming independent simplification rather than a copied normal form. All six residuals agree (both zero). `engines_agree: true`.

## Verdict justification

`clean`. I read the card, notes, and Part IV appendix before the scripts and built the claim model, then attacked each assertion. Hand-verification confirmed all six identities: static slot (`3/4 + 1/4 = 1`), pole residue (`(1/4)·Omega_Q²·(-1)/(2Omega_Q) = -Omega_Q/8`), the geometric `K0_target` (`64·243/(45·32) = 54/5`), the branch identity (`K0²/(4Omega_Q⁴) = 4·K0²/(16Omega_Q⁴)`), the sqrt form (`4^(5/2)=32` makes `9·K0^(5/2)/(32Omega_Q⁵)/K0^(3/2) = 9K0/(32Omega_Q⁵)`), and the `N_Q` normalization (`(2/5)N_Q G/c⁵ ÷ (2G/5c⁵) = N_Q`, requiring the `9/32`, `64/45`, and `2/5` constants to align — non-tautological). `K0_sym` is a free symbol, so A4–A6 cannot be passed by an `N_Q` recipe; A6 in particular would fail under any wrong constant. Every emitted value matches the appendix verbatim (see reconciliation). The two card checks the script doesn't model (l-orthogonality, minimal-module hypothesis) are explicitly disclaimed in the docstring as upstream carry-ins and are labelled carry-forward by the card itself, so this is not `script_missing_paper_claim`. Outputs are fresh (both `.txt` mtimes 2026-05-27 14:28/14:30 postdate both script mtimes 11:16). Both engines present, substantive, and agreeing. No findings.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 7 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` | py:46 / wl:33; sympy out L5, math out L5 | notes:20 `Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`; tex:13; appendix eq:app-part04-kbar-cons (L219-227) | MATCH |
| static slot `Yhat_Q^cons(0) = 1` | py:50-51 / wl:37; out L6 | implied by notes:20 / appendix L219-227 (`3/4+1/4=1`) | MATCH (structural consequence) |
| pole residue `-Omega_Q/8` ⇔ `c_pole=1/4` | py:56-57 / wl:38-39; out L7 | tex:22 (`c_pole=1/4`); appendix eq:app-part04-cpole-dynamic-geometry (L138-153), -geometry-static-split `K_pole=K0/4` (L119-123) | MATCH |
| `K0_target = 64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` | py:59-61 / wl:41-43; out L8 | notes:52-57; appendix eq:app-part04-kbar0-target (L236-245) | MATCH |
| `K2 = K0/(4 Omega_Q^2)`, `K4 = K0/(4 Omega_Q^4)` (and branch identity `K0 K4 = 4 K2^2`) | py:66-67,71 / wl:49,52; out L9,L12-13 | tex:9 (branch identity); appendix eq:app-part04-kbar-even-moments (L229-234), -minimal-branch-identity (L113-114) | MATCH |
| `Gamma5 = 9 K0/(32 Omega_Q^5)` (sqrt-form identity) | py:68,74-75,84 / wl:50,53-54,60; out L10,L14 | appendix eq:app-part04-gamma5-chiN (L278-281, chi_Q=1) | MATCH |
| `Gammabar_5/(2G/(5c^5)) = N_Q` (chi_Q=1) | py:79-80 / wl:55-56; out L11 | appendix eq:app-part04-factorized-defect-again (L284-288, chi_Q=1) | MATCH |

INTERNAL (scaffolding, no finding): `G,c,c_s,a,Omega_Q,omega,N_Q,K0_sym` symbol decls; per-check residual prints (all `= 0`); `PASS:` flags; `STAGE 099 AUDIT PASSED` sentinel. Note: notes-level carry-forward statements `rho_alpha=4/3`, `zeta_req=1/3`, and the `N_Q=1` gate are NOT emitted by either script (they are carry-ins / the open-problem statement), so they are not in scope for the script→doc reconciliation; they appear in notes:25-26,61 and tex:27.

## Self-test notes

Checked variable-independence (no `diff`/`D` in this script — all checks are algebraic residuals on closed forms, so the zero-derivative trap does not apply). Checked the residue/parity claim by hand: `Residue` of the simple pole `(1/4)Omega_Q²/((Omega_Q-omega)(Omega_Q+omega))` at `omega=Omega_Q` is `-Omega_Q/8` — correct and the assertion can genuinely fail for a wrong residue. Trivial-case/constant check: A3, A5, A6 each require a specific rational combination of literals to vanish (`64·243/(45·32)=54/5`; `4^(5/2)=32`; `9·64/(32·45)=2/5` then `÷ 2/5 = N_Q`), so none collapse to `0==0` and `K0_sym`/`N_Q` being free symbols means they are non-tautological. Paper round-trip: every script literal traces to an appendix/notes equation with no new constant introduced.
