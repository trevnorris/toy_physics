---
unit_id: 073
batch: III.4
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage073_family1_geometry_map.md]
  paper_appendix: present
---

# Audit unit 073 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_073.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage073_family1_geometry_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 124, 339 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.txt` — **PRESENT and committed** (the audit prompt's reading list labeled this "(missing)", but the file exists, is git-tracked since commit `9246e90`, and is fresh; see Value Reconciliation note below)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.txt` — present and fresh

## What the paper claims

Stage 073 is a **reference-branch numerical freeze**, not a new theorem. It fixes the geometry ratios of the first explicit Family-1 throat branch. `\stagefield{Output}` reads: *"The Family--1 geometry values \eqref{eq:app-stage073-Lambda-eta}."* The two boxed equations are the deliverables: (i) `L/a = 37/20`, `ell/a = 1/20`; (ii) `Lambda_ell = L/ell = 37`, `eta = 37`. The notes add the provenance: the wall width comes from `epsilon_r = 0.05 = 1/20` with the thin-layer identification `ell = epsilon_r a` (so `ell/a = 1/20`); the carried preferred aspect ratio `L/a ≈ 1.85` is frozen exactly as `Lambda_* := L/a = 37/20` (and `37/20 = 1.85` exactly); combining gives `Lambda_ell = (L/a)/(ell/a) = 37`; and under the Stage-071 Robin closure `K_m = T_X/ell`, the Robin variable `eta := K_m L / T_X = L/ell = Lambda_ell = 37`. The appendix row 124 restates the same four values verbatim. (The appendix line 339 `chi_s=37/2`, `kappa=12321/5` belong to Stage 074, per card line 31 and appendix row 126 — out of 073's scope.)

## What the script claims to verify

Both scripts assert, in order: (1) the symbolic combination law `Lambda_ell = (L/a)/(ell/a)` reduces to `L/ell` (the `a` cancels); (2) with the pinned rationals `L/a = 37/20` and `ell/a = epsilon_r = 1/20`, the freeze gives `Lambda_ell = 37`; (3) under the Robin closure substitution `K_m → T_X/ell`, `eta = K_m L / T_X` collapses to `L/ell` (the `T_X` cancels); and (4) substituting `L/ell → 37` confirms `eta(reference) = 37`. The SymPy docstring and the Mathematica banner both state the intent as "the explicit Family-1 geometry map." These are exactly the paper's four deliverables.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `ell/a = 1/20` (from `epsilon_r = 1/20`, `ell = epsilon_r a`) | py L44/L47 `epsilon_r = Rational(1,20)`, `ell_over_a = epsilon_r`; wl L37/L39 | match |
| `L/a = 37/20` (carried freeze) | py L45 `Lambda_star = Rational(37,20)`; wl L38 | match |
| `Lambda_ell = L/ell = 37` | py L42 (symbolic `=L/ell`) + L48/L55 (`Lambda_ell - 37 == 0`); wl L34/L40/L46 | match |
| `eta = K_m L/T_X = L/ell = 37` under `K_m = T_X/ell` | py L61–L65 (`eta - L/ell == 0`, `eta(ref) - 37 == 0`); wl L53–L57 | match |

Every paper-side deliverable has a faithful script-side check, and the script tests nothing the paper does not state. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero(Lambda_ell_sym - L/ell)` → `0` | combination law `Lambda_ell=L/ell` (box ii) | partial (light algebraic cancellation of `a`, but it IS the stage's content) |
| A2 | sympy | 55 | `expect_zero(Lambda_ell - 37)` from `(37/20)/(1/20)` | `Lambda_ell = 37` (box ii) | yes (substantive numeric-freeze check) |
| A3 | sympy | 64 | `expect_zero(eta - L/ell)` after `K_m→T_X/ell` | `eta = L/ell` (box ii) | partial (exercises the closure substitution / `T_X` cancellation) |
| A4 | sympy | 65 | `expect_zero(eta.subs(L/ell→37) - 37)` | `eta(ref) = 37` (box ii) | no (trivial substitution echo; redundant with A2) |
| A5 | math | 34 | `expectZero[lambdaEllSym - lSym/ellSym]` | combination law | partial (mirror of A1) |
| A6 | math | 46 | `expectZero[lambdaEll - 37]` | `Lambda_ell = 37` | yes (mirror of A2) |
| A7 | math | 56 | `expectZero[eta - len/ell]` | `eta = L/ell` | partial (mirror of A3) |
| A8 | math | 57 | `expectZero[(eta /. len/ell→37) - 37]` | `eta(ref) = 37` | no (trivial echo, mirror of A4) |

A4/A8 are trivial substitution echoes (substituting `L/ell→37` into the expression `L/ell` and subtracting 37 is zero by construction), but they are **not load-bearing**: the substantive freeze is carried by A2/A6, which derive 37 from the two independently-pinned rationals `37/20` and `1/20` — those rationals are the paper's stated inputs, not the answer pre-baked. A1/A3 are simple but genuine algebraic-cancellation identities (the `a` cancels; the `T_X` cancels) that constitute the actual content of notes §3 and §4. For a stage whose entire job is a numeric freeze of carried ratios, this assertion set is adequate, and no assertion is false or misleadingly-passing. No finding rises to a reportable category.

## Findings

None. (No `tautological_check`: A4/A8 are trivial but non-load-bearing redundant echoes, while the load-bearing checks A2/A6 derive 37 from independently-pinned inputs. No `hardcoded_result`: the `37/20` and `1/20` literals are the paper's stated inputs, and 37 is derived from them, not asserted against itself. No `symbol_assumption_error`: all of `L, a, ell, T_X, K_m` are physical lengths/scales correctly declared `positive`. No `missing_branch`/`insufficient_verification`: the stage is single-branch by construction — it *is* the chosen reference branch — and all four deliverables are exercised. No `paper_misalignment`: card box, notes, and appendix row all carry the identical four values the scripts verify. No `stale_output`: both outputs are newer than their scripts. No `missing_verification_script`: both engines present.)

## Independent-derivation check (Mathematica)

The `.wl` is a close structural mirror of the `.py` (same four checks, same `epsilonR=1/20`, `lambdaStar=37/20`, same Robin substitution `km -> tx/ell`). However, this stage's "derivation" is nothing more than arithmetic on two pinned rationals plus a one-symbol cancellation — there is essentially only one way to express `(37/20)/(1/20) = 37` and `K_m L/T_X` under `K_m = T_X/ell`. There is no richer independent route an honest second engine could take here, so the structural parallelism does **not** rise to a `mathematica_transliteration` finding: the policy targets cases where the second engine echoes a *non-trivial* algebraic choreography it should have re-derived. Here the "algebra" is a single fraction and a single cancellation; both engines independently arrive at the same forced result. I tried to find a place where the `.wl` could plausibly have made an independent choice (a different elimination order, a different normalization) and there is none available. Not flagged.

## Engine cross-check

Both engines agree exactly. Side by side (committed outputs):

- SymPy `..._sympy_audit_refresh.txt`: `Lambda_ell - L/ell (symbolic) = 0`; `Lambda_ell = L/ell = 37`; `Lambda_ell - 37 = 0`; `eta under K_m = T_X/ell -> L/ell`; `eta - L/ell = 0`; `eta(reference) - 37 = 0`; ledger `Lambda_ell = 37`, `eta = 37`.
- Mathematica `..._mathematica_audit.txt`: identical residuals (`Lambda_ell - 37 = 0`, `eta - L/ell = 0`, `eta(reference) - 37 = 0`), all four `PASS:` lines, `L/a = 37/20`, `ell/a = 1/20`, "Stage 073 Mathematica audit passed."

`engines_agree: true`.

## Verdict justification

Clean. The paper card, notes, and appendix row all state the identical four-value freeze (`L/a=37/20`, `ell/a=1/20`, `Lambda_ell=eta=37`), and both engines verify exactly those values via checks that, while algebraically light (this stage is by design a numeric freeze of two carried ratios), are correctly anchored: 37 is *derived* from the independently-pinned `37/20` and `1/20`, not asserted against itself, and the `eta` check genuinely exercises the Robin closure substitution `K_m=T_X/ell` (the `T_X` cancellation). Attacks tried that failed: (1) hunted for a hardcoded `37` that the check confirms against itself — found instead that 37 is produced from the two pinned rationals; (2) checked whether `37/20` equals the notes' `≈1.85` — it does, exactly; (3) checked whether the `chi_s=37/2`/`kappa=12321/5` values in appendix line 339 belong to this stage — they are Stage-074 outputs and out of scope; (4) looked for a symbol-domain error in the `positive` declarations — none, all are physical positive lengths/scales; (5) probed the `.wl` for transliteration — the algebra is too trivial to admit an independent route, so parallelism is not a defect here. I confirm I read the paper card, the notes, and the appendix rows, and the script's verified claim matches the paper's claim exactly.

## Value Reconciliation (pass-2 augmentation)

Both committed outputs were read directly (the SymPy `..._sympy_audit_refresh.txt` is present, git-tracked, and fresh — newer than its `.py`; the audit prompt's "(missing)" label for it is incorrect, so the reconciliation rests on the actual committed output, not source-only reasoning). No `stale_output` signal applies (output mtimes 02:44:04 / 02:44:34 both postdate script mtimes 02:39:28 / 02:39:30).

Deliverable-level reconciliation:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `epsilon_r = 1/20` | py L44, wl L37; out "epsilon_r = 1/20" | notes §1 ("`epsilon_r = 0.05 = 1/20`") | MATCH |
| `ell/a = 1/20` | py L47, wl L39; out "ell/a = 1/20" | tex L19 `\frac{\ell}{a}=\frac{1}{20}`; notes §1; appendix L124 `\ell/a=1/20` | MATCH |
| `L/a = 37/20` | py L45, wl L38; out "L/a = 37/20" | tex L18 `\frac{L}{a}=\frac{37}{20}`; notes §2 (`Lambda_* := L/a = 37/20`); appendix L124 `L/a=37/20` | MATCH |
| `Lambda_ell = 37` | py L48/L55, wl L40/L46; out "Lambda_ell = L/ell = 37" | tex L26 `\Lambda_\ell=\frac{L}{\ell}=37`; notes §3; appendix L124 `\Lambda_\ell=\eta=37` | MATCH |
| `eta = 37` | py L65, wl L57; out "eta = 37" (ledger) | tex L27 `\eta=37`; notes §4; appendix L124 `\Lambda_\ell=\eta=37` | MATCH |

Internal scaffolding (accounted for, no finding): the symbolic residuals `Lambda_ell - L/ell (symbolic) = 0`, `Lambda_ell - 37 = 0`, `eta - L/ell = 0`, `eta(reference) - 37 = 0`, and the intermediate symbolic forms `eta -> L/ell` / `len/ell` — these are verification residuals and intermediate expressions that exist only to drive the `expect_zero`/`expectZero` assertions, not labeled deliverable values, so they are correctly absent from the prose.

reconciliation: complete; 5 values checked, 0 misaligned.

## Self-test notes

Traps checked: (1) Variable independence — no `sp.diff`/`D[]` anywhere, so the zero-derivative trap is N/A. (2) Symmetry/parity — no integrals, N/A. (3) Trivial-case pre-check — confirmed `(37/20)/(1/20)=37` and `(T_X/ell)·L/T_X = L/ell` (T_X cancels) → `L/ell→37` gives 37, so every `expect_zero` genuinely reduces to 0. I also confirmed the only "could-this-fail" risks (A4/A8 trivial echoes) are non-load-bearing and shadowed by the substantive freeze checks A2/A6. No directive written (zero findings).
