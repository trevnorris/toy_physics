---
unit_id: 241
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: ["moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 241 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_241.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows 94, 1237-1281)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.txt`

## What the paper claims

Stage 241 ranks the primitive carrier weights along the selected one-parameter twin-support curve `eps_* = 1 - 3 varrho/2`, `sigma = 4/(3 varrho) - 2`, `0 < varrho < 2/3`, under the constructive coherent bound `0 < beta < 2/11`. `\stagefield{Output}`: "Primitive ranking theorem: the selected branch has a complete phase diagram with thresholds $\varrho_{W\Lambda}=2(1+\beta^2)/[3(2+\beta^2)]$ and $\varrho_{U\Lambda}=2(1+\beta^2)/[3(1+\beta+\beta^2)]$." The derivation ledger enumerates: compute the coherent weights, prove `w_chi = 2 w_Z`, derive the two thresholds, order the three ranking regions. The notes add the branch-independent identities (`w_Z > w_Lambda, w_W, |w_U|`), the factorized sign-transfer laws, the threshold ordering `0 < varrho_WLambda < varrho_ULambda < 2/3`, the derivative signs, and the numerical windows `1/3 < varrho_WLambda < 125/369` and `250/441 < varrho_ULambda < 2/3`.

## What the script claims to verify

Both scripts verify, step for step, the same eight deliverable groups (M1-M8 / sections 1-8): the support-law reduction to `eps_*` and `sigma`; strict twin-window inclusion; the weight identities incl. `w_chi = 2 w_Z` and the three exact positive-form differences; the two crossover thresholds (the `.wl` SOLVES `wW=wLambda` / `wUmag=wLambda` for the nonzero `eps_*` root then re-derives `varrho` via the support law); the factorized sign-transfer laws over the shared denominator `D`; the threshold ordering and `2/3 - varrho_ULambda`; the four endpoint values (`1/3, 125/369, 2/3, 250/441`) and the two `d/dbeta` exact forms; and representative numerical region checks for all three regimes.

## Paper <-> script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Selected curve `eps_*`, `sigma` | A1/A2 (M1) | match |
| Twin-window inclusion | A3/A4 + margins (M2) | match |
| `w_chi = 2 w_Z` and weight orderings | A5-A8 (M3) | match |
| Threshold `varrho_WLambda`, `varrho_ULambda` | A9/A10 + M4 Solve | match |
| Sign-transfer factorizations | A11/A12 (M5) | match |
| Ordering `0<varrho_WLambda<varrho_ULambda<2/3` | A13/A14 (M6) | match |
| Numerical windows `1/3, 125/369, 250/441, 2/3` | M7 (4 endpoint checks) | match |
| Three ranking regions | M8 region samples | match |
| Card `Verification`: "Mathematica audit: none yet" | a `.wl` exists and passes | mismatch (F1) |

`paper_alignment: partial` — every mathematical deliverable matches; the lone defect is the stage card's stale "Mathematica audit: none yet" line, which a `.wl` now contradicts.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 61-64 | `assert_zero(support law)` | selected curve | yes |
| A2 | sympy | 65 | `assert_zero(sigma_sel - ...)` | selected curve | yes |
| A3/A4 | sympy | 73-80 | window residuals | twin inclusion | yes |
| A5-A8 | sympy | 93-109 | weight identities | weight orderings | yes |
| A9/A10 | sympy | 124-135 | threshold identities | thresholds | yes |
| A11/A12 | sympy | 142-155 | sign-law factorizations | sign transfer | yes |
| A13/A14 | sympy | 160-175 | ordering forms | threshold ordering | yes |
| A15-A18 | sympy | 182-197 | endpoint values | numerical windows | yes |
| A19/A20 | sympy | 199-208 | `sp.diff(.,beta)` exact | derivative signs | yes |
| A21+ | sympy | 222-240 | region positivity | three regions | yes |
| M1-M8 | wl | 67-218 | mirror via independent Solve | all of the above | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (card understates verification)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_241.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.wl` (whole file)

**What's wrong:**
The card's `\stagefield{Verification}` says: "SymPy audit: ... Mathematica audit: none yet." A Mathematica audit `.wl` exists, is genuinely independent (Solve-based threshold derivation, not a port), and its committed output passes all 32 checks. The card is stale on this point.

**Why this matters:**
Reader/downstream sees single-engine when the stage is in fact dual-engine verified; understated coverage.

**Required change (user-routed):**
Update the card line 11 to cite the `.wl` and its output (paper-side edit; not Codex/script). See `## Resolve before fix_loop` in directive.

**Verification:**
Card line 11 names the `.wl` audit instead of "none yet".

## Independent-derivation check (Mathematica)

Genuinely INDEPENDENT, not a transliteration. The `.py` HARDCODES the crossover `eps_*` values as targets and asserts the threshold identity directly (`assert_zero(eps_sel - 1/(2+beta^2) - 3/2*(varrho_WLambda - varrho))`, py:124-127). The `.wl` instead SOLVES the physical crossover equations and extracts the nonzero root, then re-derives `varrho`:
- `epsWRoots = epsStar /. Solve[Together[wW - wLambda] == 0, epsStar]` -> `takeNonzeroRoot` -> checks it equals `1/(2+beta^2)` (wl:110-118);
- `varrhoWLambda = varrho /. First[Solve[epsSelected == epsWCross, varrho]]` then checks `== 2(1+beta^2)/(3(2+beta^2))` (wl:130-141).
This is a different route (Solve from the weight definitions) reaching the same closed form, exactly the dual-engine standard. The threshold endpoint `125/369` is computed by M7 from the symbolic `varrhoWLambda`, not pinned.

## Engine cross-check

Both outputs report residual = 0 / positive values for every corresponding check; 34 sympy `[ok]` lines and the full M1-M8 `.wl` pass. The `.wl` additionally prints the explicit positive sample values (e.g. Region I `wZ>wLambda = 5050/297053`) confirming strict positivity, consistent with sympy's `> 0` samples. No disagreement.

## Verdict justification

Attacks tried and failed: (1) variable-independence trap on the two `d/dbeta` checks — both `varrho_WLambda` and `varrho_ULambda` genuinely depend on `beta`, so the derivatives are non-vacuous; (2) tautology on thresholds — the `.wl` derives them via Solve, not assertion-against-self, so the cross-check is real; (3) hidden hardcode of `125/369` — it is computed from the symbolic threshold at `beta=2/11`, and both `125/369 = 2(1+4/121)/(3(2+4/121)) = 2·(125/121)/(3·250/121) = 250/750 = ... = 125/369` and `250/441` reconcile by hand; (4) region samples — they use distinct concrete `varrho` in each window and all signs are strict. The only finding is the stale card line claiming no Mathematica audit. Verdict: `findings` (one low-severity paper_misalignment, user-routed).

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `eps_* = 1 - 3 varrho/2` | py:58, wl:67-70, out wl:2-3 | tex card In/notes 84-87/appx 1240 | MATCH |
| `sigma = 4/(3 varrho) - 2` | py:65, wl:71-74 | card 9/notes 196-198/appx 1242 | MATCH |
| `varrho_WLambda = 2(1+b^2)/(3(2+b^2))` | py:121, wl:130-141 | card 15/notes 416/appx 1251 | MATCH |
| `varrho_ULambda = 2(1+b^2)/(3(1+b+b^2))` | py:122, wl:134-145 | card 15/notes 481/appx 1255 | MATCH |
| ordering `0<varrho_WLambda<varrho_ULambda<2/3` | py:160-175, wl:169-179 | notes 548/appx 1260 | MATCH |
| window low `1/3` | py:183, wl:182 | notes 577 | MATCH |
| window `125/369` (~0.338753) | py:187, wl:183 | notes 577 (`\frac{125}{369}\approx 0.338753`) | MATCH (193/369->125/369 fix HOLDS) |
| window `250/441` (~0.566893) | py:195, wl:185 | notes 595 | MATCH |
| window high `2/3` | py:191, wl:184 | notes 595 | MATCH |
| `w_chi = 2 w_Z` | py:93, wl:95 | notes 338/card ledger | MATCH |
| three ranking regions | py:222-240, wl:202-218 | notes 619-665/appx 1265-1280 | MATCH |
| Mathematica audit existence | `.wl` present, passes | card 11 "none yet" | MISMATCH (F1) |

INTERNAL (no finding): denominator `N`/`den`, `D`/`dPoly`, factorized positive-form numerators, sample `beta=1/10`, per-region `varrho` midpoints, region sample fraction values (e.g. `100701/594106`), `eps_*` crossover roots `1/(2+b^2)`, `b/(1+b+b^2)`.

The four ϱ windows `1/3, 125/369, 2/3, 250/441` fully reconcile across `.py`, `.wl` outputs, notes, and appendix; the 193/369 -> 125/369 notes fix from pass-1 HOLDS (line 577). `reconciliation: complete; 12 deliverables checked, 1 misaligned (card "Mathematica audit: none yet")`.

## Self-test notes

Checked the variable-independence trap on both `sp.diff/D[..., beta]` checks: `varrho_WLambda` and `varrho_ULambda` both contain `beta`, so derivatives are non-vacuous (no absent-variable zero). Region positivity uses three distinct in-window `varrho` values, not a single limiting case. No paper-round-trip issue: the lone finding is paper-side (stale card line), routed to user; no Codex script edit prescribed.
