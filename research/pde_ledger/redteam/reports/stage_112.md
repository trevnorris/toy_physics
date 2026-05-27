---
unit_id: 112
batch: IV.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage112_hybrid_robin_mixed_compensation.md"]
  paper_appendix: present
---

# Audit unit 112 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_112.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage112_hybrid_robin_mixed_compensation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_112}` at line 1258; no separate prose row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.txt`

## What the paper claims

Stage 112 (rendered "Stage 129" in the section heading; LaTeX label `stage:112`) is the "Exact Robin–Mixed Compensation Law" in the outlet-deformation block. The card quotes verbatim: *"Nontrivial compensated branch has \(\rho_R=4\sigma_W\), \(\kappa_W=1/3\), and preserves odd normalization iff \(\gamma_W=1/9\)."* The accompanying notes give the hybrid outlet model `Lambda_2^hyb = Lambda_2^out + rho_R - sigma_W/(1 - kappa_W z^2 - i gamma_W z^5) + O(z^6)`, the canonical-even system (`-L2/L0 = 1/9`, `L2^2/L0^2 - L4/L0 = 4/81`), the resulting two branches `(rho_R=sigma_W, kappa_W=0)` (trivial) and `(rho_R=4 sigma_W, kappa_W=1/3)` (nontrivial), the explicit `chi_Q^hyb = (1 - 9 sigma_W gamma_W)/(1 - sigma_W)` on the nontrivial branch, the preservation condition `gamma_W = 1/9`, and the collapse identity `Lambda_2^hyb = (1 - sigma_W) Lambda_2^out + O(z^6)` at gamma_W=1/9. The notes also offer a Stage-92 cross-check that the same gamma_W=1/9 follows from `(b, a_0, a_5) = (0, 3 sigma_W, -sigma_W gamma_W)` and the linear constraint `a_0/3 + 9 a_5 = 0`; this is presented as a redundant secondary derivation of the same gamma_W=1/9.

## What the script claims to verify

Both scripts compute the same series expansion of `Lambda_hyb` through z^5, extract L0, L2, L4, L5 (the L5 coefficient is normalized by `/I` to recover the real-valued odd coefficient), pose the two canonical-even constraints exactly as written in the notes, and verify that the solve returns exactly the two branches `{kappa=0, rho=sigma}` and `{kappa=1/3, rho=4 sigma}`. They then evaluate `(-L5/L0)/(1/27)` on each branch and assert `chi_A = 1 - 9 sigma gamma` and `chi_B = (1 - 9 sigma gamma)/(1 - sigma)`. They further check that on branch B with gamma=1/9, chi_Q → 1, and that the full `Lambda_hyb` collapses to `(1 - sigma) Lambda_out` on branch B with gamma=1/9. The docstring header on the SymPy file says "Stage 95 SymPy audit", the SymPy final print says `stage95: PASS`, and the Mathematica banner says "STAGE 095 — EXACT ROBIN-MIXED COMPENSATION LAW"; only the Mathematica trailing print says "Stage 112 Mathematica audit passed."

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Hybrid outlet form `Lambda_hyb = Lambda_out + rho_R - sigma_W/(1 - kappa_W z^2 - i gamma_W z^5)` | `Lambda_hyb = ... series(... rho - sigma/(1 - kappa*z^2 - I*gamma*z^5) ...)` (sympy:24; math:34) | match |
| Two-branch canonical-even solve: A=(sigma, 0), B=(4 sigma, 1/3) | `assert sols == [{kappa: 0, rho: sigma}, {kappa: 1/3, rho: 4*sigma}]` (sympy:38) and per-branch expectZero checks (math:50–53) | match |
| chi_Q^hyb = (1 - 9 sigma gamma)/(1 - sigma) on branch B | `assert sp.simplify(chi_B - (1 - 9*sigma*gamma)/(1 - sigma)) == 0` (sympy:48; math:62) | match |
| Preservation iff gamma_W = 1/9 | `assert sp.simplify(chi_B.subs(gamma, 1/9) - 1) == 0` (sympy:49; math:63) | match |
| Collapse identity `Lambda_hyb = (1-sigma) Lambda_out` at branch B, gamma=1/9 | `assert sp.simplify(scaled_identity) == 0` (sympy:53; math:67) | match |
| Stage-92 deformation data `(b, a_0, a_5) = (0, 3 sigma_W, -sigma_W gamma_W)` and `a_0/3 + 9 a_5 = 0` as redundant linear check | (none) | missing (secondary cross-check only; equivalent to verified gamma=1/9) |

Set `paper_alignment: partial` because the script-side label/docstring contradict the paper card (mislabeled as Stage 95/095) even though the math itself matches. The notes' secondary `(b, a_0, a_5)` linear cross-check is not exercised but is mathematically redundant with verified content; I do not file it as `script_missing_paper_claim` since the same gamma_W=1/9 is independently shown via chi_Q and the scaled identity.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 38 | `assert sols == [{kappa: 0, rho: sigma}, {kappa: 1/3, rho: 4*sigma}]` | two-branch canonical-even result | yes |
| A2 | sympy | 43 | `assert sp.simplify(chi_A - (1 - 9*sigma*gamma)) == 0` | chi_Q on trivial branch A | yes |
| A3 | sympy | 48 | `assert sp.simplify(chi_B - (1 - 9*sigma*gamma)/(1 - sigma)) == 0` | chi_Q^hyb formula on branch B | yes |
| A4 | sympy | 49 | `assert sp.simplify(chi_B.subs(gamma, 1/9) - 1) == 0` | preservation iff gamma=1/9 | yes |
| A5 | sympy | 53 | `assert sp.simplify(scaled_identity) == 0` | collapse to (1-sigma) Lambda_out | yes |
| A6 | mathematica | 50 | `expectZero["branch A rho - sigma", (rho /. solA) - sigma]` | branch A rho identity | yes |
| A7 | mathematica | 51 | `expectZero["branch A kappa", kappa /. solA]` | branch A kappa identity | yes |
| A8 | mathematica | 52 | `expectZero["branch B rho - 4 sigma", (rho /. solB) - 4*sigma]` | branch B rho identity | yes |
| A9 | mathematica | 53 | `expectZero["branch B kappa - 1/3", (kappa /. solB) - 1/3]` | branch B kappa identity | yes |
| A10 | mathematica | 61 | `expectZero["chi_Q branch A - (1 - 9 sigma gamma)", chiA - (1 - 9*sigma*gamma)]` | chi_Q branch A | yes |
| A11 | mathematica | 62 | `expectZero["chi_Q branch B - (1 - 9 sigma gamma)/(1 - sigma)", chiB - (1 - 9*sigma*gamma)/(1 - sigma)]` | chi_Q branch B formula | yes |
| A12 | mathematica | 63 | `expectZero["chi_Q branch B at gamma=1/9", (chiB /. gamma -> 1/9) - 1]` | preservation iff gamma=1/9 | yes |
| A13 | mathematica | 67 | `expectZero["scaled identity on branch B", scaledIdentity]` | collapse identity | yes |

All assertions trace to specific paper claims and are non-tautological (each branch substitution drops the symbolic dependence on the substituted variables and verifies the surviving form matches an a-priori-different target expression).

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py:22-53`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:33-67`

**What's wrong:**
The Mathematica script is a near-perfect line-for-line transliteration of the SymPy script's algebraic choreography rather than an independent re-derivation. Both scripts:

1. Build `Lambda_hyb` via `series(..., z, 0, 6)` / `Series[..., {z, 0, 5}]` of the *same* algebraic form (sympy:24; math:34):
   - sympy: `Lambda_hyb = sp.expand(sp.series(Lambda_out + rho - sigma/(1 - kappa*z**2 - I*gamma*z**5), z, 0, 6).removeO())`
   - math: `lambdaHyb = Expand[Normal[Series[lambdaOut + rho - sigma/(1 - kappa*z^2 - I*gamma*z^5), {z, 0, 5}]]];`
2. Extract coefficients in the same order with the same `/I` normalization for L5 (sympy:26-29; math:36-39):
   - sympy: `L5 = sp.simplify(Lambda_hyb.coeff(z, 5) / I)`
   - math: `l5 = FullSimplify[Coefficient[lambdaHyb, z, 5]/I, Assumptions -> $Assumptions];`
3. Pose the identical two-equation canonical-even system and call the engine's generic solver (sympy:31-34; math:41).
4. Compute `chi` as `(-L5/L0)/(1/27)` substituted with the same branch and assert the same final reduced forms (sympy:41-48; math:55-62).
5. Compute the "scaled identity" using the same name and the same substitution `(Lambda_hyb on branch B - (1-sigma) Lambda_out)` at `gamma=1/9` (sympy:51-53; math:65-67).

No independent algebraic path is taken. A genuinely independent Mathematica derivation would, for example, use `SeriesCoefficient`, factor `Lambda_hyb - (1-sigma) Lambda_out` directly via `Factor`/`Reduce`, or derive `gamma_W=1/9` from the linearized Stage-92 `a_0/3 + 9 a_5 = 0` condition that the notes provide as an independent cross-check (the notes explicitly give a second derivation path that the script ignores). Either approach would constitute genuine engine independence.

**Why this matters:**
The second-engine policy requires both engines to derive the result independently from the physical premises. A transliteration only catches typos and language-port bugs; it does not catch shared algebraic mistakes (e.g., an off-by-one in the series order, a wrong sign in the normalization, a miscopied formula from the notes), because both scripts will inherit the same bug.

**Required change:**
Refactor the Mathematica script so that it derives `chi_Q^hyb` and the `gamma_W=1/9` condition by an independent algebraic path. Recommended approach: re-derive `gamma_W = 1/9` via the Stage-92 linearized cross-check explicitly described in the notes ("Branch-selection data" section: `(b, a_0, a_5) = (0, 3 sigma_W, -sigma_W gamma_W)`, condition `a_0/3 + 9 a_5 = 0`), thereby using a *different* algebraic identity to land the same gamma. Cross-engine agreement then becomes a real test rather than a literal echo.

**Verification:**
After refactor, the Mathematica script must (i) compute `gamma_W` via the linearized `(b, a_0, a_5)` cross-check independently of the chi_Q closed-form path, (ii) assert `gamma_W = 1/9` via that route, (iii) still confirm the final scaled-identity collapse to `(1-sigma) Lambda_out`, and (iv) continue to exit 0 on `redteam exec-mathematica 112`.

### F2 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_112.tex:1` quote: `\section[Stage~129]{Stage~129: Exact Robin–Mixed Compensation Law}` (LaTeX label is `stage:112`).
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py:3` quote: `Stage 95 SymPy audit.`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py:54` quote: `print('stage95: PASS')`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl:26` quote: `banner["STAGE 095 — EXACT ROBIN-MIXED COMPENSATION LAW"];`

**What's wrong:**
The SymPy docstring and the SymPy final print label this unit as "Stage 95", and the Mathematica banner says "STAGE 095". The paper card heading prints "Stage 129" with LaTeX label `stage:112`. None of the three labels agree. The actual mathematics in the scripts matches the paper card; the discrepancy is purely in the human-readable identifying strings, but it is real and confusing for any operator reading the script transcripts side-by-side with the paper.

**Why this matters:**
Mislabeled audit transcripts make manual cross-referencing error-prone: a reader scanning the Mathematica output for "Stage 112" will not find it. It also raises the question of whether the script was forked from an earlier "Stage 95" file without checking that the algebra was updated to the current stage's claim. (The algebra here matches, so the answer is no — but the label confusion legitimately invites that doubt.)

**Required change:**
Update the script docstring/banner/final print strings to reference unit 112 (or "Stage 129" per the paper's display section number, whichever convention the user maintains). Specifically:
- sympy:3 docstring: change "Stage 95 SymPy audit." to "Stage 112 SymPy audit." (or "Stage 129" if the project uses display numbers).
- sympy:54: change `print('stage95: PASS')` to `print('stage112: PASS')`.
- math:26: change `banner["STAGE 095 — EXACT ROBIN-MIXED COMPENSATION LAW"];` to `banner["STAGE 112 — EXACT ROBIN-MIXED COMPENSATION LAW"];` (or "STAGE 129"). The math:70 trailing print already says "Stage 112 Mathematica audit passed."; align the banner to that.

**Verification:**
Re-running the scripts produces output transcripts whose headers reference Stage 112 consistently across both engines and across all banner/intermediate/final strings.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration — see F1 for quoted side-by-side excerpts. The same algebraic choreography appears in both files: identical series expansion, identical coefficient extraction with `/I` normalization, identical solve target, identical chi formula `(-L5/L0)/(1/27)` with identical branch substitution, identical "scaled identity" construction. The notes contain a Stage-92 cross-check (`a_0/3 + 9 a_5 = 0`) that the Mathematica script could use to derive `gamma_W=1/9` independently; this independent path is not exercised by either script.

## Engine cross-check

Both engines produce the same final outputs (engines agree):
- SymPy: `chi_Q branch A = -9*gamma*sigma + 1`, `chi_Q branch B = (9*gamma*sigma - 1)/(sigma - 1)`, `scaled identity on branch B = 0`.
- Mathematica: `chi_Q branch A = 1 - 9*gamma*sigma`, `chi_Q branch B = (-1 + 9*gamma*sigma)/(-1 + sigma)`, `scaled identity on branch B = 0`.

The forms are algebraically identical (`(9 g s - 1)/(s - 1) = (1 - 9 g s)/(1 - s)`); both `expectZero` / `assert` checks pass. No engine_disagreement.

## Verdict justification

The math itself matches the paper card and notes line-for-line: two-branch canonical-even solve, branch-B chi_Q formula, preservation at gamma=1/9, collapse to `(1 - sigma) Lambda_out`. All five paper deliverables have anchored, non-tautological assertions in both engines. Attempted attacks: (a) checked whether the L5 / I normalization could mask a sign in the imaginary coefficient — it does not, since L5 is then asserted against an explicit symbolic form `1/9 - sigma*gamma` (via chi_Q); (b) checked whether `simplify` under `sigma != 1` assumption could hide a removable singularity — sigma=1 is genuinely excluded by the geometric setup since branch B collapse formula `Lambda_hyb = (1-sigma) Lambda_out` becomes degenerate there, so the assumption is physically justified; (c) checked whether the series truncation order matches between engines — sympy `series(..., z, 0, 6)` and math `Series[..., {z, 0, 5}]` both retain through z^5 inclusive, so coefficients are captured consistently; (d) checked branch ordering — Mathematica explicitly sorts by kappa to ensure A=trivial, B=nontrivial, which aligns with the SymPy convention. The two findings are (F1) the Mathematica script transliterates rather than independently re-derives, and (F2) the script labels misidentify the unit as Stage 95/095. Neither blocks the result; both warrant correction.

## Self-test notes

I checked: (1) variable independence — `chi_B.subs(gamma, 1/9)` substitutes for the only remaining symbol that prevents reduction; `chi_B` is not identically free of `gamma` (i.e., `gamma` genuinely appears) so the substitution test is meaningful, confirmed by inspecting `chi_Q branch B = (9*gamma*sigma - 1)/(sigma - 1)`. (2) Symmetry/parity — N/A here, no symmetric-domain integrals. (3) Trivial-case pre-check — substituting `sigma=0, gamma=0` into `chi_B - (1 - 9*sigma*gamma)/(1 - sigma)` yields `1 - 1 = 0`, consistent with assertion. (4) Path specifications — N/A, no missing-script findings. (5) Paper round-trip — F1's proposed independent path (Stage-92 linearized cross-check) is explicitly named in the notes, so no new paper_misalignment is introduced by recommending it; F2's relabeling is cosmetic and does not touch the math.
