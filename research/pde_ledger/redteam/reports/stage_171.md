---
unit_id: 171
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage171_microscopic_grouped_obstructions.md]
  paper_appendix: present
---

# Audit unit 171 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_171.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage171_microscopic_grouped_obstructions.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 73, 417-459, 1463)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py`
- mathematica: `/var/projects/toy_projects/.../mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl` — actual path `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` states verbatim: "Decomposes \(\mathcal K_A\) and \(\mathcal G_A\) into wall, BdG, Maxwell/mixed, and outgoing-transfer defect bundles." The notes (authoritative on detail) enumerate the deliverables: (1) the regrouping `K_A = W_A - B_A - Z_A` with `W_A = 1/9 δK_A - δM_A`, `B_A = δB_{A,2} + 1/9 δB_{A,0}`, `Z_A = δZ_{A,2} + 1/9 δZ_{A,0}`; (2) the regrouping `G_A = -P_0 δK_A + P_0 δB_{A,0} + N_A` with `N_A = δN_{A,0} + P_0 δZ_{A,0}`; (3) exact BdG linearizations `δB_0`, `δB_2`, and the bundled `B_A` from `B_0 = c²/ϖ²`, `B_2 = c²/ϖ⁴`; (4) exact Maxwell/mixed linearizations `δZ_0`, `δZ_2`, the bundle `Z_A`, the static-transfer `δN_0`, and the bundle `N_A`, from `Z_0 = Q/Δ`, `Z_2 = (QS - GΔ)/Δ²`, `N_0 = P²/Δ²`; (5) the five primitive port variations `δΔ, δS, δG, δP, δQ` in terms of `δΩ_U², δΩ_W², δR, δg_U, δg_W`; (6) the weak-axisymmetric collapse to the scalar pair `(𝔎_1, 𝔊_1)` under the lane signature `(λ_{20}, λ_{21}, λ_{22}) = (1, 1/2, -1)`. The appendix (lines 434-459) restates (1) and (2) verbatim with the same coefficients.

## What the script claims to verify

The SymPy script runs 21 `expect_zero` checks: two regrouping identities for `K_A` and `G_A` (lines 39-45); the two BdG linearizations and their bundle `B_A` derived via `sp.diff` of `c²/w²` and `c²/w⁴` and compared to the notes' boxed formulas (lines 56-69); the five primitive port variations `δΔ, δS, δG, δP, δQ` derived via `sp.diff` of the port invariants and compared to the boxed forms (lines 87-107); the `δZ_0, δZ_2, δN_0` linearizations and the `Z_A`/`N_A` bundles derived via `sp.diff` of `Q/Δ`, `(QS-GΔ)/Δ²`, `P²/Δ²` (lines 115-148); and the weak-axisymmetric regrouping of `(𝔎_1, 𝔊_1)` over all three lanes `λ = 1, 1/2, -1` (lines 158-164). The Mathematica `.wl` mirrors all 21 checks in the same order. This is exactly the deliverable set the notes enumerate.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `K_A = W_A - B_A - Z_A` regrouping | sympy 39-41 / wl 31-33 (`K_A split`) | match |
| (2) `G_A = -P_0 δK_A + P_0 δB_0 + N_A` regrouping | sympy 43-45 / wl 35-37 (`G_A split`) | match |
| (3) `δB_0`, `δB_2`, `B_A` bundle | sympy 56-69 / wl 45-58 | match |
| (4) `δZ_0`, `δZ_2`, `δN_0`, `Z_A`, `N_A` | sympy 115-148 / wl 99-127 | match |
| (5) primitive `δΔ, δS, δG, δP, δQ` | sympy 87-107 / wl 70-89 | match |
| (6) weak-axisymmetric `(𝔎_1, 𝔊_1)` over 3 lanes | sympy 158-164 / wl 132-147 | match |

Every paper-side deliverable has a faithful script-side counterpart; coefficients (the recurring `1/9`, the `2`/`4` BdG slopes, the `-2`/`-2` and `Δ`-power factors, the `(1, 1/2, -1)` signature) all match. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41 | `expect_zero(K_exact - K_split)` | (1) | partial (linear regrouping; non-tautological by parenthesization, weak) |
| A2 | sympy | 45 | `expect_zero(G_exact - G_split)` | (2) | partial (same character as A1) |
| A3 | sympy | 61 | `dB0_exact - dB0_formula` (diff of c²/w²) | (3) | yes |
| A4 | sympy | 62 | `dB2_exact - dB2_formula` (diff of c²/w⁴) | (3) | yes |
| A5 | sympy | 69 | `Bcomb_exact - Bcomb_formula` | (3) | yes |
| A6 | sympy | 96 | `dDelta_expr - (W dU + U dW - 2R dR)` | (5) | yes |
| A7 | sympy | 97 | `dS_expr - (dU + dW)` | (5) | yes |
| A8 | sympy | 98 | `dG_expr - (2gu dgu + 2gw dgw)` | (5) | yes |
| A9 | sympy | 99 | `dP_expr - (gw dU + gu dR + R dgu + U dgw)` | (5) | yes |
| A10 | sympy | 100-107 | `dQ_expr - (...)` | (5) | yes |
| A11 | sympy | 122 | `dZ0_exact - (Δ dQ - Q dΔ)/Δ²` | (4) | yes |
| A12 | sympy | 123-131 | `dZ2_exact - (...)` | (4) | yes |
| A13 | sympy | 132 | `dN0_exact - (2P dP/Δ² - 2P² dΔ/Δ³)` | (4) | yes |
| A14 | sympy | 141 | `Zcomb_exact - Zcomb_formula` | (4) | yes |
| A15 | sympy | 148 | `Ncomb_exact - Ncomb_formula` | (4) | yes |
| A16-A21 | sympy | 163-164 | weak-axisymmetric K/G regrouping, λ∈{1,1/2,-1} | (6) | partial (regrouping over scalar lanes) |
| B1-B21 | mathematica | 31-147 | `expectZero[...]` mirror of A1-A21 | same | mirror of sympy |

The load-bearing checks (A3-A15) are genuine: each runs symbolic differentiation of the defining expression and compares to an independently hand-written boxed formula, so a sign or factor error in the notes would surface. The regrouping checks (A1, A2, A16-A21) are weak — both sides are the same linear combination written with different parenthesization — but they still correspond to actual paper-boxed equations and would catch a transcribed sign flip, so they are not pure tautologies.

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:1-158`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py:1-170`

**What's wrong:**
The `.wl` is a near line-by-line port of the `.py`: identical variable choreography, identical ordering of all 21 checks, identical inline comments, and — crucially — the same hand-written target formulas. Three corresponding pairs:

- sympy 64-68 `Bcomb_formula = ... 2*c*(1/w**4 + 1/(9*w**2))*dc - 2*c**2*(2/w**5 + 1/(9*w**3))*dw` vs. wl 54-57 `bCombFormula = ... 2*c*(1/w^4 + 1/(9*w^2))*dc - 2*c^2*(2/w^5 + 1/(9*w^3))*dw` — same literal target.
- sympy 134-140 `Zcomb_formula = (S/Delta**2 + Rational(1,9)/Delta)*dQ + ...` vs. wl 113-119 `zCombFormula = (s/delta^2 + 1/(9*delta))*dQ + ...` — same literal target.
- sympy 158-164 `for lam in (l20,l21,l22): K_micro = eps*lam*(...)` vs. wl 132-147 `Do[... kMicro = eps*lamVal*(...) ..., {lam, {1, 1/2, -1}}]` — same loop and same regrouped expressions.

The only genuine cross-engine signal is that each CAS runs its own `sp.diff` / `D[]` on the defining expressions (A3-A15). The right-hand "formula" targets are not re-derived independently in the second engine — they are copied. This weakens the second-engine guarantee: a wrong hand-written target formula copied identically into both scripts would pass in both.

**Why this matters:**
The second-engine policy requires both engines to derive the result independently, so a transcription error in the target formula cannot survive cross-check. Here a copied wrong target would pass in both engines, defeating the purpose. Severity is low because the load-bearing differentiations are genuinely independent per engine, so the most error-prone step (the derivative slopes) is still cross-checked.

**Required change:**
In the Mathematica script, for at least the Z/N bundles (the most factor-heavy formulas), re-derive the comparison target from the engine's own differentiation rather than restating the hand-written literal. Concretely, replace the literal `zCombFormula` (wl 113-119) and `nCombFormula` (wl 122-126) with forms assembled from `D[z2, ...]`/`D[z0, ...]`/`D[n0Expr, ...]` outputs combined per the bundle definition, so the `.wl` independently reconstructs the bundle rather than echoing the SymPy literal. (Do not change the SymPy file; the asymmetry is what restores independence.)

**Verification:**
The new `.wl` lines should build `zCombFormula`/`nCombFormula` from `D[...]` calls; `expectZero["Z obstruction bundle", ...]` and `expectZero["N obstruction bundle", ...]` must still print `= 0` and exit 0. The captured `.txt` lines 39-42 should remain `= 0` / `PASS`.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py:27`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py:5`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:26`

**What's wrong:**
The unit is Stage 171, but the banner string is `"STAGE 154 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"` in both engines (sympy line 27, wl line 26), and the SymPy docstring (line 5) says "from Stage 170." Neither label is load-bearing on an assertion — the math is correct and matches Stage 171's notes — but the mislabel propagates into the saved transcripts (`sympy .txt:11`, `mathematica .txt:11` both print "STAGE 154"), which harms traceability between the captured output and the audited unit. This is an anchoring defect, not a math error.

**Why this matters:**
A reviewer cross-referencing the transcript against the stage card sees "STAGE 154" where the card is Stage 171, which obscures whether the saved output actually belongs to this unit. It does not invalidate any check.

**Required change:**
Change the banner string in `sympy:27` and `wl:26` from `"STAGE 154 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"` to `"STAGE 171 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"`, and update the SymPy docstring at line 5 to reference Stage 171 (the notes attribute the upstream `(K_A, G_A)` pair to the renumbered upstream stage, but this card and its appendix row are Stage 171). The Mathematica final message at wl:156 already correctly says "Stage 171."

**Verification:**
Re-run capture; `sympy .txt:11` and `mathematica .txt:11` should both print "STAGE 171 ...". No assertion line changes; all `= 0` / `PASS` lines remain.

## Independent-derivation check (Mathematica)

The `.wl` derives the load-bearing slopes independently (it runs its own `D[c^2/w^2, c]`, `D[(q*s - gSym*delta)/delta^2, ...]`, etc.), which is genuine. But the comparison targets (`bCombFormula`, `zCombFormula`, `nCombFormula`, and the regrouped `kMicro`/`gMicro`) are copied literals identical to the SymPy file rather than re-derived, and the overall structure is a line-by-line transliteration. See F1. Net: partially independent — the differentiation layer is independent, the target-formula layer is not.

## Engine cross-check

Both engines emit identical residuals: every one of the 21 checks prints `= 0` in both transcripts (sympy `.txt:13-33`, mathematica `.txt:13-54`), both exit 0, both `Status: PASS`. The Mathematica transcript additionally prints `PASS:` per check. No numeric or symbolic disagreement. `engines_agree: true`.

## Verdict justification

Verdict is `findings` with two low-severity issues, neither blocking. I attacked the math on multiple fronts and it held: I re-derived `δZ_2 = ∂Z_2/∂Q·dQ + ∂Z_2/∂S·dS + ∂Z_2/∂G·dG + ∂Z_2/∂Δ·dΔ` from `Z_2 = (QS-GΔ)/Δ²` and confirmed the `G/Δ² - 2QS/Δ³` Δ-derivative term and the `-1/Δ` G-term; I confirmed `δN_0 = 2P dP/Δ² - 2P² dΔ/Δ³` from `P²/Δ²`; I confirmed `N_A = (P_0/Δ)dQ + (2P/Δ²)dP - (P_0 Q/Δ² + 2P²/Δ³)dΔ` by substituting `δZ_0 = dQ/Δ - Q dΔ/Δ²`; I checked every primitive variation (`δΔ, δS, δG, δP, δQ`) against the notes' boxed forms (exact match, including that `S` and `G` correctly omit the variables they do not depend on, and that `P` correctly omits `dW`); and I confirmed the `(1, 1/2, -1)` signature and the `1/9` weights match the notes and appendix. The variable-independence trap is clean (every `sp.diff(EXPR, VAR)` has EXPR genuinely depending on VAR). The symbol-domain assumptions (`nonzero=True` on `P0, c, w, U, W, R, Delta, P`) are justified by the physical setup (frequencies and determinants are nonzero). Outputs are fresh (both `.txt` mtimes postdate their scripts). The only defects are the structural transliteration of the second engine (F1) and the stale "STAGE 154"/"Stage 170" labels (F2). Paper alignment is exact on all six deliverables, so no `paper_misalignment` and no stop-cold.

## Self-test notes

I checked the variable-independence trap (all `sp.diff`/`D` targets genuinely depend on the differentiation variable; `S`, `G`, `P` correctly omit their non-dependencies), the factor/sign trap on every Δ-power term in `δZ_2`, `δN_0`, `Z_A`, `N_A` (all match the notes), and the symmetry/domain trap (no integrals; nonzero-denominator assumptions are physically justified). For the F1 directive I confirmed the proposed Mathematica re-derivation of `zCombFormula`/`nCombFormula` from `D[...]` would still reduce to 0 because it reconstructs the same bundle the notes box; for F2 I confirmed the banner/docstring edits touch no assertion and do not introduce a paper_misalignment (Stage 171 is the correct unit per card and appendix).
