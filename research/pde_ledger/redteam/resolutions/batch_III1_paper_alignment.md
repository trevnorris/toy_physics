---
batch_id: III.1-v2
created: 2026-05-26
auditor_version: v2 (paper-grounded)
total_user_gate_items: 4
codex_role: math-authority + cross-stage research; provide direction (a/b/c/skip) + rationale per Q
user_role: review Codex recommendations, approve or redirect before any apply session
---

# Batch III.1 v2 — user-gate item consolidation

The v2 audit of stages 037–048 (Continuum kernel, generalized branch, rank-2) produced **3 paper_misalignment items** (043 F2, 045 F4, 046 F1) plus **1 insufficient_verification item with user gate** (045 F3). Script-side findings (covered separately via fix_loop) are NOT in this file.

For each item Codex must fill in the `## Recommendation` block with `direction:` (a/b/c/skip) plus rationale citing file:line evidence. Codex MUST NOT apply any change in this session — recommendations only.

---

## Q1 — Stage 043 F2 (paper_misalignment / notes_contradicts_script, low)

**Item:** Mathematica computes `D_phi = Det[{{y0, y1}, {kappa0, kappa1}}] = y0*kappa1 - y1*kappa0`. The notes (`notes/stages/moving_throat_pde_stage043_support_direction_extraction.md:119-121`) and the SymPy script (`scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:77-78`) define `D_phi := kappa_0 y_1 - kappa_1 y_0` — opposite sign. The Mathematica `dPhiExpected` is also sign-flipped to match, so the internal check passes; but the printed transcript value has the OPPOSITE sign from the notes' definition and from SymPy's transcript. `engines_agree: false` (sign only).

**Audit finding location:** `redteam/reports/stage_043.md` § F2

**Context:** The vanishing condition `D_phi = 0` (the actual physical claim — source-tied iff `sigma_0 = 0` or `delta_U = 0`) is unaffected by sign. The math is the same; only the printed/named quantity differs. The notes/sympy definition is `D_phi := kappa_0 y_1 - kappa_1 y_0` — this puts `(kappa_0, kappa_1)` in the first row of the 2×2 determinant.

**Options:**
- **(a) Fix Mathematica** to use the notes/SymPy convention. Change `Det[{{y0, y1}, {kappa0, kappa1}}]` → `Det[{{kappa0, kappa1}, {y0, y1}}]` (rows swapped) and remove the sign-flip in `dPhiExpected`. Printed transcript value then matches SymPy and notes.
- **(b) Fix notes + SymPy** to use the Mathematica convention. Change notes' `D_phi := kappa_0 y_1 - kappa_1 y_0` → `D_phi := y_0 kappa_1 - y_1 kappa_0` (equivalent to `Det[{{y0, y1}, {kappa0, kappa1}}]`); fix SymPy to match.
- **(c) Acknowledge with documentation** — leave both as-is, add a one-line script comment + notes paragraph explaining the sign-flip convention difference.
- **skip**

## Recommendation

direction: a
rationale: The downstream interface from Stage 043 is the replacement tuple `q=tR_U`, `r=tR_phi`, not the determinant name: Stage 043 says Stage 044 uses those replacements (`paper/stages/stage_043.tex:22-36`), and Stage 044 declares the same inputs and formulas only in `R_U,R_phi` (`paper/stages/stage_044.tex:7`, `paper/stages/stage_044.tex:23-29`). Stage 045 then consumes the tracking surface `R_phi=R_U=R_tr` (`paper/stages/stage_045.tex:19-23`), while Stage 046 consumes only `F_tr`, `G_tr`, and `R_tr<1` (`paper/stages/stage_046.tex:7`). With no downstream `D_phi` sign convention visible in those consumer surfaces, the established kappa-first direction-splitting convention should control: Stage 039 defines `D_dir := kappa_0 z_1 - kappa_1 z_0` (`notes/stages/moving_throat_pde_stage039_split_u_sector.md:168-172`) and computes it that way (`scripts/moving_throat_pde_stage039_split_u_sector_sympy_audit.py:119-120`), and Stage 043 notes/SymPy mirror that convention for `D_phi` (`notes/stages/moving_throat_pde_stage043_support_direction_extraction.md:119-121`, `scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:76-78`). The Mathematica file is the outlier because it swaps the determinant rows, `Det[{{y0, y1}, {kappa0, kappa1}}]`, and expects the positive sign (`mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:52-53`).

---

## Q2 — Stage 045 F3 (insufficient_verification, medium, user-gate)

**Item:** Notes section 6 of `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md` claims: *"Because the coherent local kernel lands on the tracking surface, the Stage-27 normalization residual also collapses to the Stage-23 tracking form: `R_target = F_tr(xi, delta; R_tr)`"*. The script's current check at `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:172-189` and `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:125-136` verifies only that algebraic substitution `lam0 → 2/9` in a hand-written generic-`lam0` `F_track` expression gives the notes' `F_tr` closed form — a pure tautology by construction (the generic form was written so it would specialize this way). The substantive claim — that the Stage-27 normalization residual on the tracking branch reduces to `F_tr` — is not exercised.

**Audit finding location:** `redteam/reports/stage_045.md` § F3 (`## Resolve before fix_loop` block in directive)

**Cross-stage research needed:** "Stage 27" in the notes' language is **Stage 044** in current paper numbering ("Continuum-Selected Rank-2 Closure"). The relevant residual is the normalization residual on the rank-2 closed branch — `paper/stages/stage_044.tex` (`eq:app-stage044-Rtgt-Ftr`?) and `scripts/moving_throat_pde_stage044_continuum_selected_rank2_closure_sympy_audit.py` should be the source. Codex needs to identify the canonical Stage-044 residual expression and confirm that substituting tracking conditions (`R_phi = R_U`, `lam_0 = 2/9`) into it yields `F_tr` as stated in 045's notes.

**Options:**
- **(a) Import Stage-044's `R_target` residual into stage 045**: lift the canonical residual expression from `scripts/moving_throat_pde_stage044_..._sympy_audit.py`, substitute tracking conditions, verify it reduces to `F_tr` per notes. This is the math-substantive direction.
- **(b) Keep the hand-written F_track + relabel the assertion as a self-consistency check** rather than a tracking-collapse theorem. Drop the load-bearing claim from this stage's deliverables.
- **(c) Move the tracking-collapse theorem to a downstream stage** that has the Stage-044 residual in scope (likely 046 or 047). Stage 045 retains only the closed-form `F_tr` form check.
- **skip**

## Recommendation

direction: a
rationale: The canonical Stage 044 normalization gate is `R_target=F_cont(xi_phys)` in the paper (`paper/stages/stage_044.tex:31-35`), and the notes define the residual as `R_target - F_cont(xi_phys) = 0` after displaying `F_cont` (`notes/stages/moving_throat_pde_stage044_continuum_selected_rank2_closure.md:147-166`). The Stage 044 SymPy audit defines the same `F_cont` expression (`scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90`) and already substitutes `R_phi=R_U` to assert the generic tracking collapse (`scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:140-146`). Setting `lambda_0=2/9` in that tracked expression gives the Stage 045 notes' displayed `F_tr` form (`notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:218-220`): the two numerator factors each acquire a `1/9` scale and the squared denominator acquires a `1/81` scale, leaving the displayed overall `1/81`. Therefore (a) is mechanically applicable: import Stage 044's `F_cont/R_target` residual into Stage 045 and assert the tracking plus D/N specialization, instead of only checking the hand-written `F_track` specialization (`scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:178-188`).

---

## Q3 — Stage 045 F4 (paper_misalignment / notes_contradicts_script, low — cosmetic)

**Item:** SymPy docstring/banner/final-print at `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3,31,191` still says "Stage 28" / "STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT" / "All Stage-28 symbolic checks passed.". Mathematica banner at `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26` says "STAGE 028 — COHERENT LOCAL TRACKING" while line 139 footer correctly says "Stage 045". Paper card and filenames are "Stage 045".

**Audit finding location:** `redteam/reports/stage_045.md` § F4

**Context:** Same kind of legacy-label drift as I.2 v2's stage 022 (deemed cosmetic-only, not paper_misalignment), I.1 v2's stage 028 (same handling), and II.1 v2's stage 029 F1 (Q1=(a) cosmetic relabel approved by user). Auditor explicitly noted "direction of fix is unambiguous (script side updates to match paper); directive applies four single-line cosmetic edits without user gate."

**Options:**
- **(a) Re-categorize as script-side cosmetic** — relabel docstring/banner/final-print under fix_loop without paper change. (Matches I.2/II.1 v2 precedent.)
- **(b) Relabel scripts AND propagate "Stage 045" into the notes file** if notes also reference "Stage 28". (Notes are read-only for Codex but user could approve.)
- **(c) Acknowledge as known limitation** — leave intact, add a one-line comment noting the rename.
- **skip**

## Recommendation

direction: b
rationale: Treat this as cosmetic renumbering drift, consistent with precedent: Stage 022 explicitly called stale labels documentation drift rather than `paper_misalignment` (`redteam/reports/stage_022.md:129`), Stage 028 called the same kind of label mismatch a renumbering artifact rather than a math/paper disagreement (`redteam/reports/stage_028.md:164`), and the Stage 029 F1 resolution recommended script-side cosmetic relabeling (`redteam/resolutions/batch_II1_paper_alignment.md:34-37`). Here the stale label is not script-only: SymPy says Stage 28 in its docstring, banner, and final print (`scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3`, `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:31`, `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:191`), Mathematica says `STAGE 028` in its banner while its footer already says Stage 045 (`mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26`, `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:139`), and the notes heading also says "after Stage 28" (`notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:232`). Since option (b) exactly covers relabeling scripts plus this notes occurrence, do (b) as a cosmetic relabel only, with no math or paper-card change.

---

## Q4 — Stage 046 F1 (paper_misalignment / notes_contradicts_script, medium)

**Item:** Notes file `notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md` has 5 coefficient typos in the auxiliary positivity polynomials `P_R`, `P_1`, `P_2`:

| Polynomial | Notes value | Script + both engines |
|---|---|---|
| `P_R`: `R delta^3` term | 230 | 162 |
| `P_R`: `R delta xi^2` term | 230 | 162 |
| `P_1`: `R delta^2 xi` term | 248 | 180 |
| `P_1`: `delta^3` term | 230 | 162 |
| `P_2`: `R^2 xi^4` term | 237 | 220 |

Notes locations: `notes/stages/moving_throat_pde_stage046_..._md:87-91, 131-133, 136-139`. Both SymPy and Mathematica independently compute these from `D[F_tr, R]` and `F_flat - F_tr` and arrive at the script-asserted values. The paper card `stage_046.tex` does NOT quote these coefficients, so the disagreement is notes ↔ scripts only. Positivity argument works on either side (all coefficients positive), so the qualitative claim is unaffected — but coefficient values must agree across math sources.

**Audit finding location:** `redteam/reports/stage_046.md` § F1

**Context:** Same character as II.1 v2's Q3 (stage 035 paper polynomial coefficient fix). Independent quotient-rule / `Factor[F_flat - F_tr]` derivation should confirm which side is right.

**Options:**
- **(a) Fix notes to match scripts** (math-evident: both engines independently agree with the scripts). Update notes lines 87-91, 131-133, 136-139 to use the script-side coefficients. Scripts unchanged.
- **(b) Reverse the scripts** if the notes turn out to be derived from a different (e.g., grouped or rescaled) form of `F_tr`. Auditor notes this is "very unlikely" given both engines independently derive the script values.
- **(c) Add a paper-side note** acknowledging the discrepancy. Not advisable.
- **skip**

## Recommendation

direction: a
rationale: The script's source formula is `F_tr` (`scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:44-48`), and its expected derivative factor uses `P_R` with `162 R delta^3` and `162 R delta xi^2` (`scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:67-78`), while the notes print `230` for both coefficients (`notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:87-91`). For the coefficient check, write `A=9 delta+9 xi+2 R^2 xi`, `B=9 delta+9 xi+2 R xi`, and `C=9 delta^2+18 delta xi+(9+2 R^2)xi^2`; differentiating `A^2 B^2/C^2` gives `4 xi A B [(A+2 R B)C - 2 R xi A B]/C^3`, so `P_R=(A+2 R B)C - 2 R xi A B` after the common factors in the script derivative (`scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:82-84`). In `A+2 R B = 9 delta+9 xi+18 R delta+18 R xi+6 R^2 xi`, the only contribution to the monomial `R delta^3` is `(18 R delta)(9 delta^2)=162 R delta^3`, and the subtraction term has an explicit `xi`, so it cannot contribute to `R delta^3`. Mathematica independently prints the same `P_R` coefficients, including `162*delta^3*r` and `162*delta*r*xi^2`, in `dF_tr/dR` (`mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt:14`), and its `F_flat - F_tr` factor line confirms the other disputed script coefficients `180`, `162`, and `220` (`mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt:20`). The paper card states monotonicity and the residual comparison but does not quote these auxiliary polynomial coefficients (`paper/stages/stage_046.tex:13-36`), so fix the notes to match the engines.

---

## End of questions

After Codex fills in `## Recommendation` blocks above, halt and present to user for review. Apply session is a separate Codex invocation with explicit per-question scope.

## Apply log

### Q1 applied
- direction: a
- files modified: `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:52-53`
- destination_verified: yes — `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:52-53`
- post-edit checks: `math -script mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl` exit 0
- notes: The checked file had a positive expected value before the row swap; after switching to the kappa-first determinant, the expected value must be negative to match the notes/SymPy invariant. Mathematica verified `D_phi - expected = 0`.

### Q2 applied
- direction: a
- files modified: `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:178-207`; `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:125-156`
- destination_verified: yes — Stage 044 source grep confirmed `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:83,87,140,146` and `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:81,85,130,138`
- post-edit checks: `python3 scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py` exit 0; `math -script mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl` exit 0
- notes: Kept a generic Stage-044 tracking-collapse assertion as a subsidiary anchor, then made the D/N `F_tr` assertion use the imported Stage-044 `F_cont` expression.

### Q3 applied
- direction: b
- files modified: `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3,31,209`; `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26`; `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:232`
- destination_verified: yes — `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3,31,209`; `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26`; `notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md:232`
- post-edit checks: `python3 scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py` exit 0; `math -script mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl` exit 0
- notes: Only metadata labels and the approved notes heading were relabeled; content references to the legacy Stage-27 theorem were left unchanged.

### Q4 applied
- direction: a
- files modified: `notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:90,132-133,137`
- destination_verified: yes — `notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:90,132-133,137`
- post-edit checks: n/a — notes only
- notes: Scripts and paper card were not edited; the notes coefficients now match the already-verified script/engine values.
