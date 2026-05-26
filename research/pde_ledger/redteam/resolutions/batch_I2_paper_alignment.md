---
batch: I.2
batch_label: "Part I.2 — Maxwell bridge, parent throat action, reduced one-port"
stage_range: "013-023"
items_count: 6
created: 2026-05-25
status: awaiting_codex_recommendations
---

# Batch I.2 paper-misalignment resolutions

Six v2 audits surfaced `paper_misalignment` findings. Each is described below with the relevant file:line evidence and (a/b/c/skip) options. Codex: please read the cited paper and script files, then fill the empty `## Recommendation` block under each question with `direction: a|b|c|skip` and a one-paragraph rationale citing file:line. If you do not have a high-confidence answer for a question, write `direction: skip`.

This batch's pattern (like I.1) is dominated by `paper_missing_script_claim`: the EM-projected scripts (stages 013-020) were ported from `notes/em_projected/`, but the paper cards under `paper/stages/` were written later as compact summaries. The natural directions are (a) expand paper card to enumerate the extra script checks; (b) trim scripts to paper scope and migrate extras to whatever stage actually owns them; or (c) add a light paper acknowledgement.

---

## Q1 — Stage 013 (paper_missing_script_claim, low)

**Files:**
- paper: `paper/stages/stage_013.tex` (whole card; n_0 only; no `deltaP_2/deltaP_4/n_2/n_4` mention)
- sympy: `scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:106-151`
- mathematica: `mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl:126-166`
- appendix row: `paper/appendices/stage_appendix_part01.tex` row for 013 ("Mouth-local Taylor map for projected source and flux corrections" — too broad to anchor specific deliverables)

**What the script verifies that the paper does not declare:**

1. `deltaP_2` and `deltaP_4` defined and tested via substantive identities `∂(deltaP_2)/∂G_w' = -2P/(D_0*Delta^2)` and a nonzero `∂(deltaP_4)/∂G_w'`. Neither symbol appears in the paper card.
2. "Mechanism sieve": 2x2 linear-solve argument that `(K_1, H_even)|_{Sx=Hx=0}` and `(K_1, H_even)|_{Qx=Dx=0}` admit only the trivial pair. No paper-side mention.
3. `n_2, n_4` derived in Mathematica and asserted in SymPy. Paper card mentions only `n_0`.

**Options:**

- (a) Expand `paper/stages/stage_013.tex` Output to declare `deltaP_2, deltaP_4, n_2, n_4` and the mechanism sieve as deliverables of stage 013.
- (b) Trim the scripts to the four equation-numbered deliverables (z_0, z_2, z_4, n_0) and move the `deltaP` / sieve / `n_2 / n_4` machinery to stage 014 (where K_1, H_even are formally defined).
- (c) Add 1-2 sentences to the paper card explicitly acknowledging that the audit also exercises the secondary Taylor data that feeds stage 014.
- (skip) Codex does not have a high-confidence answer.

## Recommendation

direction: b
rationale: |
  Cross-stage check confirms duplication, not "secondary downstream data": stage 010's paper card now explicitly exports `δP_2` (`paper/stages/stage_010.tex:42-48` eq:stage010-dP2) and `δP_4` (eq:stage010-dP4) after the I.1 v2 paper sweep; stage 014's script (`scripts/moving_throat_pde_stage014_*_sympy_audit.py:77-89,129-138`) verifies the identical mechanism sieve (`qd_matrix`, `sh_matrix`, `qd_only=[{Qx:0,Dx:0}]`, `sh_only=[{Sx:0,Hx:0}]`). Codex initially recommended (c) acknowledgement, but per the "each stage builds on prior" principle, the cleanest fix is to trim the duplicates: stage 013 keeps only its Taylor map (z_0, z_2, z_4, n_0) + Xi_load (eq:stage013-xi). The δP_2/δP_4 and sieve assertions are redundant with stages 010 and 014 respectively. Keep n_2, n_4 in stage 013 if they support the Taylor map narrative (no paper card explicitly owns them yet).

## Apply: done
files_changed:
  - "scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:100-106: removed δP2/δP4 and mechanism-sieve assertions; kept Xi_load and Taylor coefficient map."
  - "mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl:113-119: removed δP2/δP4 and sieve blocks; kept n0/n2/n4 chain-rule checks and Xi_load d/dPprime."
sanity_check: "sympy + mathematica exit 0"
destination_verified: "yes — scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:68-92 asserts δP2/δP4; scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:60-115 asserts the mechanism sieve and compensation denominators."
notes: "No paper edits needed."

---

## Q2 — Stage 014 (paper_missing_script_claim, low)

**Files:**
- paper: `paper/stages/stage_014.tex:30-33` (Output paragraph)
- sympy: `scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:67,71-75,93-95,126-128`
- mathematica: `mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl:113,117-126,138-140`
- appendix row: "Gate conditions for carrying mouth-local projected data into the grouped bundle"

**What the script verifies that the paper does not declare:**

The paper card Output says only: *"Stage 014 exports the mechanism sieve: projected EM mouth data can tune the even gates only away from the degeneracy eq:stage007-compensation-denominator."* The card body exhibits exactly two equations (`K_1`, `H_even`) and one denominator identity. The scripts additionally test:

- `Xi_load = n_0/N_0 + z_0/D_0` and its `∂/∂P' = 2/P` identity
- `δP_2` and `δP_4` bundle-level transport definitions and their `∂/∂G_w' = -2P/(D_0 Δ^2)` and nonzero coefficient assertions
- Isotropic compatibility (`Compat`) printed-only check
- "5PN bottleneck" / "constant-prefactor transport" framing in the docstring with no paper-side counterpart

**Options:**

- (a) Expand the paper card Output to enumerate `Xi_load, δP_2, δP_4` as further deliverables of stage 014.
- (b) Trim the scripts to test only the even-gate sieve content the paper actually states; migrate transport-ledger assertions to a separate stage.
- (c) Add 1-2 sentences to the paper card explicitly acknowledging the audit also exercises constant-prefactor transport coefficients used downstream.
- (skip) Codex does not have a high-confidence answer.

## Recommendation

direction: b
rationale: |
  Same principle as Q1. Codex correctly identified that stage 010 owns `δP_2`/`δP_4` (`paper/stages/stage_010.tex` eq:stage010-dP2/dP4) and stage 013 owns `Xi_load` (`paper/stages/stage_013.tex:35-37` eq:stage013-xi) — but acknowledgement (c) keeps the duplicated assertions cluttering stage 014. Cleaner fix: trim `Xi_load` (line 67), `δP_2`/`δP_4` (74-75, 126-128), and `Compat` (printed) from stage 014's scripts; keep only the K_1/H_even mechanism sieve and compensation denominator that stage 014 actually owns per its paper card (`paper/stages/stage_014.tex:12-33`).

## Apply: done
files_changed:
  - "scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:56-115: trimmed Xi_load, δP2/δP4, Compat, and prefactor readouts; kept K1/H_even sieve plus denominator checks."
  - "mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl:39-130: trimmed XiLoad, δP2/δP4, and N-prefactor transport; kept K1/He solves, Jacobians, and denominator mutation."
sanity_check: "sympy + mathematica exit 0"
destination_verified: "yes — scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:68-139 asserts δP2/δP4 and compatibility transport; scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:100-103 asserts Xi_load d/dPprime."
notes: "The stage 014 SymPy printout now reports only the conservative gate bridge."

---

## Q3 — Stage 015 (paper_missing_script_claim, high)

**Files:**
- paper: `paper/stages/stage_015.tex:44-46` (Output paragraph), and whole card
- appendix row: `paper/appendices/stage_appendix_part01.tex:52` ("Parent throat-action packet and projection/reduction status boundary")
- sympy: `scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:103-208`
- mathematica: `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:104-196`
- sympy docstring (line 2) references nonexistent `step_13_parent_throat_action_master_notes.md`

**What the script verifies that the paper does not declare:**

Paper card Output says only the parent throat action promotion and the K_eta exact-quadratic-recovery formula. Yet ~6 of 12 sympy asserts (A9-A19) and ~21 of 28 mathematica checks (M8-M28) cover:

1. Wall-only specializations of full even-channel gates `K_1 = -dM + dK/9` and `H_even = (2/3) dM - dK/27`, with substitution `{B_{01}, B_{21}, B_{41}, Z_{01}, Z_{21}, Z_{41}} → 0` (sympy:123 / math:109), Jacobian determinant `1/27` (sympy:174 / math:151), Gaussian overlap integrals `dMoverlap = √(π/3)` and `dKoverlap = 23√π/(3√3)` (math:126-127).
2. Real-Y20 overlap ratios `{1, 1/2, -1}` for `m = 0, 1, 2` (sympy:195-200 / math:175-184).
3. Grouped trace identities `xbar = x0`, `bx = 3*ax` with `lam20=1, lam21=1/2, lam22=-1` (sympy:202-208 / math:186-196).

**Severity is "high"** because this is roughly half of each script — large fraction of script effort with no paper anchor.

**Options:**

- (a) Expand `paper/stages/stage_015.tex` to declare wall-only gates, Y20 ratios, and grouped trace as additional deliverables of stage 015, and restore (or rewrite) `notes/stages/moving_throat_pde_stage015_*.md` to back them.
- (b) Trim scripts to A1-A8 / M1-M7 (only the K_eta exact-quadratic block, which is the deliverable the paper actually claims); migrate the wall-only/Y20/grouped blocks to whatever stage(s) genuinely own them.
- (skip) Codex does not have a high-confidence answer.

## Recommendation

direction: b
rationale: |
  The wall-only gate obstruction and the real-\(Y_{20}\) lane ratios are already explicit stage 017 outputs (`paper/stages/stage_017.tex:13-37`), while stage 015 only claims the parent action promotion and \(K_\eta\) quadratic recovery (`paper/stages/stage_015.tex:32-46`).  The stage 015 wall-only/Y20/grouped blocks (`scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:103-208`) therefore duplicate later-stage ownership and should be trimmed or migrated rather than promoted into stage 015.

## Apply: done
files_changed:
  - "scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:1-89: removed wall-only, Gaussian-overlap gate, Y20, and grouped-trace blocks; kept IBP, boundary, and K_eta checks."
  - "mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:83-105: removed M4-M9 wall-only/Y20/grouped blocks; kept M1-M3 parent-action checks."
sanity_check: "sympy + mathematica exit 0"
destination_verified: "yes — scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py:40-110 asserts Y20 lane ratios, b=3a, and wall-only obstruction; scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:151-164 also asserts Y20 lane ratios."
notes: "Removed now-unused SymPy gaunt/grouped helpers."

---

## Q4 — Stage 018 (paper_missing_script_claim, medium)

**Files:**
- paper: `paper/stages/stage_018.tex:26-28` (Output paragraph), body 14-22 (only M_Σ, K_Σ branch integrals)
- appendix row 58: "Bundle-level parent-action identities used by the projected electromagnetic response" (plural "identities" hints at more, never enumerated)
- sympy: `scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:1-2, 25-88`
- mathematica: `mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl:1-15, 24-122`
- sympy docstring references nonexistent `step_16_parent_throat_action_bundle_master_notes.md`

**What the script verifies that the paper does not declare:**

Paper Output says only M_Σ, K_Σ branch integrals. Scripts additionally verify:

1. An algebraic identity about coefficients of a quartic pole expansion (one-pole numerator).
2. Two `KSigma` closure conditions and their compatibility (one-pole vs. normalization closures).
3. A 2x2 even-gate determinant equal to `1/27` and closed-form wall-stiffness/wall-inertia slopes.
4. A residual `Xi1` after applying those slopes.

Only the trailing Gaussian integrals A18/A19 / M8 directly correspond to the paper's M_Σ, K_Σ deliverables, and they do so only for a specific concrete profile (`μ_η = T_w = K_η + 6 T_Ω = 1`, `β = exp(-w²/2)`).

**Options:**

- (a) Expand the paper card and restore `notes/stages/moving_throat_pde_stage018_*.md` to declare the four extra claim families (one-pole closure, gate determinant, wall slopes, Xi_1 residual) as 018 deliverables.
- (b) Move the four extra claim families to a different stage (candidates: 015-017 or 021 — "Reduced Maxwell/mixed one-port normal form"); script becomes just the bridge identity check.
- (c) Prune scripts to just the M_Σ, K_Σ bridge identity (drop families 1-4 entirely).
- (skip) Codex does not have a high-confidence answer.

## Recommendation

direction: b
rationale: |
  Stage 018 should own the parent-wall bridge integrals \(M_\Sigma,K_\Sigma\) (`paper/stages/stage_018.tex:12-28`), while the one-pole/normalization compatibility is already stage 019 content (`paper/stages/stage_019.tex:13-46`) and the even-gate slope solve plus residual \(\Xi_1\) is stage 020 content (`paper/stages/stage_020.tex:38-55`).  The extra stage 018 script families (`scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:25-88`) should move to those owning stages or be removed from stage 018 once confirmed there.

## Apply: done
files_changed:
  - "scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:20-30: trimmed one-pole, compatibility, even-gate, slope-solve, and Xi1 blocks; kept Gaussian MΣ/KΣ branch integral checks."
  - "mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl:14-30: trimmed M1-M7 closure/slope/Xi checks; kept Gaussian inertia/stiffness integral checks."
sanity_check: "sympy + mathematica exit 0"
destination_verified: "yes — scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py:35-56 asserts one-pole and compatibility closure; scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:35-53 asserts gate determinant, wall-slope solve, and Xi1 residual."
notes: "No receiving-script changes were needed."

---

## Q5 — Stage 020 (paper_missing_script_claim, medium)

**Files:**
- paper: `paper/stages/stage_020.tex:1-56` (no Y20 / Gaunt / lane ratio mention anywhere in card)
- appendix row: `paper/appendices/stage_appendix_part01.tex:62` ("Weak-axisymmetric packet exported from the parent throat-action bundle" — no Y20 mention)
- sympy: `scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:45-50`
- mathematica: `mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl:34-67`

**What the script verifies that the paper does not declare:**

Both engines compute and assert Y20-squared lane overlap ratios `λ_{20} = 1, λ_{21} = 1/2, λ_{22} = -1` and vanishing of same-sign m-cross terms. The paper card and appendix row make no reference to spherical harmonics, Y20 overlaps, Gaunt integrals, or these lane ratios. Paper card content is purely the 1D wall-slope solve in the `β_2(w)` profile and identities `eq:stage020-wall-slopes` through `eq:stage020-residual-xi`.

The "weak-axisymmetric" content of the stage title suggests the Y20 ratios are conceptually load-bearing somewhere — but the paper does not say where.

**Options:**

- (a) Amend `paper/stages/stage_020.tex` to add a paragraph or equation citing the Y20 lane ratios as part of the weak-axisymmetric deliverable.
- (b) Remove the Y20 block from both scripts (sympy:45-50, wl:34-67) — they belong to a different stage.
- (skip) Codex does not have a high-confidence answer.

## Recommendation

direction: b
rationale: |
  The \(Y_{20}\) ratios are not used by the stage 020 wall-slope solve after they are asserted (`scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:45-70`), and stage 020's paper deliverable is the algebraic weak-axisymmetric packet from \(\delta M_\Sigma,\delta K_\Sigma\) through the compensated \(\Xi_1\) (`paper/stages/stage_020.tex:12-55`).  The angular lane law is already paper-owned upstream by stage 017 (`paper/stages/stage_017.tex:13-37`) and elsewhere by the broader grouped-bundle/overlap stages, so the duplicate Y20 block should be removed from stage 020.

## Apply: done
files_changed:
  - "scripts/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_sympy_audit.py:10-53: removed Y20 helper/import and lane assertions; kept determinant, wall-slope solve, and Xi1 checks."
  - "mathematica/moving_throat_pde_stage020_parent_throat_action_weak_axisym_packet_mathematica_audit.wl:6-104: removed GauntIntegral helper and angular M1 block; kept wall-slope packet checks."
sanity_check: "sympy + mathematica exit 0"
destination_verified: "yes — scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:151-164 and scripts/moving_throat_pde_stage017_parent_throat_action_weak_axisym_sympy_audit.py:40-45 assert the Y20 lane ratios."
notes: "No receiving-script changes were needed."

---

## Q6 — Stage 021 (script_missing_paper_claim, medium)

**This is the inverse direction**: paper claims more than scripts assert.

**Files:**
- paper: `paper/stages/stage_021.tex:77-81` (Output paragraph)
- paper: `paper/stages/stage_021.tex:71-75` (`eq:app-stage021-wall-odd`)
- sympy: `scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py:238-244` (`Dcorr` computed and printed, no assertion)
- mathematica: `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:135, 140` (`dCorr` computed and printed, no assertion)

**What the paper claims that the script does not verify:**

The paper's Output paragraph explicitly enumerates three exports: (1) reduced one-port self-energy, (2) transfer factor, (3) **wall-level outgoing quadrupole coefficient** `δD_2^(odd)(ω) = -i N_2(0) (a^5/(27 c_s^5)) ω^5 + O(ω^7)`. The numeric coefficient `a^5/(27 c_s^5)` is substituted in for `Γ_5^port`.

In both scripts, Section III computes:
```
Dcorr = -I * Gamma_port * omega^5 * N0   (sympy)
dCorr = -I * gammaPort * omega^5 * n0    (mathematica)
```
`Gamma_port` is a free symbolic constant; it is never specialized to `a^5/(27 c_s^5)` and `Dcorr` is never asserted against the paper's closed form. Sections III and IV separately verify N_l(0) and Γ_5^port; the composition is asserted only in prose, not in code.

The fix is algebraically trivial (one substitution + assertion in each script).

**Options:**

- (a) Add the composed assertion to both scripts (sympy and mathematica) — paper is right, scripts have a verification gap. **Recommended:** the auditor flagged paper-says-three-the-scripts-prove-two-and-a-half as the natural read.
- (b) Trim the paper Output to two exports (drop the wall-level odd-coefficient claim) — the composition is trivial enough that the paper can describe it in prose without claiming it as a deliverable.
- (skip) Codex does not have a high-confidence answer.

## Recommendation

direction: a
rationale: |
  The paper and the stage note both make the composed wall-level odd quadrupole coefficient a headline export (`paper/stages/stage_021.tex:60-81`, `notes/stages/moving_throat_pde_stage021_reduced_one_port_normal_form.md:229-244`).  The scripts separately verify \(N_l(0)\) and \(\Gamma_5^{\rm port}=a^5/(27c_s^5)\), but only print the unspecialized `Dcorr`/`dCorr` (`scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py:227-244`), so adding the composed assertion is the correct low-risk fix.

## Apply: done
files_changed:
  - "scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py:238-245: added composed δD_wall^(odd) assertion after substituting Γ_port = a^5/(27 c_s^5)."
  - "mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:123-139: added matching dCorr composed coefficient assertion with radius/cS."
sanity_check: "sympy + mathematica exit 0"
destination_verified: "yes — no destination check required; this is a script-side paper-claim assertion add."
notes: "No paper edits needed."

---

summary: 6 answered, 0 skipped

## Apply notes

Overall summary: 6 applied, 0 revised, 0 blocked.

Per-stage script exit status: stage 013 sympy + mathematica exit 0; stage 014 sympy + mathematica exit 0; stage 015 sympy + mathematica exit 0; stage 018 sympy + mathematica exit 0; stage 020 sympy + mathematica exit 0; stage 021 sympy + mathematica exit 0.

Blocked items for Trevor: none.
