# Batch I.1 v2 audit — paper alignment resolutions

This file holds the 7 questions raised by the v2 (paper-grounded) red-team audit pass of batch I.1 (stages 001–012). Each question describes a script ↔ paper disagreement and asks for a direction.

## How to answer

Each question has an empty `## Recommendation` block. Fill it in by writing:

```yaml
direction: a    # or b, c, skip
rationale: |
  1-3 sentences explaining the choice.
```

If you do not have a high-confidence answer for a question, write `direction: skip` with a one-line note saying why (e.g., "need to check upstream convention in stage X first" or "depends on physical interpretation I'm not sure about"). **Do not guess.** Skipped questions will be reviewed manually.

**Only edit this file.** Do not touch any script, paper, or notes file. Recommendations are applied by a separate workflow after review.

When done, write `summary: X answered, Y skipped` at the end of the file (X+Y = 7).

---

## Q1 — Stage 001 F1: modal-wall PDE source sign

**Files to read for context:**
- Paper: `paper/stages/stage_001.tex` (~line 161)
- Notes: `notes/stages/moving_throat_pde_stage001_geometry_lift.md` (~line 361)
- SymPy: `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:189,191`
- Mathematica: `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:188,191`
- Full audit report: `redteam/reports/stage_001.md`

**The disagreement:**
Paper and notes write `+(S_lm^(psi,A) + f_lm^ext)` on the RHS of the modal-wall PDE. The SymPy and Mathematica scripts encode the source as `-q · source_total` in the Lagrangian, which produces `-(S_lm + f_lm^ext)` on the equation RHS. The script's summary line claims it verifies the paper's positive-RHS form, but the assertion exercises the opposite sign.

**Options:**
- (a) Paper is correct → flip script source-coupling sign at the 4 line numbers above
- (b) Script is correct → flip paper text to `= -(S_lm^(psi,A) + f_lm^ext)`
- (c) Both derived from a third source that contradicts both → audit stages 002+003 to identify upstream convention author

## Recommendation

direction: a
rationale: |
  The paper and source notes agree on the positive RHS source convention in the modal wall PDE (`paper/stages/stage_001.tex:158-162`, `notes/stages/moving_throat_pde_stage001_geometry_lift.md:357-361`).  With the scripts' current `-q*source_total` coupling, the Euler equation is the opposite-sign source equation, so the scripts should flip the source coupling and target sign.

## Apply: done
files_changed:
  - scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:189-191: flipped the modal source coupling and target to the positive RHS convention.
  - mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:188-191: flipped the modal source coupling and target to the positive RHS convention.
sanity_check: "sympy + mathematica exit 0"
notes: Applied with Q2; stage 001 was run once after both edits.

---

## Q2 — Stage 001 F2: gauge-fix term sign (metric signature)

**Files to read for context:**
- Paper: `paper/stages/stage_001.tex` (Section "Localized Maxwell" / gauge-fix term)
- Notes: `notes/stages/moving_throat_pde_stage001_geometry_lift.md`
- SymPy: `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:211,225,226`
- Mathematica: `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:205,210,211`
- Full audit report: `redteam/reports/stage_001.md`

**The disagreement:**
The gauge-fix term sign in the linearized Maxwell sector disagrees between script and paper. The disagreement may be a metric-signature convention difference that the paper card never states.

**Options:**
- (a) Mostly-minus signature `(+,-,-,-,-)` → script is correct as-is; add a one-line metric-signature statement to the paper card or script docstring
- (b) Mostly-plus signature `(-,+,+,+,+)` → script's gauge-fix term sign is wrong; flip it at the 6 line numbers above
- (c) Non-standard convention (or `Z(w)` absorbs a non-trivial wall-direction metric factor) → flag for deeper review against the parent-theory definition document

## Recommendation

direction: b
rationale: |
  The project-wide parent convention is mostly-plus, with `eta_{MN}=diag(-1,+1,+1,+1,+1)` stated in `paper/frontmatter/03_notation_firewall.tex:88` and `paper/parts/part01_parent_geometry.tex:132`.  Under that convention the paper's `+(1/xi) partial^N(...)` gauge-fix term has the positive spatial-coordinate sign, while the stage 001 scripts assert the negative sign.

## Apply: done
files_changed:
  - scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:209-226: flipped the representative gauge-fix Lagrangian and Maxwell targets to the mostly-plus sign.
  - mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:205-211: flipped the representative gauge-fix Lagrangian and Maxwell targets to the mostly-plus sign.
sanity_check: "sympy + mathematica exit 0"
notes: Applied with Q1; stage 001 was run once after both edits.

---

## Q3 — Stage 006: gauge-driver placeholder in scripts not in paper

**Files to read for context:**
- Paper: `paper/stages/stage_006.tex`
- Upstream paper: `paper/stages/stage_005.tex` (defines `L_mix`)
- SymPy: `scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:7-12,57,66-69,135,142-144,150-152,166-167,182,193-195`
- Mathematica: `mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl:104,117,119,121-122`
- Full audit report: `redteam/reports/stage_006.md`

**The disagreement:**
The paper card and the upstream stage 005 name only `L_mix` (transverse leakage) as the open-system inhomogeneity of the projected Maxwell laws. The SymPy and Mathematica scripts additionally carry symbolic `Gauge_μ` / `gauge[nu]` placeholder terms in both Gauss and Ampère equations (advertised in the SymPy docstring and Section 4 summary as "gauge-driver" contributions). These placeholders are never exercised by any assertion — they cancel trivially between LHS and RHS.

**Options:**
- (a) Paper correct — drop the `Gauge_μ` placeholders from the scripts (auditor's recommended direction: never-exercised placeholder, minimum risk)
- (b) Scripts correct — paper card adds a gauge-driver derivation paragraph
- (c) Future hook — drop placeholders now; reintroduce when a concrete gauge-fixing realization exists (deferred (a))

## Recommendation

direction: c
rationale: |
  The symbolic `Gauge_mu` placeholders in stage 006 are not exercised by the assertions and are not part of the stage 005/006 paper claims, which name only the transverse leakage term (`paper/stages/stage_006.tex:38-43`).  The genuine weighted gauge-fixing channel is handled concretely in stage 008 (`paper/stages/stage_008.tex:13-22`), so drop the placeholders here and reintroduce only through a concrete `H(w)` realization.

## Apply: done
files_changed:
  - scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:8-11: removed the gauge-driver claim from the scope summary.
  - scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:57-65: removed the symbolic Gauge0..Gauge3 placeholders.
  - scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:130-163: rewrote Gauss/Ampere targets and summaries as leakage-only inhomogeneous laws.
  - mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl:100-124: removed the gauge vector placeholders from the projected inhomogeneous rearrangement.
sanity_check: "sympy + mathematica exit 0"
notes: none

---

## Q4 — Stage 007: paper exports `mu_eff^proj` AND `1/xi_eff^proj`; script only verifies `mu`

**Files to read for context:**
- Paper: `paper/stages/stage_007.tex` (equations (1) and (2), Output paragraph)
- Part appendix: `paper/appendices/stage_appendix_part01.tex:36`
- SymPy: `scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py` (docstring explicitly disclaims gauge-driver channel)
- Mathematica: `mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl`
- Full audit report: `redteam/reports/stage_007.md`

**The disagreement:**
The paper card's equation (2) exports two observer-dependent projected coefficients (`mu_eff^proj = mu0 · I_WS/I_WZ` AND `1/xi_eff^proj = I_WH/(xi · I_WZ)`), and equation (1) defines three projected integrals (`I_WZ`, `I_WH`, `I_WS`). Both scripts verify only the `mu_eff^proj` half. The script docstring openly disclaims the gauge-driver channel; the paper's `\stagefield{Output}` still promises both.

**Options:**
- (a) Paper is correct → extend scripts: add an `H(w)` profile (e.g., `H(w) = exp(-w^2/rho^2)` with new positive parameter `rho`), compute `I_WH = ∫ W H dw` for the smooth Gaussian `W`, matched-observer `W_match = Z/Z_int`, and regulated `W_eps`; assert closed forms in both engines; assert `xi_eff^proj = xi · I_WZ / I_WH` against an analogous reduction-first reference value; drop docstring lines that disclaim gauge-driver
- (b) Script is correct → trim paper card so equation (1) only defines `I_WZ` and `I_WS`, equation (2) only carries the `mu_eff^proj` line, Output paragraph references only the `mu` channel; possibly re-word the part appendix row
- (c) Both are sourced from a deeper note (the original em_projected source notes were never committed); route the gauge-driver channel to a later stage with a forward pointer

## Recommendation

direction: a
rationale: |
  Stage 007 explicitly exports both observer-dependent coefficients, including `1/xi_eff^proj = I_WH/(xi I_WZ)` (`paper/stages/stage_007.tex:16-28`, `paper/stages/stage_007.tex:37-40`).  The gauge-weighted formula is mathematically part of the same zero-mode projection law, and stage 008 already confirms this channel is real rather than decorative, so the stage 007 scripts should verify the missing `I_WH`/`xi_eff` half.

## Apply: done
files_changed:
  - scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py:7-23: replaced the gauge-driver disclaimer with the zero-mode H(w) channel and xi_eff formula.
  - scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py:48-53: added H(w)=exp(-w^2/rho^2) and H_int.
  - scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py:84-105: added smooth I_WH and xi_eff checks.
  - scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py:142-174: added matched-observer I_WH, xi_eff, and reduction-reference checks.
  - scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py:179-205: added regulated sharp-observer I_WH and xi_eff checks.
  - mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl:8-18: added rho, xi, Hprofile, and Hint.
  - mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl:36-82: added Gaussian H area, smooth I_WH, and smooth xi_eff checks.
  - mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl:119-168: added matched-observer I_WH, xi_eff, and H=Z alignment checks.
  - mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl:170-243: added regulated I_WH and sharp xi_eff limit checks.
  - paper/appendices/stage_appendix_part01.tex:36: updated the Stage 007 appendix row to include gauge-parameter mismatch.
sanity_check: "sympy + mathematica exit 0"
notes: Used a new positive rho parameter for H(w), while preserving the Stage 008 H=Z specialization as a checked subcase.

---

## Q5 — Stage 010 F1: `δu_n` vs `δP_n` (distinct objects)

**Files to read for context:**
- Paper: `paper/stages/stage_010.tex` (`\stagefield{Output}` block and eq:stage010-du2, eq:stage010-du4)
- SymPy: `scripts/moving_throat_pde_stage010_*_sympy_audit.py:60-67,71-82`
- Mathematica: `mathematica/moving_throat_pde_stage010_*_mathematica_audit.wl`
- Full audit report: `redteam/reports/stage_010.md`

**The disagreement:**
The paper card asserts three first-variation identities as the stage's `\stagefield{Output}`:
- `δu_2 = (D_0 z_2 - D_2 z_0) / D_0^2`
- `δu_4 = (D_0^2 z_4 - D_0(2 D_2 z_2 + D_4 z_0) + 2 D_2^2 z_0) / D_0^3`
- `δP_0 = (D_0 n_0 + N_0 z_0) / D_0^2`

The scripts verify `δP_0` exactly but do not verify `δu_2` or `δu_4`; instead they verify `δP_2` and `δP_4` (full numerator/denominator ratio variations including `n_n` shifts), which are distinct objects.

**Options:**
- (a) Paper is correct — `δu_2`/`δu_4` are the intended export; add explicit `assert_zero` checks for those identities to both engines (existing `δP_2`/`δP_4` blocks may stay as cross-checks or be removed)
- (b) Script is correct — `δP_2`/`δP_4` (full ratio variations) are the actual stage 010 exports; rewrite paper card `eq:stage010-du2` and `eq:stage010-du4` to display the explicit `δP_2`/`δP_4` closed forms; rename symbol `u_n → P_n`; update downstream paper text that references those eqrefs
- (c) Both `u_n` (pure D-side) and `P_n` (full ratio) are physically distinct exports — paper gets a `δP_n` paragraph, scripts get `δu_n` assertions (both sides grow)

## Recommendation

direction: c
rationale: |
  The `u_n` and `P_n` quantities are distinct and both are legitimate: `u_2,u_4` are pure denominator-inversion coefficients, while the scripts' `P_2,P_4` are full numerator/denominator prefactor variations (`scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:50-83`).  Since the paper declares shifts for both `Z_n` and `N_n` but only displays `delta u_2`, `delta u_4`, and `delta P_0` (`paper/stages/stage_010.tex:17-37`), the paper should add the `delta P_n` exports and the scripts should add the missing `delta u_n` assertions.

## Apply: done
files_changed:
  - paper/stages/stage_010.tex:39-58: added the distinct full-prefactor delta P2 and delta P4 exports with labels eq:stage010-dP2 and eq:stage010-dP4.
  - scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:50-65: added denominator-only u2/u4 definitions and assertions.
  - mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl:31-57: added denominator-only u2/u4 definitions and assertions.
sanity_check: "sympy + mathematica exit 0"
notes: Applied with Q6; stage 010 was run once after both edits.

---

## Q6 — Stage 010 F2: 7 script clusters with no paper anchor

**Files to read for context:**
- Paper: `paper/stages/stage_010.tex`
- SymPy: `scripts/moving_throat_pde_stage010_*_sympy_audit.py:2,85-169`
- Mathematica: `mathematica/moving_throat_pde_stage010_*_mathematica_audit.wl:65-220`
- Possibly relevant: stages 011 (P_2 bridge), 012 (primitive bridge), 013–014 (mouth-Taylor)
- Full audit report: `redteam/reports/stage_010.md`

**The disagreement:**
The scripts certify seven content clusters that the paper card never mentions and that have no anchoring notes file in the repo:
- K-surface compatibility
- Transported-target compatibility
- Gaunt Y_{20} lane multipliers
- Weak-axisymmetric trace anomaly
- Primitive static Xi
- Two sign-flip mutation guards

(The SymPy docstring references `step_08_projected_maxwell_push_bundle_master_notes.md`, an EM-projected step file that was renumbered to stage 010 and has no surviving notes file in the repo.)

**Options:**
- (a) Paper card is incomplete — these clusters are load-bearing for downstream stages (likely 011, 012, 013–014); `stage_010.tex` needs new paragraphs describing each cluster with explicit equation labels; optionally enrich the part appendix Verification-output column
- (b) Scripts over-verify — these are scaffolding from an earlier draft; remove the script blocks at the line ranges above, leaving only the `δP_0` check (and any `δu_n` content from Q5) as the audit body
- (c) File the missing notes — author publishes `notes/em_projected/step_08_projected_maxwell_push_bundle_master_notes.md` (or equivalent under `notes/stages/`) to anchor the seven clusters as legitimate stage 010 derivations; the stage card optionally cross-references via `\StageFile`

## Recommendation

direction: a
rationale: |
  These clusters are load-bearing downstream rather than stale scaffolding: stage 011 publishes the compatibility bridge, stage 012 explicitly uses fixed and transported target compatibility, and stages 013-014 push the primitive/Xi and constant-prefactor transport into mouth-Taylor gates (`paper/stages/stage_011.tex:13-39`, `paper/stages/stage_012.tex:23-30`, `paper/stages/stage_013.tex:33-41`, `scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:18`).  Stage 010 is the master slot-transport stage, so its paper card should grow enough anchors for the script's verified transport clusters.

## Apply: done
files_changed:
  - paper/stages/stage_010.tex:60-94: added fixed-target and transported-target compatibility transport anchors.
  - paper/stages/stage_010.tex:96-112: added the Y20 lane multiplier and weak-axisymmetric trace/anomaly anchors.
  - paper/stages/stage_010.tex:114-135: added the primitive static Xi and mutation-guard anchors.
  - paper/stages/stage_010.tex:137-144: expanded the Output field to cite the new Stage 010 anchors.
  - paper/appendices/stage_appendix_part01.tex:42: updated the Stage 010 appendix row to cover the expanded paper-card scope.
  - scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:2: replaced the stale missing-notes docstring with a Stage 010 audit description.
  - scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:180-181: updated the printed audit summary to Stage 010 and the expanded u/P scope.
sanity_check: "sympy + mathematica exit 0"
notes: Renamed the local shifts label from eq:stage006-projected-shifts to eq:stage010-projected-shifts; rg found no external refs.

---

## Q7 — Stage 011: 5 script clusters with no paper anchor

**Files to read for context:**
- Paper: `paper/stages/stage_011.tex` (Output block declares only K-eliminated compatibility surface)
- SymPy: `scripts/moving_throat_pde_stage011_*_sympy_audit.py:132-155,168-176,179-211`
- Mathematica: `mathematica/moving_throat_pde_stage011_*_mathematica_audit.wl:97-153`
- Possibly relevant: stages 012 (primitive bridge), 022 (grouped normalization bridge), 023 (full grouped bundle)
- Full audit report: `redteam/reports/stage_011.md`

**The disagreement:**
The script verifies more than the stage card declares as Output. The card promises only:
1. K-eliminated compatibility surface `C = N_0/P_{0,target} - 3 S^2/T`
2. Its first variation `δC = n_0/P_{0,target} - 6 S z_2/T + 3 S^2 z_4/T^2`
3. The `z_0`-cancellation property

The script additionally verifies:
- Transported-target variant `compat_transport = D_{0,target} - 3 S^2/T` with its `z_0`-independent first variation
- Constant-prefactor branch factorization `P_2 = (N_2 - 2 D_2 N_0/D_0)/D_0` and `P_4 = (N_4 - N_4^*)/D_0`
- Real-Y20 grouped lane signature `(λ_{20}, λ_{21}, λ_{22}) = (1, 1/2, -1)` with `xbar = x^(0)` and `b = 3a`
- Static Xi1 prefactor slope `Xi1^(proj,static) = na/N_0 + za/D_0`
- Per-lane `u_2` slope `(D_0 z_{2a} - D_2 z_a)/D_0^2`

**Options:**
- (a) Paper correct as scoped — trim the script (remove the line ranges above); verify no downstream stage script imports or quotes these results from stage 011 before deleting (if so, move the block into that downstream stage)
- (b) Script is correct — extend the paper card: add Output block (or supplemental equations) enumerating the five extra clusters; optionally add per-stage notes
- (c) Split across stages — keep only the core K-eliminated compatibility surface in stage 011; move the Y20 / Xi1 / u_2 / constant-prefactor work into the stage(s) that actually consume them (likely 012 "primitive bridge", 022 "grouped normalization bridge", or 023 "full grouped bundle"); script-side relocation + paper-side anchoring at destination

## Recommendation

direction: c
rationale: |
  The stage 011 paper is correctly scoped to the K-eliminated compatibility surface and `z_0` cancellation (`paper/stages/stage_011.tex:13-39`), while the extra Y20, Xi1, `u_2`, and constant-prefactor material has natural anchors in later stages.  In particular, stages 022-023 publish the `u_n`, `P_n`, constant-prefactor, and anisotropy transport formulas, and stage 024 publishes the weak-axisymmetric `(1,1/2,-1)` / `b=3a` signature (`paper/stages/stage_023.tex:238-290`, `paper/stages/stage_024.tex:225-235`).

## Apply: done
files_changed:
  - scripts/moving_throat_pde_stage011_projected_maxwell_p2_bridge_sympy_audit.py:1-85: trimmed the audit to the one-pole/fixed-target K surfaces, K-eliminated compatibility variation, and z0 cancellation.
  - mathematica/moving_throat_pde_stage011_projected_maxwell_p2_bridge_mathematica_audit.wl:1-63: trimmed the audit to the one-pole/fixed-target K surfaces, K-eliminated compatibility variation, and z0 cancellation.
sanity_check: "sympy + mathematica exit 0"
notes: No paper move was needed; stage 012, stage 023, and stage 024 already verify/publish the removed destination material.

---

## Summary

summary: 7 answered, 0 skipped

## Apply notes

Overall summary: 7 applied, 0 revised, 0 blocked.

Touched files outside the original question file:line lists: `paper/appendices/stage_appendix_part01.tex` rows 007 and 010 were updated to match the paper-card scope changes.  No notes files, audit pipeline files, build artifacts, or `.claude` files were edited.

Concerns for Trevor: the old Stage 010 local label `eq:stage006-projected-shifts` is now `eq:stage010-projected-shifts`; `rg` found no downstream references, but future citations should use the Stage 010 label.
