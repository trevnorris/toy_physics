---
batch: IV.6
created_at: 2026-05-27
clusters_resolved: 3
mechanical_groups: 6
status: directions_set
---

# Red-team batch IV.6 — paper-alignment resolution

13 audits landed (range 151-163). 19 findings raised across 10 dirty stages; 3 stages clean (153 status-only, 159 and 162 substantive-clean). Three user-gated clusters resolved on 2026-05-27 as `(Recommended)`; six mechanical groups applied without a gate.

## Cluster A — Mass renumbering (mechanical, 20 edits)

**Pattern:** Two uniform offsets, same shape as IV.4/IV.5 Cluster A:

- **Banners in both engines** of 9 stages have `STAGE N` where N = actual − 17:
  154→137, 155→138, 156→139, 158→141, 159→142, 160→143, 161→144, 162→145, 163→146.
- **Notes H1** in 2 files have `# Moving-Throat PDE — Stage M:` where M = actual + 85:
  159→244, 160→245.

Stages with no banner offset issue:
- 151, 152 use a title-only banner convention (no STAGE prefix).
- 157 banner already correctly says "STAGE 157".

Stages with no notes-H1 issue: 151, 152, 153, 154, 155, 156, 157, 158, 161, 162, 163 (all already correct).

Canonical numbering elsewhere (paper card filenames, appendix rows, MANIFEST, file paths) is 151-163.

**Direction:** mass-fix all 20 in place via `/tmp/iv6_mass_renumber.py`.

- 18 banner edits across 9 stages × 2 engines: `banner("STAGE <N-17> — …")` / `banner["STAGE <N-17> — …"]` → STAGE <N>.
- 2 notes H1 edits: `# Moving-Throat PDE — Stage <N+85>:` → `# Moving-Throat PDE — Stage <N>:`.

No script content other than the banner string changes; no notes content other than the H1 line.

## Cluster B — Body-text forward-stage citation re-attribution (53 edits across 13 notes files)

**Pattern:** Notes contain inline references to stages 188-250 — pre-renumber IDs. Three offsets surfaced, disambiguated by content cross-check:

- **−51 offset** for 188-199 range (IV.4 / IV.5 references — matches IV.5's own body offset):
  - 188-189 → 137-138 (mouth-gain formulas — appears in stage 158, 159 notes)
  - 197-198 → 146-147 (positive deformation expansion + first-order rigidity kernel — appears in stage 151)

- **−85 offset** for 239-248 range (IV.6 internal cross-references — NEW offset, specific to IV.6):
  - 239 → 154 (coevolving_core_mouth — appears in stages 155, 158, 163)
  - 240 → 155 (frozen_traction_fixedpoint — stage 157)
  - 241 → 156 (renormalized_canonical_branch — stages 157, 158, 159)
  - 242 → 157 (core_mouth_coevolution_status — stages 154, 158)
  - 243 → 158 (linear_defect_transport — stage 159)
  - 244 → 159 (hybrid_outlet_projection — stages 159, 160, 161, 162, 163)
  - 245 → 160 (bare_mixed_port_slippage — stages 160, 161)
  - 246 → 161 (dn_similarity_slippage — stages 161, 162)
  - 247 → 162 (parent_compensation_rigidity — stages 162, 163)
  - 248 → 163 (off_family_normal_coordinate — stage 163)

- **−102 offset** for 221 and 249-250 (IV.3 + IV.5 references — matches IV.5 body offset):
  - 221 → 119 (parent_balance_family — stage 162; content confirms vs. stage 119's "compensated throat-core condition in normalized parent ratios")
  - 249 → 147 (first_order_rigidity_kernel — stages 151, 153)
  - 250 → 148 (representative_positive_families — stage 152)

**Direction:** re-attribute all 53 inline citations via `/tmp/iv6_reattribute.py`.

The new −85 offset is content-disambiguated: any "Stage NNN" in 239-248 that could ambiguously map to either an IV.5 stage (−102) or an IV.6 stage (−85) is resolved by reading the cited content. E.g., stage 158's "Stage 241, the exact co-evolving compensated Family-1 point" matches stage 156 (Renormalized Canonical Branch Under Full Core–Mouth Co-Evolution), not stage 139 (family1_actual_mouth_gains) — so −85 applies.

## Cluster C — Stage 158 paper-card Checks downgrade

**Pattern:** Stage 158's `\stagefield{Checks}` lists three items but only item 1 ("deviations about the renormalized canonical point") is verified by stage 158's own scripts. Items 2 and 3 describe the broader transport program:

- Item 2 (even-preservation constraints) is verified downstream in stage 159 (`hybrid_outlet_projection`'s canonical-even gate).
- Item 3 (tangent motion δ⊥ = 0) is verified downstream in stages 162/163 (`parent_compensation_rigidity` + `off_family_normal_coordinate`).

**Direction:** rewrite items 2 and 3 as forward-carry citations:

```latex
\item Even-preservation of the canonical gate is verified downstream: see \ref{stage:159}.
\item Tangent motion on the parent compensation family giving \(\delta_\perp=0\) is verified downstream: see \ref{stage:162} and \ref{stage:163}.
```

**Distinction from IV.4/IV.5 Cluster C:** Both IV.4 stage 134 and IV.5 stage 144 downgraded items as **upstream** carry-forward citations (`\ref{stage:135}`, `\ref{stage:140}`, etc.). IV.6 stage 158 is the **first forward (downstream) carry-forward** because the program-level checks land downstream of the linear-transport step.

No script changes for stage 158 F1. The orchestrator edits the paper card directly per the user-chosen direction (Codex never edits paper.tex).

## Mechanical groups M1-M6 (applied without gate)

### M1 — `.wl` banner mass-rename (18 sites)

Same as Cluster A. Applied via `/tmp/iv6_mass_renumber.py`.

### M2 — `.py` banner mass-rename (none for IV.6)

Only stages 151, 152, 157 have `.py` banners and none had the −17 offset issue. Empty group for IV.6.

### M3 — Notes H1 mass-rename (2 sites)

Stage 159 and 160 H1 → `# Moving-Throat PDE — Stage 159:` / `# Moving-Throat PDE — Stage 160:`.

### M4 — Directive substance per stage (18 script-side findings)

Each non-paper_misalignment finding's "Required change" block applied directly by the orchestrator (Codex bypassed per [[feedback-codex-iterates-until-clean]] / III.4 stall lesson). 18 substantive edits across SymPy and Mathematica scripts. Three rework-loop catches required orchestrator-side corrections (stage 151 mpmath rewrite, stage 154 single-eps Series, stage 161 `depsg_branch` typo fix).

### M5 — Rework loop reruns

After substantive edits, all 12 SymPy + 12 Mathematica scripts re-run. Initial sweep had 2 failures (stage 161 SymPy + stage 154 Mathematica); both fixed and re-run clean.

### M6 — Stale output refresh

All output files refreshed by direct `python3` and `math -script` calls (not via `$RT exec-*` per [[feedback-no-parallel-exec-sympy]]). mtimes are now newer than corresponding script mtimes for all 12 dual-engine stages.

## Closing

19 findings resolved (0 blocked_legitimate). 13 stages verified. Material change in 1 stage (151, mpmath rewrite — downstream impact zero since the script verifies algebraic identities, not numeric carry-forwards). Twelfth consecutive zero-redirection batch.
