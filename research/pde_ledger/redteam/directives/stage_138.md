---
unit_id: 138
batch: IV.4
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 138

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_138.tex:1` quote: `\section[Stage 138]{Stage 138: Normalized Mouth-Gain Family and Compensation Ratio}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_138.tex:2` quote: `\label{stage:138}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage138_normalized_mouth_gain_family.md:1` quote: `# Moving-Throat PDE — Stage 138: Normalized Mouth-Gain Family and Compensation Ratio`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage138_normalized_mouth_gain_family_mathematica_audit.wl:26` quote: `banner["STAGE 121 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO"];`

## Resolve before fix_loop

The Mathematica audit script's banner identifies the audit as `STAGE 121`, but the paper card, notes file, and filename all identify the audit unit as `Stage 138`. The math content of the script is correct for Stage 138; only the banner string is wrong. Which direction should the orchestrator take?

Possible directions (the user picks one):
- (a) Banner is a typo — update Mathematica L26 to `banner["STAGE 138 — NORMALIZED MOUTH-GAIN FAMILY AND COMPENSATION RATIO"];` and re-run `redteam exec-mathematica 138` to refresh the transcript. (Recommended; the assertions, filename, and paper card all already say 138.)
- (b) Banner is canonical and the file lives under the wrong unit ID — flag for deeper review of how this stage maps to the unit-138 slot in the ledger.
- (c) The Mathematica file was copy-pasted from an unrelated Stage 121 audit and never re-derived for Stage 138 — flag for full re-derivation rather than a banner patch.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
