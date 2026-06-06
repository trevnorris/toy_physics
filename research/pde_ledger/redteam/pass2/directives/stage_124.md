---
unit_id: 124
batch: IV.3
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: deferred
needs_user_resolution: false
---

## RESOLVED (2026-06-06) — DEFERRED to the numbering pass

The `\stagefield{Purpose}` "Stage~141" (=124+17) self-label is the known +17 EM-realignment drift. It is
content-keyed paper-prose label drift → DEFERRED to the dedicated numbering pass (never offset-sweep), NOT
the red-team script loop. It is one instance of a DISCOVERED class — `\stagefield{Purpose}` card self-labels
the numbering reconciliation MISSED (its scan keyed on `\section`/`\label`, not `\stagefield{Purpose}`): also
117 ("Stage 134"), 119 ("Stage 136"), 120 ("Stage 137"). No script change (status-only unit). Logged PAPER_CLEANUP **P5-12**.


# Codex directive — unit 124

This unit's only finding is a paper-side `paper_misalignment` (stale stage-number self-label in a `.tex` card). Codex applies nothing — the orchestrator holds for user resolution. Do not edit `paper.tex`, `notes/`, or scripts. There are no script-side findings for this status-only unit (both engines are null by design; all carried-forward values reconcile against upstream stages 121/122/125 — see report).

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script (stale self-label; `+17` numbering-drift artifact)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_124.tex:7` quote: "\stagefield{Purpose}{Stage~141 is a core outlet realization ledger step.  Its audit target is the verification output quoted below.}"
- corroborating (genuinely different stage): `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_141.tex:1,7` — titled "Stage 141: Mouth-Gain Status Update", Purpose "Stage~141 is a coupled mouth fixed point and gain selection ledger step." So 141 is a distinct card; 124's "Stage~141" is a drifted self-label, not a cross-reference.
- corroborating (content key): `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:28` assigns "core realization, parent extraction, and geometric core selection" to stages 114–124 — i.e. the "core outlet realization" descriptor belongs to 124, not 141.

**Script side:**
- (none — status-only unit, no scripts)

## Resolve before fix_loop

The stage-124 card's Purpose field opens with "Stage~141" but the card's title, `\label{stage:124}`, and anchor are all 124, and the descriptor "core outlet realization ledger step" content-matches stage 124 (114–124 core block), not the genuinely-distinct stage 141 ("Mouth-Gain Status Update"). This is the known `+17` EM-realignment drift (124 + 17 = 141). Should `paper/stages/stage_124.tex:7` "Stage~141" be corrected to "Stage~124"?

Possible directions (the user picks one):
- (a) Yes, stale self-label → change `paper/stages/stage_124.tex:7` "Stage~141" → "Stage~124" (paper-side edit; Codex does not own paper/, so this is applied via the published-paper edit path, not the red-team script loop). No script change.
- (b) The "141" is intentional (a genuine pointer to stage 141) → leave as-is; but the surrounding text ("core outlet realization") then mismatches 141's actual subject, so this is unlikely.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No script-side fix is pending.
