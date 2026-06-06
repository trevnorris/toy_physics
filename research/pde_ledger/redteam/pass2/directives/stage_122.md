---
unit_id: 122
batch: IV.3
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: resolved
needs_user_resolution: false
---

## RESOLVED (2026-06-06) — direction (a)

The user authorized "get everywhere in the documents to match what the scripts derived." The r_F1
surd `√(4107−117π²)`→`√(4107−100π²)` was corrected in the published paper by an out-of-band Codex
session (019e9dea, scripts-only red-team loop forbids paper edits per codex.md), Claude-reviewed:
`paper/appendices/stage_appendix_part04.tex:562` + `paper/parts/part04_geometry_retarded_mouth.tex:576`.
The stage-122 SCRIPTS are correct and UNTOUCHED (clean). Logged PAPER_CLEANUP **P5-12**.


# Codex directive — unit 122

The single finding is a `paper_misalignment` (value_mismatch) requiring user resolution.
Codex applies NOTHING on this unit until the user picks a direction. Do not edit
`paper/`, `notes/`, or the scripts — the scripts are correct and must not be touched.

## F1 — paper_misalignment (value_mismatch)

**Subtype:** value_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:560-563` quote:
  `\mathfrak r_{F1} = \frac{\sqrt{4107-117\pi^2}}{10\pi} \approx 1.77799353547498`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:24-25` quote:
  `R = sp.Rational(37,20)` ; `rF = sp.sqrt(12*R**2/sp.pi**2 - 1)` → `√(4107-100π²)/(10π)`
- (also) `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage122_mouth_source_compensation_test.md:49` quote:
  `\mathfrak r_{F1}=\frac{\sqrt{4107-100\pi^2}}{10\pi}`

## Resolve before fix_loop

The appendix states the symbolic Family-1 branch value as `√(4107-117π²)/(10π)`, but the
scripts and the stage notes use `√(4107-100π²)/(10π)`. Arithmetic check: `12·(37/20)² = 41.07`,
so `r_F1² = 41.07/π² − 1 = (4107 − 100π²)/(100π²)` — the script/notes `100π²` form is correct.
The appendix's own quoted numeric `1.77799353547498` (and its downstream `g_-^F1≈0.758035`,
`g_+^F1≈2.79795` at L571/L573) reproduce the `100π²` form, not `117π²`, so `117` is an internal
typo. Which side is correct?

Possible directions (the user picks one):
- (a) Script/notes are correct (overwhelmingly likely) → fix appendix `stage_appendix_part04.tex:562`
  `4107-117\pi^2` → `4107-100\pi^2`. No script change. (Paper-side edit; per file-ownership policy
  this is a Codex-applied paper edit ONLY after the user authorizes it in a follow-up directive.)
- (b) Appendix `117` is intended → then the appendix's own numeric and the scripts/notes are all
  wrong; requires deeper review of the geometric input `L/a=37/20` and re-derivation. Unlikely.
- (c) Defer to broader numbering/paper-cleanup pass if this is part of a known appendix-prose drift.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
