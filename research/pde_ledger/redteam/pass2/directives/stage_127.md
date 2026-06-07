---
unit_id: 127
batch: IV.4
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T23:02:03Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 127 (USER-RESOLVED: correct notes → script)

The single finding is a `paper_misalignment` (value_mismatch) that the user has
RESOLVED in favor of the SCRIPT side (the engines are the computed truth; the notes
carry a trailing-digit transcription typo). This follow-up directive EXPLICITLY
AUTHORIZES the notes-side edit below. Apply ONLY this edit.

This is a NOTES-ONLY edit. Do NOT run any script (no `python3`, no `math -script`)
— no script changes. Do NOT touch `paper/`, any `.py`, or any `.wl`. After applying,
append an `## Applied: F1` block with `files_changed`, `summary`, and `deviation`
(or "none").

## F1 — paper_misalignment (value_mismatch) — USER-AUTHORIZED notes edit

**Subtype:** value_mismatch (resolved: script is correct)

**File / line:** `notes/stages/moving_throat_pde_stage127_penetration_families.md:78`

**Before:**
```
x_*^{\exp}\approx 0.662765402623160.
```

**After:**
```
x_*^{\exp}\approx 0.662765402623161.
```

**Rationale (for the Applied block, do not add to the file):** both engines
independently compute `x*_exp = 0.6627654026231614025…`; rounded to 15 dp this is
`0.662765402623161`. The notes' 15th decimal `0` is a transcription typo. The
companion slab depth in the same box (`0.797839360904564`) already matches the
engines to 15 dp, fixing the intended precision at 15 dp. Scripts are unchanged
(already correct).

**Verification:** the orchestrator confirms line 78 now reads
`x_*^{\exp}\approx 0.662765402623161.` and that no `.py`/`.wl`/`paper/` file changed.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage127_penetration_families.md`
- summary: Corrected the exponential penetration threshold's final decimal to match the independently computed 15 dp value.
- deviation: none
