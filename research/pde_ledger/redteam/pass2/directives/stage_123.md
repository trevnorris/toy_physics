---
unit_id: 123
batch: IV.3
created_at: 2026-06-06T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: resolved
needs_user_resolution: false
---

## RESOLVED (2026-06-06)

**F1 (r_F1 surd) — direction (a):** user authorized matching the documents to the scripts; the two
paper surds `√(4107−117π²)`→`√(4107−100π²)` (`paper/appendices/stage_appendix_part04.tex:562` +
`paper/parts/part04_geometry_retarded_mouth.tex:576`) were corrected by an out-of-band Codex session
(019e9dea), Claude-reviewed. Stage-123 SCRIPTS correct and UNTOUCHED (clean). Logged PAPER_CLEANUP **P5-12**.
**F2 (card "Mathematica audit: none yet") — direction (b), DEFERRED** to the batched paper-card sweep
(status understatement; the retro-sweep `.wl` exists and passes). Logged in PAPER_CLEANUP **P5-12** deferred list.


# Codex directive — unit 123

Both findings are `paper_misalignment` items pending USER resolution. Codex applies
NOTHING on this unit until the user picks a direction. Do NOT edit `paper/`, `notes/`,
or the scripts to "fix" either item. The scripts are already correct and exit 0 on both
engines; no script change is authorized here.

## F1 — paper_misalignment (value_mismatch)

**Subtype:** value_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:562` quote: `\frac{\sqrt{4107-117\pi^2}}{10\pi}` (with `\approx1.77799353547498`)
- `/var/projects/toy_physics/research/pde_ledger/paper/parts/part04_geometry_retarded_mouth.tex:576` quote: `\frac{\sqrt{4107-117\pi^2}}{10\pi}` (with `\approx1.77799353547498`)

**Script side (correct):**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py:50-51` quote: `R = sp.Rational(37,20); rF = sp.sqrt(12*R**2/sp.pi**2 - 1)` → output `sqrt(4107 - 100*pi**2)`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl:97-98` quote: `R = 37/20; rF1 = reduceExact[Sqrt[12*R^2/Pi^2 - 1]]`

## Resolve before fix_loop

The imported `r_{F1}` definition prints the radicand as `4107-117π²` in the paper
(appendix:562 and part04:576) but `4107-100π²` in the stage-123 scripts and in the
origin notes (stage121:69, stage122:49,56,88, stage126:100, stage148:132). Algebra
forces `12·(37/20)²/π² − 1 = (4107−100π²)/(100π²)`, and only `100π²` reproduces the
paper's own stated `≈1.77799353547498` (the `117π²` form gives ≈1.7295). Which is
correct, `100π²` or `117π²`?

Possible directions (the user picks one):
- (a) Scripts/notes are correct (`100π²`) → change the two paper surds to
  `\sqrt{4107-100\pi^2}`, NO script change. (This is an upstream stage-121 anchor
  `eq:...-rF1`; the fix belongs to the stage-121 paper pass, not a stage-123 script edit.)
- (b) `117π²` is correct → then the scripts, the stated decimal, and five notes files
  are all wrong → flag for deeper review (very unlikely; contradicts the printed numeric).

The orchestrator will not invoke Codex on this unit until the user has chosen. Note:
this surd typo is already logged in `PAPER_CLEANUP_TRACKER.md` as a known paper-side
issue (the `100π²` script form is documented as algebraically correct); this finding
records that the appendix/part04 `117π²` surd itself does not yet appear corrected.

## F2 — paper_misalignment (notes_contradicts_script / documentation status)

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_123.tex:11` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl` exists (retro-sweep 2026-06-01) and passes all six checks; output transcript ends `Stage 123 Mathematica audit passed.`

## Resolve before fix_loop

The stage-123 card says "Mathematica audit: none yet", but a passing Mathematica `.wl`
exists and is committed. Update the card to cite the `.wl`?

Possible directions (the user picks one):
- (a) Card is stale → update `stage_123.tex:11` to name the `.wl`
  (`Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl}`),
  NO script change. (Paper-side edit — orchestrator/Codex per file-ownership rules, after user OK.)
- (b) Leave the card as-is for a later batched paper-card sweep.

The orchestrator will not invoke Codex on this unit until the user has chosen.
