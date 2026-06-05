---
unit_id: 024
batch: II.1
created_at: 2026-06-04T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-04T22:16:11-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 024

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F1 was a `paper_misalignment` that has now been RESOLVED (see the `## RESOLVED — F1` block below; Claude+Codex math-coverage resolution, non-conceptual, notes-only, published card unaffected). The RESOLVED block explicitly authorizes the single notes-side edit — apply exactly that, nothing more. Make no OTHER edit to paper.tex or notes/.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** value_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage024_overlap_isotropy.md:213` quote: "`int dOmega n_i n_j n_k n_l n_m n_n = (4 pi / 122) sum_pairings delta delta delta`"
- corroborating notes inconsistency: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage024_overlap_isotropy.md:217,221` quote: "`M^(20) = kappa_* diag(1, 1/2, 1/2, -1, -1)`" with "`kappa_* = sqrt(5) / (7 sqrt(pi))`"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:128` quote: "`return sp.simplify(4 * pi * s / 105)`"
- (Mathematica computes the moment by direct surface integration, `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:42-46`, and reproduces `kappa_* = sqrt5/(7 sqrt(pi))`.)

## Resolve before fix_loop

The stage-024 notes state the sixth moment of the unit sphere with prefactor `4 pi / 122` (notes:213), but the SymPy script uses `4 pi / 105` (py:128) and both engines reproduce the boxed result `M^(20) = (sqrt5/(7 sqrt pi)) diag(1,1/2,1/2,-1,-1)`. The number `122` is incompatible with both the script and the notes' OWN downstream `kappa_* = sqrt5/(7 sqrt pi)` (a `122` prefactor would rescale `kappa_*` by `105/122`). Which value is correct?

Possible directions (the user picks one):
- (a) Notes typo → change notes:213 `4 pi / 122` to `4 pi / 105`; no script change (the script and the notes' own `kappa_*` already imply `105`). This is the strongly indicated direction.
- (b) Script is wrong → would require the boxed `M^(20)`/`kappa_*` to change, contradicting the notes' line 221 and the `.tex` card line 164; not consistent with the rest of the stage, so unlikely.
- (c) Both derive from a third source that contradicts both → flag for deeper review.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## RESOLVED — F1 (Claude+Codex; non-conceptual; published card UNAFFECTED)

**Direction chosen: (a) notes-side typo.** Correct the notes prefactor; no script change.

**Rationale (math-determined, cross-corroborated):** the sixth moment of the unit 2-sphere is
`∫ dΩ n_i n_j n_k n_l n_m n_n = (4π/105) Σ_pairings δδδ` — the prefactor is `1/(d(d+2)(d+4)) = 1/(3·5·7) = 1/105`
for d=3. `122` is a transcription typo. It is contradicted by (i) the SymPy script `4*pi/105` (py:128), (ii) the
Mathematica direct surface integral (wl:42-46), and (iii) the notes' OWN downstream `κ_* = √5/(7√π)` (notes:221),
which only follows from `105` (a `122` prefactor would rescale κ_* by 105/122). The published `.tex` card is silent
on this prefactor (eq:app-stage024-fourth-moment gives 4π/15; the sixth-moment paragraph states no number), so the
PUBLISHED CARD IS UNAFFECTED — this is a notes-only edit. Non-conceptual (corrects a typo to the already-verified
value); resolved by Claude+Codex per the math-coverage delegation ([[feedback-claude-codex-resolve-math]]), not
escalated.

**Authorized edit (Codex applies — notes are Codex-applied per [[feedback-codex-is-fix-applier]]):**
- `notes/stages/moving_throat_pde_stage024_overlap_isotropy.md:213` — change `(4 pi / 122)` → `(4 pi / 105)`.
  Value-typo only; change nothing else on the line, and make NO other notes/paper edit. After editing, re-run both
  engines (there is no script change) to confirm exit 0 and that κ_* / M^(20) are unchanged.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage024_overlap_isotropy.md`
- summary: Corrected the notes-only sixth-moment prefactor typo from `4 pi / 122` to `4 pi / 105`.
- deviation: none

## F2 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt`

**Issue:** Both saved transcripts (mtime 2026-05-26) predate their scripts (mtime 2026-06-03), so they are stale. Their content still matches the current scripts (all checks read `= 0` / `PASS`), so this is informational. This is resolved by the orchestrator's standard independent re-run, not by a Codex source edit.

**Required change:**
No source edit. The orchestrator re-runs `python3 scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py` and `math -script mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl` to refresh the committed transcripts. If Codex is invoked for F1 after user resolution, re-run both scripts then and confirm they exit 0.

**Verification command:**
After re-run, the `.txt` mtimes exceed their script mtimes and content is unchanged (all PASS / `= 0`).

## Applied: F2

- files_changed:
  - `scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt`
- summary: Re-ran both audit engines under `timeout 600` and refreshed the saved transcripts.
- deviation: Mathematica transcript content changed only by updating the stale banner from `STAGE 007` to the current script's `STAGE 024`; all checks still pass.
