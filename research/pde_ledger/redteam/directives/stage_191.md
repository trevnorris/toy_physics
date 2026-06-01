---
unit_id: 191
batch: V.3
created_at: 2026-06-01T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-01T11:33:48-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 191

This directive has two findings: F1 (`paper_misalignment`, stale stage-number label) — now RESOLVED via the settled canonical-stage-number convention (see `## RESOLVED — F1` below); apply the authorized banner relabel. F2 (`missing_mathematica`) — create the independent-route `.wl`. Do NOT edit paper.tex or notes/.

After editing/creating, RUN the affected scripts and iterate until they exit 0 (SymPy: `timeout 600 python3 .../stage191...sympy_audit.py`; Mathematica: `timeout 600 math -script .../stage191...mathematica_audit.wl`). A timeout (124) is a failure — reformulate, never raise the cap.

## F1 — paper_misalignment

**Subtype:** target_mismatch (cosmetic stage-label mismatch only; no verified identity, constant, or deliverable differs)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_191.tex:1` quote: `\section[Stage 191]{Stage 191: Minimal PDE data packet}`
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex:113` quote: `191 & Minimal PDE data packet & \StatusExactClosure{} & Defines Packet A, Packet B, \(\Delta_{\rm branch}\), \(\Delta_{\rm orbit}\), and the exact two-packet home-stretch theorem. \\`
- (context) notes `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage191_minimal_pde_data_packet.md:311` quote: `\boxed{\textbf{Theorem (Stage 242 home-stretch theorem).}}` and `:421` quote: `moving_throat_pde_stage242_minimal_pde_data_packet_sympy_audit.py`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py:65` quote: `banner("STAGE 174 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM")`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py:284` quote: `banner("STAGE 174 LEDGER")`

## Resolve before fix_loop

The script's two transcript banner strings label this unit "STAGE 174", and the source notes call the theorem the "Stage 242 home-stretch theorem" and reference `..._stage242_...` files, whereas the paper card and the part-05 appendix consistently number this unit **Stage 191**. This is a stale legacy-renumbering label only: all verified mathematics, constants, and deliverables in the script match the Stage 191 card exactly. The discrepancy is confined to (a) two `print` banner labels in the script and (b) legacy stage numbers in the notes prose. Which way should this be reconciled?

Possible directions (the user picks one):
- (a) Paper card (Stage 191) is canonical → in a follow-up directive, authorize Codex to change ONLY the two script banner string literals at lines 65 and 284 from "STAGE 174 ..." / "STAGE 174 LEDGER" to "STAGE 191 ..." / "STAGE 191 LEDGER" (string-only edit, no math change), then re-run SymPy to refresh the output header. Leave the notes prose alone (out of red-team scope; flag separately if the notes' 242/174 references should be reconciled).
- (b) The legacy numbers are intentionally preserved for cross-revision provenance → no edit; record that the banner mislabel is accepted and the report's F1 is waived.
- (c) Deeper review needed if the stage's identity (191 vs 242 vs 174) is genuinely ambiguous in the broader ledger.

## RESOLVED — F1 (settled canonical-stage-number convention)

Resolution = direction (a). Per the project's settled convention (stale script "STAGE NN" labels → the canonical internal stage number; Codex-CONCUR, established in the remediation batches), the canonical number for this unit is **191** (paper card + part-05 appendix). Authorized script-side edit — apply it:
- Relabel the two banner string literals only: line 65 `banner("STAGE 174 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM")` → `banner("STAGE 191 — MINIMAL PDE DATA PACKET AND THE EXACT HOME-STRETCH THEOREM")`; line 284 `banner("STAGE 174 LEDGER")` → `banner("STAGE 191 LEDGER")`. String-only; no math change.
- The notes-prose legacy numbers (242 home-stretch theorem; `..._stage242_...` filename references) are OUT OF red-team scope (notes are paper files) → to be logged to PAPER_CLEANUP_TRACKER, not edited here.

`needs_user_resolution` is cleared.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.txt`
- summary: Relabeled the two stale Stage 174 transcript banners to Stage 191 and regenerated the SymPy transcript.
- deviation: none

## F2 — missing_mathematica

**Issue:** Stage 191 is dual-engine-capable (exact Taylor-coefficient matches, log/exp inverse round-trips, projector algebra over a fixed `diag(1,2,2)` metric, guided substitution onto a defining surface) but has no Mathematica `.wl`. Under the dual-engine rule, an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage191_minimal_pde_data_packet_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the SymPy script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage191_minimal_pde_data_packet_sympy_audit.py` (the SymPy audit found these checks substantive and correct — verdict clean apart from the label). Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_191.tex` and the notes file are the source of truth (Packet A → Delta_branch, Packet B → Delta_orbit, and the two-packet home-stretch theorem: reduced closure ⟺ both packets vanish).
- Use Mathematica-NATIVE primitives (`Series`+`Coefficient`, `Solve`/`Reduce`, matrix products / projector algebra, `D[...]`) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port with the same variable names and step order. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (`expectZero`/`expectZeroMatrix`, `$Assumptions` positivity, `stripCE`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion must match the SymPy result.

**Anti-transliteration:** a `.wl` that merely re-types the SymPy closed forms and subtracts them is a transliteration and will be REJECTED at verification. Design a genuinely independent route.

**Verification command:** the verifier runs `redteam exec-mathematica 191`, confirms exit 0 with all PASS lines, and reviews that the `.wl` is a genuinely independent route whose conclusions agree with the SymPy engine.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage191_minimal_pde_data_packet_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage191_minimal_pde_data_packet_mathematica_audit.txt`
- summary: Added an independent Mathematica audit covering Taylor compilers, weighted grouped/projector algebra, Packet A to Delta_branch, Packet B to Delta_orbit, and the two-packet zero-set split.
- deviation: none
