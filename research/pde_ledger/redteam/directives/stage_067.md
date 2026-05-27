---
unit_id: 067
batch: III.3
created_at: 2026-05-26T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T19:10:00Z
verification_status: pending
needs_user_resolution: false
orchestrator_applied_directly:
  reason: "Pure banner relabel with unambiguous direction; auditor report stated direction (a) does not require user resolution. Orchestrator applied directly without re-invoking codex."
---

# Codex directive — unit 067 (second pass)

The first-pass directive (this file's previous content) recorded that F1-F4 from the first pass were applied at 2026-05-22T19:54:59. The current scripts and outputs reflect those fixes. The second-pass audit identified one additional finding, which is a `paper_misalignment` (subtype: stale label) and therefore routed to the user, not to Codex.

Codex must NOT edit any file under `paper/` or `notes/`. Codex must NOT edit the scripts to resolve the finding below until the user has chosen a direction in a follow-up directive.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script (script banners/docstring still use the pre-renumbering "Stage 50/050" label; paper card, filenames, and the mathematica closing print all use "067")

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_067.tex:1` quote: `\section[Stage 067]{Stage 067: Exact Sech--Gaussian Coherence Resonance Benchmark}`
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex:252` quote: `\input{stages/stage_067}`
- (notes-side, out of red-team scope) `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md:2` quote: `# Moving-Throat PDE — Stage 50: Exact Sech–Gaussian Coherence Resonance Benchmark`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:4` quote: `moving_throat_pde_stage50_sech_gaussian_sympy_audit.py`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:53` quote: `banner("STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl:38` quote: `banner["STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"];`
- (for comparison, the same .wl already uses "067" at line 174 quote: `Print["Stage 067 Mathematica audit passed."];` — so the file is internally inconsistent.)

## Resolve before fix_loop

The script docstring and printed banners still carry the pre-renumbering label "Stage 50" / "STAGE 050", while the paper card, file paths, appendix `\input`, and the mathematica closing print all use "Stage 067". The notes file's title is also still "Stage 50". Which label is canonical?

Possible directions (the user picks one):
- (a) "Stage 067" is canonical (default; paper card and filenames agree) → update the script labels:
   1. `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:4`: change `moving_throat_pde_stage50_sech_gaussian_sympy_audit.py` to `moving_throat_pde_stage067_sech_gaussian_sympy_audit.py`.
   2. `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:53`: change the banner argument from `"STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK"` to `"STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"`.
   3. `mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl:38`: change the banner argument from `"STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"` to `"STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"`.
   The notes file's "Stage 50" header would need a separate, user-initiated edit (out of red-team scope).
- (b) "Stage 50" is the historical label and should be preserved as a cross-reference (e.g., because downstream references still use it) → update the paper card's section heading and the appendix `\input` to use "stage_050" instead. This is a substantially larger change touching paper/, which the red-team will not initiate without explicit user direction; it is unlikely to be the intended direction since the filenames are already `stage_067`.
- (c) Keep both — append a parenthetical "(legacy: Stage 50)" or "(formerly Stage 50)" to the script banners and the paper card heading so both labels coexist clearly.

The orchestrator will halt and not invoke Codex on this unit until the user picks (a), (b), or (c) (or another direction). If the user chooses (a), a follow-up directive can be auto-generated using the exact three before/after edits above; no Codex re-audit needed for math, only `redteam exec-sympy 067` and `redteam exec-mathematica 067` to confirm both scripts still exit 0 with the new banner text.

**Verification (post-resolution):**
- For direction (a): the sympy transcript's `STAGE 50 — ...` banner line becomes `STAGE 067 — ...`; the mathematica transcript's `STAGE 050 — ...` banner line becomes `STAGE 067 — ...`; both scripts still exit 0; no assertion text or numeric output changes.
