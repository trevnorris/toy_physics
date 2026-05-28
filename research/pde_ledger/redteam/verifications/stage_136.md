---
unit_id: 136
batch: IV.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 136

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved

**What changed:**

Stage 136 was held for user resolution per the directive's "Resolve before fix_loop" block (no Codex auto-apply). The user resolved it through Clusters A and B of `redteam/resolutions/batch_IV4_paper_alignment.md`:

- Cluster A — Notes H1 mass-fix: `notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:2` now reads `# Moving-Throat PDE — Stage 136: Mouth-Layer Fixed-Point Status After the Coupled Solve`. Previously this line was `Stage 238` (matching the `actual + 102` offset for the affected notes group 127-136). The H1 now matches the filename, the paper card, and the manifest index.
- Cluster B — Status-only carry-forward re-attribution: `notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:6` now reads `After Stages 133–135, the mouth-source selection problem has narrowed again.` Previously this referenced "Stages 184–186" (downstream, impossible as upstream). The new citation 133-135 names the actual IV.4 upstream stages where the coupled fixed-point law (133), the Family-1 first-explicit reduction with `kappa_s=0, kappa_q=pi/2` (134), and `Σ_m^*` under the 4:−1 outlet-consistent closure (135) are derived — matching the resolution document's identification of where `(M_s, M_q)` reduction laws and `Σ_m^* ≈ 0.451485277739090` originate within batch IV.4.

**Assessment:**

Both edits land exactly where the directive's "Resolve before fix_loop" block flagged the inconsistency (H1 line and the body's upstream citation). Direction (a) from the directive was chosen: rewrite the notes to point to the actual IV.4 upstream stages. The carry-forward chain is now internally consistent — filename (136), H1 (Stage 136), and body citation (Stages 133-135 upstream of 136) all agree, and the cited stages do derive the load-bearing constants per the resolution document's batch-internal provenance table. The numerical constants on lines 30 and 40 (`M_s ≈ 1.50882951349316 − 0.658075937605429·M_q` and `Σ_m^* ≈ 0.451485277739090`) are unchanged and now have a coherent batch-internal anchor. No `paper/stages/stage_136.tex` edit was required for this resolution since the card's `\stagefield{Inputs}` (L9) describes content abstractly and the notes-level provenance fix suffices.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script exists for stage 136 (`is_status_only_candidate: True`). The paper card honestly self-declares `SymPy audit: none yet.` This is consistent with the original auditor's status-only carve-out.

**Mathematica:** exit=n/a. No Mathematica script exists for stage 136. Card self-declares `Mathematica audit: none yet.` Consistent.

**Output freshness:** No `.txt` outputs exist (status-only unit). N/A.

## Material-change assessment

`material_change`: false.

The edits are pure prose corrections to a notes file — H1 numeric label and an upstream stage-index citation. No derived results, no numerical constants, no formulas, and no script outputs changed. Downstream consumers (stages 146-153 per card L27) inherit the same numerical constants as before; only the provenance pointer they trace through this status entry is now coherent. No downstream re-audit is implied by this edit alone.

## Side observations (non-blocking)

- The paper card's `\stagefield{Inputs}` (L9 in `paper/stages/stage_136.tex`) still describes carry-forward content abstractly and does not name explicit upstream stage indices. The Cluster B resolution intentionally fixes only the notes (resolution doc L53: "The 132/136 paper cards themselves are clean"), so this is consistent with the chosen scope. No verifier action.
- The auditor's F1 framed the inconsistency as both an H1 mislabel and a body-citation error; both halves were resolved by Clusters A and B respectively, so the single finding is fully addressed even though the resolution split it across two cluster groups in the batch-level resolution doc.

## Verdict justification

The single paper_misalignment finding (F1) is resolved. The notes H1 on line 2 now reads "Stage 136" (Cluster A), and the upstream citation on line 6 now reads "Stages 133–135" (Cluster B) — exactly the lines and indices the directive's "Resolve before fix_loop" block called out as inconsistent. The carry-forward chain is internally coherent and traces to real upstream IV.4 derivations. No scripts are involved (status-only unit), so no exec logs to assess. Verdict: `verified`.
