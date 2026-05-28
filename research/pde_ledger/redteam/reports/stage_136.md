---
unit_id: 136
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md"]
  paper_appendix: present
---

# Audit unit 136 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_136.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_136}` line at L1306; no separate appendix row narrative for this stage)
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

The paper card (`paper/stages/stage_136.tex`) explicitly self-labels this as a status / ledger entry: `\stagefield{Verification}{SymPy audit: none yet. Mathematica audit: none yet.}` (L11). Its `\stagefield{Purpose}` (L7) names it "a coupled mouth fixed point and gain selection ledger step" whose "audit target is the verification output quoted below." The body block (L15–17) states: "The mouth-layer problem is reduced to a gain pair, or one gain under outlet consistency." The checks list (L21–25) requires (i) the gain pair `(M_s, M_q)` against outlet consistency, (ii) the self-matched susceptibility closure before using the one-scalar branch law, and (iii) that numerical fixed points are recorded as numerically located, not closed-form constants. `\claimstatus{\StatusExactClosure{} / \StatusOpen{}}` (L5) marks the stage as part-closed / part-open. The notes file deepens this with three carry-forward formulas: (a) the coupled fixed-point law `Pi = M_+ S(Pi, kappa_+) + M_- S(Pi, kappa_-)`, (b) Family-1 reduction `Pi = M_s + M_q S_q(Pi)` with `kappa_s = 0`, `kappa_q = pi/2`, (c) canonical compensation line `M_s ≈ 1.50882951349316 − 0.658075937605429 M_q`, and (d) Σ_m^* ≈ 0.451485277739090 under the 4:−1 weighting.

## What the script claims to verify

There is no SymPy script (`scripts/moving_throat_pde_stage136_coupled_mouth_status_sympy_audit.py` is absent) and no Mathematica script (`mathematica/moving_throat_pde_stage136_coupled_mouth_status_mathematica_audit.wl` is absent). No saved outputs exist. The manifest entry has `is_status_only_candidate: True`, and the paper card itself states no audit script exists yet, which is consistent with the missing files.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Coupled fixed-point law `Pi = M_+ S(Pi,kappa_+) + M_- S(Pi,kappa_-)` (notes item 1) | none in this unit | missing (status-only carry-forward expected) |
| Family-1 reduction `Pi = M_s + M_q S_q(Pi)`, `kappa_s=0`, `kappa_q=pi/2` (notes item 2) | none in this unit | missing (status-only carry-forward expected) |
| Canonical compensation line `M_s ≈ 1.50882951349316 − 0.658075937605429 M_q` (notes item 3) | none in this unit | missing (status-only carry-forward expected) |
| 4:−1 weighting closure `Pi = Sigma_m [4 − S_q(Pi)]`, `Sigma_m^* ≈ 0.451485277739090` (notes item 4) | none in this unit | missing (status-only carry-forward expected) |
| Card check 1: `(M_s, M_q)` against outlet consistency | none | missing |
| Card check 2: self-matched susceptibility closure precondition | none | missing |
| Card check 3: numerical FPs recorded as numerically located | bookkeeping; satisfied implicitly by notes phrasing | n/a |

Per `is_status_only_candidate: True`, the absence of a script in this unit is not by itself a finding so long as upstream scripts verify the carry-forward content. The notes attribute the narrowing to "Stages 184–186" (notes L7), but those indices are *downstream* of stage 136 in the linear ordering and so cannot be the carry-forward source. The card's `\stagefield{Inputs}` (L9) describes the carry-forward content abstractly ("shell/mixed core, the mouth source law, outlet consistency, core-to-mouth gain maps, and self-matched susceptibility closure") without naming explicit upstream stage indices. As a status-only auditor with a restricted reading list (paper card + notes for this unit only; no other-unit script reads permitted), I cannot positively confirm that an upstream unit's script verifies the canonical compensation line `M_s ≈ 1.50882951349316 − 0.658075937605429 M_q` or `Sigma_m^* ≈ 0.451485277739090`. This is the substance of finding F1.

## Assertion inventory

No scripts to inventory. No assertions exist.

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script (closest fit; the script side is empty, but the notes assert *upstream provenance* the card does not record and that is internally inconsistent on the face of it)

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:7`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:1` (header mislabel)
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_136.tex:9` (Inputs do not name upstream stages)

**What's wrong:**

The notes file's header (L1) reads `# Moving-Throat PDE — Stage 238: Mouth-Layer Fixed-Point Status After the Coupled Solve`, but the filename, the audit unit index, and the paper card all label this as stage 136. The notes body (L7) says "After Stages 184–186, the mouth-source selection problem has narrowed again." Stages 184–186 are numerically downstream of stage 136 in the part IV ordering, so they cannot be the source of carry-forward content for a stage 136 ledger entry; either the stage numbers in the notes are wrong, or this file was repurposed without re-numbering its provenance references, or the file was originally drafted for a different stage (e.g., 238) and was renamed without rewriting the body.

The paper card's `\stagefield{Inputs}` (L9) describes the imported content only abstractly ("shell/mixed core, the mouth source law, outlet consistency, core-to-mouth gain maps, and self-matched susceptibility closure") and does not name a specific upstream stage that derived `M_s ≈ 1.50882951349316 − 0.658075937605429 M_q` or `Sigma_m^* ≈ 0.451485277739090`. With the notes' explicit provenance pointers (184–186, 238) all internally inconsistent, the chain of evidence supporting the load-bearing numerical constants in the notes is not anchored anywhere a reader (or a downstream auditor) can verify.

**Why this matters:**

Status-only units are exempt from carrying their own script *only when* their carry-forward chain is unambiguous. Here the carry-forward chain is internally contradictory (header says 238, body says "After Stages 184–186", filename says 136). Downstream stages 146–153 plus the finite mouth-profile correction (card L27) inherit the compensation line and Σ_m^* value through this status entry; if no upstream stage actually derives them, the downstream chain rests on unverified numerics.

**Required change:**

Resolution belongs to the user, not Codex — this is a paper-side / notes-side correction, and the user must pick the direction. See `## Resolve before fix_loop` block in the directive.

**Verification:**

After user picks a direction: (a) if the notes were mislabeled, the H1 line, the "Stages 184–186" reference, and possibly an explicit `Inputs` upgrade in the card all settle on consistent upstream indices; or (b) if this is a genuine status carry-forward whose upstream derivation lives in a different earlier stage, that stage's index is named explicitly in the card's `\stagefield{Inputs}` and the notes header/body are realigned to stage 136. Verifier confirms by re-reading the notes header, the notes L7 sentence, and the card `Inputs` line and finding them consistent.

## Independent-derivation check (Mathematica)

Not applicable — no Mathematica script exists.

## Engine cross-check

Not applicable — no scripts.

## Verdict justification

Per the status-only carve-out, the absence of `.py` / `.wl` scripts is not by itself a finding; the paper card honestly self-declares "SymPy audit: none yet. Mathematica audit: none yet." However, the notes file that is supposed to anchor the carry-forward content contains internally inconsistent provenance (header L1 names "Stage 238", body L7 cites "Stages 184–186" as upstream of a stage labeled 136), and the paper card's `Inputs` line does not name a specific upstream stage that derived the load-bearing numerical constants (`M_s ≈ 1.50882951349316 − 0.658075937605429 M_q`, `Sigma_m^* ≈ 0.451485277739090`). That is a paper_misalignment finding requiring user resolution before Codex can do anything mechanical with it. Verdict: `findings`, `stop_cold: null`, `paper_alignment: partial`.

## Self-test notes

Variable independence and parity traps are not applicable (no scripts). I confirmed the load-bearing constants in the notes are exact numerical literals (1.50882951349316, 0.658075937605429, 0.451485277739090), which a future script would need to derive or cite an upstream derivation for; the paper card permits these as "numerically located, not closed-form constants" (card L24), so the script-side bar is "numerically locate the same fixed point from the same coupled law", not "derive the literal in closed form." Path-specification is not applicable; no missing-script directive is being written (only a `Resolve before fix_loop`).
