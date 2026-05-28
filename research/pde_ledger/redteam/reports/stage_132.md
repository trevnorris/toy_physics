---
unit_id: 132
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
  notes_stage_files: [moving_throat_pde_stage132_mouth_boundary_layer_status.md]
  paper_appendix: present
---

# Audit unit 132 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_132.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the row including this stage at line 1298)
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 132 is presented as a "Mouth Boundary-Layer Status After Explicit Source-Law Extraction" card with `\claimstatus{\StatusExactClosure{} / \StatusOpen{}.}` and `\stagefield{Verification}{SymPy audit: none yet. Mathematica audit: none yet.}` — the card itself declares it carries no script-side verification. The bottom-line prose claim is the boxed line in the card:

> "The source-shape ambiguity is reduced to the selected bias \(\Pi_m\)."

The card lists three intent-level checks (positivity of mouth source; zero-flux and boundary-layer normalizations in the GNLS/localized-Maxwell reduction; Family-1 compensation point against the lower branch) but states no closed-form equations of its own; downstream use feeds Stages 133–145. The notes file (the same `stage132_mouth_boundary_layer_status` slug) supplies the concrete carry-forward content: an explicit source family \(\sigma_\Pi(z)=\Pi e^{-\Pi z/L}/[L(1-e^{-\Pi})]\), the Family-1 bias \(\mathfrak g_\Pi=2\Pi(2\Pi e^\Pi+\pi)/[(4\Pi^2+\pi^2)(e^\Pi-1)]\), monotone behavior \(\mathfrak g_\Pi:2/\pi\to1\) as \(\Pi:0^+\to+\infty\), and the canonical compensation point \(\Pi_*\approx 1.50882951349316\), with the parent threshold given as \(\partial_z\delta V_{\rm conf}|_{\rm m}-q_*\partial_z A_0|_{\rm m}=1.50882951349316\,\Theta_\sigma/L\).

## What the script claims to verify

No scripts exist for unit 132 (`scripts/moving_throat_pde_stage132_*` and `mathematica/moving_throat_pde_stage132_*` are absent). The paper card explicitly states `Verification: SymPy audit: none yet. Mathematica audit: none yet.` There is no script-side claim to test.

## Paper ↔ script cross-check

| Deliverable (paper / notes) | Script-side check | Status |
|---|---|---|
| Boxed prose claim: source-shape ambiguity reduced to selected bias \(\Pi_m\) | — | missing (status only) |
| Note 1: explicit source family \(\sigma_\Pi\) (carry-forward) | — | missing (status only) |
| Note 2: Family-1 bias \(\mathfrak g_\Pi\) closed form | — | missing (status only) |
| Note 3: monotonicity \(\mathfrak g_\Pi:2/\pi\to1\) | — | missing (status only) |
| Note 4: \(\Pi_*\approx 1.50882951349316\) | — | missing (status only) |
| Note 5: parent threshold identity | — | missing (status only) |
| Check (a) positivity of mouth source | — | missing (status only) |
| Check (b) zero-flux / boundary-layer normalizations | — | missing (status only) |
| Check (c) Family-1 compensation against lower branch | — | missing (status only) |

`paper_alignment` is set to `partial`: the paper card's own `Verification` field truthfully declares no verification yet, so there is no false claim of coverage. However, the notes file attributes the listed results to upstream "Stages 180–182" — and 180–182 lie numerically after this unit, so the asserted carry-forward source cannot logically precede stage 132 (see F1).

## Assertion inventory

No assertions exist (no scripts). Table omitted; inventory is empty by construction.

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script (here: notes contradict the unit's own ordering and identifier; there is no script, but the notes' carry-forward chain is the de facto "verification" surface and it is broken)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md:2` (title) and `:4` (carry-forward attribution)

**What's wrong:**

The notes file's H1 title reads `# Moving-Throat PDE — Stage 234: Mouth Boundary-Layer Status After Explicit Source-Law Extraction` while the filename, paper card label, and appendix input row all say Stage 132 (`paper/stages/stage_132.tex:1` → `\section[Stage 132]{Stage 132: Mouth Boundary-Layer Status After Explicit Source-Law Extraction}`; `paper/appendices/stage_appendix_part04.tex:1298` → `\input{stages/stage_132}`). Inside the body, the notes attribute the source-law results to upstream work:

> "After Stages 180–182, the mouth-source side is no longer an abstract profile problem." (`notes/.../stage132_mouth_boundary_layer_status.md:4`)

This says the explicit GNLS+localized-Maxwell source family, the Family-1 bias \(\mathfrak g_\Pi\), the monotonicity, the canonical compensation point \(\Pi_*\approx 1.50882951349316\), and the parent threshold identity are carried forward from Stages 180–182. But 180, 181, 182 are numerically *after* 132 in the unit ordering and cannot be upstream of a stage 132 status card. The notes title further suggests this file may have been authored for a stage numbered 234 (an even later unit) and then re-anchored to the 132 slug without revising the carry-forward chain. For a status-only unit the carry-forward chain *is* the verification, so a broken chain means none of the listed numeric results (\(\Pi_*\approx 1.50882951349316\) in particular) has an upstream anchor that genuinely precedes this stage.

**Why this matters:**

Per the status-only rule in the prompt, `missing_sympy`/`missing_mathematica` is *only* valid if the carry-forward chain does not actually support the claim. Here the carry-forward chain that the notes assert (180–182) cannot support stage 132 because it is in the wrong direction. The card's prose claim, "source-shape ambiguity is reduced to the selected bias \(\Pi_m\)," is then unsupported within the visible chain: neither this stage's scripts nor an actual upstream-preceding stage's scripts (visible to this auditor through the notes) verify the numeric \(\Pi_*\) or the bias \(\mathfrak g_\Pi\) closed form quoted in the notes. Downstream Stages 133–145 then propagate a numerical anchor with no validated source.

**Required change:**

This is a paper_misalignment requiring user resolution. See the `## Resolve before fix_loop` block in the directive. Possible directions include: (a) re-source the notes against the actual upstream stages where the GNLS/localized-Maxwell source-law extraction lives (the appendix tree under Part IV may have them at numerically lower indices — the auditor cannot read other unit files to confirm); (b) re-number this stage to its intended slot (e.g., 234) and update `paper/stages/`, `paper/appendices/stage_appendix_part04.tex`, and the notes filename consistently; or (c) declare the stage's prose claim as conjectural until the upstream content lands and revise the `\claimstatus` line accordingly.

**Verification:**

After user resolution: either the notes' "After Stages X–Y" line refers to upstream stage numbers strictly less than 132 (and the title H1 reads "Stage 132"), or the unit has been renumbered consistently across `paper/stages/stage_NNN.tex`, `paper/appendices/stage_appendix_part04.tex`, and `notes/stages/moving_throat_pde_stageNNN_*.md`.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists.

## Engine cross-check

Not applicable — neither engine present.

## Verdict justification

`verdict: findings` with one finding. The card itself is truthful about its lack of script-side verification (`SymPy audit: none yet. Mathematica audit: none yet.`), so under the status-only carve-out a missing-script finding is not directly warranted. However, the carry-forward chain — which is what the status-only rule asks me to check — is broken: the notes attribute the load-bearing numerical results to Stages 180–182, which lie *after* 132 and therefore cannot be upstream sources. That is a genuine `paper_misalignment` requiring user resolution (it could be a stage-numbering mismatch — the notes H1 says "Stage 234" — or a misattribution; either way the auditor cannot pick the direction). `stop_cold` is null because this is not downstream-mathematically-propagating in the script sense (no scripts depend on it yet); the downstream prose card (Stage 132's own `Downstream use` line) does cite Stages 133–145 as consumers, so the user should consider the propagation when choosing a direction.

## Self-test notes

Traps checked:
1. Status-only gating: confirmed no scripts on disk; confirmed the card itself declares `Verification: none yet` (so missing-script is not by itself a finding). Carry-forward chain check: the notes' explicit upstream attribution ("After Stages 180–182") points to higher-numbered stages, which is invalid under any sensible ordering — this is the basis for the single paper_misalignment finding.
2. Variable-independence / parity / trivial-case substitution: not applicable (no script, no assertions to mentally execute).
3. Path specification for the directive: the directive does not propose a new script (would be premature pending user resolution on the numbering/sourcing question); it only contains a `## Resolve before fix_loop` block. No file:line edits are prescribed for Codex.
4. Paper round-trip: the finding does not propose any paper-side edit; it surfaces a question for the user. No new misalignment is introduced.
