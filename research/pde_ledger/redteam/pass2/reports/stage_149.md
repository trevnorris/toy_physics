---
unit_id: 149
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage149_mouth_rigidity_status.md]
  paper_appendix: present
---

# Audit unit 149 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_149.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage149_mouth_rigidity_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only an `\input{stages/stage_149}` aggregator line at :1332 — no separate prose row)
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 149 is explicitly a status/ledger card in the "Finite mouth-profile corrections" block. The card's `\stagefield{Purpose}` (`stage_149.tex:7`) states "Stage~149 is a finite mouth-profile corrections ledger step. Its audit target is the verification output quoted below," and the quoted output (`stage_149.tex:15-17`) is the qualitative claim "Mouth corrections are controlled by a finite rigidity kernel, not by branch ambiguity." The card carries no numeric `\stagefield{Output}` and `\stagefield{Verification}` (`stage_149.tex:11`) reads verbatim: "SymPy audit: none yet. Mathematica audit: none yet." The notes (`moving_throat_pde_stage149_mouth_rigidity_status.md`) elaborate the perturbative-rigidity summary: any positive normalized mouth deformation near the canonical exponential branch is written `Sigma_eps = (1-eps)Sigma_* + eps*varsigma` (:8-16), the first-order shift enters only through two source moments `bar_g_varsigma`, `bar_S_varsigma` (:19-23), and the traction shift `delta_T_m = eps[A_T(bar_g - g_*) + B_T(bar_S - S_*)]` with `A_T ≈ -4.27263956256927`, `B_T ≈ 0.134875005736706` (:42-46) and the overlap-dominance ratio `|A_T|/B_T ≈ 31.68` (:51-52). The notes are explicit (`:54-60, :73-74`) that these values are imported "at the same first-order deformation level used in Stages 146-148," i.e. carried forward, not derived at stage 149. The card itself declares the stage "a derivation ledger entry, not an unconditional actual-branch theorem" (`stage_149.tex:27`).

## What the script claims to verify

Nothing — there is no SymPy script, no Mathematica script, and no saved output for unit 149. The MANIFEST entry (`redteam/pass2/MANIFEST.yaml:5059-5083`) confirms `is_status_only_candidate: true`, `is_checkpoint: false`, and all four file paths (`sympy`, `mathematica`, `sympy_output`, `mathematica_output`) `path: null` / `exists: false`. A broad filename search (`find scripts mathematica -iname '*149*'`) and a body-content grep for "stage149 / Stage 149" across `scripts/` and `mathematica/` both return empty, so no engine asserts anything for this unit.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Qualitative claim: mouth corrections controlled by a finite rigidity kernel, not branch ambiguity (`stage_149.tex:16`) | none (no script) | n/a — status-only |
| Numeric `A_T`, `B_T`, `|A_T|/B_T` (notes :42-52) | none at stage 149; carried forward from Stages 146-148 (notes :54-60) | n/a — carry-forward, not a 149 deliverable |
| `\stagefield{Verification}`: "SymPy audit: none yet. Mathematica audit: none yet." (`stage_149.tex:11`) | absence of scripts | match — card and manifest agree no engine exists |

The card does not claim an executed verification, so the absence of scripts is exactly what the paper states. Under the status-only rules, `missing_sympy`/`missing_mathematica` would only be valid if a result here referenced an upstream check that no upstream unit's scripts verify. The only numbers (`A_T`, `B_T`) are explicitly inherited from the same-level Stages 146-148 deformation calculation (notes :54-60), so their verification (if any) lives upstream, not at 149. No carry-forward at 149 references a result that this card claims to have independently established. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| — | (none) | — | (no scripts exist for unit 149) | — | — |

No assertions exist to inventory. The card self-declares "none yet" for both engines (`stage_149.tex:11`), consistent with the empty MANIFEST file entries.

## Findings

None. This is a legitimately status-only / no-engine consolidation card. The paper card does not assert any executed verification (it states "none yet" for both engines), so the absence of scripts is not a `missing_verification_script` finding. The numeric values present in the notes are carried-forward from Stages 146-148, not deliverables stage 149 claims to derive, so there is nothing for a 149-local script to verify and no `script_missing_paper_claim` gap.

## Independent-derivation check (Mathematica)

N/A — no `.wl` exists.

## Engine cross-check

N/A — neither engine is present; nothing to cross-check.

## Verdict justification

Verdict is `clean`. I read the stage card, the notes, and the relevant appendix line, then attacked the only available angle for a no-script stage: whether the card claims a verification it does not actually have, or whether it carries a result no upstream check supports. It does neither. The `\stagefield{Verification}` line explicitly says "none yet" for both engines (`stage_149.tex:11`), which matches the MANIFEST's `is_status_only_candidate: true` with all paths null (`MANIFEST.yaml:5059-5083`). The numeric `A_T ≈ -4.27263956256927`, `B_T ≈ 0.134875005736706`, and `|A_T|/B_T ≈ 31.68` are internally consistent (4.27263956.../0.13487500... ≈ 31.679 → 31.68) and are stated by the notes (:54-60) to be inherited from Stages 146-148, not produced at 149. Because the card claims no executed audit and the stage is a status/ledger consolidation, there is no script-side or paper-side discrepancy to raise. No directive is written.

## Self-test notes

I checked: (1) script-existence under multiple naming conventions (`find -iname '*149*'` and body grep) — genuinely none; (2) whether the "none yet" verification line was a stale claim hiding an existing script — it is not, MANIFEST corroborates; (3) internal consistency of the only emitted numbers (`|A_T|/B_T`: 4.27264/0.134875 ≈ 31.68 — consistent); (4) whether any 149 carry-forward references an unverified upstream result — the notes attribute `A_T`/`B_T` to the same-level Stages 146-148, so verification responsibility is upstream, not a 149 gap. Variable-independence / parity / trivial-case traps are inapplicable (no derivatives, integrals, or assertions exist).

## Value Reconciliation (pass-2 augmentation)

No scripts and no saved outputs exist for unit 149, so there are no SCRIPT-emitted values to reconcile. Per the augmentation's status-only clause, the present-engine emission set is empty; an absent engine is not a finding. For completeness I reconcile the numeric values that appear in the notes (which would be the natural carrier had a script existed) against the card and against their own internal consistency. All such values are explicitly carried-forward from Stages 146-148 (notes :54-60), not stage-149 deliverables, and all reside correctly in the notes; the terse card legitimately omits them.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `A_T ≈ -4.27263956256927` | no script (carried-forward, no engine output) | notes `..._stage149_...md:44` | MATCH (lives in notes; carry-forward from 146-148, not a 149 deliverable; card terse, no finding) |
| `B_T ≈ 0.134875005736706` | no script (carried-forward) | notes `..._stage149_...md:45` | MATCH (notes) |
| `|A_T|/B_T ≈ 31.68` | no script (derived from A_T,B_T) | notes `..._stage149_...md:52` | MATCH (notes; internally consistent: 4.27264/0.134875 ≈ 31.68) |
| qualitative result "finite rigidity kernel controls corrections, not branch ambiguity" | no script | card `stage_149.tex:16` and notes :64-71 | MATCH (card + notes) |

INTERNAL / scaffolding items (no finding, name-only): none — there is no script, hence no pass/fail flags, residuals, tolerances, or intermediate driving quantities.

reconciliation: complete; 4 values checked, 0 misaligned
