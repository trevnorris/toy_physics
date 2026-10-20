# S12 build, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_build_review_r0_lookups.md`, beside this file
(generator `_scratch/s12/gen/build_r0_lookups.sh`). "B-n" cites a section of that file.

## The round

- **Artifact.** The Codex-written S12 build, uncommitted at review on top of `a48ce830`: both engines
  (`scripts/S12_order_source.py`, `mathematica/S12_order_source.wl`), their transcripts in `scripts/out/` and
  `mathematica/out/`, and the run record `S12_build_run_record.md`. It is preserved unchanged in the commit that adds
  this file. All five baseline hashes matched after both legs, and the working tree was unchanged apart from this
  file and the lookups file (B1).
- **Prompt and legs.** `_scratch/s12/build_review/s12_build_review_prompt_r0.md`, sha256 `b6377445…` (B1). Both
  legs received it identically and ran one after the other, because both ablate Wolfram.
  - **Claude:** a fresh opus agent. Report `/var/projects/_review_scratch/s12-build-r0-claude/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-build-r0-grok/grok_review.txt`.
- **Verdicts.** Claude: NEEDS REVISION, one finding. Grok: CLEAR, no findings (B2).
- **What both legs report as sound.** Each leg derived the object in its own script before opening either engine,
  and found the engines' printed region and face objects equal to its own, at the generality of fields it used. Both
  report:
  - the supplied balances and split transcribed exactly;
  - every field the list names live, faces live in the in-plane coordinates and time;
  - the continuous and the jumping order both covered;
  - item 2 general in both offsets and labelled as a specialization, with no limit;
  - dimensions computed from a stated assignment;
  - no conclusion token and no value-dependent emission;
  - a Wolfram engine that imports nothing and constructs differently.

  Every FORM ablation and one-sided corruption either leg ran moved the printed output. Guarded reruns by both legs
  reproduced both transcripts byte for byte.

## Findings

**F1. The order flux is held continuous across both faces.** (Claude 1, Medium.) **CONFIRMED.**
- **The list.**
  - The object has `J_χ` "live" (B3, line 75).
  - Under the model point: "Freezes: none intended. `J_χ` stays live" (line 104).
  - Under Retained OPEN: "the form of `J_χ`" (line 113).
- **The engines.** Each builds one order flux, differentiable across the whole domain (`flux = …` at
  `S12_order_source.py:149` and `S12_order_source.wl:77`; B7). Each states that singularities of the bulk fields are
  outside its calculus (`.py:9–12`, `.wl:7–9`) (B4).
- **The gap.** How `J_χ` behaves at a face is part of its form, which the list keeps OPEN. A flux built continuous
  across each face has its normal jump there set to zero. The face-concentrated part of the computed source then
  carries only the order's jump. The leg's own script shows the extra face term when the flux may jump, and an engine
  copy with a face-confined flux gains it (B5). Both engines make the same choice, so a later comparison of the two
  would read agreement on it (M1, M3).
- **What enabled it.** Item 1 names one source of face-concentrated content, "if the order jumps there", and is
  silent on the other fields' face behaviour. The engines covered the case it names.
- **What must be true.** The face objects cover every face-concentrated part the committed law produces with `J_χ`'s
  face behaviour left open. That includes its jump across each face. Nothing about that behaviour is chosen: no
  continuity, no vanishing outside, no boundary condition. Every other object, field and label stays as it is.
- **Grok's CLEAR** does not bear on this. Its own derivation builds a single flux over the whole domain, the same
  generality as the engines (B6).

## Notes from the Claude leg for the record (no change to what is computed)

- SymPy computes each region object by substituting that region's order. Its region objects therefore do not move
  under a SymPy support ablation, and they are no test of its support encoding. Their content is correct, by the
  leg's restriction check. Wolfram's region objects are restrictions, and they move under its support ablation.
- Both engines solve `J_χ`'s units from the order balance. The homogeneity result is therefore definitional for
  `J_χ` and informative for the other terms. The engines print that inference as such.
- SymPy's face traces carry a step in `w₊ − w₋`, which the supplied strict face order sets to one. Wolfram's do
  not. The difference is representational.

## Outcome

One finding, confirmed. It changes what is computed at the faces and what may be claimed, so the build is not
clear. Round 0 is preserved as is. The repair goes to the same builder: this is its first repair, so no change of
author is due (G4). The repaired build gets its own two-leg review.
