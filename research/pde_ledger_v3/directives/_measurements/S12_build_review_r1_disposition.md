# S12 build, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_build_review_r1_lookups.md`, beside this file
(generator `_scratch/s12/gen/build_r1_lookups.sh`). "C-n" cites a section of that file.

## The round

- **Artifact.** The S12 build after repair 1, uncommitted at review on top of `63d99406` (round 0, preserved). It is
  preserved unchanged in the commit that adds this file. All five baseline hashes matched after both legs, and the
  working tree was unchanged apart from this file and the lookups file (C1).
- **Prompt and legs.** `_scratch/s12/build_review/s12_build_review_prompt_r1.md`, sha256 `531de660…` (C1). Item 0
  asked for the repair; items 1–8 asked for the whole build again, with a new question: is any field's behaviour at a
  face chosen where the list leaves it open? Both legs received the prompt identically and ran one after the other.
  - **Claude:** a fresh opus agent. Report `/var/projects/_review_scratch/s12-build-r1-claude/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-build-r1-grok/grok_stdout.txt`.
- **Verdicts.** Claude: NEEDS REVISION, two findings. Grok: CLEAR, no findings (C2).
- **What both legs report as sound.**
  - The round-0 defect is repaired. Each engine gives `J_χ` an independent restriction and free one-sided traces in
    each region, and each face coefficient carries `J_χ`'s jump with nothing chosen about it.
  - Setting the regional fluxes equal and the added face flux (F2) to zero gives back every round-0 object.
  - Each leg's own derivation agrees with both engines' region and face objects.
  - Dimensions, supplied inputs, emission, independence and reproducibility hold as in round 0.
  - Every control either leg ran moved the printed output.

## Findings

**F1. Density and velocity are held differentiable across both faces.** (Claude F1, Low.) **CONFIRMED, as a limit
on the claim.**
- **The list.** `n` and `u` are among the "live functions" (line 75–76). Their behaviour at a face is not stated,
  and the flow profile and the slab's material law stay OPEN (lines 117–118) (C3).
- **The engines.** Both state that `n` and `u` keep their differentiable class (`.py:9`, `.wl:7–8`) (C4).
- **What it changes.** The leg's derivation lets every field jump at each face (C6):
  - with `n` and `u` read as interior-side traces in the face terms, no printed `n Γ_B` form changes (line 47);
  - the material balance gains face content of its own, which vanishes identically when `n` and `u` are
    differentiable (line 23).

  So the printed forms stand. But the object is computed only for `n` and `u` that do not jump at a face, and the
  material balance's own face content is not printed. The round-0 Claude leg asked the same question and found no
  change to the face coefficient (C7). This leg weighs what may be claimed, and is right to: a jump in density at a
  slab face is not excluded by anything the list says.

**F2. The repair adds an order flux concentrated on each face, and the headline face objects include it.** (Claude
F2, Medium.) **CONFIRMED.**
- **The list.** `J_χ` is one of the "live functions"; its form stays OPEN (lines 75, 113) (C3).
- **The engines.** Each adds, for each face, a flux the scripts call "an arbitrary distribution supported at the
  face" (`.py:9–13`, `.wl:7–12`) (C4).
  - **Its construction.** The added fields are functions of all five coordinates. They are added to `J_χ` with no
    localising factor: `.py:209–214`; in Wolfram the markers are set to 1, `.wl:117–120`.
  - **Its support.** It appears only as a printed pairing of each field with its face's equation (transcript
    line 25 in each engine). Nothing computes it.
- **What it changes.** The divergence of the added flux enters each `face/k/source_contribution`, the object a
  record would quote. A flux concentrated on a face can supply any face source: the leg's check gives a free
  coefficient at each order of `δ` (C6). So the headline face objects constrain nothing at a face. The list's own face
  object is intact and printed separately as `face/k/trace_coefficient`.
- **What enabled it.** The repair wording, mine (C7): "every face-concentrated part the committed law produces when
  `J_χ`'s behaviour at each face is left open". The builder read "behaviour at each face left open" as admitting a
  flux concentrated on the face. Such a flux would carry order along the face. That is a new mechanism, not a field
  the list names.
- **Grok** describes the addition (C5) and does not count it as a defect.

## The pattern

All three findings in rounds 0 and 1 answer one question: what class of function each live field belongs to at a
face. The list never says, and each build or repair chose a class. Round 0 chose differentiable fields. Repair 1
opened `J_χ` and went further than the list, to face-concentrated flux. That is a question for the decision list,
not for a repair brief: the class belongs in the list's object, stated once for every field, and reviewed.

## Outcome

Two findings, confirmed. Round 1 is preserved as is. A further repair would make round 3 of this build's review, so
it is raised with the user before it runs, together with the one concept it raises: whether order may be carried
along a face.
