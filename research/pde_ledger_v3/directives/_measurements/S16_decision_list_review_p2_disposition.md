# S16 decision list, review pass 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_decision_list_review_p2_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s16_dl_p2_lookups.sh`). "O-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout.

## The pass

- **Artifact.** Revision 2 of `directives/S16_decision_list.md`, committed at `ef82fb3b` (the G2 fold), with its
  setup lookups and the pass 1 dispositions. All three matched the baseline after the pass. The working tree was
  unchanged apart from the lookups file (O1).
- **Role.** The spec content of a pre-builder decision list, reviewed until clear (spec row). G2's one pass and one
  fold are done.
- **Prompt and legs.** `_scratch/s16/dl/dl_review_prompt_p2.md` (O1). Both legs received it identically, in
  parallel:
  - Codex (`gpt-6.1-sol`, xhigh, workspace-write, working directory its scratch):
    `/var/projects/_review_scratch/s16-dl-p2-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s16-dl-p2-grok/REVIEW.md`, exit 0.
- **Verdicts.** Codex: FINDINGS, 2. Grok: FINDINGS, 1 (O2).
- **The fold, as the legs read it.** Grok found D1, D2, D4, D5 and D6 holding under the user's two later decisions,
  and D3 not holding (E1 below). Codex found the removals and the domain wording present, with D3's separation
  affected by E2. Both found the later option texts and the material balance retained faithfully.

## Findings

**E1. Item 3 calls the received momentum the rate of `P`, which includes the motion term.** (Grok 1.) **CONFIRMED.**
- **The list** (`:111–113`): "This is the rate of `P`, from item 1: the stresses and the momentum carried across both
  surfaces, plus the integrated force density. It is printed as its own object, apart from item 2's kinematic terms
  and from the parts due to the surfaces' motion" (O3).
- **The sources** (O4).
  - The user's decision: "The build prints the momentum the region receives: stresses and momentum carried across both
    surfaces, plus the force density. It prints this apart from the terms caused only by the surfaces moving."
  - D3: "The terms that come from which material the moving surfaces select are kept apart from it, as kinematic
    terms."
- **The evidence.** Grok: on a shell with live radii, constant momentum density and zero flux and force density,
  `dP/dt = 16*pi*c*(t + 1)` while the sum without the motion part is `0`. On a fixed-length moving interval, the
  motion part accounts for the whole rate (O5, E1).
- **Effect.** The sentence defines one object two ways. Two engines can each follow it and print different objects.
- **Cause.** The fold's wording of item 3.
- **What must be true.** The received momentum is the object the user's decision names, kept apart from the terms due
  only to the surfaces' motion. The rate of `P` is a different object, and the list does not equate them.

**E2. The projected inputs do not determine the momentum exchanged at each slab face, or its split by face motion.**
(Codex 1.) **CONFIRMED.**
- **The list** (O3).
  - `:78`: "`f^i` includes every exchange of momentum with the faces and the bulk."
  - `:99–102`: "For each rate, the build prints separately the face terms and the integrated force density … Each
    surface term is printed in two parts".
- **The source.** The part 1 record has "two live graph faces" and "independent regional fields and one-sided traces"
  (`steps/S14a_1_projection.md:46–51`) (O4). The material face terms are supplied there; no momentum counterpart is.
- **The evidence.** Codex: three native momentum currents on the same moving slab have identical material data and
  identical projected `G`, `Π` and `f`. They have different lower and upper face receipts, `(0,0)`, `(0,−1)` and
  `(1,−1)`, and different splits between face transfer at rest and face-motion selection (O5, E2).
- **Effect.** For the momentum rate, "the face terms" and their motion split are not objects the supplied inputs
  determine. Any face-motion selection stays inside `f`, so removing "all" surface-motion terms would leave the faces'
  part hidden.
- **What must be true.**
  - The list states which surfaces the at-rest/motion split covers.
  - For momentum, it states that the face exchanges stay inside `f`, undivided.
  - It keeps the face-resolved momentum accounting OPEN, without supplying a law.

**E3. Compactness is filed as a property of the surroundings field.** (Codex 2.) **CONFIRMED.**
- **The list** (`:142–143`): "Compactness belongs here too, as a statement about that field" (O3).
- **The sources** (O4).
  - The source lists two separate assumptions, "Compact support / small body" and "Smooth external field across the
    support" (`4d_1pn_full.tex:697–704`).
  - The plan: "The worldtube reduction's premise is that the defect is **compact**" (`V3_STEP_PLAN.md:688`).
- **The evidence.** Codex: two material completions with identical data on the annulus and the same surroundings field
  have body-size to field-scale ratios `1/50` and `2` (O5, E3).
- **Cause.** The fold's new sentence.
- **What must be true.** The record files the body's smallness as its own requirement, apart from the regularity and
  normalization of a field of the surroundings.

## Agreement between the legs

Grok found E1. Codex found E2 and E3. They do not conflict. E1 and E2 both concern which motion terms the received
momentum excludes: E1 for the in-plane control surfaces, E2 for the slab faces.

## Rounds

The spec content has now had two reviews: pass 1, and pass 2 of the fold. E1 and E3 are defects in text the fold
wrote. The repairs are technical. They name no new object, choose no law and change no decision. The next review
would be the third round, so the user decides how to finish before it runs.
