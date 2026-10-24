# S15 decision list, review pass 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_decision_list_review_p2_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_dl_p2_lookups.sh`). "Q-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout; the legs' scripts were not re-run.

## The pass

- **Artifact.** Revision 2 of `directives/S15_decision_list.md` (sha256 `d83e69e5…`), committed at `cd853e16`, with
  its setup lookups and the pass 1 disposition and lookups. All four matched the baseline after the pass. The working
  tree was unchanged apart from the lookups file (Q0).
- **Role.** The list's spec content, reviewed until clear (spec row). The decision-list fold was pass 1's (G2).
- **Prompt and legs.** `_scratch/s15/dl/dl_review_prompt_p2.md`, sha256 `8b395540…` (Q0). Both legs received it
  identically:
  - Codex (`gpt-6.1-sol`, xhigh): `/var/projects/_review_scratch/s15-dl-p2-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-dl-p2-grok/REVIEW.md`, exit 0.
- **Verdicts.** Codex: FINDINGS, 1. Grok: FINDINGS, 1. Both legs report that revision 2 makes D1–D5 true, and that
  nothing else in the fold changes what is computed.

## Findings

**D6. The plan's step title, which names an outcome, is on the builder's side.** (Codex 1, MEDIUM.) **CONFIRMED.**
- **The list** (`:14`): "S15, "Two drains attract" (`V3_STEP_PLAN.md:652–669`)". It is the title's only
  occurrence, and it lies before the record section, which starts at `:142` (Q1).
- **The sources.** The plan's heading: "S15 · Two drains attract" (`V3_STEP_PLAN.md:652`). M2: "The spec says what to
  **compute** — ⛔ never what anything equals, is expected, or was measured" (`CLAUDE.md:81–82`) (Q1).
- **The legs differ.** Grok read the title as the plan's name, followed by the list's statement that the step prints
  no force law (Grok, "1. Object, not recipe or answer"). The title still states an expected outcome, and its being
  the plan's name does not change that. The orchestrator adjudicates with Codex.
- **What must be true.** No text outside the record section states or implies an outcome of the step. The record
  section's inventory already points at the plan's section, title included.

**D7. The record's conditions omit the stress law's live quantities.** (Grok 1, HIGH.) **CONFIRMED.**
- **The list** (`:145–148`): "which conditions on the profiles, the remainders, `Π_o` and `f` give a named printed
  piece, or a named sum of printed pieces, a stated dependence on the separation … Each such dependence is stated
  with the conditions on the surfaces under which it holds" (Q2).
- **The sources.** The list keeps them live: "`p_br` is a general function of `ρ_br`. `η` and `ζ` are general live
  fields (P6)" (`:64`). S9b: "`p_br` is a general function of `ρ_br` only …"; "The shear viscosity `η` and bulk
  viscosity `ζ` are live general profiles, including their gradients" (`directives/S9b_SHARED_PHYSICS.md:370–371`).
  `p_br`, `η` and `ζ` each have 0 matches in that bullet (Q2).
- **The evidence.** Grok:
  - With the same profiles and the same outer surface, the pressure piece's x-component is `-8*pi*R**3*d/3` (exact)
    for `p = ρ²`, and `-4*pi*R**3*(-R**2 + 5*d**2)/(15*d**2)` for `p = ρ`, through order `s³` in `s = R/d` (Grok
    notes it is a truncation).
  - For one smooth remainder, the viscous stress's divergence is `4*exp(-sx)/3` for `η = 1, ζ = 0`, `exp(-sx)` for
    `η = 0, ζ = 1`, and `0` with both off. The bounds on the shell integral are positive (Q2).
- **Effect.** A dependence on the separation could be stated as following from the profiles and the surfaces when it
  also rests on the pressure law or the viscosities. A flow requirement filed that way would leave out part of what
  it rests on.
- **What must be true.** Every stated dependence of a printed piece, and the flow requirement filed from it, carries
  its conditions on every live quantity the piece contains. These include `p_br`, `η` and `ζ`.

## The repair
Both repairs are to text. Neither changes what the build computes. D7 changes what the record may claim.
