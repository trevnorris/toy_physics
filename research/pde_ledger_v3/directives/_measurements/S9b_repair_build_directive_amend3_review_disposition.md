# S9b repair build directive, amendment 3: review dispositions (orchestrator)

**Artifact:** amendment 3 to `directives/S9b_repair_build_directive.md`, on top of the accepted `d14838f3`. It
prepares the build's repair round 1 after review r0 (`S9b_repair_build_r0_review_disposition.md`):
- item 4 supplies a shared vocabulary;
- item 13 requires general profiles (the user's choice, 2026-10-09);
- Parts 2 and 3 list each engine's repairs.

The user approved this repair round on 2026-10-09. The amendment is orchestrator-written, so the legs are Codex
and Grok (G1). It changes what the engines compute and how their outputs are paired, so it is reviewed until
clear.

Commands and literal output for each round are in the `_lookups.md` files beside this one. Generators:
`_scratch/s9b_build/gen/s9b_amend3_r*_lookups.sh`.

## Round 0

Both legs used the identical prompt `_scratch/s9b_build/s9b_repair_build_directive_amend3_review_prompt_r0.md` on
the uncommitted version with sha `16f908a2…`. A frozen copy is at
`_scratch/s9b_build/S9b_repair_build_directive_amend3_reviewed_r0.md`.
- **Codex** (gpt-6.1-sol, xhigh): NEEDS REVISION, 1 finding. Final report
  `_scratch/s9b_build/s9b_amend3_review_r0_codex_final.txt`; evidence `…_r0_codex_evidence/`.
- **Grok** (grok-4.7): NEEDS REVISION, 3 findings. Report `_scratch/s9b_build/s9b_amend3_review_r0_grok.txt`.

Both reported before I adjudicated.

| # | Finding | Disposition | Repair |
|---|---|---|---|
| C1 | `A_NONRECIPROCAL_G<abc>` and `B_GAMMA_DIFFERENCE` name differences without fixing their order. The spec says "half the difference" (L174) and "the difference between the two `γ`s" (L225). One name could then pair opposite signs across engines. | **ACCEPT.** The spec and vocabulary read as quoted (lookups). | The vocabulary fixes each order: `A_NONRECIPROCAL_G<abc>` is half of the emitter-to-reflector time minus the reflector-to-emitter time, and `B_GAMMA_DIFFERENCE` is the deflection `γ` minus the radar `γ`. This is a naming convention, not a value. |
| G1 | Item 13's "every shared-vocabulary object is computed for general radial profiles `δ`, `V`, `ξ_w`" collides with item 10's restrictions and forward stages, which set profiles to zero or to the solved `V` (L128, L150–152), and with Part C's constant response. Nothing says which wins. | **ACCEPT.** L182–183 against L128 and L150–152, as quoted. | The general-profile rule applies to the profiles an object leaves live. Item 10's substitutions and Part C's responses apply first. An item 10 condition that cannot be reduced is not an item 12 stop (item 10, L159). It is omitted and reported under item 17. |
| G2 | The `F_` prefixes reach every `A_` name, including `A_NONRECIPROCAL_PATH_DEPENDENCE`. Item 10's forward print list (L132–134) does not ask for it, so two unrequested names are created. | **ACCEPT.** L98 "every `A_` name", L84, and L132–134, as quoted. | The `F_` prefixes apply to the graded `A_` names and the radar-log grades only. |
| G3 | Part 3's sentence "each effective `γ` is solved on every stratum of its observable's coefficient … use that solution" widens the r0 G3 outcome, which concerns the radar `γ` only. "Δθ's coefficient" is not a spec object, and "that solution" is singular. | **ACCEPT.** L290–291 against the disposition's G3 outcome and spec L217–218, as quoted. My rewording widened it. | The sentence returns to G3's scope. The radar `γ` is solved on every stratum of the `ln(1/b²)` coefficient. The γ difference, and each restriction and forward copy of the radar `γ`, use that solution. |

**Next.** The reviewed version is committed unchanged as a preserved baseline, not accepted. The repair goes to
round 1 with fresh legs on one prompt.
