# S9b repair build directive, amendment 1: review dispositions (orchestrator)

**Artifact:** amendment 1 to `directives/S9b_repair_build_directive.md`, on top of the accepted `45daca1f`.
- **Amendment.** At the user's request (2026-10-08), item 10's forward case prints Part B's condition, where it
  printed none before. The user's rule: never change the model's concepts to fit GR, but values within a fixed
  concept may be fitted to it.
- **Memory limit.** Item 11's limit rises to 16 GiB. The user approved this on 2026-10-08, after the SymPy harness
  was killed at 8 GiB.
- **Authorship and review.** The amendment is orchestrator-written, so the legs are Codex and Grok (G1). It changes
  what the engines print, so it is reviewed as physics until clear.

Commands and literal output for each round are in the `_lookups.md` files beside this one. Generators:
`_scratch/s9b_build/gen/s9b_amend1_r*_lookups.sh`.

## Round 0

Both legs used the identical prompt `_scratch/s9b_build/s9b_repair_build_directive_amend1_review_prompt_r0.md` on
the uncommitted version with sha `18b00243…`. A frozen copy is at
`_scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r0.md`.
- **Codex** (gpt-6.1-sol, xhigh): NEEDS REVISION, 1 finding. Final report
  `_scratch/s9b_build/s9b_amend1_review_r0_codex_final.txt`; evidence `…_r0_codex_evidence/`.
- **Grok** (grok-4.7): CLEAR. Report `_scratch/s9b_build/s9b_amend1_review_r0_grok.txt`.

Both reported before I adjudicated.

| # | Finding | Disposition | Repair |
|---|---|---|---|
| C1 | The forward case prints Part B's condition, but not its Part C rewrites. The spec's Part C rows keep `V` live (L226–227), and item 10's Part C addition is restricted to `V ≡ 0`, `ξ_w ≡ 0` (L102–104). So no output covers Part C under the forward premise. | **ACCEPT.** The directive and spec read as quoted (lookups). The leg's script prints 12 forward Part C rows, and each carries `Φ`. The premise therefore changes those conditions. Grok's CLEAR says Part C "stays on that live-`V` branch". That describes the text; it does not show that the forward rows are not needed. | The second-stage forward condition is also printed under each of Part C's three responses. Each is reduced as in item 6, with domains built as in item 9 and the implied-`j_n` exception. |

The leg's evidence includes computed forward-case conditions. These are reviewer measurements. They stay out of
builder packets and leg prompts.

**Next.** The reviewed version is committed unchanged as a preserved baseline, not accepted. The repair is then
reviewed in round 1 by fresh legs on one prompt.

## Round 1 (acceptance)

The round-0 version is preserved at `8eaa17aa`. Both legs used the identical prompt
`_scratch/s9b_build/s9b_repair_build_directive_amend1_review_prompt_r1.md` on the uncommitted repair (sha
`50f3c6ae…`; frozen copy `_scratch/s9b_build/S9b_repair_build_directive_amend1_reviewed_r1.md`).
- **Codex:** CLEAR. Report `_scratch/s9b_build/s9b_amend1_review_r1_codex_final.txt`; evidence
  `…_r1_codex_evidence/`.
- **Grok:** CLEAR. Report `_scratch/s9b_build/s9b_amend1_review_r1_grok.txt`.

Neither leg reports a finding. Lookups: `S9b_repair_build_directive_amend1_r1_review_lookups.md`. The amendment is
clear and the directive is accepted for the build at this version. Findings by round: 1, 0.
