# Gravity prior-art section, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `gravity_prior_art_review_r1_lookups.md`, beside this file
(generator `_scratch/prior_art_gravity/gen_r1_lookups.sh`). "P-n" cites a section of that file.

## The round

- **Artifact.** The section after repair 1, which the authoring subagent made against the round 0 dispositions. The
  round 0 version is preserved at `eb5e95db`. The repair is `_scratch/prior_art_gravity/doc_edit_r1.diff`, equal to
  the diff of both files against `eb5e95db` (P3). Both files matched the baseline after the round, and the working
  tree was unchanged apart from the lookups file (P1).
- **Prompt and legs.** `_scratch/prior_art_gravity/review_prompt_r1.md`: the round 0 prompt, with an item 0 that tests
  R1–R7 and asks whether anything else changed (P1). Codex (`gpt-6.1-sol`, xhigh) and Grok (`grok-4.7`, high), in
  parallel, with web access. Both exited 0 (P1).
- **Verdicts.** Codex: one finding. Grok: CLEAR (P2).
- **The repairs.** Both legs found R1–R7 holding. Each checked that the repair changed nothing outside the section
  and the sources file.

## Finding

**S1. G7 omits Wang's assumption of no external body forces.** (Codex 1, Medium.) **CONFIRMED.**
- The section lists Wang's conditions: "irrotational, barotropic, a uniform constant density, an unbounded fluid at
  rest at infinity, a point sink, and a control volume between two spheres…" "Body force" does not occur in the
  section (P4).
- Wang's Theorem 2: "(4) there are no external body forces exerted on the fluid". Theorem 3 and Corollary 4 assume
  the same presuppositions (P4). The sources file already quotes that condition (P4).
- S16 keeps a general force density `f^i` in its momentum balance. So the omitted condition is one that differs from
  S16's, and the section's rule is to state those.
- **What must be true.** G7 states that Wang's results (Theorems 2 and 3, and Corollary 4 where G9 cites it) assume no
  external body forces.

## Not findings

- **Reading not completed.** Both legs record which works they could open only as abstracts or bibliographic records.
  Codex could not recover full text for Reynolds, the cited Lamb edition, Lagally, Papapetrou or Dixon. Grok opened
  Reynolds, Lamb, Lagally and Riemann in full text. The section marks what it did not open. This limits the legs'
  clearance and is carried, as in round 0.
- **Line 173** still stands (both legs).

## Outcome

One finding, confirmed. The round 1 version is preserved with these dispositions. The repair adds one condition. A
review of it would be round 2, the third round, so the orchestrator asks the user how to finish.
