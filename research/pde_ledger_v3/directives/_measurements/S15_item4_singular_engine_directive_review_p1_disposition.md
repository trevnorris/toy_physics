# S15 item-4 Singular engine directive, review pass 1: disposition

## The gate record

- **Artifact.** `research/pde_ledger_v3/directives/S15_item4_singular_engine_directive.md`, revision 1, sha256
  `1df009e9…`, with its lookups file (`S15_item4_singular_engine_directive_lookups.md`, sha256 `b8b6b043…`) and two new
  `.gitignore` lines. The orchestrator wrote them. The directive is physics-bearing (it names the objects a second
  item-4 engine computes), so it is reviewed until clear.
- **Prompt.** `/var/projects/_review_scratch/s15-singular-dir-prompt/singular_engine_directive_review_prompt_p1.md`,
  sha256 `a4c55228…`, identical for both legs.
- **Legs.** Codex (`gpt-6.1-sol`, xhigh) and Grok (`grok-4.7`, high), in parallel, 2026-11-09. Reports:
  `/var/projects/_review_scratch/s15-singular-dir-p1-codex/REVIEW.md` (sha256 `b10eade6…`, CLEAR) and
  `/var/projects/_review_scratch/s15-singular-dir-p1-grok/REVIEW.md` (sha256 `c82d1dc4…`, three findings).
- **Baseline.** Revision 1 is committed with this file, before the fold.

## Findings (Grok), verified against the sources

**G1 (confirmed; repair).** Revision 1 paraphrased item 4 and dropped the reading decision list's sentence "A condition
is on the exponents and, where they enter, on the amplitudes, the radii and their rates, and the velocities"
(`S15_reading_decision_list.md:80–81`). Without it, "the named quantities" in revision 1 has no referent. Repair:
revision 2 quotes item 4 verbatim and keeps "complete".

**G2 (confirmed; repair).** Item 4 covers each `Q` whose item-3 object is in closed form
(`S15_reading_decision_list.md:78`). Revision 1 asked for item-4 conditions for every `Q` that passes item 3's
integrand test, and its label rule covered only conditions, so an item-3 object the engine did not obtain had no
label. Repair: item 4 is asked for where the engine obtains the item-3 object exactly, and the label rule covers an
item-3 object it does not obtain.

**G3 (confirmed; repair).** The guard kills a job that runs out of memory (`scripts/s11c_guarded_run.py`,
`OOMPolicy=kill`), so the engine itself may print nothing after the kill. Round 4 of the reading build measured the same
pattern: partial blocks left inline and unmarked (`S15_reading_build_review_r4_disposition.md`, "Recorded, no repair
required"). Revision 1 did not require a label for an object whose attempt ended that way. It also did not restrict an
out-of-memory end to a run at the guard's largest admitted cap with no other guarded job running, and its report
listed only labelled objects. Repair, in revision 2:
- every object in scope ends in the final transcript with its objects or a label, and an attempt that ended otherwise
  gets a label beside what it printed, saying how it ended;
- an out-of-memory end counts only at the largest admitted cap with no other guarded job running throughout;
- the run record and the report list every object in scope and how it ended.

## Checked and found sound (Codex)

Codex found the object, the supplied premises, independence, claims and labels, model point, leakage and operation
sound. Its inventory check matched item 3's integrand test to the accepted Wolfram build's closed-form statuses, 144
and 144, with none on only one side.

## Next

Revision 2 gets pass 2 from both legs, Codex and Grok, reviewed until clear.
