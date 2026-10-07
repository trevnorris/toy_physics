# Independent review: light-leakage scoping inventory (Codex-written v0–v1; repaired by Claude authors in v2 and v3)

## Artifact
`_scratch/light_leakage/light_leakage_scoping.md` in `/var/projects/toy_physics`, written by Codex and then
repaired by Claude authors, to the directive `_scratch/light_leakage/light_leakage_scoping_directive.md`. It is an inventory for a later
worst-case/best-case leakage bracket. It contains no calculations.

## What to check
1. **The directive's requirements.** Is each one met: the columns, the statuses, the seeds, the yardsticks and
   the survival rule? Does each column answer the question the directive asks? For example, does the
   observational-constraint column name observations that bound the loss in question?
2. **Record fidelity.** Does every cited record line (under `research/pde_ledger_v3/`) say what the inventory
   claims, including the statuses it assigns (ESTABLISHED, CONDITIONAL, OPEN, UNRESOLVED)?
3. **Prior-art fidelity.** Do the cited sources say what is quoted? Are unverified items marked? Is any
   prior-art result used as a v3 bound or premise (`CLAUDE.md` M3)?
4. **Classification and destination.** Are the model-mapping and energy-destination entries physically right
   for each mechanism and route?
5. **Completeness.** Is any route a record identifies missing? Is any yardstick the bracket will need missing?
   Is any freeze unnamed?
6. **No verdicts.** Does it state any outcome, suppression requirement, or new physics?

## Required method
- Read the directive and the records first, then the artifact.
- Quote both sides for every finding.
- Each claim about a file or source needs the command or link and the literal text, or it is discarded. Write
  scratch files only under `/tmp`.
- Do not modify any file in the repository. Do not read other review outputs.

## Physics filter
Report a finding only if it changes what the later bracket would compute or claim from this inventory: a route,
a status, a mapping, a destination, a freeze, a validity domain, a yardstick, an oracle correspondence, or a
citation that would mislead the bracket. Do not report wording or formatting that changes none of these.

## Output
- **Verdict:** CLEAR or NEEDS REVISION.
- **Findings,** numbered. For each: the location, what is wrong, the evidence, and the fix.
- A short list of what you checked and found sound.

Keep it under 1,500 words, not counting command output.
