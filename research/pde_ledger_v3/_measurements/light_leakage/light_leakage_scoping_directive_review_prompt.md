# Independent review: a scoping directive (light leakage)

## Artifact
`_scratch/light_leakage/light_leakage_scoping_directive.md` in `/var/projects/toy_physics`. The orchestrator
wrote it for Codex. This is a one-pass review of the directive before Codex starts. It is not a review of any
output.

## What to check
1. **The record citations in the directive.** Do they say what the directive says they do? Check them under
   `research/pde_ledger_v3/`.
2. **The object.** Is the requested inventory the right preparation for a worst-case/best-case leakage
   bracket? Is anything needed for that bracket missing, such as a route, a parameter class, or a yardstick?
3. **Leakage and premises** (`CLAUDE.md` M2 and M3). Does the directive state or imply any expected outcome?
   Does it let prior art become a premise? Does it hide a freeze?
4. **Feasibility.** Can Codex do this from the records plus web search, without computing anything?

## Required method
- Quote both sides for every finding.
- Each claim about what a file says needs the command and its literal output, or it is discarded. Write
  scratch files only under `/tmp`.
- Do not modify any file in the repository.

## Output
- **Verdict:** SOUND or NEEDS CHANGE.
- **Findings,** numbered. For each: the location, what is wrong, the evidence, and the fix.
- A short list of what you checked and found sound.

Keep it under 800 words, not counting command output.
