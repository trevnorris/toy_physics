# S16 record review, round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_record_r0_review_lookups.md`, beside this file (generator
`_scratch/s16/gen/rec_r0_lookups.sh`). "Y-n" cites a section of that file.

## The round

- **Artifact.** These were written by Codex (`gpt-6.1-sol`, xhigh), starting from the record decision list at
  `4d4d22cb`:
  - the record `steps/S16_moving_region.md`;
  - its measurements file `steps/_measurements/S16_record_measurements.md` and generator
    `scripts/S16_record_measurements.py`;
  - the S16 pass and other edits in `SUBSTRATE_REQUIREMENTS.md`;
  - the edits to `steps/S14a_2_slow_limit.md` and `../../STATUS.md`.
- **Prompt and legs.** `_scratch/s16/rec/s16_record_review_prompt_r0.md` (Y0). Codex wrote the artifact, so the legs
  are a fresh Claude agent (Opus) and Grok (`grok-4.7`, high) (G1). They had the identical prompt and ran in parallel.
  Grok exited 0.
- **Baseline.** All six baseline files matched after the round. The working tree was unchanged (Y0).
- **Verdicts.** Claude: NEEDS REVISION, four findings (one minor). Grok: CLEAR (Y0).
- **Both legs' generator reruns** reproduced the filed measurements file byte for byte.

## Findings

**Y1. The plan's S19 gate and charter feedback are not handed on.** (Claude 1.) **CONFIRMED.**
- The plan (Y1): "S16 must be run BEFORE S19 is trusted, and its result feeds back into the charter. … Say so rather
  than preserving the boundary."
- The record says that the regional objects "do not deliver the plan's compact-body response approximation". It says
  that "S19 owns the C10 tail question" and that "No broader charter scope or later-step result is earned by the
  regional identities" (Y1). It does not say what this means for the plan's S19 gate or for the charter feedback, or
  who disposes of either. `STATUS.md`'s edit mentions neither (Y1).
- The record list, item 1: the record "states what of the plan's S16 rows this step prints and what stays with later
  steps, as its sources state them; where they do not say, the record says so."
- **What must be true.** The record and `STATUS.md` hand on the plan's S19 gate and its charter feedback beside what
  this step delivers. Each is handed on with its owner as the sources state it, or "no owner named" where they state
  none. Nothing about either is decided here.

**Y2. `R-S1-01` gains no S16 source, unlike the earlier steps on the same coordinates.** (Claude 2.) **CONFIRMED.**
- `R-S1-01`'s sources list S10, S11, S9b, S12 and S14a part 1, each as "coordinate dimension only" (Y2).
- S16 is built on S14a part 1's projected domain: "in-plane coordinates `x^i` (three) and time `t`" (the step list's
  supplied equations, Y2).
- The pass rules out an S16 source on the ground that "Control shape, moment truncation and independent dimension
  units are labelled supplied/construction choices" (Y2). The coordinate dimension is none of those, and the pass
  does not test it.
- `R-S1-02` lists S14a part 1 for its "support identification only" (Y2). S16 is built on that same projected
  balance. The pass does not test this either.
- **What must be true.** The pass applies the register's rest-on test to the coordinate dimension, and to the support
  identification S16 inherits through S14a part 1's balance, as the earlier passes did. The census, the record and
  `STATUS.md` follow from the result.

**Y3. Requirements the user did not file are attributed to the user.** (Claude 3.) **CONFIRMED. Cause: the step
list's heading.**
- The user's option texts (Y3) send two things to the register:
  - what a passive factor needs: "a named field of the surroundings, a rule separating it from the body's own field,
    and its normalization";
  - the inertial law: "The law it needs goes to the register, as O2 keeps it (S8)".
- The step list's own filing adds the body's compactness (its pass-2 finding E3), the faces' momentum exchange (E2)
  and the slab's density response (O2 `:138`). The orchestrator wrote that filing. The list's heading over all of
  them reads "Requirements the user's decisions file" (Y3). That heading is the orchestrator's paraphrase, and the
  record followed it.
- The register, the record and `STATUS.md` call these "S16 user requirement" (`R-S16-02`, `-03`, `-05`), "an
  explicit user-filed requirement" (`R-S16-02`), "the user's explicit future requirements" (the pass), and "the
  user's … compactness, face-momentum" (`STATUS.md`) (Y3).
- The pass rests its admission of the entries on that attribution (Y3).
- **What must be true.**
  - Each entry's source names who filed it: the user's decision, quoted, or the step list's filing, with the line
    or finding it comes from.
  - The disagreement between the step list's heading and the user's option texts is reported, with both sides
    quoted.
  - Each entry's admission ground rests on what its own source says.

**Y4. Twenty-two compared bindings get no cross-engine status.** (Claude 4.) **CONFIRMED; minor.**
- The comparator compares the supplied equations (4 bindings), the integral definitions (3) and the geometry (15).
  All 22 have an exact residual of `0` as computed, and every sample residual is `0` (Y4).
- The record names none of these classes (Y4).
- The record list, item 2: "The record covers every row class the comparator prints".
- **What must be true.** The record states these comparisons' cross-engine status, with their counts and domain.

## Not findings

- Both legs found the following sound:
  - the comparison counts;
  - the unreconciled domain difference;
  - the carried build and comparator notes;
  - the premise quotes, the inventories and the OPEN owners;
  - the prior-art citations, with the passive-factor row kept apart from `R-S16-01`;
  - the S14a part 2 status lines and the `R-S14a-01` row;
  - the account of edits against `git diff 4d4d22cb`.
- Grok reported no finding. Y1–Y3 are claims about what the sources assign and who filed what. Each is grounded in
  the verbatim lines in Y1–Y3.

## Outcome

Four findings confirmed. Y1–Y3 change what the record, the register and `STATUS.md` hand on. Y4 is an under-claim.
The reviewed version is preserved with these dispositions. The record's author repairs it to what must be true, and
the repair gets a round 1 review by a fresh Claude agent and Grok.
