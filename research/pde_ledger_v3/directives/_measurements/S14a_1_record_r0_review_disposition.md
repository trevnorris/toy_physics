# S14a part 1 record, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_record_r0_review_lookups.md`, beside this file
(generator `_scratch/s14a/gen/rec_r0_lookups.sh`). "K-n" cites a section of that file. Every computed result quoted
here is a leg's saved stdout.

## The round

- **Artifact.** The record, its measurements file and generator, the S14a part 1 register pass and the status and
  provenance edits, written by Codex (`gpt-6.1-sol`, effort xhigh, session `01a16048…`) to the record decision list
  at `f5d283c1`. It exited 0 after 239,180 tokens (K1). Reviewed as the seven files of the round-0 baseline:
  - `steps/S14a_1_projection.md`;
  - `steps/_measurements/S14a_1_record_measurements.md` and `scripts/S14a_1_record_measurements.py`;
  - `SUBSTRATE_REQUIREMENTS.md`, `steps/S12_order_source.md`, `directives/_measurements/S14a_1_build_runs.md` and
    `../../STATUS.md`.

  All seven matched the baseline after the round. The working tree was unchanged apart from the lookups file itself
  (K1).
- **Prompt and legs.** `_scratch/s14a/rec/s14a1_record_review_prompt_r0.md`, sha256 `eeaeebef…` (K1). The record is
  Codex-written, so the legs are a fresh Claude agent (Opus) and Grok (`grok-4.7`, effort high) (G1), run in parallel
  with the identical prompt. Reports: `/var/projects/_review_scratch/s14a1-record-r0-claude/REVIEW.md` and
  `/var/projects/_review_scratch/s14a1-record-r0-grok/REVIEW.md`.
- **Verdicts.** Claude: NEEDS REVISION, 4 findings. Grok: CLEAR (K1). Both legs reran the generator into scratch and
  got the filed measurements file byte for byte (K1).

## Findings

**T1. The seven unresolved exact rows do not carry which samples can see a derivative at another point.** (Claude 1.)
**CONFIRMED.**
- **The record.** b50, b72, b128 and b148 each "has 54 stored sample residuals `0`"; b171, b173 and b174 likewise
  (l.100–105, l.122–125). In §5 they are open exact differences checked by finite samples (l.197–199). The controls
  paragraph ends "its coupled quadratic probes detect derivative-point changes where the unresolved exact residual
  cannot resolve them; affine/cubic forms are unchanged" (l.276–280, K2).
- **The sources.** The comparator's round 1 carries to this record: "A derivative moved to another point moves these
  samples through the quadratic form only" (l.96). Its sentence "The affine and cubic forms are unchanged" (l.44)
  says the repair left those forms alone. Round 0 found those forms additively separable, so that a face derivative
  at another point "changes no sample" (l.60–64) (K2).
- **The evidence.** Only the quadratic form has a cross-argument product (396 of 396 functions; none in the affine or
  cubic form). Each of the seven rows splits its samples equally across the three forms: 18/18/18, and 36/36/36 for
  b173 (K2).
- **Effect.** A reader counts 54 samples per row toward a check that only 18 of them can make. The record's reuse of
  "affine/cubic forms are unchanged" does not say this, and read as a statement about the samples it says the
  opposite.
- **What must be true.** Wherever the record states the sample check for these seven rows, it states that only the
  quadratic-form samples can move when a derivative is taken at another point, and how many those are. No sentence
  presents the affine or cubic samples as making that check.

**T2. The U1 interior bound is given a dependence that no computed object has.** (Claude 2.) **CONFIRMED, minor.**
- **The record.** The premises table row "U1 interior `0≤f≤1`, permitted jumps and no concentrated input fields"
  lists "The regional/traced function class of the same objects" as what depends on it (l.58, K3).
- **The evidence.** Neither engine nor either transcript has a line where the order field carries an inequality or
  an assumption. The only declarations are the positive radii and an angle range (K3).
- **The list.** Item 5: "One that enters no computed object is listed as such" (K3).
- **What must be true.** The bound is listed as entering no computed object. The rest of the row keeps the dependence
  its sources give.

**T3. An instruction is attributed to the user.** (Claude 3.) **CONFIRMED, minor.**
- **The record.** "the current user expressly says **"Commit nothing."**" (l.352–355, K4).
- **The source.** The sentence is line 5 of the orchestrator's brief to the record author. It is not in the
  measurements file (count 0, K4). The user did not say it.
- **What must be true.** No sentence attributes it to the user. If the record cites it, it names its actual source.

**T4. The headline calls the channels "material-loss".** (Claude 4.) **CONFIRMED, minor.**
- **The record.** "the two material-loss channels" (l.10). It is the only occurrence (K5).
- **The sources.** The step's decision list: "both channels by which material leaves the slab around a throat"
  (l.54). Stage 006: "Total material is conserved" (l.87) (K5).
- **Effect.** "Loss" reads as material destroyed, which the committed law does not do. That is the conservation
  conflict the plan names as the reason for S14a.
- **What must be true.** The channels are named as the step's decision list names them: material leaving the slab.

## Checked and found sound by both legs

The enumeration of the comparator output against every count the record gives; the objects and their notation; the
labelled specialization, the conditioned velocity objects and the plan's standing for S14, S15 and S20; the premise
standings apart from T2; the measurements' retrieval-only commands and their byte-for-byte regeneration; no printed
difference reconciled; the build's and comparator's notes carried; the inventories and OPEN items; the register pass
(18 → 18, no new entry, merges into existing entries, no means named); the status and provenance lines, J1 and J2;
and the account of edits against `git diff f5d283c1`.

## Outcome

T1 changes what may be claimed from the sample check. T2–T4 each correct one claim: a standing, an attribution and a
name. The reviewed version is preserved before the repair. The repair goes to the same author, in a resumed session,
with these dispositions. Then review round 1, in full.

**For the orchestrator's own record.** The sentence the record reused in T1 is the orchestrator's, from the
comparator's round-1 disposition (l.44), where its referent is the repair. Dispositions that hand limits on should
state the limit, not the repair's history.
