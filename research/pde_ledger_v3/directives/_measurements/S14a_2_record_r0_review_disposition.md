# S14a part 2 record review, round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_2_record_r0_review_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s14a2_rec_r0_lookups.sh`). "P-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout.

## The round

- **Artifact.** Codex-written (session `01a160d3…`), uncommitted on `978c359c`:
  - the record `steps/S14a_2_slow_limit.md` (`cbd65642…`);
  - its measurements `steps/_measurements/S14a_2_record_measurements.md` (`d63a36d4…`) and generator
    `scripts/S14a_2_record_measurements.py` (`b9cf1509…`);
  - the edits to `SUBSTRATE_REQUIREMENTS.md`, `steps/S14a_1_projection.md` and `STATUS.md`.

  All six files are preserved unchanged in the commit that adds this file. All six matched the baseline after the
  round. The working tree was unchanged apart from the lookups file itself (P1).
- **Role.** A physics-bearing step record and register pass, reviewed until clear (record row).
- **Prompt and legs.** `_scratch/s14a2/rec/s14a2_record_review_prompt_r0.md`, sha256 `2e5b2b5a…` (P1). Both legs
  received it identically, in parallel:
  - a fresh Claude agent (Opus): `/var/projects/_review_scratch/s14a2-record-r0-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s14a2-record-r0-grok/REVIEW.md`, exit 0.
- **Verdicts.** Claude: NEEDS REVISION, 3 findings. Grok: CLEAR (P2).
- **What both legs reproduced.**
  - Each reran the seven filed witnesses under `timeout 600`. Every stdout was byte-identical to the filed one.
  - Each regenerated the measurements file into its scratch. Both outputs were byte-identical to the filed file
    (P2).

## Findings

**V1. An instruction is attributed to the user again.** (Claude 1.) **CONFIRMED.**
- **The record and the register.**
  - Record §6: "The user's current instruction is **"Commit nothing."**", cited as "Retrieve M01; current user
    instruction" (l.245–248).
  - The register's S14a part 2 pass: "resolved by the user's current **"Commit nothing"** instruction" (l.1678)
    (P3).
- **The source.** The phrase is line 5 of the orchestrator's brief to the author. It is in neither the record's
  decision list, nor its lookups, nor M01. The list says "No commit and no push" (l.96) (P3).
- **The precedent.** This is the part 1 record's T3, confirmed and repaired: "The user did not say it" (P3).
- **Cause.** The orchestrator's brief repeated the sentence that caused T3.
- **What must be true.** No sentence attributes it to the user. The list's own "No commit and no push" is the
  source the record cites.

**V2. The deliverer of the selection is a source disagreement, and the record does not report it.** (Claude 2.)
**CONFIRMED.**
- **The two sides.**
  - The option text the user selected (orchestrator-written): "the slab's own dynamics **(the open S12 forms)**
    must pick the slow limit and the far behaviour". The list's item 4: "Its target follows the sources' owners (the
    S12 record for the forms, under "Keep S12")" (P4).
  - The accepted S12 record: "Interior-order law/history/strain; faces' law; source/flux forms | No separate owner
    named for the first two; source/flux forms remain S12's. Slab material completion remains separately S8/S22",
    and "Slab's material law; puncture interior | S8/S22 for the former" (P4).
- **The evidence the requirement rests on.** The filed witnesses vary the material velocity, the material density
  or the face position. In family A, `ω` enters only through `u_r`, with `gr = gd = 1/(r**2*h)`. In family B it
  enters through `n`. The pass-1 slow family keeps "the source densities … exactly the same" while `n_eps` and
  `ur_eps` vary (P4).
- **The record.** §5, §7, the register entry and `STATUS.md` name S12 as the target. They state S8/S22 and the
  unowned laws only as separate qualifications. §6, the section on source comparisons and disagreements, contains
  none of the owner terms: count 0 (P4). The schema's target is "the step that must deliver it", one step (P4).
- **Effect.** The register could hand S12 a requirement whose evidence varies fields whose laws the S12 record
  assigns elsewhere, or to no one. The disagreement is not reported.
- **Cause.** The orchestrator's option text glossed "the slab's own dynamics" as "the open S12 forms". List item 4
  then fixed S12 as the target.
- **Disposition.** Reporting the disagreement with both sides quoted is a repair. Who must deliver the selection is
  a concept decision, and it goes to the user before the repair is written.

**V3. The record and the register state different scopes for the required object.** (Claude 3.) **CONFIRMED.**
- **The record.** "the slow limit and far behaviour selected by the slab's own dynamics**, including the
  controlled return and profile-dependent source-map obligations" (l.186–188) (P5).
- **The register.**
  - "for the controlled `ω→0` return law and the profile-dependent `J→Q` completion" (l.653).
  - "This requirement delivers none of the forms, the slow return law or the source map" (l.666–668).
  - The pass: "for S14a's still-underived controlled return and profile-dependent `J→Q` obligations" (P5).
- **Effect.** "including" would make the target owe the return law and the source map, a new owner split. The
  record's own §7 says "No new split chosen".
- **What must be true.** The record and the register state the same object.

## Notes below the filter, handled with the repair

- **The account of edits (§8).** List item 7 requires the account to match `git diff` change by change. It omits
  these changes (P6):
  - three `STATUS.md` sentences that were added, and one that was moved;
  - the `R-S12-02` rewording "The objects are" → "The printed S12 objects are".

  Grok read §8 as matching. The diff lines are retrieved in P6. **CONFIRMED, minor.**
- **Part 1 record l.366.** "The durable generator is supplied in the working tree, uncommitted". The generator is in
  `4a4dfaf5` (P6). This statement presents part 1's deliverable as not yet committed. It falls under the list's item
  6 property, whose census does not claim to be complete. **CONFIRMED, minor.**

## Notes, no change

- **Interpreter.** Reruns under the default `python3`, with SymPy 1.14.0 in the user site, match the filed stdout
  byte for byte. The filed evidence records no interpreter. Under `python3 -I`, the interpreter loads the system
  SymPy 1.12, and two `p2-grok_dc_limit_objects` stdout lines then differ (P6). No record claim uses those lines, and the record claims no rerun.
- **Hand-typed values in two witnesses (Claude, below the filter).** Codex's pass-2 family B takes a half-angle
  period mean of `1/n` typed by hand. The spreading witness types its material jump coefficients as `(0, 0)`. The
  record states each witness by its own source and standing, a review leg's instrument. Neither changes R1, which
  also rests on family A.

## Agreement between the legs

Grok found nothing outstanding. Claude found V1–V3. They do not conflict: Grok's checks of the witnesses, the
census, the merge and the status lines agree with Claude's "checked and sound" list. V1–V3 are claims Grok's report
does not address.

## Next

V1, V3 and the two minor notes are technical repairs. V2's deliverer question goes to the user first. Then the
author repairs, with no "Commit nothing" sentence in the orchestrator's brief, and round 1 reviews the whole record
again.

## The user's decision on V2 (2026-10-21)

Asked who the register should say must deliver the selection, the user selected **"All the slab's laws
(Recommended)"**. The option text is the orchestrator's: "The selection is owed jointly by S12's source and flux
forms, the slab's material law (S8/S22), and the face and interior laws, which have no owner yet. The register says
so and states the missing owner as missing; it invents none. The record reports the disagreement with both sides
quoted."

This decision replaces the list's item 4 sentence "Its target follows the sources' owners (the S12 record for the
forms, under "Keep S12")". Every other requirement of item 4 still holds: the register's schema, rest-on test and
method, both routes included; the object named, not a route, law, material, premise or candidate; the merge into an
existing entry where one covers the same requirement; and the register winning any disagreement with the list, which
is then reported. The earlier selected option text, "(the open S12 forms)", is quoted as it was given, together with
this decision.
