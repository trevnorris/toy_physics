# S9b record, review round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_record_r1_review_lookups.md` beside this file (generator
`_scratch/s9b_build/gen/s9b_record_r1_lookups.sh`).

## The round
- **Artifacts.** The record, its measurements file and generator, and the S9b register pass, as repaired in round 1
  by the same author (brief `_scratch/s9b_build/record/resume_record_r1.md`, against the round-0 disposition's "What
  must be true"). The reviewed previous version is `36f59ff1`. Four files in the baseline
  `_scratch/s9b_build/record/s9b_record_review_baseline_r1.sha256`, all OK at review.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/record/s9b_record_review_prompt_r1.md` (sha256 `956a83f5…`):
  round 0's prompt with the baseline and previous version advanced. The legs were a fresh Claude agent (opus) and
  Grok (`grok-4.7`, effort high), launched identically, with scratch under
  `/var/projects/_review_scratch/s9b-record-r1-{claude,grok}/`. Reports:
  `_scratch/s9b_build/record/s9b_record_review_r1_claude.md` and `…_r1_grok.txt`. Grok returned **CLEAR**; Claude
  returned **NEEDS REVISION** with four findings.
- **Leg note.** The Claude leg reports that the harness saved one oversized output of its own `cat` of repository
  disposition files under `~/.claude/projects/…/tool-results/`, and that it read that file there and nothing else in
  that directory. The file held repository content only.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **The deflection/radar separate-test count is handed on as OPEN, though the decision list's counting rule settles it from the printed objects.** The record and register say the Wolfram separate-test count "remains OPEN" (record 10–11, 183–186, 256–259, 339, 581–582; register 575, 838–839). Decision-list item 5 counts the two observables as separate tests "only where the printed objects show their conditions differ". The γ-difference object prints zero on `GM ≠ 0` in all six B-family prefixes in both engines. The differing Wolfram fields are `InverseWeight`, `ReferenceInverse` and the velocity reduction, and the record itself says differing serialization does not establish independent tests. | **Confirmed.** The lookups show the γ-difference literal in both transcripts and in all twelve prefixes (the leg's census), and the record's "OPEN" wording at each location. Under item 5 the count follows from what the printed objects show. Whether the differing Wolfram fields state the same condition stays open, but it is a separate question. This is material round 1 changed. My round-0 instruction ("What the record and the register conclude about counting … follows from that") did not name item 5's rule, and it led the author to the open count. Grok read the open count as acceptable. |
| 2 | Claude 2 | **The flow-only restriction is reported only in Wolfram's unevaluated form.** §3 (189–216) quotes Wolfram's `Y = -r Inactive[Integrate][…]`, `V = ±√Y` and its exchange, then says "PY separately prints its signed-component and zero-locus exchange objects". SymPy's printed restriction (measurements 10090–10100, 10159–10161) is a squared velocity `2·GM/r`, with velocity `∓√2·√GM·r^(−1/2)` on each component and its implied `j_n`. Neither the record nor `R-S21-01` gives it. | **Confirmed.** The lookups show the SymPy objects in the filed measurements and in the transcript, and the §3 text that omits them. Decision-list item 4 requires each mechanism alone to be reported as printed, with its domain, implied `j_n` and cross-engine status. As written, the record reports one engine's form of this mechanism. |
| 3 | Claude 3 | **The engines' same-named `BRANCH_EXISTENCE` tags bind different conjunctions.** Wolfram emits `propagationGate = rayConstruction && Inactive[ForAll][…]` (engine 262–264, 275). SymPy emits `real_domain`, the along-ray condition (696, 720), and puts ray-root nonemptiness in `travel_domain` (697). On the forward-flow copy, SymPy prints `PY_S9B_F_FLOW_BRANCH_EXISTENCE: true`, and Wolfram prints a conditional over `rTurnFlyby`. The record's gate description (83–89) is accurate for Wolfram only, and the difference is not listed. | **Confirmed.** The lookups show both engines' emit lines and the two forward-flow payloads. The comparator leaves the row unformed ("unresolved coverage"). Under M1 a printed difference is a finding, and the record must state it rather than describe one composition for both engines. |
| 4 | Claude 4 | **S9b is not entered as a source on `R-S8-03`, the sign of the transverse stiffness.** Every A–C observable is gated on a real positive optical ratio `μ⊥/ρbr` (record 83–84) and prints `NOT_ESTABLISHED` outside. The pass adds S9b to `R-S8-01` (the stiffness form) and `R-S8-06` only (register 824). Its table row says "No new uniform-stiffness sign result or expanded R-S8-03 domain is inferred" (850). | **Confirmed.** `R-S8-03` asks for the sign of the physical transverse stiffness, and S9 and S10 are its sources on the same basis (318–324). Its own note records that an earlier pass "consolidated one and dropped the other" because `R-S8-01` never mentions sign (330–334). The pass's reason is a delivery test ("no new sign result"). The rest-on test asks only whether a kept S9b result rests on the requirement, and these results do. |

**Sound, both legs, each with its own scripts and stdout.** Both legs reran the generator, and it reproduces the
filed measurements file byte for byte (sha256 `26d7090f…`). Both legs' comparator enumerations match the record's
row classes, sample counts, `not_formed` reasons, symbolic classes and Appendices B and C, with nothing dropped or
merged. Round 0's repairs hold: the 14 SymPy deflection/radar pairs are identical and the 14 Wolfram pairs differ at
exactly the paths M13 tabulates; the K11 claim limit is quoted; the census reads pinned revisions; `F_MASS_SOLUTION`'s
fields are attributed by engine. The forward case is labelled conditional; the premises, D5's obligation and the
dependences are carried; no expected value, sign or density-link effect is stated.

**Not clear.** Findings 1–4 change what the record and the register claim or hand on. Finding 1 is in material
round 1 changed, and my round-0 wording caused it. Findings 2–4 are in material both round-0 legs passed. This is
the second repair, by the same author. The reviewed round-1 record is preserved before the repair overwrites it.

## What must be true after the repair
1. The record and the register count the deflection and the radar slope as separate tests by decision-list item 5's
   rule, applied to the printed objects, including both engines' printed γ-difference objects and their printed
   strata. Whether the differing Wolfram fields state the same condition stays open, as its own question.
2. Each single-mechanism restriction is reported as each engine prints it, with its domain, integration constants
   and their domains, implied `j_n` and cross-engine status, in the record and in `R-S21-01`.
3. The record states what each engine's `BRANCH_EXISTENCE` and `PATH_TRAVERSAL` tags bind, from the
   engine sources, and lists each printed difference between same-named rows as an open difference. No same-named
   tag is presented as one cross-engine object.
4. The S9b pass enters S9b as a source on each register entry whose requirement a kept S9b result rests on, by the
   register's rest-on test, including `R-S8-03`. The census and the "gaining a source" list follow.
5. Nothing else changes. The measurements file is regenerated by the final generator.
