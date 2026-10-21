# S14a part 1 record decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_record_decision_list_review_p1_lookups.md`, beside this
file (generator `_scratch/s14a/gen/recdl_p1_lookups.sh`). "D-n" cites a section of that file.

## The pass

- **Artifact.** `directives/S14a_1_record_decision_list.md` (sha256 `a7fbcf25…`), with its setup lookups
  `S14a_1_record_setup_lookups.md`. Orchestrator-written, so the legs are Codex and Grok (G1). This is the single
  G2 pass; the list is folded once after it. Both files matched the baseline after the pass. The working tree was
  unchanged apart from the lookups file itself (D1).
- **Prompt and legs.** `_scratch/s14a/rec/recdl_review_prompt_p1.md`, sha256 `5cf11516…` (D1), identical for both,
  run in parallel. Codex (effort xhigh, 174,817 tokens) wrote `/var/projects/_review_scratch/s14a1-recdl-p1-codex/REVIEW.md`.
  Grok (`grok-4.7`, effort high) wrote `/var/projects/_review_scratch/s14a1-recdl-p1-grok/REVIEW.md`. Both exited 0
  (D1).
- **Verdicts.** Codex: CLEAR. Grok: FINDINGS, 1 (medium) (D1).

## Findings

**K1. Item 9 would leave the accepted S12 pass saying it awaits review.** (Grok 1.) **CONFIRMED.**
- **The list.** Item 9: "The stale status text and the known items J1 and J2 are updated as the S12 record's round-3
  disposition states them (L8). Nothing else in the S12 record or in the register's S12 pass changes." (D2)
- **The source.** The S12 round-3 disposition names "the record header, the register's status and S12 pass heading,
  and `STATUS.md`" (D2).
- **The gap.** The register's S12 pass heading (l.1298) is followed by its source sentence, which ends "this is the
  round-3 repair of reviewed baseline `fa5d1396`, awaiting the narrow check of H1/H2 (M19)" (l.1300–1302, D2). That
  sentence is not the heading. Item 9's last sentence bars changing it, so the accepted pass would keep saying it
  awaits a check. The record header's matching clause (l.5) is inside the header the disposition names (D2).
- **What must be true.** Every line in the S12 record, the register (its status and its S12 pass) and `STATUS.md`
  that describes S12's record or register pass as written for, under or awaiting review states the acceptance at
  `ac909aaf`, whether or not it is a heading.

**K2. Item 9 does not say whether J1 and J2 are repaired.** (Orchestrator, from Grok's reading.) **CONFIRMED.**
- **The list.** "the known items J1 and J2 are updated as the S12 record's round-3 disposition states them" (D2).
- **The reading.** Grok reads this as J1 and J2 "carried as known items" (D3). The disposition states them as defects
  and says the next record list "picks up J1 and J2" (D2), but it gives no "What must be true" for either (D3). So the
  list leaves the author to choose between repairing them and carrying them again.
- **What was intended.** Both are repaired: J1 in the S12 record's and the S12 pass's account of edits, J2 in the
  "18 → 18" sentence in both places (record l.340–343, register l.1315–1317, D3).
- **What must be true.** The account of edits names the two status-only edits J1 lists, and the "18 → 18" sentence
  for the round-3 repair cites the section that holds its count, in both the record and the S12 pass.

## Checked and found sound by both legs

Items 1–8 and 10: the object and its standing, cross-engine evidence and coverage, the instruments' limits, the
inventories, premises and OPEN items, the register pass, and the account of edits. Grok's census of the comparator's
output found the ten bindings whose exact residual is not literally `0` or `False` to be the set the comparator's
round 1 carries to the record. Neither leg found an outcome stated in the list.

## Outcome

K1 and K2 are confirmed. Both are in item 9, and neither changes what the engines or the comparison may be cited
for. The reviewed list is preserved before the fold. The one fold rewrites item 9's first bullet to the two "What must
be true" statements above. Then the record author starts.
