# O2 comparator build r0: review dispositions (orchestrator)

**Artifacts:** the comparator, its tests, the builder's full guarded run and the builder report, at
`_scratch/s9b_build/o2_comparator_build_review_baseline_r0.sha256`. The baseline check prints `OK` for all eight
(lookups). Codex wrote them to `directives/O2_comparator_build_directive.md` (folded `fdaa575b`).

**Legs (Codex-written → fresh Claude + Grok, identical prompt `o2_comparator_build_review_prompt_r0.md`).** This
is comparator build review round 1. Both legs reported before any adjudication.
- **Fresh Claude (opus): NEEDS REVISION**, 4 findings (`_scratch/s9b_build/o2_comparator_build_review_r0_claude.md`,
  condensed by me; leg scratch `/tmp/o2cmp_review_leg_claude_r0/`).
- **Grok: NEEDS REVISION**, 2 findings (`_scratch/s9b_build/o2_comparator_build_review_r0_grok.txt`).

**Agreed sound by both:**
- Coverage: every parsed leaf is accounted for, and every object is joined or unjoined with a specific reason.
- Both tables are injective.
- Each leg independently recomputed exact rows from the raw transcripts and matched the comparator on each one.
- No applied function is collapsed, operands are printed before residuals, and text or name pairs are never zero.
- No verdict tokens; exit 0; no ablation triple joined; the 26 tests pass unmodified.

Claude also reproduced `comparison.jsonl` and `accounting.jsonl` byte for byte.

Each verification below is a mechanical lookup. Commands and literal output are in
`O2_comparator_build_r0_review_disposition_lookups.md`.

| # | Finding (legs) | Disposition | What must be true after repair |
|---|---|---|---|
| R1 | Join row `mass_loss_orientation` (comparator L424) pairs the SymPy `outward_loss` (`jn`, audit L377–378) with the Wolfram `RHS` (`massRHS = -profiles["j_n"]`, WL L244). These are different objects: one is the outward loss, the other the right-hand side of the mass law. The equation itself is already joined at `mass_equation` (L384). The row manufactures the only nonzero exact residual in the output. (Claude 1, Grok 1) | **ACCEPT.** Lookups R1. | Each join row pairs occurrences of the same spec §9 object in the same role. An occurrence with no same-role counterpart in the other stream is unjoined with that reason. No printed residual arises from pairing two different objects. |
| R2 | `signed_actions` (L672–691) resets the sign to +1 inside every applied head. A minus on a held aggregate such as WL's `-Inactive[Total][Inactive[Map][… OpenAction[NativeToCoordinateDensity, …]]]` therefore never reaches the OPEN actions inside it. Both engines enter the face support at −1 (PY audit L259–260; WL L231), yet the census prints PY −1 against WL +1. The census manufactures a sign difference. (Claude 2, Grok 2) | **ACCEPT.** Lookups R2. | The printed orientation of each balance term is the sign that term carries in its engine's balance, including terms wrapped in held aggregates. Argument occurrences inside an OPEN action stay distinguished from signed balance terms. No orientation difference is printed that the two engines' entries do not carry. |
| R3 | The controls pass on operand text alone. `changes()` (test L35–38) compares the whole `comparison` dict, which includes `py_value`/`wl_value` (L756). Any operand mutation therefore "moves" it whether or not the residual or census moves. Claude's comparator FORM ablations leave all 26 tests green: A1 collapse `X(args)→X`, A3 sign never flips, A4 transpose layout dropped. A5 (residual hard-wired to 0) is caught only by the repoint test. Grok's `left-right → left+right` copy also leaves the suite green. (Claude 3; Grok observed the same suite result) | **ACCEPT.** Lookups R3. The legs **disagree on reading**: Grok took the green suite as admissible because directive item 6 asks each mutation to "change the printed output". My item 6 wording allowed that reading. The intent, recorded in U3 (`O2_comparator_directive_review_disposition.md`), is a control that fails when its defect goes undetected. Operand text always changes under an operand mutation, so the test cannot detect a comparison defect: a check that audits its own input. The census defect R2 went unflagged by this suite, which shows the gap is load-bearing. | Each control fails when the comparator's own comparison misses its defect. The assertion is on what the comparison produces (the outcome and residual, and the structure census for OPEN controls), never on operand text alone, and still names no value. Controls exist for argument content on the exact-subtraction path, held-aggregate orientation and the transpose layout. Each listed FORM ablation of the comparator's own logic (collapse, sign never flipping, layout dropped, residual not formed by subtraction) makes at least one test fail. |
| R4 | The name table omits four OPEN-operand pairs that the comparator's own joins and unjoined reasons treat as the same objects: `material_action_compatibility` ↔ `UnresolvedStressInertiaNormalIdentifications` (join L430), `energy_accounting_overlap` ↔ `UnresolvedEnergyOccurrenceIdentifications` (L431), `S12_reaction_system` ↔ `OPENReactionSystem` (L804), `face_support_partition` ↔ `UnresolvedSupportPartition` (L808). The census prints each as a name only one engine has. (Claude 4) | **ACCEPT.** Lookups R4: none of the eight names occurs in the comparator; all occur in their engine sources. | The name table, the join rows and the unjoined reasons are consistent. Every OPEN operand that a join or reason identifies across the engines is bound in the name table with its spec citation. If the spec does not identify them, neither the joins nor the reasons do. No census prints a one-engine name for an operand the comparator elsewhere identifies. |

**My own process defect (R3):** directive item 6 said "fails unless its mutation changes the printed output". The
printed output includes the operands, which permitted a vacuous control. The repair brief states the intent
instead of that phrase.

## Outcome

**NOT ACCEPTED.** Four findings are accepted, and all bear on what a record may claim from the output: R1 and R2
each manufacture a disagreement, and R3 and R4 leave the instrument unable to catch, or inconsistent about, the
same class of error. The reviewed baseline is preserved in a commit before the repair overwrites it (G4).
Repair round 1 goes to the same builder; it is the first repair, so no author change is called for. Then fresh
Claude + Grok legs review, until clear.
