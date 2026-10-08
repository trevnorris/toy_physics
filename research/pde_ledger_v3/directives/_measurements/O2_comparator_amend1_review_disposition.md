# O2 comparator directive, amendment 1: review dispositions (orchestrator)

**Artifact:** amendment 1 to `directives/O2_comparator_build_directive.md`. It follows the user's finite-contract
decision (2026-10-07), made after four comparator build-review rounds (`O2_comparator_build_r3_review_disposition.md`).
It rewrites items 4 and 6 and adds one exclusion. The reviewed version is sha `4db4ce57…`, frozen at
`_scratch/s9b_build/O2_comparator_build_directive_amend1_reviewed_v0.md`. It is orchestrator-written and a changed
decision list, so it gets Codex + Grok (G2): one pass, findings verified, folded once.

**Legs (identical prompt `_scratch/s9b_build/o2_comparator_amend1_review_prompt.md`).** Both reported before
adjudication.
- **Codex:** NEEDS REVISION, 2 findings (`_scratch/s9b_build/o2_comparator_amend1_review_codex_final.txt`).
- **Grok:** NEEDS REVISION, 2 findings (`_scratch/s9b_build/o2_comparator_amend1_review_grok.txt`).
- **Agreed sound by both:** the amendment states no value, sign or outcome, freezes nothing, and keeps the user's
  exclusion. Grok also finds the exclusion compatible with spec §3.2, which gives no closed argument list.

Each verification is a mechanical lookup. Commands and literal output are in
`O2_comparator_amend1_review_disposition_lookups.md`.

| # | Finding (legs) | Disposition | What must be true after the fold |
|---|---|---|---|
| A1 | The printed limit is smaller than the blindness item 4 creates. A list of objects that recur across slots does not cover a swap of two once-each objects, a change in a multiplicity that stays at least two, a partial freeze inside one nested argument, or a change to argument content outside the inventories. The engines package arguments differently, so "slot" is not defined across them: SymPy passes `CarriedMomentumW`'s operands as top-level arguments (audit L85, L249–252), while Wolfram nests `Profiles`/`Velocity`/`Metric` inside one `section` argument (WL L23, L130–131, L202–205). A record could read an empty difference as agreement beyond what is compared. (Codex 1, Grok 1) | **ACCEPT.** Lookups A1. | No slot notion is used. The comparison covers exactly the declared inventories, and an empty difference means only that those inventories match. The output states what is not compared: placement among arguments, multiplicity, and content outside the inventories. For each paired OPEN occurrence and engine, it computes and prints the live objects occurring more than once and the count of parsed leaves outside the compared inventories. The user's exclusion stays. |
| A2 | The finite list misses controls for paths items 1–5 require. Missing: that a residual is a subtraction (an addition ablation passes); operands printed before the residual (item 3); a computed object against an OPEN action for the same role printed as a difference (item 4); a different head still joining (item 1); verdict tokens (item 5). Its "or" bullets can be met by one mutation. (Codex 2, Grok 2) | **ACCEPT.** Lookups A2. Item 6 had none of these; the "or" bullets are at v0:98–99. Grok's other uncompared-slot cases (a swap, a multiplicity change) are not compared under the user's scope, so they become declared limits (A1), not controls. | Each listed behaviour is a separate control. The contract includes subtraction (the residual reconstructs the left operand, and an addition ablation fails), operands before the residual, the computed-against-OPEN difference, a join that survives a different head, the two printed-limit quantities, and the absence of verdict tokens. |

**Fold.** A1 and A2 are folded once into items 4 and 6 and the exclusion. The folded version is sha `c3ebe39e…`
(lookups, "Fold applied"). Under G2 there is no second pass. The builder may resume.
