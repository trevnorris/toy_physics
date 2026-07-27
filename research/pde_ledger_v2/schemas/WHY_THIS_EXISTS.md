# Why the dimension-survey schemas and validator exist

> **Status: PARKED (survey-era).** The survey-first plan was abandoned in favor
> of the per-stage dimension rewrite. These committed schemas, validator, and
> fixtures remain as parked work; this file preserves their design rationale
> and verification state, not a live instruction to resume the survey.

The parked artifacts are the record schema, verification-output schema,
fixtures, and `validate_dimension_survey.py` in this directory.

## Verified state when parked

- Validator suite: **122 PASS / 0 FAIL**.
- **43/43 rules** were load-bearing under per-rule ablation.
- **105/105 minimal-violation rejects** were rejected.
- **82/91 call sites** were individually protected.
- **116/116 real parameter-register rows** were recordable.

## Survey-era design decisions

The numbering below intentionally preserves the original decision numbers.

1. **Schema-first split (USER).** The directive owns **semantics**; committed
   schemas plus the validator own **shape**. Prose cannot pin a structured
   record for 44 parallel agents.
2. **Two schema artifacts.** Keep a record schema and a verification-output
   schema.
3. **`ownership` and `uses` are orthogonal.** `declaration_required` derives
   from `ownership` alone. The earlier single `role` enum made
   `ASSERTED_TARGET` imply a declaration, which would have minted eight
   spurious canonical declarations from stage038's `EXPECTED_UNIT_STATE`
   alone.
4. **§3.11 is narrowed to dimension-valued bindings.** Producers such as
   `__mul__`, `dim_add`, and `dim_of` belong to §3.10 contract features. This
   makes the recovery-to-bindings correspondence a stated one-to-one
   relationship.
5. **`orders[].axes` contains interpreted names**, not source spellings.
6. **Anchor applicability is conditional.** Roughly 40 `PRESENT` stages have
   no §4.4 anchor and positively declare `NONE`; corpus anchors and ownership
   anchors are disjoint sets.
7. **§4.8 case (ii) requires full ownership coverage**, not a sample. This
   defeats shrinking the registry as well as zeroing it.
8. **Grok is the default tertiary; GLM is only for a further opinion** (USER).
11. **Semantic dimension adjudication is descoped from the survey.** Both
    channels record verbatim evidence plus a string-decidable status;
    monomial normalization and semantic adjudication belong to a later
    computed pass. This is deferred, not dropped.
12. **Pin `model: opus` on every survey-era agent dispatch** (USER). `fork`
    inherits by contract; `general-purpose` does not.
