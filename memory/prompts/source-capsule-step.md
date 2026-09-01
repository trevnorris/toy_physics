# Step-family source-capsule task

Apply `00-snapshot-contract.md`. This task accepts one configured
`step_family` or `in_progress_step_family` and writes its declared source
capsule as a whole-body `ai_generated` page.

## Inputs

Use the prepared unit manifest to distinguish scope/specification, interpretive
step record, paper card, engines, premises, dependencies/exports, comparators,
controls, adjudications, and literal measurements. Do not infer those roles
from filenames. Never open an identity-only output.

## Extraction

Build a compact audit trail:

1. name the step's object, scope, branch/regime, assumptions, and exclusion;
2. state what the interpretive record says was computed and what remains open;
3. describe each engine by the objects it constructs and prints;
4. describe comparator coverage and controls separately from the engines;
5. map exact recorded commands to literal output paths/tags when supplied;
6. identify shared premises, actions, exports, or fixtures that limit claims of
   independence;
7. preserve historical/rejected substeps and explicit scoped supersession;
8. state the next gate or unresolved question exactly as the source does.

Never infer an output value from code or an identity-only transcript. Never
promote an exit code, assertion, or `PASS` label into a scientific conclusion.
For `measured`, cite the exact recorded invocation, readable literal stdout or
tracked transcript with operands and residual, and the separate
step/adjudication that interprets it. Code without its invocation is
insufficient. If a large output is identity-only, report its identity and the
semantic record's attributed interpretation separately; do not claim to have
verified its contents.

Paired engines are independent only to the scope explicitly documented. A
comparator measures agreement only over its actual compared keys/objects. An
in-progress unit remains open or provisional where the source has not recorded
completion, even if an engine or output artifact exists.

When an interpretive record names failed joins, unparsed categories, missing
route comparisons, or common-mode inputs, preserve those negative coverage
facts alongside positive comparator totals. When it explicitly marks an older
coverage rule, registry, harness, or provenance path as superseded, absent, or
historically unreproducible, include that scoped history even if the present
headline result is unchanged. Do not call a confirmed preregistered prediction
"corrected" merely because later controls narrow its interpretation.

## Output contract

Write the complete schema-version-2 source capsule with
`source_kind: step_family`, `content_owner: ai_generated`, and
`memory_review: ai_draft`. Copy lifecycle and all machine identity fields
exactly from the transaction.

Use the standard capsule sections. `Source-unit map` must list roles and read
modes. `Computed evidence represented by the source` must be a concise chain of
command/engine -> literal output -> comparator/control -> interpretation, with
each missing link named. `Revision and supersession relationships` must retain
the former statement and exact replacement scope; otherwise record the
difference for candidate-conflict synthesis.
Every substantive artifact-behavior claim, result, limitation, open gate,
historical position, and revision in any section uses the shared statement
block. Obey the capsule budget (default ceiling: 1,800 words and 12 key
statements).

This is a structural requirement: `Purpose and scope`, `Computed evidence
represented by the source`, `Assumptions, exclusions, and open questions`, and
`Revision and supersession relationships` may not contain uncited substantive
prose outside statement blocks. If there is no source-declared revision, leave
that section empty rather than asserting absence without a block. Prefer a
present function/tag/heading anchor over `anchor-unavailable`; in particular,
cite a Python `main` or named runner when the readable source supplies it.
