# Software and result-family source-capsule task

Apply `00-snapshot-contract.md`. This task accepts one configured
`software_project`, `software_result_family`, or `verification_project` and
writes its declared source capsule as a whole-body `ai_generated` page.

## Inputs

Use member roles from the transaction manifest. Keep architecture/workflow
documentation, entry points, engines, comparators, controls, result reports,
adjudications, and measurements distinct. Do not discover helpers, generated
artifacts, or output directories outside the resolved unit.

## Extraction

Record only retrieval-useful information:

1. the project's or result family's purpose and declared scope;
2. meaningful executable entry points and the named objects each computes,
   builds, validates, compares, renders, or reports;
3. inputs, arguments, fixtures, imported exports, and output paths/tags;
4. guards, invariants, comparator coverage, and controls;
5. the result/verdict source that owns any scientific interpretation;
6. architecture or commands that remain current despite a scoped status
   revision;
7. open limitations, known failures, and explicit supersession.

Source code supports a description of program behavior but not an unexecuted
scientific result. A result report may be summarized only at its stated scope.
Use `measured` for a substantive result only when the packet supplies the
exact recorded invocation, readable literal stdout or tracked transcript with
the relevant operands and residual, and the separate interpretive
report/adjudication. Code without its invocation is insufficient. A prose
measurement report qualifies only when it embeds that literal chain. For
identity-only result data, cite its identity but do not state unobserved values.

Preserve any source-declared authority hierarchy between a readable diagnostic
report and an authoritative packet. If the authoritative packet is
identity-only, say so, cite that exact identity with `anchor-unavailable`, and
keep the diagnostic conclusion provisional; do not imply that the packet's
contents or success state were inspected.

Do not let a verdict supersede unrelated architecture, interface, or command
documentation. Apply only the exact source-declared replacement scope. Keep
visualizations separate from verification claims.

## Output contract

Write the complete schema-version-2 source capsule with
`content_owner: ai_generated`, `memory_review: ai_draft`, and a normalized
`source_kind` copied from the task or mapped as follows:

- `verification_project` -> `verification`
- `software_project` -> `software_project`
- `software_result_family` -> `result`

Copy lifecycle, commit, extractor version, unit digest, and all resolved member
identities exactly. Use the standard capsule sections. In the source-unit map,
name meaningful entry points and supporting members rather than flattening the
bundle. In computed evidence, give the exact available chain and say `not
supplied` for missing invocation, output, comparator, or interpretation.
Every substantive behavior, result, limitation, open status, and revision in
any section uses the shared statement block. Obey the capsule budget (default
ceiling: 1,800 words and 12 key statements).

This is a structural requirement: `Purpose and scope`, `Computed evidence
represented by the source`, `Assumptions, exclusions, and open questions`, and
`Revision and supersession relationships` may not contain uncited substantive
prose outside statement blocks. If there is no source-declared revision, leave
that section empty rather than asserting absence without a block. A statement
about identity-only artifacts must stay inside a cited provisional block and
must not imply their unread contents were checked.
