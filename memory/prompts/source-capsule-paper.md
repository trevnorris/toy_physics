# Paper source-capsule task

Apply `00-snapshot-contract.md`. This task accepts one configured unit whose
kind is `paper` or `componentized_paper` and writes its declared source capsule
as a whole-body `ai_generated` page.

## Inputs

Use only the unit record, resolved member metadata, and prepared member content
listed by the transaction task. For a componentized paper, use the manifest's
resolved bundle as authoritative; do not follow `\input`, bibliography, or
other paths beyond that set. Required machine values come from the manifest.
When the task supplies `editorial_focus`, treat each item as a required
retrieval surface, not as evidence: include it only to the extent the sealed
source supports it, with exact anchors and conservative status.

## Extraction

Identify only what is useful for later retrieval, within the task's capsule
budget (default ceiling: 1,800 words and 12 key statements):

1. the paper's named object, question, scope, regime, and exclusions;
2. its main explicit results, with assumptions and exact TeX labels;
3. derivations or theorems the paper actually supplies;
4. statements the paper itself calls provisional, open, disputed, corrected,
   or superseded;
5. explicit status boundaries and invalid-inference guards whose omission could
   make readers conflate distinct operations, sectors, observables, or regimes;
6. commands, audits, outputs, or verification records explicitly cited by the
   paper and present in the packet;
7. explicit relations to earlier or later treatments;
8. a few related topic IDs or script domains already supplied by the task.

The statement ceiling applies to every `### <stable statement ID> — ...` block
in the complete body, including `Purpose and scope`, computed evidence, open
questions, and revision relationships—not only blocks under `Key statements`.
When required retrieval surfaces would otherwise exceed the ceiling, keep the
configured editorial-focus surfaces and omit the least useful secondary result.

Do not treat an abstract, conclusion, or equation display as a derivation
unless the source explicitly supplies the reasoning. Do not calculate from an
equation. Preserve definitions and convention boundaries when omitting them
would change the meaning of a result. A paper date is descriptive metadata, not
a precedence decision.

Every material clause must be supported by at least one anchor listed in that
statement's `Sources` block. Split a sentence or add the exact supporting anchor
when one citation does not support every clause; do not attach a nearby section
that merely discusses the same broad topic.

Do not combine source-declared ontology, a phenomenological closure, and an
analytic identity under one `evidence=derived` status. Split statements when
their evidence classes differ. Preserve applicability conditions such as
excluded singularities, loop conditions, gauges, asymptotic assumptions, and
closure choices whenever they bound the statement.

## Output contract

Write the complete source-capsule template from schema version 2 to the staged
output path. Use:

- `type: source_capsule`, `content_owner: ai_generated`, and
  `memory_review: ai_draft`;
- `source_kind: paper`;
- page lifecycle copied from the transaction unit;
- `generated_from_commit`, extractor version, unit digest, and every resolved
  member identity copied exactly from the manifest;
- a sorted, duplicate-free `sources` list containing all direct members used.

The body has exactly these retrieval sections:

1. `Purpose and scope`
2. `Source-unit map`
3. `Key statements`
4. `Computed evidence represented by the source`
5. `Assumptions, exclusions, and open questions`
6. `Revision and supersession relationships`
7. `Related topics and scripts`

In `Source-unit map`, distinguish the entry point and component roles. In
computed evidence, write `not supplied in the prepared unit` when the packet
does not contain the command/output/interpretation chain; use that literal
sentinel, not a paraphrase. Do not create a
candidate conflict inside the capsule; report only the paper's own position
and flag possible cross-source comparison for topic/conflict synthesis.
Every substantive item in computed evidence, assumptions/open questions, and
revision sections uses the shared statement block; neutral member metadata may
remain a list.

For an unnumbered TeX heading, cite the exact source line as the heading text
(for example, heading `\section*{Acknowledgments}`), not only its rendered
title.

This last rule is structural, not optional: those three sections may not contain
uncited substantive prose before or after their statement blocks. A negative
claim that no formal supersession is established is itself a status judgment
and needs a provisional statement with citations to any correction, extension,
or earlier-treatment language that motivated the check.
