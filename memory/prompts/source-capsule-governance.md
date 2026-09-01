# Research-governance source-capsule task

Apply `00-snapshot-contract.md`. This task accepts one configured
`research_governance` unit and writes its declared source capsule as a
whole-body `ai_generated` page.

## Inputs

Use only the prepared governance members and metadata assigned to the unit.
Treat charters, plans, requirement lists, defect registers, and handoffs as
different roles even when bundled. Do not read the governed papers or steps
unless the transaction explicitly includes them in this task.

## Extraction

Extract a compact operational map:

1. program scope and the objects the governance documents authorize;
2. current frontier, stage/step states, and gates exactly as declared;
3. status firewalls and explicitly forbidden inferences;
4. substrate or acceptance requirements, without supplying expected results;
5. defects, blockers, deferred work, and ownership boundaries;
6. explicit revision/supersession relationships among governance statements;
7. links to configured paper, step, topic, and script IDs supplied by the task.

A plan or charter is evidence of program intent and declared status, not
evidence that the proposed physics was derived or measured. A requirement that
an engine compute an object does not mean it computed it. A checked box, later
step number, or `PASS` string does not close a gate unless an interpretive
governance source explicitly records the closure and its evidence. Retain
separate current, completed, rejected, and open scopes.

Do not convert a defect hypothesis into a confirmed defect, or a future step
into an implemented result. If sources disagree about status and provide no
explicit resolution, preserve both statements and flag a candidate conflict.

## Output contract

Write the complete schema-version-2 source capsule with
`source_kind: governance`, `content_owner: ai_generated`, and all generated
wording `memory_review: ai_draft`. Copy lifecycle, target commit, unit digest,
extractor version, and complete member identities from the manifest.

Use the standard capsule sections. Within `Key statements`, prefer a small set
of stable statements covering current frontier, open gates, and status
firewalls. In `Computed evidence represented by the source`, distinguish
evidence cited by governance from work merely requested; use `not supplied`
where the packet does not carry the command/output/interpretation chain.
Every substantive status, limitation, defect, gate, and revision in any
section uses the shared statement block. Obey the capsule budget (default
ceiling: 1,800 words and 12 key statements).
