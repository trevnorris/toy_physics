# Immutable transaction contract

This contract is mandatory for every semantic extraction or synthesis task.
The task-specific prompt may narrow it but may not weaken it.

## Required input packet

The caller supplies:

- `transaction_manifest`: the frozen prepared transaction manifest;
- `task`: a sealed, self-contained unit or derived-page work item containing
  its task kind, prompt paths, page ID/title/type, normalized source kind,
  source-unit shape/entry point/extractor version, complete member identities,
  allowed citations, related pages, stable-ID registry, desired lifecycle,
  and prior page/member history when applicable;
- `output_path`: the only staged page this task may write;
- `target_commit`: the full committed snapshot SHA;
- `refresh_date`: the date to use for memory metadata;
- `schema_snapshot`: an immutable schema-version-2 copy embedded in or listed
  by the transaction manifest;
- `base_page`, when the target is mixed, copied into the transaction with its
  generated-region ownership declared;
- all semantic source material, deterministic excerpts, identity-only member
  metadata, and staged memory dependencies explicitly listed for the task.
- `output_budget`, with applicable hard `max_words`, `max_key_statements`,
  `max_entries`, and/or `max_entry_words` ceilings;
- freshness metadata for every staged memory dependency: page hash,
  `generated_from_commit`, lifecycle, unit digest when applicable, policy
  digests, and whether it is fresh at the task target;
- a pre-run isolation contract naming the absent live workspace, read-only
  packet mount, sole writable output, cleared environment policy, and selected
  narrow runtime profile. The trusted launcher creates the signed/digested
  post-run attestation after the writer exits; finalization requires it.

Treat the frozen prompt bundle and prepared packet as the complete input
universe. If a required item or applicable budget is absent, stop and report
its name, the affected output, and why the omission prevents a faithful page.
Do not fill the gap from memory.

## Read boundary

The writer must run in a Gitless isolated filesystem that mounts only its
sealed packet read-only and its one declared output directory read-write. Live
repository roots must be absent, not merely forbidden by this prompt. If the
packet does not contain the launcher's pre-run isolation contract, stop before
reading semantic input. The writer does not create or self-certify the
post-run attestation.

Read only files explicitly listed for this task in `transaction_manifest`.
Those files must be inside the prepared transaction snapshot or sealed staged
dependencies.

Never read live paths under `research/` or `software/`, even when a displayed
citation names one. Never use `atlas/` or `graph/` as evidence. Do not use Git,
the raw filesystem, network search, or an unlisted memory page to expand the
packet. A visible repository path is a citation target mapped by the manifest,
not permission to open the live file.

Honor each member's `read_mode`:

- `semantic`: the supplied committed blob may be read in full.
- `excerpt`: only the deterministic supplied excerpt may be read. State when a
  broader conclusion cannot be supported from that excerpt.
- `identity_only`: only path, role, Git mode/type, blob ID, and size may be
  used. Never claim its contents, output values, tags, or success state.

Do not reconstruct an omitted identity-only transcript from an engine, result
record, prior knowledge, or another output. A semantic interpretive record may
say what it reports about that artifact, but attribute that interpretation to
the record and identify the underlying contents as unread.

## Evidence boundary

Summarize what the prepared sources state. Do not derive missing physics,
execute scripts, repair calculations, or strengthen an argument.

- Use `measured` only when the prepared packet supplies the exact recorded
  invocation, readable literal stdout or a tracked transcript containing the
  relevant operands and residual rather than a bare asserted zero, and a
  separate source that owns the interpretation. Cite each role. Code without
  its invocation is insufficient. A prose measurement report counts as output
  evidence only when it embeds the literal chain. If any required part is
  absent, use `provisional` or `open` as the source warrants and name it.
- Use `derived` only when a prepared source explicitly gives the derivation or
  proof within stated assumptions. Your own inference is not a derivation.
- Use `provisional` for a source assertion or synthesis whose required support
  is not present or is unclear.
- Use `open` only when a source explicitly leaves the object unresolved,
  deferred, or dependent on missing work.
- Use `disputed` only for retained, apparently incompatible positions with all
  live positions cited.

A script constructs, prints, compares, or guards named objects. Its `PASS`,
exit code, or asserted residual does not by itself establish a scientific
conclusion. Cross-engine agreement covers only the compared objects and is not
automatically independent evidence when premises, exports, or actions are
shared.

Call out missing commands, outputs, anchors, interpretations, premises, or
coverage. Do not silently omit a limitation that changes what may be claimed.

## Lifecycle, precedence, and conflict boundary

Set all generated wording to `memory_review: ai_draft`. Never emit
`memory_review: reviewed` or review provenance. Source review history does not
review the memory wording.

Do not select a position because its file, date, commit, step number, or prose
is newer, later, longer, or more confident. Do not infer scientific precedence
from a configured retrieval order or a successful run. Apply a preference only
when the prepared packet supplies one of these bases:

1. an explicit correction, revision, supersession, or retraction in a source;
2. a recorded reviewed adjudication of the same object and scope;
3. configured canonical lineage as a retrieval default only;
4. an explicit date tie-break between sources already shown comparable in
   role, evidence, scope, lineage, and adjudication.

Preserve the older or losing position. Record the exact replacement scope and
basis. Partial supersession applies to statements, not an entire capsule,
unless a source explicitly replaces the whole unit. When no basis resolves a
difference, retain all positions as a candidate conflict. First test object,
notation, sign convention, gauge/coordinate, approximation order, regime,
branch, source role, and historical/current scope.

## Citations and statement form

Every substantive statement must cite one or more original repository paths
mapped by the manifest. Never display a transaction or staging path. Cite the
responsible bundle member rather than only the unit or capsule.

Use the strongest supplied stable anchor in this order: exact TeX label; exact
Markdown heading or explicit anchor; distinctive step ID, emitted tag,
function, theorem, or result key. A line number may be supplemental but never
the only anchor. Do not invent an anchor. If none is present, cite the path with
the exact sentinel `anchor-unavailable`. A statement using that sentinel
cannot be `measured` or `derived`; keep it `provisional`, or `open` only when
the source explicitly says so.

Use this form for every substantive statement: scientific conclusions,
artifact-behavior claims, status judgments, limitations, open questions,
revision/supersession claims, preferences, conflict positions, and conflict
resolutions. Only neutral navigation and member metadata may omit it.

```markdown
### <stable-kebab-case-statement-id> — <short name>

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

<Conservative statement with object, scope, assumptions, and limits.>

Sources:

- `research/or/software/path` — `<exact stable anchor>`
```

Copy existing IDs from the task registry. Allocate a new ID only through the
task's deterministic rule, namespaced as `<page-id>--<semantic-key>`. Do not
derive identity from mutable wording or reuse an ID for a different statement.

## Output and ownership

Write only `output_path`. Produce Markdown only, with no wrapper fence or
process commentary. Copy manifest-owned values exactly: commit, unit ID,
digest, resolved members, roles, read modes, Git modes/types, blob IDs, sizes,
and output target. Do not fabricate or recompute them.

For a whole source capsule, replace the entire AI-generated body and use
`content_owner: ai_generated`; do not add generated-region markers. For a
mixed page, start from `base_page`, replace only the declared generated region,
and preserve every unmarked body byte. Never nest, rename, add, or remove
region markers. In frontmatter, change only machine-owned fields the task
explicitly delegates and copy their values from the manifest; preserve all
other fields. A substantive generated change leaves the page roll-up
`memory_review: ai_draft`.

Obey the task's hard output budget. When sources contain more material than
fits, prioritize benchmark/retrieval-critical statements and explicit status
boundaries; do not exceed the ceiling. Keep purpose first, omit long recap and
unneeded equations, and never duplicate a statement across sections. Use `not
supplied` or a precise limitation instead of speculation.

Before returning, verify that the page cites only mapped original paths, has no
unsupported evidence promotion, retains scoped historical/conflicting
positions, respects read modes and ownership, and contains no live-source or
legacy dependency.
