# S11 SymPy export — decision list

**Codex-authored, 2026-08-11.** Build the SymPy engine from the physics in
`git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`; this list decides only its
merged export and the shared namespace it enters. The three earlier PY lists and `F8` are blocked records.

## Decisions

**D1 · A key is an opaque row locator, not an object name.** Consumers may not infer physics, scope, or
identity from it. Existing keys remain unchanged as locators; no upstream row is renamed or re-pointed.
This explicitly **SUPERSEDES F1**: repeated attempts to put identity in the key made wrong bindings pass
their controls. F2–F7 remain in force, with same-object meetings broadened from key collisions to the
claim join below.

**D2 · Every row carries a typed `object_claim`.** A copied key, tag, slug, package label, action
descriptor, or payload is not a claim. For a tag-derived row the claim is the specification's named
mathematical role, its mathematical scope, and references to exactly the named defining inputs and
premises on which that role depends. Those references are to other claims, never to their computed
values; producer, engine, package, and step labels are excluded unless the physics makes one part of the
object. For a traversal-manufactured symbol row the claim is the exact SymPy `Symbol`, assumptions
included. A claim is created with the live object, never assigned later by a rename map.

Claim-tree equality is total and independent of payload type: equal trees mean the same claimed object;
unequal trees mean different claimed objects; a missing, malformed, or specification-undetermined tree is
`UNDECIDED`. `UNDECIDED` is reported and may neither merge, corroborate, overwrite, nor publish. Equal
payloads never establish identity, and observed payload movement never defines it.

**D3 · The merge joins claims across the whole import and the whole candidate export.** One claim has one
merged row, regardless of locator. Every new derivation of an imported claim is compared with that row.
Different claims remain different rows even when their payloads or human-readable names coincide. A
locator collision between different claims is emitted as a finding and cannot overwrite either row.

**D4 · The row contract is typed but extension-preserving.** Every row has `object_claim`, `display`,
`value`, `value_kind`, `class`, and `step`. Every inter-row reference carries both its locator and the
target claim; this includes a `dimension_key` and its companion `dimension_claim`. A re-derived row also has
`comparisons`; `corroborated_steps` may be present only when the corresponding comparison evidence is
present. Unknown imported fields are retained unchanged. Missing mandatory fields, a reference without
its target claim, an unsupported field shape, or a false corroboration is a schema finding: emit the row,
the contract, and their discrepancy, then do not publish. No default, field inference, or repair is
allowed.

**D5 · Identity and agreement are separate decisions.** Rows with the same claim compare their exact
reconstructed values under the claim's premises and domain. Each `value_kind` has a structure-preserving
residual; subtraction is used only for objects that admit it, while sequences, matrices, relations, sets,
predicates, records, symbols, and text retain their own type and structure. The result is `AGREE`,
`DISAGREE`, or `UNDECIDED`; agreement requires the residual to establish no difference, and lack of a
decision is never agreement. Each `comparisons` entry has `imported_claim`, `derived_claim`,
`imported_value`, `derived_value`, `premises`, `residual`, `status`, and the two producing steps. Thus a
consumer of only the merged export can recompute it, as F3 requires.

**D6 · Imported physics is append-only.** Before computation, independently reconstruct the complete
import. In the merged export every imported field and reference is unchanged. A verified re-derivation
may append its D5 comparison record; it does not replace the imported value or redirect any reference.
The preservation guard compares the imported reconstruction with the imported-row projection of the
merged reconstruction, over every field, using the independently formed claim join rather than the
writer's touched-key set. Unsupported predecessor corroboration prevents publication rather than being
deleted or forwarded.

**D7 · Export eligibility is total over both row populations.** Tag-derived rows are every object the
physics specification requires `MAIN` to emit at every declared `MAIN` dimension, including an object
whose tag name the builder must supply and report, plus the corresponding premise inventories and the
observed run record. Control-package rows, solver-internal rows, local-name inventories, and standalone
copies of the previous-step comparison are not exported; re-derivation evidence lives in `comparisons`.
Symbol rows cover every distinct SymPy `Symbol` atom needed to reconstruct or interpret any
physics-bearing field of those rows, including claims, premises, references, and comparison evidence.
Same-spelled symbols with different assumptions remain distinct claims and rows.

**D8 · Omission is checked from an independent required population.** The physics specification and the
declared `MAIN` sweep supply the required claim population, independently of the emitter and writer;
required ordering objects supply only the indices whose scopes they define. Compare, with multiplicity,
required claims against eligible emitted claim–value associations, those associations against candidate
rows, and candidate rows against reconstructed written rows, in both directions. Each comparison emits
its two populations and residual before guarding. The symbol population is checked the same way. An
object absent before emission, assigned another object's value, omitted by the writer, duplicated, or
spuriously added therefore remains visible.

**D9 · Publication means current, complete, and reconstructable.** Publication is permitted only after
the observed run record accounts for every declared `MAIN` cell; all claims and references are decided;
all same-claim comparisons have status `AGREE`; D6 and D8 have completed; and the whole written artifact
round-trips field-for-field. Once a run begins, no earlier S11 artifact may be accepted as that run's
output. If the run is incomplete or any publication condition is unresolved, no artifact representing
that run is importable as complete; refusing to overwrite an older file is insufficient.

**D10 · S10 must be regenerated before S11 can publish.** Its current rows lack the claim and F3 evidence
contract. Each producer adds claims only to rows it derives; inherited rows must arrive claimed from their
producer, so the regeneration reaches upstream as far as the imported population requires and no
downstream step guesses identity. Every existing key remains a locator and every existing key-to-row
field and reference is exactly reconstructed, adding only claims, reference claims, and supported
comparison evidence. The regenerated claim–value associations are checked against the live producer
objects and the original rows, with both operands and residuals. There is no migration map or destination
change to re-point. Until that regeneration has passed its own two legs, S11 may compute and report but
may not publish a merged export.

**D11 · Every guard emits both operands and its residual, then guards.** No guard audits an artifact
against its own input, and no asserted-only residual is evidence. A physics disagreement or undecided
comparison is reported without changing the computation to remove it.

**D12 · REPORTING IS SUCCESS.** An ambiguity, missing claim dependency, unsupported comparison kind,
schema incompatibility, incomplete run, disagreement, or undecided result is emitted and reported; it is
never closed by an invented name, identity, default, or comparison. Such a finding can withhold
publication without converting the physics finding into an operational failure.
