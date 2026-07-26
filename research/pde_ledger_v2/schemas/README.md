# Dimension-survey schemas

These files implement the document shapes fixed by `SURVEY_DIRECTIVE.md` r17:

- `dimension_survey_report_v1.yaml` governs one per-script survey report.
- `dimension_survey_verification_v1.yaml` governs one independent verifier output.

Both are YAML-encoded JSON Schema 2020-12 documents. Cross-field applicability,
item-ID, reference, ownership, cardinality, and verdict rules are named under
`x-directive-rules` and enforced by the companion validator.

Every quantity records `register_evidence`; every report records
`wl_evidence`. Neither channel carries an asserted dimension-equivalence
verdict. `TEXT_IDENTICAL` and `TEXT_DIFFERS` compare characters modulo
whitespace only. Semantic adjudication is deferred to directive §4.9.

Register evidence is bound to the committed
`notes/parameter_register.md`; a synthetic file with the same basename cannot
satisfy it. A `register_locus` selects one master-table row, never a range or
`R##` edge row, and `register_cell_verbatim` records the whole dimension cell,
including multi-parameter cells, `[X]=` prefixes, annotations, and
parentheticals. The validator does not parse it into one dimension. Lowercase
annotated `pending (...)` remains distinct from an uppercase `PENDING`
reduction-route flag. `ABSENT_FROM_REGISTER` is valid bare or with a reason;
`NOT_CHECKED` requires a reason. Neither requires a consultation receipt.

Each `wl_evidence.pairs` item has its own `item_id` and a `quantity_ref` to the
report's quantity inventory. `ABSENT_FROM_PY` alone uses
`NO_SURVEYED_QUANTITY`, with a required reason. Existing-side texts must occur
at loci in the report's Python or Mathematica identity path. The two text
statuses are checked by the same whitespace-only comparison as the register
channel.

An explicit `dim_source_text` must occur at its cited source and carry
`dim_text_form` plus a `basis_locus` from the report's declared basis.
`NAMED_AXIS` requires a standalone source or interpreted axis spelling (or
exactly `1` for dimensionless); identifier fragments such as the `time` in
`run_time` do not qualify. `POSITIONAL` binds an ordered component sequence to
the cited basis and requires one component per declared axis. There is no
length floor: honest `M`, `T⁻²`, and `L⁻¹` values are valid.

`qualified_scope` records a lexical scope, with an optional terminal binding
or mapping-key spelling. For stage042's `dims` keys, all of
`dimension_guard.B`, `dimension_guard.dims.B`, and
`dimension_guard.dims["B"]` normalize to lexical scope `dimension_guard`;
the `dims` container is not treated as a nested Python lexical scope. The same
normalization applies to `C`, `K`, and `cE`.

## Published validator command

Run this from the repository root, substituting either schema and the YAML file
to validate:

```sh
python3 research/pde_ledger_v2/schemas/validate_dimension_survey.py \
  --schema research/pde_ledger_v2/schemas/<schema>.yaml \
  --document <document>.yaml
```

The command exits `0` only for a valid document. On failure it exits non-zero
and prints one or more lines in this form:

```text
RULE_CODE $.path.to.node: explanation
```

The verification validator reads the referenced report YAML so it can check
its digest, classification, report-side item references, report cardinalities,
case-(i)/(ii) ownership coverage, and the recorded mechanical integrity result.
Case-(i) anchors are keyed by script/stage, enclosing lexical scope, binding
name, and branch state; absolute source lines are not anchor identity.
For `locus_integrity.status: CHECKED`, only values reached through
schema-declared locus types are resolved and checked; free-text strings are
never scanned by shape. The report's source paths, digests, and script line
count are rechecked in the same static leg. It never imports or executes an
audit script.
Only committed records under `schemas/examples/` may use
`status: SYNTHETIC_EXEMPT`.

## Committed example gate

The committed gate covers every §3.0 classification, both r17 ownership
coverage cases, all register and Mathematica evidence states, a valid verifier
`overall_verdict: FAIL`, and at least one minimal reject fixture for every
declared rule:

```sh
python3 research/pde_ledger_v2/schemas/validate_dimension_survey.py --test-examples
```

The expected outcomes and exact reject-rule assertions are in
`example_expectations.yaml`. Compact files under `examples/reject_fixtures/`
declare an accepted base document plus the minimal mutation that violates one
named rule; the published validator materializes these descriptors directly.
Removing the mutation restores the accepted base. Most evidence paths use a
`synthetic/` namespace; the dedicated locus fixture self-cites committed YAML
lines and fixed source-identity fixtures so both mechanical checks can run
without coupling examples to production audit files.

`dimension_survey_accept_fixture_v1` descriptors use the same base-plus-
mutation materialization for compact accepted constructions. The checked
case-(i) fixture uses deliberately non-production line numbers, so the gate
exercises lexical anchor identity without coupling to live audit-script lines.
