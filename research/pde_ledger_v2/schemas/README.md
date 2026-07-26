# Dimension-survey schemas

These files implement the document shapes fixed by `SURVEY_DIRECTIVE.md` r18:

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
`dim_text_form`. **`basis_locus` is required exactly on `EXPLICIT_* +
POSITIONAL`, and forbidden on `EXPLICIT_* + NAMED_AXIS`, `EXPLICIT_* +
DIMENSIONLESS_CONSTRUCTOR`, and every non-`EXPLICIT_*` entry.**
`NAMED_AXIS` parses as a named-axis expression and requires a standalone source
or interpreted axis spelling (or exactly `1` for dimensionless); identifier
fragments such as the `time` in `run_time` do not qualify. Both the bare named
dimension and a fuller verbatim named expression are legal.

`POSITIONAL` must parse in its entirety as an ordered component sequence, binds
that sequence to the cited basis, and requires one component per declared axis.
For an arg-bearing constructor or factory, record only the narrow verbatim
parenthesized arguments—not the assignment or call wrapper. At
`scripts/ledger_stage005_sound_speed_light_ratio_sympy_audit.py:218`, record
`LENGTH = Dim(1, 0, 0)` as `dim_source_text: '(1, 0, 0)'`. The lowercase
factories have opposite orders: stage002's definition at
`scripts/ledger_stage002_matter_stress_force_assembly_sympy_audit.py:69` is
`def dim(l_power: int, t_power: int, m_power: int) -> Dim:`, so its `:199`
call uses `(L,T,M)` and records `'(1, 0, 0)'`; stage003's definition at
`scripts/ledger_stage003_transverse_photons_stray_longitudinal_sympy_audit.py:87`
is `def dim(m_power: int, l_power: int, t_power: int) -> Dim:`, so its `:866`
call uses `(M,L,T)` and records `'(1, -1, -2)'`. Read the local `def` line
before assigning axes.

`DIMENSIONLESS_CONSTRUCTOR` records the complete zero-argument call, for
example `dim_source_text: 'Dim()'` at
`scripts/ledger_stage004_gnls_action_dimensional_foundation_sympy_audit.py:143`.
It accepts only when the bare-name module-scope callee has exactly one lexical
binding event in the whole module: a preceding in-file `class`/`def` cited by
this report as its dimension machinery. Qualified, conditional, imported,
parameter-bound, and non-module forms reject as unsupported; an unmodelled
binding construct rejects fail-closed. Those diagnostics say the form is not
recordable and direct the surveyor to `UNDETERMINED` plus a finding rather than
inventing a binding. There is no length floor: honest `M`, `T⁻²`, and `L⁻¹`
values are valid.

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

The committed gate covers every §3.0 classification, both r18 ownership
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
