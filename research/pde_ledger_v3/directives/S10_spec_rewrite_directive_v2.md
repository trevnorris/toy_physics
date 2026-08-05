# Authoring contract for the next S10 shared specification

## Assignment boundary

Replace the contents of
`research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md` with a new computation contract. The two CAS
implementations will be built from that artifact without consulting one another, so the contract itself
must carry the redundancy and provenance needed to expose a bad implementation.

Edit that file only. Do not commit. Do not edit either engine, either committed output, or anything under
`reduction/`, `steps/`, or `paper/`. The engine and harness changes named below are obligations for their
future builders, not changes authorised in this task.

This assignment has no dependency on `S10_spec_rewrite_directive.md`; do not consult or reuse it. Read
these sources instead:

- `S10_SHARED_PHYSICS.md`, to inventory the supplied inputs and requested computations;
- `S10_SPEC_CONSISTENCY_FINDINGS.md`, to locate defects only;
- both S10 engine sources, as read-only text;
- both committed S10 `.out` files, reading only the tag name to the left of the first colon;
- `reduction/checks_S10.yaml` and `S10_py_repair3_directive.md`, for the live naming dependencies.

The findings file is unsanitised. It contains computed ratios, counts, and a residual judgement. Keep
those facts in working memory only: do not copy, paraphrase, negate, or warn about any of them in the new
specification. A prohibition can disclose a result just as readily as an assertion.

Do not execute either engine. In particular, do not invoke Mathematica, `wolframscript`, or a Wolfram
kernel. Wrap every script or shell pipeline in `timeout 600`; never increase that limit.

## The artifact is a closed data contract

The legacy `.md` suffix remains, but the replacement file must be one canonical JSON document. It must
contain no Markdown framing, comments, narrative text, examples, history, rationale, review findings, or
free-form string fields. This is deliberate: a second sentence cannot contradict an equation when the
format has no place for a second sentence.

Use this closed model:

- Top-level keys are exactly `contract_version`, `symbols`, `nodes`, `packages`, `q7`, `classifiers`,
  `emissions`, `quantity_registry`, and `report_contract`.
- Every object rejects unknown keys. Decode JSON with duplicate-key rejection.
- Every definitional or computational object has one globally unique `id`; an `id` has exactly one
  owner. All reuse is by `["ref","<id>"]`. Reject dangling references and definition cycles.
- Expressions are recursive prefix-array ASTs, never prose or CAS source strings. The only scalar leaves
  are `["integer",n]`, a reduced `["rational",p,q]`, `["symbol","<registry-id>"]`,
  `["enum","<allowed-token>"]`, and `["ref","<id>"]`. A composite is
  `["<operator>",<argument>...]`. The authoring-time validator owns closed operator and enum allowlists;
  the target may use only their members.
- Human labels, descriptions, notes, expectations, measured values, and result assertions are not schema
  fields. String values are identifiers, tag tokens, symbol names, or closed enums and must match their
  declared lexical grammar.
- A computation node records only its operation, inputs, scope, dependencies, outputs, and status domain.
  It never records an answer, an expected answer, a comparison outcome, or a success verdict.
- Numeric literals may occur only in supplied premises, supplied definitions, index/range declarations,
  and algorithm parameters. An emission node contains a reference and a payload type, never a literal
  result.

Define the JSON schema in the authoring-time validator, not inside the artifact being validated. A schema
stored beside the data would let the author legalise their own contradiction. Include the complete
validator as an inline command in the report; do not create another repository file.

## Immutable physics inputs

Encode the field, symbolic dimension, joint premises, real plane-wave ansatz, phase average, and
structural compute-from-action rule as typed nodes. Preserve the requested Q1--Q8 computations, but
translate them into definitions and operations rather than copying existing prose. The computation graph
must cover Euler--Lagrange equations, both matrix routes, the spectrum and loci, rank-based mode data,
wavevector scaling, dimension analysis, the three-dimensional reference calculation, and rank-drop
strata.

Use these as the author-owned spatial formulas:

```text
G_ij       := partial_i u_j
S_anti     := (1/2) sum_{i=1..D} sum_{j=1..D} (G_ij - G_ji)^2
S_full     := sum_{i=1..D} sum_{j=1..D} G_ij^2
S_div      := (sum_{i=1..D} G_ii)^2
T_iso      := (rho_br/2) sum_{j=1..D} (partial_t u_j)^2
T_aniso    := (rho_br/2) (s_rho (partial_t u_1)^2
                          + sum_{j=2..D} (partial_t u_j)^2)
```

The package table must contain the complete action, the complete spatial weight, and the density report
for every row. These are supplied definitions, not computed results. The abbreviations in this table are
for readability in this directive; in the target AST, instantiate the action's spatial expression and the
density-report expression under separate owners. The action node must not reference the weight-report or
density-report node.

| package | dimensions | action | spatial weight | spatial density |
|---|---|---|---|---|
| `MAIN` | `2,3,4,5` | `T_iso - (mu_R/2) S_anti` | `-mu_R/2` | `S_anti` |
| `XFORM_FULLGRAD` | `3,4` | `T_iso - (mu_R/2) S_full` | `-mu_R/2` | `S_full` |
| `XFORM_DIVONLY` | `3,4` | `T_iso - (mu_R/2) S_div` | `-mu_R/2` | `S_div` |
| `XFORM_SIGNFLIP` | `3,4` | `T_iso + (mu_R/2) S_anti` | `mu_R/2` | `S_anti` |
| `XFORM_ANISO` | `3,4` | `T_aniso - (mu_R/2) S_anti` | `-mu_R/2` | `S_anti` |
| `XCOEF_SCALE` | `3` | `T_iso - (s mu_R/2) S_anti` | `-s mu_R/2` | `S_anti` |

Keep the positivity, reality, nonzero-wavevector, positive-integer-dimension, rest-background,
non-dissipative, quadratic-response, and declared-dimensionless inputs as premises. Do not add a premise
to force a solver to decide.

Pin Route A rather than asking a reader to remove an unspecified common factor. With `E(a)` the vector of
Euler--Lagrange expressions after the real ansatz is substituted, define `C(a)` by exact coefficient
extraction of the supplied phase factor, then define

```text
M_A := -(1/2) Jacobian_a(C)
M_B := Hessian_a(phase_average(L))
```

Emit both matrices and their oriented residual as referenced computation objects. Do not attach a claim
about the residual's value.

Encode the dimension inputs explicitly:

```text
[energy]    := (2,-2,1)
[u]         := (1,0,0)
[L]         := [energy] + (-D,0,0)
[partial_i] := (-1,0,0)
[partial_t] := (0,-1,0)
```

For package `p`, let `C_p` be the set of atomic symbols in its full action whose declared role is inertial
or stiffness coefficient, excluding every symbol declared dimensionless. Define the scalar unknown set

```text
U_p := { dim_slot(c,slot) : c in C_p, slot in {length,time,mass} }
N_unknown(p) := cardinality(U_p)
```

Form the dimension equations from the actual action AST. Emit the equations, `U_p`, its cardinality, the
independent-equation count, their difference, and the solve status. For the linear dimension system
`A_p x = b_p`, define `N_equations(p) := rank(A_p)`, define consistency by
`rank(A_p) = rank([A_p|b_p])`, and define the requested difference as
`N_equations(p) - N_unknown(p)`. Classify over/exact/under-determination only for a consistent system; use
a separate inconsistent-system status otherwise. A coefficient expression is not an atomic coefficient
symbol, and a declared-dimensionless factor contributes no unknown slot.

## Q7 certificate: three independent reports

Adopt a new action-constructor interface. Each future engine must return
`{lagrangian, spatial_weight, spatial_density}` for every package. This is an engine requirement caused by
this specification rewrite.

The three returned objects must be independent reports:

- build `lagrangian` directly from the package's full-action AST;
- build `spatial_weight` directly from the package's weight node;
- build `spatial_density` directly from the package's density node;
- the three routes may share primitive symbols and derivative leaves, but no route may reuse, divide,
  factor, inspect, or solve another returned composite object;
- none may be defined through the Levi-Civita reference.

This duplicates a small declaration deliberately. It is the redundancy that makes the extraction check
capable of exposing a one-route error. If an implementation instead assembles `lagrangian` by multiplying
the two objects it also returns, it has removed the independent route and does not conform.

The target contract's `q7` object must canonicalise exactly to the following author-owned oracle. The
rewrite author may not edit, grade, or select among these decisions.

<!-- AUTHOR_Q7_ORACLE_BEGIN -->
```json
{
  "purpose": {
    "check": "form_and_normalisation",
    "coefficient_role": "wave_speed",
    "dimension_reason": "ordinary_cross_product_domain",
    "evidence_exclusion": "mode_count"
  },
  "applicability": {
    "predicate": ["eq", ["ref", "D"], ["integer", 3]],
    "emit_same_quantity_family_at_every_declared_dimension": true,
    "inactive_status": "not_applicable"
  },
  "constructor_report": {
    "fields": ["lagrangian", "spatial_weight", "spatial_density"],
    "composite_routes": ["full_action_ast", "weight_node", "density_node"],
    "shared_dependencies": ["primitive_symbols", "derivative_leaves"],
    "forbidden_cross_route_derivation": true,
    "levi_civita_dependency": false
  },
  "operations": [
    {
      "id": "q7_static_lagrangian",
      "operator": "substitute",
      "arguments": [
        ["ref", "package.lagrangian"],
        ["all", ["enum", "velocity_components"], ["integer", 0]],
        ["map", ["enum", "spatial_derivatives"], ["enum", "independent_g_ij"]]
      ]
    },
    {
      "id": "q7_weight_zero_region",
      "operator": "real_satisfiability_status",
      "arguments": [
        ["and", ["ref", "P_eff"], ["eq", ["ref", "package.spatial_weight"], ["integer", 0]]]
      ]
    },
    {
      "id": "q7_extracted_density",
      "operator": "guarded_exact_quotient",
      "guard": ["eq", ["ref", "q7_weight_zero_region"], ["enum", "EMPTY"]],
      "arguments": [
        ["ref", "q7_static_lagrangian"],
        ["ref", "package.spatial_weight"]
      ],
      "pre_guard_simplification": false
    },
    {
      "id": "q7_extraction_residual",
      "operator": "subtract",
      "arguments": [
        ["ref", "q7_extracted_density"],
        ["ref", "package.spatial_density"]
      ]
    },
    {
      "id": "q7_levi_civita_vector",
      "operator": "indexed_sum",
      "arguments": [
        ["levi_civita", ["symbol", "i"], ["symbol", "j"], ["symbol", "k"]],
        ["ref", "g_jk"],
        ["range", ["symbol", "j"], ["integer", 1], ["integer", 3]],
        ["range", ["symbol", "k"], ["integer", 1], ["integer", 3]]
      ],
      "expanded_literal_allowed": false
    },
    {
      "id": "q7_levi_civita_norm",
      "operator": "dot",
      "arguments": [
        ["ref", "q7_levi_civita_vector"],
        ["ref", "q7_levi_civita_vector"]
      ]
    },
    {
      "id": "q7_reference_residual",
      "operator": "subtract",
      "arguments": [
        ["ref", "package.spatial_density"],
        ["ref", "q7_levi_civita_norm"]
      ]
    }
  ],
  "forbidden_residual_shape": [
    "subtract",
    ["ref", "dividend"],
    ["multiply", ["exact_quotient", ["ref", "dividend"], ["ref", "divisor"]], ["ref", "divisor"]]
  ],
  "residual_orientation": "package_object_minus_reference"
}
```
<!-- AUTHOR_Q7_ORACLE_END -->

Zero every velocity component before extracting the static object; do not retain or discard any other
term. Map every spatial derivative to an independent `g_ij`. Query the zero region of the constructor's
returned weight under the effective premises before forming or simplifying a quotient. If the zero region
is not proved empty, emit the proof status and explicit undefined payloads for the extracted density and
extraction residual.

The extraction residual and reference residual are separate diagnostics. The first compares an extraction
from the assembled action with the independently declared density. The second compares that declared
density with the independently computed Levi-Civita construction. Do not replace them with one residual
and do not reconstruct the action from a quotient. The reference route must call the CAS Levi-Civita
operator and indexed summation; an expanded polynomial literal is forbidden.

This construction can fail: with `L0`, `w`, and `S` treated as algebraically independent reports,
`L0/w - S` is not an identity. A mutation confined to any one report changes at least one diagnostic.
The reference residual has no dependency on the weight, so a coefficient-only mutation is not silently
relabelled as a form computation.

Emit this exact Q7 quantity family for every declared `(package,D)` pair:

```text
Q7_APPLICABILITY
Q7_GRADIENT_SYMBOLS
Q7_STATIC_LAGRANGIAN
Q7_SPATIAL_WEIGHT
Q7_WEIGHT_DOMAIN_STATUS
Q7_SPATIAL_DENSITY_DECLARED
Q7_SPATIAL_DENSITY_EXTRACTED
Q7_EXTRACTION_RESIDUAL
Q7_LEVI_CIVITA_VECTOR
Q7_LEVI_CIVITA_NORM
Q7_REFERENCE_RESIDUAL
```

At dimensions outside the Q7 predicate, every member remains present with the common not-applicable
payload envelope. Presence must not encode applicability.

## A deterministic sign classifier

Route every sign call, including calls inside stratum recomputation, through one classifier. Let `P_eff`
be the conjunction of joint premises, solver conditions, and active stratum conditions. A solver query
returns exactly one of `SAT`, `EMPTY`, or `UNKNOWN` and carries a witness for every `SAT` result and the
unsettled expression for every `UNKNOWN` result.

Stage 0 queries satisfiability of `P_eff` over the reals. `EMPTY` maps to
`INCONSISTENT_PREMISES`; `UNKNOWN` maps to `CAS_UNDECIDED_PREMISES`; only `SAT` continues.

Stage 1 partitions `P_eff` into these mutually exclusive regions:

```text
DEFINED_REAL     := root_is_defined AND Im(root) = 0
DEFINED_NONREAL  := root_is_defined AND Im(root) != 0
UNDEFINED        := NOT root_is_defined
```

`root_is_defined` must be computed from the expression's actual domain conditions, including denominators
and branch/function domains. Do not infer it from the absence of a CAS exception.

Stage 2 is reached only when Stage 1 determines `DEFINED_REAL`. It partitions that region into
`POSITIVE`, `ZERO`, and `NEGATIVE` using the real root itself.

Use this one total partition reducer at both stages:

```text
if at least two regions are SAT:
    PREMISES_INSUFFICIENT(regions, witnesses)
else if exactly one region is SAT and every other region is EMPTY:
    DETERMINED(region, witness)
else if any region is UNKNOWN:
    CAS_UNDECIDED(unknown_regions, unsettled_expressions)
else:
    CLASSIFIER_ERROR(all_regions_empty)
```

For a determined Stage 1 region, continue only for `DEFINED_REAL`; map the other two region names to
`ROOT_NOT_REAL` and `ROOT_UNDEFINED`. Because non-real and undefined are disjoint regions reduced together,
the classifier has no label-precedence choice when the premises admit both: the partition reducer reports
premise insufficiency with both witnesses.

Emit, for every sign call, the effective-premise status, the three domain-region statuses and witnesses,
the three sign-region statuses and witnesses or their explicit not-reached envelopes, and the final
classification. Use one registered suffix family everywhere; do not emit a bare CAS `Sign` result as a
substitute.

## Names, scopes, and payloads

Build the quantity registry from scratch from the tag names in both committed outputs. Parse only the
substring before the first colon. For each legacy tag, record exactly one canonical quantity or classify
it as engine-local. No tag may be unclassified, multiply classified, or assigned by payload inspection.
Record legacy aliases for migration, but prescribe only canonical names for future output.

Each registry entry must declare its scope shape, payload AST type, status envelope, legacy aliases, and
classification `shared` or `engine_local`; there is no third classification. Retire the unsuffixed
`Q3_DETERMINANT` spelling as a future emission name and map each legacy occurrence to either the canonical
raw-determinant entry or the canonical factored-determinant entry. Determine that mapping from engine
source, not from a payload comparison.

Use this complete tag grammar:

```text
<ENGINE>_S10_<LOCALITY><SCOPE>_<QUANTITY>

<LOCALITY>  ::= "" | "LOCAL_"
<SCOPE>     ::= "RUN" | <PACKAGE> "_D" <n> <STRATUM> <ROOT_SCOPE>
<STRATUM>   ::= "" | "_STRATUM" <positive-index>
<ROOT_SCOPE>::= "" | "_ROOT" <positive-index>
                   | "_ROOT" <lower-index> "_ROOT" <higher-index>
```

`LOCAL_` therefore occurs immediately after `S10_`. Root pairs use increasing indices. Define root and
stratum indices by sorting their canonical serialised defining expressions, not by native CAS return order.
Use fixed collection tags as well as scoped detail records wherever a computed collection may be empty, so
absence never carries the result.

Choose one payload serialisation for both engines: canonical JSON over the contract's expression AST.
Define encodings for every payload kind appearing in the registry, including integers, reduced rationals,
symbols through a shared symbol lexicon, booleans, status tokens, sums, products, powers, equations,
logical formulas, lists, sets, matrices, piecewise expressions, conditional expressions, solver records,
and explicit undefined/not-applicable envelopes. Define canonical operand ordering and set ordering. Raw
Wolfram `InputForm`, Python `str`, and Python `srepr` are not cross-engine serialisers and may not be
selected. This is another explicit future-engine requirement.

The live `checks_S10.yaml` rows below currently carry the two Q7 spellings and must be replaced by rows
whose `wl` and `py` values both use the canonical suffixes above:

```text
main_d3_q7_stiffness
main_d3_q7_difference
xcoef_scale_d3_q7_stiffness
xcoef_scale_d3_q7_difference
xform_aniso_d3_q7_stiffness
xform_aniso_d3_q7_difference
xform_divonly_d3_q7_stiffness
xform_divonly_d3_q7_difference
xform_fullgrad_d3_q7_stiffness
xform_fullgrad_d3_q7_difference
xform_signflip_d3_q7_stiffness
xform_signflip_d3_q7_difference
```

The replacement names for each package are `q7_spatial_density_declared` and
`q7_reference_residual`. A future harness edit must also add comparison rows for the other Q7 suffixes and
for dimensions outside the applicability predicate. Generate the exact companion manifest as the Cartesian
product of the package-dimension table and the eleven Q7 suffixes, then require set equality between that
manifest and the Q7 rows in the future config. Store the fully expanded manifest in
`emissions.q7_companion_rows`; report only whether it is complete. Do not edit the config in this task.
This companion change supersedes the earlier engine-only instruction to preserve existing tag names.

## Result-leak barrier

The new contract describes inputs and operations. It may define a supplied equation, but it may not carry
a computed equality, expected/observed/measured field, answer-shaped hint, control outcome, comparison
outcome, or success condition. Delete such content rather than relocating it.

The closed JSON format is part of this barrier. Additionally validate that:

- no node kind or key contains `result`, `expected`, `measured`, `observed`, `verdict`, `pass`, or `fail`;
- no emission tag encodes a value, sign outcome, count value, ratio value, or comparison outcome;
- no output node contains a literal payload;
- no numeric token or expression copied from the findings file occurs in a computation-output position;
- no historical or explanatory record exists anywhere in the document.

Do not include a worked counterexample. Test the data-flow properties with fresh abstract symbols only and
report mutation detection as changed/unchanged, never with an expression value.

## Acceptance protocol

Acceptance is conjunctive: any rejection, unclassified legacy tag, missing required instance, or
undecided schema check means continue editing. There is no acceptable `violates` state.

1. Strictly parse the entire target as JSON with duplicate-key rejection. Validate it against the closed
   authoring-time schema, identifier grammars, operator/enum allowlists, unique-owner rule, reference
   integrity, and acyclic definition graph.
2. Extract `q7` from the target, canonicalise it with sorted object keys and compact separators, and compare
   it with `AUTHOR_Q7_ORACLE` extracted from this directive. The comparison must be byte-empty. Apply the
   same AST comparison to the action, Route A, and dimension definitions printed above.
3. Perform graph reachability checks for the three Q7 constructor routes. Reject cross-route composite
   dependencies, a Levi-Civita dependency in a package report, division before the weight-domain guard, a
   typed expanded reference polynomial, or either forbidden quotient/reconstruction pattern.
4. With fresh algebraically independent symbols, test the extraction and reference residual DAGs. Mutate
   `L0`, `w`, `S`, and the Levi-Civita reference one at a time. Each mutation must reach its designated
   diagnostic, and the reference residual must have no path from `w`.
5. Enumerate all status triples in `{SAT, EMPTY, UNKNOWN}^3` through the partition reducer. Require exactly
   one verdict for every triple. Verify the Stage 1 predicates are pairwise disjoint and exhaustive under
   `P_eff`, and verify all sign-call nodes reference this classifier.
6. Build the legacy-name inventory independently from each `.out`, before the first colon only. Require
   every name to map exactly once, every shared canonical registry entry to have both engine roles or an
   explicit future-migration status, and every engine-local entry to carry `LOCAL_` in the grammar.
7. Expand the declared package-dimension matrix against all fixed-scope registry entries. Compare it by set
   equality with required emissions. Separately expand the Q7 companion manifest and print the future config
   rows that must replace or be added.
8. Validate every payload schema against the single canonical JSON-AST serialiser. Reject a registry entry
   that permits engine-native text or has no payload type/status envelope.
9. Run the leak barrier over every JSON key, enum, identifier, node kind, and literal role. Manually inspect
   only supplied-definition literals; do not inspect committed output payloads.
10. Try to game the contract. At minimum: append a prose sentence, insert an unknown JSON member containing
    a redefinition, duplicate a definition key, add a second owner for `package.spatial_weight`, alter one
    Q7 operand, and add a result-assertion node. The strict parser/schema, unique-owner check, or external
    oracle must reject every candidate. Paste only the rejection category for each attack.

All acceptance commands must be read-only except for the target specification edit and must be wrapped in
`timeout 600`. Do not run either CAS engine.

## Handoff, at most 25 lines

Report the target path and line count; the Q7 independence/ability-to-fail result; the sign-partition
totality result; registry coverage as complete/incomplete without tag or payload totals; the exact Q7
config row identifiers requiring replacement plus whether the companion manifest was produced; and the
six adversarial rejection categories. List every future engine or harness obligation introduced by the
contract. State any premise or requested computation you believe is ill-posed. Do not report any engine
value, any direction of an engine comparison, or any computed physics conclusion.
