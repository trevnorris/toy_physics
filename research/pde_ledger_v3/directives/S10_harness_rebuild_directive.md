# Rebuild the cross-engine harness and make it compare the action

Files in scope:

- `research/pde_ledger_v3/reduction/engine_output_checks.py`
- `research/pde_ledger_v3/reduction/checks_S10.yaml`
- `research/pde_ledger_v3/reduction/checks_S9.yaml`
- `research/pde_ledger_v3/reduction/test_engine_output_checks.py`
- `research/pde_ledger_v3/reduction/harness_ablation.py`, which does not exist on HEAD and is to be written

Do not commit. Do not edit `steps/`, `paper/`, `REBUILD_HANDOFF.md`, either engine, or any committed engine
output.

## Baseline and evidence rule

The baseline is HEAD, currently `c9114773`. Commit `92461853`, tagged
`wip-2026-08-05-unreviewed`, was reverted by `ab77f25d`. It may be read as prior art, but it is neither the
baseline nor known-good code.

Every diagnostic in this directive was measured against the HEAD harness, the two HEAD configs, the
engine source, or these committed outputs:

- `mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`
- `scripts/out/S10_brane_mode_spectrum_sympy_audit.out`
- the corresponding committed S9 outputs

Do not replace these diagnostics with an expected post-rebuild number. Run the rebuilt instrument and
report what it reads. If another factual premise cannot be reproduced from HEAD or committed output,
report the contradiction instead of coding around it.

## Purpose and result interpretation

The governing failure mode is a layer that reports success without having compared anything. Empty,
missing, mis-keyed, duplicate, filtered-out, unparsed, shape-incompatible, and unwalked inputs must be
distinguishable from a clean comparison.

`AGREE` and `DISAGREE` are comparison verdicts. `MISSING`, `UNPARSED`, `SHAPE_MISMATCH`, an empty or
cardinality-invalid operand, and a declaration gap are operational failures to compare. A physics
`DISAGREE` is not an operational failure and must not be made into one.

Part A requires the action to be compared. It does not require or predict agreement. An honest persistent
`DISAGREE` is a valid finding and is not a reason to keep changing normalization until the result agrees.

## Part A: compare the action

### A1. Declare the action objects exactly

The strings `LAGRANGIAN` occur in 78 WL tag names and 130 PY tag names. Those are substring counts, not
action-object counts; most of the hits are Q6 metadata about a Lagrangian.

The action Lagrangian is emitted once by each engine in each of these 13 package-dimension cells:

- `MAIN`: D2, D3, D4, D5
- `XFORM_FULLGRAD`: D3, D4
- `XFORM_DIVONLY`: D3, D4
- `XFORM_SIGNFLIP`: D3, D4
- `XFORM_ANISO`: D3, D4
- `XCOEF_SCALE`: D3

For every declared cell `(PACKAGE, D)`, the action pair is exactly:

```text
WL_S10_<PACKAGE>_D<D>_Q1_LAGRANGIAN
PY_S10_<PACKAGE>_D<D>_Q1_LAGRANGIAN_EXPANDED
```

Declare all 13 pairs as nonempty scalar expressions with cardinality one. A Q6 dimension, term-dimension,
homogeneity, or other companion tag is not the action and cannot satisfy this declaration.

The separate averaged-Lagrangian name family also has 13 cell-level primary hits per engine: WL emits
`*_Q6_AVERAGED_LAGRANGIAN_DIMENSIONS`, while PY emits `*_Q2_PERIOD_AVERAGED_LAGRANGIAN`. These are not
declared as semantic counterparts here. The averaged family is out of scope for this rebuild, and none of
its tags counts as an action Lagrangian or as action coverage. If it is admitted in later work, it
requires a separate declaration with separate WL/PY tag pairs and cardinalities.

The Euler-Lagrange pair in each cell is exactly:

```text
WL_S10_<PACKAGE>_D<D>_Q1_EULER_LAGRANGE_SYSTEM
PY_S10_<PACKAGE>_D<D>_Q1_EULER_LAGRANGE_SYSTEM
```

Declare each as a nonempty ordered vector whose outer cardinality is the declared dimension D. These rows
and the scalar Lagrangian rows are two separate required action families.

### A2. Compare the ordered variational residual vector

For an Euler-Lagrange row, the compared object is the ordered variational residual vector `δL/δu_i`.
Convert an emitted equation to `lhs - rhs`; where an engine already emits a
residual, use that residual directly. Preserve component order and overall scale. Emit the normalized
operands with the verdict.

Do not use a zero-locus or solution-set convention. Such a convention makes `x` and `2*x` equal because
they have the same zero set, erasing a wrong overall factor in the action. The acceptance test is an
ordered, same-shape unequal residual pair that must reach `DISAGREE`. Reordering unequal components must
also remain detectable.

HEAD normalization holds a WL `Derivative[...][x1,x2,x3,t]` operand in the `OpaqueDerivative` class and a
PY derivative as `Derivative(u(t,x1,...))`. On the committed `MAIN D3` row, equation-to-`lhs - rhs`
conversion makes the outer shapes compatible, but the three residual components still compare unequal.
Representation alignment is additional work. Any such alignment must preserve derivative order, arity,
function, argument list, and argument order and must not hardcode a coordinate list or ordering.

Do not promise that representation alignment will produce `AGREE`. If the operands remain honestly
different after the required structural normalization, report `DISAGREE` and the operands.

### A3. Report action verdict coverage

`SHAPE_MISMATCH` is already an operational failure in `HarnessReport.operational_failures` on HEAD. The
silent defect is different: the 13 configured Euler-Lagrange rows currently produce only `UNPARSED` and
`SHAPE_MISMATCH`, and no Lagrangian row is configured, so no configured action row reaches either
comparison verdict.

For each of the two declared action families, report:

- distinct configured rows that reached `AGREE` or `DISAGREE`;
- distinct configured rows that did not reach a verdict, grouped as `MISSING`, `UNPARSED`,
  `SHAPE_MISMATCH`, empty, or cardinality-invalid, with row names;
- the declaration-derived denominator.

It is an operational failure when either required action family has no row that reached a verdict.
Partial action coverage remains a reported gap under the coverage rule below.

## Part B: rebuild the other silent layers

### B1. Construct declaration-oriented coverage

HEAD has no coverage layer. Build one; do not describe this as repairing code from `92461853`.

For cross-engine comparisons, use this formula:

```text
cross_engine_coverage =
  distinct declared tag-pair identities with an AGREE or DISAGREE verdict
  -----------------------------------------------------------------------
  distinct declared tag-pair identities
```

A tag-pair identity is the ordered engine-to-tag mapping in the declaration. The denominator is fixed by
the config before output is read or rows are filtered. Quantity labels do not create additional pair
identities. Missing, unparsed, shape-incompatible, empty, and cardinality-invalid pairs remain in the
denominator and do not enter the numerator.

Print the numerator, denominator, formula, and named gaps. Removing output must decrease coverage or add a
named operational gap; it can never improve coverage. An absent, empty, or fully filtered comparison
declaration is an operational declaration failure. No guard may depend on the presence of an optional
config field.

The action coverage in A3 uses the same formula over each declared action-family subset.

### B2. Derive control identity from declarations

`CONTROL_TAG_PATTERN` requires `_X<digits>_`. It matches 1380 of 1559 S9 WL tags and no S10 WL tag. S10
uses package names such as `XFORM_FULLGRAD` and `XCOEF_SCALE`, so HEAD reports zero control comparisons
and zero parity rows for S10 without an operational failure.

Each config must declare package, dimension or other cell identity, required tag suffixes, and main/control
role. S10 must also declare stiffness form. Do not infer any of those from a tag prefix or regex. The S10
package-dimension cells are the 13 cells in A1. A control or parity layer that reaches no declared cell is
an operational failure.

The stiffness-form declarations are:

- `MAIN`: `curl`
- `XFORM_FULLGRAD`: `fullgrad`
- `XFORM_DIVONLY`: `divonly`
- `XFORM_SIGNFLIP`: `curl`
- `XFORM_ANISO`: `curl`
- `XCOEF_SCALE`: `curl`

These values are the third `Package` argument in
`scripts/S10_brane_mode_spectrum_sympy_audit.py`; the `XFORM_` prefix is not a form classification.

### B3. Run control and dimension checks per engine

HEAD selects `default_values = outputs[default_engine]` and runs control and dimension checks only on that
mapping. An inert or invalid non-default engine is therefore unobserved.

Run both layers separately for every declared engine. Every control, parity, dimension, coverage, and
operational result must name its engine. If only one engine is ablated, the report must show the finding
for that engine without borrowing data from or changing the result for the other engine.

### B4. Parse output lines without absorbing diagnostics

On HEAD, parsing `TAG1: 1`, then an untagged line, then `TAG2: 2` appends the untagged line to `TAG1`'s
payload. The committed S10 WL output also contains ten `Solve::svars:` diagnostic lines; they are treated
as duplicate `Solve` tags and reduce the CLI result to one error line.

Define a tag grammar that distinguishes emitted tags from CAS diagnostics. A line outside that grammar
must be reported, must not become a tag or part of a neighboring payload, and must not erase the rest of
the report. Print the count and the first few lines verbatim. A duplicate line that really satisfies the
tag grammar is an operational finding, but it also must not erase all other findings.

### B5. Declare dimension sources by engine, package, symbol, and shape

The two S10 `derived_dimensions` entries on HEAD are WL tags under `MAIN`; neither tag exists in the PY
output. A per-engine dimension run using that schema would silently build no PY derived table.

Replace the symbol-to-one-tag schema with declarations keyed by engine and package. Each symbol entry must
name its exact source tag or source-tag family and its shape. The shape vocabulary is:

- `vector`: exactly three scalar components;
- `family`: a nonempty sequence of three-component vectors that collapses only when every member has the
  same vector.

Within an engine-package table, the shorthand `symbol: TAG` means `shape: vector`. Do not infer `family`
from payload shape. A live S10 family can be a square list whose rows happen to agree, which is
indistinguishable by shape alone from a matrix.

Never take the first family member. A family with unequal members is ill-posed for a one-symbol mapping.
A declared source that is absent, empty, unparsed, or inconsistent with its declared shape is an
operational failure and leaves that engine-package table unassessable. It must not default to primitives,
another engine, or another package.

### B6. Use each package's own symbolic dimension table

The configured S10 sources on HEAD are `_SPECIALIZED` values at D3 and are applied across the sweep. Use
symbolic sources whose components retain `braneDimension` or `D`, not a vector specialized to one numeric
dimension.

Build and apply one table for each declared `(engine, package)`. If valid package tables disagree, report
the package names and vectors and continue to assess each package with its own vector. If a package lacks
a valid own vector, mark that package unassessable. Reporting a cross-package disagreement and then using
the other package's table is not a disposition and is forbidden.

### B7. Walk mappings and partition dimension accounting

`MAPPING` is absent from `DIMENSIONFUL_KINDS` on HEAD, and `walk_container` returns immediately on a
mapping. The payload `<|a -> (rhoBr + omega)|>` therefore reports no dimension error and no checked tag.

Walk dimension-bearing values in supported containers, including mappings. If any recognized container
or member cannot be walked, report it as `unwalked`; do not count the tag as compared, homogeneous, or
silently non-dimensional.

Assign exactly one primary dimension status to every emitted tag and print this arithmetic per engine:

```text
total_tags = compared + no_comparison + not_applicable + unwalked + unassessable + unparsed
```

The statuses mean:

- `compared`: at least one homogeneity proposition was evaluated for the tag;
- `no_comparison`: a dimension-bearing payload was walked but offered no evaluable homogeneity
  proposition;
- `not_applicable`: the parsed payload kind has no dimensional semantics by a stated rule;
- `unwalked`: a recognized container path was unsupported;
- `unassessable`: walking or a required dimension lookup failed visibly;
- `unparsed`: normalization failed with a reason.

`homogeneous` is a subset of `compared`, not another partition bucket. A residual bucket defined only as
the arithmetic remainder is not acceptable; each status must be assigned from the tag's observed path.

### B8. Count a tag as checked only after a comparison

`self.comparisons` occurs zero times on HEAD. The live defect is that `checked` increments after any
successful walk and `homogeneous` increments when that walk has no recorded issue. Measured probes `2`,
`x`, and `Element[x, Reals]` each report `checked=1 homogeneous=1` even though no homogeneity comparison
was performed.

A tag on which no homogeneity proposition was evaluated belongs in `no_comparison`, not `compared` or
`homogeneous`. Report the number of evaluated propositions by site kind so the content of `compared` can
be audited. A proposition true only by syntax, such as the dimensionlessness of a literal integer
exponent, is not a homogeneity comparison.

### B9. Parse the measured WL constructs without discarding structure

HEAD leaves 220 of 2983 S10 WL payloads unparsed and none of 4227 S10 PY payloads unparsed. The measured
WL substring populations overlap: 42 contain `Derivative`, 74 contain `Element`, 18 contain `Inequality`,
and 37 contain `ConditionalExpression`.

There are no top-level bare-string payloads. There are 71 unparsed containers containing quoted strings.
There is one unparsed real-literal payload, `WL_S10_RUNTIME_SECONDS: 60.62794`.

Implement these measured structures:

1. Parse a derivative of arbitrary arity applied to an arbitrary function and argument list. Arity,
   derivative orders, function, argument list, and argument order are identity, and no coordinate names or
   ordering may be hardcoded.
2. Parse set membership over `Integers` as well as `Reals`, including membership inside boolean chains,
   while preserving the set and member.
3. Parse associations and other live containers with quoted keys and quoted or list-valued members. Keep
   strings as strings and preserve container structure; do not coerce a container into a bare string.
4. Parse the named marker heads present in committed output when they have list arguments, preserving the
   head and argument structure. Do not install a rule that treats every head with a list argument as the
   same object.
5. Parse `Piecewise` with every branch, condition, and default preserved. Different branches or conditions
   must remain unequal.
6. Parse list-valued relations and represent `Inequality` as boolean structure with its operands and
   operators preserved, rather than as an opaque scalar function.
7. Parse the real-literal syntax used by `WL_S10_RUNTIME_SECONDS` and preserve its numeric value.

`ConditionalExpression` is not a parser defect on HEAD: its expression and condition are preserved,
different conditions compare unequal, and the dimension walk raises a visible error when it cannot handle
the condition. Do not claim to fix it by reproducing behavior already present. It may be changed only if a
measured remaining defect and its acceptance test are reported.

Where a construct genuinely cannot be parsed, retain the raw payload and reason as `UNPARSED`. Never drop
branches, conditions, keys, members, derivative arguments, or operators to make a payload parse.

### B10. Observe boolean truthiness exposure

Do not change `symbolic_equal` equality semantics in this rebuild. Its measured defect is truthiness
matching between Python booleans and numbers: `symbolic_equal(True, 2)` is true while
`symbolic_equal(False, 2)` is false. It also makes `symbolic_equal({a: True}, {a: 1})` true.

For this observer, a boolean value is exactly the population recognized by `_is_boolean_value`: Python or
SymPy booleans, objects whose `is_Boolean` is true, applied predicates, relations, and `Contains`. B9 will
expose more structured relations. For A2 rows, count the normalized residual tree sent to the comparator,
not the source equation wrapper that normalization removed.

Report these two row counts separately over distinct configured cross-engine rows after selection and
action normalization:

- `selected_boolean_rows`: rows in which at least one complete selected operand is a boolean value;
- `tree_boolean_rows`: rows in which at least one boolean value occurs anywhere in either selected
  operand's recursively walked tree, including sequence members, matrix entries, mapping keys and values,
  and rule sides.

These are row populations, not counts of boolean nodes. Name every row in either population. A nonempty
population is an operational failure because the deferred equality defect is exposed.

The acceptance fixture must be nested-only: compare `{a: True}` with `{a: 1}`. It must leave the selected
operand count at zero, make the tree count nonzero, and raise the observer's operational finding. A
top-level boolean fixture is insufficient because shape gating can prevent it from reaching the false
`AGREE` path.

### B11. Reject duplicate and empty comparisons

On HEAD, three identical configured rows produce three `AGREE` rows from one tag pair. Empty lists, empty
mappings, and empty multisets also report `AGREE` without comparing an element.

Count distinct tag-pair identities, not config rows. Duplicate declarations are an operational config
error, are named, and do not increase a verdict count or coverage denominator.

Every cross-engine object declaration must state a nonempty structural cardinality: scalar cardinality
one, exact positive outer cardinality for a sequence, exact positive entry count for a mapping, exact
positive element count for a multiset, or exact nonempty matrix shape. A selected operand that violates
the declaration gets an operational empty/cardinality status and can never get `AGREE`.

### B12. Give registry residuals declaration coverage

`registry_residual` is configured with no rows on both HEAD configs: S9 declares an empty list and S10
omits the key and defaults to an empty list. This is a declaration gap, not a clean residual check.

Apply the B1 rule to registry residuals. An absent, empty, duplicate-only, or fully filtered declaration is
an operational failure. Use this formula:

```text
registry_residual_coverage =
  distinct declared registry identities with a ZERO or NONZERO verdict
  ---------------------------------------------------------------------
  distinct declared registry identities
```

The denominator comes from the config declaration. Missing, unparsed, invalid-shape, or unsubstitutable
rows remain named gaps. `NONZERO` is an honest comparison verdict, not an operational failure.

### B13. Write the ablation instrument and select form controls correctly

`reduction/harness_ablation.py` does not exist on HEAD. Write it as a deterministic instrument that takes
the harness, config, and outputs as arguments and stays within the project script timeout. The reverted
version is only prior art.

The engines and committed output are frozen, so ablations operate on scratch copies or in-memory payload
substitutions. Do not describe payload substitution as an engine mutation.

For the form-response ablation, select packages by the declared stiffness-form value from the third
`Package` argument. `XFORM_FULLGRAD` and `XFORM_DIVONLY` are the genuine stiffness-functional controls;
their values are `fullgrad` and `divonly`, while `MAIN` is `curl`. `XFORM_SIGNFLIP`, `XFORM_ANISO`, and
`XCOEF_SCALE` all retain `curl` and alter a sign, inertial coefficient, or stiffness coefficient. The
`XFORM_` prefix is a trap and must never select the form controls.

The action comparison must show a changed normalized operand difference or changed verdict when a genuine
different-form payload is substituted at a shared declared cell. A pre-existing `DISAGREE` may remain
`DISAGREE`; in that case the changed operand difference is the required response. A coefficient-only
substitution may be detected by equality, as scale preservation requires, but it cannot create a
`FORM_RESPONSE` finding or satisfy the form ablation. Form-control classification must remain unchanged
when only coefficients change.

## Acceptance battery

Run every item against the candidate harness and paste literal stdout. Each item below names its oracle;
movement in an unrelated total or kind bucket does not satisfy it.

1. Run both committed configs to a complete report. The report must include per-engine control and
   dimension sections, both coverage formulas and declared denominators, both action-family verdict/gap
   sections, registry residual coverage, bucket arithmetic, ignored-line diagnostics, and operational
   findings. An exception or one-line fatal report fails this item.
2. For every declared scalar Lagrangian and Euler-Lagrange pair, print the row status. Print normalized
   operands for a representative row from each family. `AGREE` and `DISAGREE` both satisfy comparison;
   the other statuses do not.
3. Compare an equal synthetic ordered residual pair, then change one side to a same-shape unequal vector.
   The target row must move from `AGREE` to `DISAGREE`. Repeat with component order changed and show that
   order is not erased.
4. Delete every declared control tag at one cell from one engine only. The target engine's declared-tag
   coverage must decrease and its named missing-cell findings must increase; the other engine's control
   and dimension results must remain unchanged.
5. Delete all but one emitted main-package tag without changing the declaration. Cross-engine coverage
   must decrease and named declared-pair gaps must increase. Paste the formula and unchanged denominator.
6. Rename the main package in a scratch output, then separately rename one control package. In each run,
   matched declared cells must decrease and a named missing-package operational finding must appear.
7. Empty `cross_engine`, then separately empty `registry_residual`. The respective declaration-gap
   operational finding must appear. Neither `configured=0` report is a pass.
8. Duplicate one configured cross-engine tag pair under new quantity names. The distinct-pair denominator
   and verdict count must not increase, and a named duplicate-declaration finding must appear. Feed empty
   list, mapping, and multiset operands; each target row must get an empty/cardinality operational status
   and none may report `AGREE`.
9. Start with a control row reported `INVARIANT` and change one same-shape operand. That named row must move
   to `RESPONSIVE`; movement only in a total-tag or kind counter does not satisfy this item.
10. Use the nested-only boolean fixture `{a: True}` versus `{a: 1}`. The selected-operand count must be
    zero, the tree count must be nonzero, the row must be named by the observer, and its operational
    finding must appear.
11. Remove a declared PY dimension source while leaving WL intact. The named PY package must become
    unassessable and WL results must remain unchanged. Repeat in the other direction. This is the
    per-engine-table oracle.
12. Give one package a valid own vector that differs from another package's vector. A cross-package
    disagreement must be printed, each package must use its own vector, and unrelated package verdicts
    must remain unchanged. Then remove the target package's source and show that the target becomes
    unassessable rather than adopting another package's vector.
13. Replace a configured `vector` source and a configured `family` source, one at a time, with a bare
    integer. Each named source must get a shape operational failure and its engine-package table must
    become unassessable; no default table may be applied.
14. Feed the dimension layer `2`, `x`, and `Element[x, Reals]` as separate synthetic tags. Each addition
    must increase `no_comparison` without increasing `compared` or `homogeneous`. Then add a payload with
    an evaluable unequal-dimension addition; `compared` and the named non-homogeneous findings must
    increase. Print the evaluated-proposition site breakdown.
15. Feed a mapping whose value contains the same evaluable unequal-dimension addition. `compared` and the
    named non-homogeneous findings must increase, `unwalked` must not absorb it, and the printed partition
    equation must still hold.
16. Exercise each numbered B9 construct with structural equality and inequality fixtures. Each live form
    must parse; changing a derivative arity or argument order, membership set, quoted member, marker head,
    piecewise branch or condition, inequality operand or operator, or real value must remain observable in
    the normalized structure or verdict.
17. Feed a file with a stray solver diagnostic and a file with a stray line between two tags. The
    ignored-line count must increase, the lines must be printed, the neighboring parsed operands must
    remain unchanged, and the rest of the report must still be present.
18. Run the form-response ablation at a shared dimension for `MAIN` versus `XFORM_FULLGRAD` and for `MAIN`
    versus `XFORM_DIVONLY`, selecting solely by the declared stiffness-form value. Each must produce a
    `FORM_RESPONSE` with a changed normalized action operand difference or verdict. Run coefficient-only
    substitutions from `XFORM_SIGNFLIP`, `XFORM_ANISO`, and `XCOEF_SCALE`; they must not produce
    `FORM_RESPONSE` or satisfy the form-ablation status. Print the selected packages and forms.
19. Run the unit suite without invoking Mathematica. Give a one-line cause for every failure: regression,
    stale test, or test intentionally updated by this rebuild. Do not silently rebaseline a frozen
    expectation or delete a test.

## Constraints on conclusions

Do not make physics disagreement exit nonzero; only operational failure does. Do not trim a declaration
to match emitted output. Do not change `symbolic_equal` semantics in this scope. Do not hardcode expected
physics values, coordinates, coordinate order, or a post-rebuild report count.

## Report back

Keep the summary under 40 lines, followed by the literal acceptance output. Give one line for A1-A3 and
B1-B13 with file and line references, then list any partial or unimplemented item, the ablation oracle
results, and anything in this directive contradicted by measurement. Do not report an expected physics
value; report verdicts and operational findings read by the rebuilt instrument.
