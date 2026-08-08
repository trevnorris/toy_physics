# Reduction registry v2

This directory is the phase-1 semantic input for reduction.  It contains the
finite scalar registry grown from the seven-scalar medium block, including the
S9 transverse-brane constants and their derived cone relation.  It is not a
corpus census.  The versioned contract is
`registry_schema.yaml`; `quantities.yaml` and `relations.yaml` are its two data
documents.

## Why `prefix-v1`

`residual` is a YAML nested-list prefix expression, not Python, SymPy,
Mathematica, or infix text.  There is one stored expression.  YAML loaders in
both engines return the same tree of exact integers, ASCII strings, and lists;
each reader maps that tree to its own expression objects.  Operator position and
arity determine grouping, so there is no precedence or implicit-multiplication
decision on which the engines can disagree.

No floating-point literal is valid.  Use an integer or `[Rat, p, q]`.  A
quantity leaf is `[Q, canonical-qid]`; engines must key identity by QID, never by
the display `symbol`.

| Prefix node | Arity | SymPy construction | Mathematica construction |
|---|---:|---|---|
| integer | 0 | `Integer[n]` | exact integer `n` |
| `[Q,id]` | 1 | QID-symbol table lookup | QID-symbol table lookup |
| `[Rat,p,q]` | 2 | `Rational[p,q]` | `Rational[p,q]` |
| `[Add,x,...]` | 2+ | `Add[...]` | `Plus[...]` |
| `[Mul,x,...]` | 2+ | `Mul[...]` | `Times[...]` |
| `[Sub,x,y]` | 2 | `x-y` | `Subtract[x,y]` |
| `[Div,x,y]` | 2 | `x/y` | `Divide[x,y]` |
| `[Pow,x,y]` | 2 | `Pow[x,y]` | `Power[x,y]` |
| `[Neg,x]` | 1 | `-x` | `Minus[x]` |
| `[Sqrt,x]` | 1 | `sqrt(x)` | `Sqrt[x]` |

A Mathematica parser can therefore be written directly as:

```wl
parse[n_Integer] := n;
parse[{"Q", id_String}] := qidSymbols[id];
parse[{"Rat", p_Integer, q_Integer}] := Rational[p, q];
parse[{"Add", xs__}] := Plus @@ (parse /@ {xs});
parse[{"Mul", xs__}] := Times @@ (parse /@ {xs});
parse[{"Sub", x_, y_}] := parse[x] - parse[y];
parse[{"Div", x_, y_}] := parse[x]/parse[y];
parse[{"Pow", x_, y_}] := parse[x]^parse[y];
parse[{"Neg", x_}] := -parse[x];
parse[{"Sqrt", x_}] := Sqrt[parse[x]];
```

Anything outside that grammar is an error.  In particular, neither reader may
call a general-purpose evaluator on registry text.  `applied_functions` is
empty in v1; adding `Apply`, derivatives, integrals, or profiles requires a
schema/language version bump and the DESIGN §4 closure checks.

## Quantity document

Every quantity record has exactly these required fields:

- `qid`: stable identity, matching `Q.<ASCII path>`;
- `symbol`: display spelling only;
- `kind`: one closed schema value among `parameter`, `field`, `coordinate`,
  `boundary-datum`, `observable`, `intermediate`, `control`,
  `discrete-choice`, or `function-profile`;
- `scope` and `regime`: nonempty string lists;
- `state`: `live` or `retired`;
- `counting_axis`: `continuous-model`, `convention-orbit`, or
  `discrete-structural`;
- `dimension`: a constant `[L,T,M]` vector or a `BoundDimensionLaw` declaration
  under `LTM-exponent-vector-or-bound-law-v2`, plus provenance;
- `aliases`: alternate input spellings; and
- `source_loci`: one or more `{path,line_start,line_end}` mappings relative to
  `research/pde_ledger_v2/`.

QID, symbol, and aliases form one global alias table.  Any collision or unknown
reference is invalid.  The current counting contract is the document's
`active_regime`; phase 1 includes every live `continuous-model` quantity in
that regime in the finite scalar ambient space.  Live `convention-orbit`
coordinates and `discrete-structural` choices are reported separately.
The throat radius `a` is absent: it is a defect-sector quantity, not a
medium-level input or output, and no defining relation for it is admitted here.
`Q.medium.n_eos` records the discrete structural EOS choice and is
not a continuous knob.  A quantity may also carry an exact `value` (integer or
`[Rat,p,q]`); the seed declares `Q.medium.n_eos.value: 5`.  Declared values are
not automatic numeric inputs.  A caller may opt into their use during forward
propagation with `allow_declared_defaults=True`; every default consumed is
reported with `DeclaredValueDefaultWarning`.  The opt-in preserves the
missing-input guard and prevents a newly declared value from silently changing
an existing caller's dataflow.

A bound law uses `dimension-prefix-v1` components, explicit local-name to
structural-QID `bindings`, and `reference_values`.  Its retained `exponents`
triple is the checked evaluation at those reference values; it is not the law.
`Ref`, `Add`, `Sub`, `Mul`, and `Neg` are the closed law grammar.  An unknown
binding target, a non-structural target, an unbound reference, or a
reference-vector mismatch fails closed.  The loader returns the unspecialised
object through `dimension_declaration` or `dimension_law`;
`dimension_at(qid, structural_values)` requires an explicit integer value for
every binding needed by that law.  For example,
`r.dimension_at("rho_br", {"D_brane": 4}) == (-4, 0, 1)`.

`Quantity.dimension` remains a compatibility view for legacy engine consumers:
it is the checked reference evaluation and therefore remains iterable at the
registry's declared structural value.  Registry dimension validation uses
`Quantity.dimension_declaration`, not that compatibility view.  This preserves
the unspecialised symbolic object for the later engine-symbolic comparison.

The law's `reference_values` must also equal the declared integer values of its
bound structural QIDs.  This ties the retained reference triple, the law's
reference evaluation, and the structural quantity used to name that reference
point; absence, non-integrality, or disagreement fails closed.

`dimension_law_check.py` contains no second declaration of the dimension laws.
It reads the symbolic operands already present in the committed S9/S10 outputs
and S11-Python output, normalizes their dimension symbol to the registry binding,
then prints each engine operand, registry-law operand, and residual before
guarding it.  The fixed transcript population catches a single-sided registry
coefficient change.  `rho_br` and `mu_R` have S9, S10, and S11-Python operands;
`B_comp` is corroborated by S11-Python alone, which is single-engine coverage
but is not circular.  Because the `mu_R` and `B_comp` laws are identical,
swapping those two S11 coefficient tags leaves every comparison unchanged;
attribution between them carries no independent content in this pin.

The production homogeneity gate still does not anchor an absolute `D`
coefficient: it checks relation homogeneity at resolved bindings.  A common-mode
shift applied to the registry and every selected transcript together leaves the
pin, the gate, and `dimension_law_check.py` green.  That one-wrong-shared-premise
class still needs independent re-derivation, and the medium-sector structural
laws described below remain deferred.  The pin printer reports S11-Python's
symbolic-minus-numeric-reference residual and its value at the declared binding;
interpretation is outside the production script.

The medium sector has a known unrepresented dependency.  `Q.medium.rho0` and
`Q.medium.K` have structural dimension dependence but remain bare constant
vectors.  `K` also depends on `Q.medium.n_eos`.  The shipped relation dimension
grammar cannot express that linkage: a `Pow` exponent in `relations.yaml` is a
bare integer and is treated as dimensionless.  Extending the medium sector
would require a declared medium spatial-dimension choice, bound laws for
`rho0`/`K`, and a versioned relation grammar plus dimension algebra in which a
structural exponent can be bound to `n_eos`; this round does none of those.

## Relation document

Every relation record has exactly the required fields named in
`registry_schema.yaml`.  `residual` is one `prefix-v1` expression whose value is
zero, or `null` for a prose-only candidate.  For an explicit definition it must
have the exact transport shape
`[Sub,[Q,designated_output],rhs]`.  The output must occur only in that left leaf,
and `input_qids` must equal the alias-canonicalized QID leaves of `rhs`.  This is
the dataflow check that prevents independently freezing an alleged output.
Multi-output source rows are normalized into one scalar record per output; R2
is consequently `R2.xi_h` and `R2.h0`.

`domain_measure_bcs` is explicit even for a finite scalar relation.  A reader
must evaluate every `denominator_guards` expression and reject zero or
undetermined denominators before numeric evaluation.  `assumptions` are typed
predicates over QIDs.  `literal_consistency` assertions name an integer node in
the transported residual by path and require it to equal a valued quantity plus
an integer offset.  R1 uses these assertions to require its coefficient and
power to be `n_eos` and `n_eos - 1`; this is a structural check, not a symbolic
exponent in the residual.  `source_locus` identifies the claim and
`execution_locus` identifies executable evidence.  Paths are ledger-root
relative and line ranges must exist.  `benchmark_refs` is reserved for external
anchors and is empty in this seed.

`provenance_status` is this closed enum:

```text
DERIVED-EXECUTED, PENDING, CALIBRATED, CALIBRATED-UNCOMMITTED,
UNCOMMITTED, CONVENTIONAL, IMPOSED, STRUCTURAL, PROSE-ONLY,
CLOSED-NEGATIVE, CONTROL-ONLY
```

Only the exact value `DERIVED-EXECUTED` is eligible for the earned constraint
set.  A parseable equation with any other status remains a registry candidate
but is refused.  Eligibility also requires active/live QIDs, compatible regime,
known operators/functions, nonrecursive explicit output, and transitive input
closure.  In this seed, scalar terminals are live `parameter`,
`boundary-datum`, `control`, or `discrete-choice` QIDs; an `intermediate` input
must itself be the output of an admitted relation.

The residue is every live active-regime QID that is not a designated output of
an admitted relation on the `continuous-model` axis.  Finite-block dimension is
the Krull dimension of the exact polynomialized constraints, computed from the
leading-monomial ideal of a grevlex Gröbner basis.  Zero residuals are removed
and nonzero constant residuals denote an empty locus.  This algebraic dimension
is the maximum component dimension; on a reducible locus it need not equal the
dimension of the positive real locus.  A positive-real claim therefore has the
explicit precondition that `certify_positive_real_dimension` finds an exact
positive smooth witness whose local Jacobian dimension equals the algebraic
dimension.  Witness ranks certify that precondition but never select the
dimension.  Unlike the historical fixture helper, `constraint_dimension`
rejects a constraint containing symbols outside its declared dimension
variables instead of treating them as coefficient parameters (for example,
`a*x - 1` over `(x,)`).  This fail-closed divergence prevents an accidentally
omitted ambient variable from silently changing the dimension problem.

## Dimensional-homogeneity gate

`dimensional_homogeneity_gate.py` audits every non-null residual directly from
the transport documents.  It reports `HOMOGENEOUS`, `INHOMOGENEOUS`, and
`UNDETERMINED` populations separately.  Additive terms must agree; multiplication
and division add/subtract vectors; `Sqrt` halves a vector; and `Pow` accepts only
a bare integer exponent.  Missing dimension data, a convention mismatch,
an unresolved or non-integer law binding, unsupported syntax, or a non-integer
exponent is `UNDETERMINED`, never a pass.

An inhomogeneous relation is a finding, not an exception: the report names the
relation source/execution loci and the dimension-declaration loci as the two
candidate culprit classes.  The per-QID provenance output also records whether
the dimension-owning stage uses `scripts/ledger_dimensions.py`; `false` records
migration state, not an error.  The gate exits `0` only when every residual is
homogeneous and every QID has complete dimension provenance; otherwise it exits
`1`.

The relation gate constrains differences only.  With the speed declarations
held fixed, adding any shared symbolic dimension-vector function (including an
arbitrary `D_brane`-dependent one, not merely a constant shift) to `rho_br`,
`mu_R`, and `B_comp` is invisible.  If speed declarations may move too, the
full invisible kernel is wider: perturbations satisfying
`delta(mu_R)-delta(rho_br)=2 delta(c_gamma)` and
`delta(B_comp)-delta(rho_br)=2 delta(c_L)` remain homogeneous.  Consequently,
a green relation gate is not evidence about any absolute brane declaration.

The committed-output pin compares each registry law with selected symbolic
engine emissions and therefore polices a registry-only `D_brane` coefficient
change against those fixed operands.  This is narrower than a fresh independent
derivation: a common-mode change to the registry and every transcript survives
the pin, the relation gate, and the law check.  Coverage is also asymmetric:
`B_comp` has only the S11-Python operand described above.

## Registry-dimension witnesses

`registry_dimension_witness.py` compares the main-cell dimension vectors named
by the existing engine configs with the declarations in `quantities.yaml`.
`registry_dimension_witness.yaml` classifies each selected source as `DERIVED`
or `ECHOED`, binds symbolic substitutions to that source, and declares the
power whose dimension the tag carries.  Template branches use the selected
branch dimension, not `D_brane`'s registry value; a differing branch is
`BRANCH_DIMENSION_MISMATCH`.  `ECHOED` is excluded from agreement, with a
nonzero echo residual reported separately as guarded `ECHOED_MISMATCH`.

`declared_emitted_convention` is an input.  It is printed as
`emitted_convention_input`, never described as an engine measurement.  Where
an engine emits axis-labelled solution components, the witness reads them in
the registry's `ordered_bases` order and checks that reconstructed vector
against the bare vector.  Sources without labelled components are
`UNDETERMINED`, even when their numerical residual is zero.

Every selected QID prints its declared vector, declared multiplier,
`multiplier * declared`, raw and specialised emitted vectors, axis evidence,
and `emitted - multiplier * declared` residual before the status guard.  Every
selected `(artifact, quantity)` pair has a row; zero-tag template expansions
are `NO_TAGS`, and artifacts with no rows are `NO_ROWS`.  Both are guarded.
`NOT_EMITTED` and the aggregate counts print the artifact set to which they
apply, so they are coverage accounting rather than project-wide claims.

## Forward propagation

`load_registry().evaluate_output(output, numeric_inputs)` resolves aliases and
recursively evaluates admitted explicit definitions from primitive/residue
inputs.  It checks denominator guards, assumptions, and the final zero
residual.  Supplying a derived QID on the dependency path is an error: the
reader recomputes it instead of accepting an independently frozen value.
Missing inputs remain errors even when a quantity declares a `value`; passing
`allow_declared_defaults=True` opts into those defaults and visibly warns with
the QID and value of each one used.

Example:

```python
from registry_read import load_registry
r = load_registry()
assert r.evaluate_output(
    "h0", {"K": 3, "rho0": 2, "mass": 7}
) == 60
assert r.dimension_at("mu_R", {"D_brane": 4}) == (-2, -2, 1)
```

A Mathematica reader must use the same algorithm: build
`designated_output -> relation`, recurse through `input_qids`, parse and
evaluate the RHS (the third child of explicit `Sub`), refuse a supplied derived
node, then verify the parsed residual is zero.

## Runnable checks

From this directory, all commands are bounded externally with `timeout 600` —
⚠ **except `registry_import_fence.py --verify`**, which regenerates the
discovered engines and needs roughly **35 minutes** (its own per-engine bounds
are 3600 s for Python, 910 s for Mathematica).
⛔ **Measured: it writes nothing to a redirected stdout until it finishes**,
because Python block-buffers to a pipe and nothing flushes — so under a harness
or log capture it shows **no observable progress**. ⭐ Run it with `python -u`.

```text
python registry_read.py
python acceptance_check.py
python able_to_fail.py
python able_to_fail.py --case vacuous
python able_to_fail.py --demonstrate-crash
python able_to_fail.py --demonstrate-empty
python able_to_fail.py --demonstrate-spoof
python dimensional_homogeneity_gate.py
python engine_dimension_pin.py
python dimension_law_check.py
python dimension_law_able_to_fail.py
python w3_acceptance_ablations.py
python w3_duplicate_pin_ablation.py
python w4_weaker_implementation_runs.py
python w4_shifted_registry_printer.py
python w4_pin_completeness_runs.py
python registry_import_fence.py --list
python registry_import_fence.py --verify
python dimension_law_able_to_fail.py --case wrong-law
python dimension_law_able_to_fail.py --case wrong-binding
python dimension_law_able_to_fail.py --case unbound
python dimension_law_able_to_fail.py --case unresolvable
python registry_dimension_witness.py
python registry_dimension_witness_able_to_fail.py --case inherited-config
python registry_dimension_witness_able_to_fail.py --case multiplier
python registry_dimension_witness_able_to_fail.py --case echoed
python registry_dimension_witness_able_to_fail.py --case branch-dimension
python show_reduced.py
python -m pytest -q
python w3_runner_check.py
```

`dimension_law_check.py` exits nonzero for a missing engine operand, a nonzero
engine-minus-registry symbolic residual, or a declaration-form failure.
`w4_shifted_registry_printer.py` is an adverse driver and intentionally exits 1
after the shifted registry makes the production pin fail.
`w3_duplicate_pin_ablation.py` relies on the computed baseline/mutated behavior
flip; it does not infer absence of a duplicate from a spelling denylist.
`w4_pin_completeness_runs.py` builds and tests weakenings that remove the
covered-QID, population-count, error-set, and unmapped-symbol guards.
`registry_import_fence.py --list` computes the engine rerun population from
`registry_read` import syntax; `--verify` runs those engines sequentially,
captures stderr separately, and diffs stdout against committed output.
`generate_rows.py` is explicitly retired and exits 2 without writing: its
historical engine inputs contain only numeric triples, so it cannot construct
required v2 bound laws without recreating a bare `D=3` declaration.

`acceptance_check.py` computes the registry payload first and compares last
against the medium portion of its in-module literal fixture; disagreement exits
1 and is not reconciled.  `able_to_fail.py` runs five child mutations, including
separate syntactic-duplicate and semantic-entailment controls.  A caught child
must both exit 1 and print its case-specific `ABLE_TO_FAIL_CAUGHT` marker
without stderr.  Exit 0 is an escaped mutation; any other outcome, including an
exit-1 traceback without the marker, is an error.  `--demonstrate-crash`
exercises that distinction and intentionally exits 1 with
`ABLE_TO_FAIL_HARNESS: ERROR`.  The harness also refuses any case count other
than the expected five.  Every echoed child line is framed by case and stream,
and the parent-only verdict token is escaped in child text.  The empty and spoof
demonstrations exercise these controls.  `show_reduced.py` traverses only
admitted explicit definitions, detects cycles or stalled substitutions, and
reports the actual terminal QIDs left in each fully reduced expression.

`python -m pytest -q` is the intended full-population runner.  Plain
`python -m unittest discover` does not collect the pytest population and now
fails an explicit runner-contract test with `WRONG_RUNNER` instead of ending in
`OK`.  `w3_runner_check.py` prints that wrong-runner output and guards the
nonzero exit/message.  `dimension_law_able_to_fail.py` runs its complete W3
mutation population by default; every child must exit 1 with its computed
`CAUGHT` marker, while any escape makes the aggregate exit nonzero.  Its escape
demonstration runs a weaker child and computes the observed exit, residual, and
classification from that process rather than printing a fixture.
