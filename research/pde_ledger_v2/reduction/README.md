# Reduction registry v1

This directory is the phase-1 semantic input for reduction.  It contains only
the ten-scalar medium block used by `scripts/midway_knob_audit_codimension_sympy.py`.
It is not a corpus census.  The versioned contract is
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
- `dimension`: `{convention: LTM-exponent-vector-v1, exponents: [L,T,M],
  provenance: {stage_id, stage_uses_shared_dimensions_module, source_locus}}`;
- `aliases`: alternate input spellings; and
- `source_loci`: one or more `{path,line_start,line_end}` mappings relative to
  `research/pde_ledger_v2/`.

QID, symbol, and aliases form one global alias table.  Any collision or unknown
reference is invalid.  The current counting contract is the document's
`active_regime`; phase 1 includes every live `continuous-model` quantity in
that regime in the finite scalar ambient space.  Live `convention-orbit`
coordinates and `discrete-structural` choices are reported separately.
`Q.medium.a_pin` is a continuous-model output of an admitted explicit
definition; `Q.medium.n_eos` records the discrete structural EOS choice and is
not a continuous knob.

## Relation document

Every relation record has exactly the required fields named in
`registry_schema.yaml`.  `residual` is one `prefix-v1` expression whose value is
zero, or `null` for a prose-only candidate.  For an explicit definition it must
have the exact transport shape
`[Sub,[Q,designated_output],rhs]`.  The output must occur only in that left leaf,
and `input_qids` must equal the alias-canonicalized QID leaves of `rhs`.  This is
the dataflow check that prevents independently freezing an alleged output.
Multi-output source rows are normalized into one scalar record per output; R2
is consequently `R2.xi_h`, `R2.a_pin`, and `R2.h0`.

`domain_measure_bcs` is explicit even for a finite scalar relation.  A reader
must evaluate every `denominator_guards` expression and reject zero or
undetermined denominators before numeric evaluation.  `assumptions` are typed
predicates over QIDs.  `source_locus` identifies the claim and
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
ambient QIDs minus the largest Jacobian rank observed across exact positive
constraint-satisfying solver branches and exact parameter witnesses.  The
helper never evaluates rank off the constraint locus.  Zero residuals are
removed and nonzero constant residuals denote an empty locus.

## Dimensional-homogeneity gate

`dimensional_homogeneity_gate.py` audits every non-null residual directly from
the transport documents.  It reports `HOMOGENEOUS`, `INHOMOGENEOUS`, and
`UNDETERMINED` populations separately.  Additive terms must agree; multiplication
and division add/subtract vectors; `Sqrt` halves a vector; and `Pow` accepts only
a bare integer exponent.  Missing dimension data, a convention mismatch,
unsupported syntax, or a non-integer exponent is `UNDETERMINED`, never a pass.

An inhomogeneous relation is a finding, not an exception: the report names the
relation source/execution loci and the dimension-declaration loci as the two
candidate culprit classes.  The per-QID provenance output also records whether
the dimension-owning stage uses `scripts/ledger_dimensions.py`; `false` records
migration state, not an error.  The gate exits `0` only when every residual is
homogeneous and every QID has complete dimension provenance; otherwise it exits
`1`.

## Forward propagation

`load_registry().evaluate_output(output, numeric_inputs)` resolves aliases and
recursively evaluates admitted explicit definitions from primitive/residue
inputs.  It checks denominator guards, assumptions, and the final zero
residual.  Supplying a derived QID on the dependency path is an error: the
reader recomputes it instead of accepting an independently frozen value.

Example:

```python
from registry_read import load_registry
r = load_registry()
assert r.evaluate_output(
    "lambda_gamma", {"K": 1, "rho0": 1, "mass": 5, "c_gamma": 1}
) == 1
```

A Mathematica reader must use the same algorithm: build
`designated_output -> relation`, recurse through `input_qids`, parse and
evaluate the RHS (the third child of explicit `Sub`), refuse a supplied derived
node, then verify the parsed residual is zero.

## Runnable checks

From this directory, all commands are bounded externally with `timeout 600`:

```text
python registry_read.py
python acceptance_check.py
python able_to_fail.py
python able_to_fail.py --case vacuous
python dimensional_homogeneity_gate.py
python show_reduced.py
```

`acceptance_check.py` computes the registry payload first and compares last
against the medium portion of the existing literal fixture; disagreement exits
1 and is not reconciled.  `able_to_fail.py` runs five child mutations, including
separate syntactic-duplicate and semantic-entailment controls.  Each child
deliberately presents a forbidden expectation and must exit 1; the parent
prints each failing output/exit code and exits 0 only when all five failures
are observed.  `show_reduced.py` traverses only admitted explicit definitions,
detects cycles or stalled substitutions, and reports the actual terminal QIDs
left in each fully reduced expression.
