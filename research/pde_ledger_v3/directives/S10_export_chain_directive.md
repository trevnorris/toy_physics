# S10 export-chain directive

## Authority and scope

Apply `S9_REWRITE_PLAN.md` to the S10 SymPy engine, using the committed S9 engine and its generated
`S9_exports.py` as the implementation reference. Change only the S10 SymPy engine, its generated export,
its regenerated stdout record, and this directive. Do not change S9, a Wolfram engine, S11, or any
registry, YAML, runner, comparator, or committed test machinery.

The Stage 1 stdout capture is the value baseline. Preserve every exported `MAIN` emission by its full tag
name. Compare those name-bound payloads outside the engine so the comparison remains independent of the
export writer. Report any movement outside that export boundary separately; do not force a control-package
readout to stay fixed when the imported assumptions change its computation.

## Imports and declarations

Import S9's flat `LEDGER`. Bind every S9 object that S10 consumes from its record; do not construct a
look-alike symbol or unit object under the same name. In particular, bind the material knobs, the
structural dimension symbol, the unit-basis objects, and the supplied unit premises used by S10. Use the
single module-level `omegaSquared` object at every construction site; a bare same-named symbol is not the
spectral coordinate. Objects such as neutral, wave-number, and spectral unit expressions that follow from
the imports are derived in S10. Scan every value and dimension in the merged ledger for same-named symbols
with different assumptions; this is the full defect class that import binding is meant to remove.

Every declaration site must carry exactly one trailing annotation of the form
`# TAG · English description`, where `TAG` is one of `KNOB`, `STRUCTURAL`, `COORDINATE`, `CONTROL`,
`PREMISE`, or `DERIVED`. The description names the role of the declaration and never states a result.
Supplied premise expressions remain inputs; do not replace them with typed conclusions.

## Export boundary and names

Export every non-local emission from each `MAIN` construction. Control-package emissions remain ablation
evidence and are not exports. Preserve all imported S9 records, then add the S10 records to the same flat
mapping.

An exported key is an authored object name. A solver root, coefficient position, stratum position, or any
other runtime index must not become any part of a key. Collect every indexed emission into the fixed
`indexed_derivations` object, with the authored field name, runtime indices, payload, and available unit
readout stored inside the collection. Emit the collection even when its varying payload is empty.

The key suffix records the component count of the construction path. Pass that construction metadata to
the key writer explicitly. Do not infer it from free symbols, payload shape, payload content, or a solved
answer. Dimension-solve inputs constructed in a fixed-component action carry the suffix. Imported unit
premises and symbolic dimension-solve outputs live in fixed slot collections under unsuffixed keys.

## Upstream overwrite protocol

S10 independently re-solves the inertia-coefficient dimension, stiffness-coefficient dimension, and their
difference. Those authored names overwrite the corresponding S9 records. For every S10 record whose key
already exists in the imported ledger, emit one collection row containing the key name, the committed S9
value, the freshly computed S10 value, and an exact reconstruction residual. Guard only after all three
objects have been emitted. Apply the same mechanism generically so an accidental collision also reaches
the residual instead of silently replacing upstream data.

Do not add an in-run chain-integrity comparison. Untouched-chain integrity is established by diffing the
committed generated modules.

## Generated module

The generated module is a deterministic function of the run. Every export-path collection with observable
iteration order is ordered explicitly; no set or hash iteration may decide record order.

At the start of a run, remove any prior `scripts/S10_exports.py`. Build the new module source in memory,
reconstruct it from that source, compare exact Python type and exact `srepr` against the live merged
records, emit the writer readouts, and apply the existing guards. Write the module only after every guard
passes, so a failed run leaves no generated module that a consumer could mistake for a validated result.

Use S9's readable record shape: `display`, exact reconstructible `value`, available computed `dim`,
`class`, and last-setting `step`. Preserve upstream record metadata exactly unless S10 overwrites that
key. Reconstruction must preserve unevaluated relational operands as well as their relational type; it
must not evaluate a stored numeric relation into a Boolean verdict. Emit the writer readouts and guards
that S9 emits; add no further writer self-check.

## Recorded limits and movement

The deleted registry comparison was the engine's only comparison between coefficient dimensions derived
from the imported field-dimension premise and independently declared coefficient dimensions in the
registry. It thereby policed the premise that the displacement field carries the supplied unit dimension.
Nothing in this export-chain build replaces that comparison. Changing the imported premise while leaving
the upstream coefficient records unchanged is caught by the overwrite residual, but changing the premise
and those upstream records consistently leaves every guard green and can write coefficient dimensions
wrong by two powers of length. The field-dimension premise is therefore unfalsifiable within this engine.

`coefficient_dimension_difference` is computed from
`stiffness_coefficient_dimension - inertia_coefficient_dimension`, and the field-dimension premise cancels
in that subtraction. The difference record adds no falsification power over the two records from which it
is built.

The control readouts that moved when S10 adopted the imported assumptions are:

- `PY_S10_XFORM_ANISO_D3_Q8_STRATUM1_ROOT2_Q3_SIGN`: before
  `undecided_under_joint_assumptions`; after `1`.
- `PY_S10_XFORM_ANISO_D4_Q8_STRATUM1_ROOT2_Q3_SIGN`: before
  `undecided_under_joint_assumptions`; after `1`.

## Acceptance evidence

Keep all acceptance scripts and their literal stdout outside the repository.

- Run the unmodified final engine first and show that the harness reproduces the committed stdout.
- Compare surviving Stage 1 emissions to the final run by full tag name and payload.
- Mutate the construction component count across runs and show that fixed-construction keys follow it,
  while symbolic-output keys do not acquire a suffix. State that this separates live construction metadata
  from a literal only across counts; it does not decide physical dependence on the structural symbol or
  police the legacy stdout tag spelling.
- Redirect an S10 record to an existing S9 key with a wrong value and show the emitted upstream operand,
  S10 operand, nonzero residual, and failing guard.
- For each export guard, apply the weakest mutation in its stated rejection domain and retain the failing
  stdout or traceback. State what each guard does not establish.
- Run the same unchanged engine under several Python hash seeds and compare the generated modules byte for
  byte. Repeat with the ordered symbolic-slot collection weakened to hash iteration and retain the
  divergent module diff.
- Reintroduce a bare `omegaSquared` at only one construction site and show that the namespace-wide symbol
  scan rejects it; also exercise determinant substitution with the stored solve variable and the module
  coordinate. Retain the clean scan from the unmodified final engine.
- Force an existing export guard to fail after pre-seeding an output module and show that the final engine
  leaves no module. Repeat with only the write moved behind validation but prior-output removal omitted;
  the stale module must make that weaker change fail acceptance.
- Run a construction that actually produces an allowed stratum. Compare the live and reconstructed
  relational objects, including their operands, type, and `srepr`. Repeat with plain `srepr` to `eval`
  restoration and retain its rejection.
- For each record-only item above, remove the smallest required clause or movement entry and show that the
  directive audit rejects the weakened text.
- Perturb every declared input, including an action-form perturbation. Report, for every exported
  `DERIVED` entry, which perturbations move it and explicitly list entries moved by none.

Do not write measured counts, tallies, partition sizes, dimensions, ranks, multiplicities, or other
computed results into this directive. The engine emits them for the reader.
