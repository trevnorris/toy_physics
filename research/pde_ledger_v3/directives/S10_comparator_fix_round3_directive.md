# S10 cross-engine comparator fix round 3 directive

## Authority and boundary

Repair only the existing cross-engine comparator and regenerate only its comparator-side stdout records
from the committed PY and WL transcripts. Treat both transcripts as read-only inputs. Do not modify or run
either engine, do not start `wolframscript`, and do not modify any Wolfram file, step record, or TeX file.
Add no configuration, runner, authored name correspondence, expected tally, or acceptance population.

The live repository changes for this round are limited to:

```text
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_comparator_fix_round3_directive.md
```

Every tally remains a computed emission. Neither this directive nor the comparator may state its measured
value as an expected result.

## Preserve constants and shadow namespace objects

A bare payload identifier must be classified from the object the parser would otherwise resolve, not from
an authored list of spellings. If it resolves to a closed SymPy symbolic object, retain that object: a
mathematical constant continues to denote itself. If it instead resolves to a callable, class, registry,
namespace, or other non-symbolic library object, supply a local free symbol so the library object cannot
capture a bare payload name. Callable heads and attribute namespaces remain syntax and retain their library
meaning.

This comparator assumes that a name denoting a constant has the same meaning in both engines. It does not
derive that assumption. A future pair of engines that assign different constants to the same name would
therefore be unpoliced by this instrument and requires an explicit comparison rule.

## Refuse transcript-refuted rename proposals

Keep row-local symbol substitution as diagnostic evidence only; it does not alter a residual, guard, or
disagreement classification. Before placing an inferred pair on the D12 worklist, inspect the other parsed
non-route payloads from each transcript. If both engines exhibit each proposed normalized spelling as a
distinct free symbol, whether in the same payload or separate payloads, the transcripts refute the claim
that they are two spellings of one object. Emit the rejected pair and the candidate and per-engine,
per-spelling distinction evidence, but do not propose it on the worklist.

For every pair not refuted by that evidence, retain it on the worklist and explicitly mark it `UNREFUTED`.
That label means only that this transcript pair supplies no contradiction; it is not acceptance of the
rename and does not create a comparator correspondence. Route word tokens remain content and supply neither
rename candidates nor distinction evidence.

## Decompose content divergences

Partition every computed content divergence into mutually exclusive diagnostic populations:

- a residual containing a shape, type, length, relational, or other structural mismatch;
- a route-tag row whose residual is a route word token rather than a structural mismatch;
- a remaining numeric or algebraic residual.

Emit each population and the residual between their sum and the computed content-divergence population.
Classification remains diagnostic and cannot make a failing row pass. Include the decomposition residual
in the final guard without encoding an expected population.

## Recorded load-bearing distinctions

The nullspace span residual and matrix-membership residuals catch different mutations. A rank loss can move
only the span residual, while a perturbed matrix against unchanged spans can move only membership; the latter
includes a re-pointed wavevector component that no other comparator residual catches. Retain both halves.

The currently unparsed shared rows are the Wolfram-payload parser refusing a Boolean that contains a
domain-membership assertion. Record that parser boundary; do not repair it in this round.

## Acceptance and saved evidence

Use constructed in-memory records and disposable transcripts outside the repository. Save every ablation
script and its literal stdout at named absolute paths outside the repository. The checks must reject weaker
repairs and preserve the earlier controls:

- Show a same-constant complex expression that the former captured-name list reports as divergent while
  the object-based parser compares it equal. Show separately that a constant is not proposed as a rename to
  a real payload symbol, and that an unlisted callable or namespace collision is still captured as a free
  payload symbol. Retain callable-head and attribute-namespace syntax.
- Construct a row-local spelling substitution for objects that both transcripts distinguish elsewhere.
  Show that row-local inference alone proposes it, while the repaired worklist emits it as transcript-refuted
  with evidence and does not emit it as an unrefuted proposal. Also show a genuine spelling candidate with
  no transcript contradiction remains present and is explicitly marked `UNREFUTED`.
- Check that the emitted content decomposition accounts exactly for its computed parent population on both
  constructed representatives and the committed transcripts.
- Re-run the round-two derivative-dependence, namespace capture, mapping-order, duplicate-accounting, and
  computed-nullspace-classification controls under the repaired comparator.
- On real committed nullspace payloads, show the complete residual moves for a wrong subspace and remains
  zero under a symbolic invertible change of basis. Also demonstrate independently that rank loss moves the
  span half and an unchanged-span matrix perturbation moves the membership half.

Regenerate the comparator and existing re-point ablation stdout without running either engine. Save the
complete before/after diff for every regenerated `.out` outside the repository and account for every changed
line. Literal stdout is the measurement; do not transcribe its measured tallies into this directive.
