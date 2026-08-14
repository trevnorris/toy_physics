# S11 WL engine — fix round 1 brief

Target: `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` at `46ba77c2`.
`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` wins every conflict; the build directive
`research/pde_ledger_v3/directives/S11_wl_build_directive.md`, including its run discipline, continues
to bind. Change nothing beyond what the four items require.

## Item 1 — the §Q6 dimension walk must compute on every exponent its own streams produce

At `:907` the branch predicate references `RationalQ`, which is defined nowhere in the file and is not
a language builtin. On a radical exponent the condition stays unevaluated and the emitted payload is an
inert conditional carrying module temporaries instead of the §Q6 object; the derived-homogeneity family
then computes its test object over those poisoned operands.

What must become TRUE: for every expression the engine's own cells send through the §Q6 families, the
dimension walk reaches a computed result — every branch predicate is a defined, evaluating test — and
the `_DIM_HOMOGENEITY_DERIVED` test object is computed from computed operands. The spec's §Q6 declares
the payload shapes; add no shape it does not declare.

## Item 2 — a §Q8b status token reports an attempted decision

At `:643` the component-count status is assigned by structure alone —
`If[freeParameters === {}, CONSTANT, UNDECIDED]` — and the certificate is the literal `True` (`:646`).

What must become TRUE: each status token is the outcome of an attempted decision on the actual
component object; the mapping from the attempt's outcomes to the spec's declared tokens is total — no
token excluded by the shape of the code rather than by the attempt's result; the certificate, when not
`NOT_APPLICABLE`, is the attempt's own exact CAS object over the component's free parameters; the
change-locus objects are produced by the attempt when its outcome calls for them; `VALUE` reports what
the computation obtained.

## Item 3 — a failed exact-point extraction must be observable, never a silent drop

At `:1153`, `Select[strata, ListQ[#["Point"]] &]` silently removes an admissible stratum from the
emitted ordering when exact-point extraction fails. The spec declares an exact `STRATUM<s>_POINT` and
no failure payload, so the honest outcomes are exactly two: an exact point, or an operational failure
of that extraction surfaced through the engine's operational reporting.

What must become TRUE: an admissible component whose point extraction does not return a point produces
an observable, attributable operational outcome — never a silent value-conditioned drop, and never an
invented payload shape or fabricated point. Whether a tag appears depends only on package, dimension
and quantity (spec §5 corollary 4).

## Item 4 — a solver non-answer must not become a proved claim

At `:198-207` with `truthStatus` at `:91-95`: `inconsistentTest = SameQ[unrestrictedReduction, False]`
maps an unevaluated `Reduce` result to the decided boolean `False`, which `truthStatus` then emits as
`PROVED_FALSE` — a proved claim produced from a solver non-answer. The spec pins the payload
(`:253-260`): `STATUS_TOKEN` exactly one of `PROVED_TRUE · PROVED_FALSE · UNDECIDED`, `TEST_OBJECT`
the live boolean-valued object. Live streams do reach this site with radical systems (14 of 79
`_INCONSISTENT` families in the preserved `XKIN_ANISO` D2 capture carry `Sqrt`), so the branch is not
unreachable.

What must become TRUE: a proved token is emitted only when the underlying reduction actually decided;
a reduction that returns without deciding maps to `UNDECIDED`. The same property holds for every site
that maps a test object to a proved token.

## Acceptance — executable, no expected values

Probes may extract engine functions into a `/tmp` harness (the file runs its top-level on load); every
probe names its script and literal stdout.

1. Walker identity residuals: on symbolic atoms with declared dimension vectors, emit — operands and
   residual, then guard — the §Q6 algebra identities for the walker: `dim(x·y) − dim(x) − dim(y)` and
   `dim(x^r) − r·dim(x)` with `r` a symbolic rational exponent including `1/2`. Every residual is a
   computed object; perturbing one atom's declared dimension moves the operands. (A walker with a
   fixed offset or a broken rational branch fails the residuals; no expected dimension value appears.)
2. On a demonstration cell whose roots carry a radical, every `_DIM_` family payload parses as the
   spec-declared object: no unevaluated conditional heads, no module temporaries. The `XKIN_ANISO` D2
   cell emits the family early enough to observe within the run-discipline budget.
3. On the demonstration cells, every emitted count record's certificate that is not `NOT_APPLICABLE`
   carries as `TEST_OBJECT` the attempt's computed object: on a `/tmp` copy, ablating the component's
   defining data moves it — a `TEST_OBJECT` byte-identical under that ablation fails, whatever its
   value. The emitted token equals the total mapping applied to the attempt's outcome; where an
   attempt outcome cannot be produced honestly on the demonstration cells, show the attempt's literal
   returned object and record why — never a fabricated demonstration.
4. On a `/tmp` copy with only the `:1149` extraction call corrupted to return a non-list, the
   demonstrated admissible stratum remains present in the emitted record and the failure is observable
   in the run's output; the corrupted run's output must differ from the baseline run's (an empty diff
   fails — it means the extraction was not on the emission path).
5. Unit probe on the item-4 mapping: applied to an unevaluated (inactive) reduction object, the
   emitted token is the undecided one; applied to the decided booleans, the proved tokens. Run through
   the engine's own mapping function on a `/tmp` harness, operands and outcome printed.

Regression scope — which tags may move and which must not — is owed to the script review legs after
the build, judged from the diff; it is deliberately not a builder-facing criterion.
