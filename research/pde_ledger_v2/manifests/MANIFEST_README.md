# Stage Manifest System v2.1

The v2.1 schema (`stage_manifest_schema_v2.json`) defines the shape of a stage
interface. This file defines its semantics: how expressions and dimensions are
interpreted, what the composite build checks, and how manifests are produced
without becoming a second source of truth.

The per-stage audits are unit tests; the composite build is the integration
suite. The current ledger has 44 stages. Partial extraction coverage is normal,
but unresolved material edges must remain visible.

## ⭐⭐ SCOPE — THE SEMANTIC CORE (user decision, 2026-07-29/30; corrected 2026-07-30)

The manifest system **continues**. The **semantic core** — required, extended, and allowed to block — is:

1. **quantity identity** — `C1`
2. **citation integrity** — `C2`
3. **export and lifecycle enforcement** — `C3`
4. **dimensional relations** — `DIMENSIONAL_CONSISTENCY`
5. **the range-aware lifecycle census** — `C5`
6. **the dependency graph and its cycles** — `C6`
7. **mutation** — `C7`
8. **genesis completeness** — `C8`
9. **consumption completeness** — `IMPORT-COMPLETENESS`

⭐⭐ **C2/C3/C5/C8 ARE NOT BOOKKEEPING — corrected 2026-07-30, after a round that trimmed them.** Each was
verified to catch a way the *physics* could be wrong, which is the governing test
(`docs/development_pipeline.md`, *THE POSTURE*):
- **`C2` detects a CHANGED CONSUMED EQUATION** — it re-proves each `consumes.ref` against the producer's actual
  exported relation, so a downstream stage that cites an upstream result it no longer matches fails closed.
  That is a wrong derivation, not a wrong citation string.
- **`C3` detects CONSUMPTION OF RETIRED PHYSICS** — a `RETIRED` claim, or one carrying `discharged_by`, is not
  an operative export and cannot be newly consumed. Without it a stage can build on a result the ledger has
  already withdrawn.
- **`C5` catches a WRONG IRREDUCIBLE-COUNT RANGE** — it recomputes both convention branches from the register
  and compares them to stage043's typed `record_range` and stage044's `[41,50]`. The irreducible count is a
  headline claim of Part VII; a mis-partitioned or double-counted knob makes it wrong.
- **`C8` catches CALIBRATED / TARGET-MATCHED GENESIS** — it is the fit-vs-derive detector: it refuses an
  unsupported `origin: independent` and forces `coordinated` / `target_matched` to name what they were
  coordinated with. ⚠ Its own honest limit still stands (below): `origin: independent` is an *asserted*
  historical claim no `genesis` record can establish — the same concession `DIMENSION_REWRITE.md` §1 makes
  about engine independence. It rejects unsupported independence; it does not prove independence.
⇒ **These four are in scope, extend normally, and block normally.** ⛔ The earlier trim was wrong on the
merits and is withdrawn; do not re-trim them without a measurement that they catch only bookkeeping.
⛔ **Also retired as false:** the premise that dimension recovery covering "~16 of 43 scripts" is what
blocks the 44-stage fanout. Measured 2026-07-29/30: recovery reaches **10 of 43** (7 `class Dim` +
3 digest-pinned bare-tuple orders). ⚠ Converting a stage **lowers** that count only for the **seven** stages
that carry an exact `class Dim` the AST walk reads (005, 006, 007, 008, 009, 030, 031); everywhere else it
neither raises nor lowers it, because those scripts were never recoverable to begin with. Either way
conversion does not *raise* recovery, since the shared module exports `Dimension`, not `Dim`. See
`DIMENSION_REWRITE.md` §2.

## Composite build checks

**C1 — Quantity identity, ownership, and symbol aliases.** The ledger-wide key is
`quantity_id`, not a printed symbol name. Every occurrence of one `quantity_id`
must have the same named `dim`, assumptions, definition, and callable signature.
Each quantity has exactly one owning stage: the earliest stage whose referenced
claim actually binds the stage-local parse alias. A binding claim is a
`defines`/`eq` relation with that alias on the lhs, a `spectrum` whose
`kernel_symbol` is that alias, or an `operator_identity` whose `acts_on` is that
alias. Every later appearance must point `definition_ref` at the owner's
defining claim and declare the matching consume. A later `here/...` reference is
`FAIL FALSE_LOCAL_DEFINITION`; a reference to a claim that does not bind the
symbol is `FAIL DEFINITION_NOT_BINDING`.
`parse_alias` is the stage-local SymPy spelling. If one printed `name` denotes
two different quantity ids, every colliding symbol must carry an explicit,
consistent `alias` declaration listing the distinct ids; otherwise C1 fails.
This makes the stage031 `kappa` and stage044 `ell` collisions detectable without
conflating distinct quantities.

**IMPORT-COMPLETENESS — undeclared nonlocal use.** This is a separate mandatory
check before citation comparison. Walk every expression in every typed claim
payload and resolve each symbol through its `quantity_id` and `definition_ref`.
If a used symbol has a nonlocal `definition_ref`, the stage must have a
`consumes` entry whose `ref` resolves to that definition. Absence is
`FAIL UNDECLARED_IMPORT`. This is how a missing `consumes` entry becomes
machine-detectable; declared graph edges alone cannot reveal an omitted edge.
Ownership is established before this walk, so an extractor cannot suppress an
edge by co-declaring a false local definition.

**C2 — Citation integrity, fail closed.** ⭐ **IN THE SEMANTIC CORE — it detects a changed consumed equation; see §SCOPE.** Resolve each `consumes.ref`, apply only
its structured `substitutions`, and dispatch by `check`:

- `cas_equivalence`: relation to relation; prove the two normalized relations
  equivalent by exact CAS under producer assumptions and declared
  substitutions.
- `implication`: prove the producer relation implies the relation as consumed;
  direction matters.
- `specialization`: a deliberate narrowing. It requires
  `specialization: true` and a structured `domain`.
- `value_equal`: compare exact numeric values; float atoms are forbidden.
- `dim_equal`: compare named dimensions only. It uses `as_consumed_dim` and may
  not treat a formula as value-equivalent.
- `token_match`: compare true status/verdict string-token payloads exactly.
  It is never legal for a quantity.
- `spectrum_match`: require spectrum↔spectrum and compare `operator` and
  `kernel_symbol` exactly, `kernel`/`eigenvalue`/`gap` by exact CAS, and
  multiplicity as an integer.
- `range_match`: require record-range↔record-range and compare low, high,
  spread, convention axes, and components exactly.
- `adjudication_match`: require adjudication↔adjudication and compare
  `outcome_token`, `domain_cardinality`, and `oracle_digest`.
- `set_match`: require set-cardinality↔set-cardinality and compare count and
  elements (as sets).
- `opaque_quantity_match`: cite a quantity without re-expanding it, but pin
  both its `quantity_id` and the producer export `source_digest`.

Every substitution must be backed by an exported producer claim of
`kind: convention`. Prose in `note` is never interpreted as a substitution.
The mode is selected from the producer payload kind, not consumer preference.
In particular `cas_equivalence` on a non-relation producer is `UNSUPPORTED`.
An ineligible producer payload kind, unsupported expression, or inconclusive
proof is `UNSUPPORTED`, never a silent pass.

**C3 — Export and lifecycle enforcement.** ⭐ **IN THE SEMANTIC CORE — it detects consumption of retired physics; see §SCOPE.** Every export `claim_id` must be one
of its manifest's own claim ids. A consume of a claim not exported by the
producer fails. A `RETIRED` claim, or a claim carrying `discharged_by`, is not
an operative export and cannot be newly consumed. The checker also verifies the
producer/consumer C7 facet metadata is present.

**DIMENSIONAL_CONSISTENCY — manifest-internal dimensional algebra plus source
certificate.** Verifies manifest-internal dimensional algebra: relation
homogeneity, declared-vs-recovered agreement for symbols whose dimensions are
recoverable from the stage's audit script, and cross-stage agreement of shared
quantities. It does NOT independently certify that a stage's dimensions are
physically correct — that is owned by the stage's dual-engine unit audit.

For relations, both sides must agree; when `expected_dim` is present, lhs and
rhs must both equal it. Claims of `kind: dimensional` require `expected_dim`.
Callable functions and operators use their declared domain/codomain signatures.
Integer and exact rational powers, including `sqrt`, are supported. Derivatives
and integrals adjust dimensions by their measures; transcendental arguments
must be dimensionless.

Each new stage declares `dimension_basis` as a stable basis `id` and an ordered,
unique `axes` array. A dimension is a sparse map from any of those named axes
to exact rational strings; omitted axes are zero, and JSON numbers/floats are
not accepted. `dim_source_order` is an axis array for arbitrary bases. The
legacy compact permutations such as `"LTM"` remain valid, and omission of
`dimension_basis` means `{id: "LMT", axes: ["L", "M", "T"]}` so existing
manifests do not migrate. Cross-stage dimensional comparisons proceed only
when both the basis id and axes declaration agree; otherwise the checker emits
`CROSS_BASIS_DIMENSION_COMPARISON` instead of comparing exponent vectors.

The check also reads each symbol's `dim_source`, recovers the positional order
from the script's `Dim(...)` structure or an independently digest-pinned
bare-tuple order, and reads the raw tuple. It fails if the recovered order
differs from `dim_source_order`, if the raw tuple differs from
`dim_source_tuple`, or if transposition by the recovered order differs from the
named dimension map. Missing or unreadable source certificates and unsupported
function rules are `UNSUPPORTED` and fail the build. The named map is
authoritative for interchange, but it is still a manifest claim rather than an
independent physical-correctness proof.

**C5 — Range-aware lifecycle census.** ⭐ **IN THE SEMANTIC CORE — it catches a wrong irreducible-count range; see §SCOPE.** The composite config pins
`notes/parameter_register.md` by path and SHA-256; a changed digest fails before
the census. Treat each `knobs` item as a lifecycle
event keyed by stable `knob_id` and `registry_row`. Partition the complete
parameter-register universe exactly once into
`base|route_ful_debt|route_less_structure|force_mag|discrete_postulate`.
For every event, add both endpoints of `count_effect`; pending debt remains
counted, and departures contribute `{low:0, high:0}`. An `inherited` event does
not double count its origin. A `discharged` event requires
`discharge_evidence` resolving to a DERIVED claim with executed engine
evidence.

C5 computes both convention branches and compares the result to stage043's
typed `record_range` and stage044's `[41,50]` sensitivity. Every register row
must be classified. An unclassified or multiply classified row, invalid
discharge, or wrong endpoint is a failure.

`record_range.convention_axes` has a deterministic rule: sum all component
low/high endpoints; enumerate the Cartesian product of axis choices; add each
choice's `low_delta` and `high_delta`; then report the minimum branch low, the
maximum branch high, and `spread = high-low`. `range_match` requires exact axis
structure equality.

**C6 — Dependency DAG.** Build the stage-level graph from `consumes` and require
it to be acyclic. A cycle means two stages each treat the other as upstream.

**C7 — Cross-stage executable binding (separate, elevated check).**
`c7_binding` and `c7_expect` are optional because fictitious bindings are worse
than visible missing coverage. When present, an export binding declares a
producing primitive, mutation environment, mutation command, and semantic
facet; a consume expectation declares an injection point, facet used, and
expected first failing predicate. That predicate must be an existing
`verification.teeth[].predicate` in the consumer. In a composite session the
harness mutates the producer facet and requires the declared downstream tooth
to fire.

- `DECORATIVE_DEPENDENCY(S,E)`: the consumer stays green when its declared
  producer facet is mutated.
- `UNDECLARED_DEPENDENCY(S,E)`: the consumer fires for a producer facet it did
  not declare.

**C7 mutation-harness execution contract.** An edge is covered only when the
producer export has
`c7_binding{producing_primitive,mutation_env,mutation_command,exported_facet}`,
the consumer consume has
`c7_expect{injection_point,facet_used,expected_first_failure}`, and
`facet_used == exported_facet`. The harness runs `mutation_command` as a
subprocess with `env[mutation_env] = exported_facet`, then parses the **last**
stdout line as JSON of the form
`{"consumer_stage": "tooth_or_PASS"}`. The verdicts are:

- consumer returns `"PASS"` → `DECORATIVE_DEPENDENCY` failure;
- consumer returns a tooth other than `expected_first_failure` →
  `WRONG_FIRST_FAILURE` failure;
- consumer returns `expected_first_failure` → C7 pass; and
- any other loaded stage fires → `UNDECLARED_DEPENDENCY` failure.

`expected_first_failure` must name a predicate in the consumer's
`verification.teeth`. The mutator must honestly perturb the exported facet and
compute the downstream outcome; it must never hardcode the expected tooth. The
real 030→031 proof reads the exported facet from `C7_FACET`, changes the
zero-mode profile power from 2 to 1, computes
`-diff(f0,w,2) + V_H*f0` in SymPy, and emits the stage031 tooth only when that
residual is nonzero. Its `--decorative` variant ignores the facet and emits
`PASS`, demonstrating the `DECORATIVE_DEPENDENCY` failure path.

**What the pilot caught.** The 030→031→032 pilot exposed a spurious
`stage030/kernel_determinant_star` citation: the real cited result is
`stage030/kernel_stability`, with `D*=7/4`.

C7 is the derivation-independence test; the prose/script comparison in the
production workflow is not.

C7 runs on a declared causal slice. A slice containing consume edges cannot
report headline PASS until C7 has run on every such edge; without complete C7
coverage its best outcome is PARTIAL. The v2.1 build acceptance uses synthetic
decorative-dependency and undeclared-dependency fixtures. The real 030→031
zero-mode mutation is also populated and passing; its deliberately decorative
variant demonstrates the expected `DECORATIVE_DEPENDENCY` failure.

**C8 — Genesis completeness.** ⭐ **IN THE SEMANTIC CORE — it catches calibrated / target-matched genesis, i.e. fit-vs-derive; see §SCOPE.**
⚠ Its honest limit, which bounds the claim without dropping the check: `origin: independent` is an *asserted*
historical claim that no `genesis` record can establish — the same thing §1 of `DIMENSION_REWRITE.md` concedes
about engine independence. C8 **rejects an unsupported independence assertion**; it does not establish a
supported one. POSTULATED, CONV, and CALIBRATED claims require
`genesis`. `origin: independent` is a positive historical assertion and
requires dated genesis evidence containing a path+commit record reference,
record span, and the later claim it predates. `coordinated` and
`target_matched` require `coordinated_with` refs. `unknown` is the safe default.
C8 rejects unsupported independence assertions and verifies the direction and
dates of discharge references.

## Coverage and outcomes

Each check reports exactly one of `PASS`, `FAIL`, `PARTIAL`, or `UNSUPPORTED`.
`PARTIAL` is distinct from success: it is used when material edges remain
unresolved at current extraction coverage. `UNSUPPORTED` means the checker
could not perform a required proof and is a build failure, not a skip.

Every report includes a coverage matrix with:

- resolved/total citations;
- checked/total claims;
- unresolved-producer count; and
- covered causal closure.

A headline never reports PASS while material edges are unresolved. Under full
44-stage coverage, or for a declared closed slice, an absent producer is FAIL.

## Expression and payload conventions

- Every claim has a non-empty typed `payload`; `note` is never the only
  machine-readable statement.
- Payload kinds are `relation`, `operator_identity`, `spectrum`,
  `adjudication`, `record_range`, `set_cardinality`, and `token`.
- Canonical expressions are `sympy` strings. LaTeX is a human mirror.
- Parse with `sympy.parsing.sympy_parser.parse_expr`; use the merged
  `parse_alias` table as `local_dict`; keep implicit multiplication off.
- Exact arithmetic only: use `Rational(8,3)`, never `2.6667`.
- Function whitelist:
  `exp log sqrt sin cos tan sinh cosh tanh atan Abs Derivative Integral Sum
  Function Eq Rational oo pi`. Write `sech(x)` as `1/cosh(x)`.
- `holds_within` is a recursive conditions AST with only
  `all|any|not|xor|domain`; its leaves are claim refs or structured domain
  predicates.
- Named dimensions omit zero exponents. For example,
  `a_B = L^-2 M^1 T^-2` is `{"L":"-2","M":"1","T":"-2"}` in the legacy
  basis; arbitrary declared axis names use the same sparse exact-string form.
- Claim kind and payload kind are coupled: identity/inequality/dimensional
  claims use relations, operator identities use operator payloads, spectral
  claims use spectra, adjudication claims use adjudications, range claims use
  record ranges, set claims use set cardinalities, reductions use relations or
  sets, nogos use relations/adjudications/tokens, and conventions use relations
  or true tokens.
- Adjudication `bucket_counts`, when present, sum to `domain_cardinality`;
  adjudication `axes`, when present, have cardinality product equal to it.

## Evidence and extraction records

Every symbol, claim, knob, departure, and tooth carries evidence with a source
path, stable locus, SHA digest, engine, and method. `method: prose_only` makes
an asserted-but-not-script-verified claim visible; it must never be silently
promoted to engine evidence and therefore requires `engine: prose`. The
composite checker re-hashes every evidence `source_path`; stale digests fail.
Each manifest persists its extraction report and all source digests in
`extraction`.

Future enhancement, not part of v2.1: add `--emit-registry` to audit scripts so
they emit verified identities and dimensions directly for comparison.

## Production workflow

1. **Extract.** Read the stage prose and both audit engines. Transpose the
   script's declared source order into named dimensions, emit typed payloads,
   and persist per-item evidence and source digests. Genesis comes only from
   dated records.
2. **Validate locally.** Validate against the v2.1 schema, parse exact
   expressions, compare engine-grounded claims and dimension certificates, and
   run import-completeness.
3. **Prose/script consistency review.** A second extractor is blinded to both
   the scripts and the primary manifest and works from prose alone. Compare
   normalized semantic payloads and a claim/source coverage matrix. Persist
   discrepancies and adjudications. This tests prose/script consistency; it is
   explicitly not proof of derivation independence.
4. **Composite build.** Merge manifests, run IMPORT, evidence/adjudication
   integrity, C1–C6, then C7 on the selected slice, and emit the coverage
   matrix. C7 runs after the cheap structural and CAS checks.

The v2.1 pilot slices are `006→030→031→032` and `043→044`.

## Versioning and retrofit

`schema_version: "2.1"` is the folded v2 contract. A stage amendment creates a new
manifest revision and, where appropriate, a new amending claim; it never erases
the historical dependency direction. Keep the v1 schema and v1 example intact
for provenance.

For the 44-stage retrofit, extract dense consume graphs first, run C1–C6 plus
IMPORT-COMPLETENESS over every wave, then populate real C7 bindings. Causal
slices remain PARTIAL until their C7 edges are exercised.
