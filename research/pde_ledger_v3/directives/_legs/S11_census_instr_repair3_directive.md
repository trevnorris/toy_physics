# Repair directive — S11 census instruments, round 3 (folded after two legs)

You are repairing the instrument modules under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` at commit `89ed80c9`, governed by
the folded brief (`S11_locus_census_instrument_brief.md`, `fa8c58b3`) — its eight obligations and
population section remain the acceptance surface. Two review legs on the round-2 repair, two
legs on this directive, and orchestrator re-computation verified the six defects below in the
round-2 census (the measured commands and literal outputs are in this directive's
`../_measurements/` twin — evidence lives there, not here). Repair them all, extend the
calibration, and re-run the full census. Everything the round-2 legs certified live must stay
live — the round-3 calibration must still catch a regression of every round-2 class.

## Verified defects

1. **One sibling branch poisons union coverage** (violates obligation 3). A candidate is covered
   when it lies in the UNION of the emitted branches. The current test multiplies one constraint
   from EVERY emitted branch (`s11_census_math.py:683-715`), so a branch undefined at the
   candidate makes the product undefined; the sampling fallback discards any sample where ANY
   sibling is undefined (`:639-680`); and `COVERAGE_UNDECIDED` is collapsed to `NOT_COVERED`
   (`:713`) and promoted to `OMITTED_BRANCH` (`:773-774`). Measured: candidates satisfying an
   emitted branch identically, and a candidate covered piecewise by the union of two branches
   (each covering one sign of a real symbol) while no single branch covers it, are all verdicted
   `OMITTED_BRANCH`. Required: the union test is computed over the branches DEFINED at the
   candidate — a branch undefined there is non-covering and is excluded from the union
   computation, never poison; the union decision object must cover the piecewise case (the
   product of one substituted constraint per defined branch vanishing identically on the
   candidate, or an equivalent union-variety membership — a per-branch containment OR is NOT
   sufficient); `COVERAGE_UNDECIDED` keeps its own token and is never collapsed to a decided
   one.
2. **A pole of the as-written equation is treated as a candidate component** (violates
   obligation 3's artifact exclusion). Candidates at which the as-written equation is undefined
   are kept and verdicted `OMITTED_BRANCH` (`s11_census_math.py:717-718` drops only `Set`
   objects). Measured: most of the round-2 WL omitted records carry only such candidates (the
   two legs' counts and the per-record evidence are in the twin). Required: a candidate where
   the as-written equation is undefined is an excluded artifact — excluded for that entailment,
   with the candidate and the undefined substitution printed on its per-record line; the
   record's completeness verdict is then computed over the surviving candidates (no surviving
   missing candidate ⇒ the covered token, with the exclusions itemized), never left
   `COVERAGE_UNDECIDED` for the exclusion alone. Where a record's equation is defined but every
   emitted branch is undefined at the candidate, defect 1's rule applies (all non-covering);
   where the equation itself is undefined, this exclusion wins.
3. **The exact sampler is non-generic: symbols in the same assumption pool collide** (violates
   obligation 5's refutation direction). `exact_sample_assignments`
   (`s11_census_math.py:794-831`) assigns `positive[(serial + index) % len(positive)]` with a
   pool of length 4, so two positive symbols at sorted-index distance ≡ 0 (mod 4) receive
   IDENTICAL values in EVERY sample; residuals carrying the coincidence factor of two such
   symbols sample to zero everywhere and are reported `UNDECIDED_ZERO_SAMPLES` although one
   generic point refutes them (measured population in the twin). Required: every sample
   assignment is generic — pairwise-distinct values within each assumption pool for all symbols
   simultaneously; the refutation direction must be able to fire on a residual vanishing only
   on a coincidence locus of ANY two sampled symbols.
4. **Witness premise semantics: classification must happen AFTER substitution** (supersedes the
   round-2 directive's class 5 wording; this round's directive originally partitioned by
   pre-substitution symbol binding, and both directive legs refuted that with real records).
   Required: substitute the witness into each premise conjunct FIRST (the WL `PREMISES` `And`
   is split into conjuncts — never evaluated as one opaque object), then classify each
   substituted conjunct: identically TRUE, identically FALSE, or contingent (free symbols
   remain and neither truth value is forced). Any FALSE conjunct ⇒ `WITNESS_FAILURE`. All
   conjuncts TRUE and membership contained on a coherent sheet ⇒ `WITNESS_VALIDATED`.
   Otherwise the contingent conjuncts are printed verbatim on the per-record line and the
   verdict is `WITNESS_UNDECIDED`. Assumption atoms on concrete numbers decide (a rational
   substituted into a realness/positivity atom yields TRUE or FALSE, in both dialects — the PY
   `Q.real(<rational>)` case must decide TRUE, not remain undecided). Membership-driven
   failures (witness point failing its own locus equations on all sheets) are premise-
   independent and must be unchanged. The round-2 premise-violating plant must still fail.
5. **The reducer counts sheet-level progress lines as objects** (violates obligation 6's count
   semantics). `CONTAINMENT_BRANCH_SHEET` / `CONTAINMENT_WITNESS_SHEET` lines carrying a
   limitation token are counted alongside the branch/witness verdict lines for the same objects
   (`s11_acceptance_reducer.py:257-273`; measured split in the twin). Required: reducer buckets
   count OBJECT-level verdict lines once; sheet-level lines are evidence for their parent
   verdict, not separately-counted objects. `failures` and `findings` semantics are already
   object-level and must not move.
6. **Premise evaluation is branch-unsound before substitution** (orchestrator-measured; the
   mechanism behind two false `premise_truth="FALSE"` witness failures). The parser expands
   `Element[X, Reals]` into re/im components of complex symbols BEFORE the witness is
   substituted (`parse_payload` with `assumption_free=True`); for radical/Abs-bearing `X` the
   expansion takes a non-principal branch, and two real witnesses are verdicted FALSE on an
   atom that is identically TRUE at the witness point (the `Sqrt` argument is exactly 0;
   records and commands in the twin). Required: premise atoms are evaluated by substituting the
   witness FIRST and then evaluating the principal value of the substituted atom; no
   pre-substitution re/im expansion may decide a premise's truth.

## Calibration extension (obligation 2 binds: production path, byte-shaped, must fail)

Keep every round-2 plant. Add, each demonstrated ABLE TO FAIL before the census runs:

- a candidate covered only by the UNION of two emitted branches (each covering one sign of a
  real symbol) with an additional sibling branch undefined at the candidate — must NOT be
  omitted (fails under per-branch OR and under the current product);
- a pole candidate zeroing the equation's denominator — must be excluded as an artifact, with
  the record's completeness verdict computed over the survivors;
- a residual vanishing exactly on the coincidence locus of two same-pool symbols — must be
  refuted by sample;
- a witness with a partially-bound premise that becomes identically FALSE after substitution —
  must FAIL (byte-shaped; a fully-bound violator does not exercise the partition);
- a byte-shaped witness whose premises include atoms still contingent after substitution
  (realness of unbound field symbols) alongside all-TRUE decided atoms — must be UNDECIDED
  with the contingent atoms printed;
- a byte-shaped witness on its locus whose every premise conjunct decides TRUE after
  substitution — must VALIDATE (fails at `89ed80c9`, where it lands UNDECIDED);
- a radical/Abs-bearing realness atom identically TRUE at the witness — must NOT fail (fails
  at `89ed80c9` via the pre-substitution expansion);
- a sheet-line/object-line counting case — object counted once;
- the `Abs` plant's argument replaced with one the parser cannot pre-fold, so the binding can
  fail.

## Order of work

1. Repairs + calibration extension. 2. Calibration FIRST, literal output (stop if any planted
case cannot fail). 3. Full census, both committed records (WL `a4cf6539`, PY `19591194`),
complete stdout preserved under `/home/trevnorris/.s11_build/census_build3/` (create it).
4. Reducer. 5. Build report: every measured claim carries its command and literal output, plus a
verdict-delta table against the round-2 summaries (the literal terminal lines in
`/home/trevnorris/.s11_build/census_build2/build_report.md` §5–§9).

## Constraints

As rounds 1–2: census findings are reported, never repaired; no engine, spec, committed-record,
or register modification; no commits; no Wolfram kernel; absolute paths; long + silent is
failure. Do not weaken anything the round-2 legs verified: the seven round-2 repair classes,
reducer failure/finding arithmetic, population reconciliation, the NONZERO refutation direction
where samples are generic, the `PROVED_CONSISTENT` decisions, and the two-route residual
comparisons all stay.
