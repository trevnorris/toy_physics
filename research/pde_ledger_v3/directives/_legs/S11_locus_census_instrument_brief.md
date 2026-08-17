# S11 locus-census instrument round — decision list (folded once after two legs)

Targets: the obligation-4 census instruments (round-2 brief, obligation 4), rebuilt as repo-homed
tools under `research/pde_ledger_v3/reduction/`: an **independent containment/soundness census**
over locus records of a committed record file, and the **undecided-record probe census**.
Companion measurements (generated): `../_measurements/S11_locus_census_instrument_brief.md`.
No engine, spec, or committed record is modified this round; the seven DEFECT_REGISTER.md entries
stay untouched — this round builds the instrument that measures them.

## The measured wall (every claim twin-backed)

- The build-round instruments have never executed on real records: the SymPy containment harness
  parse-fails on all 484 locus pairs of the final D2/D4 records; the WL probe census dies on
  `ToExpression::sntxi` before its first probe and still exits 0. Both were calibrated only on
  planted synthetic records that never took the production parse path.
- The two committed records use different payload dialects (WL `InputForm`, PY `srepr`), and the
  undecided token appears in several forms per dialect; the twin's per-class counts are the
  population. A WL-shaped scan reports zero PY undecideds against hundreds present — this exact
  mistake is in the twin's own first generation, caught by both legs.
- Sheet handling is a measured trap in both directions: a naive evaluation misjudged a radical
  branch as a non-solution that sheet-consistent evaluation decides is a solution (16/16 minors on
  both global sheets); a different branch of the same list genuinely fails 1 minor on both sheets.
  Solution branches themselves introduce radicals the equations do not carry.
- The PY producer has no resource budgets anywhere; the WL producer's `canonicalLocus` route is an
  unbounded `GroebnerBasis`. "At least the emitting budget" is undefined for these classes.
- `"NOT_APPLICABLE"` solution slots pair with `"NOT_APPLICABLE"` equations everywhere (392/392 on
  the WL canonical; zero live-equations/NA-solution pairs on any record) — the spec's Q8b
  non-`VARIES` form.

## The census population — defined here, not by the implementation

- **Containment/soundness census**: every §5 locus pair on BOTH committed records — a
  `_SOLUTION`-suffixed tag with its same-stem `_EQUATIONS` sibling — including live
  `*_CHANGE_LOCUS_*` and `KW_ZERO_LOCUS_*` pairs, and excluding exactly two named non-§5 classes
  (`C1_*`, `STRATUM<i>_DEFINING_*`) and the `DIM_*` dimension-solve pairs. PLUS every
  `_REAL_WITNESS` whose sibling `_REAL_STATUS` is `PROVED_NONEMPTY`: the witness point is
  substituted into its `_EQUATIONS` and the §3 premises, exactly, and the residuals are census
  output. The instrument enumerates every raw `_SOLUTION`/`_EQUATIONS`/`_REAL_WITNESS`-suffixed
  line of both files and reconciles each one: in-population, exact paired-`NOT_APPLICABLE`
  sentinel, or one of the named excluded classes — an unreconciled line is an instrument failure,
  printed with its tag.
- **Probe census**: every `UNDECIDED` occurrence in any payload of both records, in every dialect
  form (WL: `STATUS_TOKEN`-keyed, `SIGN_TOKEN`-keyed, and bare; PY: the `srepr` `Str(...)`
  equivalents of all three, and nested forms). The dialect-aware enumeration prints per-class
  counts and reconciles them against a raw whole-file `UNDECIDED` count; an unexplained gap is an
  instrument failure.

## What must be true after the build

1. **Parse 100% of the population.** Only the exact paired-`NOT_APPLICABLE` sentinel bypasses
   semantic parsing; every other in-population payload is parsed to objects in both dialects —
   including attempt-provenance-bearing solution lists (trailing additive fields are stripped for
   the membership computation, never grounds to classify a record opaque), multi-megabyte PY
   lines, and `ConditionSet`/`FiniteSet`/`EmptySet` returns. A parse failure on any in-population
   line is an instrument failure, printed per-tag with the failing text.
2. **One path, then calibration first.** Calibration and census share the same entrypoint,
   parser, enumerator, branch and sheet machinery, solver route, and verdict reducer. Planted
   records are byte-shaped copies of real lines from BOTH dialects — among them an
   attempt-provenance solution, a `ConditionSet` line, and a paired-`NOT_APPLICABLE` pair — with
   one defect introduced per planted case: a planted omitted branch, a planted spurious branch,
   and a planted decidable-but-undecided record. Every planted defect must be detected, before
   any census, through the production path; literal calibration output precedes census output.
3. **The containment solve is independent of the producer's route.** The PY producer solves with
   `sp.solve` (engine:1005); the census containment route must be structurally different (its
   choice is the builder's, documented in code), and containment/soundness is decided up to
   algebraic equivalence. Finite generic sampling is admissible only as a fallback and its
   verdict token names itself as sampled.
4. **Live operands and residual recompute (parent obligation restored).** Undecided-record probes
   attack the record's own emitted operands; emitted residual-class payloads are recomputed from
   those operands and compared; an absent, unparseable, or stale operand field is a census
   finding with its own token, never a skip.
5. **Branch × sheet granularity.** For every exposed branch of every in-population `_SOLUTION`:
   substitute the branch as written (radicals inside the branch are part of the emitted object
   and are never sign-flipped), then evaluate every equation under every global sheet assignment
   of the remaining radicands, recording per-equation residual and definedness status, then the
   per-branch verdict. No any-branch or whole-list reduction may hide a per-branch verdict. Where
   independent radicands or `Abs` make full enumeration exceed a stated bound, the branch carries
   an explicit sheet-incompleteness token with the count of untested assignments.
6. **Verdicts are computed objects, and a named reducer turns them into the round verdict.**
   Per-record lines carry operands, computed memberships/residuals, then the verdict token;
   nothing is asserted before it is printed; instrument exit codes reflect execution only. A
   separate named acceptance reducer reads the full census stdout and computes the round verdict:
   any calibration miss, unreconciled line, in-population parse failure, decided
   undecided-record, or residual mismatch fails the round; census verdicts on real loci
   (spurious/omitted branches, witness failures) are FINDINGS the reducer counts and lists,
   never silently folded into a pass/fail bit without their per-record lines. The reducer's
   literal output is part of the deliverable.
7. **Budgets are explicit, floored, and honest.** Each audit operation's time and memory budget
   is set in code where used: at least the emitting engine's primary budget where one exists
   (twin lists the WL constants), and at least 60 s / 512 MB per record for classes whose
   emitting route carries no budget (all PY routes; the WL `GroebnerBasis` route). An expiry is
   an explicit per-record token, never a crash, skip, or exit; per-class expiry counts appear in
   the summary.
8. **The full census runs and is recorded**: both instruments over both engines' committed
   records (WL at `a4cf6539`, PY at `19591194`), complete literal stdout preserved, per-record
   and per-class summary counts computed from the per-record lines, and the reducer verdict.

## Out of scope

Any repair of what the census finds (the registered entries' round follows, taking these verdicts
as measured facts — subject to obligation 6's rule that findings are carried as per-record lines,
not summary bits). Any change to `S11_SHARED_PHYSICS.md` or either engine. The count-payload
provenance question (register entry 3) — trailing provenance fields are engine-local additions,
per the register.

## Builder operational constraints

As round 2 (absolute paths, scratch outside the repository, observable progress). The containment
census must not require a Wolfram kernel; if the probe census needs one, one kernel at a time,
never in parallel with anything, and never re-running `XKIN_ANISO` cells — the census reads
committed records, it does not run cells.
