# Repair directive — S11 census instruments, round 3

You are repairing the instrument modules under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` at commit `89ed80c9`, governed by
the folded brief (`S11_locus_census_instrument_brief.md`, `fa8c58b3`) — its eight obligations and
population section remain the acceptance surface. Two independent review legs plus orchestrator
re-computation verified the round-2 repairs (all seven classes live on the production path, one
ablation diff each) and measured the five defects below in the round-2 census. Repair them all,
extend the calibration, and re-run the full census. Everything the legs certified live must stay
live — the round-3 calibration must still catch a regression of every round-2 class.

## Verified defects

1. **The union-membership test lets one sibling branch poison coverage** (violates obligation 3).
   A candidate is covered when it lies in the union of emitted branches — i.e. when it satisfies
   the FULL constraint set of AT LEAST ONE emitted branch. The current test multiplies one
   constraint from EVERY emitted branch (`s11_census_math.py:683-715`), so a sibling branch that
   is undefined at the candidate (a rational chart with a pole there) makes the product
   undefined; the sampling fallback discards any sample where ANY sibling is undefined
   (`:639-680`); and `COVERAGE_UNDECIDED` is then collapsed to `NOT_COVERED` (`:713`), which the
   completeness pass promotes to `OMITTED_BRANCH` (`:773-774`). Measured: a candidate satisfying
   an emitted branch's constraint identically (simplified residual 0) is verdicted
   `OMITTED_BRANCH`; 7 of the 70 round-2 PY omitted lines carry this shape. Required: membership
   is decided per emitted branch — a branch undefined at the candidate is non-covering for that
   branch, never poison for the union; the union verdict is the OR over branches;
   `COVERAGE_UNDECIDED` keeps its own token and is never collapsed to a decided one.
2. **A pole of the as-written equation is treated as a candidate component** (violates
   obligation 3's artifact exclusion). Candidates at which the equation as written is undefined
   (the cleared-denominator points, e.g. a candidate zeroing a factor of the equation's
   denominator) are currently kept and verdicted `OMITTED_BRANCH`
   (`s11_census_math.py:717-718` drops only `Set` objects). Measured: 4 of the 13 round-2 WL
   omitted records are exactly this (candidate zeroes the `_EQUATIONS` denominator; substitution
   yields an undefined value, not a solution). Required: a candidate where the as-written
   equation is undefined is an excluded artifact — excluded for that entailment, with the
   candidate and the undefined substitution printed on its per-record line.
3. **The exact sampler is non-generic: symbols in the same assumption pool collide** (violates
   obligation 5's refutation direction). `exact_sample_assignments`
   (`s11_census_math.py:794-831`) assigns `positive[(serial + index) % len(positive)]` with a
   pool of length 4, so two positive symbols whose sorted-index distance is ≡ 0 (mod 4) receive
   IDENTICAL values in EVERY sample. Measured: for the free-symbol set
   `(bComp, k1, k2, k3, muR)` every sample lies on `bComp = muR`, so a residual carrying the
   factor `(bComp − muR)` samples to zero in all samples and is reported
   `UNDECIDED_ZERO_SAMPLES` although one generic point refutes it; 28 WL probe records are
   affected. Required: every sample assignment is generic — pairwise-distinct values within each
   assumption pool for all symbols simultaneously (a pool at least as large as the symbol count,
   or per-symbol distinct constructions); the refutation direction must be able to fire on a
   residual vanishing only on a coincidence locus of ANY two sampled symbols.
4. **Witness premise semantics leave the validated direction unreachable** (round-2 directive
   defect, superseding its class 5 wording). Every real `_REAL_WITNESS` binds only its solve
   variables, so premises over unbound symbols (`k1`, `k2`, `A`, `a1`, `a2`, …) evaluate
   UNDECIDED and every witness lands in `WITNESS_UNDECIDED`: round-2 emitted ZERO
   `WITNESS_VALIDATED` across both records (round 1: 107 + 99). Required semantics: partition
   the §3 premises by whether their free symbols are all bound by the witness; evaluate the
   bound ones conjuncted at the witness point on assumption-free symbols; print the unbound
   premises verbatim on the per-record line as unevaluated. A witness VALIDATES when its
   equations vanish on a coherent sheet AND every bound premise is TRUE; it FAILS when it fails
   membership on all sheets or any bound premise is FALSE; it is UNDECIDED only when membership
   or a bound premise is genuinely undecidable. The round-2 failure direction (premise-violating
   plant fails; the 172 measured witness failures) must be unchanged.
5. **The reducer counts sheet-level progress lines as objects** (violates obligation 6's
   count semantics). `CONTAINMENT_BRANCH_SHEET` / `CONTAINMENT_WITNESS_SHEET` lines carrying a
   limitation token are counted alongside the branch/witness verdict lines for the same objects
   (`s11_acceptance_reducer.py:257-273`): round-2 `limitations=885` counts 78 sheet-progress
   lines whose parent verdict lines are already counted (object-level count: 807). Required:
   reducer buckets count OBJECT-level verdict lines once; sheet-level lines are evidence for
   their parent verdict, not separately-counted objects. `failures` and `findings` semantics are
   already object-level and must not move.

## Calibration extension (obligation 2 binds: production path, byte-shaped, must fail)

Keep every round-2 plant. Add, each demonstrated ABLE TO FAIL before the census runs: a covered
chart candidate with a polar sibling branch (must NOT be omitted); a pole candidate zeroing the
equation denominator (must be excluded as an artifact, not omitted); a residual vanishing
exactly on the coincidence locus of two same-pool symbols (must be refuted by sample); a
premise-satisfying witness on its locus (must VALIDATE) alongside the round-2 premise-violating
plant (must still FAIL); a sheet-line/object-line counting case (object counted once). Replace
the `Abs` plant's constant argument with one the parser cannot pre-fold, so the binding can
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
