# Repair directive — S11 census instruments, round 2

You are repairing the six instrument modules under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` at commit `90ab5e2d`, governed by
the same folded brief (`S11_locus_census_instrument_brief.md`, commit `fa8c58b3`) — its eight
obligations and population section still ARE the acceptance surface. Two independent review legs
plus orchestrator re-computation measured the defects below in the first real census; every one
was verified by literal output. Repair them all, extend the calibration to cover each class, and
re-run the full census.

## Verified defects (each violates a numbered brief obligation)

1. **Sheet semantics manufacture SPURIOUS verdicts** (violates obligation 5's "substitute the
   branch as written"). The instrument masks branch-internal radicals, flips equation-side
   occurrences of the SAME radicand independently, and lets any failing sheet win
   (`s11_census_math.py:349-393`, `s11_containment_census.py:139-151`). Measured: KW_ZERO
   branches whose as-written residual is identically zero are verdicted SPURIOUS (72 mixed-sheet
   cases; 10/10 sampled spurious verdicts failed independent re-decision). Required: substitute
   the branch as written FIRST, simplify, and only radicands still present in the simplified
   residual are sheet-swept; a branch-embedded radical is never flipped; the per-branch verdict
   follows from the per-sheet lines it prints. Numeric radicands (e.g. `Sqrt[2]`) are constants,
   not sheets.
2. **Coverage is exact-match against single branches, not variety containment** (violates
   obligation 3's "up to algebraic equivalence"). Chart-equivalent candidates (`k1 = ±Sqrt[-k2^2]`
   vs emitted `k2 = ∓I k1` — same variety), `Abs`-restricted rewrites, and `EmptySet` solver
   artifacts are all verdicted OMITTED (`s11_census_math.py:524-537, 559-567`; ~294 of 355
   omitted verdicts match these false patterns). Required: a candidate is covered when it lies in
   the UNION of emitted branches up to algebraic equivalence; a non-point artifact is not a
   candidate; sampling fallback keeps its explicit token.
3. **The probe treats a solver's empty return as a proof** (violates the spec's own opaque-object
   rule, which the brief binds). `nonlinsolve → EmptySet` becomes `PROVED_INCONSISTENT`
   (`s11_undecided_probe_census.py:258-261`); orchestrator-verified counterexample: a
   PROVED_INCONSISTENT system with exact solution point `k1=0, k2=0` (all residuals zero).
   Required: an empty solver return without an independent emptiness certificate leaves the
   record UNDECIDED-CONFIRMED, with the return recorded.
4. **Finite zero-sampling promoted to a proof of identical vanishing.** 28 WL identity records
   were decided by 4-point zero samples (`s11_undecided_probe_census.py:135-136`,
   `s11_census_math.py:435-436`). A sample refutes; it never proves. Required: zero-sample
   agreement leaves the record UNDECIDED-CONFIRMED unless a symbolic proof decides it; the
   NONZERO direction (refutation by sample) stays.
5. **The PY witness premise check is unreachable** (violates the population section's "substituted
   into its `_EQUATIONS` and the §3 premises, exactly"). Premises arrive as a Tuple never
   conjuncted, and positivity premises auto-evaluate True because parsed symbols carry
   `positive=True` (`s11_containment_census.py:203-211, 283-284`); a planted premise-violating
   witness (`B_comp = -1` under `B_comp > 0`) passes WITNESS_VALIDATED. Required: premises
   conjuncted and evaluated at the witness point on assumption-free symbols; the planted case
   must fail.
6. **Parser gaps leave whole sectors unmeasured** (violates obligation 1): the empty branch
   `{{}}` (62 WL pairs); Boolean premises containing `Unequal` (246 + 446 records — every XKIN
   witness and most WL admissibility probes); `Abs` parsed as an inert Function so `Abs[0]` is
   classified NONZERO (≥2 of 46 witness failures false); 160 PY probe operands rejected by a
   single-assignment branch-unwrapping bug (`s11_undecided_probe_census.py:77-95`). Repair all
   four; `parse_wl` must bind `Abs` (and audit its whole head-binding table the same way).
7. **The reducer's verdict-token taxonomy is open** (violates obligation 6). Tokens
   `NON_VERDICT_TEXT`, `BRANCH_MEMBERSHIP_UNDECIDED`, `COMPLETENESS_UNDECIDED`,
   `WITNESS_UNDECIDED`, `WITNESS_SHEET_INCOMPLETE`, `WITNESS_VALIDATED_SAMPLED` fall into no
   failure/finding/limitation bucket (`s11_acceptance_reducer.py:23-38, 170-179`). Required: a
   closed taxonomy — every token the censuses can emit maps to exactly one bucket, and an
   unrecognized token is itself an acceptance failure.

## Calibration extension (obligation 2 still binds: production path, byte-shaped, must fail)

Add planted cases covering every class above, each demonstrated to FAIL before the census runs:
a radical-bearing branch whose as-written residual is zero (must NOT be spurious) and one that
genuinely fails all coherent sheets (must be spurious); a chart-rewritten covered candidate (must
NOT be omitted) and a genuinely omitted component (must be); an EmptySet-return inconsistency
(must stay undecided-confirmed); a premise-violating witness (must fail); a `{{}}` branch, an
`Unequal`-bearing premise, an `Abs`-bearing residual, and a single-assignment PY branch (all must
parse and verdict correctly).

## Order of work

1. Repairs + calibration extension. 2. Calibration FIRST, literal output (stop if any planted
case cannot fail). 3. Full census, both committed records (WL `a4cf6539`, PY `19591194`),
complete stdout preserved under `/home/trevnorris/.s11_build/census_build2/` (create it).
4. Reducer. 5. Build report: every measured claim carries its command and literal output, plus a
verdict-delta table against the round-1 census summaries (the six summary lines in
`/home/trevnorris/.s11_build/census_build/build_report.md` §4–§8).

## Constraints

As the build directive (`07ab921d`): census findings are reported, never repaired; no engine,
spec, committed-record, or register modification; no commits; no Wolfram kernel for containment;
absolute paths; long + silent is failure. Do not weaken any verdict logic the legs verified
correct: reducer arithmetic, population reconciliation, the NONZERO refutation direction, the
PROVED_CONSISTENT decisions, and the two-route residual comparisons all stay.
