# Repair directive — S11 census instruments, round 4 (single defect class)

You are repairing the instrument modules under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` at commit `fd9a5835`, governed by
the folded brief (`S11_locus_census_instrument_brief.md`, `fa8c58b3`). Rounds 2–3 repaired
thirteen verified defect classes; both round-3 review legs certified all six round-3 repairs
live and every other verdict change sound, and converged on exactly ONE surviving defect,
orchestrator-reproduced (commands and literal outputs in this directive's `../_measurements/`
twin). Repair it, extend the calibration, and re-run the full census. Everything the round-3
legs certified must stay live.

## The verified defect

**A constant that fails to simplify is classified NONZERO, and sampled coverage treats that as
a refutation** (violates obligation 5's "a sample refutes, never proves", relocated into the
zero-status decision itself). `simplify_residual` (`s11_census_math.py:544-578`) classifies a
free-symbol-less expression NONZERO by the structural test `simplified != 0` (`:550`, `:570`).
Measured: an exact algebraic zero (a nested-radical constant arising when a covering branch's
constraint is evaluated at a sampler point; 50-digit value ≈ 1e-161) is classified NONZERO, so
`_sample_union_coverage` (`:689-717`) sees the covering branch as non-covering, the point as
uncovered, and emits `NOT_COVERED_SAMPLED`; completeness promotes it to `OMITTED_BRANCH`
(`:829`). One committed-record verdict is false by this route (the record named in the twin —
both its candidates lie identically in the union of its emitted branches; the symbolic union
product exceeds the simplify size cap, so the sampled fallback decided). Required:

1. **Zero-status of a constant is decided by proof, never by structure — and the exact route
   is mandatory, not one option.** For a free-symbol-less algebraic expression the
   classification is: ZERO only with an exact zero certificate (the minimal-polynomial /
   exact equals-zero decision on algebraic numbers); NONZERO only with an exact nonzero
   certificate or a rigorous numeric enclosure that excludes zero; UNDECIDED only after the
   exact route has been attempted and failed (resource-guarded, with the guard's expiry
   recorded) — never NONZERO by failure-to-simplify, and never UNDECIDED merely because a
   fixed working precision could not separate a nonzero value from zero when the exact route
   decides it. This binds every consumer of the status (coverage, membership, spurious
   refutation, probe refutation, definedness). Both directive legs measured the
   under-specified version: two compliant implementations (fixed-precision interval vs
   minimal-polynomial) returned different verdicts on the same branch — the exact-route
   mandate is what removes that freedom.
2. **Sampled coverage refutes only on a certified-uncovered point, under per-branch AND
   semantics.** At a candidate sample point: a DEFINED branch covers the point iff EVERY one
   of its substituted constraints is certified ZERO under rule 1; the point is
   certified-uncovered iff EVERY defined branch has at least one constraint certified NONZERO
   there (equivalently: some one-constraint-per-defined-branch combination product is
   certified NONZERO); any other status pattern leaves the point undecided. ⛔ Constraints
   belonging to ONE branch are never multiplied together — that turns the branch's AND into
   an OR and manufactures coverage (both directive legs measured a real candidate this
   mis-covers). A sample point that is certified-uncovered refutes union coverage
   (`NOT_COVERED_SAMPLED`, the witness point printed); if every sampled point is covered or
   undecided, the coverage verdict is `COVERAGE_UNDECIDED` (its own token, never collapsed) —
   a sample refutes coverage, it never proves non-coverage.

The genuine refutation direction must stay: a point where the union product is certified
nonzero is a witness of non-coverage and must still be verdicted `NOT_COVERED_SAMPLED`; sampled
NONZERO classifications backed by certified-nonzero values (the round-3 spurious and probe
refutations both legs re-decided and upheld) must not change.

## Calibration extension (obligation 2 binds: production path, byte-shaped, must fail)

Keep every round-2/3 plant. Add, each demonstrated ABLE TO FAIL before the census runs:

- a free-symbol-less nested-radical exact zero routed through the production status
  classification — must NOT be NONZERO (fails at `fd9a5835`);
- a free-symbol-less algebraic constant that is nonzero but smaller than any fixed working
  precision (shape: `Sqrt[10^200 + 1] - 10^100`) — must be certified NONZERO by the exact
  route (fails under a fixed-precision-only implementation);
- a solution/equation pair whose candidates are inverse charts of the emitted branches with a
  union product too large for the symbolic simplify route — must NOT be omitted (fails at
  `fd9a5835`);
- a genuinely uncovered candidate whose union product is certified nonzero at a sampled
  point — must still be `NOT_COVERED_SAMPLED` (guards the refutation direction).

## Order of work

1. Repair + calibration extension. 2. Calibration FIRST, literal output (stop if any planted
case cannot fail). 3. Full census, both committed records (WL `a4cf6539`, PY `19591194`),
complete stdout preserved under `/home/trevnorris/.s11_build/census_build4/` (create it).
4. Reducer. 5. Build report: every measured claim carries its command and literal output, plus a
verdict-delta table against the round-3 summaries (the literal terminal lines in
`/home/trevnorris/.s11_build/census_build3/build_report.md` §5–§9).

## Constraints

As rounds 1–3: census findings are reported, never repaired; no engine, spec, committed-record,
or register modification; no commits; no Wolfram kernel; absolute paths; long + silent is
failure. Do not weaken anything the round-3 legs certified: the six round-3 classes, the seven
round-2 classes, reducer object-level arithmetic, population reconciliation, certified-nonzero
refutations, `PROVED_CONSISTENT`, and the two-route residual comparisons all stay.
