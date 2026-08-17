# S11 locus-census instrument round — decision list

Targets: the two obligation-4 census instruments (round-2 brief, obligation 4), rebuilt as
repo-homed tools under `research/pde_ledger_v3/reduction/`:
an **independent containment/soundness census** (non-Wolfram engine) over every
`_EQUATIONS`/`_SOLUTION` locus pair of a committed record file, and the **undecided-record probe
census** over every undecided-class record. Companion measurements (generated):
`../_measurements/S11_locus_census_instrument_brief.md`. No engine, spec, or committed record is
modified this round; the seven DEFECT_REGISTER.md entries stay untouched — this round builds the
instrument that measures them.

## The measured wall (every claim twin-backed)

- The build-round instruments have never executed on real records: the SymPy containment harness
  parse-fails on all 484 locus pairs of the final D2/D4 records; the WL probe census dies on
  `ToExpression::sntxi` before its first probe and still exits 0. Both were calibrated only on
  planted synthetic records.
- Sheet handling is a measured trap: one review leg misjudged a radical branch as a non-solution
  (0/16 minors) that sheet-consistent evaluation decides is a solution (16/16 on both global
  sheets); a different branch genuinely fails 1/16 on both sheets. Naive per-minor evaluation
  gets radical loci wrong in both directions.
- `"NOT_APPLICABLE"` solution slots pair with `"NOT_APPLICABLE"` equations everywhere (178/178
  on D4; zero live-equations/NA-solution pairs on any record) — the spec's Q8b non-`VARIES` form.

## What must be true after the build

1. **Both instruments parse 100% of the record classes they audit** on the real committed records
   of BOTH engines — the WL canonical at `a4cf6539` and the SymPy record at `19591194`. A parse
   failure on any in-class line is an instrument failure, printed per-tag with the failing text;
   spec-sanctioned `NOT_APPLICABLE` pairs are classified, not errors.
2. **Calibration runs first and can fail.** Before any census: the probe decides a planted
   decidable-but-undecided record; the containment census flags a planted omitted branch AND a
   planted spurious branch in a planted solution list. Literal calibration output precedes census
   output in the deliverable.
3. **Verdicts are computed objects.** Per-locus lines carry the operands and the computed
   memberships/residuals, then the verdict token; nothing is asserted before it is printed;
   summary counts are computed from the per-locus lines. An instrument's exit code reflects
   execution, never verdict content.
4. **Radical loci are decided sheet-consistently**: every radicand-bearing system is evaluated
   under each global sheet assignment consistently across the whole system, and the per-locus
   line records which sheets were tested.
5. **Budgets expire honestly.** Per-locus budgets at least as large as the emitting engine's
   primary budget for that class; an expiry is an explicit per-locus token, never a crash, a
   silent skip, or an exit.
6. **The full census runs and is recorded**: both instruments over both engines' committed
   records, complete literal stdout preserved, with per-record and per-class summary counts.

## Out of scope

Any repair of what the census finds (the registered entries' round follows this one, with these
verdicts as its measured facts). Any change to `S11_SHARED_PHYSICS.md`. The count-payload
provenance question (register entry 3) — the census treats trailing provenance fields as
engine-local additions, per the register.

## Builder operational constraints

As round 2 (absolute paths, scratch outside the repo, observable progress). The containment
census must not require a Wolfram kernel; if the probe census needs one, one kernel at a time,
never in parallel with anything, and never on `XKIN_ANISO` 2 or 4 cell runs — the census reads
committed records, it does not re-run cells.
