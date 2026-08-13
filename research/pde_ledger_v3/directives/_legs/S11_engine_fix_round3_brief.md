# S11 SymPy engine — fix round 3. ⭐ Three items. ⛔ Two block; the third unblocks the publish.

## Authority and boundary

Edit `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` **only**. ⛔ Change no other
file. `CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` is the sole physics
authority. `research/pde_ledger_v3/directives/S11_sympy_build_directive.md` is **closed**.

⭐ **Baseline is `4d5ff0f6`.** Round 2 succeeded at what it targeted: `MAIN` D2–D5 completes (~295 s), the
export digest path resolves, a publish failure is attributable to the publish, the nullspace basis is a
computed object again, and the `ISSUES` truncation is gone. ⛔ Do not revisit any of that.
⚠ One round-2 fix left a residue that is this round's item 1: the hard-coded size gates were removed and
**the class survived** — the refusal moved from a count to another attempt-free pre-check.

Every claim below was measured by the orchestrator at HEAD `784ec815` and re-measured independently by
two review legs; commands and literal output: `directives/_measurements/S11_engine_fix_round3_brief.md`.
Probes: `~/.s11_build/round3_probes/`, `/tmp/s11r3_leg_Wlem/`, `/tmp/s11r3_leg_xfHB/`.

## ⛔⛔ 1 · BLOCKING — "MEASURED UNAVAILABLE" IS RECORDED WITHOUT A MEASUREMENT

`:820-828` decides that `ROOT_COINCIDENCE_*_COEFF` solves are unavailable from a pre-check of the input
(`compound_radical_present`, `:763-770`), emits a `ConditionSet`, and records *"measured unavailable"* —
⛔ **no CAS call was attempted**. Measured at `XKIN_ANISO` D2: the pre-check answers in <0.001 s; an
independent solve of the engine's own recorded system **terminates with real branches in ~9–10 s**
(three independent measurements). `XKIN_ANISO` D3 emits **three** sibling refusals of the same class
before its wall.
⛔ Against the physics authority: `S11_SHARED_PHYSICS.md:245` requires `_SOLUTION` be *the solution set
exactly as your CAS returns it*.
⇒ ⭐ **The only route to an "unavailable" record is an attempt** — the record carries what the attempt
returned or raised, under §10. ⛔ A property of the input is not evidence about what an attempt will do
(round 2 item 7 stands: if some attempted solve genuinely does not return, that is *measured*, then
reported under §10). ⭐ The property holds at **every** emission site of the class, in every package at
every dimension — ⛔ not at a sampled one. ⚠ This lengthens `XKIN_ANISO` cells; item 3 exists so the
publish does not hang on it — ⛔ cost is never a reason to narrow this (CLAUDE.md rule 11).

## ⛔⛔ 2 · BLOCKING — STATUS TOKENS DECIDED BY NOTHING, AND ONE TEST READ OFF THE WRONG OBJECT

**(a)** `evaluate_premise` (`:692-711`) returns `None` — hence `UNDECIDED` — for premises **decided
under the branch, in both directions and by both routes**: a relational whose substitution
auto-evaluates to a Boolean (measured: `True` and `False` both come back `None`), and an
`AppliedPredicate` whose CAS attribute decides it **false** (measured: `Q.real(k1)` at `k1=I` has
`is_real=False` and still returns `None`). A `_REAL_ADMISSIBLE` entry whose `TEST_OBJECT` carries a
literal `False` sits beside `STATUS_TOKEN: UNDECIDED`; distribution at `MAIN` D2/D3/D4 is
`{UNDECIDED: 22/35/39}` with **zero ADMISSIBLE and zero EXCLUDED anywhere measured**, `XKIN_ANISO` D2
included. `S11_SHARED_PHYSICS.md:248` and `:266` define the entry; `:657` defines a **stratum** by that
token — the §Q8b build consumes what this defect currently empties.
⇒ ⭐ **A premise decided under the branch carries its decision's token, whichever route decided it and
in either direction**; `UNDECIDED` is reserved for a premise the CAS genuinely cannot decide.

**(b)** `:870` builds `_INCONSISTENT`'s test object as `sp.Eq(solution_payload, sp.Tuple(),
evaluate=False)` — ⛔ **read off the solver's return**. `S11_SHARED_PHYSICS.md:250-251` carries the full
obligation: the test is *computed from `_EQUATIONS`, ⛔ never read off the solver's empty token* — and
`:246-247` define both it and `_IDENTICALLY_SATISFIED`. Two measured consequences at HEAD: the token can
never reach `PROVED_TRUE` (the unevaluated `Equality` matches nothing in `symbolic_truth_status`
`:640-647`); and the weakest agreement-only edit — letting that `Eq` evaluate — types a locus
`PROVED_TRUE` (no solution at all) **where an exact premise-satisfying assignment zeros every equation**,
because a CAS solve returning `[]` is not a proof of emptiness. ⇒ a wrong degenerate-case split feeding
`:657` strata.
⇒ ⭐ **The `_INCONSISTENT` and `_IDENTICALLY_SATISFIED` tokens are decided from `_EQUATIONS`** — the
`TEST_OBJECT`/`OPERANDS` they carry are objects derived from the equations, ⛔ never the `_SOLUTION`
payload or its shape — **and the token must agree with what that object evaluates to.** An exact
assignment satisfying the premises and zeroing every equation is dispositive against `PROVED_TRUE`.
⛔ Do not specify the routine (rule 3).

## ⛔ 3 · THE PUBLISH WAITS ON CELLS ITS OWN GATE DOES NOT READ

`write_exports` runs once, after the full `PACKAGE_ORDER` loop (`:1866-1870`), while its gate reads only
`MAIN` (`:1816-1818`) and `MAIN` runs first (`:136-147`). Measured: `MAIN` D2/D3/D4 = 2.6/7.1/30.1 s;
`XKIN_ANISO` D2 = 139–141 s; `XKIN_ANISO` D3 **exceeds a 600 s cap at HEAD** (leg-measured this round:
probe killed at exactly 600 s), scaling 4–5×/dim — and item 1 adds the solves the gate was suppressing.
⇒ a sweep must survive every package before `MAIN`'s own results can publish.
⚠ **The spec pins the post-sweep record**: `S11_SHARED_PHYSICS.md:1038` — the run record is *observed,
not declared*; `RUN_PAIRS` / `SKIPPED_PAIRS` are emitted **after** the sweep with
`SKIPPED = declared ∖ completed`. ⇒ those two objects cannot move mid-loop. What a mid-loop publish may
carry is a **publish-time record that cannot be mistaken for them**, in which a cell not yet attempted
is stated as such — ⛔ never claimed `skipped`, failed, or run. The 18-false-rows export round 2 removed
is the measured instance of getting this wrong; the *"cell skipped after"* issue text (`:1859`) for an
attempted-and-failed cell is the conflation's other half.
⇒ ⭐ **A publish whose gate is satisfied does not wait on cells the gate does not read**, and ⭐ **every
record the publish emits is true at the moment it is emitted**: publish-time state distinguishable from
the spec's post-sweep objects, no row claiming a run record it does not have. ⭐ Both properties, or
neither: ⛔ if they cannot both hold, keep end-of-run publication and record the conflict under §10.
⭐ Round-2 properties survive any placement change: a publish failure is attributable to the publish
(`:1867-1877`), and a failed run still emits its §10 report (`:1901-1912`).

## Boundaries

- ⛔ No memory cap, no timeout, no handler that swallows a failure to make a run finish.
- ⛔ Do not change `PACKAGE_DIMS`, the `D` range, or any package.
- ⛔⛔ **§Q8b machinery is out of scope** — `stratum_candidates`, the `STRATUM_ORDERING` emission, and
  every §Q8b evidence object belong to their own build directive. Items 1–2 repair what feeds them;
  ⛔ do not start the build here.
- ⭐ The objects §§4–8 name must not change. **What the engine computes for them must change where items
  1–2 say it is currently refused or mistyped.**
- ⛔ **No expected value and no acceptance criterion referencing one** (rule 5). ⛔ Do not treat any prior
  output as a reference — ⭐ a changed value is a **finding to report under §10**, ⛔ never something to
  tune away.
- ⛔ Do not add a registry, completeness map or exit policy. ⚠ This run's own observed bookkeeping
  (`completed_pairs` and the like) is not a registry; item 3's publish-time record reads it.

## Acceptance

⚠ Checked against the items: with any item unfixed, the matching demonstration below cannot be produced.
1. `/home/trevnorris/.s11_build/repro_d5.py` runs `run_cell('MAIN', 5)` to completion; literal stdout.
2. **Item 1:** run the single `XKIN_ANISO` D2 cell (~140 s at baseline) **and** the single `XKIN_ANISO`
   D3 cell under `timeout 600`, reporting each literal stream (D3's partial stream, if the cap hits, is
   the report). ⛔ **The criterion is class-wide over every stream this acceptance produces**: any
   `ROOT_COINCIDENCE_*_COEFF` record of unavailability produced without an attempt — whatever its token,
   in any package at any dimension — fails this item. An emission carries the provenance of an attempt:
   what the CAS returned, or the raised failure under §10.
3. **Item 2:** from the same runs, report for `MAIN` D2 every `_REAL_ADMISSIBLE` entry's `STATUS_TOKEN`
   beside its `TEST_OBJECT`, and every `_INCONSISTENT` and `_IDENTICALLY_SATISFIED` `STATUS_TOKEN`
   beside its `TEST_OBJECT` **and `OPERANDS`**. ⛔ Fails on any of: a premise decided under the branch —
   by either route, in either direction — beside `UNDECIDED`; a test object that cannot evaluate by
   construction; a `_INCONSISTENT`/`_IDENTICALLY_SATISFIED` token whose `TEST_OBJECT` or `OPERANDS`
   derive from the `_SOLUTION` payload rather than from `_EQUATIONS`.
4. **Item 3, closed path:** on a reduced `PACKAGE_DIMS` copy under `/tmp` — ⛔ never by editing the
   engine's own declarations — demonstrate with literal stdout: (i) a run whose declared `MAIN` completes
   publishes **before later packages run**, and the publish-time record states not-yet-attempted cells as
   such, distinguishable from the spec's post-sweep objects; (ii) the same run still emits `RUN_PAIRS` /
   `SKIPPED_PAIRS` **after** the sweep per `S11_SHARED_PHYSICS.md:1038`, and a later cell that raises
   appears there and under §10 — with no published row contradicted by it; (iii) a run whose `MAIN` does
   not complete does not publish and §10 says why; (iv) a publish failure is recorded under the publish,
   not under any cell. **Hatch path:** if the two item-3 properties cannot both hold, the acceptance for
   this item is instead the §10 entry recording the conflict, with end-of-run publication retained — ⛔
   and then no publish-time demonstration is owed, but (iii) and (iv) still are.
5. ⛔ Do not run the full package loop; it has OOM-killed this machine. ⭐ Any engine-cell run happens
   with `/tmp/s11_watchdog.sh` armed (capture its pid; kill that pid when the run ends — ⛔ never
   `pkill` by name).

## Deliverable

The edited script, and a note saying per item what changed, plus every value that moved anywhere and
where you reported it.
