# Measurements behind `S11_engine_fix_round3_brief.md`

Commands and literal output. CLAUDE.md rule 2. Regenerate; do not transcribe.
⛔ NO S11 RESULT VALUE APPEARS HERE — rule 5. This file carries defects, locations, timings and
status-token distributions, because a builder reads it. One elision is marked where a probe printed
a solution payload; the full stdout is at `~/.s11_build/round3_probes/`.
All measurements taken at HEAD `784ec815` (engine state `4d5ff0f6`), 2026-08-13, by the orchestrator —
⛔ not transcribed from a leg. Probe scripts: `~/.s11_build/round3_probes/probe_r3_pure.py`,
`~/.s11_build/round3_probes/probe_r3_cells.py` (interceptors wrap `emit_quantity`/`emit_locus` and
always call the originals unchanged; per-cell tag streams beside the probes).

## 1 · D1 — the gate refuses a solve that terminates (BLOCKING)

The gate:
```
$ sed -n '819,828p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    solve_skipped = False
    if (base_quantity.startswith("ROOT_COINCIDENCE") and base_quantity.endswith("_COEFF")
            and compound_radical_present(residual_exprs)):
        solution = locus_conditionset(equations, variables)
        solve_skipped = True
        record_issue(
            f"{base_quantity}: multivariate radical coefficient solve measured unavailable; "
            "emitted exact ConditionSet for the unresolved locus",
            package=package, n=n, root=root, stratum=stratum,
        )
```
`compound_radical_present` (`:763-770`) is a pre-order syntactic walk; no CAS call precedes the
"measured unavailable" record.

Where it fires, and the refused solve attempted independently on the engine's own recorded system:
```
$ cd ~/.s11_build/round3_probes && python3 probe_r3_cells.py \
    /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py XKIN_ANISO 2
== CELL XKIN_ANISO D2 -> /home/trevnorris/.s11_build/round3_probes/r3_XKIN_ANISO_D2.out
CELL_RETURNED True  CELL_SECONDS 141.2
REAL_ADMISSIBLE_STATUS_TOKEN_COUNTS {'UNDECIDED': 19}
ROOT_COINCIDENCE_COEFF_CALLS 1
-- ROOT_COINCIDENCE_R1_R2_COEFF  vars=['B_comp', 'mu_R', 'rho_br', 's_rho']  compound_radical_present=True  radical_check_seconds=0.00
   GATE_PATH locus_conditionset seconds=0.00  result_type=ConditionSet
   ATTEMPTING sp.solve on the same system
   SOLVE seconds=9.73  branches=2
   SOLVE_RESULT [ ...ELIDED, rule 5 — an S11 solution payload; see r3_XKIN_ANISO_D2 stdout... ]
```
The gate does NOT fire on `MAIN` D2/D3/D4 (`compound_radical_present=False` on the single
`ROOT_COINCIDENCE_R1_R2_COEFF` call per cell — §2 transcript below). Carried, leg-measured, not
re-measured here: at `XKIN_ANISO` D3 two refused pair-solves returned in 3.8 s / 3.9 s (STATUS top
block at `1d187a25`; verified by the orchestrator in the prior session).

## 2 · D2 — both truth directions of a premise fall through to UNDECIDED (BLOCKING)

Mechanism — substitution auto-evaluates a relational to a Boolean, which matches no branch of
`evaluate_premise` (`:692-711`):
```
$ python3 probe_r3_pure.py .../S11_stray_longitudinal_sympy_audit.py   # SECTION A
PREMISE=x > 0  BRANCH={x: 1}  SUBSTITUTED=True  SUBSTITUTED_TYPE=BooleanTrue  EVALUATE_PREMISE_RETURNS=None
PREMISE=x > 0  BRANCH={x: -1}  SUBSTITUTED=False  SUBSTITUTED_TYPE=BooleanFalse  EVALUATE_PREMISE_RETURNS=None
PREMISE=Eq(x, 1)  BRANCH={x: 1}  SUBSTITUTED=True  SUBSTITUTED_TYPE=BooleanTrue  EVALUATE_PREMISE_RETURNS=None
PREMISE=Eq(x, 1)  BRANCH={x: 2}  SUBSTITUTED=False  SUBSTITUTED_TYPE=BooleanFalse  EVALUATE_PREMISE_RETURNS=None
```
On the engine's own premises (SECTION B): every relational premise of `MAIN` D2/D3
(`rho_br > 0`, `mu_R > 0`, `B_comp > 0`, `k1**2 + ... > 0`) substitutes to `BooleanTrue` under the
engine's `candidate_witness` and returns `None`; only `AppliedPredicate` premises return `True`.

A branch carrying a provably violated premise still types UNDECIDED (SECTION C — `TEST_OBJECT`
contains a literal `False` for `k1**2 + k2**2 > 0` at `k1=k2=0`):
```
ADMISSIBILITY_ENTRY=(... (STATUS_TOKEN, UNDECIDED), (TEST_OBJECT, (True, True, True, False, ...)) ...)
```

Scale, re-measured (same probe_r3_cells run, `MAIN` cells; full transcripts `r3_MAIN_D*.out`):
```
== CELL MAIN D2 ...  CELL_SECONDS 2.6   REAL_ADMISSIBLE_STATUS_TOKEN_COUNTS {'UNDECIDED': 22}
== CELL MAIN D3 ...  CELL_SECONDS 7.1   REAL_ADMISSIBLE_STATUS_TOKEN_COUNTS {'UNDECIDED': 35}
== CELL MAIN D4 ...  CELL_SECONDS 30.1  REAL_ADMISSIBLE_STATUS_TOKEN_COUNTS {'UNDECIDED': 39}
```
Zero ADMISSIBLE, zero EXCLUDED anywhere, `XKIN_ANISO` D2 included ({'UNDECIDED': 19}, §1 transcript).

## 3 · Same class — `_INCONSISTENT`'s test object is built so it cannot evaluate

`:870` builds `sp.Eq(solution_payload, sp.Tuple(), evaluate=False)`; `symbolic_truth_status`
(`:640-647`) types by `== sp.true` / `== sp.false`. Probe (SECTION D), on a toy system, not S11:
```
SYSTEM=[y + 1, y + 2]  SOLVE_RETURNS=[]  PAYLOAD=()
INCONSISTENT_OBJECT=Eq((), ())  TYPE=Equality
SYMBOLIC_TRUTH_STATUS=((STATUS_TOKEN, UNDECIDED), (TEST_OBJECT, Eq((), ())), (OPERANDS, ((),)))
```
A solve that returned the empty set — the inconsistent case §5 names — types UNDECIDED. The only
route to PROVED_TRUE is `constant_false` (`:862-864`), which requires an equation with no free
symbols.

## 4 · Publish placement at HEAD

```
$ sed -n '1866,1870p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    if main_is_complete(completed_pairs):
        CURRENT_ISSUE_CONTEXT = "PUBLISH"
        try:
            ledger = merged_export(main_dim_data, run_pairs_payload, skipped_pairs_payload)
            write_exports(ledger)
```
This sits after the full `PACKAGE_ORDER` loop (`:1835-1861`); the gate reads only `MAIN`
(`main_is_complete`, `:1816-1818`), and `MAIN` is first in `PACKAGE_ORDER` (`:136-147`). So a
publish whose gate is satisfied ~5 min in waits on every other package. Cell timings measured above:
`MAIN` D2/D3/D4 = 2.6/7.1/30.1 s; `XKIN_ANISO` D2 = 141.2 s. Carried, leg-measured: `XKIN_ANISO`
exceeds 600 s at D3, scaling 4–5×/dim (STATUS at `1d187a25`). ⛔ Not re-measured — re-measuring it
means running the wall.

The run-record vocabulary conflates attempted-and-failed with never-reached:
`run_record_payloads` (`:1808-1813`) computes `skipped = declared − completed`, and the cell
exception path records `"cell skipped after ..."` (`:1859`). At end-of-run every declared pair was
attempted, so the conflation is latent; any mid-loop publish inherits it — the 18-false-rows
version round 2 removed is the measured instance.

Round-2 properties that must survive a placement change (both verified present at HEAD): a publish
failure is attributable to the publish (`PUBLISH` context, `:1867-1877`), and a failed run still
emits its §10 report (`:1901-1912`, `publish_error` re-raised only after the report).

---

# Fold-time verifications (after both brief legs reported)

Leg reports: `~/.s11_build/fix3_brief_leg_grok.log` (probes `/tmp/s11r3_leg_Wlem/`),
`~/.s11_build/fix3_brief_leg_codex.log` (probes `/tmp/s11r3_leg_xfHB/`). Every finding folded below was
re-verified by the orchestrator's own command first (rule 13).

## 5 · The decided-false predicate route (Grok finding 2, verified)

```
$ python3 - <<'EOF'   # imports engine via importlib, then:
premise = sp.Q.real(k1); branch = {k1: sp.I}; m.evaluate_premise(premise, branch)
EOF
PREMISE=Q.real(k1)  BRANCH={k1: I}  SUBSTITUTED=Q.real(I)  TYPE=AppliedPredicate
argument.is_real=False
EVALUATE_PREMISE_RETURNS=None
ADMISSIBILITY_ENTRY_TOKEN=(STATUS_TOKEN, UNDECIDED)
```
`:701-710`: only `is_* is True` returns; a decided-false attribute falls to `return None`.

## 6 · The spec's full `_INCONSISTENT` obligation (both legs, verified)

```
$ sed -n '250,251p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
⭐ `_IDENTICALLY_SATISFIED` and `_INCONSISTENT` are what separate the three degenerate cases, and they are
computed from `_EQUATIONS`, ⛔ never read off the solver's empty token.
```
The first brief draft cited `:250` but paraphrased only the first clause.

## 7 · The nonempty locus a solver-empty test would call inconsistent (Codex finding 1, re-run)

```
$ timeout 600 python3 /tmp/s11r3_leg_xfHB/probe_item2_physics_closure.py \
    research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
ROOT 1 SOLVE_RESULT []
ROOT 1 WITNESS_RESIDUALS [0, 0, 0, 0]
ROOT 1 ROOT_MATRIX_WITNESS_RANK 0
ROOT 1 WEAKEST_AGREEING_INCONSISTENT_STATUS ((STATUS_TOKEN, PROVED_TRUE), (TEST_OBJECT, True), ...)
ROOT 2 [same shape]
TOTAL_SECONDS 1.934449
```
Witness assignment and locus expressions ELIDED (rule 5 — computed S11 facts a builder could anchor
on); full stdout: `~/.s11_build/round3_probes/reverify_item2_physics_closure.out`. `sp.solve` returning
`[]` while an exact premise-satisfying assignment zeros every residual ⇒ emptiness of the solver return
is not a proof of inconsistency.

## 8 · The spec pins the post-sweep run record (Codex finding 2, verified)

```
$ sed -n '1038,1040p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
⛔⛔ **THE RUN RECORD MUST BE OBSERVED, ⛔ NOT DECLARED.** ⭐ Accumulate a `(package, D)` pair **only after
that package has finished emitting**, and emit `RUN_PAIRS` / `SKIPPED_PAIRS` **after** the sweep, with
`SKIPPED = declared ∖ completed`.
```
⇒ a mid-loop publish may not move those objects; its publish-time record must be distinguishable.

## 9 · Leg re-measurements folded as leg-produced (not orchestrator-re-run)

- `XKIN_ANISO` D3 at HEAD: probe exit 124 at exactly 600.00 s (Codex, `probe_cell_XKIN_ANISO_D3_timed.stdout`);
  three `ROOT_COINCIDENCE_*_COEFF` refusals emitted before the wall (`independent_tags_XKIN_ANISO_D3.stdout`).
- Un-refused `XKIN_ANISO` D2 cell completes in 153.5 s vs ~139 s gated (Codex ablation).
- Independent solve of the D2 refused system: 9.46 s / 11.72 s (Codex), ~9.5 s (Grok), 9.73 s (orchestrator, §1).

## 10 · Leg remedies NOT adopted verbatim, and why

- Grok: *"drop the escape hatch"* — kept, restructured: an honest builder exit must exist (a rule that
  forces fabrication is worse than none); the acceptance now makes the hatch a §10-recorded path instead
  of a silent one.
- Codex: *"acceptance must exercise an actual solver-empty but equation-nonempty S11 locus"* — the
  property is adopted; naming the locus in the brief is not (it hands the builder a computed S11 fact,
  rule 5). The script legs carry that check instead.
