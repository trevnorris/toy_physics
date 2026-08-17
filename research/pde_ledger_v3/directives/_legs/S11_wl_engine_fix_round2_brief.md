# S11 WL engine fix round 2 — decision list (strata memory wall)

Target: `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` at HEAD
of this commit. Companion measurements (generated): `../_measurements/S11_wl_engine_fix_round2_brief.md`.
The round-1 brief `S11_wl_engine_fix_round1_brief.md` is a governing document for this round: its
obligations 1, 5, 6, 7 and its **Post-legs fold** section bind this round verbatim — re-read them
there; they are not restated here.

## The measured wall (every claim here is twin-backed)

- **D4**: guard-killed 3× fully deterministically at exactly emission 999, last tag
  `..._STRATUM3_ROOT2_N2_RANK_CHANGE_LOCUS_EQUATIONS` — the `emitCell` at wl:311. The kernel died
  inside wl:312, `solution = Quiet[Solve[And @@ equations, variables]]` (the following `_SOLUTION`
  never printed), on 16 radical-bearing minor equations in 7 unknowns, after ~193 s of terminal
  silence. The call carries no resource bound of any kind.
- **D2**: guard-killed 2× fully deterministically at exactly emission 1,411; the last tag closes a
  completed record block, so the killer left no trace of its own. The next expressions in code
  order are wl:1076–1077: `stacked = Join[matrixAtRoot, {wavevector}]; stackedRank =
  assumedRank[stacked, assumptions]` — `MatrixRank` with a `FullSimplify` zero test (wl:72–78),
  unbounded. The identical call on the sign-twin stratum completed in 327.6 s earlier in the same
  run; the dying instance was killed only 22.3 s in, with available memory already at the floor.
  **The D2 wall is accumulation plus a terminal call**: the two preceding strata account for
  ~1,223 s of the run's last ~1,240 s, in unbounded `FullSimplify`-chain calls of 100–330 s each
  on nested-radical/Abs/I entries.
- **D3 is the control**: same engine, same package, every stratum, 255 s / 2,501 emissions; its
  corresponding operands are radical-free and small.
- **Feasibility** (diagnosis: two independent analysts, all scripts + literal stdout archived):
  independent exact routes on the independently reconstructed death-site operands decide the D4
  object in well under 1 s and certify the D2 object in under 2 s where the completed twin needed
  327.6 s. Both walls are route/presentation costs, not object costs.
- **No cross-engine oracle exists for strata objects**: the SymPy sibling emits an empty
  `STRATUM_ORDERING` on every cell and never attempts either dying computation. Acceptance oracles
  for this round are therefore intra-record (exact specialization against emitted operands), as in
  round-1 obligation 6.

## Out of scope

The two DEFECT_REGISTER.md entries dated 2026-08-16 (pointwise `_IDENTICALLY_SATISFIED` wl:203;
joint `ROOT_COINCIDENCE` wl:887–890): untouched, exactly as in round 1. The SymPy engine, the
spec, and the comparator: untouched.

## What must be true after the fix

1. **Round-1 obligations 1, 5, 6, 7 and the Post-legs fold rulings hold verbatim** (objects and
   spec strength unchanged; the 19-cell byte regression run for real against the baseline of this
   commit; manifest census for newly completing cells; generated twin; primary-named slots carry
   the primary's actual outcome and multi-attempt routes emit their attempt sequences).
2. **The two measured killer classes are bounded.** No invocation of the wl:312 locus `Solve` and
   no invocation of the `assumedRank` / `engineSimplify` family may run without an explicit
   resource bound — time AND memory (the D2 killer was 22.3 s old when the machine hit the floor;
   a clock alone cannot catch it). A bound lives in the shared helper and applies at every call
   site of that helper — ⛔ never selected by package, dimension, root index, stratum index, or
   tag identity. If the builder bounds any further call class, it is bounded the same uniform way.
   The diff is the evidence; a bound reachable only for some cells' identity is a build failure.
3. **Expiry is per-call and recoverable.** An attempt exceeding its budget yields a measured
   failure object; the route then runs its bounded fallback attempts and, only after all have run
   and failed, an undecided token carrying the attempt sequence per the round-1 fold. An expiry
   must never terminate the cell and must never reach the machine-level guard.
4. **Honest incompleteness, only after a real attempt** — round-1 obligation 3 with its acceptance
   probe, widened: the bounded, independent, exact-rational decision probe re-attacks every
   undecided-class record of every kind (locus, rank, nullity, count, certificate) in each newly
   completing cell, on that record's own emitted operands; any such record the probe decides is a
   build failure.
5. **Accumulation is part of the wall.** `XKIN_ANISO` D2 and D4 run to cell completion under the
   committed guard contract — `~/.s11_build/fix1_build/run_guarded_cell.sh`, fresh kernel per
   cell, 1 GiB available-memory floor, 14,400 s outer wall, no silent inter-emission gap over
   1,800 s. A guard death is a measurement, not a papering-over target — ⛔ and its death record
   must be analyzed: a death inside any call class this fix claims to have bounded means the fix
   did not apply, and the build fails.
6. **Committed partial records stay at least as strong.** Every record present in the committed
   D2/D4 partials is present in the new run. An undecided-class token may become decided.
   ⛔ A decided token that changes value, or becomes undecided, is a HALT-AND-REPORT — it is a
   finding about the routes and is never silently accepted or repaired by the builder.

## Builder operational constraints

As round 1 (one kernel at a time, two-seat licence, per-cell only, absolute paths, transcript
outside the repository, observable progress — long + silent is the failure mode under repair),
plus: ⛔ `XKIN_ANISO 2` and `XKIN_ANISO 4` are run ONLY via
`~/.s11_build/fix1_build/run_guarded_cell.sh`; never unguarded, never in parallel with anything.
