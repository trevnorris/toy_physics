# S11 WL engine fix round 2 — decision list (strata memory wall)

Target: `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` at HEAD
of this commit. Companion measurements (generated): `../_measurements/S11_wl_engine_fix_round2_brief.md`.
The round-1 brief `S11_wl_engine_fix_round1_brief.md` is a governing document for this round: its
obligations 1, 5, 6, 7 and its **Post-legs fold** section bind this round verbatim — re-read them
there; they are not restated here.

## The measured wall (every claim here is twin-backed)

- **D4**: three runs all died at exactly emission 999 with the same last tag (two memory-guard
  kills plus one external harness kill at the same frontier), tag
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
  **The D2 wall is accumulation plus a terminal call**: from the first STRATUM5 emission to the
  guard kill is ~1,223 s of the 1,315 s run, spent in unbounded `FullSimplify`-chain calls of
  100–330 s each on nested-radical/Abs/I entries.
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

The DEFECT_REGISTER.md entries dated 2026-08-16 (pointwise `_IDENTICALLY_SATISFIED`; joint
`ROOT_COINCIDENCE`; the additive count-payload provenance extension): untouched, exactly as in
round 1. ⚠ The register's engine line numbers predate the round-1 repair and have shifted — the
register's TEXT, never a line number, identifies each defect. The SymPy engine, the spec, and the
comparator: untouched. Provenance fields remain ADDITIVE-ONLY: appended after the spec's pinned
payload fields, never replacing or reordering them.

## What must be true after the fix

1. **Round-1 obligations 1, 5, 6, 7 and the Post-legs fold rulings hold verbatim** (objects and
   spec strength unchanged; the 19-cell byte regression run for real against the baseline of this
   commit; manifest census for newly completing cells; generated twin; primary-named slots carry
   the primary's actual outcome and multi-attempt routes emit their attempt sequences).
2. **Every CAS call class on the live cell path is bounded.** No invocation of the locus `Solve`
   (wl:312), the spectrum `Solve`s (wl:768, wl:982), the minor construction
   `allMaximalMinors`/`Factor[Det[…]]`, or the `assumedRank` / `assumedNullSpace` /
   `engineSimplify` / `unrestrictedSimplify` family may run without an explicit resource bound —
   time AND memory (the D2 killer was 22.3 s old when the machine hit the floor; a clock alone
   cannot catch it). A bound lives in the shared helper and applies at every call site of that
   helper — ⛔ never selected by package, dimension, root index, stratum index, or tag identity.
   Any further call class the builder bounds is bounded the same uniform way. The diff is the
   evidence; a bound reachable only for some cells' identity is a build failure.
3. **Expiry is per-call and recoverable, and the route is uniform all the way down.** An attempt
   exceeding its budget yields a measured failure object; the route then runs its bounded
   fallback attempts and, only after all have run and failed, an undecided token carrying the
   attempt sequence per the round-1 fold. An expiry must never terminate the cell and must never
   reach the machine-level guard. Round-1 obligation 2's uniformity binds every multi-attempt
   route this round builds, generalized: selection of attempts, fallbacks, and post-processing
   depends ONLY on measured expiry/failure of prior attempts — never on package, dimension, root
   index, stratum index, tag identity, or the names of symbols. A structural dispatch (on radical
   content, degree, operand class) must be a total rule applied identically at every call site,
   and every fallback branch must be exercised by the acceptance harness on at least one operand
   OUTSIDE the two previously-dying cells (synthetic if necessary) — a branch exercised in
   practice by exactly one cell is an identity gate wearing a structural predicate, and the build
   fails.
   **Returns are never coerced.** A primary or fallback return that is not the solver's genuine
   branch list (aborted, expired, a failure object) is emitted AS that object — spec:245 ("the
   solution set exactly as your CAS returns it") and the opaque-object rule at spec:271–272 bind
   every new route; coercing a non-answer to an empty or truncated list is a build failure. A
   successful fallback's record carries the fallback's ENTIRE return — filtering or truncating
   branches between the solver and the emit is a build failure.
4. **Honest incompleteness, only after a real attempt** — round-1 obligation 3 with its
   acceptance probe, widened four ways:
   (a) every undecided-class rank / nullity / count record carries its live operands (the matrix
   and the premises it was computed under), appended additively after the spec's pinned payload
   fields, so the probe has something to attack;
   (b) the probe re-attacks every undecided-class record of every kind in each newly completing
   cell on that record's own emitted operands, with per-record budgets at least as large as the
   engine's own primary budget for that class — a probe budgeted to lose is not a probe, and the
   probe harness must demonstrate it can fail by deciding a planted decidable-but-undecided
   record before the real census runs;
   (c) any such record the probe decides is a build failure;
   (d) emitted residual-class payloads in newly completing cells are recomputed from their own
   record operands by the probe; a mismatch is a build failure.
   **Completeness has an oracle**: for every newly emitted locus `_SOLUTION` in the two
   previously-dying cells, an independent engine (SymPy) solves the same emitted `_EQUATIONS`
   under the same premises within a fixed budget; where the independent route completes, every
   branch/point it finds must lie in the emitted `_SOLUTION` up to algebraic equivalence — an
   omitted branch is a build failure. Membership specialization alone (round-1 obligation 6)
   does not test completeness; this does.
5. **Accumulation is part of the wall, and completion is the acceptance.** `XKIN_ANISO` D2 and D4
   run to cell completion — rc=0, `WL_S11_LOCAL_TAG_NAMES` emitted, guard=NONE — under the
   committed guard contract: `~/.s11_build/fix1_build/run_guarded_cell.sh` (sha256 pinned in the
   twin; the runner is not modified), fresh kernel per cell, 1 GiB available-memory floor,
   14,400 s outer wall, no silent inter-emission gap over 1,800 s. ⛔ ANY guard death or
   incompleteness on either cell is a build failure — there is no measurement exception this
   round; the builder STOPS at the first death, preserves the death record, and reports — never
   iterates past it. Kernel memory is instrumented as evidence: a sidecar samples the kernel's
   RSS on a fixed interval to a separate file for both cells, and the build report includes the
   per-stratum high-water profile — accumulation is measured, not inferred from survival.
6. **Committed partial records stay at least as strong.** Every record present in the committed
   D2/D4 partials is present in the new run. An undecided-class token may become decided.
   ⛔ A decided token that changes value, or becomes undecided, is a HALT-AND-REPORT — it is a
   finding about the routes and is never silently accepted or repaired by the builder.

## Builder operational constraints

As round 1 (one kernel at a time, two-seat licence, per-cell only, absolute paths, transcript
outside the repository, observable progress — long + silent is the failure mode under repair),
plus: ⛔ `XKIN_ANISO 2` and `XKIN_ANISO 4` are run ONLY via
`~/.s11_build/fix1_build/run_guarded_cell.sh`; never unguarded, never in parallel with anything.
