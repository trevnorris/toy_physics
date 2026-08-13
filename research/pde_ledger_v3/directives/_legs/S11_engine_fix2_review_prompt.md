# Independent review — the S11 SymPy engine after fix round 2

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`
(Codex-written, uncommitted). Built to
`research/pde_ledger_v3/directives/_legs/S11_engine_fix_round2_brief.md`, which had two legs.
Physics authority: `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`. `CLAUDE.md` binds.
Previous baseline: `19d4a44a`. Round-2 builder log: `~/.s11_build/codex_fix_round2.log`.

Round 1 cleared a `MemoryError` wall and its rank construction was verified by two legs as the exact
algebraic rank. Round 2 fixed a broken export path, removed a mid-loop publish, restored a regressed
nullspace basis, and **removed hard-coded size gates that were refusing computations before calling the
CAS**.

## ⛔⛔ A FORM ABLATION IS MANDATORY. ⛔ NOT OPTIONAL.

Change the **structure** of a load-bearing object — flip a sign *and* an off-diagonal, collapse two
independent symbols into one — re-run, report the **literal** diff. ⭐ A **coefficient** rescale tests
arithmetic; only a **FORM** change tests physics.
⚠ Measured 2026-08-04: a script emitted a physics conclusion as a typed sentence and the ablation produced
**byte-identical output**. ⛔ Eight fidelity legs missed it by reading.
⇒ ⭐ Ask of every claim: **which line computed this?** Give the line number or report it uncomputed.

## What to check

1. ⛔⛔ **THE SIZE GATES ARE GONE AND THE LOCI NOW COMPUTE. ⭐ ARE THE NEW OBJECTS RIGHT?** `ISSUES` at
   `MAIN` D5 went from four *"solve unavailable"* entries to **empty**. Loci that were tokens are now
   solution sets. ⭐ **Verify independently that what is emitted is the solution set of the stated system**,
   not a Gröbner artifact, a partial branch, or a set over the wrong variables. ⚠ The builder now computes a
   Gröbner basis once and solves from it — ⭐ check that this is the same object as solving the original
   system.
2. ⛔ **DID ANY VALUE MOVE AT D2/D3/D4?** Round 2 changed what is computed for the loci. ⚠ A changed value
   is a **finding to report**, ⛔ never a regression to tune away — ⛔ do not treat prior output as a
   reference. Say precisely what moved and whether it is now more or less correct.
3. ⛔⛔ **THE NULLSPACE RESTORATION.** The builder accepts an *undecided* source minor **only if the live
   matrix residual validates to zero**. ⭐ Is that sound? ⚠ Can it accept a vector that is not in the
   nullspace, or miss one that is? Check `XKIN_ANISO` D3 roots 2 and 3 and report the residual you compute
   yourself.
4. ⭐ **THE PUBLISH PATH.** A publish failure must be attributable to the publish, ⛔ never to a physics
   cell; §10 must still emit; stale exports must be deleted; and no row may claim a run record it does not
   have. ⭐ Verify on a reduced-`PACKAGE_DIMS` copy under `/tmp` — ⛔ never by editing the engine.
5. ⭐ **§10 ATTRIBUTION.** All `ISSUES` are now emitted with cell/publish attribution. ⭐ Is every entry
   attributable, and can the report still be truncated by anything?
6. ⭐ **ITEM 7's WALL.** `XKIN_ANISO` D3/D4 radical root-coincidence loci now emit `ConditionSet` plus a
   §10 entry. ⭐ Is that an honest unavailable construction, or a computation degraded to a token — the
   defect round 2 removed elsewhere? ⚠ `S11_SHARED_PHYSICS.md:245` and `:285` bound what may be substituted.
7. ⭐ **WOULD A FULL RUN NOW PUBLISH?** ⛔ Do not run the full loop. Reason from what you can measure per
   cell, and say what would still block it.
8. **Anything else**, including whether either round's fix is unsound.

## Method

- **Every claim carries the command that produced it and its literal output** (rule 2). A prose
  re-derivation is discarded. Save probes to named absolute paths and report them.
- Report a finding only if it catches a way the physics could be wrong, or a wrong result could survive.

## Operational constraints — these bind both legs identically

- ⛔ **Do not run the full package loop.** It has OOM-killed this machine. Use
  `/home/trevnorris/.s11_build/repro_d5.py` (one cell, ~2 min) or per-cell probes.
- ⛔ Copy the script to `/tmp/<yourname>/` and ablate the **copy**. ⛔ Never modify the working tree.
- ⛔ Wrap every run in `timeout 600`. A 600 s hit is a failed ablation — report it and move on.
- ⚠ A watchdog kills any S11 python over 6 GB RSS (`~/.s11_build/watchdog.log`).

## Deliverable

A verdict per numbered item, the literal output of your form ablation, every defect with file:line and the
command that found it.
