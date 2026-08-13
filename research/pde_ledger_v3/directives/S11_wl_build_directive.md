# S11 Wolfram engine — build directive

## Authority and boundary

Rewrite `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` in full. Its
product is the flushed stdout tag stream; that is its only write.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` is the sole physics
authority and the sole physics input: every equation, ansatz, premise, package and tag obligation enters
from it alone, and it wins every conflict with this directive. Implement its §§1–10 directly. Point to
rather than duplicate §4; §5, including its corollaries, no-verdict rule and locus protocol; §6 Q1–Q11;
§7; §8; §9; and §10. Add no expected value or acceptance criterion (`CLAUDE.md` rule 5).

The build's inputs besides the spec carry no S11 physics: this directive, `CLAUDE.md`, and the two
precedent citations in the next section. The engine this build writes reads no file at run time — it
re-derives from the spec's equations alone. `S9_export_chain_rebuild_directive.md:16-18`: this
independence is the design's one blindness control, and nothing else may be built pretending to be one.
The spec's §8 engine-local scoping (`:1097-1105`) applies as written; this engine neither imports a
`LEDGER` nor writes one.

Do not add parallel machinery for the spec's contracts: no per-cell completion registry and no
directive-local exit policy. §5 corollary 4 applies as a property, without named exceptions.

## Mechanical precedent

`research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl` is mechanical
precedent, not authority. Its emit site (`:25-30` at HEAD) shows the working shape: an internal tag
mapped to a standard emitted-name string at the emit site, a duplicate emitted name a hard stop, every
payload printed as `NAME: value` in `InputForm` as computed. Its physics is another step's. §8
standardises emitted-name strings only; symbol spellings inside the engine are its own.

That precedent does not demonstrate the spec's flush obligation (`:1044-1047`): `Print` alone does not
establish flushing on a redirected stream. Satisfying that obligation under redirected stdout is this
build's to establish; acceptance item 2 observes it.

## Run discipline — binds the build's own demonstration runs

- The engine runs as `wolframscript -file <absolute path> [PACKAGE D]`: two optional positional
  arguments, a §7 package name and an integer `D`, select a single cell; with no arguments it runs §7's
  declared sweep. No other argument shape.
- The build never runs the declared sweep to completion. Demonstrations are single cells. The full
  sweep, and with it the replacement of
  `mathematica/out/S11_stray_longitudinal_mathematica_audit.out`, belongs to the orchestrator after
  review. Demonstration stdout goes under the build's own scratch paths, never under `mathematica/out/`.
- The licence has two seats: run at most one kernel at any time.
- Kill criteria for a demonstration kernel, both by PID, both recorded with the cell they interrupted
  and reported, never answered by narrowing or cheapening the requested object: 600 seconds with **no
  new output** (the spec's runtime rule, `:1042-1054`, is about observable progress, not elapsed time —
  a visibly streaming cell may run long), or resident memory past 6 GB.

## Acceptance — executable, no expected values

1. Copy the finished `.wl` alone into an empty scratch directory. Run the same single-cell
   demonstration there and in-repo: both runs exit 0, the streams are non-empty and byte-identical, and
   after the run the scratch directory contains exactly the one file that was copied in. Any run-time
   dependence on a repository file, and any write other than stdout, breaks this.
2. One completed single-cell demonstration per §7 package, at a `D` of the builder's choosing from that
   package's declared sweep, each with stdout redirected to a file. On at least one demonstrated cell,
   the capture file must be observed growing while the kernel is alive. If every demonstrated cell
   completes too quickly to observe that, record the fact in the report — do not claim the flush
   obligation was demonstrated.
3. Sweep-path probe: start a no-argument run, kill it by PID only after tags from at least two distinct
   declared `(package, D)` cells have appeared, and compare: the first cell's complete tag block must be
   byte-identical to that same cell's single-selection run. A no-argument branch that emits one cell and
   exits, or emits a prefix and stalls, fails this.
4. After all demonstrations: `git status --porcelain research/pde_ledger_v3/mathematica/out/` prints
   nothing.

## Conflicts

Use the spec's §10 to report any conflict, ambiguity or unavailable construction. Do not fill such a gap
with new physics, an expected result, or a self-review mechanism.
