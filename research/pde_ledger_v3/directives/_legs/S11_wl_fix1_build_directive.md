# Build directive — S11 WL engine fix round 1

You are repairing `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
at commit `c891769d`, in the working tree at `/var/projects/toy_physics`.

## The governing decision list — it binds every choice you make

Read `/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round1_brief.md`
and its measurements twin in `../_measurements/` first. Its seven numbered obligations ARE the
acceptance surface; nothing below overrides them. The spec for every emitted object is
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (the locus
protocol is §5). Two known committed-record defects listed in the brief's out-of-scope section
must not be touched.

## The one structural rule for any expression you author

The ONLY place physical symbols may be combined by hand is in constructing the action and the
ansatz — and this repair authors NEITHER. Every expression the repair emits must be REACHED BY
COMPUTATION from objects the engine already builds. Scripts print computed objects; they never
state conclusions; residuals are printed, then guarded — never asserted before emission.

## Order of work

1. Design and implement the uniform bounded route at the three QE sites the brief names. The
   route's deciding attempts and budgets are your choice, subject to obligations 1–3; document
   each attempt's budget in the code where it is set.
2. Run the 19-cell regression (obligation 5): every completed cell via
   `wolframscript -file <abs .wl> <PACKAGE> <D>` — sequentially, ONE kernel at a time — and
   byte-compare all emitted `WL_S11_*` tags per cell against the committed `.out` at `c891769d`.
   Report the literal comparison output per cell.
3. Run XKIN_ANISO 3, then 4, then 2 (obligation 4), each in a fresh kernel, each wrapped in a
   memory guard that kills the cell's own wrapper and kernel PIDs (never a process group) when
   available system memory drops below 1024 MB, logging the guard event. Timestamp every
   emission so silent-gap and wall measurements exist. Report last-emitted-tag and the guard/exit
   record for any cell that does not complete.
4. Run the obligation-3 decision probe and the obligation-6 manifest census against every newly
   completing cell, and report their literal outputs.
5. Write the build report with, for every measured claim, the command and its literal output.

## Operational constraints

- Mathematica licence has TWO seats and one must stay free: never more than one kernel at any
  moment, including your own test evaluations.
- Absolute paths in every command. All scratch and logs go under `/home/trevnorris/.s11_build/fix1_build/`
  (create it), NOT under the repository and NOT under /tmp.
- Never run the full sweep in one kernel. Never modify anything outside the target `.wl` and your
  scratch directory. Do not commit; the orchestrator commits after independent review.
- If a run is long it must be visibly emitting; long + silent is the failure mode you are here to
  remove. Do not leave any kernel running when you finish.
