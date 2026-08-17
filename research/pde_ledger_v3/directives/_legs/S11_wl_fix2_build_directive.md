# Build directive — S11 WL engine fix round 2 (strata memory wall)

You are repairing `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
at commit `a8f26909`, in the working tree at `/var/projects/toy_physics`.

## The governing decision list — it binds every choice you make

Read `/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round2_brief.md`
and its measurements twin in `../_measurements/` first. Its six numbered obligations ARE the
acceptance surface; nothing below overrides them, and the round-1 brief sections it binds
verbatim (`S11_wl_engine_fix_round1_brief.md`: obligations 1, 5, 6, 7 and the Post-legs fold)
bind you too. The spec for every emitted object is
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` (the locus
protocol is §5). Three registered committed-record defects listed in the brief's out-of-scope
section must not be touched.

## The one structural rule for any expression you author

The ONLY place physical symbols may be combined by hand is in constructing the action and the
ansatz — and this repair authors NEITHER. Every expression the repair emits must be REACHED BY
COMPUTATION from objects the engine already builds. Scripts print computed objects; they never
state conclusions; residuals are printed, then guarded — never asserted before emission.

## Order of work

1. Design and implement the uniform bounded routes for every call class the brief's obligation 2
   names. The deciding attempts, fallback algebra, and budgets are your choice, subject to
   obligations 1–4; document each budget in the code where it is set. Remember obligation 3's
   exercise requirement: every fallback branch must be demonstrably reachable on at least one
   operand outside XKIN_ANISO D2/D4 — build the synthetic-operand exercise into your acceptance
   harness and report its literal output.
2. Build the obligation-4 probe harness (independent exact-rational decision probe with
   no-starvation budgets, the planted-record calibration, the residual recompute, and the SymPy
   completeness containment check). Run the calibration FIRST and report its literal output —
   a probe that cannot fail the planted record is a broken instrument and the build stops there.
3. Run the 19-cell regression (round-1 obligation 5): every completed cell via
   `wolframscript -file <abs .wl> <PACKAGE> <D>` — sequentially, ONE kernel at a time — and
   byte-compare all emitted `WL_S11_*` tags per cell against the committed `.out` at `a8f26909`.
   Report the literal comparison output per cell.
4. Run XKIN_ANISO 3, then 4, then 2 — each ONLY via
   `/home/trevnorris/.s11_build/fix1_build/run_guarded_cell.sh <D>` (sha256 pinned in the twin;
   do not modify it), fresh kernel per cell. Arm the RSS sidecar sampler (a separate script
   writing `<stem>.rss_samples.tsv` beside the runner's outputs; do not modify the runner
   itself) for D4 and D2. Obligation 5 is absolute: ANY guard death or incompleteness on D2 or
   D4 is a build failure — STOP at the first death, preserve the death record, report it, and
   do not iterate past it.
5. Run the obligation-4 probe census and the round-1 obligation-6 manifest census against both
   newly completing cells, and the obligation-6 (round-2) partial-record comparison — every
   committed partial record present, undecided-class tokens may only strengthen, any decided
   token that changes value is a HALT-AND-REPORT. Report literal outputs for all three.
6. Write the build report with, for every measured claim, the command and its literal output,
   including the per-stratum RSS high-water profile for D2 and D4.

## Operational constraints

- Mathematica licence has TWO seats and one must stay free: never more than one kernel at any
  moment, including your own test evaluations.
- Absolute paths in every command. All scratch and logs go under
  `/home/trevnorris/.s11_build/fix2_build/` (create it), NOT under the repository and NOT under
  /tmp — except the runner's own outputs, which land where the runner puts them.
- Never run the full sweep in one kernel. Never run XKIN_ANISO 2 or 4 unguarded. Never modify
  anything outside the target `.wl` and your scratch directory. Do not commit; the orchestrator
  commits after independent review.
- If a run is long it must be visibly emitting; long + silent is the failure mode you are here
  to remove. Do not leave any kernel running when you finish.
