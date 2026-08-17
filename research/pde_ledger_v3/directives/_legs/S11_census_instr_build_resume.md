# Resume — S11 locus-census instrument build

You are resuming the build governed by
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_census_instr_build_directive.md`
at commit `07ab921d` (read it and its brief first; they bind everything). A prior run completed
order-of-work step 1 — the five instrument modules exist under
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/` and both records' population
reconciliation printed `POPULATION_RECONCILED` with `unreconciled_count=0` — and then stopped
honestly because the scratch directory was read-only in its sandbox. That is fixed; verify with
`touch /home/trevnorris/.s11_build/census_build/.probe && rm` before anything else.

Continue from order-of-work step 2 exactly as the directive states: production-path planted
calibration FIRST (report its literal output; stop if it cannot fail), then the full census over
both committed records with complete stdout preserved under
`/home/trevnorris/.s11_build/census_build/`, then the named acceptance reducer, then the build
report. Re-run the step-1 reconciliation as part of the deliverable so the record is one
contiguous run. Every constraint in the directive still binds: census findings are reported,
never repaired; no engine, spec, committed record, or register modification; no commits.
