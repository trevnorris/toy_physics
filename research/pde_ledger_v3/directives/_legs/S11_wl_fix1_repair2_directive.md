# Build directive — S11 WL fix round 1, post-legs repair (small, scoped)

Target: the WORKING-TREE
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl`
(carries the uncommitted round-1 repair; baseline commit `19b607ab`). Working directory
`/var/projects/toy_physics`.

Governing decision list:
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_wl_engine_fix_round1_brief.md`
— read it first, especially the "Post-legs fold" section. Implement EXACTLY the first two fold
rulings; the third is a scope ruling requiring no code change. Nothing else in the engine may
change.

## The two changes

1. Wherever `boundedQERoute`'s outcome fills an operand slot that names the primary's return
   (e.g. `UNRESTRICTED_REDUCTION`), bind the PRIMARY attempt's actual outcome to that slot —
   including its budget-expiry `Failure` object when it expired — and append the route's full
   attempt sequence to the record's operands whenever more than one attempt ran, regardless of
   whether the outcome decided. The status token and test object continue to come from the
   route's outcome as they do now.
2. `emitComponentCount` (the count site): when the count's decision route ran more than one
   attempt, the emitted record carries that attempt sequence. Today the site builds
   `QE_DECISION_ATTEMPTS` and never emits it.

## Acceptance you must run and report (commands + literal output)

1. Byte-identity spot-run: `MAIN 2` and `XKIN_ANISO 3` fresh runs, all `WL_S11_*` tags
   byte-compared against the committed `.out` at the baseline: identical (their primaries never
   expire, so these changes must be invisible there).
2. A tiny-primary-budget probe run of `MAIN 2` (budget forced small in a /tmp COPY, never the
   working tree): show the affected records now carrying the primary's expiry object and the
   attempt sequence, with unchanged status tokens versus the same probe before your change.
3. Parse check (`wolframscript -file` syntax load) passes.

## Constraints

- One kernel at a time; every kernel run wrapped in `timeout 600` except the two spot-runs
  (XKIN_ANISO 3 needs ~300 s; give it `timeout 900`).
- ⛔ Do not run `XKIN_ANISO 4` or `XKIN_ANISO 2`.
- Scratch and logs under `/home/trevnorris/.s11_build/fix1_build/repair2/` (create it). Absolute
  paths. Do not commit. Leave no kernel running.
