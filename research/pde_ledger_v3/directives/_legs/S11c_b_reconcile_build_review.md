# Independent review — S11c-b hand-coded reconcile layer (a build leg; the instrument, not the physics)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_handcoded_comparison.py` (Codex-built).
It imports the committed comparator `S11c_b_cross_engine_comparator.py`, applies a HAND-VERIFIED enumerated
WL→PY rename map (80 entries), re-checks each core-family residual to zero, and prints per-family
`MATCH` / `FLAG (residual shown)` / `COVERAGE` / `NAMESPACE_INCOMPLETE`. It decides NOTHING about physics
(no PASS/FAIL/VERDICT/target). Default run reported 2 MATCH, 14 FLAG.

## The ONE failure mode that matters
A reconcile layer earns its keep only if it **cannot manufacture a false MATCH** — i.e. cannot fold two
genuinely different objects into agreement, and cannot silently reduce away a real difference. A false FLAG
(surfacing a difference that is actually representational) is safe — the orchestrator adjudicates those. A
**false MATCH hides a cross-engine disagreement** and is the defect this leg exists to catch. Report a
finding only if it catches one of:
1. **A rename that maps two DIFFERENT physical quantities to one name**, creating a false MATCH or spuriously
   cancelling a residual. (The governing sin — never blanket-collapse.)
2. **A blanket `camelCase→snake_case` transform** (rather than an enumerated, source-cited, per-entry map).
3. **A representative/scale collapse that WAS applied** where it must be FLAGGED: the ENERGY_BASIS_* non-unique
   quotient (a variable-coefficient IBP first-background-jet term is PHYSICS — must not be folded), the
   coupling-kernel adjointness residual (already IBP-reduced by the comparator — no further fold), `bRho`, or
   the quotient-specific `DivGrad` gamma representatives (the builder says it FLAGGED these — verify it did).
4. **A MATCH that is not real** — the two default MATCH families (`ADMISSIBILITY_SUPPORT_OPERAND` 20/20,
   `ENERGY_BASIS_COUNT` 2/2): confirm each is a genuine cross-engine agreement, not an artifact of an
   over-broad rename or a trivial/degenerate operand.
5. **The layer deciding physics** — any PASS/FAIL/VERDICT/target, or a residual asserted rather than printed.

## The completeness question (report your assessment even if no defect)
Is the 80-entry rename map COMPLETE for the in-scope core families? After the renames, do the FLAG residuals
still contain **un-bridged cross-engine spelling pairs** — a WL camelCase atom with an obvious snake_case PY
sibling that is simply missing from the map — which would inflate the FLAGs with naming noise? Or do the FLAG
residuals contain only genuine same-namespace differences (different coefficients, scales, representatives,
extra/missing terms)? Sample several FLAG residuals (across COUPLING_KERNEL, SLAB_OPERATOR, MU_THETA_OPERATOR,
ADMISSIBILITY_*) and state which. This decides whether the FLAGs are ready for physics adjudication or the map
needs another pass. ⛔ Do NOT propose the physics verdict on any FLAG (that is the orchestrator's call) — only
whether the residual is naming-noise vs a genuine (representational-or-real) difference.

## What you are handed
- The reconcile layer + the comparator it imports.
- The two committed transcripts (`scripts/out/S11c_b_brane_operator_sympy_audit.out`,
  `mathematica/out/S11c_b_brane_operator_mathematica_audit.out`) and the two engine SOURCES
  (`scripts/S11c_b_brane_operator_sympy_audit.py`, `mathematica/S11c_b_brane_operator_mathematica_audit.wl`)
  and the spec `directives/S11c_b_SHARED_PHYSICS.md` — use the sources to verify a rename is the SAME variable.
- ⛔ You are NOT handed a "correct" set of MATCH/FLAG verdicts. A genuine FLAG is expected and is not a bug.

## Required method — ABLATE, do not merely read
1. **`--drop-rename <WLname>` ablation (the load-bearing test):** for a SAMPLE of the 80 renames spanning
   families, drop each and re-run; a GENUINE rename (same variable) must turn a MATCH into a FLAG (or move a
   residual). A rename whose removal changes NOTHING is dead or suspicious — report it. Report literal diffs.
2. **False-MATCH probe:** for the 2 MATCH families, read the operands from the transcripts and independently
   confirm they agree for the SAME physical reason (not because a rename mapped two different atoms together).
   Try to BREAK a MATCH by removing the rename(s) its cases depend on; if nothing breaks it, dig into why it
   matches.
3. **Non-fold probe:** grep the layer for any energy-basis / adjointness / `bRho` / `DivGrad`-gamma
   substitution; confirm those are NOT in the applied rename map (they must remain FLAGGED). If one IS folded,
   that is a finding.
4. **Source-verify a sample of renames:** pick ~8 map entries (e.g. `WZero→W_0`, `muR→mu_R`, a gamma mapping)
   and confirm from BOTH engine sources that the two names denote the same physical quantity (cite lines). A
   rename mapping near-but-different quantities is the critical finding.
5. For any claim, name the line. Report any `assert` on a measured payload.

## Physics filter
Report a finding only if it catches a way the reconcile could manufacture a false MATCH, hide a real
difference, blanket-collapse, or decide physics. Do NOT report "engine A and engine B differ here" (a real
FLAG is expected) or "would be wrong on different input".

## Ablation sandbox & operational constraints (identical for both legs)
- Copy the layer to `/tmp` and ablate the COPY (or use its `--drop-rename` flag against the committed file
  read-only). ⛔ NEVER modify the working tree.
- Pure Python — no Mathematica kernel, no licence seat. Run everything FOREGROUND; ⛔ no background job, ⛔ no
  monitor/poll loop; report directly. (A background-monitor leg has stalled twice on this task.)
- A full default run is ~13 min. Run it at most once or twice; do not loop it. `--drop-rename` re-runs are the
  same cost — sample, don't exhaust.
- Save every ablation script/command AND its literal stdout to named /tmp paths; report those paths. A prose
  "looks fine" with no command + stdout is discarded.

## Report
Per probe: the ablation, its literal diff, verdict (sound / manufactures-false-MATCH / folds-what-must-be-
flagged / decides-physics). Then your completeness assessment (map complete vs naming-noise-remains, with the
sampled FLAG residuals you based it on). Then any unaddressed definition-of-done item. Findings, not narration.
