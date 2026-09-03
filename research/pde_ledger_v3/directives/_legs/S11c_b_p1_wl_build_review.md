# Build review — S11c-b P1-WL `S11CB_PRIMARIES_ONLY` WL gate

## Artifact
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` — WORKING TREE (uncommitted). The
change is `git diff HEAD -- research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (+17
lines, insertion-only, a `S11CB_PRIMARIES_ONLY` emit gate). It was written by Codex; you are one of two independent
build legs. Derive and ablate from first principles — do NOT trust Codex's self-report.

## What the change must do (verify each; the directive is `directives/S11c_b_p1_wl_residual_emit_directive.md`)
When `S11CB_PRIMARIES_ONLY` is SET, the engine must emit ONLY the primary families — for the FULL case set (all
`branches × densities`, route EULERIAN) — and SKIP every control / tower-depth variant / ablation, then exit. When
UNSET, the emit must be byte-identical to HEAD. Required primary families (verify NONE are missing — a missing one
makes the downstream `row_residual` RAISE): `ENERGY_BASIS_VARIABLE/COUNT/NEW_INVARIANTS/OMISSIONS` (branch-scoped),
`SLAB_OPERATOR`, `SLAB_OPERATOR_TERM_ORIGINS`, `MU_THETA_OPERATOR`, `COUPLING_KERNEL`, `COUPLING_KERNEL_TERM_ORIGINS`,
`ADMISSIBILITY_OPERATOR_OPERAND`, `ADMISSIBILITY_SUPPORT_OPERAND`, `ADMISSIBILITY_RESIDUAL`, and the term-origin
SUM_RESIDUAL integrity records.

## What to check (ablate a /tmp COPY; ground every claim in literal stdout)
1. **Set-mode content (E1/E2):** run the engine with `S11CB_PRIMARIES_ONLY=1` at a TRACTABLE reduced scale (reduce
   the basis / cases as needed to fit time — state your reduction). Confirm the emitted `.out` carries EXACTLY the
   primary families above for all cases and ZERO control/tower tags (no `CONTROL_*`, no `DEPTH_TRUNCATED/EXTENDED`,
   no tractability/MATERIAL/frozen/uniform-limit tags). List the tags emitted.
2. **Extractor pass (E4):** the set-mode `.out` must pass the COMMITTED comparator extractors
   (`scripts/S11c_b_cross_engine_comparator.py`: `extract_slab`, `extract_coupling`, `extract_mu`, `extract_energy`,
   the admissibility/term-origin extractors) with NO exception and the correct `(BRANCH,DENSITY)`/branch scope.
   Show the literal parse result.
3. **Faithfulness (E7) — MANDATORY:** confirm each primary tag's payload in set-mode is BYTE-IDENTICAL to the same
   tag from an UNSET (HEAD-path) run at the SAME reduced scale — proving the payload came through the engine's own
   `emitShared`/`modelRecord`/`kernelRecord`/origin maps, NOT a hand-assembled association. Show the `cmp`/diff.
4. **Default-path invariance (E3):** with the env UNSET, the engine's emit is byte-identical to HEAD (at reduced
   scale). Also inspect the diff STRUCTURALLY: every changed line is inside a branch guarded by the env; no shared
   iterator/helper/emit-ordering on the unset path is altered. ⚠ Codex reports the dynamic
   `S11CB_MAX_MEMORY_USED_BYTES` trailer differs run-to-run (a self-observing counter) — CONFIRM that is the ONLY
   whole-file difference and that it is a benign runtime counter, not a physics change.
5. **FORM ablation — MANDATORY (this is the only test that catches the worst defect):** structurally CORRUPT the
   gate on your copy (e.g. make set-mode ALSO emit a control tag, or SKIP one required primary family, or exit
   BEFORE `COUPLING_KERNEL`) and show the check that BITES (the tag list / extractor / E7 byte-compare must change).
   A gate that is byte-identical under such a corruption is unimplemented.
6. **Footprint sample (E5):** report the set-mode reduced run's peak kernel RSS + min MemAvailable. (The FULL 4-case
   production footprint is a SEPARATE orchestrator step — you are NOT asked to run it.)

## Method + constraints (both legs identical)
⛔ Copy the artifact to /tmp and ablate the COPY. ⛔ NEVER modify the working tree.
⛔ Wrap EVERY `wolframscript` kernel run in `timeout 600`. A 600s hit is a FAILED ablation — report it and move on.
⛔ NEVER raise the timeout; ⛔ run ONE kernel at a time (2-seat licence). Kill an orphaned kernel by exact pid; if a
run dies with a healthy log, check `free -h` FIRST.
⭐ Save every ablation script AND its literal stdout to named absolute paths and report them.
Physics filter: report a finding only if it catches a way the gate could produce a WRONG or incomplete `.out` (a
missing family, an altered payload, a broken default path, a control leak); do not report "would differ on a
different input." A clean pass is weak evidence — say what you ablated and what bit.
