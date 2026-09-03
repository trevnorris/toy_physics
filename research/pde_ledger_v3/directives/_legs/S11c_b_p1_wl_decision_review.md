# Decision-list review — S11c-b P1-WL residual-mode single-case emit (WL engine)

## Artifact
`research/pde_ledger_v3/directives/S11c_b_p1_wl_residual_emit_directive.md` — an orchestrator-written DECISION LIST
for a WL-engine change. You are ONE of two independent decision legs (rule 7). The builder will TRUST this list:
everything downstream is checked twice, the list itself is checked once, here. Find its defects now — a productive
review FINDS things; "looks fine" is weak evidence.

## What the change is for
Produce a fresh single-case WL `.out` — the CURRENT un-frozen slab operator (#89b) + the primary COUPLING
`FINAL_KERNEL` (#90-era), for ONE case — that the cross-engine `row_residual` can PARSE, on this 30 GB box. The
committed WL `.out` is stale (frozen operator), so it cannot be used. STEP 0 measured that the primary single-case
operator + kernel fit at ~8 GB, while the per-case loop's tower-depth control variants force ~16 GB; residual-mode
must emit the primaries for one case and skip those variants, WITHOUT changing the default full-run emit.

## Context you are handed (read the CODE, cite file:line)
- WL engine `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`: the per-case emit loop
  (~`2186`-`2320`); the tower-variant build (`operatorTruncated/Extended`, `kernelTruncated/Extended`, tower-spool,
  tractability — ~`2204`-`2258`); `compactModel = KeyDrop[processed,…]` + `AssociateTo[mainModels,…]` (~`2259`-`2261`);
  the post-loop primary emits `emitShared["SLAB_OPERATOR", Map[modelRecord, mainModels]]` (~`2268`), term-origins
  (~`2274`-`2318`), `MU_THETA_OPERATOR` (~`2287`-`2297`), `COUPLING_KERNEL` (~`2312`-`2313`); `modelRecord` (~`1570`);
  `extractCouplingData`/`FINAL_KERNEL` (~`1744`); `emitShared` (~`35`); `ENERGY_BASIS_NEW_INVARIANTS` emit.
- Consumer `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py` (`load_wl`, `extract_slab`,
  `extract_coupling`, `extract_mu`) and `research/pde_ledger_v3/scripts/S11c_b_row_residual.py` (its `_family_cases`
  family list) — the residual-mode `.out` must carry every WL family these JOIN on.
- STEP 0 record `research/pde_ledger_v3/directives/_measurements/S11c_b_step0_residual_scope.md`.

## What to check (find defects; ground every claim in code, cite file:line)
1. **Required-tag COMPLETENESS (item E1).** Independently GREP the comparator + `row_residual` for EVERY WL family
   they read in the SLAB / COUPLING / MU_THETA / §3a-basis join. Is E1's list complete, or would residual-mode risk
   OMITTING a family the comparator needs (⇒ the residual silently drops that object)? List the exact families with
   file:line. Flag any family that is source-scoped or run-scoped (not per-case) and could be lost when the loop is
   restricted to one case.
2. **DEFAULT-PATH INVARIANCE (item E3).** Is the invariance actually VERIFIABLE as written, or unfalsifiable? A gate
   that touches shared code before the env check could perturb the default run. State the concrete check the build
   legs must run to prove the default `.out` is byte-identical.
3. **FOOTPRINT (items E2/E5).** Does skipping ONLY the ~2204-2258 region actually hold the run at ~8 GB — or does a
   PRIMARY emit (the post-loop SLAB/COUPLING payload, term-origins, or `ENERGY_BASIS_NEW_INVARIANTS`) DEPEND on an
   object built inside the skipped region, so skipping it either breaks the emit or fails to save the memory? Inspect
   what `2204-2258` feeds downstream.
4. **FAITHFULNESS (governing invariant 2).** Could a builder satisfy E1 by REBUILDING a tag payload by hand (drifting
   from the engine's real emit) instead of calling `modelRecord`/`kernelRecord`/`emitShared`? Does the directive
   close that hole?
5. Any ambiguity a builder could satisfy WRONGLY: the case-selector format/validation (E6), source-scoped vs
   case-scoped tags, the term-origins sum residuals, the licence/timeout constraints. Any missing property.

## Required method
This is a DOCUMENT review — do NOT modify the tree. Ground every claim in the code (file:line); where a claim is
checkable by a command (grep the comparator for its WL families; inspect the 2204-2258 dependencies), RUN it and
quote the LITERAL output — a prose assertion with no command behind it is discarded. Report defects as a numbered
list; for each, name the directive item (E#) it fixes and the file:line evidence. If you believe the list is
complete and correct, say what you checked and why nothing survives — but a clean pass is weak evidence.
