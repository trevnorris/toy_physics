# P1-WL residual-mode emit — decision-leg review record (2 legs, computation-backed, convergent → v2 redesign)

Directive `directives/S11c_b_p1_wl_residual_emit_directive.md` (orchestrator-written ⇒ Codex + Grok decision legs,
rule 7). Prompt `directives/_legs/S11c_b_p1_wl_decision_review.md`. Raw logs
`~/.s11_build/S11c_b_p1_wl/decision_{codex,grok}.log`. **v1 REJECTED by BOTH legs, independently, computation-backed
(file:line + literal greps).** v2 folds both once (rule 7, no re-leg).

## The convergent architectural finding — v1's single-CASE framing was wrong
Both legs, independently:
- Committed `row_residual` UNIONS the case set and RAISES on any non-aligned key — `ADMISSIBILITY_OPERATOR_OPERAND`
  at `row_residual.py:595-603`, coupling at `:589-592`, slab likewise. ⇒ a 1-case WL `.out` cannot be consumed
  against a 4-case PY `.out`. (Orchestrator rule-13 verified: read `:595-603` and `:585-592` — both raise
  `InputError` on `aligned.py is None or aligned.wl is None`.)
- `ENERGY_BASIS_*` are BRANCH-scoped (pre-loop, `…audit.wl:1884/2131`, extractor axis BRANCH only
  `comparator.py:736`); `ADMISSIBILITY_SUPPORT_OPERAND` tables over ALL branches×densities (`…:2384`). A singleton
  main-loop restriction yields inconsistent scopes (1 operator case vs 4 support/energy cases).
⇒ **The correct shape is a `S11CB_PRIMARIES_ONLY` gate SYMMETRIC with the PY engine's existing one**
(`sympy_audit.py:41` `PRIMARY_TASKS`; `:5471` skip of `CONTROL_TASKS`): emit the PRIMARY families for the FULL case
set, skip the controls. This dissolves the singleton mismatch (WL case-set matches PY), the branch-scope
inconsistency, and the case-selector ambiguity. v2 is this redesign.

## Convergent defects folded into v2
1. **E1 omitted `ADMISSIBILITY_OPERATOR_OPERAND`** (both; `row_residual:595` hard-RAISES) + the other shared-primary
   families (Codex: 4 `ENERGY_BASIS_*`, 3 `ADMISSIBILITY_*`; E1 had 1 energy + 0 admissibility). v2 E1 lists all
   families the comparator (`comparator.py:84-95,1084-1096`) + `row_residual` (`:566-602`) join, over all cases.
2. **E2 "skip 2204-2258" was a naive line-range** (both): (a) the tower-spool built at `:2233` is CONSUMED at
   `:2300-2311` (`CONTROL_BACKGROUND_JET_TOWER_OPERANDS`, between `MU_THETA_OPERATOR` `:2297` and `COUPLING_KERNEL`
   `:2313`) — skipping only the build corrupts the stream (rule-13 verified: `:2297/2301/2305/2313`); (b) many LATER
   control loops re-run `evaluatedModel`/`extractCoupling` (MATERIAL `:2438`, corrupted `:2466`, frozen `:2522`,
   form-ablation `:2742/2772`, uniform-limit `:2885`) — so a mere region-skip does NOT save the memory. v2 E2 names
   the emit-scope as OBJECTS (all WL analogs of the PY `CONTROL_TASKS`) with an allowlisted exit.
3. **~8 GB not established for the full transcript path** (Codex #5): the primary emit also builds
   `mainKernelOrigins` (`:2199` → `kernelFromOrigin`→`extractCoupling` `:1806/1824`) + `ADMISSIBILITY_SUPPORT` table.
   STEP 0 measured only `evaluatedModel` + one `FINAL_KERNEL` ⇒ a LOWER bound. v2 E5 MEASURES the real footprint,
   with a documented single-case + P2-intersection fallback if it exceeds the box.
4. **E3 default-path invariance unfalsifiable** (both): the "demonstrate it's guarded" disjunct is a typed sentence;
   the full-`.out` byte disjunct is ≥64 GB work. v2 E3 = STRUCTURAL git-diff (changed lines inside a gate branch
   false when unset; unset iterators/helpers unchanged) + PAIRED REDUCED-SCALE byte-identity (HEAD vs patched,
   mechanically-identical scale edit, `cmp`, both exit 0).
5. **Faithfulness had no acceptance test** (both): a builder could `emitShared["SLAB_OPERATOR", <hand association>]`
   omitting `DIVERGENCE_FORM_SOURCE`/`FACE_SHAPE_SUBSTRATE`/`VIRTUAL_CONSTRAINT` that `extract_slab` reads → P3's
   μ_θ vanishes. v2 E7 = per-tag byte-equality vs HEAD's payload for the same case (proves the engine's own
   `modelRecord`/`kernelRecord`/origin maps produced it).
6. **E4 only checked `load_wl` tag-presence** (both): `load_wl` never parses the payload (`comparator.py:177-209`).
   v2 E4 requires every required family to pass its COMMITTED extractor with no exception + correct scope.
7. **E6 case-selector ambiguous** (both): env name not pinned, `route` not a production axis (loop hard-codes
   `EULERIAN` `:2192`; `MATERIAL` is a control route), `caseKey` is `branch|density` not slash-triple, empty≠unset.
   v2 E6 DISSOLVES the selector — a boolean `S11CB_PRIMARIES_ONLY` mirroring PY; production route stays EULERIAN.
8. **Term-origin SUM_RESIDUAL + operational policy** (both): `SLAB/COUPLING_KERNEL_TERM_ORIGINS_SUM_RESIDUAL`
   (`:2278/2321`) classified as required integrity records (v2 E1); E5 pins timeout/floor/exit-0/one-kernel-lock.

## Status
v1 rejected; v2 REDESIGNED (single-case → symmetric `S11CB_PRIMARIES_ONLY`) folding both legs once (rule 7). ⇒ Codex
build (WL) → 2 build legs (fresh Claude + Grok, Mathematica serialized). ⚠ The footprint of the all-case primaries
run is the one genuine open measurement (E5); STEP 0 established only the per-case lower bound.
