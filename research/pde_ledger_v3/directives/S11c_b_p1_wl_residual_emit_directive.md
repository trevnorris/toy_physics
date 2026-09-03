# P1-WL — WL `S11CB_PRIMARIES_ONLY` emit gate (fresh cross-engine `.out`) — decision list v2

⚠ **v2 — REDESIGNED after 2 decision legs (Codex + Grok, both computation-backed, convergent).** v1 proposed a
single-CASE selector; both legs showed that is the wrong architecture (branch-scoped `ENERGY_BASIS_*` and
all-branch `ADMISSIBILITY_SUPPORT` make a singleton inconsistent, and committed `row_residual` unions the case set
and RAISES on any non-aligned key, so a 1-case WL `.out` cannot be consumed against a 4-case PY `.out`). The correct
shape is a **`S11CB_PRIMARIES_ONLY` gate SYMMETRIC with the PY engine's existing one** — emit the PRIMARY families
for the FULL case set, skip all build-heavy controls. Folded once (rule 7), no re-leg. Decision-review record:
`directives/_measurements/S11c_b_p1_wl_residual_emit_directive.md`.

## Objective
Add an ADDITIVE, env-gated `S11CB_PRIMARIES_ONLY` mode to the WL engine
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`, mirroring the PY engine's existing
`S11CB_PRIMARIES_ONLY` (`scripts/S11c_b_brane_operator_sympy_audit.py:41,5471`): when set, emit ONLY the PRIMARY
families — for the full case set (all `branches × densities`) — and SKIP every build-heavy control / tower-depth
variant / ablation / uniform-limit path (the objects the cross-engine residual does not compare and that force the
~16 GB/case peak). When UNSET, the emit is byte-identical to HEAD. This produces a fresh 4-case primaries WL `.out`
(the CURRENT un-frozen operator #89b + #90-era kernel) whose case-set MATCHES the PY `S11CB_PRIMARIES_ONLY` `.out`,
so the committed comparator + `row_residual` consume both directly. This is the WL half of integration-pass P1 (PY
half = run the existing PY `S11CB_PRIMARIES_ONLY`, no build).

## Why (context — do not re-litigate)
The committed WL `.out` is stale (`d4adbd99`, pre-#89b → FROZEN operator); the residual must run on the CURRENT
engine. STEP 0 (`directives/_measurements/S11c_b_step0_residual_scope.md`) measured a single case's primary operator
+ `FINAL_KERNEL` at 7.95 GB RSS / 0.99 GB in-kernel — establishing the PER-CASE primary footprint; the per-case
streaming (`KeyDrop`/spool) frees between cases, so a primaries-only all-case run is EXPECTED to hold ~8 GB peak.
⚠ That is a LOWER BOUND, NOT established for the full transcript path — the primary emit also builds
`mainKernelOrigins` (`…:2199`, which re-invokes `extractCoupling` via `kernelFromOrigin` `…:1806/1824`) and the
`ADMISSIBILITY_SUPPORT_OPERAND` table (`…:2384`, all branches×densities). E5 MEASURES the real footprint.

## Governing invariants
- ⛔ ADDITIVE ONLY. With `S11CB_PRIMARIES_ONLY` UNSET, the engine's emit is byte-identical to HEAD (default-path
  invariance, verified per E3). The gate must not mutate shared iterators/helpers/cleanup/emit-ordering on the unset
  path.
- ⛔ SAME OBJECTS via the ENGINE'S OWN emit sites (`emitShared`/`modelRecord`/`kernelRecord`/origin helpers/the
  `ADMISSIBILITY_*` maps) on `mainModels`/`mainKernels` populated ONLY from `evaluatedModel["EULERIAN",branch,density]`
  (`…:2192`) + `extractCouplingData[…]["FINAL_KERNEL"]` (`…:2195`). ⛔ No hand-assembled tag payloads (E7).
- ⛔ No physics value is changed, asserted, or leaked (rule 5): this is an emit-scope gate, not a computation.

## What must be TRUE
E1. **Emit EVERY primary family the comparator/`row_residual` join** (verified against
    `scripts/S11c_b_cross_engine_comparator.py:84-95,1084-1096` and `scripts/S11c_b_row_residual.py:566-602`), for
    the FULL case set, so the WL case-set matches the PY `S11CB_PRIMARIES_ONLY` `.out`:
    `ENERGY_BASIS_VARIABLE`, `ENERGY_BASIS_COUNT`, `ENERGY_BASIS_NEW_INVARIANTS`, `ENERGY_BASIS_OMISSIONS`
    (branch-scoped, pre-loop — must stay emitted over all branches, NOT reduced); `SLAB_OPERATOR`,
    `SLAB_OPERATOR_TERM_ORIGINS`, `MU_THETA_OPERATOR`, `COUPLING_KERNEL`, `COUPLING_KERNEL_TERM_ORIGINS`;
    `ADMISSIBILITY_OPERATOR_OPERAND` (`row_residual:595` RAISES if unaligned — the fatal v1 omission),
    `ADMISSIBILITY_SUPPORT_OPERAND`, `ADMISSIBILITY_RESIDUAL`. The term-origin integrity records
    `SLAB_OPERATOR_TERM_ORIGINS_SUM_RESIDUAL` (`…:2278`) and `COUPLING_KERNEL_TERM_ORIGINS_SUM_RESIDUAL` (`…:2321`)
    are on the primary path and are RETAINED (classified as required integrity records, not controls).
E2. **Skip every non-primary path** (the memory hogs) — named as OBJECTS, not a line range: the tower-depth
    variants + tower-spool BUILD (`operatorTruncated/Extended`, `kernelTruncated/Extended`, `towerCase`, `…:2204-2258`)
    AND its CONSUMER emit `CONTROL_BACKGROUND_JET_TOWER_OPERANDS` (`…:2300-2311`); the tractability equivalence
    (`…:2241-2257`); the operator-freeze rank diagnostic (`…:2124-2128`); and every later control/ablation loop that
    re-runs `evaluatedModel`/`extractCoupling`/`frozenEvaluatedModel` — MATERIAL (`…:2438`), corrupted (`…:2466`),
    frozen (`…:2522`), form-ablation (`…:2742/2772`), uniform-limit (`…:2885`), rep-invariance/independence/
    homogeneity — i.e. the WL analogs of the PY `CONTROL_TASKS` (`sympy_audit.py:42-52`). After the last required
    primary emission the gate takes an allowlisted exit; no later block that dereferences a per-case control object
    (e.g. `towerSpoolPaths[key]`) may run.
E3. **Default-path invariance is VERIFIED, two ways** (⛔ not the v1 "demonstrate it's guarded" typed sentence, and
    ⛔ not a full 30 GB byte run which is ≥64 GB work): (a) STRUCTURAL — `git diff HEAD` shows every changed line
    inside a gate branch whose predicate is false when `S11CB_PRIMARIES_ONLY` is unset; the unset path's `Do`
    iterators stay `{branch, branches},{density, densities}` (`…:245-246` literals) and no shared helper
    (`emitShared`/`modelRecord`/`kernelRecord`/`evaluatedModel`/`extractCouplingData`) gains a gate argument. (b)
    PAIRED REDUCED-SCALE BYTE-IDENTITY — a HEAD copy and the patched copy, reduced by a MECHANICALLY IDENTICAL scale
    edit in both, each run with the gate UNSET, `cmp`-identical and both exit 0.
E4. **The residual-mode `.out` passes the COMMITTED extractors**, not merely `load_wl` tag-presence (`load_wl`
    accepts any `WL_S11CB_*:` line and never parses the payload, `comparator.py:177-209`). For each required family,
    the committed `extract_slab`/`extract_coupling`/`extract_mu`/`extract_energy`/`extract_admissibility`/
    `extract_term_origins` must produce cases with NO exception, the correct `(BRANCH,DENSITY)` (or branch) scope,
    and the required object/DOF leaves present.
E5. **The real footprint + wall time are MEASURED** by a memory-watched full primaries-only run (not asserted from
    STEP 0): record peak RSS, min MemAvailable, and exit status; require exit 0 (⛔ not 124/137), a predeclared
    MemAvailable floor of 2.5 GB, ONE kernel at a time (2-seat licence), and a concrete `timeout` set ABOVE the
    expected wall time (STEP 0 was ~26 min/case; four cases + origins ⇒ budget generously, e.g. ≥4 h — ⛔ do NOT
    inherit the review-legs 600 s, which is for structural ablations only). ⚠ If the measured footprint exceeds the
    box, the DOCUMENTED FALLBACK is a single-case emit + a P2 comparator case-intersection (do NOT silently ship a
    partial `.out`; the runner writes a temp artifact and promotes it only after exit 0 + E4 extractor validation).
E6. **The gate is a BOOLEAN env `S11CB_PRIMARIES_ONLY`** (mirroring PY; ⛔ no route/branch/density selector — that
    was the v1 mistake). Unset (WL `$Failed`) → full default run. Set to any string (per PY's `os.environ.get`
    truthiness; match the WL `StringQ[Environment[…]]` convention already used for `S11CB_SKIP_HEAVY_CONTROLS`
    `…:2177`) → primaries-only all-cases. Production route stays hard-coded `EULERIAN` (`…:2192`); no case is
    selected or dropped.
E7. **Faithfulness has an ACCEPTANCE TEST** (governing invariant 2): for each required tag, the residual-mode
    payload for a case EQUALS the HEAD full-run payload for the same case (from the E3 paired reduced-scale run),
    byte-identical per tag — proving the payload was produced through the existing `modelRecord`/`kernelRecord`/
    origin/`ADMISSIBILITY` maps on `mainModels`/`mainKernels`, NOT a hand-assembled association.

## Build execution note (operational — not a folded decision)
The BUILD implements the gate and verifies its CORRECTNESS at a TRACTABLE scale: the E3 paired reduced-scale
byte-identity, an E4 extractor pass and E7 per-tag equality on a REDUCED / smoke emit (e.g. a reduced basis or a
gate-plus-manual single-branch invocation that still exercises the real emit sites), and a memory sample of that
smoke. ⛔ Do NOT run the full ~4 h all-case production emit inside the build. The full-scale E5 production run — the
actual fresh 4-case WL `.out` generation, memory-watched — is a SEPARATE orchestrator step AFTER the build legs
clear and the engine is committed. State any check performed at reduced scale as reduced-scale.

## Legs
Orchestrator-written directive ⇒ this v2 folded both decision legs (no re-leg, rule 7). Codex build (WL,
`--sandbox danger-full-access`, xhigh) ⇒ 2 build legs (fresh Claude agent + Grok; ⛔ SERIALIZE — Mathematica 2-seat
licence; wrap structural-ablation kernel runs in `timeout 600`, the production emit run in the E5 budget; kill an
orphaned kernel by exact pid; check `free -h` first if a job dies). Ablations that must BITE: unset-env output
byte-identical to HEAD at reduced scale (E3); set-env output carries exactly the E1 primary families for all cases
and NO control/tower payloads (E2); every required tag passes its committed extractor (E4); a required-tag payload
equals HEAD's for that case (E7).
