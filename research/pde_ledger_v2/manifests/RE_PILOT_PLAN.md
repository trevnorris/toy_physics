# Integration-test re-pilot — RESUME ANCHOR (next session)

> **Status: the framework is BUILT + COMMITTED (`e849e303`).** This doc is the resume anchor for the NEXT
> step: re-extract the two pilot slices under schema v2.1 and run the real checker. READ-FIRST next session,
> then `manifests/EXTRACTION_PROTOCOL.md` + `manifests/MANIFEST_README.md` + `manifests/stage_manifest_schema_v2.json`.

## ⭐ RE-PILOT PROGRESS (2026-07-24, this session — read this first)
**MILESTONE ACHIEVED + independently verified:** stages 030, 031, 032 extracted to v2.1; the 3-stage
`{030,031,032}` composite build is CLEAN (exit 0, headline PARTIAL — all findings honest edges into unextracted
stages). `citations=17/20` resolve; `claims=106/106`; **ADJUDICATION PASS** (stage032 23040-cell grid: buckets
sum + axes `5·6·6·4·2·2·2·2·2` both = 23040); the 030→031 spectrum handoff resolves (`spectrum_match` on the
corrected `O_perp` operator, able-to-fail confirmed). 030 + 031 fidelity-audited FAITHFUL (031 audit pending on
032). The checker was HARDENED (func/composite dim recovery) → Grok-reviewed (4 holes found+closed) → refined
(over-correction on real stage030) → CONVERGED at **49 self-test fixtures** (SHA `faa7e8f1de1d3755`, frozen for
extraction). Extraction pivoted from Codex → **agent-as-coder** (Codex CLI stalled twice on API hangs; agents run
in-harness, reliable — this is the likely fanout default, see FANOUT_PLAN.md).

**The pilot EARNED ITS KEEP** — real integration errors it caught: 031's v1 `det_m` axis handling, the spurious
`stage030/kernel_determinant_star` citation (→ real `kernel_stability` D*=7/4), the drifted 031-claim renames that
032's v1 consumes still pointed at.

**SYSTEMATIC FINDINGS (fold into #27 / the fanout approach — USER input wanted on #2):**
1. **Checker gaps to fix in ONE consolidated round (#27) before fanout:** (a) `dimension_of` doesn't treat sympy
   `pi`/NumberSymbol as dimensionless → π in a dim expr = `DIMENSION_RULE_UNSUPPORTED` (stage031 worked around with
   a `declare_pi` dimensionless local symbol); (b) dim-recovery requires a live `Dim` class → stage032's script uses
   BARE TUPLES (`dim_E=(2,-2,1)`) so its C4 certificate had to BORROW stage031's script (values right, mis-anchored).
   Fix both + able-to-fail fixtures; keep all 49 fixtures + real 030/031/032 passing; then re-anchor 032 dims + drop
   the declare_pi workaround.
2. **Systematic UNDER-EXPORT:** stages don't export every claim a downstream stage cites (030 lacked `O_perp` — added;
   031 lacks `S_gg`/`declare_sgg` etc. → 032 fell back to opaque citations + an `A_V=m_gg·C_V` opaque form vs the note's
   `m_gg·φ²/S_gg²`). The fanout likely needs an EXPORT-COMPLETENESS pass (a stage exports every claim downstream cites).
   ⚠ This is a design choice affecting all 44 stages — confirm the approach with the user at the fanout opt-in.

**REMAINING pilot work:** the real 030→031→032 **C7 mutation** (task #23, the capstone anti-rederivation proof —
manifests currently have `c7_edges=0/17`, no C7 metadata; author an HONEST mutator + add c7_binding/c7_expect to one
clean edge, e.g. 030→031 `zero_mode`);

> **⭐ C7 EXECUTION RECIPE (contract read from `composite_build.py:1870-1946`, committed `30a432e9`):** The harness
> covers an edge when the producer export has `c7_binding{producing_primitive, mutation_env, mutation_command,
> exported_facet}` AND the consumer consume has `c7_expect{injection_point, facet_used, expected_first_failure}` AND
> `facet_used == exported_facet`. It runs `mutation_command` as a subprocess with `env[mutation_env]=exported_facet`,
> and reads the LAST stdout line as JSON `{consumer_stage: tooth_or_"PASS"}`. Verdicts: consumer `"PASS"` →
> `DECORATIVE_DEPENDENCY` FAIL; `!= expected_first_failure` → `WRONG_FIRST_FAILURE` FAIL; `== expected_first_failure`
> → PASS. Any OTHER loaded stage firing → `UNDECLARED_DEPENDENCY` FAIL. `expected_first_failure` MUST be a tooth in the
> consumer's `verification.teeth`. **HONESTY (non-negotiable):** the mutator must ACTUALLY perturb the facet + COMPUTE
> the downstream tooth outcome, never hardcode the tooth. Concrete plan for the 030→031 `zero_mode` edge: mutator
> reads `C7_FACET`, perturbs f0's power 2→1 (`f0=1/(ell*cosh(w/ell))`), computes the residual
> `-diff(f0,w,2)+V_H*f0` in sympy (≠0 ⇒ the zero-mode property genuinely breaks), and emits
> `{"stage031": "<the 031 tooth whose assertion consumes the zero mode>"}` iff the residual is nonzero, else PASS.
> Then add the c7 metadata to 030's `zero_mode` export + 031's `zero_mode` consume, run the checker → that edge's C7
> goes PASS; also prove a `--decorative` variant (ignores the facet, prints PASS) → `DECORATIVE_DEPENDENCY` FAIL, per
> the self-test fixtures at `composite_build.py:3536-3548`. Do this in a FRESH context (it's the subtlest piece). 043 (record_range path — the one distinct checker path still untested) + 006/044
(lower marginal value — candidates for early fanout instead); the #27 checker round; commit. cleanup list: stage031
`declare_pi` + its `det_m` extraction-report narrative (v1 was `[L,M,T]` order, not an error).

## What exists (committed `e849e303`, all under `research/pde_ledger_v2/manifests/`)
- `stage_manifest_schema_v2.json` — schema **v2.1** (`schema_version` const "2.1"). The per-stage machine interface.
- `composite_build.py` (~100KB) — the checker: IMPORT-COMPLETENESS + C1–C6 + the **C7 mutation harness**.
  Run: `python3 manifests/composite_build.py` (on `manifests/stages/*.json`) or `--self-test`.
- `EXTRACTION_PROTOCOL.md` — the reusable per-stage extraction protocol (v2). **Used by every extraction.**
- `MANIFEST_README.md` — C1–C8 check semantics + conventions.
- `DIM_ORDER_DECISION.md` — named-dict `{L,M,T}` dims + per-stage source-order table (scripts are MIXED; 030/031/
  032/044 = LTM, 027 = LMT).
- `composite_config.json` — pins the `parameter_register.md` path+digest for C5.
- `examples/stage_manifest_v2_example.json` — a schema-v2.1-valid worked example.
- `stages/` — EMPTY (`.gitkeep`); v2.1 manifests land here.
- **Verified:** 34/34 able-to-fail self-tests (independently reproduced, exit 0); genuineness audit `YES`
  (8 novel author-unscripted probes: real CAS, real file-re-read for the C4 certificate, real subprocess C7);
  the schema-invalid reporting bug fixed (skipped checks report `SKIPPED`, never `PASS`).

## Reference (NOT committed — under `_scratch/pilot_v1_manifests/`)
The four **v1-schema** pilot manifests (030/031/032/044) + the stale v1 example + the v1-pilot run report. These
informed the design. ⭐ The re-extraction can START from these (upgrade v1→v2.1) rather than from scratch — faster
and they already captured the hard content (032's 23040-cell adjudication, 044's roster + S_hold map, etc.).

## THE RE-PILOT (next session's work)
Two slices (Codex's recommendation — they exercise the hardest paths):
- **Slice A `006→030→031→032`** — convention/citation/operator/spectrum/R1 behavior.
- **Slice B `043→044`** — census/range-object/verdict/structural + count sensitivity.
So re-extract **6 stages**: 006, 030, 031, 032, 043, 044 (030/031/032/044 = upgrade the v1 manifests; 006/043 = new).

Steps:
1. Re-extract all 6 under `manifests/EXTRACTION_PROTOCOL.md` → `manifests/stages/stageNNN.json` (schema v2.1). Parallel
   agents; each self-validates against the v2.1 schema. Honor the v2.1 rules: named-dict dims + `dim_source_order`
   + `dim_source_tuple`; `quantity_id` OWNERSHIP (a quantity's first-defining stage owns it; downstream references
   it, never re-declares `here`); typed payloads; structured `consumes` + `substitutions`; `c7_binding`/`c7_expect`;
   per-item `evidence` {path,locus,digest,engine,method}; genesis evidence for `independent`.
2. Run `python3 manifests/composite_build.py` on the 6 → the FIRST genuine integration report on ledger content.
   Expect PARTIAL on edges into not-yet-extracted stages (that's honest, not a fail). Triage every real FAIL —
   these are genuine cross-stage findings (symbol collisions, citation drift, dim breaks, miscounts). Fix/record.
3. **One real C7 mutation on `030→031→032`:** author an HONEST mutator (`mutation_command` that perturbs behavior
   based on `C7_FACET` — see the mutator-honesty note in EXTRACTION_PROTOCOL/CHECKER directive) for a 030 exported
   facet; confirm 031/032's declared `expected_first_failure` teeth fire (and that a deliberately-decorative
   consumer is caught as `DECORATIVE_DEPENDENCY`). This proves the anti-rederivation guarantee end-to-end on real
   stages.
4. Commit the v2.1 pilot manifests + the integration report + any fixes.

## After the re-pilot
- **The 44-stage fanout** — extract the remaining ~38 manifests. ⭐ This is the point to greenlight a WORKFLOW
  (user opt-in): a pipeline per stage (extract → self-validate → composite-check). ~40 stages × (primary + blind
  consistency rebuild). Then run the full composite build → the ledger-wide integration report.
- **C8 + the prose/script consistency review** (the reframed "blind rebuild"), then wire `composite_build.py --self-test`
  + a partial composite run into the per-commit flow so integration stays live.

## Process wins to carry (validated this build)
- Reviews drive correctness: **Codex → Grok → Codex** on foundational directives; the ACCEPTANCE is MECHANICAL
  (self-tests proving each evasion FAILS), not subjective — "truth in output, not prose."
- ⛔ NO can't-fail conjuncts / no hollow green: every check must be able-to-fail (planted-defect fixtures) + a
  genuineness audit that drives the checker with NOVEL mutations (not just the built-in fixtures).
- Codex builds/codes; Claude reviews + owns scaffolding/directives; independently re-run every script (don't trust
  self-reports). setsid-detached launches; short waiters re-armed (long background waiters get reaped ~10–15 min).
- Feed reviewers LEAN inputs (the v1 4672-line log blew a Codex budget — point at distilled issues, not raw logs).
- `/var/projects/toy_physics` NOT `toy_projects`.

## ⭐ PAUSED ledger threads (do NOT lose — resume after the integration-test work reaches a milestone)
The integration-test framework was a (correct, user-directed) DETOUR. Two ledger threads are paused:
- **044-v2 (redo stage 044 with a dynamical-Σ sleeve — un-freeze `S_hold`):** see `notes/stage044_v2_unfreeze_prep.md`.
  User DECIDED: redo 044 (not a downstream patch), commit the bending/anchoring knobs `κ_bend/κ_anchor/collar-tension`
  (Tier-2), un-freeze at the action level (dynamical free boundary; full free-boundary solve stays sim-deferred). The
  reopen is **044-local** (locality confirmed by fresh reads AND mechanically by the stage044 manifest's S_hold-map:
  S_hold-dependent = `s_hold`/`action_roster`(6→5)/`field_set_union`(λ_Σ removed, pinned mode → physical DOF)/dim-
  homogeneity term count; everything else — dedup `r_B≡χ_B`, `Z_χ`, wave speeds, drain, P-retirement — survives).
- **045 (VII-2b, the non-variational drain/return block):** `notes/stage045_nonvariational_block_prep.md`. Drain =
  the DYNAMICAL `Γ_B` (frozen-wall ruled out); re-host the G0 §6 return/BC/force-partition STRUCTURE on it.
The manifest system now SUPPORTS these (mechanical locality + knob census). Likely order next: finish integration
re-pilot (+ maybe fanout) → 044-v2 → 045. USER to confirm sequencing.
