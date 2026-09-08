# S11c-c2 N6 ablation-harness knife-list directive — G2 review record (review-until-clear, CLEARED)

Directive: `directives/S11c_c2_N6_ablation_harness_directive.md` — the fixed four-engine knife list for the N6
ablation harnesses (the durable per-engine cert of the diagnostic/reconcile/covariance SymPy engines + the blind WL
N6 engine `e11f2f82`). Physics-bearing ⇒ review-until-clear.

## Round 1 — orchestrator-written draft → BOTH legs NOT-SOUND (Codex-sol + Grok)
Reports: `scratchpad/n6_dir_{grok,codex}_report.md` (session scratch). Convergent NOT-SOUND, ~9 grounded findings —
systematic engine-internal SITE / data-flow errors + leak language:
- Wrong sites / crashes: K_circular (`predicted_amplitude` scope can't reach the target); K_rank (`background_depth`
  is the wrong object; a BFS-only truncate trips `if uncovered: raise`); K_source_route (type crash — `build_increment`
  returns a tuple, needs unpacking); K_advection (at production `t=1` the `shift` term is identically 0 → moves nothing).
- FORM-vs-COEFFICIENT mislabels: K_ewsign (W→−W = ×−1) and K_operand_swap (ms↔es = −ds) were labeled FORM.
- Cross-file two-sidedness: K_normal's global S11c-a patch also moves `ms` via `face_velocity_raw`.
- Gaps: nothing corrupted the WL SPLIT_CHECK identity; covariance omitted the rank-certification tags; K_carrier's
  DEAD set wrongly included `R_COV_INCREMENT` / `ACTUAL_CONTROL_PARAMETERS`.
- Leak language ("MUST be 0", "must leave byte-identical", "print that the checker fires").

## Orchestrator adjudications on the 3 divergent/judgment points (E2 FORM-vs-COEFFICIENT)
1. **K_ewsign is COEFFICIENT** (sign-flip = scalar ×−1) — Codex right, Grok wrong. The diagnostic FORM knife becomes a
   genuine row-family removal (`K_EW_rowdrop`); sign-flip demoted to the coefficient companion.
2. **K_operand_swap is COEFFICIENT** (ms↔es = −ds) — replaced by a structural `es−es` operand collapse (211-212).
3. **K_advection DROPPED** (0 at production `t=1`; it only nullifies the engine's own N4 control) — certify instead via
   the engine's own `CONTROL_INDEPENDENCE {BASE,CORRUPTED,RESIDUAL}` triples for TILT + N4.

## Re-author — delegated to Codex-sol (rule 15 / the N6 route-2 precedent)
Given round 1 exposed systematic engine-internal site errors in the hand-authored draft, the re-author was delegated
to Codex-sol (with both leg reports + the adjudication key + the 3 adjudications above). Codex re-authored the
directive in place AND verified the four high-risk mutations (row-drop, rank-image, slot-rescale, scoped-normal)
actually run on temporary engine copies before finalizing ("all pinned Python mutations parse and the high-risk forms
completed on temporary engine copies"). Codex-written ⇒ the round-2 legs are a fresh Claude agent + Grok.

## Round 2 — re-authored directive → BOTH legs SOUND (fresh Claude Sonnet + Grok)
Reports: `scratchpad/n6_r2_grok_report.md` (Grok, with check script `/tmp/n6_r2_review/knife_checks.py`); the fresh
Claude Sonnet leg verdict + per-knife confirmations (below). Both legs traced every knife against the live engine
code and confirmed:
- **Every fragment unique** (Grok's script: count 1 for all 11 sites); **every patch parses/runs** (PARSE_OK for
  K_EW_pop / K_rank_form / K_normal_wrapper / K_source_unpack / K_slotdrop_coeff); **every API exists**
  (`a.grad_W`, `a.dot`, `n.replace`=`dataclasses.replace`, `a.build_material_face_source`, `a._FACE_CACHE`).
- **No crash**: `rows.pop('E_W')` is safe (build_increment's output-domain completion is name-independent; `slots`
  unused in its body); K_rank's `image_of` patch keeps coverage (no new `uncovered` raise).
- **FORM-vs-COEFFICIENT correct** throughout; K_EW_rowdrop is a real row-family removal, K_operand_swap the `es−es`
  collapse (not the ×−1 swap), K_split_route genuinely admits the affine sig-0 family into the source channel only.
- **Cones correct**: each FORM knife reaches its certified object; every DEAD object genuinely independent
  (K_carrier DEAD correctly excludes `R_COV_INCREMENT`/`ACTUAL_CONTROL_PARAMETERS`).
- **K_normal carrier-only**: `m_v = build_material_velocity(...)` runs (rec:204) BEFORE the wrapper installs
  (rec:205), so `ms`/`es`/the source bridge use the unpatched builder; wrapper + cache restored in `finally`;
  Eulerian (a:850) / global / `material_inverse_transpose` (a:694) sites untouched.
- **No leak** (0 hits for PASS/FAIL/must-vanish/non-trivial/verdict); builder bounds clean (build→run→report→stop,
  no other AI, no commit).

Fresh Claude leg verdict: **SOUND**, 2 non-blocking observations. Grok verdict: **SOUND**, no finding survives the
physics filter. (The Claude leg disclosed it incidentally saw the adjudication-key filename in a directory listing
but did not use it — reviewers are not the blindness target; it derived everything independently.)

## Verdict — CLEARED
Both round-2 legs SOUND; nothing outstanding changes what is computed or may be claimed. One cosmetic prose nit
folded (K_rank rank≥3 recursion wording; code was already correct per both legs). The WL K_junk case-scoping
observation needs no fix (H1 runs all four cases). The directive is build-clear. Adjudication key
(`_measurements/S11c_c2_N6_ablation_harness_adjudication_key.md`, ORCHESTRATOR-ONLY) updated to the cleared knives.
Committed as the reviewed baseline; astra builds the four harnesses against it (handed the directive alone).
