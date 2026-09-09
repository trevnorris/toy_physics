# S11c-c2 N6 cross-engine comparator — SCRIPT review record (review-until-clear → SOUND)

## Artifact / role
`scripts/S11c_c2_N6_cross_engine_comparator.py` (+ `test_…py`) — astra-built (`gpt-6-astra` high) against the cleared
build directive (`1228e357`). Cross-engine N6 measurement instrument (joins the WL N6 engine `.out` vs the 3 SymPy N6
`.out` streams by object name; three-valued; PIT/digests sealed; PRINTS/decides nothing).

## Legs (Codex/astra-written → fresh Opus agent + Grok; identical prompt; parallel — pure Python, no kernel)
- **Fresh Opus agent** — report `directives/_legs/S11c_c2_N6_comparator_review_opus.md`.
- **Grok** (`grok-4.6` high) — report `directives/_legs/S11c_c2_N6_comparator_review_grok.md`.
- Prompt `directives/_legs/S11c_c2_N6_comparator_review_prompt.md`. Both derived an independent join contract from the
  engines' committed streams + constructors BEFORE opening the comparator, then ablated /tmp copies (mandatory FORM).

## Verdict: BOTH SOUND — convergent, ablation-backed. No finding survives the physics filter.
The two legs independently confirmed every load-bearing decision:
- **PIT SEAL holds:** forced cross-engine PIT subtraction in a /tmp copy DOES emit `PIT_A_MINUS_B`; the working
  comparator emits 0 (Opus: 868 pit_sealed / 0 residual; Grok: forced-subtraction ablation). Support derives from the
  boolean witness, never enters `residual()`.
- **SCHEMA BRIDGE (the load-bearing B adjudication) is correct + correctly implemented:** block-name / kernel-family /
  wave / component keys decode to DISJOINT axis sets → `axis_set_mismatch`/`unmatched_key` + pairing-table rows, ⛔ never
  auto-joined. FORM-collapsing those axes in a /tmp copy MANUFACTURES a join (proving the check bites); the working
  comparator has no such path. The conservative rule-5-safe pre-registration (name/CAS identities only: grade
  `{1,η,σ}↔(η,σ)`, face `SUM↔0`, jet/U-origin spelling) + pairing table was the right call.
- **Residual genuinely computed (not an echo):** Opus spurion FORM ablation moved all 160 matched R_COV off zero;
  Grok two-slot C_E ablation moved both residuals; one-sided corruption moved ONLY the corrupted slot; repoint moved
  the residual. The all-zero R_COV joins are real algebraic cancellations.
- **Carrier-rep computes a residual (DO-NOT-FOLD, ⛔ not sealed):** C_E 320 joins → 40 genuinely nonzero symbolic
  residuals surfaced.
- **grade fold is not an M3 freeze:** the dropped leading element is constant=1 across all 288 real WL columns (verified).
- **_NODES/formal_jet round-trip** on the real committed DAG; missing `ARITHMETIC_DAG` → `parse_failed`.
- **Join-set complete + no leak:** the 4 added covariance families present; SymPy-only → `sympy_only`/"no WL sibling";
  `residual_target=none`; no prime literals; `ACTUAL_CONTROL_PARAMETERS` field-name mismatch (kappa_a vs ADVECTION) is
  NOT aliased (unmatched = the measurement).
- **DoD teeth:** zero-extract from two nonempty shared containers → exit 2; `join>0` never the exit; all-deferred → exit 2.
- **No banned verdict token, no assert** on measured payloads. Tests 23/23.

## Orchestrator G4 verification (mechanical, own)
Ran the committed test suite: **23 passed**. Grepped: 0 printed `PASS/FAIL/AGREE/DISAGREE/COVARIANT/VERDICT` tokens; 0
asserts in the comparator (asserts live only in the test file). Concur with both legs.

## Non-blocking observations (correct behavior — no change required)
`dimension_record` always `KEY_DISAGREE` (different engine schemas; the signal is in the sibling `dimension_vectors`);
`addition_consistency` always `BooleanNotResidualable` (required rejection of a native bool). Both honestly surfaced.

## Builder overstep (quarantined — not used)
astra ran its OWN grok/opus "reviews" + a "review_adjudication" (`_measurements/S11c_c2_N6_comparator/`), a recurrence of
the c1 build overstep ([[feedback_builder_directive_no_orchestrator_process]]). Quarantined out of the tree to
`scratchpad/c2_n6_comparator_builder_overstep/` BEFORE the real legs ran (⛔ not read, ⛔ not used); the real legs are the
fresh independent Opus + Grok above. ⚠ Future comparator directives should carry the explicit build→run→report→stop
fence (the builder ran its own reviews despite the build-skill guardrail).

## Disposition
**SOUND → the N6 cross-engine comparator is the CLEARED instrument.** Stopping rule (G4): nothing outstanding changes
what is computed or may be claimed. Commit the reviewed baseline (script + test + reports + this record). NEXT = run the
comparator on the committed streams → the staged representational reconcile (surface, ⛔ not pre-adjudicate, the carried
cross-engine questions) → the c2 step record.
