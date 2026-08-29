# canon_jet_name time-order fix — BUILD review, two legs, consolidated (SOUND, NO REGRESSION)

**Artifact:** the `canon_jet_name` fix in `scripts/S11c_a_cross_engine_comparator.py` (Boolean `has_time` →
`time_order` counter; `tt` added to the derivative-token set) + pinned fixtures in
`scripts/test_S11c_b_cross_engine_comparator.py`. **Legs (Codex-written ⇒ fresh Claude agent + Grok):** prompt
`_legs/S11c_b_canon_jet_time_order_build_review.md`. Raw: Grok `~/.s11_build/S11c_b_jetfix_build_grok.txt`;
Claude-agent scratch `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-.../scratchpad/opus_*`.
**Both legs: SOUND, NO BLOCKING FINDINGS.** Both ran the ablation + the full S11c-a regression in the
foreground.

## Convergent verifications
- **Fix correctness + FORM ablation (load-bearing):** both legs reverted the counter to the Boolean `has_time`
  form (ablated /tmp copy) and the pinned fixtures FAIL (`u_1_t_t` → `u_1_t` ≠ `u_1_tt`); post-fix WL `*_t_t*` and
  PY `*_tt*` canonicalize to the identical string. The counter is the load-bearing change (a coefficient tweak
  would not do it). Claude: `opus_form_ablation_counter.out`; Grok: `03_run_pinned_fixtures.out`.
- **No spatial / order-≤1 regression:** every spatial-only and order-1 name byte-identical pre→post; the fix is
  injective (0 collisions) and UN-MERGES the pre-fix false collapse (`A_T_1_t_d2d3` order-1 vs
  `A_T_1_tt_d2d3` order-2 were merged pre-fix, now distinct). Claude: `opus_injectivity.out`; Grok:
  `04_vocabulary_regression.out`.
- **S11c-a regression (the key check — neither leg trusted the builder's "4 diffs"):** both ran the S11c-a
  comparator PRE (from `git show HEAD`) and POST to completion (~2500s each) and diffed. Both found ~28
  differing lines and independently computed each residual to ZERO (Claude 28/28 canonically zero, 0 nonzero, 1
  `runtime_seconds`-only; Grok same). `RUN_ACCOUNTING`/`MEASUREMENT_SCOPE` byte-identical. Root cause of the
  28-vs-4 discrepancy (both legs): the comparator's `--residual-leaf-budget 0.1s` fallback makes residual
  SERIALIZATION timing-dependent under CPU contention — pre-existing, orthogonal to the fix. Corroborated
  structurally: S11c-a carries ZERO order-2 time content (all 199 PY symbols identical pre/post canon; WL has no
  4-slot time-order≥2), so the fix CANNOT change any S11c-a canonical name. NO REGRESSION.
- **Coupling asymmetry exposed, not massaged:** the WL `Derivative[0,1,1,2][transverseTrialPotentialOne]`
  (864×) now canonicalizes to `A_T_s11cb_1_tt_d2d3` (was `_t_d2d3`); the regenerated reconcile carries
  `u_i_tt`/`e_W_tt` in `SLAB_OPERATOR` (8 each). PY zeros the dynamical `u_tt`/`e_tt` in `sector_substitution`
  (`S11c_b_brane_operator_sympy_audit.py:1659`, symbols at 708-709). The fix is an injective relabel applied
  identically to both engines, so it CANNOT manufacture false agreement — it corrects the pre-fix false merge.

## Non-blocking observations (carry to Phase 2 adjudication)
1. ⚠ **The transverse-trial ∂²ₜ terms do NOT appear in the printed reconcile `A_minus_B` residual** (absent under
   both PRE and POST — PY carries `A_T_s11cb_i` only as base symbols; the ∂²ₜ may legitimately cancel). Both
   legs: this is NOT the fix silently massaging the asymmetry; the instrument now represents the true order, and
   whether the coupling ∂²ₜ is a real cross-engine finding or a legitimate cancellation is the ADJUDICATION's
   call. ⛔ Do not pre-claim it a finding. The SLAB inertia (`u_i_tt`/`e_W_tt`) DOES surface and is in scope.
2. ⚠ **The comparator's `--residual-leaf-budget 0.1s` timing fallback makes residual serialization
   non-deterministic run-to-run on heavy families** — pre-existing, orthogonal to this fix, but relevant to
   adjudication REPRODUCIBILITY (run the adjudication with a generous or disabled leaf budget, or expect
   serialization churn on the heavy residual families).

## Verdict
Both build legs SOUND, no blocking findings. Fix repairs the asymmetric time-order collapse correctly and
minimally; S11c-a not regressed; S11c-b inertia exposed rather than collapsed. Commit.
