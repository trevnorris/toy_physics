# Independent review — canon_jet_name time-order fix directive (decision list, pre-build)

## Artifact
`research/pde_ledger_v3/directives/S11c_b_canon_jet_time_order_fix_directive.md` — orchestrator-written brief to
repair a time-order collapse in the SHARED comparator `scripts/S11c_a_cross_engine_comparator.py`. This is a
DECISION-LIST review BEFORE the builder launches (no code yet). The fix changes a canonicalizer used by BOTH
S11c-a and S11c-b, so a wrong fix can silently corrupt cross-engine comparison. Find defects now.

## Context (read all)
- The directive above.
- `scripts/S11c_a_cross_engine_comparator.py`: `canon_jet_name` (~L785-812) and `jet_suffix_from` (~L610-641).
- PY second-time-derivative spellings `scripts/S11c_b_brane_operator_sympy_audit.py:214,217,248-252`.
- WL time derivatives `mathematica/S11c_b_brane_operator_mathematica_audit.wl:837-838,1658-1660`.
- The run showing the artifact: `~/.s11_build/S11c_b_reconcile_run.out` (grep `e_W_tt`).

## Required checks — report a finding only if it changes what gets built or what may be claimed
1. **Defect diagnosis.** Confirm from source that `canon_jet_name` collapses N time tokens to a single `_t`
   (Boolean `has_time`), while `jet_suffix_from` correctly emits N `t` codes. Confirm PY emits order-2 time as a
   single `tt` token and order-1 as `t`, and that `canon_jet_name` currently leaves PY `u_1_tt` as base (so WL
   `u_1_t_t`→`u_1_t` ≠ PY `u_1_tt`). If the diagnosis is wrong, say so.
2. **Named canonical form.** The directive says canonicalize order-2 time to `_tt` (PY's spelling) so WL
   `u_1_t_t` ≡ PY `u_1_tt`. Verify `_tt`/`_t` are exactly PY's spellings (not `_t2`, not `_t_t`). Check for a
   spelling where the two engines would STILL disagree after the fix (e.g. a base with an existing `t` in its
   name, a field whose PY form is not `tt`, a time+space mixed jet whose ordering differs).
3. **Scope / regression.** The fix must touch only the time-order branch. Check nothing in the proposed scope
   changes a spatial-only jet (`w1_profile_d1`, `mu_R_bg_d1_d2`, `theta_d1`) or an order-≤1 time jet. Is the
   S11c-a regression check (re-run S11c-a comparator, diff) adequate to catch a shared-canonicalizer break, and
   does S11c-a actually carry (or not) multi-time jets — check its engine source if reachable.
4. **DoD able to fail + no leak.** Is the synthetic fixture decisive (would a no-op or a wrong-form fix fail
   it)? Is the DoD value-free (no expected residual value for the inertia families; the COUPLING_KERNEL "tt-free
   / unchanged" statement is a verified structural invariant, not a leaked measurement — confirm it is not a
   physics-value leak)? Flag any leaked expected value or any DoD clause a broken fix would pass.

## Method
Read `canon_jet_name`, `jet_suffix_from`, and the two engine spellings FIRST; trace by hand what
`canon_jet_name` returns for `u_1_t_t`, `u_1_tt`, `u_1_t`, and a mixed jet, before and after the proposed fix.
Show those traces. A prose "looks fine" is worth nothing. Numbered findings (blocking vs non-blocking), file+
line, concrete correction. If sound, say so and name what you verified. ⛔ Reading + tracing only; do not spawn
Mathematica.
