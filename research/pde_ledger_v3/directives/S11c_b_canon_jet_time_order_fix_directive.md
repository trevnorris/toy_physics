# Build brief — repair `canon_jet_name` time-order collapse in the shared comparator

## The defect (mechanism, verified)
`canon_jet_name` in `research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py` (~L800-811) records
time differentiation as a BOOLEAN and emits a single `_t` regardless of order:
```
has_time = False
for token in derivatives:
    if token == "t": has_time = True
suffix = "_t" if has_time else ""
```
`jet_suffix_from` (~L610-641) correctly preserves the count — WL `D[…,{time,2}]` produces two `t` codes
(`…_t_t`) — but `canon_jet_name` then collapses `u_1_t_t` → `u_1_t` (order 2 → order 1). PY emits second time
derivatives as a single `tt` token (`u_{a}_tt`, `e_W_tt`;
`scripts/S11c_b_brane_operator_sympy_audit.py:249,252`), which `canon_jet_name` leaves as base (`tt` is not a
recognized derivative token) so PY `u_1_tt` stays `u_1_tt`. Result: WL order-2 time jets canonicalize to order
1 while PY's stay order 2 — an ASYMMETRIC loss that manufactures false cross-engine differences on every inertia
(second-time-derivative) term. Confirmed in the current reconcile log: `e_W_tt`/`u_i_tt` residuals in
`SLAB_OPERATOR` (THICKNESS_ROW / U_MOMENTUM_ROWS). ⚠ The defect ALSO reaches `COUPLING_KERNEL`: the WL transcript
carries second-time mixed derivatives of the transverse trial potentials
(`Derivative[0, 1, 1, 2][transverseTrialPotentialOne]` and siblings, 864 each;
`mathematica/out/S11c_b_brane_operator_mathematica_audit.out`), which the defect collapses to first order — so
the coupling residual currently shows `A_T_s11cb_*_t_*` where the WL payload actually held `*_t_t_*`. The fix
WILL change the coupling-kernel comparison and may EXPOSE a genuine cross-engine asymmetry (WL keeps ∂²ₜ of the
transverse trial; PY zeros `u_tt`/`e_tt` in its sector restriction, `S11c_b_brane_operator_sympy_audit.py:696`).
That exposure is the correct behaviour — a disagreement is a finding, never something this fix may suppress.

## Object
Repair `canon_jet_name` (and only what that requires) so the canonical jet name PRESERVES time-derivative
multiplicity and is SYMMETRIC across the two engines' spellings: a field differentiated N times in time
canonicalizes to the SAME name PY uses for N time derivatives — `_t` for order 1, `_tt` for order 2 (PY's
spelling) — so WL `u_1_t_t` and PY `u_1_tt` map to the identical symbol, and likewise `e_W_t_t` ≡ `e_W_tt`.
Recognize PY's single `tt` token as order-2 time as well, so both spellings converge. ⛔ Leave the spatial-jet
handling (`dN`, `dw`, direction sorting) and the non-time base parsing UNCHANGED — this fix touches only the
time-order branch.

## Constraints
- The change is to the SHARED comparator used by S11c-a and S11c-b. It must not alter any canonical name that
  has no time derivative, and must not alter order-1 time jets (`u_1_t`, `e_W_t`, `theta_t` stay `…_t`).
- This is an instrument fix, not an emit script. The engine RUN outputs still PRINT and never assert on measured
  payloads. Assertions are permitted ONLY in the synthetic canonicalizer fixture below (synthetic inputs, not a
  measured payload).

## Definition of done (value-free; the build legs check empirically)
1. **Synthetic canonicalizer fixture** (add to `scripts/test_S11c_b_cross_engine_comparator.py` or a sibling
   test): assert `canon_jet_name` maps to the EXACT canonical string (pairwise equality alone is not decisive —
   a boolean-plus-`tt` fix makes pairs equal at the wrong collapsed form and would pass). Pin the values:
   - `canon_jet_name("u_1_t_t") == canon_jet_name("u_1_tt") == "u_1_tt"`
   - `canon_jet_name("e_W_t_t") == canon_jet_name("e_W_tt") == "e_W_tt"`
   - `canon_jet_name("u_1_t_t_d1") == canon_jet_name("u_1_tt_d1") == "u_1_tt_d1"` (mixed space+time keeps both)
   - the coupling-shaped case: `canon_jet_name("A_T_s11cb_1_t_t_d2_d3") == "A_T_s11cb_1_tt_d2d3"`
   - order-1 unchanged: `canon_jet_name("u_1_t") == "u_1_t"`
   - spatial-only equal to the PRE-FIX `canon_jet_name` output (NOT raw input):
     `canon_jet_name("mu_R_bg_d1_d2") == "mu_R_bg_d1d2"`, `w1_profile_d1`, `theta_d1` unchanged.
   A no-op fails the pins; a `_t_t`- or `_t2`-form fails them; a boolean+`tt` form fails `== "u_1_tt"`.
2. **No spatial regression:** every spatial-only and order-≤1 jet name canonicalizes exactly as before the fix.
3. **No S11c-a regression:** re-run the S11c-a comparator on its committed transcripts
   (`out/S11c_a_interface_geometry_sympy_audit.out` + the WL default) and diff its output against the
   pre-fix output — it must be unchanged except for any intended time-order canonicalization (report the diff;
   if S11c-a carries no multi-time jets, it is byte-identical).
4. **Regenerate S11c-b artifacts:** re-run the v1 reconcile
   (`scripts/S11c_b_handcoded_comparison.py`) against the committed transcripts to
   `~/.s11_build/S11c_b_reconcile_run.out`; it completes at exit 0. Report the per-family MATCH/FLAG/
   NAMESPACE tally and a diff of the newly-visible order-2 time jets vs the pre-fix log (which families/keys now
   carry `_tt` where they carried `_t`) — including COUPLING_KERNEL, where the transverse-trial ∂²ₜ terms will
   surface. ⛔ Do NOT suppress or "reconcile away" any resulting change, and do NOT state an expected residual
   value or tally for any family — the adjudication measures those; this step only exposes the true jets.
5. All existing comparator fixtures still pass.

## Builder report (≤20 lines)
The `canon_jet_name` diff (the time branch, before → after); the fixture results (which pairs now canonicalize
equal); the S11c-a regression diff (empty, or only time-order lines); the regenerated reconcile tally + runtime;
explicit confirmation that spatial jets and order-≤1 time jets are byte-identical to pre-fix. No expected
residual value; no `PASS`/`FAIL`/`VERDICT`/target on any measured payload.
