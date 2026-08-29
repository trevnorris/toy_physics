# Independent build review — canon_jet_name time-order fix (SCRIPT change to the shared comparator)

## Artifact
The uncommitted change to `research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py` — `canon_jet_name`
now counts time-derivative tokens (`time_order`) and emits `"t"*time_order` (order-2 → `_tt`), and recognizes
PY's single `tt` token — plus the new pinned fixtures in
`research/pde_ledger_v3/scripts/test_S11c_b_cross_engine_comparator.py`. Inspect the working-tree diff
(`git diff research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py`). The fix repairs an asymmetric
time-order collapse (WL `D[…,{time,2}]`→`u_t_t`→was `u_t`; PY `u_tt`) in a canonicalizer SHARED by S11c-a and
S11c-b. A wrong fix silently corrupts cross-engine comparison of every inertia term.

## What to check — report a finding only if it catches a way the physics/comparison could be wrong

1. **Fix correctness + MANDATORY form ablation.** Trace `canon_jet_name` for `u_1_t_t`, `u_1_tt`, `u_1_t`,
   `u_1_t_t_d1`, `A_T_s11cb_1_t_t_d2_d3`, and a spatial-only name, before and after the fix; show the literal
   outputs. Then ABLATE: revert the time handling to the Boolean `has_time` form (in a /tmp copy — never the
   working tree) and re-run the pinned fixtures + confirm WL order-2 time jets re-collapse to order 1 (the
   fixtures must FAIL). A coefficient-style tweak is not enough — change the FORM (the counter) and show the
   literal diff. Confirm the fix makes WL `*_t_t*` and PY `*_tt*` canonicalize to the identical string.

2. **No spatial / order-≤1 regression.** Confirm every spatial-only jet (`w1_profile_d1`, `mu_R_bg_d1_d2`,
   `theta_d1`) and every order-1 time jet (`u_1_t`, `e_W_t`, `theta_t`) canonicalizes byte-identically to the
   pre-fix output. Enumerate the jet vocabulary and check for ANY name whose canonical form moved that should
   not have (e.g. a base name containing a literal `t`, or a `tt` substring that is not a time jet).

3. **S11c-a REGRESSION — the key check; do NOT trust the builder's self-report.** Re-run the S11c-a comparator
   (`python3 S11c_a_cross_engine_comparator.py` with its committed default transcripts) BEFORE and AFTER, and
   diff the two outputs yourself. ⛔ Do NOT `git stash` or mutate the working tree (another leg runs
   concurrently) — get the pre-fix comparator via
   `git show HEAD:research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py > /tmp/<you>_prefix_cmp.py`
   and run that copy against the same committed transcripts. The builder claims 4
   "serialization-only" differences that are canonically zero and an identical symbol multiset. Verify that
   INDEPENDENTLY: for each differing line, compute the actual residual (`simplify`/`expand`/`together`/`cancel`,
   and integral linearity if bound integrals appear) and confirm it is genuinely zero — a "canonically zero"
   claim that is actually nonzero, or any change beyond serialization order, is a blocking regression.

4. **The coupling asymmetry must SURFACE, not be massaged.** The fix is expected to EXPOSE second-time
   derivatives of the transverse trial in COUPLING_KERNEL (WL had `Derivative[0,1,1,2][transverseTrialPotential*]`
   collapsed to `_t`; now `_tt`). Confirm the regenerated reconcile (`~/.s11_build/S11c_b_reconcile_run.out`) now
   carries those `A_T_s11cb_*_tt_*` jets on the WL side, that PY has no `_tt` counterpart in its coupling
   (it zeros `u_tt`/`e_tt` in restriction, `S11c_b_brane_operator_sympy_audit.py:696`), and that the fix did NOT
   suppress, zero, or "reconcile away" this difference. (You are confirming the instrument exposes the
   disagreement; you are NOT judging whether the disagreement is a real finding — that is the later
   adjudication.)

5. **Script hygiene.** The regenerated reconcile still PRINTS operands/residuals and asserts nothing on a
   measured payload; the fixtures assert only on synthetic canonical strings.

## Method + constraints
- Derive/trace independently; SAVE every ablation script and its literal stdout to named absolute paths under
  your scratch dir and report those paths. A prose "I checked" is discarded.
- ⛔ Copy anything you ablate to /tmp and ablate the COPY; NEVER modify the working tree.
- ⭐ Run your ablations in the FOREGROUND; do NOT background them or spawn a monitor loop.
- This artifact is pure Python (reads committed `.out` transcripts); ⛔ do NOT spawn Mathematica/wolframscript.
- Physics filter: report a finding only if it catches a way the comparison could be wrong, not "it would be
  wrong on a different input."
- Report: numbered findings (blocking vs non-blocking), each with the file+line, the ablation command + literal
  output, and the concrete correction. If sound, say so and name what you verified with which command.
