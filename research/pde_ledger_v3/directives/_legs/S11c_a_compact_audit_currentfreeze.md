# Compact-prep doc audit — S11c-a T7 current-freezing fix round

You are auditing whether this session's DOC UPDATES accurately and honestly reflect the actual repository
state, before a context compaction. Report ONLY discrepancies (a doc claim that contradicts the real
git/engine state, a stale/misleading statement, a leaked expected value, an inaccurate command/sha). If a
doc is accurate on a point, say so briefly. Do NOT fix anything; just report. Run commands from repo root
`/var/projects/toy_physics`.

## Ground truth to check against (run these yourself)
- `git log --oneline -6` — confirm the current-freezing fix is committed as `49b5c525` on `ledger-v3-rebuild`,
  atop `ff634ef6`, `6fae82b8`, etc.
- `git show 49b5c525 --stat` and `git show 49b5c525 -- research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
  — confirm the fix changed `projection_terms` + `uniform_projection_reference` so the normal current is a
  `w`-dependent field (`delta_j_bulk_4`/`normal_current` as a `Function(w)`), NOT the bare constant `j_bulk[3]`,
  and that it did NOT add a face axis or a second `∫Ω·∂_w δj` channel.
- In the current engine, inspect `projection_terms` (~lines 1150–1170): the `WINDOW_NORMAL_CURRENT` origin
  should use the `w`-dependent current, and `current_divergence` should still be in-plane-only (`range(3)`).

## Docs to audit (read them in the WORKING TREE — updated, some uncommitted)
1. `STATUS.md` — the top "CURRENT FRONT" section + the new `#3 CURRENT-FREEZING FINDING` block + the
   `COMPARATOR DEFERRED` block + the `NET`/`NEXT` bullets. Verify: (a) it states THREE §3c-class findings all
   fixed, with the correct shas (`c36beac4`, `6fae82b8`, `49b5c525`); (b) it says the projection integrand was
   NOT mechanical — it was finding #3 (PY froze the perturbation current; `WINDOW_NORMAL_CURRENT` was
   identically 0); (c) it honestly RETRACTS the §5c uniform-limit "corroboration" as a non-simplified-zero
   false alarm (it must NOT claim the uniform-limit confirmed #3); (d) it says the ULTIMATE cross-engine
   projection AGREEMENT is NOT yet computed (deferred to the comparator). Flag any leftover "mechanical, not a
   finding" / "two findings" / "BOTH §3c findings" language that is now wrong.
2. `research/pde_ledger_v3/directives/S11c_a_comparator_reemit_plan.md` — the Step 5 `NET`/`NEXT`, the Step-5
   PLAN (should say DEFERRED with the ~14-defect rev-1 list incl. the 2 SMUGGLING folds), and the `## State`
   line (shas; PY tag stream now `~/.s11_build/S11c_a_py_fixed_run2.out`). Same checks; confirm no stale
   "two findings / mechanical bridge / NEXT = projection IBP bridge" survives.
3. The fix directive `research/pde_ledger_v3/directives/S11c_a_py_projection_current_fix_directive.md` (rev-2,
   committed) — confirm it names the OBJECT (∂_wδj_w must enter), does NOT prescribe the face-trace recipe
   (`affine_bulk_perturbation`/`dw_delta_j_bulk[face]`), states the no-double-count constraint, and leaks no
   expected residual value (rule 5).

## Specifically flag
- Any sha that does not match `git log` / `git show`.
- Any doc still calling the projection integrand a "mechanical bridge" or "NOT a finding".
- Any place that claims the uniform-limit CORROBORATED #3 (it was retracted — the residual is a
  non-simplified zero under integral linearity; #3 stands on the CAS consult + adjudication + the pre-fix
  `WINDOW_NORMAL_CURRENT ≡ 0` build finding).
- Any claim that the cross-engine projection is confirmed to AGREE post-fix (it is NOT yet computed — that is
  the deferred comparator's job).
- Any expected-value leak in the fix directive (rule 5).

Return a short list of discrepancies (or "in the clear" if none), each with the file, the quoted text, and
the ground-truth it contradicts.
