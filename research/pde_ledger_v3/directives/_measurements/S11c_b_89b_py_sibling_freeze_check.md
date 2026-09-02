# S11c-b #89b — PY sibling freeze-check (rule 16): does the SymPy engine have the WL re-freeze?

After the WL build legs caught (and the repair fixed) a re-freeze in the WL EMITTED operator (it reduced
`operatorLive` while the outer `Inactive[Div]` was held, dropping the mixed/higher U-row Hessian jets), we
must check whether the sibling PY engine `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
has the analog — otherwise the pre-#89b cross-engine operator agreement could have been two engines agreeing on
the same FROZEN object, and fixed-WL would now disagree with a still-frozen PY. Investigated by a fresh agent
(source reading + an isolated reduced computation; ⛔ the full PY run OOMs ~3.5 h, not run).

## VERDICT: PY IS CORRECT (activate-then-reduce). No freeze defect. No PY change needed for the freeze.

### Order of operations (emitted operator) — reduction is LAST
- `task_slab_operator` (`:3120-3124`): `operator = build_operator(...)` THEN `retained_grade(operator)`. The
  profile reduction (`retained_grade` `:1005-1022` → `first_shape_series` `:873-916`, whose DEFAULT path is
  `first_shape_series_fast` [first-order Taylor, `:888-913`] with `sp.series(...,eta_bg,0,2).removeO()` as the
  reference/fallback [`:873-886`], + `PROFILE_GRADE_SUBS`) runs only AFTER the operator is fully built.
- The build keeps the background LIVE: `operator_from_density` (`:1968-2062`) does the EL variation
  (`sp.diff(density,u[a])`) and the in-plane divergence via `total_derivative(..., background_depth=STRONG_ROW_JET_DEPTH)`
  (`:1988-2006`); `total_derivative` (`:696-737`) maps `W_bg→grad_W→background_jet_expression(source,2,…)` (the
  Hessian) while the coefficient is still symbolic. No `subs`/`series` on the background before the variation.
- The frozen derivative (`frozen_derivative` `:740-745`, background frozen at 1st order) is used ONLY in the
  diagnostic controls `committed_*` (`:2089-2175`) inside the HESSIAN_FREEZE task (`:3540-3599`) — it never
  feeds the emitted primary.

### Computational confirmation (isolated, ~24 s; `…/scratchpad/py_freeze_probe2.py`)
LAB_HELD / RHO4_CONSTANT / EULERIAN, live strong rows (emitted path, depth 2) vs frozen `committed_strong_rows`,
both after `retained_grade`. Literal stdout:
```
LIVE   2nd-order bg jets (HESSIAN): [m1_profile_d1d1 … w1_profile_d3d3]   (full symmetric Hessian)
FROZEN 2nd-order bg jets (HESSIAN): []
LIVE_vs_FROZEN residual nonzero-scalar count: 5   (= exactly the dropped Hessian terms; all 5 strong rows)
```
The Hessian atoms are σ_W¹ so they survive `retained_grade`'s first-shape truncation. ⇒ the emitted PY strong
rows retain the Hessian; only the diagnostic frozen rows drop it.

### #89 coverage
`git show f655ea65` (the #89 checkpoint = the current committed PY engine) INTRODUCED `total_derivative` with the
`background_depth` jet tower and threaded it through the entire live path (`density_pair`, `material_pullback`,
`construct_energy`, `live_basis_substitution`, `local_thickness_map`, `operator_from_density`, curl/divergence),
demoting `frozen_derivative` to the `committed_*` controls. So PY #89's "strong-rows un-freeze" genuinely covers
the emitted operator's differentiation/reduction ORDER, not merely the basis count.

## ⚠ RECONCILIATION ITEM for integration (a FLAG, not a freeze defect)

PY caps strong rows at **`STRONG_ROW_JET_DEPTH = 2`** — a single in-plane divergence yields background jets up to
2nd order (Hessian) and **NO 3rd-order** background jets in the strong rows (3rd order lives in the coupling
cascade, `COUPLING_JET_DEPTH = 3`, `:60`). The WL #89b fix now emits **3rd-order** background jets in the strong
U-momentum rows (the re-review's decisive evidence: "order-3 restored", higher-jet atoms 16 vs 6). ⇒ WL strong
rows will carry 3rd-order terms PY does not — a **grade/jet-depth convention mismatch** that `scripts/S11c_b_row_residual.py`
will surface as a strong-row disagreement. This is a SPEC question (which strong-row background-jet depth is
correct — a single divergence naturally raises to 3rd, PY truncates to 2), SEPARATE from the freeze (which PY
passes). ⇒ adjudicate at integration (likely a spec-confirm, like the other both-engine convention questions)
BEFORE reading the row_residual as a physics disagreement. Do NOT assume WL over-emits or PY under-emits — compute
which depth the retained grade requires.

Evidence scripts (scratchpad, ephemeral): `…/scratchpad/py_freeze_probe2.py` (+ stdout above).
