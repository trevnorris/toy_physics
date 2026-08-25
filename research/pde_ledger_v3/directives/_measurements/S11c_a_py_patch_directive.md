# Measurements twin — S11c_a_py_patch_directive.md

Every factual claim the directive makes about the SymPy engine's CURRENT state carries the command that
produced it. ⚠ rule-5: no output-expression value appears here; claims are structural.

CONVENTIONS: run from repo root `/var/projects/toy_physics` with
`PY=research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`,
`SPEC=research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`. Engine at HEAD = committed `9b6438fa`.

## Claim (Verdict C): virtual-work is DIAGONAL; one `dof` drives physical shape AND virtual displacement
CMD: `sed -n '854,862p' "$PY"` -> `virtual_work_cases(branch, dof, representative, ...)` — single `dof`.
CMD: `sed -n '740,756p' "$PY"` -> in `build_face_source`, `dof` sets `zeta`→`dh`→`normal_exact`/`measure_exact`
(physical shape) AND `virtual_zeta`/`virtual_vertical`→`virtual_displacement` (the virtual `δ_v`).
CMD: `sed -n '919,924p' "$PY"` -> emit loop keys `(branch, dof, representative)` (diagonal).
CMD: `grep -n 'emit_primary("VIRTUAL_WORK_SHAPE_DERIV"' "$PY"` -> line 1269, `export=True` path.
Spec: `sed -n '417,420p' "$SPEC"` (§4 T-d: "Which virtual-displacement pairings occur is part of the computation").

## Claim (Verdict BG): BACKGROUND_STATE carries the zeros but NO boundary loads; loads live in the premise
CMD: `sed -n '1320,1332p' "$PY"` -> `state` per (branch,representative) = profile_context + θ⁰=0 + V_0/J_0/A_0=0;
premise Tuple carries `(support_body, support_face_plus, support_face_minus)`.
CMD: `sed -n '174,176p' "$PY"` -> `support_body=f_hold_0`, `support_face_plus=t_hold_plus_0`,
`support_face_minus=t_hold_minus_0` (supplied PREMISE symbols).
CMD: `grep -n 'emit("BACKGROUND_STATE"' "$PY"` -> line 1335, plain `emit` (NOT export=True).
CMD: `sed -n '1649p' "$PY"` -> the `BACKGROUND_STATE` dimension tuple (must grow with the added fields).
Spec: `sed -n '246,272p' "$SPEC"` (§2d:251 `𝔅⁰ ≡ {…, boundary loads}`, `𝒮_hold⁰ ≡ {f_hold⁰, t_hold,s⁰}`).

## Claim (density map): BACKGROUND_DENSITY_MAP carries a redundant branch axis; sigma_e is branch-independent
CMD: `sed -n '897,905p' "$PY"` -> `build_background_density_raw` loops branch×representative, keys
`(branch, representative)`; `sigma_e = rho4_value * W_bg` with `rho4_value,_ = density_pair(representative, W_bg)`
(no branch dependence).
CMD: `grep -n 'emit_primary("BACKGROUND_DENSITY_MAP"' "$PY"` -> line 1232, `export=True`.
CMD: `sed -n '1611,1616p' "$PY"` -> `independent_uniform_reference` keys `(branch, representative)` (must match).
Spec: `sed -n '228,231p' "$SPEC"` (§2b: "Emit the two computed maps" = per representative, pre-anchoring y).

## Export feed (why C and the density map mutate exports; BG does not)
CMD: `sed -n '333,343p' "$PY"` -> `emit` appends a `CandidateRow` to `export_candidates` only when
`export_key is not None` (i.e. `emit_primary`, export=True). CMD: `grep -n 'export_candidates' "$PY"`.
⇒ VIRTUAL_WORK_SHAPE_DERIV + BACKGROUND_DENSITY_MAP are export=True (payload changes propagate to
`S11c_a_exports.py`); BACKGROUND_STATE is plain emit (no export change).

## Discipline
Directive relays the step-1 verdicts (C full grid, BG loads, density-map branch drop) as spec-grounded case-
structure changes, ⛔ never with an expected value and ⛔ never leaking any sibling-engine finding (e.g. it does
NOT state whether the off-diagonal pairings are redundant — that is the comparator's measurement, not a target).
