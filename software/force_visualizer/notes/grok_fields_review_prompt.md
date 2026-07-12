# Grok review — force_visualizer FIELD visualizations (READ-ONLY)

You are an adversarial expert reviewer. READ-ONLY: do not modify files. COMPUTE-verify with SymPy/numpy where useful.

## Context
`software/force_visualizer/` is a 2D phenomenology visualizer for a toy 4D-superfluid analog model (four sectors: gravity, light, charge, magnetism). The physics core (`physics/`) already passed a full fidelity review (FAITHFUL). A new increment added **per-sector spatial FIELD visualizations** (spec §10 in `software/force_visualizer/notes/build_spec.md`), computed from the existing core. Your job: check that the field REPRESENTATIONS are faithful to each sector's mechanism and honestly labeled, and that the new field checks are able-to-fail. **No force law should have changed.**

## The model's sector distinctions you must enforce (do not let the viz conflate them)
- **Gravity** = the medium's one-way inflow/drain (a real medium flow toward masses; streamlines terminate at the masses, no circulating-back-out).
- **Charge** = the throats' 4D-body interaction (an electric FIELD), explicitly **NOT a medium flow** — a flow carrying energy would be a mass current = gravity. So charge must be shown as field lines/vectors, NOT flowing tracers.
- **Magnetism** = the moving throat's 4D-body swirl, felt on the brane as a **circulating** magnetic field. The literal swirl is in the throat's 4D body (off-brane).
- **Light** = the brane transverse shear wave (already the field).

## Files
- New/updated field overlays: `software/force_visualizer/scenes/*.py`, `scenes/data.py`
- Field-direction checks: `software/force_visualizer/report.py` + `output/verification_report.txt` (the three new checks: gravity drain-field toward mass; charge E away-from-+/toward-−; magnetism B circulation)
- Tests: `software/force_visualizer/tests/*.py`
- Core (should be UNCHANGED): `software/force_visualizer/physics/*.py`

## Check and COMPUTE-verify
1. **Gravity field** — is it the inward drain/inflow field (∝ 1/r² toward masses, one-way, terminating at sinks — no closed loops)? Report shows field at (2,1) = (−0.179,−0.089), toward-mass alignment 1.00, |field|=0.200=1/r² — verify. Label must read as the medium inflow (gravity = the drain), not "force field."
2. **Charge field** — Coulomb field lines/vectors with away-from-+/toward-− geometry (report: +source E=+x away, −source E=−x toward). Confirm it is rendered as a FIELD, NOT as flowing/streaming tracers, and labeled as such (a field, not a medium flow).
3. **Magnetism field** — circulating B around a current (report: forward-current B=+y at +x, reverse=−y). Confirm the circulation (curl) structure and the sign-flip with current direction, labeled as the brane manifestation of the throat's 4D swirl.
4. **Labels/honesty** — no label conflates one sector's mechanism with another's; the charge overlay must not imply a medium flow; gravity's must not imply a static field. The "phenomenology visualizer, not emergence" honesty holds.
5. **The 3 new field-direction tests are genuinely able-to-fail** (not tautological) — a wrong field direction/circulation must fail them.
6. **No physics regression** — `physics/*.py` force laws unchanged; `physics/` still imports no rendering module.

## Output
Per new field (gravity/charge/magnetism): FAITHFUL or DISCREPANCY (with file:line + what's wrong). Then: any label that conflates/misleads; any new field test that cannot fail; confirm no physics-form change. Overall verdict: CLEAN / ISSUES (prioritized). Be specific.
