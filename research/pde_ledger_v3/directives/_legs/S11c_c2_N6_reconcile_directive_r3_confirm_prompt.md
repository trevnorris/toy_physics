# Round-3 SCOPED confirmation — the S11c-c2 N6 RECONCILE directive, after folding the 2 round-2 findings

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_directive.md` (v3). Rounds 1–2
(2 legs each) cleared the whole directive EXCEPT two items, now folded. This is a **tight confirmation**, ⛔ not a
fresh full review: verify the two folds landed correctly + faithfully, and scan ONLY for a regression they introduced.
⛔ No CAS, ⛔ no build, ⛔ no fictional-script ablation. Working dir `/var/projects/toy_physics`.

## The two folds to confirm (verify against the sources, not just that text is present)
1. **Carrier bite** (directive "Controls", the "Carrier bite:" bullet). Confirm it now corrupts the **material**
   covector/normal-map feeding `C_M` (`material_inverse_transpose`/material normal, route-2
   `_measurements/S11c_c2_N6_route2_spec_astra.md:113`) with `ms` held uncorrupted ⇒ moves
   `CARRIER_BRIDGE_RESIDUAL`(`C_E−C_M`)/`CARRIER_CHANNEL`, leaves `C_E` (imported `e_coeff` diagnostic:797), the
   Eulerian operand, and `SOURCE_BRIDGE_RESIDUAL` unchanged. Confirm the Eulerian-factory tilt (diagnostic:808,
   `tilt_coeff`) is correctly demoted to the SEPARATE reconstruction/independence probe, gated on unmodified-factory
   = imported `C_E` (baseline nonzero = reconstruction drift; route-2 :195,:209) — ⛔ not the carrier-bridge bite.
   Is the corruption target now one that actually moves `C_E−C_M`, one-sided, and non-A−A?
2. **SPLIT_CHECK verify line** (directive "Deliverable + verification"). Confirm it now says an independent residual
   node with **samplewise-zero shared-PIT numerators** (print node + per-sample numerators), ⛔ not "structurally
   zero", ⛔ not `number(0)`. Confirm no other "structurally zero" mention survives anywhere in the directive.

## Regression scan (only report a NEW defect the 2 folds introduced)
Did fixing the carrier bite or the verify line contradict any previously-cleared content — the three-way split, the
`es`/`ms` source pinning + in-process `pit()`, the jet bridge, the frozen relations, the coefficient-level carrier
bridge, the SPLIT_CHECK contract at the affine-split section, or the script clauses/corollaries?

## Sources
`_measurements/S11c_c2_N6_route2_spec_astra.md` (:113, :195, :209); `scripts/S11c_c2_N6_diagnostic_sympy.py`
(:797, :808, :809, :832-836, :137, :163).

## Output
End with **DIRECTIVE SOUND — CLEAR TO BUILD** or the exact remaining fold (with `file:line` + minimal fix). Brief.
