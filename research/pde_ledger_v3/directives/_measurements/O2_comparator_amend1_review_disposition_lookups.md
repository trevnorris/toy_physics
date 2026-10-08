# Measurements — O2 comparator directive amendment 1 review dispositions (generated 2026-10-07 18:47)

Generator: `_scratch/s9b_build/gen/o2_comparator_amend1_lookups.sh` (sed/grep/sha256sum/diff only; lines cut at 300 characters).

```
$ sha256sum _scratch/s9b_build/O2_comparator_build_directive_amend1_reviewed_v0.md research/pde_ledger_v3/directives/O2_comparator_build_directive.md
4db4ce5794457d39b7f973aded8a7eadebe90f6c95b6bdf9cc2b680b883e69ae  _scratch/s9b_build/O2_comparator_build_directive_amend1_reviewed_v0.md
c3ebe39eaca4295314efba381d720fb0f239379bc6e8633a2c640b234381e93e  research/pde_ledger_v3/directives/O2_comparator_build_directive.md
```

## A1 — the reviewed limit sentence; how each engine packages an OPEN action's arguments
```
$ grep -n 'occur in more than' _scratch/s9b_build/O2_comparator_build_directive_amend1_reviewed_v0.md
76:   paired OPEN occurrence and each engine, the comparator computes and prints which live objects occur in more than
```

```
$ sed -n '85p;249,252p' research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py
def action(role, *operands):
    carried_w = action('CarriedMomentumW', op['Pi_n'], jn, op['J_map'],
                       op['N_br_live'], graph_velocity, geom, state, source,
                       boundary, Str('outward_native_relative_mass_current'),
                       Str('premise_3_local_material_velocity_at_each_transfer'))  # KNIFE_K13
```

```
$ sed -n '23p;130,131p;202,205p' research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl
response[role_, operands_, state_] := OpenAction[role, operands, state];
section = <|"Profiles" -> profiles, "Velocity" -> vMaterial,
  "Metric" -> metric, "OtherDependence" -> UnrestrictedSection[
carriedBulk = response[OutwardCarriedBulkMomentum,
  {exchangeOperand,mapOperand,normalResponse,branch,sourceInventory,boundaryInventory,
   NativeRelativeMassCurrent,Premise3LocalMaterialVelocity,
   UnspecifiedFaceToMaterialVelocity},section]; (* SITE K13 *)
```

## A2 — the reviewed control list: subtraction, operand-before-guard, representation pair, 'or' bullets, verdict token
```
$ grep -n 'subtraction\|before the residual\|VERDICT' _scratch/s9b_build/O2_comparator_build_directive_amend1_reviewed_v0.md
55:   Wolfram operand and their residual. Use exact symbolic subtraction where feasible. Otherwise use numeric
81:   measured streams is expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`VERDICT`/`STATUS` token. It exits 0
```

```
$ grep -n ' or binder structure\| label, live object or orientation' _scratch/s9b_build/O2_comparator_build_directive_amend1_reviewed_v0.md
98:   - a changed derivative order, derivative evaluation point or binder structure;
99:   - a changed OPEN head, named operand, label, live object or orientation;
```

## Fold applied
```
$ diff _scratch/s9b_build/O2_comparator_build_directive_amend1_reviewed_v0.md research/pde_ledger_v3/directives/O2_comparator_build_directive.md | grep -c '^[<>]'
83
```

