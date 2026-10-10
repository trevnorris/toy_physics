# Measurements — S9b directive amendment 6, review round 0 (generated 2026-10-10 10:39)

Generator: `_scratch/s9b_build/gen/s9b_amend6_r0_lookups.sh` (sha256sum/grep/sed/cat/wc only). The reviewed
directive is frozen as `_scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md`.

## The reviewed version
```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md
46b942e4fac1ee8343b411ce44ee587f1a5ff476ba5f173cee93d990fde235b6  _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md
```

```
$ head -1 _scratch/s9b_build/s9b_repair_build_directive_amend6_review_baseline_r0.sha256
46b942e4fac1ee8343b411ce44ee587f1a5ff476ba5f173cee93d990fde235b6  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend6_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend6_review_r0_grok.txt
_scratch/s9b_build/s9b_amend6_review_r0_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_amend6_review_r0_grok.txt:3:Verdict: NEEDS REVISION
```

## Grok finding 1: item 14's bullet does not say which part of a logical payload is the value
```
$ sed -n 223,225p _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md
    - For a payload that carries a domain or conditions, the difference shows separately whether the value changed
      and whether the domain changed. Each part's difference is computed by the engine, so an unchanged part prints
      as zero.
```

### A Wolfram payload whose computed content is logical sits beside a Domain key
```
$ sed -n 166,171p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
emit["A_NONRECIPROCAL_PATH_DEPENDENCE", <|"OneForm" -> oddOneForm,
  "ExteriorDerivative" -> Simplify[oddExterior],
  "Closedness" -> Simplify[And @@ Thread[Flatten[oddExterior] == 0]],
  "Domain" -> (r > 0 && physicalTraversal && physicalPropagation && speedChart),
  "RadialPrimitive" -> Inactive[Integrate][oddOneForm[[1]], {r,rAnchor,rr}]|>];
(* COMPACT_K11_BOUNDARY *)
```

### The SymPy harness differences a bare boolean as one Xor
```
$ sed -n 213,214p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    if isinstance(a,sp.logic.boolalg.Boolean) and isinstance(b,sp.logic.boolalg.Boolean):
        return sp.Xor(a,b)
```

### Codex's probe paired a relation with a domain; it did not test a bare logical value
```
$ grep -n 'separated_boolean_value' _scratch/s9b_build/s9b_amend6_review_r0_codex_evidence/payload_review.wl
37:emit["separated_boolean_value",splitDifference[ConditionalExpression[x==0,e],ConditionalExpression[x==0,d]]];
```

## Grok finding 2: the Part 2 line
```
$ sed -n 309,310p _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md
- **Repair round 3 (amendment 6).** Items 5 and 6 hold for every `γ` object, and for its restriction and forward
  copies, on every stratum inside the branch domain, including `GM = 0`.
```

### 'the branch domain' occurs once, undefined
```
$ grep -n 'branch domain' _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md
310:  copies, on every stratum inside the branch domain, including `GM = 0`.
```

### Item 6's stratum sentence is about the γ difference only
```
$ sed -n 135,136p _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r0.md
   - The `γ` difference is reduced by the engine's CAS, so its printed form shows its value on each stratum. It is
     not printed as an uncombined difference of two objects.
```

### C1's required outcome, as the r2 disposition states it
```
$ grep -o 'Inside the branch domain, each .γ. object prints[^|]*' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r2_review_disposition.md
Inside the branch domain, each `γ` object prints, on every stratum including `GM = 0`, what its relation determines there, computed by the engine. `NOT_ESTABLISHED` appears only outside the branch-existence conditions. 
```

### The spec's branch-existence paragraph states two conditions
```
$ sed -n 193,198p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
along the ray and the condition under which a ray can traverse the full flyby path and both legs of the round
trip in the required directions. Outside either condition, report the branch type (growing, decaying, absent,
or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
them. State the first condition through `c_γ²` without presupposing its sign, so that each of these branch types
can follow from it.
```

## Process: Codex's orphaned probe
```
$ cat _scratch/s11c/s9b-amend6-review-r0/codex-payload-wl-01/RECONCILIATION.md
# Reconciliation of a stale reservation (orchestrator, 2026-10-10 10:38 MDT)

Unit `s11c-guard-f4020c055d31.service` ran the Codex amendment-6 review leg's auxiliary Wolfram probe. The Codex
session exited at 10:25 with its report, which says the verdict does not rely on this probe. Its guard launcher
(owner pid 426824) died with it. The probe's last output was at 10:20; at 10:37 it held ~200 MB and one kernel seat.

The orchestrator stopped the unit, then reconciled per the guard docstring (L19–22):
- terminal systemd state: `LoadState=not-found ActiveState=inactive SubState=dead`;
- no queued job for the unit (`systemctl --user list-jobs`: 0 lines);
- owner pid 426824 dead (`kill -0` fails);
- the old registry is kept as `_scratch/s11c/resource-pool/reservations.json.stale-20261010T163825Z`;
- only this unit's entry was removed, by atomic rewrite. The registry is now `{}`.
```

```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

