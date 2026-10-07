# Measurements — O2 build r2 review dispositions (generated 2026-10-07 11:43)

Generator: `_scratch/s9b_build/gen/o2_build_r2_lookups.sh` (grep/sed/git-show lookups only; output lines cut at 400 characters; this file is written by it).

## Baseline
```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_build_review_baseline_r2.sha256)
research/pde_ledger_v3/scripts/O2_exports.py: OK
research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py: OK
research/pde_ledger_v3/scripts/O2_live_balance_sympy_ablation.py: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl: OK
```

## C1 — WL bulk-normal load amplitude
```
$ grep -n 'bulkAmplitude\|bulkTraction = ' mathematica/O2_live_balance_mathematica_audit.wl
57:bulkAmplitude = operand[TBulkNormalLive[s], "3.3,5", "7"];
104:bulkTraction = bulkAmplitude faceNormal;
```

```
$ grep -n "T_bulk" scripts/O2_live_balance_sympy_audit.py
229:    amplitude = action('T_bulk_n_s_live', native_point_context)
```

```
$ sed -n '155p' directives/O2_SHARED_PHYSICS.md
| `t_bulk,s^live = 𝒯_bulk,n,s^live n̂_s`; directional restriction from **adopted premise 4**, C §7 | Bulk mechanical traction on its native face. `𝒯_bulk,n,s^live` is a general signed **OPEN** normal-load amplitude; geometry and projection are live. No pressure/affinity/DC response law or sign of the amplitude is given. |
```

```
$ sed -n '34,36p' directives/O2_SHARED_PHYSICS.md
Keep `V`, `ρ_br`, `μ_⊥`, `ξ_w`, `h`, `δ`, `j_n`, the bulk state, stress and inertial responses,
reference/strain state, loading and boundary data live, including their spatial derivatives and
material-history dependence. Eulerian steadiness does not impose material constancy. No particular
```

## C2 / G2 — WL height-chart face restriction; spec §10 rule; r1 disposition wording
```
$ grep -n 'height-chart\|Coverage\|faceHeight = ' mathematica/O2_live_balance_mathematica_audit.wl
70:(* The supplied brane graph and a CONDITIONAL native height-chart ansatz.
76:  "Coverage" -> "Restricted to that height-chart domain; faces without such a chart are outside this representation",
98:faceHeight = qFace[s][x1,x2,x3,t];
```

```
$ sed -n '623,626p' directives/O2_SHARED_PHYSICS.md
compatibility remain at their named later steps. The named untruncated conditional balance is the
finite construction target. If a more explicit action requires a physical restriction not present
here, retain the named OPEN action where possible; otherwise report the missing restriction as a
question for the user, without choosing it. A second method failure, new unnamed sub-problem or
```

```
$ grep -n '^| C2 ' directives/_measurements/O2_build_r1_review_disposition.md
30:| C2 | WL represents each native face as a graph `qFace[s][x,t]` over the far-field coordinates. That excludes faces that are not single-valued over x, and the restriction is not stated. PY uses a general immersion. (Claude) | **ACCEPT, WL.** Lookup: WL lines 81–84 and 91. The spec leaves native geometry and measures to `𝒥_map` (§5). | Every restriction an engine's face representation pla
```

## G1 — WL balance relations
```
$ grep -n '"Relation"' mathematica/O2_live_balance_mathematica_audit.wl
236:  "Relation" -> Thread[holdBalance == ConstantArray[0,4]],
297:  "Object" -> energyBalance,"Relation" -> (energyBalance == 0),
```

```
$ grep -n 'holdBalance\|energyBalance' scripts/O2_live_balance_sympy_audit.py
```

```
$ sed -n '81,83p' directives/O2_SHARED_PHYSICS.md
The target is conditional on premises 1–4 below. It is a named balance with unresolved responses,
not a derived substrate law or a solution for the profiles. This specification supplies no assembled
momentum or energy balance, expected residual, sign, cancellation, profile or compatibility outcome.
```

```
$ sed -n '272,274p' directives/O2_SHARED_PHYSICS.md
**`ℬ_E^steady` is OPEN**, with the accounting requirement fixed by **adopted premise 1** and its S21
label in §2 (C §8). It denotes the steady energy relation to be constructed, not a precomputed
residual. Its required inputs are:
```

```
$ sed -n '620,623p' directives/O2_SHARED_PHYSICS.md

Constitutive elimination, additional premise selection, S11c repair/composition, S14a bridge work,
source/response `GM` matching, physical holder selection/solve, S21 integration and optical
compatibility remain at their named later steps. The named untruncated conditional balance is the
```

## G3 — WL K10 site versus the directive
```
$ grep -n 'K10' directives/O2_build_directive.md
62:   - **K10, energy transport:** the live energy transport `𝒥_E^live` is removed from the energy accounting.
```

```
$ grep -n 'K10' mathematica/O2_live_balance_mathematica_ablation.wl
38: {"K10","energyTransport = Sum[sectionDerivative[energyCurrent[[i]],energySection,x[[i]]],{i,3}]; (* SITE K10 *)",
39:   "energyTransport = 0; (* SITE K10 *)"},
```

```
$ grep -n 'energyCurrentOperand\|energyInputs = ' mathematica/O2_live_balance_mathematica_audit.wl
61:energyCurrentOperand = operand[JELive, "6", "8"];
265:energyInputs = {energyOperand,energyCurrentOperand,relaxationPower,
```

## Provenance — present before the round-2 repair? (r0 = 93ed661d, r1 = ba1c65cb)
```
$ git show 93ed661d:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl | grep -c 'bulkAmplitude = operand\[TBulkNormalLive\[s\]'
1
```

```
$ git show 93ed661d:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl | grep -c '"Relation" -> '
2
```

```
$ git show 93ed661d:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl | grep -c 'energyInputs = {energyOperand,energyCurrentOperand'
1
```

```
$ git show 93ed661d:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl | grep -n 'K10' | cut -c1-200
38: {"K10","energyTransport = Sum[sectionDerivative[energyCurrent[[i]],energySection,x[[i]]],{i,3}]; (* SITE K10 *)",
39:   "energyTransport = 0; (* SITE K10 *)"},
```

```
$ git show ba1c65cb:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl | grep -c 'bulkAmplitude = operand\[TBulkNormalLive\[s\]'
1
```

```
$ git show ba1c65cb:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl | grep -c '"Relation" -> '
2
```

```
$ git show ba1c65cb:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl | grep -c 'energyInputs = {energyOperand,energyCurrentOperand'
1
```

```
$ git show ba1c65cb:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl | grep -n 'K10' | cut -c1-200
38: {"K10","energyTransport = Sum[sectionDerivative[energyCurrent[[i]],energySection,x[[i]]],{i,3}]; (* SITE K10 *)",
39:   "energyTransport = 0; (* SITE K10 *)"},
```

```
$ git show ba1c65cb:research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl | grep -c 'height-chart'
0
```

