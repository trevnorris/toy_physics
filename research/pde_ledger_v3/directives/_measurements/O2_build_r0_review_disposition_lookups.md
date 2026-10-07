# Measurements — O2 build r0 review dispositions (generated 2026-10-07 05:48)

Generator: `_scratch/s9b_build/gen/o2_build_r0_lookups.sh` (grep/sed lookups only; output lines cut at 400 characters; this file is written by it).

## Baseline
```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_build_review_baseline.sha256)
research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py: OK
research/pde_ledger_v3/scripts/O2_live_balance_sympy_ablation.py: OK
research/pde_ledger_v3/scripts/O2_exports.py: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl: OK
```

## F1 — free face label in the reduced load and face work (both engines)
```
$ grep -n "s = sp.Symbol('s_face')\|FaceLoadReduction_\|FaceWorkReduction" scripts/O2_live_balance_sympy_audit.py
158:    s = sp.Symbol('s_face')
171:    reduced_load = vector(action('FaceLoadReduction_' + str(a), op['J_map'],
200:    mechanical_power = action('FaceWorkReduction', op['J_map'], native_area,
```

```
$ grep -n 'qFace\[s\]\|orientation\[s\]\|THold\[s\]' mathematica/O2_live_balance_mathematica_audit.wl
56:holdOperand = operand[THold[s], "3.3,4,5", "7"];
81:faceHeight = qFace[s][x1,x2,x3,t];
84:faceNormal = orientation[s] Join[-faceSlope, {1}]/faceArea; (* SITE K7 *)
93:  "NativeFaceOrientation" -> (orientation[s]^2 == 1),
```

```
$ grep -c 'Sum(' scripts/O2_live_balance_sympy_audit.py; grep -n 'Sum\[' mathematica/O2_live_balance_mathematica_audit.wl
0
132:transport = Table[Sum[sectionDerivative[momentumCurrent[[a,i]],
135:  Sum[vPlane[[i]] D[value,x[[i]]],{i,3}]], profiles];
242:energyTransport = Sum[sectionDerivative[energyCurrent[[i]],energySection,x[[i]]],{i,3}]; (* SITE K10 *)
```

## F2 — PY face normal and area are independent OPEN actions; WL builds them from a face chart
```
$ grep -n 'NativeFaceNormal_\|NativeFaceAreaFactor\|bulk_load = ' scripts/O2_live_balance_sympy_audit.py
160:    native_normal = vector(action('NativeFaceNormal_' + str(a), face_context)
162:    native_area = action('NativeFaceAreaFactor', face_context)
164:    bulk_load = amplitude * native_normal
```

```
$ grep -n 'faceSlope = \|faceArea = \|faceNormal = ' mathematica/O2_live_balance_mathematica_audit.wl
82:faceSlope = D[faceHeight, #] & /@ x;
83:faceArea = Sqrt[1 + faceSlope.faceSlope];
84:faceNormal = orientation[s] Join[-faceSlope, {1}]/faceArea; (* SITE K7 *)
```

## F3 — PY harness difference(): matrix branch returns before the equality test
```
$ sed -n '60,78p' scripts/O2_live_balance_sympy_ablation.py
def difference(baseline, corrupted):
    """Corrupted minus baseline; structured fields retained even when zero.

    CAS tuples/matrices are compared recursively. Metadata has a structural
    edit pair if changed, zero if identical; no fictitious subtraction of text.
    """
    if isinstance(baseline, sp.MatrixBase) and isinstance(corrupted, sp.MatrixBase):
        if baseline.shape == corrupted.shape:
            return corrupted-baseline
    if isinstance(baseline, sp.Tuple) and isinstance(corrupted, sp.Tuple):
        if len(baseline) == len(corrupted):
            return sp.Tuple(*(difference(a,b) for a,b in zip(baseline,corrupted)))
    if baseline == corrupted:
        return sp.S.Zero
    if isinstance(baseline, sp.Equality) and isinstance(corrupted, sp.Equality):
        return sp.Tuple(corrupted.lhs-baseline.lhs, corrupted.rhs-baseline.rhs)
    if isinstance(baseline, sp.Expr) and isinstance(corrupted, sp.Expr):
        return corrupted-baseline
    return sp.Tuple(Str('STRUCTURAL_EDIT'), baseline, corrupted)
```

## F4 — PY K6c replaces the undifferentiated profile; the slope is built before the replacement
```
$ grep -n "'K6c'" -A1 scripts/O2_live_balance_sympy_ablation.py
41: 'K6c': ('    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6',
42-         "    momentum_profiles = sp.Tuple(vr, rho, sp.Symbol('xi_w_constant'))  # KNIFE_K6"),
```

```
$ grep -n 'tangent = embedding.jacobian\|graph_velocity = tangent\|momentum_profiles = \|momentum_velocity = ' scripts/O2_live_balance_sympy_audit.py
113:    tangent = embedding.jacobian(x)
121:    graph_velocity = tangent * v  # KNIFE_K2
128:    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6
130:    momentum_velocity = graph_velocity.xreplace(momentum_replacements).doit()
```

## F5 — PY name rho_br matches an upstream row that is a different object
```
$ grep -n "'rho_br'" scripts/O2_live_balance_sympy_audit.py | head -3
92:        'V_r', 'rho_br', 'mu_perp', 'xi_w', 'h', 'delta', 'j_n', 'f')}
232:    grades = text_tuple('rho_br', 'mu_perp', 'I_br_live', 'T_br_live', 'N_br_live',
```

```
$ grep -n -A7 "^    'rho_br': " scripts/S11c_b_exports.py
9274:    'rho_br':     {
9275-        'display': 'rho_br',
9276-        'value': _restore("Symbol('rho_br', positive=True)"),
9277-        'value_kind': 'COMPUTED_OBJECT',
9278-        'class': 'KNOB',
9279-        'step': 'S11c-b',
9280-        'description': 'uniform integrated brane density',
9281-        'f9_operands': _restore("Tuple(Symbol('rho_br', positive=True), Symbol('rho_br', positive=True))"),
```

## F6 — WL typed construction claims
```
$ grep -n 'StressOccurrences\|TransformationsUsedInMomentum\|no such transfer used' mathematica/O2_live_balance_mathematica_audit.wl
150:  "StressOccurrences" -> {fullStress}, "Origin" -> origin["3.2,4", "3,4"]|>];
216:  "TransformationsUsedInMomentum" -> {}, "Measure" -> CoordinateVolume[x],
217:  "Qualification" -> "Recorded relative O(epsilon) qualification on j_n when transferred to induced measure; no such transfer used",
```

## F7 — operand registers: R_br and M_perp
```
$ grep -n 'R_br\|M_perp' scripts/O2_live_balance_sympy_audit.py
96:             'A_rot_live R_ref_strain_live M_perp R_br E_h_live J_map '
```

```
$ grep -n '^    material = ' -A2 scripts/O2_live_balance_sympy_audit.py
123:    material = sp.Tuple(op['I_br_live'], op['P_br_cons'], op['N_br_live'],
124-                        op['A_rot_live'], op['J_map'], op['B_A13'],
125-                        op['material_action_compatibility'])
```

```
$ grep -n 'materialOperands = ' -A2 mathematica/O2_live_balance_mathematica_audit.wl
105:materialOperands = {inertia,consMomentum,fullStress,consStress,normalResponse,
106-  rotation,branch,mapOperand,densityResponse,stiffnessResponse,
107-  sourceInventory,boundaryInventory};
```

```
$ grep -n 'ℛ_br. (O7)' directives/O2_SHARED_PHYSICS.md
138:| `ℛ_br` (O7); **OPEN**, C §6 | General brane-density response with live bulk, flow, embedding, thickness/projection dependence and gradients. It is input to the density in mass/momentum/energy accounting; it is not inferred from the bulk EOS or slab factorization. | S8, substrate reduction and S22 completion. |
```

## Note — tag sets (directive item 4)
```
$ grep -o 'emit\["[A-Z_]*"' mathematica/O2_live_balance_mathematica_audit.wl | sort -u | wc -l
12
```

```
$ awk 'NR>=230 && NR<=341' scripts/O2_live_balance_sympy_audit.py | grep -c "^        '[A-Z_]*'"
25
```

