# Measurements — O2 build r1 review dispositions (generated 2026-10-07 08:21)

Generator: `_scratch/s9b_build/gen/o2_build_r1_lookups.sh` (grep/sed lookups only; output lines cut at 400 characters; this file is written by it).

## Baseline
```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_build_review_baseline_r1.sha256)
research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py: OK
research/pde_ledger_v3/scripts/O2_live_balance_sympy_ablation.py: OK
research/pde_ledger_v3/scripts/O2_exports.py: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl: OK
```

## C2 — WL native face as a graph chart over x
```
$ grep -n 'faceHeight = \|faceSlope = \|faceArea = \|faceNormal = \|NativeFaceMeasure' mathematica/O2_live_balance_mathematica_audit.wl
81:faceHeight = qFace[s][x1,x2,x3,t];
82:faceSlope = D[faceHeight, #] & /@ x;
83:faceArea = Sqrt[1 + faceSlope.faceSlope];
84:faceNormal = orientation[s] Join[-faceSlope, {1}]/faceArea; (* SITE K7 *)
92:  "DensityMeasure" -> CoordinateVolume[x], "NativeFaceMeasure" -> faceArea CoordinateVolume[x],
```

## C3 — PY material pairing
```
$ grep -n 'material_pairing = \|MaterialStressNormalWork' scripts/O2_live_balance_sympy_audit.py
258:    material_pairing = internal.dot(graph_velocity)
262:    material_work = action('MaterialStressNormalWork', material_pairing,
```

```
$ grep -n 'MaterialStressNormalRotationalWork' -A3 mathematica/O2_live_balance_mathematica_audit.wl
240:stressWork = response[MaterialStressNormalRotationalWork,
241-  {stressInputs,inertia,normalResponse,rotation},
242-  <|"ForceAction" -> internalForce,"GraphVelocity" -> vMaterial,
243-    "GeneralizedRates" -> operand[UnspecifiedRotationalNormalRates,"6","3,8"],
```

## C4 — PY load and work reductions; spec §6 map rule
```
$ grep -n "FaceLoadReduction_\|'FaceWorkReduction'" scripts/O2_live_balance_sympy_audit.py
225:    per_face_load = vector(action('FaceLoadReduction_' + str(a), op['J_map'],
255:    per_face_power = action('FaceWorkReduction', op['J_map'], native_domain,
```

```
$ grep -n 'reduceNative\[' mathematica/O2_live_balance_mathematica_audit.wl
160:reduceNative[object_] := response[NativeToCoordinateDensity,
239:faceWork = totalNativeFaces[reduceNative[nativeFaceWork]];
```

```
$ grep -n 'same geometric map' directives/O2_SHARED_PHYSICS.md
285:degree of freedom on which it acts, on the same measure (§1) and with the same geometric map as its
```

## C5 — c_gamma labelled OPEN in PY; spec supplies it
```
$ grep -n "action('c_gamma'" scripts/O2_live_balance_sympy_audit.py
286:    optical = sp.Tuple(sp.Eq(action('c_gamma', r)**2, mu/rho, evaluate=False),
287:                       sp.Eq(action('c_gamma', r), c0*(1+delta), evaluate=False))
```

```
$ grep -n 'Supplied optical-regime live identifications' directives/O2_SHARED_PHYSICS.md
111:| **Supplied optical-regime live identifications**, C §6 | `c_γ(r)² ≡ μ_⊥(r)/ρ_br(r)`, `c_γ(r) ≡ c₀[1+δ(r)]`. | `μ_⊥` is premise 1's optical elastic stiffness. The ratio supplies no full stress, steady-load stiffness or inertial law. `c₀` is the supplied asymptotic light speed. |
```

## C6 — hand-typed graph normal (both engines)
```
$ grep -n 'graph_normal = vector' scripts/O2_live_balance_sympy_audit.py
114:    graph_normal = vector((*(-slope), sp.S.One)) / sp.sqrt(1 + slope.dot(slope))
```

```
$ grep -n 'graphNormal = ' mathematica/O2_live_balance_mathematica_audit.wl
78:graphNormal = Join[-slope, {1}]/Sqrt[1 + slope.slope];
```

## C7 — K6 sites
```
$ grep -n "'K6[abc]'" -A1 scripts/O2_live_balance_sympy_ablation.py
37: 'K6a': ('    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6',
38-         "    momentum_profiles = sp.Tuple(sp.Symbol('V_r_constant'), rho, xi)  # KNIFE_K6"),
39: 'K6b': ('    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6',
40-         "    momentum_profiles = sp.Tuple(vr, sp.Symbol('rho_momentum_constant'), xi)  # KNIFE_K6"),
41: 'K6c': ('    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6',
42-         "    momentum_profiles = sp.Tuple(vr, rho, sp.Symbol('xi_w_constant'))  # KNIFE_K6"),
--
195:        if knife == 'K6c':
196-            engine.emit('REPAIR_K6C', repair_evidence(corrupted, engine), local=True)
```

```
$ grep -n '"K6' -A1 mathematica/O2_live_balance_mathematica_ablation.wl
26: {"K6a","momentumDerivativeRules = {}; (* SITE K6a K6b K6c *)",
27-   "momentumDerivativeRules = {VR -> Function[{z},VRConstant]}; (* SITE K6a K6b K6c *)"},
28: {"K6b","momentumDerivativeRules = {}; (* SITE K6a K6b K6c *)",
29-   "momentumDerivativeRules = {RhoBr -> Function[{z},RhoBrConstant]}; (* SITE K6a K6b K6c *)"},
30: {"K6c","momentumDerivativeRules = {}; (* SITE K6a K6b K6c *)",
31-   "momentumDerivativeRules = {XiW -> Function[{z},XiWConstant]}; (* SITE K6a K6b K6c *)"},
```

## G1 — PY recorded grades
```
$ sed -n '384,390p' scripts/O2_live_balance_sympy_audit.py
            epsilon=GM/(c0**2*r),
            optical_monomial_indices=sp.Tuple(*(sp.Tuple(a,b,c) for a in range(2)
                                               for b in range(3) for c in range(2))),
            recorded_grades=record(delta=sp.Symbol('epsilon'),
                                   slope_squared=sp.Symbol('epsilon'),
                                   velocity_over_c0=sp.sqrt(sp.Symbol('epsilon')))),
    }
```

```
$ grep -n 'ε(r) ≡\|δ = O(ε)\|V/c₀ = O' directives/O2_SHARED_PHYSICS.md
323:ε(r) ≡ GM/(c₀²r) ,
324:δ = O(ε) ,                  (∂ξ_w)² = O(ε) ,
325:V/c₀ = O(ε^{1/2}) ,         (V/c₀)² = O(ε) ,
```

