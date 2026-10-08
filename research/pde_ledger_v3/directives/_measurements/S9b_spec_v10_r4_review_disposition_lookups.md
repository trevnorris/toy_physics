# Measurements — S9b spec v10 review round 4 dispositions (generated 2026-10-08 14:17)

Generator: `_scratch/s9b_build/gen/s9b_spec_v10_r4_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). The reviewed files are the working copies checked against the baseline below.

```
$ sha256sum -c _scratch/s9b_build/s9b_spec_v10_review_baseline_r3.sha256
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_v10_spec_lookups.md: OK
```

```
$ grep -n -o 'Verdict:[* ]*[A-Z][A-Z ]*' _scratch/s9b_build/s9b_spec_v10_review_r4_claude.md
3:Verdict: NEEDS REVISION
```

```
$ grep -n -o 'Findings.\*\* None' _scratch/s9b_build/s9b_spec_v10_review_r4_grok.txt
5:Findings.** None
```

## F1 — the coordinate-only rule is unscoped; Part D keeps native face-area and measure content OPEN
```
$ sed -n '64,70p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`
  coordinates (O2-R §§2, 8; O2-S §§1, 3.1). `μ_⊥` in the optical ratio is on the same measure as
  `ρ_br`. This input supplies no induced-measure or finite-slab replacement mass law.
  S9b computes and claims on the coordinate `d³x` measure only. No induced-measure object,
  re-expression or correction is computed or printed. The induced-measure readings and
  O2-R L478–480's recorded qualification for those readings are outside this step's
  computation and claims (D4).
```

```
$ grep -n 'Native face-area factors' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
261:law. Native face-area factors and reductions remain explicit through O6. Coordinate `w` content and
```

```
$ grep -n -o 'native face geometry/measure' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
416:native face geometry/measure
```

```
$ grep -n -o 'chart/measure/map' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
488:chart/measure/map
```

```
$ git show f416bb0d:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | grep -n 'dvol_g'
73:    dvol_g ≡ √det(g_ij) d³x ,
74:    ρ_br^(g) dvol_g ≡ ρ_br d³x ,      j_n^(g) dvol_g ≡ j_n d³x .
```

```
$ grep -n -i -c 'dvol' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
0
```

```
$ git show 217a92e9:research/pde_ledger_v3/directives/O2_input_contract.md | sed -n '313p'
a_s^α = sqrt(1+|∇_x h_s^α|²) ,
```

```
$ git show 4680e251:research/pde_ledger_v3/directives/O2_SHARED_PHYSICS.md | sed -n '61,64p'
`c_γ² ≡ μ_⊥/ρ_br`, `μ_⊥` is referred to the same measure as `ρ_br`. Induced-metric factors, such as
`g_ij`, `g^{ij}`, `det g_ij` and the graph normal, enter explicitly in the geometric actions; no
occurrence is moved to another measure. A quantity defined per native face area enters through its
native geometric factor and `𝒥_map` (§5). This fixes the volume measure of densities only. It
```

```
$ grep -n -o 'native_area`, `native_normal`, `native_measure' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
375:native_area`, `native_normal`, `native_measure
```

```
$ grep -n 'PY-only native/chart' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
470:including all eight WL-only derivative keys at every location in §4 and all PY-only native/chart,
```

```
$ grep -n -o 'use the same material and measure' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
486:use the same material and measure
```

The repair-3 brief's first must-be-true item (the orchestrator's wording):
```
$ grep -n 'no induced-measure object' _scratch/s9b_build/s9b_spec_v10_repair3_prompt.md
15:1. The spec asks the engines to compute or print no induced-measure object, re-expression or correction.
```

The leg's script output:
```
$ grep -n 'FACE_AREA_FACTOR\|SUM_s\|SYMMETRIC_FACES\|CENTRE_GRAPH_DET' _scratch/s9b_build/s9b_spec_v10_review_r4_claude_evidence/geom_measure.stdout
4:CENTRE_GRAPH_DET (radial): Subs(Derivative(xi_w(_xi_1), _xi_1), _xi_1, sqrt(x1**2 + x2**2 + x3**2))**2 + 1
6:CENTRE_GRAPH_DET at xi_w==0: 1
7:SYMMETRIC_FACES centre graph (h_+ + h_-)/2: 0
8:FACE_AREA_FACTOR a_+ : sqrt(Derivative(Hs(x1, x2, x3), x1)**2 + Derivative(Hs(x1, x2, x3), x2)**2 + Derivative(Hs(x1, x2, x3), x3)**2 + 1)
9:FACE_AREA_FACTOR a_- : sqrt(Derivative(Hs(x1, x2, x3), x1)**2 + Derivative(Hs(x1, x2, x3), x2)**2 + Derivative(Hs(x1, x2, x3), x3)**2 + 1)
11:SUM_s T n_s (components x1,x2,x3,w): [-2*T_bulk_n*Derivative(Hs(x1, x2, x3), x1)/sqrt(Derivative(Hs(x1, x2, x3), x1)**2 + Derivative(Hs(x1, x2, x3), x2)**2 + Derivative(Hs(x1, x2, x3), x3)**2 + 1), -2*T_bulk_n*Derivative(Hs(x1, x2, x3), x2)/sqrt(Derivative(Hs(x1, x2, x3), x1)**2 + Derivative(Hs(x1, x2, x3), x2)**2 + Derivativ
12:SUM_s a_s T n_s (per coordinate area): [-2*T_bulk_n*Derivative(Hs(x1, x2, x3), x1), -2*T_bulk_n*Derivative(Hs(x1, x2, x3), x2), -2*T_bulk_n*Derivative(Hs(x1, x2, x3), x3), 0]
13:SUM_s T n_s radial face profile Hr(r): [-2*T_bulk_n*x1*Subs(Derivative(Hr(_xi_1), _xi_1), _xi_1, sqrt(x1**2 + x2**2 + x3**2))/sqrt(x1**2*Subs(Derivative(Hr(_xi_1), _xi_1), _xi_1, sqrt(x1**2 + x2**2 + x3**2))**2 + x1**2 + x2**2*Subs(Derivative(Hr(_xi_1), _xi_1), _xi_1, sqrt(x1**2 + x2**2 + x3**2))**2 + x2**2 + x3**2*Subs(Deriv
```

## Notes — version line and provenance
```
$ sed -n '3p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Version:** v10, repair 2 (2026-10-08).
```

```
$ grep -n 'outside this step' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | head -3
69:  O2-R L478–480's recorded qualification for those readings are outside this step's
152:  outside this step, the `v_dr` freeze is open, and `v_dr` is distinct from both the live in-plane `V` and
```

