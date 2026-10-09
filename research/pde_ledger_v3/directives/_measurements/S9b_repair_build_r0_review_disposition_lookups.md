# Measurements — S9b repair build review r0 (generated 2026-10-09 04:19)

Generator: `_scratch/s9b_build/gen/s9b_build_r0_lookups.sh` (sha256sum/grep/sed/sort/comm/count only).
Baseline streams: the builders' final guarded baseline runs (PY engine sha recorded in the run's inventory).

## Reviewed artifacts and verdicts
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r0.sha256
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/scripts/S9b_exports.py: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_staged_candidates.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_not_evaluated_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_not_established_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_executable_checks.json: OK
```

```
$ grep -o '"engineSHA256": "[0-9a-f]*"' research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json | head -1
"engineSHA256": "16b00f137623c897a5f890c37fa5ea0082e8c1156c9907409b2884603a6290d9"
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r0_claude.md _scratch/s9b_build/s9b_repair_build_review_r0_grok.txt
_scratch/s9b_build/s9b_repair_build_review_r0_claude.md:5:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_build_review_r0_grok.txt:5:Verdict: NEEDS REVISION
```

## C1/G2: quantity names per engine (prefix stripped), and their intersection
```
$ grep -o '^PY_S9B_[A-Z0-9_]*' /var/projects/toy_physics-s9b-r2-py/_scratch/s11c/s9b-amend1-baseline-02/stdout | sed 's/^PY_S9B_//' | sort -u | wc -l
213
```

```
$ grep -o '^WL_S9B_[A-Z0-9_]*' /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-02-finished-repository/stdout | sed 's/^WL_S9B_//' | sort -u | wc -l
181
```

```
$ comm -12 <(grep -o '^PY_S9B_[A-Z0-9_]*' /var/projects/toy_physics-s9b-r2-py/_scratch/s11c/s9b-amend1-baseline-02/stdout | sed 's/^PY_S9B_//' | sort -u) <(grep -o '^WL_S9B_[A-Z0-9_]*' /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-02-finished-repository/stdout | sed 's/^WL_S9B_//' | sort -u)
B_DEFLECTION_CONDITION
B_RADAR_CONDITION
BRANCH_EXISTENCE
BRANCH_TYPES
FORWARD_NO_FAR_ZONE_LOSS_MASS_SOLUTION
MASS_BALANCE
OBSERVABLE_DOMAIN
PATH_TRAVERSAL
RAY_DOMAIN
SUPPLIED_DEPENDENCIES
```

```
$ grep -n 'same set of' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
53:   the same set of `<QUANTITY>` names, with one tag per named object and one line per tag: `TAG: <payload>`,
```

## C2: SymPy dispersion symbol and the item-8 one-form construction
```
$ sed -n '443p;618,627p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    kinetic = (omega-U*kr)**2  # K1: kinetic construction
    radius = sp.sqrt(sum(a*a for a in xyz))
    beta_general = sum(betaj.c.values(),sp.S.Zero).subs(r, radius)
    azimuthal_velocity = sp.zeros(3, 1)  # K11: advection construction
    advected_velocity = sp.Matrix([velocity.subs(r, radius)*sp.diff(radius, a) for a in xyz])+azimuthal_velocity
    spatial_metric = sp.eye(3)+slope.subs(r, radius)**2*sp.Matrix([sp.diff(radius, a) for a in xyz])*sp.Matrix([sp.diff(radius, a) for a in xyz]).T
    # Randers odd form is obtained by completing the same advected quadratic.
    speed_cartesian = (c0*(1+delta.subs(r, radius)))**2
    odd_cartesian = -spatial_metric*advected_velocity/(speed_cartesian-(advected_velocity.T*spatial_metric*advected_velocity)[0])
    oneform = [beta_general*sp.diff(radius, a)+sp.factor(odd_cartesian[i]-odd_cartesian[i].subs({z: 0 for z in azimuthal_velocity.free_symbols if z.name == 's9b_azimuthal'})) for i,a in enumerate(xyz)]
    curl = sp.ImmutableMatrix([sp.simplify(sp.diff(oneform[(i+2)%3], xyz[(i+1)%3])-
```

```
$ grep -n 'kinetic' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | head -20
443:    kinetic = (omega-U*kr)**2  # K1: kinetic construction
447:    dispersion = kinetic-local_speed_squared*(covector.T*inverse_metric*covector)[0]
```

## C3: the declared profile families and the spec's Profiles item
```
$ grep -n 'PROFILE_FAMILY\|PROFILE_ANSATZ' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl | head
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py:459:    emit('PROFILE_FAMILY', (eq(sp.Function('delta')(r), delta), eq(sp.Function('V')(r), velocity),
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl:90:emit["PROFILE_ANSATZ", <|"Delta" -> delta, "V" -> velocity,
```

```
$ sed -n '137,140p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
- **Profiles.** Prefer general radial functions in Parts A–C. If an engine restricts those objects to a
  family, it says so and keeps every exponent and coefficient symbolic. Part D's unsupplied responses
  remain general live unknowns, including admissible gradients and material history; no engine-chosen
  response family or derivative/history cutoff is permitted (D3; O2-R §8).
```

## C4: the SymPy flow-only deflection tag's line carries the zeroed profiles' amplitude symbols
```
$ grep '^PY_S9B_FLOW_ONLY_DELTA_XI_ZERO_DEFLECTION' /var/projects/toy_physics-s9b-r2-py/_scratch/s11c/s9b-amend1-baseline-02/stdout | grep -o "Symbol('s9b_[dw]'" | sort | uniq -c
    108 Symbol('s9b_d'
    135 Symbol('s9b_w'
```

## G1: rayRadiusDomain uses and assignments in the Wolfram engine
```
$ grep -n 'rayRadiusDomain' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
138:  "AlongRay" -> Inactive[ForAll][r, rayRadiusDomain[r],
143:  "SubcriticalDomain" -> Inactive[ForAll][r, rayRadiusDomain[r], physicalTraversal]|>];
198:   perturbative branch of the boundary-value solution.  rayRadiusDomain
319:rayGate = Inactive[ForAll][r, rayRadiusDomain[r],
523:  forwardEmit[prefix <> "BRANCH_EXISTENCE", Inactive[ForAll][r, rayRadiusDomain[r], physicalPropagation /. stageRules]];
524:  forwardEmit[prefix <> "PATH_TRAVERSAL", Inactive[ForAll][r, rayRadiusDomain[r], physicalTraversal /. stageRules]];
```

```
$ grep -c 'rayRadiusDomain\[[^]]*\] *:\?=' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
0
```

## G3: the Wolfram radar-gamma construction and its baseline payload head
```
$ sed -n '369,373p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
referenceLog = Coefficient[Expand[referenceRadar /. Log[4 rE rR/b^2] -> logBasis], logBasis];
gammaTheta = gamma /. First[Solve[referenceTheta == thetaFirst, gamma]];
gammaRadar = gamma /. First[Solve[referenceLog == radarFirst, gamma]];
thetaResidual = thetaFirst - (referenceTheta /. gamma -> 1);
radarResidual = radarFirst - (referenceLog /. gamma -> 1);
```

```
$ grep -o '^WL_S9B_B_GAMMA_RADAR: .\{0,40\}' /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-02-finished-repository/stdout
WL_S9B_B_GAMMA_RADAR: <|"OnDomain" -> ConditionalExpression[-1
```

```
$ grep -o '^WL_S9B_[A-Z_]*GAMMA_RADAR: .\{0,40\}' /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-02-finished-repository/stdout
WL_S9B_B_GAMMA_RADAR: <|"OnDomain" -> ConditionalExpression[-1
WL_S9B_RESTRICTION_FLOW_ONLY_GAMMA_RADAR: <|"OnDomain" -> ConditionalExpression[-1
WL_S9B_RESTRICTION_SPEED_ONLY_GAMMA_RADAR: <|"OnDomain" -> ConditionalExpression[-1
WL_S9B_RESTRICTION_TILT_ONLY_GAMMA_RADAR: <|"OnDomain" -> ConditionalExpression[-1
WL_S9B_FORWARD_NO_FAR_ZONE_LOSS_FLOW_GAMMA_RADAR: <|"Premise" -> <|"NormalExchange" -> jn[
WL_S9B_FORWARD_NO_FAR_ZONE_LOSS_LIVE_OPTICS_GAMMA_RADAR: <|"Premise" -> <|"NormalExchange" -> jn[
```

