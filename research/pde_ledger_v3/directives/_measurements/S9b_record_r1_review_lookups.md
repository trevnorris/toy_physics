# Measurements: S9b record review, round 1 (generated 2026-10-19 03:25)

Generator: `_scratch/s9b_build/gen/s9b_record_r1_lookups.sh` (sha256sum/grep/sed/head/wc only). The legs'
evidence is hard-linked at `_scratch/s9b_build/record/s9b_record_review_r1_{claude,grok}_evidence/`.

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_baseline_r1.sha256 | grep -c ': OK$')
4
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_prompt_r1.md | cut -c1-64
956a83f590201681ab5440d4d70bf13ae6bd6d896e41caf1160468f641d311e0
```

## Verdicts
```
$ head -1 /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r1_claude.md
**Verdict:** NEEDS REVISION
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r1_grok.txt
Verdict: CLEAR
```

## Both legs: the generator reproduces the filed measurements
```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md | cut -c1-64
26d7090fa3990449b68e478409d1bd7f8bde51fde39919c1d713ddedf793686e
```

```
$ grep -h 'REGENERATED_IDENTICAL_TO_FILED' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r1_claude_evidence/regen/regen.stdout | head -1
REGENERATED_IDENTICAL_TO_FILED
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r1_grok_evidence/S9b_record_measurements.rerun.md | cut -c1-64
26d7090fa3990449b68e478409d1bd7f8bde51fde39919c1d713ddedf793686e
```

## Claude 1: the deflection/radar separate-test count
```
$ sed -n '109,111p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_record_directive.md

5. **The two observables.** The record states, from printed objects, which combination of inputs each observable's
   Part B condition constrains, and the printed difference between the two effective `γ`s. It counts the deflection
```

```
$ sed -n '10,11p;183,186p;256,259p;339p;581,582p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
Whether those differing fields state the same condition is OPEN. The record therefore establishes
no second independent test from the PY pairs and leaves the WL separate-test count OPEN; it does
condition is OPEN, no owner named. Consequently PY supplies no second independent profile test
from these identical pairs; WL's separate-test count remains OPEN. Matching `GM`/`ProfileEquation`
fields alone do not settle the full conditions, and differing serialization alone does not
establish independent physical tests. These limits apply to every condition comparison in
M6 counts marker-only comparisons for all 18 γ-family rows. The PY payload identity and the WL
field differences retrieved in M13 yield §3's bounded counting conclusion: no second independent
test established by PY, WL count OPEN. The printed γ relation does not settle equivalence of
the differing WL condition fields. Nothing extends this to mixed grades or full radar time.
pair (with the wrappers listed in §3); equivalence and the WL separate-test count remain OPEN.
Questions requiring work beyond retrieval remain OPEN: whether differing within-WL condition fields
state the same condition and whether the observables furnish separate WL tests; complete cross-engine equivalence or
```

```
$ sed -n '575p;838,839p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
  the WL separate-test count remains OPEN. Text differences alone do not establish independent
Whether the differing WL fields state the same condition, and the WL separate-test count,
remain OPEN. No second independent test follows from the identical PY pairs.
```

```
$ grep -h 'GM!=0 member literal present' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r1_claude_evidence/gamma_difference_census.stdout | head -12
PY_S9B_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
PY_S9B_R_FLOW_ONLY_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
PY_S9B_R_SPEED_ONLY_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
PY_S9B_R_TILT_ONLY_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
PY_S9B_F_FLOW_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
PY_S9B_F_LIVE_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
WL_S9B_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
WL_S9B_R_FLOW_ONLY_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
WL_S9B_R_SPEED_ONLY_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
WL_S9B_R_TILT_ONLY_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
WL_S9B_F_FLOW_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
WL_S9B_F_LIVE_B_GAMMA_DIFFERENCE | GM!=0 member literal present: True | occurrences 1
```

```
$ grep 'PY_S9B_B_GAMMA_DIFFERENCE' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | head -1 | grep -o "Unequality(Symbol('GM', real=True), Integer(0))), FiniteSet(Integer(0))"
Unequality(Symbol('GM', real=True), Integer(0))), FiniteSet(Integer(0))
```

```
$ grep 'WL_S9B_B_GAMMA_DIFFERENCE' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | head -1 | grep -o '(gammaDifference == 0 && GM != 0)'
(gammaDifference == 0 && GM != 0)
```

## Claude 2: the flow-only restriction as SymPy prints it
```
$ sed -n '189,216p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
**Implied exchange and domain.** Speed-only and tilt-only both print `jn=0`, with `ρbr` still
symbolic. Flow-only prints both velocity signs on each positive-squared-velocity component, the
zero intervals, differentiable zero-set joins, and excludes negative squared velocity. WL prints

```
Y = -r Inactive[Integrate][(-4 GM)/tailRadius^3, {tailRadius,r,Infinity}],
V = -√Y or +√Y,
```

with its general `r velocityIntegrationConstant` term excluded by supplied asymptotic counting,
and the retained constant domain `velocityIntegrationConstant==0`. Its two positive-component
exchange numerators are the negative and positive versions, respectively, of the literal field

```
(rhoBr[r]*(4*GM + 5*r^2*Inactive[Integrate][(-4*GM)/tailRadius^3,
 {tailRadius,r,Infinity}]) + 2*r^3*Inactive[Integrate][(-4*GM)/tailRadius^3,
 {tailRadius,r,Infinity}]*Derivative[1][rhoBr][r]) /
 (2*r^2*Sqrt[-r*Inactive[Integrate][(-4*GM)/tailRadius^3,
 {tailRadius,r,Infinity}]])
```

No integral is evaluated here. PY separately prints its signed-component and zero-locus exchange
objects. The common coordinate mass expression is printed engine-locally as
`jn=−ρbr V′−Vρbr′−2ρbr V/r`; each restriction retains its own computed substitutions, not a
constant-density shortcut. WL prints component-boundary limits and zero-interval `jn=0`; PY prints
one-sided differentiable-join conditions and the corresponding exchange. Constants with different
engine names are not identified. (M3/M4 `LOCAL_MASS_BALANCE`, `R_FLOW_ONLY_B_IMPLIED_JN_*`;
full WL restrictions M11; cross-engine bounds M6/M8.)
```

```
$ sed -n '10090,10100p;10159,10161p' /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION line=137 path=1/0/0
Str('ON_EACH_CONNECTED_COMPONENT_OF_POSITIVE_V_SQUARED')
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION line=137 path=1/0/2
Equality(Function('s9b_V', **{'complex': True, 'extended_real': True, 'finite': True, 'hermitian': True, 'imaginary': False, 'infinite': False, 'real': True})(Symbol('s9b_r', positive=True)), Mul(Integer(-1), Pow(Integer(2), Rational(1, 2)), Pow(Symbol('GM', real=True), Rational(1, 2)), Pow(Symbol('s9b_r', positive=True), Rational(-1, 2))))
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION line=137 path=1/0/3
Equality(Function('j_n', **{'complex': True, 'extended_real': True, 'finite': True, 'hermitian': True, 'imaginary': False, 'infinite': False, 'real': True})(Symbol('s9b_r', positive=True)), Add(Mul(Pow(Integer(2), Rational(1, 2)), Pow(Symbol('GM', real=True), Rational(1, 2)), Pow(Symbol('s9b_r', positive=True), Rational(-1, 2)), Derivative(Function('rho_br_profile', **{'complex': True, 'extended_n
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION line=137 path=1/1/0
Str('ON_EACH_CONNECTED_COMPONENT_OF_POSITIVE_V_SQUARED')
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION line=137 path=1/1/2
Equality(Function('s9b_V', **{'complex': True, 'extended_real': True, 'finite': True, 'hermitian': True, 'imaginary': False, 'infinite': False, 'real': True})(Symbol('s9b_r', positive=True)), Mul(Pow(Integer(2), Rational(1, 2)), Pow(Symbol('GM', real=True), Rational(1, 2)), Pow(Symbol('s9b_r', positive=True), Rational(-1, 2))))
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION line=137 path=1/1/3
Str('squared')
PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION line=137 path=2/2/1
Mul(Integer(2), Symbol('GM', real=True), Pow(Symbol('s9b_r', positive=True), Integer(-1)))
```

```
$ grep -c '2GM\|2\*GM\|2 GM\|2·GM' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
2
```

```
$ grep -n 'R_FLOW_ONLY' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | head -5
568:  `VelocityReduction/SquaredVelocityEquation` additionally differs for general B, R_FLOW_ONLY B
581:  `C_*_IMPLIED_JN_*`, and the `R_FLOW_ONLY_`, `R_SPEED_ONLY_`, `R_TILT_ONLY_`,
```

```
$ grep 'PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | head -1 | grep -o "Tuple(Str('squared'), Mul(Integer(2), Symbol('GM', real=True), Pow(Symbol('s9b_r', positive=True), Integer(-1))))"
Tuple(Str('squared'), Mul(Integer(2), Symbol('GM', real=True), Pow(Symbol('s9b_r', positive=True), Integer(-1))))
```

## Claude 3: what each engine's BRANCH_EXISTENCE tag binds
```
$ sed -n '83,89p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
The native gates require real coefficients, nonzero material density, a real positive optical
ratio `μ⊥/ρbr`, its equality to `c₀²(1+δ)²` with the nonnegative speed, and the ray's own turning
and angular-closure construction. Traversal additionally requires the printed along-ray inequality
`μ⊥/ρbr > V²(1+ξw′²)`. The turning roots, outward monotonicity and finite-endpoint closure stay
attached; a pointwise inequality alone is not the full domain. WL uses its turning-radius/impact
conditions; PY uses its constructed ray-root sets. Their equivalence is not established by the
comparator. PY declares `rho_br_profile` positive; WL's gate states real, nonzero `rhoBr` and the
```

```
$ grep -n 'BRANCH_EXISTENCE' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md | head -6
374:| BRANCH_EXISTENCE, PATH_TRAVERSAL, BRANCH_TYPE and conditional F_* copies | Stored objects and attempted predicates | No formed along-ray comparison; no equivalence of type/turning-ray constructions |
```

```
$ grep '^PY_S9B_F_FLOW_BRANCH_EXISTENCE' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | head -1
PY_S9B_F_FLOW_BRANCH_EXISTENCE: true
```

```
$ grep 'WL_S9B_F_FLOW_BRANCH_EXISTENCE' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out | head -1 | grep -c 'rTurnFlyby'
1
```

```
$ sed -n '262,266p;275p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
  physicalPropagation && speedChart && rhoBr[r] != 0 &&
  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];
traversalGate = rayConstruction && Inactive[ForAll][r, traversedRadii, physicalTraversal];
rayGate = propagationGate && traversalGate;
  emit[prefix <> "BRANCH_EXISTENCE", wrapper[propagationGate /. rules]];
```

```
$ sed -n '696,697p;720p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | cut -c1-260
    real_domain = everywhere(r,existence.subs(branch_map),traversed)
    travel_domain = sp.And(*(sp.Ne(ray['roots'],sp.EmptySet,evaluate=False) for ray in rays.values()),
    emit('BRANCH_EXISTENCE',real_domain)
```

```
$ grep -h 'F_FLOW_BRANCH_EXISTENCE' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r1_claude_evidence/value_predicates.stdout | head -3
F_FLOW_BRANCH_EXISTENCE comparison boolean_residual= ~((wl::rTurnFlyby > s9b_b_min) & (wl::rTurnRadar > s9b_b_min) & (wl::rTurnRadar <= Min(sqrt(Z_E**2 + b**2), sqrt(Z_R**2 + b**2))) & Eq(sqrt(c_0**2*wl::rTurnFlyb formed= False {'unresolved coverage': 37} py_cond= True wl_cond= (wl::rTurnFlyby > s9b_b_min) & (wl::rTurnRadar > s9b_b_min) & (wl::rTurnRadar <= Min(sqrt(Z_E**2 + b**2), sqrt(Z_R**2 + 
```

## Claude 4: R-S8-03 (the sign of the transverse stiffness)
```
$ sed -n '318,324p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
### R-S8-03 — the SIGN of the physical transverse stiffness

- **source** S9, S10; S11 (`steps/S11_stray_longitudinal.md`, unchanged transverse branch);
  unified S11b (`steps/S11b_interface_coupling_law.md`, transverse mode); S11c
  (`steps/S11c_PARTIAL_CLOSEOUT.md`, conditional uniform result;
  `steps/S11c_d_profile_conditioned_scattering.md`, uniform check) · **target** S8 · **status** OPEN
- **requirement** — the **sign** of the physical transverse stiffness in the quadratic brane
```

```
$ sed -n '334,336p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
  so the pass consolidated one and dropped the other.
  ⚠ **It is not closed by defect `B2`**, which shuts the route for deriving `μ_R` — its sign as well as
  its magnitude — from a polar substructure `P`. Other substrate routes remain open: a step that
```

```
$ sed -n '824p;850p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
S21, `R-S21-01`.** Existing entries gaining a source: **`R-S8-01`, `R-S8-06`**.
| Positive-real optical ratio/constructed ray; restricted response domains | Travel with R-S21-01. No new uniform-stiffness sign result or expanded R-S8-03 domain is inferred; PY/WL brane-density assumptions remain an open domain difference. |
```

```
$ sed -n '83,84p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
The native gates require real coefficients, nonzero material density, a real positive optical
ratio `μ⊥/ρbr`, its equality to `c₀²(1+δ)²` with the nonnegative speed, and the ray's own turning
```

```
$ grep -n 'R-S8-03' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | head -10
80:| **S8** | `R-S8-03` | the **sign** of the physical transverse stiffness | S9, S10, S11, unified S11b, S11c uniform |
318:### R-S8-03 — the SIGN of the physical transverse stiffness
616:sources (`R-S7-01`, `R-S1-01`, `R-S1.5-01`, `R-S6-02`, `R-S8-01`, `R-S8-02`, `R-S8-03`). New entries:
850:| Positive-real optical ratio/constructed ray; restricted response domains | Travel with R-S21-01. No new uniform-stiffness sign result or expanded R-S8-03 domain is inferred; PY/WL brane-density assumptions remain an open domain difference. |
```

