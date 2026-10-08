# Measurements — S9b spec v10 review round 1 dispositions (generated 2026-10-08 11:02)

Generator: `_scratch/s9b_build/gen/s9b_spec_v10_r1_lookups.sh` (sed/grep/sha256sum/git only; lines cut at 330 characters). The reviewed files are the working copies checked against the baseline below.

```
$ sha256sum -c _scratch/s9b_build/s9b_spec_v10_review_baseline.sha256
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_v10_spec_lookups.md: OK
```

```
$ grep -n -o 'Verdict:[* ]*[A-Z][A-Z ]*' _scratch/s9b_build/s9b_spec_v10_review_r1_claude.md _scratch/s9b_build/s9b_spec_v10_review_r1_grok.txt
_scratch/s9b_build/s9b_spec_v10_review_r1_claude.md:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_spec_v10_review_r1_grok.txt:1:Verdict: CLEAR
```

## C1 — P5 has no steady content without a sourced momentum transport
Both O2 engines keep storage (time derivative of the density) and transport (divergence of a separate current) as separate roles:
```
$ grep -n 'storage = density\|^    transport = \|internal = vector\|entries = ((1, storage)' research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py
197:    storage = density.applyfunc(lambda p: derivative(p, t))
198:    transport = vector(sum(derivative(flux[a, i], x[i]) for i in range(3))
211:    internal = vector(action('InternalForce_' + str(a), stress_state) for a in range(4))
259:    entries = ((1, storage), (1, transport), (1, carry), (1, partners),
```

```
$ grep -n 'storage *=\|transport *=\|internalForce *=' research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl
152:storage = sectionDerivative[#, momentumDifferentiatedSection, t] & /@ momentumDensity;
153:transport = Table[Sum[sectionDerivative[momentumCurrent[[a,i]],
168:internalForce = Table[response[InternalForce[a],stressInputs,section],{a,4}];
```

The O2 record requires any advective kinetic input to be sourced by Part D's own premise work:
```
$ sed -n '509,513p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
In particular, O2 supplies **no momentum-density map `ρ_br V`**. Its open inertia/momentum mapping
cannot be cancelled or converted into a chosen advective kinetic law using the mass relation.
The directive's folded disposition keeps the separate user premise for S9b repair at that later
decision-list gate (M1). Part D must source any additional kinetic, constitutive, projection,
profile, support, order or energy input through its own authorized premise/specification work.
```

What the spec says P5 supplies, and what it leaves OPEN:
```
$ sed -n '322,331p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**P5 — momentum density (2026-10-07; D2).** The supplied brane in-plane momentum-density map is

```
𝒫_br,inplane^live ≡ ρ_br V .
```

This supplies the previously OPEN in-plane momentum-density identification. If stressed brane
material carried additional momentum from its stress, P5 would change. Whether that enters at a
retained grade is OPEN and owned by S8. P5 supplies no total-energy law, normal-response map or
independent closure of every O2 transport/history action.
```

```
$ grep -n '^| Material momentum storage' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
383:| Material momentum storage and transport; `ℐ_br^live`, `𝒫_br^cons` | P5 supplies in-plane momentum density. | Remaining flowing/embedded inertia, normal and transport/history actions and their unresolved relations. The conservative antecedent is named, rather than added as another momentum species. |
```

```
$ grep -n '^| Internal material force' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
384:| Internal material force; `𝒯_br^live`, `𝒯_br^cons`, `𝒩_br^live`, `𝒜_rot^live`, `ℛ_ref/strain^live` | P6 supplies the steady in-plane pressure part of the full stress. | Normal material response, conservative antecedents, reference evolution and rotational/couple/frame content wherever unsupplied. Their overlap
```

```
$ sed -n '443,446p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
P5 supplies the in-plane momentum density, P6 the steady in-plane pressure stress, and the density
link the optical stiffness response, each with the dependence stated in its adopted equation.
Their use in a constructed flux/action is flagged; they are not blanket replacements of every O2
momentum-flux or energy-flux action. Keep the general OPEN objects, their full dependence declarations,
```

```
$ grep -n -i 'advect\|momentum current\|momentum flux' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
53:- **Advection (supplied).** `V(x)` is the in-plane velocity of the material that `u` displaces. Its
99:  `χ` is the inverse material map. MATERIAL_ADVECTED is not selected. LAB_HELD supplies neither a
106:- mixed `ωk` content other than advection.
```

The decision list leaves the determination to the spec author:
```
$ grep -n 'The spec author determines' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
75:  - The spec author determines, from the O2 sources, which content each premise supplies.
```

```
$ grep -n '^| \*\*P5\*\*' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
27:| **P5** (2026-10-07) | The brane's momentum density is `ρ_br V`. | If the stressed brane material carried additional momentum from its stress, P5 would change. Whether that enters at a retained grade is OPEN; S8 owns it. |
```

The leg's script output (steady storage; an OPEN remainder current absorbs the P5/P6 content; the FORM-control difference depends on the carrying velocity U_r — counted, not printed):
```
$ grep -n 'storage components\|residual B_with\|== (3) FORM control' _scratch/s9b_build/s9b_spec_v10_review_r1_claude_evidence/r1_partD_transport.stdout
2:storage components: [0, 0, 0]
11:== (3) FORM control: momentum advected at a DIFFERENT velocity U (not V) ==
16:component 0 residual B_with - B_without(F'): 0
17:component 1 residual B_with - B_without(F'): 0
18:component 2 residual B_with - B_without(F'): 0
```

```
$ grep 'difference to B_adv' _scratch/s9b_build/s9b_spec_v10_review_r1_claude_evidence/r1_partD_transport.stdout | grep -c 'U_r'
1
```

## C2 — P6's stated reason is not entailed by P1
```
$ sed -n '333,335p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**P6 — steady in-plane pressure (2026-10-08; D2).** In Part D, the supplied steady in-plane stress
is isotropic pressure, as the P1 material relaxes under steady load. In the Cartesian Cauchy-stress
convention where stress contracted with a unit normal gives traction, the premise is
```

```
$ grep -n 'because shear relaxes under steady load' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
28:| **P6** (2026-10-08) | In steady flow the brane's in-plane stress is an isotropic pressure `p_br(ρ_br)` that depends on `ρ_br` only, because shear relaxes under steady load (P1). `p_br` stays a general function. Its compressional speed `c_comp`, with `c_comp² = dp_br/dρ_br`, is live wherever `ρ_br` varies. In the optica
```

```
$ sed -n '17p' research/pde_ledger_v3/directives/O2_premise_decision_list.md
   optical shear regime and relaxes under steady load. The user chose this on 2026-10-06 over a separately offered
```

The leg's script output (steady linear Maxwell counterexample; FORM control V_r = C R):
```
$ grep -n 'Maxwell residual\|profiles with edev\|V_r = ' _scratch/s9b_build/s9b_spec_v10_review_r1_claude_evidence/r2_P6_steady_flow_stress.stdout
7:Maxwell residual with sigma'=0, [0,0] on x1-axis: 4*G*tau*(R*Derivative(V_r(R), R) - V_r(R))/(3*R)
8:profiles with edev == 0: Eq(V_r(R), C1*R)
15:V_r = A/R^2 (div-free) | edev00 = -2*A/R**3 | f_visc const eta = 0 | f_visc eta(R) = -4*A*Derivative(eta(R), R)/R**3
16:V_r = A/R^2 + B/R (div != 0) | edev00 = 2*(-3*A - 2*B*R)/(3*R**3) | f_visc const eta = -8*B*eta0/(3*R**3) | f_visc eta(R) = 4*(-3*A*Derivative(eta(R), R) - 2*B*R*Derivative(eta(R), R) - 2*B*eta(R))/(3*R**3)
17:V_r = C R (FORM control: edev=0) | edev00 = 0 | f_visc const eta = 0 | f_visc eta(R) = 0
```

## C3 — the blind Wolfram engine cannot retain operands it is never given
```
$ sed -n '393,396p;468p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
locality, finite internal-variable set, derivative order, stress split or constitutive family.
Retain every named operand and complete live-object dependence printed by either O2 engine in the
content no adopted premise supplies, including PY-only native/chart and core/material-compatibility
content. Neither engine's finite inventory, nor their union, exhausts admissible dependences.
- **Engines.** SymPy, plus a blind Wolfram engine that imports nothing. No Lean (CLAUDE.md L5).
```

```
$ sed -n '469p' research/pde_ledger_v3/steps/O2_steady_brane_balance.md
Part D must retain **every named operand and complete live-object dependence printed by either engine**,
```

```
$ for n in OpenDomain StateHistory NativeFaceChartDomain material_action_compatibility energy_accounting_overlap face_support_partition CarriedMomentumW UnfixedNormalGeneralizedRates RotationalGeneralizedWork NativeFaceReductionMap; do printf '%s spec=%s py=%s wl=%s\n' $n $(grep -c $n research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md) $(grep -c $n research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py) $(grep -c $n research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl); done
OpenDomain spec=0 py=2 wl=0
StateHistory spec=0 py=1 wl=0
NativeFaceChartDomain spec=0 py=1 wl=0
material_action_compatibility spec=0 py=2 wl=0
energy_accounting_overlap spec=0 py=2 wl=0
face_support_partition spec=0 py=2 wl=0
CarriedMomentumW spec=0 py=1 wl=0
UnfixedNormalGeneralizedRates spec=0 py=1 wl=0
RotationalGeneralizedWork spec=0 py=1 wl=0
NativeFaceReductionMap spec=0 py=1 wl=0
```

## C4 — the induced-measure rule: two readings; the D4 condition gives no relative order
```
$ sed -n '65,71p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  The divergence and both densities are on the coordinate `d³x` measure of the far-field `x^i`
  coordinates (O2-R §§2, 8; O2-S §§1, 3.1). `μ_⊥` in the optical ratio is on the same measure as
  `ρ_br`. This input supplies no induced-measure or finite-slab mass law. For an induced-measure
  interpretation, keep `∂_r[(∂ξ_w)²]` live. D4 supplies the optional additional scale condition
  `∂_r[(∂ξ_w)²] = O(ε/r)`; it is not imposed here. No order for the relative correction to `j_n`
  is supplied from the slope-amplitude grade alone.

```

```
$ grep -n 'That holds only if\|No order is attached\|or states the condition above' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
55:  correction to the implied `j_n`". That holds only if `∂_r[(∂ξ_w)²] = O(ε/r)` (v9 review, Grok).
57:  the induced measure carries the scale `∂_r[(∂ξ_w)²]` live, or states the condition above. No order is attached
```

The leg's script output (Claude R3):
```
$ grep -n 'check vs\|sqrt(g)\*div_cov\|case rho' _scratch/s9b_build/s9b_spec_v10_review_r1_claude_evidence/r3_induced_measure_jn.stdout
2:check vs -rho V_r d_R ln sqrt(g): 0
5:sqrt(g)*div_cov(rho~) - div_coord(rho): 0
9:case rho V_r R^2 = F0: j_n(coord) = 0 ; absolute difference = F0*R_s*a/(2*R**3*(R + R_s*a))
10:case rho V_r R^2 = F0 R^c: j_n(coord) = -F0*R**(c - 3)*c ; relative difference = -R_s*a/(2*c*(R + R_s*a)) ; leading as R_s/R->0: -a*u/(2*c)
```

The Grok leg's own script reports the same absolute difference with zero coordinate j_n:
```
$ grep -n -i 'j_cov\|j_coord' _scratch/s9b_build/s9b_spec_v10_review_r1_grok_evidence/s9b_v10_review_physics.stdout
7:j_cov - j_coord - (-rho*V_r*d_r ln sqrt(g)) = 0
9:j_coord = -V_r(r)*Derivative(rho(r), r) - rho(r)*Derivative(V_r(r), r) - 2*V_r(r)*rho(r)/r
10:j_cov - j_coord = -V_r(r)*rho(r)*Derivative(s(r), r)/(2*s(r) + 2)
18:j_coord(constant flux) = 0
19:j_cov(constant flux) = C*k/(2*r**3*(k + r))
21:j_flat(rho_coord) - sqrt(g)*j_cov(rho_prop) = 0
```

