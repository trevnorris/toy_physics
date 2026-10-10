# Measurements — S9b repair build, review round 1 (generated 2026-10-09 18:22)

Lookups file name: `S9b_repair_build_r1_review_disposition_lookups.md`. Generator: `_scratch/s9b_build/gen/s9b_build_r1_lookups.sh` (sha256sum/grep/sed/cat only). Paths under `/tmp` are
the Claude leg's working copies and outputs; they are its evidence, retrieved here verbatim.

## The reviewed artifacts are the working-tree files, unchanged after both legs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r1.sha256
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py: OK
research/pde_ledger_v3/scripts/S9b_exports.py: OK
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
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_checks.json: OK
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r1_claude.txt _scratch/s9b_build/s9b_repair_build_review_r1_grok.txt
_scratch/s9b_build/s9b_repair_build_review_r1_claude.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_build_review_r1_grok.txt:5:Verdict: NEEDS REVISION
```

```
$ grep -n '^\*\*Verdict' _scratch/s9b_build/s9b_repair_build_review_r1_grok.txt
```

## Claude finding 1: the implied j_n is built from the profile replacement, not from the condition
### What the documents ask for
```
$ sed -n 234p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  - For each condition, the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic.
```

```
$ sed -n 242,244p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  Here `ρ(x)` is the bulk number density at the brane and `ρ₀` is its asymptotic value. For each condition,
  print the `j_n` it implies through the supplied mass balance, with `ρ_br` symbolic. Print where `n` enters,
  if it enters at all.
```

```
$ sed -n 112,117p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
   - Each Part B and Part C condition is solved relative to `GM` for every far-zone `b`. The solving is case by
     case, over every branch the reduction produces, and each case carries its domain.
   - An unevaluated quantified set, `ConditionSet` or `ForAll` is a restatement, not a reduction.
   - For each reduced condition, print the `j_n` it implies through `∇·(ρ_br V) = −j_n`, with `ρ_br` symbolic and
     the sign of `V` symbolic.
   - The method is the builder's choice (item 13). If an engine cannot reduce a condition, item 12 applies.
```

### SymPy: the j_n payload pairs the condition with the balance under the replacement only
```
$ sed -n 600,602p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    mass_density = rho  # K8
    exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)
    emit('MASS_BALANCE',eq(jn,exchange))
```

```
$ sed -n 608,619p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
        domain=sp.And(substitute(gate,replacement),own_domain)
        for name,mass in condition_mass.items():
            local_mass=substitute(mass,replacement)
            condition=cas(('RADIAL_DIFFERENTIAL_IDENTITY',r,sp.Interval.open(bmin,sp.oo),domain,eq(GM,local_mass)))
            tag=prefix+'CONDITION_'+name
            emit(tag,condition)
            objects['s9b_'+tag.lower()]=condition
            if implied:
                tag=prefix+'IMPLIED_JN_'+name
                value=cas((condition,eq(jn,substitute(exchange,replacement))))
                emit(tag,value)
                objects['s9b_'+tag.lower()]=value
```

### Wolfram: the same construction
```
$ sed -n 380,389p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
massDensity = rhoBr[r]; (* K8 *)
massVector = (massDensity velocity cartesian/r) /. r -> cartRadius;
divergence = Simplify[Total[MapThread[D,{massVector,cartesian}]] /.
  {x1 -> r,x2 -> 0,x3 -> 0}, $Assumptions];
exchange = jn[r] /. First[Solve[massBalanceInput /. divRhoV -> divergence, jn[r]]];
emit["LOCAL_MASS_BALANCE", <|"Measure" -> "coordinate d3x", "Components" -> massVector,
  "Divergence" -> divergence, "Exchange" -> exchange|>];
impliedExchange[condition_, rules_] := <|"ReducedProfileCondition" -> condition,
  "ImpliedExchange" -> (jn[r] == (exchange /. rules)),
  "SignedVelocity" -> (velocity /. rules), "Density" -> (massDensity /. rules)|>;
```

### The leg's harness tallies: the j_n part is unchanged under every knife except K8
```
$ grep -c 'JN_EQ_IDENTICAL' /tmp/s9b_rev_r1_fc_runs/an/py_harness_summary.stdout
224
```

```
$ grep 'JN_EQ_DIFFERS' /tmp/s9b_rev_r1_fc_runs/an/py_harness_summary.stdout | awk '{print $2}' | sort | uniq -c
     10 K8
```

```
$ grep 'JN_EQ_IDENTICAL' /tmp/s9b_rev_r1_fc_runs/an/py_harness_summary.stdout | awk '{print $2}' | sort | uniq -c
     18 K1
     18 K10
     18 K2
     18 K3
     18 K4a
     18 K4b
     18 K4c
     18 K5
     18 K6
     18 K7
      8 K8
     18 K9a
     18 K9b
```

```
$ grep 'IMPLIED_EXCHANGE_DIFFERS' /tmp/s9b_rev_r1_fc_runs/an/wl_harness_summary.stdout | awk '{print $2}' | sort | uniq -c
     10 K8
```

```
$ grep 'IMPLIED_EXCHANGE_IDENTICAL' /tmp/s9b_rev_r1_fc_runs/an/wl_harness_summary.stdout | awk '{print $2}' | sort | uniq -c
     18 DEAD_PATH
     18 K1
     18 K10
     18 K2
     18 K3
     18 K4a
     18 K4b
     18 K4c
     18 K5
     18 K6
     18 K7
      8 K8
     18 K9a
     18 K9b
```

## Grok finding 1: Wolfram's Part C domains carry a cut that SymPy's do not
```
$ sed -n 413,420p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
responseDeltas = (Simplify[Normal[Series[#/c0 - 1,{fLocal,0,1}]],
  rho0 > 0 && c0 > 0 && Element[{n,s,fLocal},Reals]] &) /@ responseInputs;
responseRules = Table[{deltaProfile -> Function[{r}, Evaluate[responseDeltas[[i]] /. fLocal -> fProfile[r]]]}, {i,3}];
responseNames = {"CONSTANT","FIXED_RATIO","POWER"};
responseDomain[i_] := With[{resp = responseInputs[[i]],change = responseDeltas[[i]]},
  Element[resp,Reals] && resp > 0 && Abs[change] < 1 &&
    If[i == 2, localSoundSquared > 0 && (bulkSoundSquared /. rho -> rho0) > 0,
      If[i == 3,localBulkDensity/rho0 > 0,True]]] /. fLocal -> fProfile[r];
```

```
$ grep -n 'Abs' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
544:        bound=sum((sp.Abs(term) for term in sp.Add.make_args(sp.expand(envelope))),sp.S.Zero)
```

```
$ sed -n 122,123p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
9. **Part C domains (D6 B6).** Each Part C row's domain predicates are built from that row's own response. They are
   not built from another row, or from the amplitude before the response is applied.
```

## Claude finding 2: branch existence
```
$ sed -n 188,192p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
along the ray and the condition under which a ray can traverse the full flyby path and both legs of the round
trip in the required directions. Outside either condition, report the branch type (growing, decaying, absent,
or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
them.
```

```
$ sed -n 120p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
- **Optical smallness.** Define `δ(x) ≡ c_γ(x)/c₀ − 1`. The retained optical set is every monomial
```

```
$ sed -n 284p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
`c_γ² ≡ μ_⊥/ρ_br`. It relaxes under steady load. The named equations identifying its unsupplied
```

```
$ grep -n -m1 'PY_S9B_BRANCH_EXISTENCE:' /tmp/s9b_rev_r1_fc_runs/py_baseline.stdout
7:PY_S9B_BRANCH_EXISTENCE: Equality(ConditionSet(Symbol('s9b_r', positive=True), LessThan(Mul(Pow(Symbol('c_0', positive=True), Integer(2)), Pow(Add(Function('s9b_delta', **{'complex': True, 'extended_real': True, 'finite': True, 'hermitian': True, 'imaginary': False, 'infinite': False, 'real': True})(Symbol('s9b_r', positive=Tr
```

```
$ grep -c 'c0\*(1 + deltaProfile\[r\]) > 0' /tmp/s9b_rev_r1_fc_runs/wl_baseline.stdout
171
```

```
$ grep -o -m1 'c0^2\*(1 + deltaProfile\[r\])^2 > 0 && r^2\*(1 + Derivative\[1\]\[xiProfile\]\[r\]^2) > 0 && c0\*(1 + deltaProfile\[r\]) > 0' /tmp/s9b_rev_r1_fc_runs/wl_baseline.stdout
c0^2*(1 + deltaProfile[r])^2 > 0 && r^2*(1 + Derivative[1][xiProfile][r]^2) > 0 && c0*(1 + deltaProfile[r]) > 0
```

## Claude finding 3: the γ difference is printed as an uncombined difference
```
$ grep -n 'GAMMA_DIFFERENCE' research/pde_ledger_v3/directives/S9b_repair_build_directive.md | head -3
91:   | `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE` | Part B's effective `γ`s, and their difference: the deflection `γ` minus the radar `γ` |
99:     `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_<OBS>`, `B_CONDITION_<OBS>` and `B_IMPLIED_JN_<OBS>`.
104:     `B_GAMMA_<OBS>`, `B_GAMMA_DIFFERENCE`, `B_RESIDUAL_<OBS>` and `B_CONDITION_<OBS>`. `F_LIVE_` also prefixes
```

```
$ grep -n "'B_GAMMA_DIFFERENCE'" research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
327:    b = ['B_'+name+'_'+o for name in ('GAMMA','RESIDUAL','CONDITION','IMPLIED_JN') for o in obs]+['B_GAMMA_DIFFERENCE']
624:        emit(prefix+'B_GAMMA_DIFFERENCE',gated(substitute(gammas['DEFLECTION']-gammas['RADAR'],replacement),sp.And(domain,sp.Ne(GM,0))))
```

```
$ grep -n -m1 'WL_S9B_B_GAMMA_DIFFERENCE:' /tmp/s9b_rev_r1_fc_runs/wl_baseline.stdout
91:WL_S9B_B_GAMMA_DIFFERENCE: <|"OnDomain" -> ConditionalExpression[((2*GM)/(b*c0^2) < 0 && (((2*GM)/c0^3 < 0 && gammaDifference == (b*c0^5*((2*GM*Inactive[Integrate][(2*b*Derivative[1][deltaProfile][r])/Sqrt[-b^2 + r^2] + (b*velocityProfile[r]*(velocityProfile[r] - 2*r*Derivative[1][velocityProfile][r]))/(c0^2*r*Sqrt[-b^2 + r^2
```

```
$ grep -n -m1 'PY_S9B_B_GAMMA_DIFFERENCE:' /tmp/s9b_rev_r1_fc_runs/py_baseline.stdout
109:PY_S9B_B_GAMMA_DIFFERENCE: Tuple(Tuple(And(Unequality(Symbol('GM', real=True), Integer(0)), StrictLessThan(Symbol('s9b_r_turn', positive=True), Min(Pow(Add(Pow(Symbol('Z_E', positive=True), Integer(2)), Pow(Symbol('b', positive=True), Integer(2))), Rational(1, 2)), Pow(Add(Pow(Symbol('Z_R', positive=True), Integer(2)), Pow(S
```

