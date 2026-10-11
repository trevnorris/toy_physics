# Measurements — S9b repair build, review round 3 (generated 2026-10-10 19:19)

Generator: `_scratch/s9b_build/gen/s9b_build_r3_lookups.sh` (sha256sum/grep/sed/cat/wc, and awk for literal
counts only). `_scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/` holds Grok's scripts and stdout, copied verbatim from `/tmp/s9b_r3_grok2_*`. `_scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/py_ablation_baseline.stdout` and
`_scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/wl_ablation_BASELINE.stdout` are the two engines' harness baseline streams from Grok's run of the reviewed artifacts.

## The reviewed artifacts are the working-tree files, unchanged after both legs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r3.sha256
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
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_evaluated_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_not_established_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_executable_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_payload_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_validation.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round3_measurements.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

## Both verdicts
```
$ grep -n -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r3_claude.md _scratch/s9b_build/s9b_repair_build_review_r3_grok_cont.txt
_scratch/s9b_build/s9b_repair_build_review_r3_claude.md:8:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_build_review_r3_grok_cont.txt:1:Verdict: NEEDS REVISION
```

## The Grok leg's runs: one incomplete run, then one session that was resumed once
```
$ cat _scratch/s9b_build/build_r3_grok_incomplete.DONE _scratch/s9b_build/build_r3_grok.DONE _scratch/s9b_build/build_r3_grok_cont.DONE
EXIT=0
EXIT=0
EXIT=0
```

```
$ grep -c -i 'verdict' _scratch/s9b_build/s9b_repair_build_review_r3_grok_incomplete.txt _scratch/s9b_build/s9b_repair_build_review_r3_grok.txt
_scratch/s9b_build/s9b_repair_build_review_r3_grok_incomplete.txt:0
_scratch/s9b_build/s9b_repair_build_review_r3_grok.txt:1
```

```
$ tail -c 400 _scratch/s9b_build/s9b_repair_build_review_r3_grok.txt
nished and is the probe on record. The SymPy harness is still the pool holder, finishing the K6 engine run, with K7 through K11 still ahead.The successful probe rerun is the one already on record, at `/tmp/s9b_r3_grok2_construct/s9b_engine_probe.stdout`. The SymPy harness has finished the K6 engine (exit 0) and is differencing i
```

### The resume prompt (it carries no physics, finding or value)
```
$ cat _scratch/s9b_build/s9b_repair_build_review_r3_grok_continue_prompt.md
Continue the review task from your previous turns in this session; the original prompt still governs it.

Your previous turn ended with a progress update while the guarded SymPy harness (unit `s11c-guard-6d273f233e2f`, log
`/var/projects/toy_physics/_scratch/s11c/s9b-review-r3/grok2-sympy-harness-1`) was still running. In this mode a turn
that ends without a tool call ends the whole run: no one resumes you or answers. So keep working with tool calls until
you have produced the final report that the original prompt's "Output" section defines, and make that report your
final message.

To wait for a guarded job, use a blocking shell command that polls until the unit is no longer active, for example
`while systemctl --user is-active --quiet s11c-guard-6d273f233e2f; do sleep 60; done`, and then read its log. The
original prompt's rules on stuck runs still apply.
```

## Both accepted findings predate the r2→r3 repair (r1 baseline f4dc2d3d, r2 baseline ab19dd1d)
```
$ git show f4dc2d3d:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | grep -n 'Interval(rt,sp.oo))' | cut -c1-110
433:    real_domain = everywhere(r,existence.subs(profile_map),sp.Interval(rt,sp.oo))
434:    travel_domain = everywhere(r,traversal.subs(profile_map),sp.Interval(rt,sp.oo))
```

```
$ git show f4dc2d3d:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | grep -n "emit('A_NONRECIPROCAL_PATH_DEPENDENCE'"
490:    emit('A_NONRECIPROCAL_PATH_DEPENDENCE',gated(path['predicate'],gate))
```

```
$ git show f4dc2d3d:research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl | grep -n '"Closedness"'
152:  "Closedness" -> Simplify[And @@ Thread[Flatten[oddExterior] == 0]],
```

```
$ git show ab19dd1d:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | grep -n 'Interval(rt,sp.oo))' | cut -c1-110
537:    real_domain = everywhere(r,existence.subs(branch_map),sp.Interval(rt,sp.oo))
538:    travel_domain = everywhere(r,traversal.subs(branch_map),sp.Interval(rt,sp.oo))
610:                everywhere(r,root_domain,sp.Interval(rt,sp.oo)),
```

```
$ git show ab19dd1d:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | grep -n "emit('A_NONRECIPROCAL_PATH_DEPENDENCE'"
620:    emit('A_NONRECIPROCAL_PATH_DEPENDENCE',gated(path['predicate'],gate))
```

```
$ git show ab19dd1d:research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl | grep -n '"Closedness"'
168:  "Closedness" -> Simplify[And @@ Thread[Flatten[oddExterior] == 0]],
```

## Finding 1 (Claude F1 = Grok 1): SymPy's along-ray gates quantify over a free turning radius
### What the spec and the directive ask for
```
$ sed -n 193,198p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
along the ray and the condition under which a ray can traverse the full flyby path and both legs of the round
trip in the required directions. Outside either condition, report the branch type (growing, decaying, absent,
or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
them. State the first condition through `c_γ²` without presupposing its sign, so that each of these branch types
can follow from it.
```

```
$ sed -n 313,317p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
- **Repair round 3 (amendment 6).** The branch domain is the set where both conditions in the spec's "Branch
  existence" paragraph hold. On every stratum inside it, including `GM = 0`, each `γ` object (`B_GAMMA_DEFLECTION`,
  `B_GAMMA_RADAR`, `B_GAMMA_DIFFERENCE`), and each of its restriction and forward copies, prints what that object's
  relation determines there, computed by the engine. `NOT_ESTABLISHED` is printed only outside the branch domain
  (item 5).
```

### The constructed-domain requirement is written into the Wolfram part only
```
$ grep -n 'constructed domain' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
330:  - Every along-ray gate quantifies over a constructed domain, the radii the ray traverses. No gate quantifies
```

```
$ sed -n 330,331p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
  - Every along-ray gate quantifies over a constructed domain, the radii the ray traverses. No gate quantifies
    over an undefined head.
```

### SymPy: s9b_r_turn is a bare positive symbol; both gates quantify over [s9b_r_turn, oo)
```
$ sed -n 575p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    RE,RR,rt,bmin = sp.symbols('r_E r_R s9b_r_turn s9b_b_min',positive=True)
```

```
$ sed -n 600,601p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    real_domain = everywhere(r,existence.subs(branch_map),sp.Interval(rt,sp.oo))
    travel_domain = everywhere(r,traversal.subs(branch_map),sp.Interval(rt,sp.oo))
```

```
$ sed -n 511,517p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
def gated(value, domain, *, applicability=None):
    # NOT_ESTABLISHED belongs only to the complement of the physical branch
    # domain. Additional construction assumptions remain attached to the
    # object itself, separately from that branch classification.
    if applicability is not None:
        value=cas(('APPLICABLE_OBJECT',attached_domain(applicability),value))
    return cas(((attached_domain(domain),value),(attached_domain(sp.Not(domain)),'NOT_ESTABLISHED')))
```

### SymPy: the turning relation is built later and enters only the applicability gate
```
$ sed -n 668,676p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    turning=sp.factor(geometry['even_squared']/arr).subs(profile_map).subs(J,B/c0)
    turning_domain=sp.And(eq(turning.subs(r,rt),0),sp.Gt(radial_diff(turning,r).subs(r,rt),0),
                          everywhere(r,sp.Gt(turning,0),sp.Interval.open(rt,sp.oo)))
    observable_gate = real_domain  # K10: existence only
    branch_domain=sp.And(observable_gate,travel_domain)
    gate=sp.And(observable_gate,travel_domain,turning_domain,
                everywhere(r,root_domain,sp.Interval(rt,sp.oo)),
                sp.Lt(rt,sp.Min(RE,RR))).subs(endpoint_map).subs(B,b)
    path=nonreciprocity(r,c0,delta,velocity,xi)
```

### Wolfram: the gates quantify over the radii the constructed ray traverses
```
$ sed -n 218,232p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
radarEndpointEquation = Total[fi[radarAngleIntegrand,{r,rTurnRadar,#}] & /@ endpointRadii] == angularSeparation;
turningCondition[j_, t_] := (rayRadius /. r -> t) == j && t > bFar &&
  (D[rayRadius, r] /. r -> t) > 0;
flybyRadii = r >= rTurnFlyby;
radarRadii = rTurnRadar <= r <= Max[endpointRadii];
rayConstruction = radarEndpointEquation && turningCondition[b, rTurnFlyby] &&
  turningCondition[impactRadar, rTurnRadar] && rTurnRadar <= Min[endpointRadii] &&
  Inactive[ForAll][r, r >= rTurnFlyby, D[rayRadius, r] > 0] &&
  Inactive[ForAll][r, radarRadii, D[rayRadius, r] > 0];
traversedRadii = flybyRadii || radarRadii;
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
  physicalPropagation && speedChart && rhoBr[r] != 0 &&
  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];
traversalGate = rayConstruction && Inactive[ForAll][r, traversedRadii, physicalTraversal];
rayGate = propagationGate && traversalGate;
```

### Baseline prints: the SymPy gates' quantified interval, and the Wolfram gates' radius set
```
$ sed -n 7p _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/py_ablation_baseline.stdout | grep -o "Interval(Symbol('s9b_r_turn', positive=True), oo, false, true)" | wc -l
1
```

```
$ sed -n 8p _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/py_ablation_baseline.stdout | grep -o "Interval(Symbol('s9b_r_turn', positive=True), oo, false, true)" | wc -l
1
```

```
$ sed -n 14p _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/wl_ablation_BASELINE.stdout | grep -o 'r >= rTurnFlyby || rTurnRadar <= r <= Max\[Sqrt\[b^2 + ZE^2\], Sqrt\[b^2 + ZR^2\]\]' | wc -l
1
```

```
$ sed -n 15p _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/wl_ablation_BASELINE.stdout | grep -o 'r >= rTurnFlyby || rTurnRadar <= r <= Max\[Sqrt\[b^2 + ZE^2\], Sqrt\[b^2 + ZR^2\]\]' | wc -l
1
```

### How many of SymPy's 377 baseline tags carry s9b_r_turn
```
$ grep -c 's9b_r_turn' _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/py_ablation_baseline.stdout
325
```

```
$ wc -l < _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/py_ablation_baseline.stdout
377
```

## Finding 2 (Claude F2): the path-dependence tag pairs opposite-sense predicates, and Wolfram does not gate it
### What the spec and the directive ask for
```
$ sed -n 190,191p research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
3. **The two one-way excess times** between the same endpoints, and their **nonreciprocal part** (half the
   difference). Print whether that part depends on the path or only on the endpoints.
```

```
$ sed -n 100p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
   | `A_NONRECIPROCAL_PATH_DEPENDENCE` | item 8's object |
```

```
$ sed -n 140,141p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
8. **Path dependence (D6 B4).** Whether the nonreciprocal part depends on the path or only on the endpoints is
   computed on the nonreciprocal object Part A computed, not on a placeholder.
```

### SymPy: the predicate is curl ≠ 0, gated on the branch domain
```
$ sed -n 505,507p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    curl = sp.ImmutableMatrix([sp.simplify(radial_diff(oneform[(i+2)%3],xyz[(i+1)%3])-
                                         radial_diff(oneform[(i+1)%3],xyz[(i+2)%3])) for i in range(3)])
    return dict(oneform=oneform,curl=curl,predicate=sp.Ne(curl,sp.zeros(3,1)),
```

```
$ grep -n "emit('A_NONRECIPROCAL_PATH_DEPENDENCE'" research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
684:    emit('A_NONRECIPROCAL_PATH_DEPENDENCE',gated(path['predicate'],branch_domain,applicability=gate))
```

### Wolfram: the predicate is closedness (curl == 0), emitted without gated[]
```
$ sed -n 166,170p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
emit["A_NONRECIPROCAL_PATH_DEPENDENCE", <|"OneForm" -> oddOneForm,
  "ExteriorDerivative" -> Simplify[oddExterior],
  "Closedness" -> Simplify[And @@ Thread[Flatten[oddExterior] == 0]],
  "Domain" -> (r > 0 && physicalTraversal && physicalPropagation && speedChart),
  "RadialPrimitive" -> Inactive[Integrate][oddOneForm[[1]], {r,rAnchor,rr}]|>];
```

```
$ sed -n 237,238p research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
gated[z_, rules_:{}] := <|"OnDomain" -> ConditionalExpression[z, rayGate /. rules],
  "OutsideDomain" -> ConditionalExpression["NOT_ESTABLISHED", Not[rayGate /. rules]],
```

### Baseline prints: NOT_ESTABLISHED occurrences in the path-dependence tag
```
$ grep '^PY_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE:' _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/py_ablation_baseline.stdout | grep -o 'NOT_ESTABLISHED' | wc -l
1
```

```
$ grep '^WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE:' _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/wl_ablation_BASELINE.stdout | grep -o 'NOT_ESTABLISHED' | wc -l
0
```

```
$ grep '^WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE:' _scratch/s9b_build/s9b_repair_build_review_r3_grok_evidence/art/wl_ablation_BASELINE.stdout | grep -o '"Closedness" -> [A-Za-z]*'
"Closedness" -> True
```

## Grok 2: SymPy domain differences print as unreduced formulas
### What the directive asks for (item 14) and what K1 mutates
```
$ sed -n 222,229p research/pde_ledger_v3/directives/S9b_repair_build_directive.md
      whose difference is zero. It prints no verdict.
    - For a payload that carries a domain, the difference shows separately whether the value changed and whether
      the domain changed. The domain is the condition a payload attaches to its computed object, such as the spec's
      branch-existence conditions or the condition that a case of a reduction carries. The value is the computed
      object itself, including when that object is a relation, a predicate such as a closedness test, or a
      condition. Each part's difference is computed by the engine, so an unchanged part prints as zero. This holds
      for K11's compact differences and method residual too. K11's sample domain is coverage, not a domain in this
      sense.
```

```
$ grep -n 'K1, advection' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
234:    - **K1, advection:** in the supplied kinetic symbol, `ω − V^i k_i` is replaced by `ω`.
```

### SymPy harness: identical objects difference to zero; an undecided Boolean difference is returned as a formula
```
$ sed -n 233,254p research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
def _value_difference(a,b):
    if a == b:
        return sp.S.Zero
    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple) and len(a)==len(b):
        parts=sp.Tuple(*(_value_difference(x,y) for x,y in zip(a,b)))
        return sp.S.Zero if all(part==0 for part in parts) else parts
    if isinstance(a,sp.MatrixBase) and isinstance(b,sp.MatrixBase) and a.shape==b.shape:
        return b-a
    if isinstance(a,sp.Expr) and isinstance(b,sp.Expr):
        return b-a
    if isinstance(a,sp.logic.boolalg.Boolean) and isinstance(b,sp.logic.boolalg.Boolean):
        # Normalize arithmetic relations without evaluating functionals,
        # then compare their Boolean structure. This also removes spurious
        # changes from nonzero scalar multiples of a compatibility equation.
        result=sp.simplify_logic(sp.Xor(_normalized_predicate(a),_normalized_predicate(b)),deep=False)
        return sp.Integer(int(bool(result))) if result in (sp.true,sp.false) else result
    if isinstance(a,sp.Set) and isinstance(b,sp.Set):
        # General-profile solution sets need not have decidable membership.
        # Preserve their exact symmetric difference without asking Integral
        # membership tests to evaluate arbitrary profiles. Equality above
        # still reduces identical sets to zero.
        return sp.SymmetricDifference(a,b,evaluate=False)
```

### Grok's cited triple (K1, A_DEFLECTION_G100): s9b_V occurrences per segment, split at the role keys
Segments 1–3 are the baseline payload's slots, 4–6 the corrupted payload's, 7 the value difference, 8 the domain
difference.
```
$ sed -n 55p _scratch/s11c/s9b-review-r3/grok2-sympy-harness-1/stdout | cut -c1-60
PY_LOCAL_S9B_ABLATION_K1_PY_S9B_A_DEFLECTION_G100: Tuple(Tup
```

```
$ sed -n 55p _scratch/s11c/s9b-review-r3/grok2-sympy-harness-1/stdout | awk 'BEGIN{RS="ATTACHED_DOMAIN|VALUE_DIFFERENCE|DOMAIN_DIFFERENCE"} {n=gsub(/s9b_V/,"&"); print NR-1": length="length($0)" s9b_V="n}'
0: length=80 s9b_V=0
1: length=4508 s9b_V=2
2: length=12209 s9b_V=8
3: length=4519 s9b_V=2
4: length=2279 s9b_V=0
5: length=6580 s9b_V=0
6: length=2284 s9b_V=0
7: length=28 s9b_V=0
8: length=23756 s9b_V=12
```

### Domain differences that print as zero, of all domain differences, across Grok's SymPy harness run
```
$ grep -o "DOMAIN_DIFFERENCE'), Integer(0))" _scratch/s11c/s9b-review-r3/grok2-sympy-harness-1/stdout | wc -l
3083
```

```
$ grep -o "DOMAIN_DIFFERENCE')" _scratch/s11c/s9b-review-r3/grok2-sympy-harness-1/stdout | wc -l
4907
```

## Process: the pool registry after the Grok leg
```
$ ls _scratch/s11c/resource-pool/ | grep -E 'grok'
reservations.json.before-reconcile-grok2-construct1
reservations.json.before-reconcile-grok2-harness1
reservations.json.stale-grok-r3-1644.bak
```

```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

