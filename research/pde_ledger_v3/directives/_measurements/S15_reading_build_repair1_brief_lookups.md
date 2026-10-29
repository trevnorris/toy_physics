# Measurements: S15 reading build repair-1 brief, leak gate and citation check

Generator: `_scratch/s15/gen/s15_reading_repair1_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by the
orchestrator.

## R0. The brief

```
$ sha256sum _scratch/s15/reading_build/resume_s15_reading_repair_1.md; wc -c < _scratch/s15/reading_build/resume_s15_reading_repair_1.md; git log -1 --format='%h %s'
723c921a94b615fc42d2a5ca290864c0887ecf7cc05b63b15cba135d80d4d214  _scratch/s15/reading_build/resume_s15_reading_repair_1.md
4888
5f0fc03a Preserve the S15 reading build round 0 as reviewed (not accepted)
```

## R1. Values, exponents and methods from the review legs (fixed strings)

```
$ /usr/bin/grep -n -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= -1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= -1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
5:  - Every family integrand of a group that contains drain 1 fails that test, so each is labelled as not obtained, and no exact integration of it was attempted. These ar
count: 1
```

```
$ /usr/bin/grep -n -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
11:- V2, both engines: item 4's condition is the closed item-3 object set equal to zero, and "Q under the condition" is the family's `Q` with no restriction.
12:  - The condition is the item-3 object `== 0` (`.py:504`, `.wl:323`), not a condition reached from it. `Q_under_condition` is the family's `Q` (`.py:509`, `.wl:327`).
count: 2
```

```
$ /usr/bin/grep -n -F -- '−4' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−4' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-4' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-4' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−5' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−5' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'A_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'A_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'B_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'B_1' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'A₁' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'A₁' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'R_t' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'R_t' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
27:Then, from the final code, run each engine through the guard. Use the guard command your brief gives, with a new log directory name for each run, and one Wolfram kerne
count: 1
```

```
$ /usr/bin/grep -n -F -- 'chord' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'chord' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'azimuth' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'azimuth' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'substitut' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'substitut' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'large' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'large' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

## R2. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
37:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
1:Note from the orchestrator. Your S15 reading build was reviewed, and two defects were confirmed. The reviewed version is preserved at commit 5f0fc03a. Repair both in th
13:  - The status lines (`.py:512`, `.wl:330`) and the run record (`research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md:47`) say that t
18:    - No status line or run-record sentence states a property of a printed object that no computation in that engine produced.
29:Then update `research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md`. For each engine, it names:
count: 4
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
1:Note from the orchestrator. Your S15 reading build was reviewed, and two defects were confirmed. The reviewed version is preserved at commit 5f0fc03a. Repair both in th
6:  - Independent review reached closed forms for representative objects of each listed class by exact routes. The gaps are not limits of the objects.
count: 2
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_1.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_1.md)"
count: 0
```

## R3. Every engine line the brief cites, at the preserved commit 5f0fc03a

```
$ git show 5f0fc03a:research/pde_ledger_v3/scripts/S15_separation_departure.py | sed -n '433,434p;486p;495,499p;504p;509p;512p' | cut -c1-170
            if not angular_polynomial:
                return S.Integral(f,*expr.limits)
            status='Closed form' if not result.has(S.Integral) else 'NOT OBTAINED IN CLOSED FORM: angular integrand is outside this engine\'s exact polynomial integration
        elif result.has(d,*direction_D):
            emit(tag+'/item4/status','NOT OBTAINED EXACTLY: elimination of the universal separation and direction quantifiers is not supplied by this integration route.')
            solo(tag+'/item4/equation_before_quantifier_elimination',S.Eq(result,0,evaluate=False))
            solo(tag+'/item4/Q_on_family',family_q)
            gaps.append((name,'item4','Universal separation/direction elimination unresolved.'))
            condition=S.Eq(result,0,evaluate=False)
            qclosed=closed_route(family_q) if not has_general_field(family_q) else family_q
            emit(tag+'/item4/status','Exact implicit equation in the named parameters; its Piecewise branches include exceptional exponent values. Q is printed restricted
```

```
$ git show 5f0fc03a:research/pde_ledger_v3/mathematica/S15_separation_departure.wl | sed -n '266,269p;318,322p;323p;327p;330p' | cut -c1-170
 poly=PolynomialQ[angularInput,angularVars] && FreeQ[angularInput,theta|phi];
 answer=If[poly,
  Fold[Integrate[#1,Sequence@@{#2},Assumptions->(assumptions && Element[Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}]],Reals]),GenerateConditions->True]&,integrand,SortBy[ra
  Inactive[Integrate][integrand,ranges]];
  If[!FreeQ[familyResult,dist|directionSeparation],
   solo[readingTag<>"/item4/equation_before_quantifier_elimination",familyResult==0];
   solo[readingTag<>"/item4/Q_on_family",familyQ];
   emit[readingTag<>"/item4/status","NOT OBTAINED EXACTLY: universal separation and direction elimination is unresolved."];
   AppendTo[exactGaps,{objectName,"item4","Universal separation/direction elimination unresolved."}],
   parameterCondition=familyResult==0;
   conditionedQ=If[hasGeneral[familyQ],familyQ,angularRoute[familyQ]];
   emit[readingTag<>"/item4/status","Exact implicit equation in named parameters, retaining every Piecewise exponent branch; Q is restricted to the printed equation."]
```

```
$ git show 5f0fc03a:research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md | sed -n '47p' | cut -c120-300
d forms reached, item-4 conditions are printed as implicit equations in the named quantities with the engine’s exponent branches, alongside Q restricted to the condition.
```

```
$ git show 5f0fc03a --stat --format='%h %s' | head -1
5f0fc03a Preserve the S15 reading build round 0 as reviewed (not accepted)
```
