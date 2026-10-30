# Measurements: S15 reading build repair-2 brief, leak gate and citation check

Generator: `_scratch/s15/gen/s15_reading_repair2_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by the
orchestrator.

## R0. The brief

```
$ sha256sum _scratch/s15/reading_build/resume_s15_reading_repair_2.md; wc -c < _scratch/s15/reading_build/resume_s15_reading_repair_2.md; git log -1 --format='%h %s'
a164ee597746479933557bd0eae7a47236dfd7f32e512830bcdf894dde41868c  _scratch/s15/reading_build/resume_s15_reading_repair_2.md
5471
6c0fc66b Preserve the S15 reading build repair 1 as reviewed in round 1 (not accepted)
```

## R1. Values, exponents and methods from the review legs (fixed strings)

```
$ /usr/bin/grep -n -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= -1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= -1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
16:    - Every item-3 closed form each engine prints is exact at every real exponent value on the domain, including those where a printed divisor vanishes, or carries a n
21:- Every object of items 1, 2, 5 and 6 stays as it is, at the generality it has now, and no item-3 object changes except as W2 requires. Make no change beyond what this
count: 2
```

```
$ /usr/bin/grep -n -F -- '−4' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−4' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-4' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-4' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−5' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−5' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'A_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'A_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'B_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'B_1' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'A₁' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'A₁' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'R_t' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'R_t' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
14:  - Their status line says "including logarithmic branches". It is the same when the exceptional branch of the primitive is removed, so it is typed text, not derived f
27:Run every engine job through the guard, with the command your brief gives and a new log directory name for each run, one Wolfram kernel at a time. If the guard refuses
count: 2
```

```
$ /usr/bin/grep -n -F -- 'chord' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'chord' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'azimuth' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'azimuth' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'substitut' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'substitut' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
22:- The rules of your original brief still hold: printed objects and no conclusions; operands and residuals before any guard; names that name objects only; no emission t
count: 1
```

```
$ /usr/bin/grep -n -F -- 'large' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'large' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

## R2. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
37:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
1:Note from the orchestrator. Your repair was reviewed, and two defects were confirmed. The reviewed version is preserved at commit 6c0fc66b. Repair both in this session,
29:Finish with the final guarded run of each engine, from the final code. Then update `research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record
count: 2
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
1:Note from the orchestrator. Your repair was reviewed, and two defects were confirmed. The reviewed version is preserved at commit 6c0fc66b. Repair both in this session,
6:  - Independent review reached the conditions, and `Q` under them, for representative objects of these classes by exact routes, and checked them against a second route.
count: 2
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_2.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_2.md)"
count: 0
```

## R3. Every engine and transcript line the brief cites, at the preserved commit 6c0fc66b

```
$ git show 6c0fc66b:research/pde_ledger_v3/scripts/S15_separation_departure.py | sed -n '776p;782,788p;820p' | cut -c1-170
            solver_operand,definitions=abbreviate_identity(coefficient,slots) if coefficient.has(d) else (coefficient,{})
                result_for_coefficient=cas_call('solve_undetermined_coeffs '+str(powers),S.solve_undetermined_coeffs,solver_operand,slots,d,dict=True)
            except (NotImplementedError,S.PolynomialError) as exc:
                result_for_coefficient=type(exc).__name__+': '+str(exc)
            attempted.append((powers,result_for_coefficient))
            emit(tag+'/scalar_identity_solver_return',(powers,result_for_coefficient))
        emit(tag+'/identity_solver_return',attempted)
        complete=(not reduced.has(d,*direction_D,S.Piecewise,S.appellf1,S.polylog))
        else: solutions=[]
```

```
$ git show 6c0fc66b:research/pde_ledger_v3/mathematica/S15_separation_departure.wl | sed -n '492p;499p;501,507p;510p' | cut -c1-170
 query=Implies[dist>rout[t]/.rules,directionAnswer];
  Block[{$MessageList={}},identityReturn=Quiet[trackedSolveAlways[query,dist]];messages=$MessageList];
  If[MatchQ[identityReturn,{{(_Rule|_RuleDelayed)...}...}] && messages==={},
   identityConditions=Or@@Map[And@@(Equal@@#& /@ (List@@@#))&,identityReturn];
   solo[prefix<>"/identity_candidate_conditions",identityConditions];
   (* SolveAlways supplies generic solutions; a returned list is not a proof
      of completeness on all exceptional parameter strata. Keep that gap. *)
   answer=identityReturn,
   answer=identityReturn]
 If[!FreeQ[answer,Reduce|SolveAlways|ForAll|Exists] || !FreeQ[answer,dist] || ListQ[answer],
```

```
$ sed -n '21458,21463p;22853,22858p;24223,24228p' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-110
reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -(dist*directionSeparation[1]*Piecewis
reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand/comparison: "There is nothing to compar
WORKING_ON: "reading/region/component_1/rate_labels_{drain1} / closed geometric differentiation"
reading/region/component_1/rate_labels_{drain1}/item3/object: 2*(-(dist*directionSeparation[1]*Piecewise[{{(2*
reading/region/component_1/rate_labels_{drain1}/item3/object/comparison: "There is nothing to compare it with 
reading/region/component_1/rate_labels_{drain1}/item3/status: "Closed form from generating-polynomial angular 
reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand: -(dist*directionSeparation[2]*Piecewis
reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand/comparison: "There is nothing to compar
WORKING_ON: "reading/region/component_2/rate_labels_{drain1} / closed geometric differentiation"
reading/region/component_2/rate_labels_{drain1}/item3/object: 2*(-(dist*directionSeparation[2]*Piecewise[{{(2*
reading/region/component_2/rate_labels_{drain1}/item3/object/comparison: "There is nothing to compare it with 
reading/region/component_2/rate_labels_{drain1}/item3/status: "Closed form from generating-polynomial angular 
reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand: -(dist*directionSeparation[3]*Piecewis
reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand/comparison: "There is nothing to compar
WORKING_ON: "reading/region/component_3/rate_labels_{drain1} / closed geometric differentiation"
reading/region/component_3/rate_labels_{drain1}/item3/object: 2*(-(dist*directionSeparation[3]*Piecewise[{{(2*
reading/region/component_3/rate_labels_{drain1}/item3/object/comparison: "There is nothing to compare it with 
reading/region/component_3/rate_labels_{drain1}/item3/status: "Closed form from generating-polynomial angular 
```

```
$ for f in '(2 + aa[1] + bb[1])' '(4 + aa[1] + bb[1])'; do echo "$f in line 21461: $(sed -n '21461p' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | grep -o -F "$f" | wc -l)"; done
(2 + aa[1] + bb[1]) in line 21461: 20
(4 + aa[1] + bb[1]) in line 21461: 10
```

```
$ git show 6c0fc66b --stat --format='%h %s' | head -1; git diff --quiet 6c0fc66b -- research/pde_ledger_v3/mathematica/out/S15_separation_departure.out && echo 'Wolfram transcript unchanged since 6c0fc66b'
6c0fc66b Preserve the S15 reading build repair 1 as reviewed in round 1 (not accepted)
Wolfram transcript unchanged since 6c0fc66b
```
