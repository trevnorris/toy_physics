# Measurements: S15 reading build review, round 3 (round 4 of the build): mechanical lookups

Generator: `_scratch/s15/gen/s15_reading_r3_lookups.sh` (GNU grep, sed, awk, git, sha256sum, `LC_ALL=C`). Each block
is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes and line lengths.
Extracts are in `_scratch/s15/reading_review/r3_lookup_extracts/`. The review legs' conditions and closed forms are
not reproduced here: their files are cited by residual lines and counts only.

## M0. Artifact, reports, prompt, tree and guard state

```
$ sha256sum -c _scratch/s15/reading_review/s15_reading_build_baseline_r3.sha256
research/pde_ledger_v3/scripts/S15_separation_departure.py: OK
research/pde_ledger_v3/mathematica/S15_separation_departure.wl: OK
research/pde_ledger_v3/scripts/out/S15_separation_departure.out: OK
research/pde_ledger_v3/mathematica/out/S15_separation_departure.out: OK
research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md: OK
```

```
$ sha256sum -c _scratch/s15/reading_review/r3_reports.sha256; sha256sum _scratch/s15/reading_review/s15_reading_build_review_prompt_r3.md
/var/projects/_review_scratch/s15-reading-r3-grok/REVIEW.md: OK
/var/projects/_review_scratch/s15-reading-r3-claude/REVIEW.md: OK
f0f1bc369397c7c6cac0ba2a17c8a350244ec4e1fdd9ac84b9aa84fd9f39d95b  _scratch/s15/reading_review/s15_reading_build_review_prompt_r3.md
```

```
$ git rev-parse --short HEAD; git status --porcelain --untracked-files=all
0e51c6de
 M .claude/skills/prep-for-compact/SKILL.md
 M research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
 T research/pde_ledger_v3/scripts/out/S15_separation_departure.out
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r3_lookups.md
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; echo "active s11c-guard units: $(systemctl --user list-units --type=service --state=active --no-legend | grep -c s11c-guard)"
{}

active s11c-guard units: 0
```

```
$ cat /var/projects/_review_scratch/s15-reading-r3-grok/DONE; for d in $(ls _scratch/s11c/s15-reading-review-r3-guard); do o=_scratch/s11c/s15-reading-review-r3-guard/$d/outcome.json; if test -f $o; then echo "$d exit $(grep -o '"exitCode": *[-0-9]*' $o | grep -o '[-0-9]*$')"; else echo "$d no outcome.json"; fi; done
exit 0
claude-abl-py-A0_base exit 0
0
claude-abl-py-A1_family_u1_exponent_shift exit 0
0
claude-abl-py-A2_closedQ_times_d_over_R exit 0
0
claude-abl-py-A3_drop_special_exponent_branches exit 0
0
claude-abl-wl-A0_base no outcome.json
claude-abl-wl-A0_base-2 exit 0
0
claude-abl-wl-A1_family_u1_exponent_shift exit 0
0
claude-abl-wl-A2_closedQ_times_d_over_R exit 0
0
claude-abl-wl-A3_drop_special_exponent_branches exit 0
0
claude-c01 exit 0
0
claude-c01b exit 0
0
claude-c01c exit 0
0
claude-c02 exit 0
0
claude-c03py exit 0
0
claude-c03wl exit 0
0
claude-c04a exit 0
0
claude-c04b exit 0
0
claude-c05 exit 0
0
claude-c06 exit 0
0
claude-c07 exit 0
0
claude-d01a exit 0
0
claude-d01a2 exit 0
0
claude-d01a3 exit 0
0
claude-d01b exit 1
1
claude-d01b2 exit 1
1
claude-d01b3 exit 0
0
claude-d01c exit 1
claude-d01c2 exit 0
0
claude-d01c3 exit 0
0
claude-d01d exit 0
0
claude-d02 exit 0
0
claude-p01 exit 0
0
claude-p02 exit 0
0
claude-p03 exit 0
0
claude-rerun_py exit 0
0
claude-rerun_wl exit 0
0
grok-check-closed exit 1
grok-check-closed-2 exit 1
grok-check-num exit 1
1
grok-check-num2 exit 0
0
grok-derive-01 exit 0
0
grok-derive-02 exit 1
grok-derive-03 exit 0
0
grok-evalf exit 0
0
grok-sympy-dep exit 0
0
grok-sympy-form exit 0
0
grok-wolfram-dep exit 1
1
grok-wolfram-form no outcome.json
grok-wolfram-form2 exit 0
0
```

```
$ for d in repair3_wolfram_final_01 repair3_sympy_final_02; do echo "$d: $(grep -o '"exitCode": *[-0-9]*' _scratch/s15/guard/$d/outcome.json) $(grep -o '"wallSeconds": *[0-9.]*' _scratch/s15/guard/$d/outcome.json)"; done; stat -c '%y %n' research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/scripts/out/S15_separation_departure.out research/pde_ledger_v3/mathematica/out/S15_separation_departure.out _scratch/s15/guard/repair3_wolfram_final_01/invocation.json _scratch/s15/guard/repair3_sympy_final_02/invocation.json
repair3_wolfram_final_01: "exitCode": 0
"exitCode": 0 "wallSeconds": 653.0475128099788
"wallSeconds": 652.9804659249494
repair3_sympy_final_02: "exitCode": 0
"exitCode": 0 "wallSeconds": 5262.722074120073
"wallSeconds": 5262.656777929049
2026-10-30 08:18:09.409321032 -0600 research/pde_ledger_v3/scripts/S15_separation_departure.py
2026-10-30 06:40:44.842254395 -0600 research/pde_ledger_v3/mathematica/S15_separation_departure.wl
2026-10-30 09:46:49.784059308 -0600 research/pde_ledger_v3/scripts/out/S15_separation_departure.out
2026-10-30 06:52:28.603130439 -0600 research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
2026-10-30 06:41:36.020660712 -0600 _scratch/s15/guard/repair3_wolfram_final_01/invocation.json
2026-10-30 08:19:08.872562365 -0600 _scratch/s15/guard/repair3_sympy_final_02/invocation.json
```

## M1. The contract and X1

```
$ sed -n '76,84p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
   field besides the family's amplitudes, the region's quantities and the velocities, the integrals are evaluated in
   closed form. Otherwise the object is printed with the family substituted and left unevaluated.
4. For each `Q` whose item-3 object is in closed form:
   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
   - `Q` on the family under each condition.
5. For each drain `k`, the part of `Q` that depends on that drain's velocity: `Q` minus `Q` at `Ẋ_k = 0`.
6. Dimensional homogeneity of items 2 and 5, each change `F′` carrying the dimension of `F`.
```

```
$ sed -n '102,117p' research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r2_disposition.md
- **What must be true.**
  - For each of the 48 objects, in each engine, the conditions under which the item-3 object vanishes for every `D`
    on the domain:
    - computed exactly from the printed closed form;
    - complete, with every special exponent value treated;
    - printed as equations and inequalities in the named quantities: values of, or relations among, the exponents;
      and equations or inequalities in the amplitudes, the radii and their rates, and the velocities.
  - `Q` on the family under each condition: one per condition, computed on that condition.
  - An equation that sets the object, its coefficients, its derivatives or its jet to zero restates the question. It
    may be printed as an intermediate object, and it does not replace the conditions.
  - A not-obtained label appears only where an exact attempt that could have returned a condition did not. The label
    says what was attempted and stands beside what was printed.
  - No printed line or run-record sentence states a property that no computation in that engine produced: an
    equivalence, a bound, completeness, or "solves".
  - Each engine reaches its conditions by its own construction.
  - No object of items 1, 2, 3, 5 or 6 changes by value.
```

## M2. Wolfram item 4 for the 48 separation-dependent objects

```
$ wc -l < _scratch/s15/reading_review/r3_lookup_extracts/wl_gap48.txt; head -3 _scratch/s15/reading_review/r3_lookup_extracts/wl_gap48.txt
48
inner/carried/component_1/labels_{drain1}
inner/carried/component_1/labels_{drain1, drain2}
inner/carried/component_2/labels_{drain1}
```

```
$ echo "numbered item-4 condition lines in the transcript: $(grep -c '/item4/condition_[0-9]*: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out)"; awk '{c+=$1; e+=$2; if ($1>0) n++; if ($2>0) o++} END {print "48 objects: conditions", c, "objects with a condition", n+0, "conditions naming aa[ or bb[", e, "objects with such a condition", o+0}' _scratch/s15/reading_review/r3_lookup_extracts/wl_dep48_counts.txt
numbered item-4 condition lines in the transcript: 210
48 objects: conditions 210 objects with a condition 48 conditions naming aa[ or bb[ 0 objects with such a condition 0
```

```
$ grep '/item4/condition_[0-9]*: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | grep -c -E 'aa\[|bb\['
0
```

```
$ grep -c 'Inactive\[SolveAlways\]' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; grep -n -m1 'Inactive\[SolveAlways\]' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-220
48
1752:reading/inner/carried/component_1/labels_{drain1}/item4/identity_solver_return: Inactive[SolveAlways][Implies[dist > rout[t], (2*(2*dist^2*(Piecewise[{{-((Pi*Log[dist - rin[t]])/dist^2) + (Pi*Log[dist + rin[t]])/dis
```

```
$ grep -n -E 'SolveAlways|polynomialSlots' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-200
23:SetAttributes[trackedSolveAlways,HoldAll]; trackedSolveAlways[args___] := casStep["SolveAlways "<>ToString[HoldForm[{args}],InputForm],SolveAlways[args]];
492: coefficients,query,answer,messages,monomials,polynomialSlots,rows,monomialQuery,
522:  answer=Quiet[trackedSolveAlways[query,dist]/.SolveAlways->Inactive[SolveAlways]];messages=$MessageList];
525: emit[prefix<>"/identity_solver_scope","SolveAlways returns are candidate substitutions; no exceptional-parameter completeness is inferred from a rule list."];
530: polynomialSlots=Pick[slots,Map[!MemberQ[Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}],{rin[t],rout[t]}],#]&,physical]];
531: solo[prefix<>"/polynomial_parameters",polynomialSlots/.inverse];
535:  operands=Expand[Block[{directionDifferentiationSymbols=polynomialSlots},directionPolynomial[coefficients[[coefficientIndex]]]]];
537:  rows=List@@@CoefficientRules[operands,polynomialSlots];
541:  monomials=Join[monomials,Times@@MapThread[Power,{polynomialSlots,First[#]}]& /@ Select[rows,Last[#]=!=0&]],
546: Block[{$MessageList={}},monomialAnswer=Quiet[trackedReduce[monomialQuery,polynomialSlots,Complexes]/.Reduce->Inactive[Reduce]];messages=$MessageList];
573: emit[prefix<>"/status","NOT OBTAINED EXACTLY: complete conditions, all exceptional exponent strata, and Q under every condition were not obtained. Native SolveAlways on the separation-domain impl
574: AppendTo[exactGaps,{name,"item4","Complete parameter classification not obtained by native SolveAlways; sufficient monomial branches and their full-identity residuals are printed separately."}]
645:  (* SolveAlways is asked for the exact identity, including the separation
648:  Block[{$MessageList={}},identityReturn=Quiet[trackedSolveAlways[query,dist]];messages=$MessageList];
653:   (* SolveAlways supplies generic solutions; a returned list is not a proof
659: If[!FreeQ[answer,Reduce|SolveAlways|ForAll|Exists] || !FreeQ[answer,dist] || ListQ[answer],
660:  emit[prefix<>"/status","NOT OBTAINED EXACTLY: indicator-polynomial direction elimination was computed; SolveAlways on the separation-domain implication, or Reduce for a separation-independent co
661:  AppendTo[exactGaps,{name,"item4","Indicator-polynomial direction elimination; SolveAlways on the separation-domain implication (Reduce for separation-independent conditions) returned unresolved.
```

```
$ sed -n '515,550p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-200
 reduction=trackedPolynomialReduce[polynomial,{sphere},directionDifferentiationSymbols];
 pair[prefix<>"/direction_reconstruction",polynomial,reduction[[1]].{sphere}+reduction[[2]]];
 coefficients=(Last/@CoefficientRules[reduction[[2]],directionDifferentiationSymbols])/.Thread[coefficientSlots->values];
 solo[prefix<>"/scalar_identity_operands",coefficients/.inverse];
 query=Implies[dist>(rout[t]/.rules),And@@Thread[coefficients==0]];
 solo[prefix<>"/identity_solver_operand",query/.inverse];
 Block[{$MessageList={}},
  answer=Quiet[trackedSolveAlways[query,dist]/.SolveAlways->Inactive[SolveAlways]];messages=$MessageList];
 solo[prefix<>"/identity_solver_return",answer/.inverse];
 solo[prefix<>"/identity_solver_messages",messages];
 emit[prefix<>"/identity_solver_scope","SolveAlways returns are candidate substitutions; no exceptional-parameter completeness is inferred from a rule list."];

 (* An independent exact algebraic sufficient construction. The full analytic
    coefficients remain coefficients; only the computed physical monomials are
    annihilated. No analytic independence of the kernels is assumed. *)
 polynomialSlots=Pick[slots,Map[!MemberQ[Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}],{rin[t],rout[t]}],#]&,physical]];
 solo[prefix<>"/polynomial_parameters",polynomialSlots/.inverse];
 monomials={};
 Do[
  stepProgress["item4 polynomial monomials "<>ToString[coefficientIndex]];
  operands=Expand[Block[{directionDifferentiationSymbols=polynomialSlots},directionPolynomial[coefficients[[coefficientIndex]]]]];
  solo[prefix<>"/physical_polynomial_operand_"<>ToString[coefficientIndex],operands/.inverse];
  rows=List@@@CoefficientRules[operands,polynomialSlots];
  solo[prefix<>"/physical_polynomial_coefficients_"<>ToString[coefficientIndex],rows/.inverse];
  (* Zero coefficients are retained in the table above. Their monomials need
     no restriction; a possibly nonzero coefficient imposes its monomial. *)
  monomials=Join[monomials,Times@@MapThread[Power,{polynomialSlots,First[#]}]& /@ Select[rows,Last[#]=!=0&]],
  {coefficientIndex,Length[coefficients]}];
 monomials=DeleteDuplicates[monomials];
 monomialQuery=And@@Thread[monomials==0];
 solo[prefix<>"/sufficient_solver_operand",monomialQuery/.inverse];
 Block[{$MessageList={}},monomialAnswer=Quiet[trackedReduce[monomialQuery,polynomialSlots,Complexes]/.Reduce->Inactive[Reduce]];messages=$MessageList];
 solo[prefix<>"/sufficient_solver_return",monomialAnswer/.inverse];
 solo[prefix<>"/sufficient_solver_messages",messages];
 If[FreeQ[monomialAnswer,Reduce|Resolve|ForAll|Exists],
  branches=If[Head[LogicalExpand[monomialAnswer]]===Or,List@@LogicalExpand[monomialAnswer],{LogicalExpand[monomialAnswer]}];
```

```
$ sed -n '15649p' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-420
reading/outer/carried/component_1/labels_{drain1}/item4/status: "NOT OBTAINED EXACTLY: complete conditions, all exceptional exponent strata, and Q under every condition were not obtained. Native SolveAlways on the separation-domain implication was attempted and its exact return is printed. Polynomial-monomial Reduce and exact closed-form substitution produced the separately printed sufficient branches, where verified
```

```
$ grep -E '^reading/region/component_1/momentum_labels_\{drain1\}/item4/(condition|Q_under_condition)_[0-9]+: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-200
reading/region/component_1/momentum_labels_{drain1}/item4/condition_1: AA[1][t] == 0 && rout[t] > rin[t] > 0
reading/region/component_1/momentum_labels_{drain1}/item4/Q_under_condition_1: 0
reading/region/component_1/momentum_labels_{drain1}/item4/condition_2: BB[1][t] == 0 && rout[t] > rin[t] > 0
reading/region/component_1/momentum_labels_{drain1}/item4/Q_under_condition_2: 0
```

```
$ grep -E '^reading/outer/carried/component_1/labels_\{drain1\}/item4/(condition|Q_under_condition)_[0-9]+: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-200
reading/outer/carried/component_1/labels_{drain1}/item4/condition_1: AA[1][t] == 0 && rout[t] > rin[t] > 0
reading/outer/carried/component_1/labels_{drain1}/item4/Q_under_condition_1: 0
reading/outer/carried/component_1/labels_{drain1}/item4/condition_2: BB[1][t] == 0 && rout[t] > rin[t] > 0
reading/outer/carried/component_1/labels_{drain1}/item4/Q_under_condition_2: 0
```

## M3. SymPy item 4 for the 48 separation-dependent objects

```
$ wc -l < _scratch/s15/reading_review/r3_lookup_extracts/py_gap39.txt; wc -l < _scratch/s15/reading_review/r3_lookup_extracts/py_complete9.txt; cat _scratch/s15/reading_review/r3_lookup_extracts/py_complete9.txt
39
9
outer/carried/component_1/labels_drain1
outer/carried/component_2/labels_drain1
outer/carried/component_3/labels_drain1
outer/transport_centre/component_1/labels_drain1
outer/transport_centre/component_2/labels_drain1
outer/transport_centre/component_3/labels_drain1
outer/transport_radius/component_1/labels_drain1
outer/transport_radius/component_2/labels_drain1
outer/transport_radius/component_3/labels_drain1
```

```
$ echo "numbered item-4 condition lines in the transcript: $(grep -c '/item4/condition_[0-9]*: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out)"; awk '{c+=$1; e+=$2; if ($1>0) n++; if ($2>0) o++} END {print "48 objects: conditions", c, "objects with a condition", n+0, "conditions naming an exponent symbol", e, "objects with such a condition", o+0}' _scratch/s15/reading_review/r3_lookup_extracts/py_dep48_counts.txt; awk '$2==0 {print $3}' _scratch/s15/reading_review/r3_lookup_extracts/py_dep48_counts.txt | sed 's#/component_[123]/#/component_*/#' | sort | uniq -c
numbered item-4 condition lines in the transcript: 890
48 objects: conditions 782 objects with a condition 48 conditions naming an exponent symbol 596 objects with such a condition 18
      3 reading/inner/carried/component_*/labels_drain1+drain2
      3 reading/inner/transport_centre/component_*/labels_drain1+drain2
      3 reading/inner/transport_radius/component_*/labels_drain1+drain2
      3 reading/outer/carried/component_*/labels_drain1+drain2
      3 reading/outer/transport_centre/component_*/labels_drain1+drain2
      3 reading/outer/transport_radius/component_*/labels_drain1+drain2
      3 reading/region/component_*/momentum_labels_drain1
      3 reading/region/component_*/momentum_labels_drain1+drain2
      3 reading/region/component_*/rate_labels_drain1
      3 reading/region/component_*/rate_labels_drain1+drain2
```

```
$ grep -c "unresolved_parameter_branches: \[('all_parameters'" research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep -c "unresolved_parameter_branches: \[('branch_" research/pde_ledger_v3/scripts/out/S15_separation_departure.out
30
9
```

```
$ grep -n -E 'appellf1|def parameter_conditions_from_closed_form|solve_undetermined_coeffs' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
17:def parameter_conditions_from_closed_form(coefficients, q, distance, directions,
118:            result=call(scope+'/solve_undetermined_coeffs',S.solve_undetermined_coeffs,
155:                'solve_undetermined_coeffs of the displayed formal-coefficient identity '
222:        if atoms and len(atoms)==2 and not value.has(S.log,S.appellf1,S.polylog):
307:    if any(v.has(S.appellf1,S.polylog) for v in all_values):
886:                S.appellf1(1,-power_r,-power_affine,2,*parameters))
1089:            if (isinstance(value,S.Piecewise) or value.func in (S.appellf1,S.polylog)
1196:                result_for_coefficient=cas_call('solve_undetermined_coeffs '+str(powers),S.solve_undetermined_coeffs,solver_operand,slots,d,dict=True)
1202:        complete=(not reduced.has(d,*direction_D,S.Piecewise,S.appellf1,S.polylog))
1254:            reason='NOT OBTAINED EXACTLY: Unit-sphere monomial reduction, exact kernel abbreviation retaining all arguments, and solve_undetermined_coeffs on each scalar coefficient were attempte
```

```
$ sed -n '300,335p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
            if combined==restrictions:
                direct_attempt(transformed,scope+'/unchanged_candidate',restrictions,support);continue
            recurse(values,combined,support,scope)

    # One independent affine exponent coordinate admits an explicit complete
    # partition into the printed primitive's exceptional values and complement.
    all_values=[*coefficients,q]
    if any(v.has(S.appellf1,S.polylog) for v in all_values):
        direct_attempt(coefficients,'all_parameters');return answers,gaps
    forms=[]
    for value in all_values:
        for node in value.atoms(S.Pow):
            if node.base.has(distance) and node.exp.has(*exponent_parameters):
                constant,form=S.expand(node.exp).as_coeff_Add()
                if form not in forms:forms.append(form)
    show('exponent_forms',forms)
    if not forms:
        direct_attempt(coefficients,'all_parameters');return answers,gaps
    chosen=forms[0]
    variables=sorted(chosen.free_symbols & set(exponent_parameters),key=str)
    if not variables:
        direct_attempt(coefficients,'all_parameters');return answers,gaps
    p=S.Symbol('item4_exponent',real=True)
    coordinate=call('solve_exponent_coordinate',S.solve,S.Eq(p,chosen),variables[0],dict=True)
    show('exponent_coordinate_operand',S.Eq(p,chosen))
    if len(coordinate)!=1:
        direct_attempt(coefficients,'all_parameters');return answers,gaps
    rule=coordinate[0]
    mapped=[v.subs(rule,simultaneous=True) for v in all_values]
    remaining=set(exponent_parameters)-set(rule)
    varying_forms=[S.expand(f.subs(rule)) for f in forms]
    if any(f.has(*remaining) for f in varying_forms):
        direct_attempt(coefficients,'all_parameters');return answers,gaps
    original_parameters=parameters
    parameters=tuple([p]+[v for v in parameters if v not in rule])
    exponent_parameters=tuple([p]+[v for v in exponent_parameters if v not in rule])
```

```
$ grep -n -E '^reading/region/component_1/momentum_labels_drain1/item4/(condition|Q_under_condition)_[0-9]+: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-200
111831:reading/region/component_1/momentum_labels_drain1/item4/condition_0: Eq(A_1(t), 0) & (R(t) > r_in(t))
111833:reading/region/component_1/momentum_labels_drain1/item4/Q_under_condition_0: 0
111835:reading/region/component_1/momentum_labels_drain1/item4/condition_1: Eq(B_1(t), 0) & (R(t) > r_in(t))
111837:reading/region/component_1/momentum_labels_drain1/item4/Q_under_condition_1: 0
```

```
$ sed -n '111743p' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-260
reading/region/component_1/momentum_labels_drain1/item4/all_parameters/solve_undetermined_coeffs/return: [{A_1(t): 0}, {B_1(t): 0}]
```

```
$ grep -E '^reading/outer/carried/component_1/labels_drain1/item4/(condition|Q_under_condition)_[0-9]+: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-200
reading/outer/carried/component_1/labels_drain1/item4/condition_0: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Ne(2*a_1 + b_1, -6) & Ne(2*a_1 + b_1, -4) & Ne(2*a_1 + b_1, -2) & Ne(2*a_1 + b_1, 0)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_0: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_1: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Ne(2*a_1 + b_1, -6) & Ne(2*a_1 + b_1, -4) & Ne(2*a_1 + b_1, -2) & Ne(2*a_1 + b_1, 0)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_1: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_2: (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -1)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_2: 4*pi*N_1*A_1(t)**2*B_1(t)*R(t)**3/(3*d**2)
reading/outer/carried/component_1/labels_drain1/item4/condition_3: (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -2)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_3: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_4: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -4)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_4: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_5: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -4)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_5: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_6: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -6)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_6: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_7: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -6)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_7: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_8: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, 0)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_8: 0
reading/outer/carried/component_1/labels_drain1/item4/condition_9: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, 0)
reading/outer/carried/component_1/labels_drain1/item4/Q_under_condition_9: 0
```

```
$ grep -E '^reading/inner/carried/component_1/labels_drain1/item4/(condition_[0-9]+|unresolved_parameter_branches): ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-260
reading/inner/carried/component_1/labels_drain1/item4/condition_0: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Ne(2*a_1 + b_1, -6) & Ne(2*a_1 + b_1, -4) & Ne(2*a_1 + b_1, -2) & Ne(2*a_1 + b_1, 0)
reading/inner/carried/component_1/labels_drain1/item4/condition_1: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Ne(2*a_1 + b_1, -6) & Ne(2*a_1 + b_1, -4) & Ne(2*a_1 + b_1, -2) & Ne(2*a_1 + b_1, 0)
reading/inner/carried/component_1/labels_drain1/item4/condition_2: (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -1)
reading/inner/carried/component_1/labels_drain1/item4/condition_3: (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -2)
reading/inner/carried/component_1/labels_drain1/item4/condition_4: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -4)
reading/inner/carried/component_1/labels_drain1/item4/condition_5: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -4)
reading/inner/carried/component_1/labels_drain1/item4/condition_6: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -6)
reading/inner/carried/component_1/labels_drain1/item4/condition_7: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, -6)
reading/inner/carried/component_1/labels_drain1/item4/condition_8: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, 0)
reading/inner/carried/component_1/labels_drain1/item4/condition_9: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, 0)
reading/inner/carried/component_1/labels_drain1/item4/condition_10: Eq(B_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, 0)
reading/inner/carried/component_1/labels_drain1/item4/condition_11: Eq(A_1(t), 0) & (R(t) > r_in(t)) & Eq(2*a_1 + b_1, 0)
reading/inner/carried/component_1/labels_drain1/item4/unresolved_parameter_branches: [('branch_4/node_13/unclassified_factors', 'NOT OBTAINED EXACTLY: instantaneous-polynomial expansion and solve_undetermined_coeffs of the displayed formal-coefficient identity
```

## M4. The legs' derivation evidence, by residual and count

Claude, route-1 closed forms (`derive/d01b_objects_v2.stdout`), the rate (`derive/d01c2_rate.stdout`), the series
route (`derive/d01c3_rate_series.stdout`), the engine comparison (`compare/c04_engine_Q_vs_mine.stdout`) and the
Wolfram probes (`probe/`):

```
$ grep -c 'residual components' /var/projects/_review_scratch/s15-reading-r3-claude/derive/d01b_objects_v2.stdout; grep -o 'residual components \[[0, ]*\]' /var/projects/_review_scratch/s15-reading-r3-claude/derive/d01b_objects_v2.stdout | sort | uniq -c
19
     17 residual components [0, 0, 0]
      1 residual components [0]
```

```
$ grep -c 'E_1 residual = ' /var/projects/_review_scratch/s15-reading-r3-claude/derive/d01c2_rate.stdout; grep -c 'E_1 residual = 0 ;' /var/projects/_review_scratch/s15-reading-r3-claude/derive/d01c2_rate.stdout
15
11
```

```
$ grep -c 'structure residual at D = X(1/3,2/3,2/3): 0$' /var/projects/_review_scratch/s15-reading-r3-claude/derive/d01c3_rate_series.stdout
6
```

```
$ grep -c 'difference = ' /var/projects/_review_scratch/s15-reading-r3-claude/compare/c04_engine_Q_vs_mine.stdout; grep -c 'difference = 0$' /var/projects/_review_scratch/s15-reading-r3-claude/compare/c04_engine_Q_vs_mine.stdout
12
12
```

```
$ cat /var/projects/_review_scratch/s15-reading-r3-claude/probe/p01_solvealways.stdout
A  SolveAlways[Implies[x>1, a x^2 + b x == 0], x] -> SolveAlways[Implies[x > 1, b*x + a*x^2 == 0], x]
B  SolveAlways[a x^2 + b x == 0, x] -> {{a -> 0, b -> 0}}
C  SolveAlways[Implies[x>1, {a x^2 == 0, b x == 0}], x] -> SolveAlways[Implies[x > 1, {a*x^2 == 0, b*x == 0}], x]
D  SolveAlways[Implies[x>1, a (x+1)^(p+2) - a (x-1)^(p+2) == 0], x] -> SolveAlways[Implies[x > 1, a*(-(-1 + x)^(2 + p) + (1 + x)^(2 + p)) == 0], x]
E  SolveAlways[a ((x+1)^(p+2) - (x-1)^(p+2)) == 0, x] -> SolveAlways[a*(-(-1 + x)^(2 + p) + (1 + x)^(2 + p)) == 0, x]
F  SolveAlways[Implies[x>1, Piecewise[{{a Log[x], p == 0}}, a x^p] == 0], x] -> SolveAlways[Implies[x > 1, Piecewise[{{a*Log[x], p == 0}}, a*x^p] == 0], x]
```

```
$ grep -h 'numeric control' /var/projects/_review_scratch/s15-reading-r3-claude/probe/p02_wl_own_closed_form.stdout /var/projects/_review_scratch/s15-reading-r3-claude/probe/p03_wl_region_momentum.stdout | cut -c1-150; grep -h -c '| numeric -> 0' /var/projects/_review_scratch/s15-reading-r3-claude/probe/p02_wl_own_closed_form.stdout /var/projects/_review_scratch/s15-reading-r3-claude/probe/p03_wl_region_momentum.stdout
reading/outer/carried/component_1/labels_{drain1}/item3/object | numeric control at aa[1]=3/10 -> 1.0765625333594511894152567576417244528`20.
reading/outer/carried/component_1/labels_{drain1}/item3/object | numeric control at aa[1]=3/10 -> 1.0765625333594511894152567576417244528`20.
reading/outer/transport_radius/component_1/labels_{drain1}/item3/object | numeric control at aa[1]=3/10 -> -1.52800392654019863072101050029568338659`2
reading/region/component_1/momentum_labels_{drain1}/item3/object | numeric control at aa[1]=3/10 -> 2.`20.*(-1.05`20.*(1.85049620396310738494116148747
reading/region/component_1/momentum_labels_{drain1}/item3/object | numeric control at aa[1]=3/10 -> -1.33205274984102761993749673232471396363`20.
3
1
```

Grok (`exact_conditions.stdout`: the tags whose value is exactly `0`; `evalf_points.stdout` in full; `inventory_48.txt` head):

```
$ grep -E '^[A-Z0-9_]+: 0$' /var/projects/_review_scratch/s15-reading-r3-grok/exact_conditions.stdout
CARRIED_Q_LIMIT: 0
CARRIED_E_LIMIT: 0
CARRIED_Q_LIMIT_SERIES: 0
CARRIED_E_LIMIT_SERIES: 0
RADIUS_E_ON: 0
RADIUS_E_ON_SERIES: 0
RADIUS_Q_MINUS_PRED: 0
MOM_K_MINUS_PRED: 0
MOM_E_ON: 0
RATE_ON_EAMP_GEOM: 0
RATE_ON_ER_GEOM: 0
RATE_ON_ERIN_GEOM: 0
```

```
$ cat /var/projects/_review_scratch/s15-reading-r3-grok/evalf_points.stdout
START: evalf_points
ON_d5_Q_EVALF: -1.1728612573401894757
ON_d5_Q_MINUS_PRED_EVALF: 0E-130
ON_d5_E_EVALF: 0E-129
ON_d8_Q_EVALF: -0.45814892864851151394
ON_d8_Q_MINUS_PRED_EVALF: 0E-129
ON_d8_E_EVALF: -0E-154
OFF_d5_E_EVALF: -3036.0351404291761857
RATE_ON_V0_EVALF: 0E-129
RATE_ON_V1_EVALF: -0.46914450293607579028
RATE_OFF_V0_EVALF: -3036.0351404291761857
PRED_d5: -1.1728612573401894757
PRED_d8: -0.45814892864851151394
DONE: evalf_points
```

```
$ head -3 /var/projects/_review_scratch/s15-reading-r3-grok/inventory_48.txt
SYMPY objects 48
with an exponent symbol in some condition 18
with no exponent symbol in any condition 30
```

## M5. The run record's claims

```
$ grep -n -i -E 'complete|establish|solve|equivalen|bound|sufficient|not obtained|missing' research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md | cut -c1-240
5:This record supersedes the previous record's claims that finite-identity equations solved all requested conditions and supplied Q under every condition. Those claims were unsupported. The repaired engines print parameter branches obtained
17:- The 336 item-4 statuses and transcript end marker are present. There are 210 numbered sufficient conditions and 210 corresponding Q restrictions across the 48 separation-dependent objects; none of these entries is `Missing`.
19:Missing requested objects: complete item-4 classifications, including all exceptional exponent strata, and Q under every condition, for all 48 separation-dependent objects. These are the three Cartesian components of each drain-1 and mix
21:Attempt: direction-polynomial reduction of the printed closed form, native `SolveAlways` on the separation-domain implication, and `Reduce` of computed physical-monomial equations. Native identity-solver returns remain candidates or unev
34:- Nine separation-dependent objects have no unresolved branch entry: the three Cartesian components of the outer surface's drain-1 carried current, centre transport and radius transport. Their primitive-support partitions, necessary fact
36:Missing requested objects: complete item-4 classifications and Q under every condition for the following **39** objects (each row includes all three Cartesian components):
44:The exact unresolved parameter branches are indexed beside each object's `/item4/unresolved_parameter_branches`, with the attempted construction and its return. Final `EXACT_OBJECT_GAP` lines index all 39 objects. Conditions and Q are pr
46:Attempt: a computed affine exponent coordinate, partition by the closed primitive's exceptional supports, differential determinant equations used only as necessary candidate equations, exact factor zero sets, and substitution into the fu
48:The first full SymPy attempt (`repair3_sympy_final_01`) was stopped by the author through its live guard after a prolonged silent direct coefficient-matching call. Its unit became inactive and its reservation was released. The fallback w
52:The old finite-jet condition and reference-point Euler template were removed from item 4 in both engines. No finite multiplicity bound or equivalence is claimed. Intermediate necessary equations are labelled as intermediate candidates ra
58:SymPy: **18602** of **18608** retained lines match byte-for-byte. The six differing lines are `reading/outer/pressure/component_i/whole_minus_expansion/E_geo` and `/E_resp`, for `i = 1, 2, 3`. Exact parsing and bound-variable canonicaliz
```

## M6. Unchanged objects and reruns (the legs' summaries)

```
$ cat /var/projects/_review_scratch/s15-reading-r3-claude/compare/c02_six_lines.stdout | cut -c1-160
reading/outer/pressure/component_1/whole_minus_expansion/E_geo | len r2/r3 = 3178 3177 | expand(r2 - r3) = 0
reading/outer/pressure/component_1/whole_minus_expansion/E_resp | len r2/r3 = 5794 5793 | expand(r2 - r3) = 0
reading/outer/pressure/component_2/whole_minus_expansion/E_geo | len r2/r3 = 3178 3177 | expand(r2 - r3) = 0
reading/outer/pressure/component_2/whole_minus_expansion/E_resp | len r2/r3 = 5794 5793 | expand(r2 - r3) = 0
reading/outer/pressure/component_3/whole_minus_expansion/E_geo | len r2/r3 = 3142 3141 | expand(r2 - r3) = 0
reading/outer/pressure/component_3/whole_minus_expansion/E_resp | len r2/r3 = 5758 5757 | expand(r2 - r3) = 0
```

```
$ cat /var/projects/_review_scratch/s15-reading-r3-claude/compare/c03_rerun_wl.stdout /var/projects/_review_scratch/s15-reading-r3-claude/compare/c05_rerun_conditions.stdout | cut -c1-200
{'identical': 33736}
differing item4 tags: 0
objects with item-4 conditions: production 336 rerun 336 ; same object set: True
objects identical in bytes: 336 ; identical after condition-set canonicalisation: 0 ; differing: 0
```

```
$ grep -E 'aligned_steps' /var/projects/_review_scratch/s15-reading-r3-grok/preservation_compare.txt
sympy aligned_steps 18587 value_or_key_diffs 6
wolfram aligned_steps 20438 value_or_key_diffs 0
```

## M7. Ablations

```
$ grep -c '^=====' /var/projects/_review_scratch/s15-reading-r3-claude/ablate/wl_ablation_diffs.txt; grep -E '^=====|^[<>] ' /var/projects/_review_scratch/s15-reading-r3-claude/ablate/wl_ablation_diffs.txt | sed 's/: .*//' | cut -c1-120
3
===== diff A0_base vs A1_family_u1_exponent_shift (Wolfram)
< outer/carried/component_1/labels_{drain1}/item3/exponent_branch_inventory
> outer/carried/component_1/labels_{drain1}/item3/exponent_branch_inventory
< outer/transport_radius/component_1/labels_{drain1}/item3/exponent_branch_inventory
> outer/transport_radius/component_1/labels_{drain1}/item3/exponent_branch_inventory
< region/component_1/momentum_labels_{drain1}/item3/exponent_branch_inventory
> region/component_1/momentum_labels_{drain1}/item3/exponent_branch_inventory
===== diff A0_base vs A2_closedQ_times_d_over_R (Wolfram)
===== diff A0_base vs A3_drop_special_exponent_branches (Wolfram)
< outer/carried/component_1/labels_{drain1}/item3/exponent_branch_inventory
> outer/carried/component_1/labels_{drain1}/item3/exponent_branch_inventory
< outer/transport_radius/component_1/labels_{drain1}/item3/exponent_branch_inventory
> outer/transport_radius/component_1/labels_{drain1}/item3/exponent_branch_inventory
< region/component_1/momentum_labels_{drain1}/item3/exponent_branch_inventory
> region/component_1/momentum_labels_{drain1}/item3/exponent_branch_inventory
```

```
$ grep -c '/item4/condition' /var/projects/_review_scratch/s15-reading-r3-claude/ablate/py_ablation_diffs_A1_A3.txt; grep -c '/item4/condition' /var/projects/_review_scratch/s15-reading-r3-claude/ablate/py_ablation_diffs_A2.txt
83
52
```

```
$ cat /var/projects/_review_scratch/s15-reading-r3-grok/ablation_wl_dep_diff.txt | cut -c1-160
DIFF item4/condition_1
  BASE AA[1][t] == 0 && rout[t] > rin[t] > 0
  DEP  Missing["NotObtained"]
DIFF item4/condition_2
  BASE BB[1][t] == 0 && rout[t] > rin[t] > 0
  DEP  
DIFF item4/condition_2/comparison
  BASE "There is nothing to compare it with in this engine."
  DEP  
```
