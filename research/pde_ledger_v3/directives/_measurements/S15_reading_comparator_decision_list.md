# Source lookups for the S15 reading comparator decision list

Commands run from the repository root. The hashes and native counts identify the compared inputs. The source excerpts define the objects and conditions used for the review.

```text
$ python3 research/pde_ledger_v3/scripts/S15_reading_input_inventory.py
sympy original bytes and sha256: 24641539 4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c
sympy END_OF_TRANSCRIPT bytes and sha256: 192065031 10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6
manifest objects with selected tags: 336
item 1 same-engine lookup: {'identical payload': 309, 'new sum with all original operand tags located': 27}
item 3 status tag text: {'Exact family object retained unevaluated: integrand contains a general field, as permitted by item 3.': 192, 'Closed-form return from the exact moment, distance and finite-interval integration construction. The branch inventory is printed separately.': 144}
wolfram original bytes and sha256: 9797390 811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6
wolfram END_OF_BASE_PASS bytes and sha256: 333550018 75622128d04da2d934c2a52468505b4438fff9f17a902d7a6cd35993d769b6c4
manifest objects with selected tags: 336
item 1 same-engine lookup: {'identical payload': 309, 'new sum with all original operand tags located': 27}
item 3 status tag text: {'"Exact family integral with general fields retained, as permitted by item 3."': 192, '"Closed-form return from the angular, distance and finite-interval integration construction. The exponent-branch and retained-function inventories are printed separately."': 144}
common object names after label spelling normalization: 336
SymPy-only names: []
Wolfram-only names: []
E_geo literal-zero spelling differs across engines: 3
E_resp literal-zero spelling differs across engines: 72
sympy changed native payloads with same tag: []
sympy sum rows: 27 original operand references: 144
wolfram changed native payloads with same tag: []
wolfram sum rows: 27 original operand references: 144
```

```text
$ sed -n 65,86p research/pde_ledger_v3/directives/S15_reading_decision_list.md
**The objects `Q`:**
- every piece and group S15's item 1 prints, at each surface, for each Cartesian component (L4);
- the region's momentum and its rate, each whole and by group, and the integrated force density (item 2; L4);
- at each surface, three named sums: the rest part (carried, whole pressure, viscous, other), the transport (centre
  and radius), and the crossing (rest plus transport); and the same three sums over both surfaces.

**For each `Q`:**
1. `Q` itself.
2. `E_Q^geo` and `E_Q^resp`, for general profiles, exact on the domain under the condition each carries, with every
   live quantity kept.
3. `E_Q^geo` on the profile family, with its dependence on `d` explicit and exact. Where the integrand holds no general
   field besides the family's amplitudes, the region's quantities and the velocities, the integrals are evaluated in
   closed form. Otherwise the object is printed with the family substituted and left unevaluated.
4. For each `Q` whose item-3 object is in closed form:
   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
   - `Q` on the family under each condition.
5. For each drain `k`, the part of `Q` that depends on that drain's velocity: `Q` minus `Q` at `Ẋ_k = 0`.
6. Dimensional homogeneity of items 2 and 5, each change `F′` carrying the dimension of `F`.

Every object is exact on S15's domain under the conditions it carries: S15's pressure conditions, and for each
```

```text
$ sed -n 137,167p research/pde_ledger_v3/directives/S15_comparator_decision_list.md
   Wolfram operand and their residual.
   - **Exact first.** Use exact symbolic subtraction where feasible, and print the exact residual as computed. Where it
     is feasible, evaluation at points is added beside it, never in its place. Otherwise evaluate at several generic
     rational points (E1; numeric-PIT practice).
   - **Test functions.** Under evaluation, every live function is replaced by test functions of several different
     forms, applied identically to both operands. That covers each field of item 3, both centres and both radii, with
     their time derivatives.
     - The test functions are not required to satisfy the supplied balance. Print the test functions, the points and
       both operands' values.
     - Distinct live functions get distinct test functions in every form and stratum. That includes `Π_o^{ij}` and
       `Π_o^{ji}`, and the two drains' profiles.
     - In every form, no test function is constant in any of its arguments, and none vanishes identically on a surface
       or in the region.
     - For a function of several arguments, the test functions tell its arguments apart: at the evaluation points,
       its derivatives in different arguments differ.
     - A derivative taken at another point is a different object (item 1), and evaluation sees it. For every compared
       value whose exact residual is not resolved, a live function's derivative taken at another point, in any one
       argument, moves the samples in at least one form and stratum. That includes a derivative inside a
       substitution, a derivative at a composite argument, and a derivative inside a surface or volume integral.
     - The geometry stays inside the printed domain at every evaluation time: `|X_1 − X_2| > R > r_in > 0`.
   - **Integrals.** An integral is evaluated after the substitution, exactly where feasible. Otherwise it is evaluated
     numerically, and the working precision and an error estimate are printed. An integral that cannot be evaluated is
     printed as unresolved coverage.
   - **Strata.** The evaluation covers the generic case, and strata where, at the evaluation time, the velocity of
     `X_2` vanishes, where both radii's velocities vanish, where both vanish, and where the velocity of `X_1` vanishes.
     In each, the higher time derivatives stay generic. A stratum the comparator cannot evaluate is printed as
     unresolved coverage.
   - **No supplied equation applied.** No supplied equation is applied to an operand to make two operands meet. That
     covers the momentum balance, the forms of `G`, `Π` and `T`, and the split of the fields. Where a pair differs by a
     supplied equation, the residual as computed is printed. The residual after that equation is applied may be
     printed beside it, naming the equation and its source line, never in its place.
```

```text
$ sed -n '343,355p;685,709p' research/pde_ledger_v3/scripts/S15_separation_departure.py
    # Reading construction: every Q was constructed above from this engine's
    # own tensor. The registry records values, never imports another engine.
    for i in range(3):
        sums = []
        for j, surf in enumerate(('inner','outer')):
            entries = surface_parts[j][i]
            rest = sum(entries[k] for k in ('carried','pressure','viscous','other'))
            transport = sum(entries[k] for k in ('transport_centre','transport_radius'))
            row = dict(rest=rest,transport=transport,crossing=rest+transport)
            sums.append(row)
            for key,value in row.items(): objects[f'{surf}/sum_{key}/component_{i+1}']=value
        for key in sums[0]:
            objects[f'both_surfaces/sum_{key}/component_{i+1}']=sums[0][key]+sums[1][key]
        solo(tag+'/Q',q)
        emit(tag+'/departure_existence_condition',existence)
        geo_derivative=variation(q,False)
        full_derivative=variation(q,True)
        geo=geo_derivative+2*q
        full=full_derivative+2*q
        response=full-geo
        solo(tag+'/separation_derivative',full_derivative)
        solo(tag+'/E_geo',geo)
        solo(tag+'/E_resp',response)
        pair(tag+'/departure_partition',full,geo+response)
        for objname,obj in (('E_geo',geo),('E_resp',response)):
            dimension_audit(tag+'/'+objname,obj)
            emit(tag+'/'+objname+'/dimension_reference_operand',sorted(dim(q)))
            emit(tag+'/'+objname+'/dimension_reference_residuals',[sub(a,b) for a in dim(obj) for b in dim(q)])
        for k in range(2):
            replacement={S.diff(c,t):S.Integer(0) for c in centres[k]}
            without=q.xreplace(replacement)
            emit(tag+f'/velocity_part_{k+1}/operand_a',q)
            emit(tag+f'/velocity_part_{k+1}/operand_b',without)
            part=q-without
            solo(tag+f'/velocity_part_{k+1}/residual',part)
            dimension_audit(tag+f'/velocity_part_{k+1}',part)
            emit(tag+f'/velocity_part_{k+1}/dimension_reference_residuals',[sub(a,b) for a in dim(part) for b in dim(q)])
        family_geo=on_family(geo)
```

```text
$ sed -n '216,225p;1505,1527p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl
Do[
 sumRows=Table[
  row=surfaceResults[[surfaceIndex,componentIndex]];
  restObject=Total[Lookup[row,channelNames]];
  transportObject=Total[Lookup[row,{"transport_centre","transport_radius"}]];
  sums=<|"rest"->restObject,"transport"->transportObject,"crossing"->restObject+transportObject|>;
  KeyValueMap[(AssociateTo[qObjects,({"inner","outer"}[[surfaceIndex]]<>"/sum_"<>#1<>"/component_"<>ToString[componentIndex])->#2])&,sums];
  sums,{surfaceIndex,2}];
 Do[AssociateTo[qObjects,("both_surfaces/sum_"<>key<>"/component_"<>ToString[componentIndex])->Total[Lookup[sumRows,key]]],{key,Keys[First[sumRows]]}],
 {componentIndex,3}];
 solo[readingTag<>"/Q",qValue];
 emit[readingTag<>"/departure_existence_condition",existenceCondition];
 geoDerivative=geometricDerivative[qValue];
 respDerivative=responseDerivative[qValue];
 totalDerivative=geoDerivative+respDerivative;
 departure=totalDerivative+2 qValue;
 geoDeparture=geoDerivative+2 qValue;
 respDeparture=departure-geoDeparture;
 solo[readingTag<>"/separation_derivative",totalDerivative];
 solo[readingTag<>"/E_geo",geoDeparture];
 solo[readingTag<>"/E_resp",respDeparture];
 pair[readingTag<>"/departure_partition",departure,geoDeparture+respDeparture];
 Do[dimensionAudit[readingTag<>"/"<>part[[1]],part[[2]]];
  emit[readingTag<>"/"<>part[[1]]<>"/dimension_reference_operand",dim[qValue]];
  emit[readingTag<>"/"<>part[[1]]<>"/dimension_reference_residuals",Flatten[Table[da-db,{da,dim[part[[2]]]},{db,dim[qValue]}],1]],
  {part,{{"E_geo",geoDeparture},{"E_resp",respDeparture}}}];
 Do[
  withoutVelocity=qValue/.Thread[D[centres[[k]],t]->ConstantArray[0,3]];
  velocityPart=qValue-withoutVelocity;
  pair[readingTag<>"/velocity_part_"<>ToString[k],qValue,withoutVelocity];
  solo[readingTag<>"/velocity_part_"<>ToString[k]<>"/object",velocityPart];
  dimensionAudit[readingTag<>"/velocity_part_"<>ToString[k],velocityPart];
  emit[readingTag<>"/velocity_part_"<>ToString[k]<>"/dimension_reference_residuals",Flatten[Table[da-db,{da,dim[velocityPart]},{db,dim[qValue]}],1]],{k,2}];
```
