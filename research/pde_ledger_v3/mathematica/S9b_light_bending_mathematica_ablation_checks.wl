(* Comparator demonstrations mechanically load the final harness definitions. No engine expected-value test. *)
$HistoryLength=0;
$Messages={OutputStream["stderr",2]};
src=Import["research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl","Text"];
method=StringSplit[StringSplit[src,"(* BEGIN_DOMAIN_DIFFERENCE_METHOD"][[2]],"(* END_DOMAIN_DIFFERENCE_METHOD *)"][[1]];
method="(* BEGIN_DOMAIN_DIFFERENCE_METHOD"<>method;
Begin["S9bHarness`"];
ToExpression[Global`method];
End[];
$Assumptions=True;
show[name_,a_,b_] := Print[name,": ",ToString[S9bHarness`difference[a,b],InputForm,PageWidth->Infinity]];
show["ZERO_INTEGRAL",Inactive[Integrate][0,{r,1,Infinity}],0];
show["POLYNOMIAL_IDENTITY",(x+1)^2,x^2+2x+1];
show["CROSSING_VALUES",2x,x];
show["PROFILE_DERIVATIVE",velocityProfile'[r],0];
show["PROFILE_INTEGRAL",Inactive[Integrate][velocityProfile[r],{r,1,2}],0];
show["PROFILE_PREDICATE",velocityProfile'[r]==0,True];
show["QUANTIFIED_STATEMENT",Inactive[ForAll][r,velocityProfile[r]==0],True];
show["QUANTIFIER_EQUAL",Inactive[ForAll][r,r^2>=0],True];
show["SET_ELEMENTS",ImplicitRegion[r>1,{r}],ImplicitRegion[r>0,{r}]];
show["SET_EMPTINESS",Inactive[Exists][r,False],False];
show["ATTACHED_DOMAINS",<|"Value"->x,"Domain"->x>1|>,<|"Value"->x,"Domain"->x>0|>];
show["AGREEMENT_NOT_ZERO",unknownFunction[r],0];
$Assumptions=r>0 && c0>0 && b>bFar>0;
show["BOUND_ASSUMPTION_NOT_CAPTURED",Inactive[ForAll][r,r>0],True];
show["SET_NEGATIVE_ELEMENTS",ImplicitRegion[r< -1,{r}],ImplicitRegion[r< -2,{r}]];
show["QUANTIFIED_PROFILE_DOMAIN",<|"Value"->1,"Domain"->Inactive[ForAll][r,velocityProfile[r]==0]|>,<|"Value"->1,"Domain"->True|>];
show["ZERO_BEFORE_SEARCH",Inactive[Integrate][0,{r,1,Infinity}],0];
show["POLE_NOT_WITNESS",1/velocityProfile[r],0];

$Assumptions=True;
show["ALPHA_EQUAL_SETS",ImplicitRegion[r>0,{r}],ImplicitRegion[x>0,{x}]];
show["DISJOINT_DOMAINS",<|"Value"->x,"Domain"->x>0|>,<|"Value"->2x,"Domain"->x<0|>];
show["PROFILE_VALUE_DOMAIN",<|"Value"->velocityProfile'[r],"Domain"->r>0|>,<|"Value"->0,"Domain"->r>0|>];
show["NESTED_BINDERS",Inactive[ForAll][r,Inactive[Exists][x,x>r]],True];
show["FREE_PARAMETER_STATEMENT",Inactive[ForAll][r,r^2>=x],Inactive[ForAll][r,r^2>=x+1]];
show["AGREEING_PROFILES_NOT_ZERO",(velocityProfile[r] (velocityProfile[r]-1/3) (velocityProfile[r]-1/(5(1+r^2))) (velocityProfile[r]+1/(3(1+r^2)))),0];

show["GROUND_VALUES_NO_VARIABLES",2,1];
show["INFINITE_CHART_REJECTED",Inactive[Limit][r^(-1),r->Infinity],0];

show["PROFILE_NAME_METADATA",{velocityProfile,rhoBr,deltaProfile,xiProfile},{velocityProfile,rhoBr}];
show["PROFILE_HEADS",velocityProfile,xiProfile];
show["PURE_FUNCTIONS",Function[{q},q^2],Function[{q},1/(1+q^2)]];

show["SINGULAR_CANDIDATES_BEFORE_VALID_POINT",Inactive[Integrate][1/(r-y)^2,{r,0,1}],0];
