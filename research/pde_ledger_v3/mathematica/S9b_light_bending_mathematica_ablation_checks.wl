(* Executable syntax/structural demonstrations from the actual harness.
   No expected physics value, canonicalization or verdict. *)
$HistoryLength=0;
$Messages={OutputStream["stderr",2]};
src=Import["research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl","Text"];
syntax=Block[{$Context="S9bSyntax`",$ContextPath={"System`"}},ToExpression[src,InputForm,HoldComplete]];
Print["HARNESS_PARSE_HEAD: ",Head[syntax]];
method=StringSplit[StringSplit[src,"(* BEGIN_LITERAL_METHOD *)"][[2]],"(* END_LITERAL_METHOD *)"][[1]];
Begin["S9bHarness`"];
Block[{$Context="S9bHarness`",$ContextPath={"System`"}},ToExpression[Global`method]];
End[];
show[name_,a_HoldComplete,b_HoldComplete] := Print[name,": ",
 ToString[<|"Operands"->{a,b},"Comparison"->S9bHarness`comparison[a,b]|>,InputForm,PageWidth->Infinity]];
show["SAME_SYNTAX",HoldComplete[x+1],HoldComplete[x+1]];
show["POLYNOMIAL_FORMS",HoldComplete[(x+1)^2],HoldComplete[x^2+2 x+1]];
show["INERT_ZERO_INTEGRAL",HoldComplete[Inactive[Integrate][0,{r,1,Infinity}]],HoldComplete[0]];
show["DIFFERENT_APPROXIMATE_NUMBERS",HoldComplete[1.],HoldComplete[1.0000000000000002]];
show["DIFFERENT_PRECISIONS",HoldComplete[1.0`20],HoldComplete[1.0`30]];
show["DIFFERENT_DOMAINS",HoldComplete[ConditionalExpression[x,x>0]],HoldComplete[ConditionalExpression[x,x>1]]];
show["ALPHA_RENAMED_BINDERS",HoldComplete[Inactive[ForAll][r,r>0]],HoldComplete[Inactive[ForAll][x,x>0]]];
show["NO_EVALUATION",S9bHarness`freeze["1+1"],S9bHarness`freeze["2"]];
show["ORDERED_CASES",HoldComplete[Piecewise[{{1,x>0},{2,x>1}},3]],
 HoldComplete[Piecewise[{{2,x>1},{1,x>0}},3]]];
(* Arithmetic fixtures exercise the actual witness bookkeeping definitions.
   These are synthetic objects, not additional engine physics inputs. *)
witnessMethod=StringSplit[StringSplit[src,"(* BEGIN_WITNESS_METHOD *)"][[2]],"(* END_WITNESS_METHOD *)"][[1]];
Block[{$Context="S9bHarness`",$ContextPath={"System`"}},ToExpression[Global`witnessMethod]];
$Assumptions=True;
fixture1=<|"State"-><|"RadialMetric"->1|>|>;
fixture2=<|"State"-><|"RadialMetric"->2|>|>;
Print["COMMON_FIXED_SOLUTION_SETS: ",ToString[S9bHarness`fixedReading[{Global`fLocal},fixture1,fixture1],InputForm,PageWidth->Infinity]];
Print["DISJOINT_FIXED_SOLUTION_SETS: ",ToString[S9bHarness`fixedReading[{Global`aMetric},fixture1,fixture2],InputForm,PageWidth->Infinity]];
Print["ORDERED_CASE_SELECTION: ",ToString[S9bHarness`caseReadings[
 Piecewise[{{11,Global`n==0},{12,Global`n==1}},13],fixture1],InputForm,PageWidth->Infinity]];

(* Domain decisions and binder substitution use the real round-9 helpers. *)
Print["DOMAIN_BOOLEAN_DECISION: ",ToString[<|"Operands"->{False,True},
 "Decision"->S9bHarness`domainDifferenceDecision[False,True]|>,InputForm]];
Print["DOMAIN_UNRESOLVED_DECISION: ",ToString[<|"Operands"->{x>0,True},
 "Decision"->S9bHarness`domainDifferenceDecision[x>0,True]|>,InputForm]];
Print["BOUND_VARIABLE_SCOPE: ",ToString[S9bHarness`siteEvaluate[
 Inactive[ForAll][Global`r,Global`r>1,Global`r^2>0] &&
 Inactive[Integrate][Global`r,{Global`r,1,2}]>0,fixture1],InputForm,PageWidth->Infinity]];

Print["BOUND_ASSOCIATION_SCOPE: ",ToString[S9bHarness`siteEvaluate[
 <|"Primitive"->Inactive[Integrate][Global`r,{Global`r,1,2}],"Domain"->Global`r>1|>,fixture1],InputForm,PageWidth->Infinity]];

S9bHarness`object[data_,tag_] := ReleaseHold[data["Held"][tag]];
(* Amendment-12 fixtures print computations before any external validation.
   These synthetic objects test readings, not expected engine physics values. *)
fixtureCatalog=<|"Held"-><|"WL_S9B_B_CONDITION_DEFLECTION"->HoldComplete[
 <|"Domain"->True,"Cases"->{<|"GM"->{{Global`GM->1}},
   "VelocityReduction"-><|"Status"->"fixture status","Cases"->{
     <|"Component"->"each connected open component of squared velocity > 0","Velocity"->-Sqrt[Global`r],"Domain"->Global`r>1|>,
     <|"Component"->"each connected open component of squared velocity > 0","Velocity"->Sqrt[Global`r],"Domain"->Global`r>1|>}|>|>}|>]|>|>;
Print["SEMANTIC_BRANCH_CATALOG: ",ToString[KeyValueMap[
 Function[{id,rec},S9bHarness`publicPart[rec]],S9bHarness`catalog[fixtureCatalog]["Parts"]],InputForm,PageWidth->Infinity]];
fixtureDensity=<|"Held"-><|"WL_LOCAL_S9B_C_RESPONSES"->HoldComplete[
 <|"Responses"->{1,<|"SuppliedBulkDensityPremise"->Global`rho0>0,"DensityParts"->{
   <|"Label"->"rho > 0","LabelCondition"->Global`rho0(1+Global`fProfile[Global`r])>0,"Object"->ConditionalExpression[Global`r,Global`r>0]|>,
   <|"Label"->"rho <= 0","LabelCondition"->Global`rho0(1+Global`fProfile[Global`r])<=0,"Object"->ConditionalExpression[Global`r,Global`r>0]|>}|>,3}|>]|>|>;
Print["DENSITY_SIGN_CATALOG: ",ToString[KeyValueMap[Function[{id,rec},S9bHarness`publicPart[rec]],
 S9bHarness`catalog[fixtureDensity]["Parts"]],InputForm,PageWidth->Infinity]];
Print["TEXT_COMPARISON_DECISIONS: ",ToString[<|"Operations"->{Inactive[SameQ]["one","two"],Inactive[SameQ]["one","one"]},
 "Decisions"->{SameQ["one","two"],SameQ["one","one"]}|>,InputForm]];
Print["NEGATIVE_DENSITY_ASYMPTOTE: ",ToString[<|"Profile"->Function[{q},-9/q],
 "DensityAtPointComputation"->Inactive[ReplaceAll][1-9/Global`r,Global`r->3],
 "DensityAtPoint"->(1-9/Global`r/.Global`r->3),
 "AsymptoticComputation"->Inactive[Limit][1-9/Global`r,Global`r->Infinity],
 "AsymptoticDensity"->Limit[1-9/Global`r,Global`r->Infinity]|>,InputForm]];

(* Round-10 scoped fixtures: every computed value is printed before inspection. *)
knifeManifest=StringSplit[StringSplit[src,"knives = "][[2]],"(* BEGIN_LITERAL_METHOD *)"][[1]];
Block[{$Context="S9bHarness`",$ContextPath={"System`"}},ToExpression["knives = "<>Global`knifeManifest]];
Print["WHOLE_BRANCH_CONDITION_NEGATION: ",ToString[<|
 "ActualK10SourceReplacement"->S9bHarness`knives["K10"],
 "SyntheticWholeNegation"->Not[Inactive[ForAll][x,0<=x<=1,x>0]],
 "SyntheticPointwiseNegation"->Inactive[ForAll][x,0<=x<=1,Not[x>0]],
 "ActiveDecisions"->{Resolve[Not[ForAll[x,0<=x<=1,x>0]],Reals],
   Resolve[ForAll[x,0<=x<=1,Not[x>0]],Reals]}|>,InputForm,PageWidth->Infinity]];
fixtureProfiles=<|"Held"-><|"WL_LOCAL_S9B_ABEL_REDUCTION"->HoldComplete[
 <|"Domain"-><|"Profiles"->{Global`deltaProfile,Global`velocityProfile,Global`xiProfile,Global`rhoBr,Global`fProfile}|>|>]|>|>;
Print["LITERAL_PROFILES_READING: ",ToString[<|
 "Parts"->Values[S9bHarness`catalog[fixtureProfiles]["Parts"]],
 "ExtractedNameKeys"->S9bHarness`extractedNameKeys|>,InputForm,PageWidth->Infinity]];

(* K10's actual replacement applied to synthetic gate factors. The branch
   range differs from each constructor range, as in the engine. Print the
   operands and literal decisions; no expected physics value is supplied. *)
Block[{Global`rayConstruction,Global`traversedRadii,Global`propagationGate,
 Global`traversalGate,Global`rayGate,Global`r},
 Global`traversedRadii=(Global`r>=2 || 3<=Global`r<=4);
 Global`rayConstruction=Inactive[ForAll][Global`r,Global`r>=2,Global`r>0] && Global`fixtureEndpointEquation;
 Global`propagationGate=Global`rayConstruction &&
   Inactive[ForAll][Global`r,Global`traversedRadii,Global`fixtureBranchPredicate];
 Global`traversalGate=Global`rayConstruction &&
   Inactive[ForAll][Global`r,Global`traversedRadii,Global`fixtureTraversalPredicate];
 beforeK10={Global`rayConstruction,Global`propagationGate,Global`traversalGate};
 Block[{$Context="Global`",$ContextPath={"Global`","System`"}},
   ToExpression[S9bHarness`knives["K10"][[2]]]];
 Print["K10_GATE_FACTOR_READING: ",ToString[<|
   "OriginalFactors"->beforeK10,"ObservableGateAfterReplacement"->Global`rayGate,
   "DiagnosticFactorsAfterReplacement"->{Global`rayConstruction,Global`propagationGate,Global`traversalGate},
   "DiagnosticFactorsLiteralIdentity"->SameQ[beforeK10,{Global`rayConstruction,Global`propagationGate,Global`traversalGate}],
   "NegatedQuantifiers"->Cases[Global`rayGate,Not[q:Inactive[ForAll][___]]:>q,Infinity]|>,InputForm,PageWidth->Infinity]]];

(* Exercise K10's whole-object reading on both sides of an active domain. *)
Print["K10_CONDITIONAL_OBJECT_READING: ",ToString[S9bHarness`k10ObjectEvaluate[
 <|"OnDomain"->ConditionalExpression[False,Inactive[ForAll][x,x>0,x>1]],
   "OutsideDomain"->ConditionalExpression["NOT_ESTABLISHED",Not[Inactive[ForAll][x,x>0,x>1]]]|>,fixture1],InputForm,PageWidth->Infinity]];
