(* Amendment 11. One kernel; unchanged engine baseline; one source
   replacement per knife. Literal frozen emitted syntax is the only tag
   comparison. Witness searches are separate, one per knife. *)
Begin["S9bHarness`"];
$HistoryLength = 0;
root = Directory[];
engine = FileNameJoin[{root,"research","pde_ledger_v3","mathematica","S9b_light_bending_mathematica_audit.wl"}];
work = FileNameJoin[{root,"_scratch","s9b_wl_repair8","ablation-final"}];
If[!DirectoryQ[work],CreateDirectory[work,CreateIntermediateDirectories->True]];
source = Import[engine,"Text"];
S9bHarness`progress = OpenWrite[FileNameJoin[{work,"progress.log"}]];
note[stage_,name_] := (WriteString[S9bHarness`progress,
  ToString[{stage,name,MemoryInUse[],MaxMemoryUsed[]},InputForm]<>"\n"];Flush[S9bHarness`progress]);
print[name_,value_] := (WriteString[First[$Output],"WL_LOCAL_S9B_ABLATION_"<>name<>": "<>
  ToString[value,InputForm,PageWidth->Infinity]<>"\n"];Flush[First[$Output]]);
knives = <|
 "K1" -> {"kineticSymbol = (omega - vVector.kCov)^2;", "kineticSymbol = omega^2;"},
 "K2" -> {"metric2 = DiagonalMatrix[{aMetric, r^2}];", "metric2 = DiagonalMatrix[{1, r^2}];"},
 "K3" -> {"dispersionSpeed = cSquared;", "dispersionSpeed = c0^2;"},
 "K4a" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] /. HoldPattern[Derivative[j_][deltaProfile][arg_]] :> 0;"},
 "K4b" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] /. HoldPattern[Derivative[j_][velocityProfile][arg_]] :> 0;"},
 "K4c" -> {"xiSlope = Simplify[D[xiAnsatz, r], $Assumptions];", "xiSlope = Simplify[D[xiConstant, r], $Assumptions];"},
 "K5" -> {"returnDirection = -1;", "returnDirection = 1;"},
 "K6" -> {"radarSource = roundTrip;", "radarSource = oneWayER;"},
 "K7" -> {"quantifierMode = \"EVERY_B\";", "quantifierMode = \"FIXED_B\";"},
 "K8" -> {"massDensity = rhoBr[r];", "massDensity = rhoBrConstant;"},
 "K9a" -> {"localSoundSquared = bulkSoundSquared /. rho -> rho0 (1 + fLocal);", "localSoundSquared = bulkSoundSquared /. rho -> rho0;"},
 "K9b" -> {"localBulkDensity = rho0 (1 + fLocal);", "localBulkDensity = rho0;"},
 "K10" -> {"branchGateLocal = propagatingLocal;", "branchGateLocal = Not[propagatingLocal];"},
 "K11" -> {"vVector = {vLocal, 0};", "vVector = {vLocal, vLocal/r};"}|>;
(* BEGIN_LITERAL_METHOD *)
(* ToExpression's third argument wraps the parsed emitted syntax before
   evaluation. BinarySerialize records its literal tree, including exact
   binary approximate numbers. SameQ is applied to integer byte lists,
   never to approximate operands. No normalization or arithmetic precedes
   a mark. The original emitted text is retained for every printed operand. *)
freeze[text_String] := Block[{$Context="Global`",$ContextPath={"Global`","System`"}},
  ToExpression[text,InputForm,HoldComplete]];
literalBytes[held_HoldComplete] := Normal[BinarySerialize[held]];
literalMark[a_HoldComplete,b_HoldComplete] := If[
  SameQ[literalBytes[a],literalBytes[b]],"identical","not compared"];
numberEncodings[held_HoldComplete] := Cases[held,
  q_Real :> <|"Number"->q,"LiteralBytes"->Normal[BinarySerialize[HoldComplete[q]]]|>,Infinity];
comparison[a_HoldComplete,b_HoldComplete] := <|
  "Mark"->literalMark[a,b],
  "Comparison"->"SameQ on BinarySerialize integer bytes of HoldComplete parsed emitted InputForm; no evaluation",
  "ApproximateNumberEncodings"->{numberEncodings[a],numberEncodings[b]}|>;
(* END_LITERAL_METHOD *)

run[path_,name_] := Module[{capture,text,lines,payloads,held,state},
  ClearSystemCache[];note["engine begin",name];
  capture=OpenWrite[FileNameJoin[{work,name<>".stdout"}]];
  Block[{$Context="Global`",$ContextPath={"Global`","System`"},$Output={capture}},Get[path]];
  Close[capture];
  text=Import[FileNameJoin[{work,name<>".stdout"}],"Text"];
  lines=StringSplit[text,"\n"];
  payloads=Association[Map[Function[line,With[{pos=First[First[StringPosition[line,": "]]]},
    StringTake[line,pos-1]->StringDrop[line,pos+1]]],lines]];
  If[Length[payloads]!=Length[lines] || !And@@(StringStartsQ[#,"WL_"]&/@Keys[payloads]),Quit[95]];
  held=AssociationMap[freeze[payloads[#]]&,Keys[payloads]];
  state=<|"Assumptions"->$Assumptions,"LocalAssumptions"->Global`localAssumptions,
    "Dispersion"->Global`dispersion,"Roots"->Global`frequencyRoots,
    "SelectedRoot"->Global`positiveRoot,"Subcritical"->Global`subcritical,
    "PhysicalPropagation"->Global`physicalPropagation,"PhysicalTraversal"->Global`physicalTraversal,
    "SpeedChart"->Global`speedChart,"RadialMetric"->Global`atOne[Global`radialMetric],
    "RayConstruction"->Global`rayConstruction|>;
  note["engine end",name];<|"Text"->payloads,"Held"->held,"State"->state|>];
object[data_,tag_] := ReleaseHold[data["Held"][tag]];
mutate[name_] := Module[{pair=knives[name],target},
  print[name<>"_SITE",<|"Construction"->First[pair],"Mutation"->Last[pair],
    "Occurrences"->StringCount[source,First[pair]],"SourceLine"->sourceLine[First[pair]]|>];
  If[StringCount[source,First[pair]]!=1,Quit[92]];
  target=FileNameJoin[{work,name<>".wl"}];Export[target,StringReplace[source,Rule@@pair],"Text"];target];
sourceLine[needle_] := 1+StringCount[StringTake[source,First[First[StringPosition[source,needle]]]-1],"\n"];
emitPair[name_,a_,b_,tag_] := Module[{raw},
  raw="<|\"Baseline\" -> "<>b["Text"][tag]<> ", \"Corrupted\" -> "<>a["Text"][tag]<>
    ", \"StructuralComparison\" -> "<>ToString[comparison[a["Held"][tag],b["Held"][tag]],InputForm,PageWidth->Infinity]<>"|>";
  WriteString[First[$Output],"WL_LOCAL_S9B_ABLATION_"<>name<>"_"<>tag<>": "<>raw<>"\n"];Flush[First[$Output]]];

(* BEGIN_WITNESS_METHOD *)
(* Deliberately bounded finite search, without computation deadlines. A
   failed candidate, an unresolved equation or an agreeing evaluation is
   never a claim of change. No symbolic inequality proof substitutes for
   a concrete witness. These inputs never enter an engine execution. *)
profileHeads={Global`deltaProfile,Global`velocityProfile,Global`xiProfile,
  Global`rhoBr,Global`fProfile,Global`muPerp,Global`cGamma,Global`jn};
profiles=<|
  Global`deltaProfile->Function[{q},1/(100 q)],
  Global`velocityProfile->Function[{q},1/(10 q)],
  Global`xiProfile->Function[{q},1/(10 q)],
  Global`rhoBr->Function[{q},1+1/q],
  Global`fProfile->Function[{q},1/(100 q)],
  Global`muPerp->Function[{q},4(1+1/q)(1+1/(100 q))^2],
  Global`cGamma->Function[{q},2(1+1/(100 q))],
  Global`jn->Function[{q},-1/(10 q^2)]|>;
(* j_n above is the CAS divergence consequence of the chosen rho and V;
   it is used only where that balance is relevant. Print its verification. *)
point={Global`r->3,Global`rr->4,Global`c0->2,Global`ell->1,
  Global`b->2,Global`bFar->1,Global`ZE->4,Global`ZR->5,Global`impact->2,
  Global`GM->1,Global`n->2,System`K->1,Global`m->1,Global`s->2,Global`rho0->1,
  Global`Phi->1,Global`rAnchor->2,Global`kr->1,Global`kphi->1,
  Global`tau->1,Global`dr->1/2,Global`dphi->1/3,
  Global`xD->1,Global`xV->1,Global`xW->1,
  Global`xiConstant->0,Global`rhoBrConstant->4/3};
(* There are no binders in any candidate accepted by this method. Full ray
   candidates retain all binders and are recorded unresolved without point
   substitution. In particular no inactive quantifier is treated as true. *)
noBinderQ[z_] := FreeQ[z,Inactive|_Integrate|_Limit|_ForAll|_Exists|_ImplicitRegion|_Function|_HoldForm|_HoldComplete];
scalarSymbols[z_] := DeleteDuplicates[Cases[z,q_Symbol /;
  (Context[q]==="Global`" || q===System`K) && !MemberQ[profileHeads,q],{0,Infinity}]];
evaluate[z_] := FullSimplify[z/.Normal[profiles]/.point,Assumptions->True];
finiteExactQ[z_] := (NumericQ[z] && Precision[z]===Infinity &&
  FreeQ[z,Indeterminate|Undefined|DirectedInfinity] && TrueQ[Im[z]==0]);
part[data_,path_List] := Fold[If[AssociationQ[#1],#1[#2],Part[#1,#2]]&,object[data,First[path]],Rest[path]];
heldEquation[left_,right_] := Inactive[Equal][left,right];
profileCertificate[] := Module[{q=Unique["S9bHarness`radius"],derivatives,smooth,tail,balance},
  derivatives=Association[KeyValueMap[Function[{head,fun},head->Table[D[fun[q],{q,j}],{j,0,4}]],profiles]];
  smooth=Resolve[ForAll[q,q>0,And@@(Element[#,Reals]&/@Flatten[Values[derivatives]])],Reals];
  tail={Limit[q profiles[Global`deltaProfile][q],q->Infinity],
    Limit[Sqrt[q] profiles[Global`velocityProfile][q],q->Infinity],
    Limit[q D[profiles[Global`xiProfile][q],q]^2,q->Infinity]};
  balance=FullSimplify[D[q^2 profiles[Global`rhoBr][q] profiles[Global`velocityProfile][q],q]/q^2+
    profiles[Global`jn][q],q>0];
  <|"Profiles"->Normal[profiles],"DerivativeVariable"->q,"DerivativeRange"->q>0,
    "DerivativesThroughOrder4"->derivatives,
    "SmoothRealDecisionComputation"->Inactive[Resolve][Inactive[ForAll][q,q>0,
      And@@(Element[#,Reals]&/@Flatten[Values[derivatives]])],Reals],
    "SmoothRealDecision"->smooth,
    "OpticalCountingComputations"->{Inactive[Limit][q profiles[Global`deltaProfile][q],q->Infinity],
      Inactive[Limit][Sqrt[q] profiles[Global`velocityProfile][q],q->Infinity],
      Inactive[Limit][q D[profiles[Global`xiProfile][q],q]^2,q->Infinity]},
    "OpticalCountingEvaluations"->tail,"MassBalanceResidual"->balance|>];
certificate=.;
(* The only auxiliary local quantities are fixed by the same profile inputs
   in each copy. The copied xiSlope construction is allowed to differ. *)
localFix[data_] := {Global`aMetric->data["State"]["RadialMetric"],
  Global`cSquared->Global`muPerp[Global`r]/Global`rhoBr[Global`r],
  Global`vLocal->Global`velocityProfile[Global`r],
  Global`fLocal->Global`fProfile[Global`r]};
localEvaluate[z_,data_] := evaluate[z/.localFix[data]];
fixedRules[data_,variables_] := Map[#->localEvaluate[#,data]&,variables];
ruleCondition[rules_] := And@@(Equal@@#&/@(List@@@rules));
frequencyReading[a_,b_] := Module[{ar,br,as,bs,var=Global`omega,variables,
  fixes,sets,used,condition,intersection,common,decision,equations},
  variables={Global`aMetric,Global`cSquared,Global`vLocal,Global`omega};
  {ar,br}=Map[Function[data,localEvaluate[data["State"]["Roots"],data]],{a,b}];
  {as,bs}=Map[Function[data,localEvaluate[data["State"]["SelectedRoot"],data]],{a,b}];
  fixes=Map[fixedRules[#,Most[variables]]&,{a,b}];
  sets=MapThread[Function[{rules,roots},Append[rules,var->#]&/@roots],{fixes,{ar,br}}];
  used=MapThread[Append[#1,var->#2]&,{fixes,{as,bs}}];
  condition=(Or@@(ruleCondition/@sets[[1]]))&&(Or@@(ruleCondition/@sets[[2]]));
  intersection=Reduce[condition,variables,Reals];
  decision=FreeQ[intersection,Reduce] && (SameQ[intersection,False] || TrueQ[FullSimplify[intersection/.used[[1]]]]);
  common=TrueQ[FullSimplify[And@@MapThread[Equal,{Last/@used[[1]],Last/@used[[2]]}]]] &&
    TrueQ[FullSimplify[intersection/.used[[1]]]];
  equations=Map[Function[data,<|"Dispersion"->heldEquation[data["State"]["Dispersion"],0],
    "LocalParameterEquations"->(heldEquation[First[#],Last[#]]&/@Take[localFix[data],3]),
    "DispersionAtInputs"->heldEquation[localEvaluate[data["State"]["Dispersion"],data],0]|>],{a,b}];
  <|"FixingEquations"->equations,"SolutionVariables"->variables,
    "SolutionSets"->sets,"SolutionsUsed"->used,
    "Selection"->"engine positive comoving root: positive when vLocal is set to zero",
    "SelectionEvaluations"->Map[Function[data,localEvaluate[
      (data["State"]["SelectedRoot"]/.Global`vLocal->0)>0,data]],{a,b}],
    "SelectedSolutionEquationEvaluations"->MapThread[Function[{data,value},
      localEvaluate[data["State"]["Dispersion"]/.var->value,data]==0],{{a,b},{as,bs}}],
    "IntersectionComputation"->Inactive[Reduce][condition,variables,Reals],
    "SolutionIntersection"->intersection,"IntersectionDecided"->decision,
    "CommonSolutionUsed"->common,
    "AdmissiblePairing"->(decision&&(SameQ[intersection,False]||common))|>];
fixedReading[variables_,a_,b_] := Module[{sets,condition,intersection,used,common},
  used=Map[fixedRules[#,variables]&,{a,b}];sets=List/@used;
  condition=ruleCondition[used[[1]]]&&ruleCondition[used[[2]]];
  intersection=If[variables==={},{{}},Reduce[condition,variables,Reals]];
  common=And@@MapThread[TrueQ[FullSimplify[#1==#2]]&,{Last/@used[[1]],Last/@used[[2]]}];
  <|"FixingEquations"->Map[Function[data,
      heldEquation[#, #/.localFix[data]]&/@variables],{a,b}],
    "SolutionVariables"->variables,"SolutionSets"->sets,"SolutionsUsed"->used,
    "IntersectionComputation"->If[variables==={},Inactive[Intersection][{{}},{{}}],Inactive[Reduce][condition,variables,Reals]],
    "SolutionIntersection"->intersection,"IntersectionDecided"->FreeQ[intersection,Reduce],
    "CommonSolutionUsed"->common,
    "AdmissiblePairing"->(FreeQ[intersection,Reduce]&&(intersection===False||common)),
    "Reading"->If[variables==={},"no auxiliary quantity in this site part; common empty assignment",
      "all fixed local quantities included in the joint solution sets and intersection"]|>];
caseReadings[z_,data_] := Map[Function[piece,Module[{conditions,effective,ev,chosen},
  conditions=Last/@piece[[1]];
  effective=Append[Table[And@@Join[Not/@Take[conditions,j-1],{conditions[[j]]}],{j,Length[conditions]}],And@@(Not/@conditions)];
  ev=localEvaluate[#,data]&/@effective;
  chosen=Flatten[Position[ev,True]];
  <|"Object"->piece,"ConditionsInOrder"->conditions,"EffectiveConditionsIncludingDefault"->effective,
    "EvaluationComputation"->"localEvaluate of each effective condition, including all earlier-case negations",
    "ConditionEvaluations"->ev,"SelectedCase"->chosen,
    "SelectedCondition"->If[Length[chosen]===1,effective[[First[chosen]]],Missing["UNRESOLVED_CASE"]],
    "SelectedConditionEvaluation"->If[Length[chosen]===1,ev[[First[chosen]]],False]|>]],
  DeleteDuplicates[Cases[z,p_Piecewise:>p,{0,Infinity}]]];
specs=<|
 "K1"->{"local",{"WL_LOCAL_S9B_FERMAT_OBJECT","TimePolynomial"},"timePolynomial = Expand["},
 "K2"->{"local",{"WL_LOCAL_S9B_FERMAT_OBJECT","TimePolynomial"},"timePolynomial = Expand["},
 "K3"->{"local",{"WL_LOCAL_S9B_FERMAT_OBJECT","TimePolynomial"},"timePolynomial = Expand["},
 "K4a"->{"optical",{"WL_LOCAL_S9B_OPTICAL_COORDINATE","PushedDensity"},"opticalDensity += (-1)^j"},
 "K4b"->{"optical",{"WL_LOCAL_S9B_OPTICAL_COORDINATE","PushedDensity"},"opticalDensity += (-1)^j"},
 "K4c"->{"slope",{"WL_LOCAL_S9B_PROFILE_ANSATZ","XiDerivative"},"\"XiDerivative\" -> xiSlope"},
 "K5"->{"odd",{"WL_LOCAL_S9B_NONRECIPROCAL_CONSTRUCTION","OneForm",1},"oddOneFormLocal = D["},
 "K6"->{"ray",{"WL_S9B_A_RADAR_LOG_G100","OnDomain"},"finiteSlope = Expand["},
 "K7"->{"ray",{"WL_S9B_B_CONDITION_DEFLECTION"},"If[quantifierMode === \"FIXED_B\""},
 "K8"->{"mass",{"WL_LOCAL_S9B_MASS_BALANCE","Divergence"},"massVector = (massDensity"},
 "K9a"->{"response",{"WL_LOCAL_S9B_C_RESPONSES","Responses",2,"DensityParts",1,"Object"},"responseValues = {c0"},
 "K9b"->{"response",{"WL_LOCAL_S9B_C_RESPONSES","Responses",3,"DensityParts",1,"Object"},"responseValues = {c0"},
 "K10"->{"ray",{"WL_S9B_A_NONRECIPROCAL_PATH_DEPENDENCE","OnDomain"},"physicalPropagation = branchGateLocal"},
 "K11"->{"local",{"WL_LOCAL_S9B_FERMAT_OBJECT","TimePolynomial"},"timePolynomial = Expand["}|>;
searchWitness[name_,a_,b_] := Module[{spec=specs[name],kind,path,aa,bb,conditions,reading,
  evaluatedConditions,values,roles,symbols,fixed,assigned,unknown,entry,ceConditions={True,True},
  density,checks,location,detail,cases},
  {kind,path}=Take[spec,2];note["witness begin",name];
  aa=part[a,path];bb=part[b,path];
  location=<|"Tag"->First[path],"PayloadLocation"->Rest[path],"SitePart"->"value",
    "Source"->engine,"SourceLine"->sourceLine[Last[spec]],"EntryConstruction"->Last[spec]|>;
  If[kind==="ray",
    (* Read the full candidate and fixing equations. The finite search stops
       at the unresolved full ray boundary problem; no guessed turning radius,
       impact, endpoint cutoff or inactive-domain decision is used. *)
    entry=<|"Result"->"none found","Site"->location,"Attempt"-><|
      "CandidateParts"->{aa,bb},"ProfilesTried"->Normal[profiles],
      "ConstructionEquationsAfterProfiles"->{a["State"]["RayConstruction"],b["State"]["RayConstruction"]}/.Normal[profiles],
      "UnassignedDependentQuantities"->{Global`rTurnFlyby,Global`rTurnRadar,Global`impactRadar},
      "Reason"->"no exact solution of the full endpoint/turning construction and all ranged conditions was obtained; no scalar point or bound-variable substitution was made"|>|>;
    note["witness end",name];Return[entry]];
  If[!ValueQ[certificate],certificate=profileCertificate[]];
  If[kind==="response",
    ceConditions=Map[If[Head[#]===ConditionalExpression,#[[2]],True]&,{aa,bb}];
    {aa,bb}=Map[If[Head[#]===ConditionalExpression,#[[1]],#]&,{aa,bb}]];
  conditions=MapThread[Function[{data,ce},Join[
    {data["State"]["Assumptions"]},
    If[MemberQ[{"local","optical","odd"},kind],{data["State"]["LocalAssumptions"]},{}],
    If[MemberQ[{"optical","odd"},kind],
      {data["State"]["PhysicalPropagation"],data["State"]["PhysicalTraversal"],data["State"]["SpeedChart"]},{}],
    If[kind==="response",{ce,object[data,"WL_LOCAL_S9B_C_RESPONSES"]["Domains"][[If[name==="K9a",2,3]]],
      object[data,"WL_LOCAL_S9B_C_RESPONSES"]["Responses"][[If[name==="K9a",2,3]]]["SuppliedBulkDensityPremise"]},{}]]],
    {{a,b},ceConditions}];
  If[!noBinderQ[{aa,bb,conditions}],
    note["witness end",name];Return[<|"Result"->"none found","Site"->location,
      "Attempt"->"candidate contains a binder or unevaluated functional outside the exact local evaluator"|>]];
  cases=MapThread[caseReadings,{ {{aa,conditions[[1]]},{bb,conditions[[2]]}}, {a,b}}];
  conditions=MapThread[Join[#1,Lookup[#2,"SelectedCondition",{}]]&,{conditions,cases}];
  evaluatedConditions=MapThread[Function[{items,data},localEvaluate[#,data]&/@items],{conditions,{a,b}}];
  values={localEvaluate[aa,a],localEvaluate[bb,b]};
  symbols=scalarSymbols[{aa,bb,conditions}];
  fixed=Intersection[symbols,{Global`aMetric,Global`cSquared,Global`vLocal,Global`fLocal}];
  assigned=Complement[symbols,fixed];unknown=Complement[assigned,First/@point];
  If[MemberQ[{"local","optical","odd"},kind],fixed=Union[fixed,{Global`aMetric,Global`cSquared,Global`vLocal}]];
  reading=If[MemberQ[{"local","optical","odd"},kind],frequencyReading[a,b],fixedReading[fixed,a,b]];
  roles=Join[Association[Map[#-><|"Role"->"assigned at the point","Value"->(#/.point)|>&,First/@point]],
    Association[Map[Function[var,var-><|"Role"->"fixed by each copy's construction",
      "Equations"->Map[Function[data,heldEquation[var,var/.localFix[data]]],{a,b}],
      "Evaluations"->{localEvaluate[var,a],localEvaluate[var,b]}|>],fixed]],
    If[MemberQ[{"local","optical","odd"},kind],<|Global`omega-><|"Role"->"fixed by each copy's dispersion solution and positive comoving root selection","Solutions"->reading["SolutionsUsed"]|>|>,<||>],
    If[kind==="response",<|Global`rho-><|"Role"->"bound differentiation variable of pressure; then substituted by local or asymptotic density","Range"->Reals,"DerivativeOrder"->1|>|>,<||>],
    Association[Map[#-><|"Role"->"assigned radial profile","Value"->profiles[#]|>&,
      DeleteDuplicates[Cases[{aa,bb,conditions},q_ /; MemberQ[profileHeads,q],{0,Infinity},Heads->True]]]]];
  density=If[kind==="response",<|"DensityComputation"->Global`rho0(1+Global`fProfile[Global`r]),
    "DensityEvaluation"->evaluate[Global`rho0(1+Global`fProfile[Global`r])],
    "Part"->If[TrueQ[evaluate[Global`rho0(1+Global`fProfile[Global`r])]>0],"rho > 0","rho <= 0"],
    "PartitionLabelsAreConditions"->False|>,None];
  detail=Switch[kind,
    "slope",<|"BaselineDifferentiation"->Inactive[D][Global`xiProfile[Global`r],Global`r],
      "CorruptedDifferentiation"->Inactive[D][Global`xiConstant,Global`r],"Evaluations"->{values[[2]],values[[1]]}|>,
    "mass",<|"CoordinateDivergenceComputation"->Inactive[D][Global`r^2 Global`rhoBr[Global`r] Global`velocityProfile[Global`r],Global`r]/Global`r^2,
      "ConstantDensity"->(Global`rhoBrConstant/.point)|>,
    "response",<|"BulkPressure"->System`K Global`densityPower[Global`rho,Global`n],
      "SoundSquareComputation"->Inactive[D][System`K Global`densityPower[Global`rho,Global`n],Global`rho]/Global`m,
      "LocalFractionEquation"->heldEquation[Global`fLocal,Global`fProfile[Global`r]],
      "DensityPowerReading"-><|"Exponent"->evaluate[Global`s],
        "IsZeroExponent"->TrueQ[evaluate[Global`s]==0],
        "OrdinaryPowerAtThisInput"->If[name==="K9b",{values[[1]],values[[2]]},None]|>,
      "KnifeDefinition"->knives[name]|>,
    "optical",<|"DerivativeConstruction"->knives[name],
      "Reading"->"evaluate the emitted retained-ring coefficient sum at all counting markers equal to one; frozen derivatives are those of the copied rayDerivative, not independent fitted jets"|>,
    _,<|"Reading"->"evaluate the emitted local construction; local metric/speed/velocity inputs are fixed from the printed physical profiles"|>];
  checks=<|"AllSymbolsAssignedOrFixed"->(unknown==={}),
    "ConditionsDecidedTrue"->And@@Flatten[evaluatedConditions],
    "ExactFiniteValues"->And@@(finiteExactQ/@values),
    "ExactDifferenceDecision"->If[And@@(finiteExactQ/@values),FullSimplify[values[[1]]!=values[[2]]],False],
    "SolutionPairing"->reading["AdmissiblePairing"],
    "SelectedSolutionsAdmissible"->If[MemberQ[{"local","optical","odd"},kind],
      And@@Join[reading["SelectionEvaluations"],reading["SelectedSolutionEquationEvaluations"]],True],
    "CaseSelectionsDecided"->And@@(TrueQ/@Flatten[Lookup[#, "SelectedConditionEvaluation",{}]&/@cases]),
    "SmoothProfiles"->certificate["SmoothRealDecision"]|>;
  entry=<|"Result"->If[And@@(TrueQ/@Values[checks]),"witness","none found"],
    "OperandOrder"->{"corrupted","baseline"},"Site"->location,"IndependentPoint"->point,"SymbolRoles"->roles,
    "RangedSymbolsInPart"->{},"BoundSymbolsInPart"->{},"UnknownSymbols"->unknown,
    "AdmissibleProfiles"->certificate,"ReadingChoices"->detail,
    "FixingAndSolutions"->reading,"OrderedCaseReadings"->cases,"AttachedConditions"->conditions,
    "ConditionEvaluationComputation"->"substitute fixed radial functions and copy-local equations, then the printed point; FullSimplify with Assumptions -> True",
    "ConditionEvaluations"->evaluatedConditions,"DensityPart"->density,
    "OperandsAtInputs"->{aa/.localFix[a]/.Normal[profiles]/.point,bb/.localFix[b]/.Normal[profiles]/.point},
    "Evaluations"->values,"Arithmetic"->"exact real algebraic","ErrorBound"->0,"Checks"->checks,
    "Attempt"->"one printed exact profile configuration and scalar point; no unresolved condition certifies a witness"|>;
  note["witness end",name];entry];
(* END_WITNESS_METHOD *)
S9bHarness`compactBoundary = "(* COMPACT_K11_BOUNDARY *)";
S9bHarness`prepareCompact[path_, name_] := Module[{text, target, position, prefix},
  text = Import[path, "Text"];
  If[!StringQ[S9bHarness`compactBoundary] ||
    StringCount[text, S9bHarness`compactBoundary] != 1, Quit[93]];
  position = First[First[StringPosition[text, S9bHarness`compactBoundary]]];
  prefix = StringTake[text, position - 1];
  If[!StringQ[prefix] || StringLength[prefix] == 0 ||
    StringLength[prefix] >= StringLength[text] ||
    StringContainsQ[prefix, S9bHarness`compactBoundary], Quit[93]];
  target = FileNameJoin[{work, name <> ".wl"}];
  Export[target, prefix, "Text"];
  <|"Path" -> target, "Text" -> prefix|>];
baseline=run[engine,"BASELINE"];
print["MANIFEST",knives];
Do[copy=mutate[knife];corrupted=run[copy,knife];
  If[Keys[corrupted["Text"]]=!=Keys[baseline["Text"]],Quit[96]];
  Do[note["structural pair",knife<>"_"<>tag];emitPair[knife,corrupted,baseline,tag],{tag,Keys[baseline["Text"]]}];
  print[knife<>"_WITNESS_SEARCH",searchWitness[knife,corrupted,baseline]],
  {knife,DeleteCases[Keys[knives],"K11"]}];
copy=mutate["K11"];
S9bHarness`compactInputs={
  S9bHarness`prepareCompact[engine,"K11_COMPACT_BASELINE"],
  S9bHarness`prepareCompact[copy,"K11_COMPACT_CORRUPTED"]};
If[S9bHarness`compactInputs[[2]]["Text"]=!=
  StringReplace[S9bHarness`compactInputs[[1]]["Text"],Rule@@knives["K11"]],Quit[93]];
compactBaseline=run[S9bHarness`compactInputs[[1]]["Path"],"K11_COMPACT_BASELINE"];
compactCorrupted=run[S9bHarness`compactInputs[[2]]["Path"],"K11_COMPACT_CORRUPTED"];
coverage=<|"Method"->"exact symbolic mechanical prefix extraction",
  "Boundary"->S9bHarness`compactBoundary,"RetainedGrades"->Global`grades,
  "LocalConstructions"->"untruncated dispersion/Fermat and retained-ring nonreciprocal one-form",
  "AdditionalTruncation"->None,"SampleDomain"->"symbolic constructor domains; coverage only, not an attached condition",
  "Arithmetic"->"exact","Samples"->{},"Seeds"->{}|>;
Do[note["structural pair","K11_"<>tag];
  If[KeyExistsQ[compactBaseline["Text"],tag]&&KeyExistsQ[compactCorrupted["Text"],tag],
    raw="<|\"Baseline\" -> "<>baseline["Text"][tag]<>
      ", \"DifferencePair\" -> {"<>compactCorrupted["Text"][tag]<>", "<>compactBaseline["Text"][tag]<>"}"<>
      ", \"DifferencePairComparison\" -> "<>ToString[comparison[compactCorrupted["Held"][tag],compactBaseline["Held"][tag]],InputForm,PageWidth->Infinity]<>
      ", \"MethodResidualPair\" -> {"<>compactBaseline["Text"][tag]<>", "<>baseline["Text"][tag]<>"}"<>
      ", \"MethodResidualPairComparison\" -> "<>ToString[comparison[compactBaseline["Held"][tag],baseline["Held"][tag]],InputForm,PageWidth->Infinity]<>
      ", \"Coverage\" -> "<>ToString[coverage,InputForm,PageWidth->Infinity]<>"|>",
    raw="<|\"Baseline\" -> "<>baseline["Text"][tag]<>
      ", \"DifferencePair\" -> Missing[\"NOT_EVALUATED\"], \"DifferencePairMark\" -> \"not compared\", \"MethodResidualPair\" -> Missing[\"NOT_EVALUATED\"], \"MethodResidualPairMark\" -> \"not compared\", \"Coverage\" -> "<>
      ToString[coverage,InputForm,PageWidth->Infinity]<>"|>"];
  WriteString[First[$Output],"WL_LOCAL_S9B_ABLATION_K11_"<>tag<>": "<>raw<>"\n"];Flush[First[$Output]],
  {tag,Keys[baseline["Text"]]}];
print["K11_WITNESS_SEARCH",searchWitness["K11",compactCorrupted,compactBaseline]];
Do[note["structural pair","HARNESS_DEAD_PATH_"<>tag];emitPair["HARNESS_DEAD_PATH",baseline,baseline,tag],{tag,Keys[baseline["Text"]]}];
note["complete","ALL"];Close[S9bHarness`progress];
End[];
