(* Amendment 12. One kernel; unchanged engine baseline; one source
   replacement per knife. Literal frozen emitted syntax is the only tag
   comparison. Witness searches are separate, one per knife. *)
Begin["S9bHarness`"];
$HistoryLength = 0;
root = Directory[];
engine = FileNameJoin[{root,"research","pde_ledger_v3","mathematica","S9b_light_bending_mathematica_audit.wl"}];
work = FileNameJoin[{root,"_scratch","s9b_wl_repair10","ablation-final"}];
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
 "K10" -> {"propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,\n  physicalPropagation && speedChart && rhoBr[r] != 0 &&\n  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];", "propagationGate = Not[rayConstruction && Inactive[ForAll][r, traversedRadii,\n  physicalPropagation && speedChart && rhoBr[r] != 0 &&\n  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]]];"},
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
(* Local candidates have no binders. Functional candidates use scoped
   substitution below; an unevaluated quantified condition decides nothing. *)
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
    "ConditionEvaluations"->ev,"SelectedEffectiveConditions"->If[Length[chosen]===1,{effective[[First[chosen]]]},{}],
    "SelectedCondition"->If[Length[chosen]===1,effective[[First[chosen]]],Missing["UNRESOLVED_CASE"]],
    "SelectedConditionEvaluation"->If[Length[chosen]===1,ev[[First[chosen]]],False]|>]],
  DeleteDuplicates[Cases[z,p_Piecewise:>p,{0,Infinity}]]];
(* BEGIN_ROUND9_SITE_SEARCH *)
(* These three searches evaluate only parts reached by K6, K7 or K10.
   They do not change the comparison, another knife's search, or the engine.
   The separate boolean path is a domain comparison, not a numerical value
   comparison. An unsuccessful attempt reports evaluated objects only. *)
decidedBooleanQ[z_] := SameQ[z,True] || SameQ[z,False];
domainDifferenceDecision[a_,b_] := If[decidedBooleanQ[a] && decidedBooleanQ[b],Xor[a,b],False];

(* Rename bound variables before substituting a scalar point. Actual engine
   quantifiers, definite integrals and limits each bind one radius. Fresh
   names and their original scopes are printed with the evaluated attempt. *)
siteScope[z_] := Module[{bindings={},walk,result},
  walk[q_] := Module[{h=Head[q],v,fresh,args},Which[
    AssociationQ[q],AssociationMap[walk[q[#]]&,Keys[q]],
    AtomQ[q],q,
    MemberQ[{ForAll,Exists,Inactive[ForAll],Inactive[Exists]},h],
      v=Flatten[{q[[1]]}];fresh=Table[Unique["S9bHarness`siteBound"],{Length[v]}];
      AppendTo[bindings,<|"Role"->"bound","Head"->h,"OriginalVariables"->v,
        "Variables"->fresh,"Scope"->Rest[List@@q]|>];
      args=walk /@ (Rest[List@@q]/.Thread[v->fresh]);
      h@@Prepend[args,If[ListQ[q[[1]]],fresh,First[fresh]]],
    MemberQ[{Integrate,Inactive[Integrate]},h] && Length[q]===2 && ListQ[q[[2]]] && Length[q[[2]]]===3,
      v=q[[2,1]];fresh=Unique["S9bHarness`siteBound"];
      AppendTo[bindings,<|"Role"->"bound","Head"->h,"OriginalVariables"->{v},
        "Variables"->{fresh},"Range"->{q[[2,2]],q[[2,3]]}|>];
      h[walk[q[[1]]/.v->fresh],{fresh,walk[q[[2,2]]],walk[q[[2,3]]]}],
    MemberQ[{Limit,Inactive[Limit]},h] && Head[q[[2]]]===Rule,
      v=q[[2,1]];fresh=Unique["S9bHarness`siteBound"];
      AppendTo[bindings,<|"Role"->"bound","Head"->h,"OriginalVariables"->{v},
        "Variables"->{fresh},"Limit"->q[[2,2]]|>];
      h[walk[q[[1]]/.v->fresh],fresh->walk[q[[2,2]]]],
    True,Map[walk,q]]];
  result=walk[z];<|"Object"->result,"Bindings"->bindings|>];
siteExact[z_] := Block[{$Assumptions=True},Module[{walk},
  walk[q_Association] := AssociationMap[walk[q[#]]&,Keys[q]];
  walk[q_List] := walk/@q;
  walk[q_] := Module[{a=q},
    a=a/.HoldPattern[Inactive[Integrate][f_,range:{_,_,_}]] :>
      Integrate[f,range,GenerateConditions->True,Assumptions->True];
    a=Activate[a,ForAll|Exists|Limit];
    If[!FreeQ[a,ForAll|Exists],a=Resolve[a,Reals]];
    FullSimplify[a,Assumptions->True]];
  walk[z]]];
siteEvaluate[z_,data_] := Module[{scoped=siteScope[z],inputs},
  inputs=scoped["Object"]/.localFix[data]/.Normal[profiles]/.point;
  <|"ScopedObject"->scoped["Object"],"BoundSymbols"->scoped["Bindings"],
    "AtInputs"->inputs,"Evaluation"->siteExact[inputs],
    "Computation"->"alpha-rename binders; substitute copy-local fixing equations, specific profile functions and the point; ordinary Integrate with GenerateConditions -> True; active Resolve over Reals; FullSimplify with Assumptions -> True"|>];
siteRoles[expressions_,fixed_,reading_] := Module[{symbols,assigned,unknown,roles},
  symbols=scalarSymbols[expressions];assigned=Complement[symbols,fixed];
  unknown=Complement[assigned,First/@point];
  roles=Join[Association[Map[#-><|"Role"->"assigned at the point","Value"->(#/.point)|>&,First/@point]],
    Association[Map[#-><|"Role"->"fixed by each copy's construction","SolutionsUsed"->reading["SolutionsUsed"]|>&,fixed]],
    Association[KeyValueMap[#1-><|"Role"->"assigned radial profile","Function"->#2|>&,profiles]]];
  <|"Roles"->roles,"UnknownSymbols"->unknown|>];
siteLocation[tag_,path_,kind_,entry_] := <|"Tag"->tag,"PayloadLocation"->path,
  "SitePart"->kind,"Source"->engine,"SourceLine"->sourceLine[entry],"EntryConstruction"->entry|>;

searchReachedSite["K6",mutatedData_,baselineData_] := Module[{tag="WL_LOCAL_S9B_RADAR_SLOPE_DERIVATION",path={"IntegralKernels","100"},
  parts,conditions,conditionRecords,valueRecords,sourceRecords,reading,roles,checks,values},
  parts=Map[part[#,Prepend[path,tag]]&,{mutatedData,baselineData}];
  conditions=Map[Function[data,{data["State"]["Assumptions"],data["State"]["LocalAssumptions"],
    data["State"]["PhysicalPropagation"],data["State"]["PhysicalTraversal"],data["State"]["SpeedChart"],
    Global`r>Global`b}],{mutatedData,baselineData}];
  conditionRecords=MapThread[siteEvaluate,{conditions,{mutatedData,baselineData}}];
  valueRecords=MapThread[siteEvaluate,{parts,{mutatedData,baselineData}}];values=Lookup[valueRecords,"Evaluation"];
  sourceRecords=Map[Function[data,siteEvaluate[object[data,tag]["Source"]["100"],data]],{mutatedData,baselineData}];
  reading=frequencyReading[mutatedData,baselineData];
  roles=siteRoles[{parts,conditions},{Global`aMetric,Global`cSquared,Global`vLocal,Global`omega},reading];
  checks=<|"ConditionsDecidedTrue"->And@@(TrueQ/@Flatten[Lookup[conditionRecords,"Evaluation"]]),
    "AllSymbolsAssignedOrFixed"->(roles["UnknownSymbols"]==={}),
    "OrdinaryFiniteSourceValues"->And@@(finiteExactQ/@Lookup[sourceRecords,"Evaluation"]),
    "ExactFiniteValues"->And@@(finiteExactQ/@values),
    "ExactDifferenceDecision"->If[And@@(finiteExactQ/@values),FullSimplify[values[[1]]!=values[[2]]],False],
    "SolutionPairing"->reading["AdmissiblePairing"],
    "SelectedSolutionsAdmissible"->And@@Join[reading["SelectionEvaluations"],reading["SelectedSolutionEquationEvaluations"]],
    "SmoothProfiles"->profileCertificate[]["SmoothRealDecision"]|>;
  <|"Result"->If[And@@(TrueQ/@Values[checks]),"witness","undecided"],"OperandOrder"->{"corrupted","baseline"},
    "Site"->siteLocation[tag,path,"value","kernel = Simplify[(finiteSlope - boundary)"],
    "IndependentPoint"->point,"SymbolRoles"->roles["Roles"],"UnknownSymbols"->roles["UnknownSymbols"],
    "AdmissibleProfiles"->profileCertificate[],"FixingAndSolutions"->reading,
    "ReadingChoices"-><|"RadialArgument"->"r is the assigned argument of the emitted local kernel, not the bound radius in its source integrals",
      "ConstructionSourceOrder"->{"oneWayER","roundTrip"},"SourceEvaluationRecords"->sourceRecords,
      "AttachedDomainReading"->"the emitted local derivation container has no attached ray gate; constructor assumptions and the kernel's open integration chart r > b are evaluated explicitly"|>,
    "AttachedConditions"->conditions,"ConditionEvaluationRecords"->conditionRecords,
    "ConditionEvaluations"->Lookup[conditionRecords,"Evaluation"],"ValueEvaluationRecords"->valueRecords,
    "Evaluations"->values,"Arithmetic"->"exact","ErrorBound"->0,"Checks"->checks,
    "Attempt"->"evaluated grade-100 local radar kernels, finite-endpoint source integrals and the printed constructor conditions at one exact profile/point configuration"|>];

(* Round 10: decide the engine's full branch condition. The flat profiles are
   witness inputs only. The ray equations are read from each actual copy;
   no turning radius or radar impact is assigned independently. *)
k10Exact[z_,assumptions_:True] := Block[{$Assumptions=True},Module[{result=z},
 result=result/.HoldPattern[Inactive[Integrate][f_,range:{_,_,_}]] :>
   Integrate[f,range,GenerateConditions->True,Assumptions->assumptions];
 result=Activate[result,ForAll|Exists|Limit];
 If[!FreeQ[result,ForAll|Exists],result=Resolve[result,Reals]];
 FullSimplify[result,Assumptions->assumptions]]];
k10Construction[data_] := Module[{scope,input,variables={Global`rTurnFlyby,Global`rTurnRadar,Global`impactRadar},
 turningVariables={Global`rTurnFlyby,Global`rTurnRadar},j=Global`impactRadar,conjuncts,algebraic,
 equations,turningSolutions,turningDomains,admissibleTurning,turningRules,domain,upper,interior,boundaryInput,boundary,
 integralInput,integrated,statement,solutions,used,checks},
 scope=siteScope[data["State"]["RayConstruction"]];
 input=scope["Object"]/.localFix[data]/.Normal[profiles]/.point;
 conjuncts=If[Head[input]===And,List@@input,{input}];
 algebraic=k10Exact[And@@Select[conjuncts,FreeQ[#,Inactive[Integrate]]&]];
 equations=Select[conjuncts,Head[#]===Equal&&FreeQ[#,Inactive[Integrate]]&];
 turningSolutions=Solve[And@@equations,turningVariables,Reals];
 turningDomains=If[ListQ[turningSolutions],Reduce[algebraic/.#,j,Reals]&/@turningSolutions,{}];
 admissibleTurning=If[ListQ[turningSolutions],Pick[turningSolutions,(FreeQ[#,Reduce]&&#=!=False)&/@turningDomains,True],{}];
 If[!ListQ[turningSolutions]||!And@@(FreeQ[#,Reduce]&/@turningDomains)||Length[admissibleTurning]!=1,
  Return[<|"Decided"->False,"ScopedObject"->scope,"AtInputs"->input,
    "TurningEquations"->equations,"TurningSolutions"->turningSolutions,"TurningCandidateDomains"->turningDomains,
    "UnmetRequirements"->{"the constructor's turning equations did not decide a unique turning-rule substitution"}|>]];
 turningRules=First[admissibleTurning];domain=Reduce[algebraic/.turningRules,j,Reals];
 upper=Min@@localEvaluate[object[data,"WL_LOCAL_S9B_RAY_DOMAIN"]["Endpoints"],data];
 interior=domain&&j<upper;
 (* Integrate on the open admissible interior. The closed upper endpoint is
    evaluated separately from the same constructor, so no solution is lost. *)
 integralInput=input/.turningRules;
 integrated=k10Exact[integralInput,interior];
 boundaryInput=integralInput/.j->upper;boundary=k10Exact[boundaryInput];
 statement=And[algebraic,integrated,interior];
 statement=Reduce[statement,variables,Reals];
 solutions=Solve[statement,variables,Reals];
 used=If[ListQ[solutions]&&Length[solutions]===1,First[solutions],{}];
 checks=<|"OrdinaryIntegralsAndActiveQuantifiersDecided"->FreeQ[integrated,Inactive|Integrate|ForAll|Exists|Limit|ConditionalExpression],
   "CompleteInteriorSolutionDecision"->FreeQ[statement,Reduce],
   "ClosedUpperEndpointExcludedByConstructor"->SameQ[boundary,False],
   "UniqueFlatConnectedSolution"->(ListQ[solutions]&&Length[solutions]===1)|>;
 <|"Decided"->And@@(TrueQ/@Values[checks]),"FixingEquations"->data["State"]["RayConstruction"],
   "ScopedObject"->scope["Object"],"BoundSymbols"->scope["Bindings"],"AtInputs"->input,
   "TurningEquationExtraction"->"top-level Equal conjuncts without an integral; the original endpoint equation and all its integrals are retained in the next step",
   "TurningEquations"->equations,"TurningSolutions"->turningSolutions,"TurningCandidateDomains"->turningDomains,"TurningRulesUsed"->turningRules,
   "NonintegralConstructorConditions"->algebraic,"ImpactDomainComputation"->Inactive[Reduce][algebraic/.turningRules,j,Reals],
   "ImpactDomain"->domain,"EndpointUpperBound"->upper,"InteriorRange"->interior,
   "OrdinaryIntegralInput"->integralInput,"OrdinaryIntegralAndQuantifierEvaluation"->integrated,
   "IntegrationAssumptionsReading"->"only the impact range derived from the engine constructor; its closed upper endpoint is checked separately",
   "ClosedUpperEndpointInput"->boundaryInput,"ClosedUpperEndpointEvaluation"->boundary,
   "SolutionStatement"->statement,"SolutionVariables"->variables,"SolutionSets"->solutions,
   "SolutionsUsed"->used,"Selection"->"the engine's flat-connected branch; these witness profiles are already flat, and the complete construction has one solution",
   "CASOperations"->{"Solve turning equalities","Reduce original nonintegral conjuncts","ordinary Integrate with GenerateConditions -> True on the derived open impact range","active Resolve","evaluate the closed endpoint","Reduce full intersection of constructor constraints","Solve the resulting complete real solution set"},
   "Checks"->checks|>];

searchReachedSite["K10",mutatedData_,baselineData_] := Block[{profiles=<|
 Global`deltaProfile->Function[{q},0],Global`velocityProfile->Function[{q},0],Global`xiProfile->Function[{q},0],
 Global`rhoBr->Function[{q},1],Global`fProfile->Function[{q},0],Global`muPerp->Function[{q},4],
 Global`cGamma->Function[{q},2],Global`jn->Function[{q},0]|>,
 point=Join[DeleteCases[point,Rule[Global`GM|Global`ZE|Global`ZR,_]],{Global`GM->0,Global`ZE->2,Global`ZR->2}]},
 Module[{tag="WL_S9B_BRANCH_EXISTENCE",datasets={mutatedData,baselineData},construction,variables,
  sets,used,intersection,reading,records,conditionRecords,traversalRecords,values,roles,scopeSymbols,boundRoles,checks},
 construction=k10Construction/@datasets;variables={Global`rTurnFlyby,Global`rTurnRadar,Global`impactRadar};
 If[!And@@(TrueQ/@Lookup[construction,"Decided"]),Return[<|"Result"->"undecided",
   "Site"->siteLocation[tag,{},"domain","propagationGate = rayConstruction"],
   "IndependentPoint"->point,"AdmissibleProfiles"->profileCertificate[],"ConstructionReadings"->construction,
   "Attempt"->"actual completed branch-condition site; turning/endpoint construction remained unsettled"|>]];
 sets=Lookup[construction,"SolutionStatement"];used=Lookup[construction,"SolutionsUsed"];
 intersection=Reduce[And@@sets,variables,Reals];
 reading=<|"FixingEquations"->Lookup[construction,"FixingEquations"],"SolutionVariables"->variables,
   "SolutionSets"->Lookup[construction,"SolutionSets"],"SolutionsUsed"->used,
   "IntersectionComputation"->Inactive[Reduce][And@@sets,variables,Reals],"SolutionIntersection"->intersection,
   "IntersectionDecided"->FreeQ[intersection,Reduce],
   "CommonSolutionUsed"->And@@(TrueQ[FullSimplify[intersection/.#]]&/@used),
   "Reading"->"both full ray constructions are solved at the same independent inputs before either post-construction gate is evaluated"|>;
 records=MapThread[Function[{data,solution},siteEvaluate[object[data,tag]/.solution,data]],{datasets,used}];
 values=Lookup[records,"Evaluation"];
 conditionRecords=MapThread[Function[{data,solution},siteEvaluate[data["State"]["Assumptions"]/.solution,data]],{datasets,used}];
 traversalRecords=MapThread[Function[{data,solution},siteEvaluate[object[data,"WL_S9B_PATH_TRAVERSAL"]/.solution,data]],{datasets,used}];
 roles=siteRoles[{Lookup[records,"ScopedObject"],Lookup[conditionRecords,"ScopedObject"]},variables,reading];
 boundRoles=Flatten[Lookup[records,"BoundSymbols"],1];
 checks=<|"BothWholeConstructionsDecided"->And@@(TrueQ/@Lookup[construction,"Decided"]),
   "OtherConditionsDecidedTrue"->And@@(TrueQ/@Lookup[conditionRecords,"Evaluation"]),
   "OtherTraversalGateDecidedTrue"->And@@(TrueQ/@Lookup[traversalRecords,"Evaluation"]),
   "AllSymbolsAssignedOrFixed"->(roles["UnknownSymbols"]==={}),
   "ComparedConditionsDecided"->And@@(decidedBooleanQ/@values),
   "ExactDifferenceDecision"->domainDifferenceDecision@@values,
   "IntersectionDecidedAndCommonSolutionUsed"->(TrueQ[reading["IntersectionDecided"]]&&TrueQ[reading["CommonSolutionUsed"]]),
   "WholeEmittedObjectsEvaluated"->FreeQ[values,Inactive|Integrate|ForAll|Exists|Limit|ConditionalExpression],
   "SmoothProfiles"->profileCertificate[]["SmoothRealDecision"]|>;
 <|"Result"->If[And@@(TrueQ/@Values[checks]),"witness","undecided"],"OperandOrder"->{"corrupted","baseline"},
   "Site"->siteLocation[tag,{},"domain","propagationGate = rayConstruction"],
   "IndependentPoint"->point,"SymbolRoles"->roles["Roles"],"BoundSymbolRoles"->boundRoles,
   "UnknownSymbols"->roles["UnknownSymbols"],"AdmissibleProfiles"->profileCertificate[],
   "FixingAndSolutions"->reading,"ConstructionReadings"->construction,
   "ReadingChoices"-><|"ComparedCondition"->"the complete propagationGate is negated, after rayConstruction and the along-ray ForAll are assembled; the local predicate and the interior of that ForAll are unchanged",
     "PayloadReading"->"BRANCH_EXISTENCE directly emits this full condition at engine 275; there is no containing conditional value or additional case selection",
     "SourceLocations"->{262,263,264,266,271,275},
     "OtherGate"->"PATH_TRAVERSAL is an independently emitted condition; it is computed true in both copies, not assumed",
     "BoundVariableReading"->"integral and quantified radii are alpha-renamed before any point substitution; their whole source ranges are retained"|>,
   "AttachedConditions"->Lookup[conditionRecords,"ScopedObject"],"ConditionEvaluationRecords"->conditionRecords,
   "ConditionEvaluations"->Lookup[conditionRecords,"Evaluation"],"OtherTraversalGateRecords"->traversalRecords,
   "DomainEvaluationRecords"->records,"EmittedObjectEvaluationRecords"->records,"Evaluations"->values,
   "Arithmetic"->"exact ordinary integration and active real quantifier decisions","ErrorBound"->0,"Checks"->checks,
   "Attempt"->"both full post-construction branch conditions at one exact flat profile/point configuration; both complete turning/endpoint solution sets and their intersection are computed"|>]];

(* END_ROUND9_SITE_SEARCH *)

(* BEGIN_AMENDMENT12_CATALOG *)
(* Exposed mathematical objects, not their arithmetic subexpressions, are
   parts. Associations, density wrappers and solution-branch lists are
   containers. Tensors and equation/solution sets are mathematical objects.
   All identities below use engine keys, row names, basis or branch meaning;
   integer storage offsets are never identities across copies. *)
semanticPaths=<|
 "K1"->{"TimePolynomial"},"K2"->{"TimePolynomial"},"K3"->{"TimePolynomial"},
 "K4a"->{"PushedDensity"},"K4b"->{"PushedDensity"},"K4c"->{"XiDerivative"},
 "K5"->{"OneForm","component along dr (basis dr,dphi; engine lines 146,192)"},
 "K6"->{"IntegralKernels","100"},"K8"->{"Divergence"},
 "K9a"->{"Responses","FIXED_RATIO","DensityParts","rho > 0","Object"},
 "K9b"->{"Responses","POWER","DensityParts","rho > 0","Object"},
 "K10"->{},"K11"->{"TimePolynomial"}|>;
rowNames={"CONSTANT","FIXED_RATIO","POWER"};
sourceRanges=<|"local dispersion"->{58,188},"odd form"->{192,204},
 "optical density"->{230,242},"ray domain"->{247,281},"radar source"->{356,377},
 "profile reduction"->{395,425},"mass"->{464,472},"velocity branches"->{490,577},
 "shared B"->{584,592},"shared A"->{593,599},"response rows"->{606,636},
 "forward premise"->{651,670}|>;
sourceAnchor[tag_] := Which[
 StringStartsQ[tag,"WL_LOCAL_S9B_SUPPLIED"],106,
 StringStartsQ[tag,"WL_LOCAL_S9B_PROFILE"],122,
 StringContainsQ[tag,"BRANCH_FREQUENCIES"],177,
 StringContainsQ[tag,"LOCAL_S9B_BRANCH_TYPES"],184,
 StringContainsQ[tag,"HAMILTON"],185,StringContainsQ[tag,"FERMAT"],187,
 StringContainsQ[tag,"LOCAL_S9B_NONRECIPROCAL"],200,
 StringContainsQ[tag,"OPTICAL_COORDINATE"],241,
 StringContainsQ[tag,"LOCAL_S9B_RAY_DOMAIN"],270,
 StringContainsQ[tag,"RADAR_SLOPE_DERIVATION"],374,
 StringContainsQ[tag,"ABEL_REDUCTION"],426,
 StringContainsQ[tag,"LOCAL_S9B_MASS_BALANCE"],470,
 StringContainsQ[tag,"LOCAL_S9B_C_RESPONSES"],624,
 StringContainsQ[tag,"F_MASS_SOLUTION"],662,
 StringContainsQ[tag,"_C_"],630,
 StringContainsQ[tag,"_B_"],588,
 StringContainsQ[tag,"_A_"],596,
 StringContainsQ[tag,"BRANCH_"]||StringContainsQ[tag,"PATH_TRAVERSAL"],274,
 True,14];
(* Source data-flow audit. Dependence includes control flow and attached
   gates, and is not tested by inspecting a structural mark or cancellation.
   Substitution to a zero coefficient therefore never removes a site. *)
siteUses[knife_,tag_,path_] := Module[{p=ToString[path,InputForm],flow,optic,radar,gate,hit=False,row},
 flow=MemberQ[{"K1","K2","K3","K4c","K11"},knife];
 optic=flow||MemberQ[{"K4a","K4b"},knife];radar=optic||MemberQ[{"K5","K6"},knife];
 gate=flow||knife==="K10";
 If[StringContainsQ[tag,"LOCAL_S9B_"],Return[Which[
  StringContainsQ[tag,"SUPPLIED_INPUTS"],
    (MemberQ[path,"Dispersion"]&&MemberQ[{"K1","K2","K3","K11"},knife]) ||
    (MemberQ[path,"Metric"]&&knife==="K4c") ||
    (MemberQ[path,"Responses"]&&((knife==="K9a"&&MemberQ[path,"FIXED_RATIO"])||(knife==="K9b"&&MemberQ[path,"POWER"]))),
  StringContainsQ[tag,"PROFILE_ANSATZ"],knife==="K4c"&&MemberQ[path,"XiDerivative"],
  Or@@(StringContainsQ[tag,#]&/@{"BRANCH_FREQUENCIES","HAMILTON_VELOCITY","FERMAT_OBJECT"}),MemberQ[{"K1","K2","K3","K11"},knife],
  StringContainsQ[tag,"BRANCH_TYPES"],MemberQ[{"K1","K2","K3","K10","K11"},knife],
  StringContainsQ[tag,"NONRECIPROCAL_CONSTRUCTION"],gate||(knife==="K5"&&!MemberQ[path,"Domain"]),
  StringContainsQ[tag,"OPTICAL_COORDINATE"],If[MemberQ[path,"PushedDensity"],optic,flow],
  StringContainsQ[tag,"RAY_DOMAIN"],flow,
  StringContainsQ[tag,"RADAR_SLOPE_DERIVATION"],radar&&!MemberQ[path,"BoundedEndpointSymbols"],
  StringContainsQ[tag,"ABEL_REDUCTION"],(MemberQ[path,"DeflectionKernel"]&&optic)||(MemberQ[path,"RadarKernel"]&&radar),
  StringContainsQ[tag,"MASS_BALANCE"],knife==="K8"&&!MemberQ[path,"Measure"],
  StringContainsQ[tag,"C_RESPONSES"],(knife==="K9a"&&MemberQ[path,"FIXED_RATIO"])||(knife==="K9b"&&MemberQ[path,"POWER"]),
  True,False]]];
 If[StringStartsQ[tag,"WL_S9B_F_"],hit=knife==="K8"];
 If[MemberQ[path,"Premise"],Return[knife==="K8"&&MemberQ[path,"Flux"]]];
 If[StringContainsQ[tag,"F_MASS_SOLUTION"],Return[knife==="K8"]];
 If[StringContainsQ[tag,"_C_"],
   row=If[StringContainsQ[tag,"FIXED_RATIO"],"K9a",If[StringContainsQ[tag,"POWER"],"K9b",""]];
   Return[hit||gate||If[StringContainsQ[tag,"RADAR"]||StringContainsQ[tag,"N_DEPENDENCE"],radar,optic]||knife==="K7"||
     (knife==="K8"&&StringContainsQ[tag,"IMPLIED_JN"])||knife===row]];
 If[StringContainsQ[tag,"_B_"],Return[hit||gate||
   If[StringContainsQ[tag,"RADAR"]||StringContainsQ[tag,"DIFFERENCE"],radar,optic]||
   (knife==="K7"&&(StringContainsQ[tag,"CONDITION"]||StringContainsQ[tag,"IMPLIED_JN"]))||
   (knife==="K8"&&StringContainsQ[tag,"IMPLIED_JN"])]];
 If[StringContainsQ[tag,"_A_"],Return[hit||gate||Which[
   StringContainsQ[tag,"NONRECIPROCAL_PATH_DEPENDENCE"],flow||knife==="K5",
   StringContainsQ[tag,"RADAR_LOG"],radar,
   StringContainsQ[tag,"ROUND_TRIP"]||StringContainsQ[tag,"ONE_WAY_RE"]||StringContainsQ[tag,"NONRECIPROCAL"],optic||knife==="K5",
   True,optic]]];
 hit||If[StringContainsQ[tag,"PATH_TRAVERSAL"],flow,gate]];

extractedNameKeys={"FreeInputs","RemainingInputs"};
conditionKeys={"Domain","BranchDomain","AlongRayDomain","ConstantDomain","GeneralConstantDomain","SuppliedBulkDensityPremise"};
(* Fields with defining mathematical text retain that reading; they are not
   tried by literal string comparison. Counting/Measure are descriptions. *)
wordedKeys={"FunctionalDomain","ZeroSetJoining","AtComponentBoundary"};
textSourceDifferenceQ[knife_,left_,right_] := StringQ[left]&&StringQ[right]&&left=!=right&&
 StringContainsQ[knives[knife][[1]],ToString[left,InputForm]]&&
 StringContainsQ[knives[knife][[2]],ToString[right,InputForm]] ||
 (StringQ[left]&&StringQ[right]&&left=!=right&&
 StringContainsQ[knives[knife][[2]],ToString[left,InputForm]]&&StringContainsQ[knives[knife][[1]],ToString[right,InputForm]]);
partID[tag_,path_,role_] := ToString[{tag,path,role},InputForm,PageWidth->Infinity];
branchIdentity[branch_Association] := Which[
 KeyExistsQ[branch,"GM"],"unique linear-GM solution (linearRules; engine 226-228,405-425)",
 KeyExistsQ[branch,"Component"],{branch["Component"],If[TrueQ[Lookup[branch,"Velocity",None]===0],
    "zero root (Solve vRoot^2==0, engine 488)",
    If[MatchQ[Lookup[branch,"Velocity",None],Times[-1,_]],
      "negative square-root solution of vRoot^2==yRoot (engine 487,533-538)",
      "positive square-root solution of vRoot^2==yRoot (engine 487,533-538)"]]},
 KeyExistsQ[branch,"Velocity"],"unique prescribed or unconstrained velocity branch (engine 494-516)",
 KeyExistsQ[branch,"HomogeneousSolution"],"excluded nonzero integration-constant branch (engine 529-531)",
 True,{"branch defined by domain",Lookup[branch,"Domain",branch]}];
(* The map contains expressions for evaluation but the printed inventory
   refers back to the already printed emitted payloads, avoiding reprinting
   megabytes of the same attached ray condition at each child. *)
catalog[data_] := Module[{entries=<||>,containers={},walk,add,tag},
 add[z_,path_,role_,conditions_,reading_] := AssociateTo[entries,partID[tag,path,role]-><|
   "Tag"->tag,"Location"->path,"Role"->role,"Object"->z,
   "AttachedConditionLocations"->conditions,"Reading"->reading,"SourceLine"->sourceAnchor[tag]|>];
 walk[z_,path_,conditions_] := Module[{keys,ownConditions,conditionPath,rows,earlier={},cond,effective,branch,hasNames},
  Which[
   AssociationQ[z],
    keys=Keys[z];ownConditions=Join[conditions,({tag,Append[path,#]}&/@Intersection[keys,conditionKeys])];
    If[KeyExistsQ[z,"Premise"],AppendTo[ownConditions,{tag,Append[path,"Premise"]}]];
    AppendTo[containers,<|"Tag"->tag,"Location"->path,"Reading"->"association of named objects; conditions attach to children at every containing level","Conditions"->ownConditions|>];
    KeyValueMap[Function[{key,value},walk[value,Append[path,key],DeleteCases[ownConditions,{tag,Append[path,key]}]]],z],
   ListQ[z]&&Length[path]>0&&MemberQ[{"Responses","DeltaSeries","Domains"},Last[path]],
    AppendTo[containers,<|"Tag"->tag,"Location"->path,"Reading"->"response rows identified by CONSTANT,FIXED_RATIO,POWER; engine Table i=1..3, responseNames 611,624-625; positions only implement that named mapping"|>];
    MapThread[walk[#1,Append[path,#2],conditions]&,{z,Take[rowNames,Length[z]]}],
   ListQ[z]&&Length[path]>0&&Last[path]==="DensityParts",
    AppendTo[containers,<|"Tag"->tag,"Location"->path,"Reading"->"same shared Object in each density sign; LabelCondition partitions it, is not an attached condition; 618-623"|>];
    Scan[Function[row,walk[row,Append[path,row["Label"]],conditions]],z],
   ListQ[z]&&And@@(AssociationQ/@z)&&Length[path]>0&&MemberQ[{"Cases","ExcludedBranches"},Last[path]],
    AppendTo[containers,<|"Tag"->tag,"Location"->path,"Reading"->"branch set, not first-match; identities are the defining domain/component/solution, or the constructor's unique linear-GM branch; engine 405-425,529-553,564-577"|>];
    Scan[Function[row,walk[row,Append[path,branchIdentity[row]],conditions]],z],
   Head[z]===ConditionalExpression,
    AppendTo[containers,<|"Tag"->tag,"Location"->path,"Reading"->"conditional value and its domain are separate site parts; domain attaches to value; engine 102,271-273"|>];
    add[z[[2]],Append[path,"conditional domain"],"domain",conditions,"ConditionalExpression second operand" ];
    walk[z[[1]],Append[path,"conditional value"],Append[conditions,{tag,Append[path,"conditional domain"]}]],
   Head[z]===Piecewise,
    AppendTo[containers,<|"Tag"->tag,"Location"->path,"Reading"->"ordered first-match; effective branch predicate includes failure of every earlier predicate; default requires all to fail","OriginalConditions"->(Last/@z[[1]])|>];
    Scan[Function[row,cond=row[[2]];effective=And@@Join[Not/@earlier,{cond}];
      branch={"case selected by",cond};
      add[effective,Append[path,{branch,"selection"}],"domain",conditions,"ordered first-match effective predicate"];
      walk[row[[1]],Append[path,branch],Append[conditions,{tag,Append[path,{branch,"selection"}]}]];
      AppendTo[earlier,cond]],z[[1]]];
    add[And@@(Not/@earlier),Append[path,"default selection"],"domain",conditions,"all earlier cases fail"];
    walk[z[[2]],Append[path,"default"],Append[conditions,{tag,Append[path,"default selection"]}]],
   True,
    add[z,path,Which[Length[path]>0&&Last[path]==="Profiles","literal profile list",StringQ[z],If[AnyTrue[path,MemberQ[wordedKeys,#]&],"worded mathematical part","literal text"],
      AnyTrue[path,MemberQ[conditionKeys,#]&],"domain",True,"value"],conditions,
      Which[Length[path]>0&&Last[path]==="Profiles","literal descriptive list of profile symbols written directly at engine 401; not extracted from expressions, not a condition, and no profile-function substitution is applied to its elements",ListQ[z],"one mathematical vector/tensor, equation list or solution set; tensor basis and solution variables are those of its source constructor",
        StringQ[z],"engine literal status/label unless a defining worded mathematical condition is named",
        True,"one named mathematical object; nested Piecewise has ordered first-match semantics; arithmetic subexpressions are not exposed parts"]]]];
 Scan[Function[t,tag=t;walk[object[data,t],{},{}]],Keys[data["Held"]]];
 <|"Parts"->entries,"Containers"->containers|>];
(* K10 is introduced only at propagationGate (262-264). Direct fields
   are the emitted branch condition and ray-gate domains; a child also
   depends on it when one of those fields is an attached condition. *)
k10DirectSiteQ[tag_,path_] := Module[{p=ToString[path,InputForm]},
 If[StringContainsQ[tag,"LOCAL_S9B_"]||MemberQ[path,"Premise"],Return[False]];
 If[StringEndsQ[tag,"BRANCH_EXISTENCE"],Return[True]];
 If[StringEndsQ[tag,"PATH_TRAVERSAL"]||StringEndsQ[tag,"BRANCH_TYPE"],Return[False]];
 If[MemberQ[path,"BranchDomain"]||MemberQ[path,"AlongRayDomain"],Return[True]];
 MemberQ[path,"conditional domain"]&&(MemberQ[path,"OnDomain"]||MemberQ[path,"OutsideDomain"])];
k10PartUses[record_] := k10DirectSiteQ[record["Tag"],record["Location"]]||
 AnyTrue[record["AttachedConditionLocations"],k10DirectSiteQ[#[[1]],#[[2]]]&];

publicPart[record_] := KeyDrop[record,{"Object"}];
inventory[knife_,mutatedData_,baselineData_] := Module[{ac,bc,ids,included={},excluded={},av,bv,rec,why,reached},
 ac=catalog[mutatedData];bc=catalog[baselineData];ids=Union[Keys[ac["Parts"]],Keys[bc["Parts"]]];
 Do[av=Lookup[ac["Parts"],id,Missing["not emitted"]];bv=Lookup[bc["Parts"],id,Missing["not emitted"]];
   rec=If[AssociationQ[av],av,bv];why=None;
   If[MemberQ[rec["Location"],"Label"]||MemberQ[rec["Location"],"LabelCondition"]||
     MemberQ[rec["Location"],"SuppliedBulkDensityPremise"],
     why="the supplied density premise and physical sign partition are computed independently of the mutated construction (engine 111,618-623)"];
   If[why===None&&!If[knife==="K10",k10PartUses[rec],siteUses[knife,rec["Tag"],rec["Location"]]],why="computation and its attached conditions do not use the mutated construction; independent or before entry (source data-flow ranges printed)"];
   If[AnyTrue[rec["Location"],MemberQ[extractedNameKeys,#]&]||StringContainsQ[rec["Tag"],"ENGINE_LOCAL_NAMES"],
     why="additional exclusion 1: lists of symbol/function/row names extracted from expressions (and engine local name registry)"];
   If[AssociationQ[av]&&AssociationQ[bv]&&textSourceDifferenceQ[knife,av["Object"],bv["Object"]],
     why="additional exclusion 2: literal text difference is the knife replacement source"];
   reached=If[knife==="K11",KeyExistsQ[compactCorrupted["Held"],rec["Tag"]]&&KeyExistsQ[compactBaseline["Held"],rec["Tag"]],True];
   rec=Join[publicPart[rec],<|"ID"->id,"Copies"->{If[AssociationQ[av],"emitted","not emitted"],If[AssociationQ[bv],"emitted","not emitted"]},
     "CompactReached"->reached,"LocationIdentity"->"engine object name/response row/density sign/defining branch, never a storage position"|>];
   If[why===None,AppendTo[included,rec],AppendTo[excluded,Append[rec,"Reason"->why]]],{id,ids}];
 print[knife<>"_SITE_INVENTORY",<|"OperandOrder"->{"corrupted","unmutated"},"SiteParts"->included,
   "Exclusions"->excluded,"ContainersAndReadings"->{ac["Containers"],bc["Containers"]},
   "SourceDataFlowRanges"->sourceRanges,
   "AbsentEntryReading"->"not emitted denotes only missing storage; engine 405-425 and 490-577 have different returned fields but define no value for missing fields; absence alone proves nothing",
   "K11Reading"->"inventory includes full-baseline downstream parts; compact-unreached corrupted parts are not evaluated, and no no-change decision is made",
   "ConditionsReading"->"Domain/BranchDomain/AlongRayDomain/ConstantDomain/premises and conditional/case selections attach at every level. Density sign labels partition. Profiles/Measure/Closedness/OpticalCounting are not conditions"|>];
 <|"Parts"->included,"CorruptedCatalog"->ac,"BaselineCatalog"->bc|>];
textAttempt[knife_,mutatedData_,baselineData_] := Module[{tag="WL_S9B_B_IMPLIED_JN_DEFLECTION",texts,decision},
 texts=Map[object[#,tag]["VelocityStatus"]&,{mutatedData,baselineData}];
 decision=SameQ@@texts;
 <|"Result"->If[decision,"no change","witness"],"OperandOrder"->{"corrupted","baseline"},
   "Site"-><|"Tag"->tag,"PayloadLocation"->{"VelocityStatus"},"SitePart"->"literal status text","Source"->engine,"SourceLines"->{498,544,575}|>,
   "ReadingChoices"->"the engine copies solution[Status] into VelocityStatus (575); it supplies no mathematical parser for this literal status. The associated equations are separate parts. It is neither a name list nor the replacement text EVERY_B/FIXED_B",
   "Texts"->texts,"CASOperation"->Inactive[SameQ]@@texts,"LiteralEqualityDecision"->decision,
   "Evaluations"->texts,"Checks"-><|"EngineLiteralText"->And@@(StringQ/@texts),"LiteralDifference"->Not[decision]|>,
   "Attempt"->"Wolfram SameQ comparison of both emitted literal status texts; mathematical point/solution tests are replaced by the item-14 text rule"|>];
(* END_AMENDMENT12_CATALOG *)

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
 "K10"->{"ray",{"WL_S9B_BRANCH_EXISTENCE"},"propagationGate = rayConstruction"},
 "K11"->{"local",{"WL_LOCAL_S9B_FERMAT_OBJECT","TimePolynomial"},"timePolynomial = Expand["}|>;
primaryAttempt[name_,a_,b_] := Module[{spec=specs[name],kind,path,aa,bb,conditions,reading,
  evaluatedConditions,values,roles,symbols,fixed,assigned,unknown,entry,ceConditions={True,True},
  density,checks,location,detail,cases},
  {kind,path}=Take[spec,2];note["witness begin",name];
  aa=part[a,path];bb=part[b,path];
  location=<|"Tag"->First[path],"PayloadLocation"->semanticPaths[name],"SitePart"->"value",
    "Source"->engine,"SourceLine"->sourceLine[Last[spec]],"EntryConstruction"->Last[spec]|>;
  If[kind==="ray",entry=searchReachedSite[name,a,b];
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
    note["witness end",name];Return[<|"Result"->"undecided","Site"->location,
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
  entry=<|"Result"->If[And@@(TrueQ/@Values[checks]),"witness","undecided"],
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

(* Finite search order is printed before the first attempt. The complete union
   is enumerated first, including unreached K11 parts and both density signs.
   A witness ends the search. A failed point never becomes a no-change decision.
   There are no mathematical no-change/cannot-be-tried shortcuts in this search. *)
searchWitness[name_,mutatedData_,baselineData_] := Module[{inv,attempt,negative,positive,result,remaining,mathRecords={},record,aa,bb},
 inv=inventory[name,If[name==="K11",mutatedData,mutatedData],If[name==="K11",baseline,baselineData]];
 print[name<>"_SEARCH_PLAN",<|"Order"->If[name==="K7",{"literal VelocityStatus", "all remaining semantic site parts in printed ID order if no witness"},
   If[MemberQ[{"K9a","K9b"},name],{"one combined trial of the shared response object on rho <= 0 and rho > 0, at both printed concrete inputs","remaining parts in printed ID order if no witness"},
   {semanticPaths[name],"remaining semantic site parts in printed ID order if no witness"}]],
   "AttemptBound"->"one candidate configuration per semantic part; K9's first configuration includes both required sign trials; K11 compact-unreached parts remain not evaluated",
   "GeneralInputCoverage"->"no mathematical no-change or cannot-be-tried decision is inferred from a point; all live profiles and inputs remain admitted outside the attempted points"|>];
 Which[
  name==="K7",attempt=textAttempt[name,mutatedData,baselineData],
  MemberQ[{"K9a","K9b"},name],
    negative=Block[{certificate,profiles=Join[profiles,<|Global`fProfile->Function[{q},-9/q]|>],
       point=Join[DeleteCases[point,Rule[Global`n,_]],{Global`n->3}]},primaryAttempt[name,mutatedData,baselineData]];
    negative=Join[negative,<|"DensitySignTrial"->"rho <= 0","Site"->Join[negative["Site"],
      <|"PayloadLocation"->(semanticPaths[name]/."rho > 0"->"rho <= 0")|>],
      "DensityEvaluationsEachEmittingCopy"->Table[negative["DensityPart"]["DensityEvaluation"],{2}],
      "SignDecision"->TrueQ[negative["DensityPart"]["DensityEvaluation"]<=0],
      "AsymptoticDensityComputation"->Inactive[Limit][Global`rho0(1-9/Global`r),Global`r->Infinity],
      "AsymptoticDensityEvaluation"->Limit[1-9/Global`r,Global`r->Infinity]|>];
    positive=primaryAttempt[name,mutatedData,baselineData];
    positive=Join[positive,<|"DensitySignTrial"->"rho > 0",
      "DensityEvaluationsEachEmittingCopy"->Table[evaluate[Global`rho0(1+Global`fProfile[Global`r])],{2}],
      "SignDecision"->TrueQ[evaluate[Global`rho0(1+Global`fProfile[Global`r])]>0],
      "AsymptoticDensityComputation"->Inactive[Limit][Global`rho0(1+1/(100 Global`r)),Global`r->Infinity],
      "AsymptoticDensityEvaluation"->Limit[1+1/(100 Global`r),Global`r->Infinity]|>];
    (* This is one combined two-sign trial, so both computations are printed
       before returning the first witness from that trial. *)
    print[name<>"_DENSITY_SIGN_TRIALS",{negative,positive}];
    attempt=If[negative["Result"]==="witness",negative,positive],
  True,attempt=primaryAttempt[name,mutatedData,baselineData]];
 print[name<>"_ATTEMPT",attempt];
 If[attempt["Result"]==="witness",Return[Join[attempt,<|"InventoryCount"->Length[inv["Parts"]],
   "SearchTermination"->"first witness in the printed finite attempt order; remaining parts were enumerated, not tried and have no equality decision"|>]]];
 (* If a primary configuration is unresolved, make a finite exact attempt at
    every remaining emitted part. Evaluation failures remain undecided.
    These fallbacks are not mathematical universal decisions. *)
 Do[record=site;aa=Lookup[inv["CorruptedCatalog"]["Parts"],record["ID"],Missing["not emitted"]];
   bb=Lookup[inv["BaselineCatalog"]["Parts"],record["ID"],Missing["not emitted"]];
   result=If[name==="K11"&&!TrueQ[record["CompactReached"]],
     <|"Evaluations"->{"not evaluated","not evaluated"},"UnmetRequirements"->{"compact method does not reach this part"}|>,
     <|"Inputs"->point,"Profiles"->Normal[profiles],"CASOperations"->"scoped profile/point substitution, ordinary Integrate, active Resolve and FullSimplify; no universal agreement inferred",
       "Evaluations"->MapThread[If[AssociationQ[#1],If[#1["Role"]==="literal profile list",<|"Evaluation"->#1["Object"],"Computation"->"evaluate the source-written list of symbols literally; no profile replacement; not an extracted-name or text comparison"|>,siteEvaluate[#1["Object"],#2]],"not emitted"]&,
         {{aa,bb},{mutatedData,baselineData}}],
       "UnmetRequirements"->{"no decision of agreement over all admissible inputs; attached-condition and construction solution-set coverage not established by this fallback"}|>];
   result=Join[publicPart[record],result,<|"Result"->"undecided"|>];
   print[name<>"_REMAINING_ATTEMPT",result];AppendTo[mathRecords,KeyTake[result,{"Tag","Location","Role","UnmetRequirements"}]],{site,inv["Parts"]}];
 <|"Result"->"undecided","PrimaryAttempt"->attempt,"PartsWithoutNoChangeDecision"->mathRecords,
   "NoneFoundRule"->"none found requires CAS no-change decisions at every site part; no such universal decisions were obtained"|>];
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
