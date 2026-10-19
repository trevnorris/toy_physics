(* Native comparison instrument. The only literal physical combinations are
   in the ORACLE TRANSCRIPTION block. Other expressions are source reads,
   algebra on those reads, or explicitly printed identification substitutions.
   Run via the Python driver in the pooled guard; never directly. *)
$HistoryLength = 0;
(* Keep internal algebra threads within the pooled guard's single-CPU task cap.
   Wolfram's internal thread controls are separate from OMP/MKL environment vars.
   https://reference.wolfram.com/language/Compile/tutorial/Parallel.html *)
SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
root = "/var/projects/toy_physics/_scratch/s9b_oracle_check/";
data = Import[root <> "inputs.json", "RawJSON"];
emit[tag_, obj_] := Print[tag, ": ", ToString[obj, InputForm, PageWidth -> Infinity]];
progress[obj_] := emit["WORKING_ON", obj];
quote[key_] := emit["ORACLE_QUOTED_" <> key, data["quotes"][key]];
checks = {};
check[label_, value_] := AppendTo[checks, {label, value}];
parse[x_] := ToExpression[x, InputForm];
pget[tag_] := parse[data["py"][tag]["wl"]];
wget[tag_] := parse[data["wl"][tag]["raw"]];
node[x_] := parse[x["wl"]];
(* No heuristic Simplify search: explicit rational normalization followed by
   Refine under the operand's source domain. All residual calls go through report. *)
domainSimplify[x_, ass_] := Refine[Cancel[Together[x]], ass];
domainAtoms=<||>;
opaqueAtom[x_] := Module[{id="sourcePredicate$"<>ToString[Hash[x,"SHA256"]],sym},
 sym=Symbol[id];
 If[!KeyExistsQ[domainAtoms,id],AssociateTo[domainAtoms,id->x];emit["SOURCE_DOMAIN_PREDICATE_DEFINITION",{sym,x}]];
 sym];
opaqueQ[x_] := !FreeQ[x,sourceConditionSet|sourceInterval|sourceProductSet|sourceIntersection|sourceUnion|sourceContains|Inactive[Integrate]|Inactive[Limit]|Inactive[ForAll]];
domainLogic[x_And] := And@@(domainLogic/@(List@@x));
domainLogic[x_Or] := Or@@(domainLogic/@(List@@x));
domainLogic[x_] := If[opaqueQ[x],opaqueAtom[x],x];
makeDomain[id_,parts_Association] := Module[{full=And@@Values[parts],cas},
 emit["OPERAND_DOMAIN",<|"id"->id,"source_parts"->parts|>];
 (* Nonlocal ray existence/traversal predicates retain their exact source
    definitions as logical premises; the local residual does not solve those
    functional/integral existence problems afresh. *)
 cas=And@@KeyValueMap[If[MemberQ[{"ray_domain","along_ray"},#1],opaqueAtom[#2],domainLogic[#2]]&,parts];
 <|"id"->id,"CAS"->cas|>];
report[tag_,a_,b_,dom_] := Module[{ans},
 progress[{tag,dom["id"]}];
 emit["RESIDUAL_OPERANDS",{a,b}];
 emit["RESIDUAL_ASSUMPTIONS",dom];
 ans=domainSimplify[a-b,dom["CAS"]];emit[tag,ans];ans];
activateOperand[x_, ass_] := Block[{$Assumptions=ass}, domainSimplify[Activate[x, Integrate], ass]];
quoteMany[keys_] := Scan[quote, keys];
emit["KERNEL_THREAD_OPTIONS",{SystemOptions["ParallelOptions"],SystemOptions["InequalitySolvingOptions" -> "ThreadOr"]}];
emit["PROVENANCE", data["provenance"]];
emit["SUPPLIED_ORDER", <|"S9b" -> "Exterior flow-only restriction; Part B retained first optical order: velocity through its squared contribution, with supplied half-order velocity counting. General integration constant remains live in a separately printed, counting-excluded family.",
 "Optics" -> "Supplied leading-eikonal dispersion, not an S9b result; coefficients of the quadratic characteristic polynomial only. No amplitude transport, finite-wavelength terms, wave equation, force law or material law is compared.",
 "Oracle" -> "Exact quoted spherical profiles and cone polynomials restricted to the common exterior patch. Characteristic polynomials are quadratic in the flow; no higher optical order is inferred from exact metric inverses. No rotating, core, bubble-interior or horizon extension of S9b is taken."|>];
Scan[emit["IDENTIFICATION_MAP", #]&, data["map"]];
quoteMany[Keys[data["quotes"]]];

(* ORACLE TRANSCRIPTION: symbols refer only to the quoted oracle operands.
   gmH and gmV denote each paper's written G M, separately from S9b's GM.
   signO encodes the two printed +/- alternatives. norm53 records the unspecified
   proportionality in (53); norm55 records the omitted constant in (55). *)
hlLine = -odt^2 + (odr + betaO odt)^2 + orad^2 (odtheta^2 + Sin[otheta]^2 odphi^2);
hlProfile = betaO == Sqrt[2 gmH/orad];
hlMassProfile = betaO == Sqrt[2 massH[orad]/orad];
visMetric = (densityO/soundO) {{-(soundO^2-vxO^2-vyO^2-vzO^2),-vxO,-vyO,-vzO},
 {-vxO,1,0,0},{-vyO,0,1,0},{-vzO,0,0,1}};
visContinuity = HoldComplete[Inactive[D][densityField[ot,ox,oy,oz],ot] +
 Inactive[Div][densityField[ot,ox,oy,oz] {velocityX[ot,ox,oy,oz],velocityY[ot,ox,oy,oz],velocityZ[ot,ox,oy,oz]}, {ox,oy,oz}] == 0];
visEuler = HoldComplete[densityO materialAccelerationO == densityO (partialTimeVelocityO + advectiveDerivativeO) == forceO];
visForce = HoldComplete[forceO == -pressureGradientO-densityO gravityGradientO-densityO drivingGradientO];
visBackground7 = {speed7O == Sqrt[2 gmV/orad], soundO == 1,
 density7O == densityNormO orad^(-3/2), pressureO == pressureInfinityO + soundO^2 density7O};
visLine53 = norm53 orad^(-3/2) (-odt^2 + (odr + signO Sqrt[2 gmV/orad] odt)^2 + orad^2 (odtheta^2+Sin[otheta]^2 odphi^2));
visProfile54 = speed8O == soundO radius0O^2/orad^2;
visLine55 = -soundO^2 odt^2 + (odr + signO soundO radius0O^2/orad^2 odt)^2 + orad^2 (odtheta^2+Sin[otheta]^2 odphi^2);
visFactor4 = densityO/soundO;
visFactor53 = norm53 orad^(-3/2);
(* END ORACLE TRANSCRIPTION. norm55 is the constant multiplier omitted by (55). *)
emit["ORACLE_TRANSCRIBED_OBJECTS", <|"HL1"->hlLine,"HL2"->hlProfile,"HL4"->hlMassProfile,
 "V4"->visMetric,"V5"->visContinuity,"V6"->visEuler,"V7"->visForce,
 "V7background"->visBackground7,"V53"->visLine53,"V54"->visProfile54,"V55"->visLine55|>];

emit["CONTROL_VARIANT",data["control"]];
If[data["control"]=="visser_form",
 visContinuity=visContinuity/.Inactive[Div][Times[scalar_,vector_List],coords_List]:>Times[scalar,Inactive[Div][vector,coords]];
 emit["CONTROL_CHANGED_ORACLE_OPERATOR",visContinuity]];
If[data["control"]=="hl_spatial_form",
 hlLine=Plus[hlLine,Times[Power[hlProfile[[2]],2],Power[differentialPlaceholder,2]]/.differentialPlaceholder->odr];
 emit["CONTROL_CHANGED_HL_LINE_ELEMENT",hlLine]];

(* Hessian extraction, orthonormal spatial basis, inverse-metric eikonal
   contraction and frequency-coefficient normalization are computations. *)
differentials={odt,odr,odtheta,odphi};
lineMetric[line_] := Table[D[line,differentials[[i]],differentials[[j]]]/2,{i,4},{j,4}];
localMetric[metric_,factor_] := Module[{spatial,basis},
 spatial=Cancel[metric[[2;;4,2;;4]]/factor];
 basis=DiagonalMatrix[Join[{1,1},1/Sqrt[Diagonal[spatial][[2;;3]]]]];
 Refine[Cancel[Transpose[basis].metric.basis],orad>0 && 0<otheta<Pi]];
momenta={frequencyQ,radialQ,tangentQ,otherQ};
covector=DiagonalMatrix[{-1,1,1,1}].momenta;
cone[metric_] := Cancel[covector.Inverse[metric].covector];
normalize[polynomial_] := Cancel[polynomial/Coefficient[Expand[polynomial],frequencyQ,2]];
drift[metric_] := Cancel[-metric[[1,2]]/metric[[2,2]]];
clockPullback[metric_,rate_] := Module[{jac=IdentityMatrix[Length[metric]]},jac[[1,1]]=rate;Cancel[Transpose[jac].metric.jac]];
metricSpeedSquared[metric_] := Cancel[-Subtract[metric[[1,1]],metric[[1,2;;4]].Inverse[metric[[2;;4,2;;4]]].metric[[2;;4,1]]]/metric[[2,2]]];

baseHL=localMetric[lineMetric[hlLine],1];
base53=localMetric[lineMetric[visLine53],visFactor53];
base55=localMetric[lineMetric[visLine55],1];
models={
 <|"name"->"HL2","quotes"->{"HL1","HL2"},"metric"->(baseHL/.Apply[Rule,hlProfile]),"factor"->1,"mass"->gmH,"profile"->hlProfile[[2]],"profile_equation"->"HL2"|>,
 <|"name"->"HL4","quotes"->{"HL1","HL4"},"metric"->(baseHL/.Apply[Rule,hlMassProfile]),"factor"->1,"mass"->massH[orad],"profile"->hlMassProfile[[2]],"profile_equation"->"HL4"|>,
 <|"name"->"V7_plus","quotes"->{"V7background","V53","V86","V87"},"metric"->(base53/.signO->1),"factor"->visFactor53,"mass"->gmV,"profile"->visBackground7[[1,2]],"profile_equation"->"V7background"|>,
 <|"name"->"V7_minus","quotes"->{"V7background","V53","V86","V87"},"metric"->(base53/.signO->-1),"factor"->visFactor53,"mass"->gmV,"profile"->visBackground7[[1,2]],"profile_equation"->"V7background"|>,
 <|"name"->"V8_plus","quotes"->{"V54","V55","V86","V87"},"metric"->(norm55 base55/.signO->1),"factor"->norm55,"mass"->Missing["NoMassParameter"],"profile"->visProfile54[[2]],"profile_equation"->"V54"|>,
 <|"name"->"V8_minus","quotes"->{"V54","V55","V86","V87"},"metric"->(norm55 base55/.signO->-1),"factor"->norm55,"mass"->Missing["NoMassParameter"],"profile"->visProfile54[[2]],"profile_equation"->"V54"|>};
AppendTo[models,<|"name"->"V4","quotes"->{"V4","V85","V86","V87"},"metric"->visMetric,"factor"->visFactor4,"mass"->Missing["NoMassParameter"]|>];

pyProfiles=pget["PY_LOCAL_S9B_GENERAL_PROFILES"];
pySpeed=pget["PY_LOCAL_S9B_SPEED_IDENTIFICATION"];
wlSupplied=wget["WL_LOCAL_S9B_SUPPLIED_INPUTS"];
wlProfiles=wget["WL_LOCAL_S9B_PROFILE_ANSATZ"];
wlRestriction=wget["WL_LOCAL_S9B_RESTRICTION_FLOW_ONLY"];
wlSample=wget["WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION"];
(* This speed relation is obtained from the attached flow-only domain, not
   from rewriting the supplied dispersion with the profile restriction. *)
wlSpeedRelations=DeleteDuplicates[Cases[wlSample, e_Equal /; MemberQ[List@@e,c0] && !FreeQ[e,muPerp] :> e,Infinity]];
emit["WL_FLOW_SPEED_RELATIONS_FROM_ATTACHED_DOMAIN",wlSpeedRelations];
wlSpeedSquareRules=If[Length[wlSpeedRelations]>0,{Apply[Rule,Map[#^2&,First[wlSpeedRelations]]]},{}];
pyRestrict=Thread[pyProfiles[[{1,3}]]->0];
pyC0=py`c$0;
engineData={
 <|"name"->"PY","radius"->py`s9b$r,"c0"->pyC0,"mass"->py`GM,
 "velocity"->pyProfiles[[2]],"density"->pyProfiles[[4]],
 "speedSquared"->(pySpeed[[1,1]]/.pyRestrict),"restriction"->pyRestrict,
 "dispersion"->pget["PY_LOCAL_S9B_VECTOR_DISPERSION"],
 "waveMap"->Thread[{py`s9b$omega,py`s9b$k$r,py`s9b$k$t,py`s9b$k$z}->momenta],
 "openBindings"->{py`s9b$chi->(pySpeed[[1,1]]/.pyRestrict),py`s9b$A->lineMetric[hlLine][[2,2]],py`s9b$U->pyProfiles[[2]]},"fixedBindings"->{}|>,
 <|"name"->"WL","radius"->r,"c0"->c0,"mass"->GM,
 "velocity"->wlProfiles["V"],"density"->wlProfiles["RhoBr"],
 "speedSquared"->(wlSupplied["Speed"][[2]]/.Reverse/@wlSpeedSquareRules),"restriction"->wlRestriction,
 "dispersion"->Subtract@@(List@@wlSupplied["Dispersion"]),
 "waveMap"->Join[Thread[{omega,kr}->{frequencyQ,radialQ}],{kphi->Times[tangentQ,Sqrt[wlSupplied["Metric"][[2,2]]]]}],
 "openBindings"->{aMetric->(wlSupplied["Metric"][[1,1]]/.wlRestriction),vLocal->wlProfiles["V"]},
 "fixedBindings"->Join[{Apply[Rule,wlSupplied["Speed"]]},Reverse/@wlSpeedSquareRules]|>};
(* All source equalities above are converted directly, preserving their sides. *)

(* Source declarations are read per package, including all PY srepr assumptions.
   Opaque source set/quantifier predicates are losslessly named Boolean atoms;
   their complete definitions are printed. None is set True or dropped. *)
pyDeclared[tag_] := And@@(node/@data["py"][tag]["declared"]);
pyCommonDeclared=And@@(pyDeclared/@{"PY_LOCAL_S9B_GENERAL_PROFILES","PY_LOCAL_S9B_COUNTING","PY_LOCAL_S9B_SPEED_IDENTIFICATION","PY_LOCAL_S9B_VECTOR_DISPERSION"});
engineData[[1]]=Join[engineData[[1]],<|"DeclaredDomain"->pyCommonDeclared,"balance"->pget["PY_LOCAL_S9B_MASS_BALANCE"]|>];
engineData[[2]]=Join[engineData[[2]],<|"DeclaredDomain"->wlProfiles["Domain"],"balance"->wlSupplied["Balance"]|>];
packages=<||>;
prepare[eng_,route_] := Module[{tag,obj,meta,declared,attached,common,general,squared,cases,dom,roots,constantDomain,rayDomain,caseDomain},
 tag=eng["name"]<>"_S9B_R_FLOW_ONLY_B_IMPLIED_JN_"<>route;
 progress[{tag,"source domains"}];
 If[eng["name"]=="PY",
  obj=data["py"][tag];emit["S9B_RAW_SOURCE",obj];
  declared=And[eng["DeclaredDomain"],pyDeclared[tag]];
  attached=And@@(node/@Most[obj["attached"]]);rayDomain=node[Last[obj["attached"]]];
  constantDomain=True; (* integration constant's real declaration is already read above *)
  general=node[obj["general"]][[2]];squared=node[obj["squared"]];
  cases=Map[Function[c,<|"Velocity"->node[c["velocity"]][[2]],"ImpliedExchange"->node[c["exchange"]],
    "VelocityDerivative"->If[c["velocity_derivative"]===Null,D[node[c["velocity"]][[2]],eng["radius"]],node[c["velocity_derivative"]][[2]]],
    "EvaluationPoint"->node[c["velocity"]][[1,1]],"Domain"->node[c["domain"]],"FullCase"->c["raw"]|>],obj["cases"]];
  emit["S9B_CONSTANT_METADATA",obj["constants"]],
  obj=wget[tag];emit["S9B_RAW_SOURCE",data["wl"][tag]];
  declared=eng["DeclaredDomain"];
  meta=First[Cases[obj,a_Association /; KeyExistsQ[a,"GeneralSquaredVelocity"],Infinity]];
  attached=And[obj["ReducedProfileCondition"]["Domain"]["Radius"],obj["ReducedProfileCondition"]["Domain"]["Anchor"]];rayDomain=obj["ReducedProfileCondition"]["BranchDomain"];
  constantDomain=meta["GeneralConstantDomain"];
  general=meta["GeneralSquaredVelocity"];squared=meta["SquaredVelocity"];cases=obj["Cases"];
  emit["S9B_VELOCITY_METADATA",meta]];
 common=makeDomain[tag<>"/general",<|"declared"->declared,"attached"->attached,"ray_domain"->rayDomain,"constant"->constantDomain|>];
 general=activateOperand[general,common["CAS"]];squared=activateOperand[squared,common["CAS"]];
 roots=unknownRoot/.Solve[Power[unknownRoot,2]==general,unknownRoot];
 emit["S9B_GENERAL_SQUARED",general];emit["S9B_RETAINED_SQUARED",squared];emit["S9B_GENERAL_SIGNED_ROOTS",roots];
 emit["GENERAL_COMPONENT_DOMAIN",<|"positive_components"->Greater[general,0],"zero_locus"->Equal[general,0],"real_domain"->GreaterEqual[general,0],"qualification"->"Formal general roots retain the printed real-velocity declaration and attached domain. Positive components and differentiable zero-set joining have the source qualifications; the supplied counting excludes the nonzero integration constant."|>];
 cases=MapIndexed[Function[{cc,index},
  caseDomain=If[eng["name"]=="WL",Block[{$Assumptions=declared && attached},Activate[cc["Domain"],Integrate]],cc["Domain"]];
  emit["SOURCE_CASE_DOMAIN_EVALUATION",{cc["Domain"],caseDomain}];
  dom=makeDomain[tag<>"/case/"<>ToString[First[index]],<|"declared"->declared,"attached"->attached,
    "case"->caseDomain,"ray_domain"->rayDomain,"constant"->Lookup[cc,"ConstantDomain",True],"along_ray"->Lookup[cc,"AlongRayDomain",True]|>];
  Join[cc,<|"Assumptions"->dom,"EvaluationPoint"->Lookup[cc,"EvaluationPoint",eng["radius"]]|>]],cases];
 <|"tag"->tag,"general"->general,"squared"->squared,"roots"->roots,"cases"->cases,"domain"->common|>];

(* All source expressions, integrands and limits stay unchanged. Native algebraic
   equalities are used before evaluating boundary limits, on that same domain. *)
evaluateCase[x_,cc_,eng_] := Module[{eqs,rules,y},
 eqs=Cases[List@@cc["Domain"],e_Equal /; !FreeQ[e,eng["mass"]] && FreeQ[e,Inactive|sourceConditionSet|sourceContains] :> e,{1}];
 eqs=Select[eqs,FreeQ[#, _Derivative|_sourceConditionSet|_Inactive|_sourceContains]&];
 rules=If[FreeQ[x,Inactive[Limit]],{},If[eqs==={}, {},First[Solve[And@@eqs,eng["mass"]]]]];
 y=domainSimplify[x/.rules,cc["Assumptions"]["CAS"]];
 (* Constant limits are evaluated algebraically before calling Limit, so the
    kernel does not discard source assumptions on a dummy limit variable. *)
 y=y/.Inactive[Limit][value_,Rule[variable_,destination_],options___]/;FreeQ[value,variable]:>value;
 Block[{$Assumptions=cc["Assumptions"]["CAS"]},domainSimplify[Activate[y,Integrate|Limit],$Assumptions]]];

compareFlow[eng_,package_,oracle_,orientation_,map_] := Module[{rad=eng["radius"],dom,point,source},
 emit["FLOW_ORACLE_OPERAND",<|"orientation"->orientation,"value"->oracle|>];emit["FLOW_IDENTIFICATION_MAP",map];
 report["FLOW_SQUARED_RESIDUAL",package["general"],Power[oracle,2],package["domain"]];
 Do[report["FLOW_SIGNED_RESIDUAL",root,oracle,package["domain"]],{root,package["roots"]}];
 Do[point=cc["EvaluationPoint"];dom=cc["Assumptions"];
  source=evaluateCase[cc["Velocity"],cc,eng];
  emit["FLOW_S9B_CASE",cc];emit["FLOW_ORACLE_AT_COMPONENT_POINT",oracle/.rad->point];
  report["FLOW_SQUARED_RESIDUAL",package["squared"]/.rad->point,Power[oracle/.rad->point,2],dom];
  report["FLOW_SIGNED_RESIDUAL",source,oracle/.rad->point,dom],{cc,package["cases"]}]];

profileComparisons[eng_,package_] := Module[{rad=eng["radius"],metric,rate,nativeSpeed,orientationOperands,clock,root,coord,map,massRule},
 Do[If[model["name"]=="V4",Continue[]];
  quoteMany[model["quotes"]];
  metric=model["metric"]/.orad->rad;nativeSpeed=metricSpeedSquared[metric];
  rate=domainSimplify[Sqrt[eng["speedSquared"]/nativeSpeed],package["domain"]["CAS"]];
  Do[
   root=Times[clock[[2]],model["profile"]/.orad->rad];
   coord=drift[clockPullback[metric,clock[[2]]]];
   orientationOperands={
    <|"label"->"quoted_nonnegative_profile","equations"->{model["profile_equation"]},"clock_conversion_equations"->model["quotes"],"value"->root|>,
    <|"label"->"coordinate_drift","equations"->model["quotes"],"value"->coord|>};
   If[StringStartsQ[model["name"],"HL"],
    AppendTo[orientationOperands,<|"label"->"reversed_radial_identification","equations"->model["quotes"],"value"->Minus[coord]|>]];
   map=<|"oracle"->model["name"],"radius"->Rule[orad,rad],"clock"->clock,
    "orientation_identification"->"UNMATCHED: equation-root, coordinate drift, and reversed common-radial identification are separate operands; no orientation is selected by a residual.",
    "mass"->"UNMATCHED: independent oracle parameter"|>;
   Do[
    compareFlow[eng,package,orientation["value"],KeyDrop[orientation,"value"],map];
    If[clock[[1]]=="conditional_speed_match" && !MissingQ[model["mass"]],
     massRule=Rule[model["mass"]/.orad->rad,Divide[eng["mass"],Power[clock[[2]],2]]];
     (* Simultaneous one-pass substitution is printed, including any remaining
        parameters when a FORM control makes the proposed map implicit. *)
     compareFlow[eng,package,orientation["value"]/.massRule,KeyDrop[orientation,"value"],Join[map,<|"conditional_mass_rule"->massRule|>]]],
    {orientation,orientationOperands}],{clock,{{"quoted_units",1},{"open_clock",clockRate},{"conditional_speed_match",rate}}}],{model,models}]];

(* Evaluate the entire quoted equation (5) on the common steady radial fields.
   Div is evaluated by Wolfram in the 3D spherical chart; rho remains a function.
   The FORM control changes the operator before this specialization. *)
continuityOnFields[eng_] := Module[{eq,rad=eng["radius"],fieldMap},
 fieldMap={densityField[ot,ox,oy,oz]->eng["density"],velocityX[ot,ox,oy,oz]->eng["velocity"],velocityY[ot,ox,oy,oz]->0,velocityZ[ot,ox,oy,oz]->0};
 emit["CONTINUITY_FIELD_IDENTIFICATION",<|"source"->"Visser (5); conditional common steady, spherically symmetric three-dimensional background", "map"->fieldMap,"chart"->{rad,otheta,odphi}|>];
 eq=ReleaseHold[visContinuity]/.fieldMap;
 eq=eq/.Inactive[Div][v_,{ox,oy,oz}]:>Inactive[Div][v,{rad,otheta,odphi},"Spherical"];
 emit["CONTINUITY_ORACLE_OPERATOR_BEFORE_EVALUATION",eq];
 Activate[eq,D|Div]];

exchangeComparisons[eng_,package_] := Module[{rad=eng["radius"],rho=eng["density"],eq,sign,op,point,v,dv,ex,required,sourceResidual,operatorResidual,domain,zeroeq,sol,sourceDeriv},
 quoteMany[{"V5","V7background","V54"}];
 emit["SOURCE_BALANCE_SIGN_INPUT",eng["balance"]];
 (* These source coefficients fix only the exchange sign convention. *)
 sign=If[eng["name"]=="PY",Cancel[Coefficient[eng["balance"][[2]],D[rho,rad]]/eng["velocity"]],Coefficient[eng["balance"][[2]],jn[r]]];
 emit["SOURCE_EXCHANGE_SIGN",sign];
 eq=continuityOnFields[eng];emit["CONTINUITY_ORACLE_EVALUATED_EQUATION",eq];
 op=Subtract@@(List@@eq);
 Do[domain=cc["Assumptions"];point=cc["EvaluationPoint"];
  progress[{package["tag"],domain["id"],"continuity and exchange"}];
  emit["EXCHANGE_S9B_CASE",cc];
  ex=evaluateCase[cc["ImpliedExchange"][[2]],cc,eng];
  v=evaluateCase[cc["Velocity"],cc,eng];dv=evaluateCase[cc["VelocityDerivative"],cc,eng];
  required=Times[sign,op/.{D[eng["velocity"],rad]->dv,eng["velocity"]->v}]/.rad->point;
  emit["EXCHANGE_IDENTIFICATION_MAP",{"background_continuity_Visser","section7_EOS_and_density","section8_EOS_and_density","no_off_sheet_exchange","Native signed velocity and derivative; live native density; no oracle density substitution."}];
  emit["EXCHANGE_ORACLE_OPERATOR_REQUIRED_SOURCE",required];emit["EXCHANGE_S9B_EVALUATED_OPERAND",ex];
  operatorResidual=report["EXCHANGE_BALANCE_OPERATOR_RESIDUAL",ex,required,domain];
  (* Keep the original source-vs-zero-source comparison and its density set.
     Also print the independently evaluated operator's source requirement. *)
  sourceResidual=report["EXCHANGE_RESIDUAL",ex,Times[sign,Last[eq]],domain];
  report["EXCHANGE_ORACLE_OPERATOR_SOURCE_RESIDUAL",required,Times[sign,Last[eq]],domain];
  zeroeq=Equal[sourceResidual,0];
  emit["DENSITY_ZERO_SET",<|"source_domain"->domain,"density"->(rho/.rad->point),"equation"->zeroeq|>];
  If[!FreeQ[sourceResidual,Derivative] && !FreeQ[sourceResidual,rho],
   sol=Block[{$Assumptions=domain["CAS"]},DSolve[zeroeq,rho,rad]];emit["DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT",sol]];
  emit["HL_CONTINUITY_COMPARISON",<|"oracle"->"UNMATCHED: Hamilton-Lisle states no material continuity/exchange operand.","S9b"->ex,"residual"->"UNMATCHED"|>];
  check[{domain["id"],"exchange_evaluated"},FreeQ[{ex,required},Inactive[Integrate]|Inactive[Limit]]],{cc,package["cases"]}]];

bindingPackages[eng_] := Module[{rules=eng["openBindings"],objects,choices,normalized},
 objects=If[eng["name"]=="PY",{
  "PY_LOCAL_S9B_SPEED_IDENTIFICATION / equality LHS after flow restriction (conditional chi binding)",
  "Hamilton-Lisle (1) / coordinate spatial rr coefficient from line-element Hessian (conditional oracle binding; no S9b metric operand)",
  "PY_LOCAL_S9B_GENERAL_PROFILES / velocity field (conditional U binding)"},{
  "WL_LOCAL_S9B_SUPPLIED_INPUTS / Metric[[1,1]] after WL_LOCAL_S9B_RESTRICTION_FLOW_ONLY",
  "WL_LOCAL_S9B_PROFILE_ANSATZ / V"}];
 choices=MapThread[{{}, {<|"rule"->#1,"object"->#2|>}}&,{rules,objects}];
 If[eng["name"]=="PY",
  normalized=Divide[lineMetric[hlLine][[2,2]],lineMetric[hlLine][[2,2]]];
  AppendTo[choices[[2]],{<|"rule"->Rule[py`s9b$A,normalized],"object"->"Unit radial coefficient imposed by normalization of the oracle spatial rr coefficient; construction-fixed, not an S9b object"|>}]];
 Flatten[#,1]&/@Tuples[choices]];

opticalComparisons[eng_,package_] := Module[{rad=eng["radius"],dom,sourceDisp,disp,sourceCone,metric,rawCone,factor,nativeSpeed,rate,transformed,coneOp,binding,rules,cmp,driftRules},
 dom=package["domain"];
 emit["OPTICAL_INPUT_STATUS","SUPPLIED leading-eikonal input; not an S9b result"];
 emit["S9B_DISPERSION_AS_PRINTED",eng["dispersion"]];emit["FLOW_ONLY_RESTRICTION",eng["restriction"]];
 sourceDisp=eng["dispersion"]/.eng["restriction"];
 emit["S9B_DISPERSION_AFTER_PROFILE_RESTRICTION",sourceDisp];emit["FIXED_DISPERSION_BINDINGS",eng["fixedBindings"]];emit["WAVE_COVECTOR_MAP",eng["waveMap"]];
 Do[quoteMany[model["quotes"]];
  metric=model["metric"]/.orad->rad;factor=model["factor"]/.orad->rad;
  emit["ORACLE_METRIC_OPERAND",metric];emit["ORACLE_OVERALL_COVARIANT_FACTOR",factor];
  rawCone=cone[metric];emit["ORACLE_INVERSE_METRIC_CONE",rawCone];
  emit["ORACLE_CONFORMAL_REMOVAL_OPERATION",<|"covariant_factor"->factor,"inverse_factor"->Power[factor,-1],"frequency_coefficient_divided_out"->Coefficient[Expand[rawCone],frequencyQ,2]|>];
  nativeSpeed=metricSpeedSquared[metric];rate=domainSimplify[Sqrt[eng["speedSquared"]/nativeSpeed],dom["CAS"]];
  Do[transformed=clockPullback[metric,clock[[2]]];coneOp=normalize[cone[transformed]];
   If[eng["name"]=="WL",coneOp=coneOp/.otherQ->0];
   Do[rules=Join[eng["fixedBindings"],Lookup[binding,"rule",{}]]/.orad->rad;
    disp=Fold[ReplaceAll,sourceDisp,rules]/.eng["waveMap"];sourceCone=normalize[disp];
    emit["OPTICAL_IDENTIFICATION_MAP",<|"oracle"->model["name"],"clock"->clock,"binding_operands"->binding,
     "radius_identification"->Rule[orad,rad],"applied_binding_rules"->rules,"fixed_bindings"->eng["fixedBindings"],"unbound_symbols"->Complement[First/@eng["openBindings"],First/@Lookup[binding,"rule",{}]],
     "frame"->"UNMATCHED: conditional common coordinates only","mass"->"UNMATCHED: independent oracle parameters retained"|>];
    cmp=report["OPTICAL_CHARACTERISTIC_RESIDUAL",sourceCone,coneOp,dom];
    If[model["name"]=="V4",
     driftRules=Join[First[Solve[drift[transformed]==eng["velocity"],vxO]],{vyO->0,vzO->0}];emit["CONDITIONAL_RADIAL_VELOCITY_IDENTIFICATION",driftRules];
     report["OPTICAL_CHARACTERISTIC_RESIDUAL",sourceCone,coneOp/.driftRules,dom]];
    check[{eng["name"],model["name"],clock[[1]],"cone_polynomial"},PolynomialQ[cmp,momenta]],
    {binding,bindingPackages[eng]}],{clock,{{"open_clock",clockRate},{"conditional_speed_match",rate}}}];
  emit["LINE_ELEMENT_COMPARISON",<|"oracle"->metric,"S9b"->"UNMATCHED: no line element is printed by S9b.","residual"->"UNMATCHED"|>],{model,models}]];

Do[
 Do[package=prepare[eng,route];AssociateTo[packages,eng["name"]<>route->package];
  profileComparisons[eng,package];exchangeComparisons[eng,package],{route,{"DEFLECTION","RADAR"}}];
 opticalComparisons[eng,packages[eng["name"]<>"DEFLECTION"]],{eng,engineData}];
emit["UNAPPLIED_ORACLES",<|"Unruh1981"->"Full text unverified in the authorized scope; remains open.","Cahill2003"->"Not opened in scoping section 5; survey flags its status; remains open.","Sbitnev"->"Named in survey but not opened in scoping section 5; remains open.","second_route_history"->"Out of scope; remains open with the register."|>];
emit["COMPUTED_STRUCTURAL_GUARDS",checks];
If[!And@@checks[[All,2]],Exit[2]];
emit["INSTRUMENT_COMPLETION",Length[checks]];
Exit[0];
