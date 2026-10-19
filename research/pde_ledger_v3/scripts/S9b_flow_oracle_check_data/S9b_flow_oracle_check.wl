(* Native comparison instrument. The only literal physical combinations are
   in the ORACLE TRANSCRIPTION block. Other expressions are source reads,
   algebra on those reads, or explicitly printed identification substitutions.
   Run via the Python driver in the pooled guard; never directly. *)
$HistoryLength = 0;
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
residual[a_, b_, assumptions_:True] := Simplify[a-b, assumptions];
activateOperand[x_, ass_] := Block[{$Assumptions=ass}, Simplify[Activate[x, Integrate], ass]];
quoteMany[keys_] := Scan[quote, keys];
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

(* Hessian extraction, orthonormal spatial basis, inverse-metric eikonal
   contraction and frequency-coefficient normalization are computations. *)
differentials={odt,odr,odtheta,odphi};
lineMetric[line_] := Table[D[line,differentials[[i]],differentials[[j]]]/2,{i,4},{j,4}];
localMetric[metric_,factor_] := Module[{spatial,basis},
 spatial=Simplify[metric[[2;;4,2;;4]]/factor];
 basis=DiagonalMatrix[Prepend[1/Sqrt[Diagonal[spatial]],1]];
 Simplify[Transpose[basis].metric.basis,orad>0 && 0<otheta<Pi]];
momenta={frequencyQ,radialQ,tangentQ,otherQ};
covector=DiagonalMatrix[{-1,1,1,1}].momenta;
cone[metric_] := Simplify[covector.Inverse[metric].covector];
normalize[polynomial_] := Cancel[polynomial/Coefficient[Expand[polynomial],frequencyQ,2]];
drift[metric_] := Simplify[-metric[[1,2]]/metric[[2,2]]];
clockPullback[metric_,rate_] := Module[{jac=IdentityMatrix[Length[metric]]},jac[[1,1]]=rate;Simplify[Transpose[jac].metric.jac]];
metricSpeedSquared[metric_] := Simplify[-Subtract[metric[[1,1]],metric[[1,2;;4]].Inverse[metric[[2;;4,2;;4]]].metric[[2;;4,1]]]/metric[[2,2]]];

baseHL=localMetric[lineMetric[hlLine],1];
base53=localMetric[lineMetric[visLine53],visFactor53];
base55=localMetric[lineMetric[visLine55],1];
models={
 <|"name"->"HL2","quotes"->{"HL1","HL2"},"metric"->(baseHL/.Apply[Rule,hlProfile]),"factor"->1,"mass"->gmH|>,
 <|"name"->"HL4","quotes"->{"HL1","HL4"},"metric"->(baseHL/.Apply[Rule,hlMassProfile]),"factor"->1,"mass"->massH[orad]|>,
 <|"name"->"V7_plus","quotes"->{"V7background","V53","V86","V87"},"metric"->(base53/.signO->1),"factor"->visFactor53,"mass"->gmV|>,
 <|"name"->"V7_minus","quotes"->{"V7background","V53","V86","V87"},"metric"->(base53/.signO->-1),"factor"->visFactor53,"mass"->gmV|>,
 <|"name"->"V8_plus","quotes"->{"V54","V55","V86","V87"},"metric"->(norm55 base55/.signO->1),"factor"->norm55,"mass"->Missing["NoMassParameter"]|>,
 <|"name"->"V8_minus","quotes"->{"V54","V55","V86","V87"},"metric"->(norm55 base55/.signO->-1),"factor"->norm55,"mass"->Missing["NoMassParameter"]|>};
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
 "openBindings"->{py`s9b$chi->(pySpeed[[1,1]]/.pyRestrict),py`s9b$A->baseHL[[2,2]],py`s9b$U->pyProfiles[[2]]},"fixedBindings"->{}|>,
 <|"name"->"WL","radius"->r,"c0"->c0,"mass"->GM,
 "velocity"->wlProfiles["V"],"density"->wlProfiles["RhoBr"],
 "speedSquared"->(wlSupplied["Speed"][[2]]/.Reverse/@wlSpeedSquareRules),"restriction"->wlRestriction,
 "dispersion"->Subtract@@(List@@wlSupplied["Dispersion"]),
 "waveMap"->Join[Thread[{omega,kr}->{frequencyQ,radialQ}],{kphi->Times[tangentQ,Sqrt[wlSupplied["Metric"][[2,2]]]]}],
 "openBindings"->{aMetric->(wlSupplied["Metric"][[1,1]]/.wlRestriction),vLocal->wlProfiles["V"]},
 "fixedBindings"->Join[{Apply[Rule,wlSupplied["Speed"]]},Reverse/@wlSpeedSquareRules]|>};
(* All source equalities above are converted directly, preserving their sides. *)

profileComparisons[eng_,route_] := Module[{tag,obj,meta,general,squared,cases,definition,vel,ex,ass,rad,roots,allSigned,signed,dom,oraw,ometric,osq,rr,rate,nativeSpeed,conditionalMetric,mappedMass,calibrated,zeroeq,solutions,rho,rawex,caseIndex,domainAss,domainRules,casePoint,oracleAtPoint},
 progress[{eng["name"],route,"flow and exchange"}];
 rad=eng["radius"]; rho=eng["density"]; ass=rad>0;
 tag=eng["name"]<>"_S9B_R_FLOW_ONLY_B_IMPLIED_JN_"<>route;
 If[eng["name"]=="PY",
  obj=data["py"][tag]; emit["S9B_RAW_SOURCE",obj];
  general=node[obj["general"]][[2]]; squared=node[obj["squared"]];
  definition=node[obj["definition"]];
  roots=eng["velocity"]/.Solve[definition/.Apply[Rule,node[obj["general"]]],eng["velocity"]];
  cases=Map[Function[c,<|"Velocity"->node[c["velocity"]][[2]],"ImpliedExchange"->node[c["exchange"]],"EvaluationPoint"->node[c["velocity"]][[1,1]],"AlgebraicDomain"->(And@@(node/@c["algebraic_domain"])),"Domain"->c["domain_raw"],"FullCase"->c["raw"]|>],obj["cases"]],
  obj=wget[tag]; emit["S9B_RAW_SOURCE",data["wl"][tag]];
  meta=First[Cases[obj,a_Association /; KeyExistsQ[a,"GeneralSquaredVelocity"],Infinity]];
  (* Some transcript versions nest VelocityReduction in a case. Native key
     discovery preserves the object instead of retyping it. *)
  If[!AssociationQ[meta],meta=First[Cases[obj,a_Association /; KeyExistsQ[a,"GeneralSquaredVelocity"],Infinity]]];
  emit["S9B_VELOCITY_METADATA",meta];
  general=activateOperand[meta["GeneralSquaredVelocity"],ass];
  squared=activateOperand[meta["SquaredVelocity"],ass];
  roots=unknownRoot/.Solve[Power[unknownRoot,2]==general,unknownRoot];
  cases=obj["Cases"];
 ];
 emit["S9B_GENERAL_SQUARED",general];emit["S9B_RETAINED_SQUARED",squared];
 emit["S9B_GENERAL_SIGNED_ROOTS",roots];
 emit["GENERAL_COMPONENT_DOMAIN",<|"positive_components"->Greater[general,0],"zero_locus"->Equal[general,0],"real_domain"->GreaterEqual[general,0],"qualification"->"Each root applies on a connected positive component; source counting excludes the nonzero integration constant. No extension through a zero is asserted without the printed differentiable joining conditions."|>];
 allSigned=Join[MapIndexed[<|"name"->{"general",First[#2]},"value"->#1,"point"->rad|>&,roots],
 MapIndexed[<|"name"->{"retained",First[#2]},"value"->activateOperand[#1["Velocity"],ass],"point"->Lookup[#1,"EvaluationPoint",rad],"domain"->#1["Domain"]|>&,cases]];
 Do[
  If[model["name"]=="V4",Continue[]];
  quoteMany[model["quotes"]];
  ometric=model["metric"]/.orad->rad;
  oraw=drift[clockPullback[ometric,clockRate]];
  nativeSpeed=metricSpeedSquared[ometric];
  rate=Simplify[Sqrt[eng["speedSquared"]/nativeSpeed],ass && eng["c0"]>0 && soundO>0];
  conditionalMetric=clockPullback[ometric,rate];
  emit["FLOW_IDENTIFICATION_MAP",<|"clock"->"UNMATCHED: clockRate>0 parameterizes the common-radius local time conversion; speed-matched specialization is conditional.","radius"->Rule[orad,rad],"speed_matching_rate"->rate,"mass"->"UNMATCHED: independent oracle parameter; conditional orbital calibration printed separately when a mass operand exists.","orientation"->model["name"]|>];
  Do[
   osq=If[mode=="raw_clock",oraw,drift[conditionalMetric]];
   emit["FLOW_ORACLE_SIGNED_OPERAND",{model["name"],mode,osq}];
   emit["FLOW_ORACLE_SQUARED_OPERAND",Power[osq,2]];
   Do[emit["FLOW_S9B_SQUARED_OPERAND",branch];emit["FLOW_SQUARED_RESIDUAL",residual[branch[[2]],Power[osq,2],ass]],{branch,{{"general",general},{"retained",squared}}}];
   Do[emit["FLOW_S9B_SIGNED_OPERAND",signed];emit["FLOW_ORACLE_AT_COMPONENT_POINT",osq/.rad->signed["point"]];emit["FLOW_SIGNED_RESIDUAL",residual[signed["value"],osq/.rad->signed["point"],ass]],{signed,allSigned}];
   If[mode=="speed_matched" && !MissingQ[model["mass"]],
    (* Explicit conditional dimensional calibration: a geometric mass is
       converted with the square of the clock's length/time conversion. *)
    mappedMass=Rule[model["mass"]/.orad->rad,Divide[eng["mass"],Power[rate,2]]];
    calibrated=Simplify[osq/.mappedMass,ass];
    emit["CONDITIONAL_ORBITAL_CALIBRATION_MAP",mappedMass];
    emit["FLOW_ORACLE_SIGNED_OPERAND",calibrated];
    emit["FLOW_ORACLE_SQUARED_OPERAND",Power[calibrated,2]];
    Do[emit["FLOW_S9B_SQUARED_OPERAND",branch];emit["FLOW_SQUARED_RESIDUAL",residual[branch[[2]],Power[calibrated,2],ass]],{branch,{{"general",general},{"retained",squared}}}];
    Do[emit["FLOW_S9B_SIGNED_OPERAND",signed];emit["FLOW_ORACLE_AT_COMPONENT_POINT",calibrated/.rad->signed["point"]];emit["FLOW_SIGNED_RESIDUAL",residual[signed["value"],calibrated/.rad->signed["point"],ass]],{signed,allSigned}];
   ],{mode,{"raw_clock","speed_matched"}}],{model,models}];
 (* The zero on the oracle side is read from equation (5), not supplied as a
    new S9b law. The full functional zero set is the attached domain conjoined
    with the computed equation; DSolve is additional on open radial components. *)
 quoteMany[{"V5","V7background","V54"}];
 oracleExchange=Last[First[visContinuity]];
 caseIndex=0;
 Do[caseIndex++;
  progress[{eng["name"],route,"exchange",caseIndex}];
  emit["EXCHANGE_ORACLE_OPERAND",oracleExchange];
  emit["EXCHANGE_S9B_CASE",cc];
  rawex=cc["ImpliedExchange"];
  casePoint=Lookup[cc,"EvaluationPoint",rad];
  domainAss=Lookup[cc,"AlgebraicDomain",True];
  emit["EXCHANGE_EVALUATION_POINT",casePoint];
  emit["EXCHANGE_ALGEBRAIC_DOMAIN_FROM_SOURCE",domainAss];
  (* The source zero-set constraint may determine a parameter before a
     boundary limit is evaluated. Solve only that printed algebraic equality;
     retain every other domain condition unchanged. *)
  domainRules=If[!FreeQ[rawex,Inactive[Limit]],Solve[domainAss,eng["mass"]],{{}}];
  emit["EXCHANGE_DOMAIN_PARAMETER_RULES",domainRules];
  ex=activateOperand[rawex[[2]],ass];
  If[!FreeQ[ex,Inactive[Limit]],
   ex=Block[{$Assumptions=ass && domainAss},Simplify[Activate[ex/.First[domainRules],Limit],ass && domainAss]]];
  emit["EXCHANGE_S9B_EVALUATED_OPERAND",ex];
  emit["EXCHANGE_IDENTIFICATION_MAP",{"background_continuity_Visser","section7_EOS_and_density","section8_EOS_and_density","no_off_sheet_exchange","S9b density remains live; only zero source is compared"}];
  rr=residual[ex,oracleExchange,ass];emit["EXCHANGE_RESIDUAL",rr];
  zeroeq=Equal[rr,0];
  emit["DENSITY_ZERO_SET",<|"source_domain"->cc["Domain"],"density"->(rho/.rad->casePoint),"equation"->zeroeq,"definition"->"All density functions satisfying this equation AND every attached source-domain and boundary condition. No oracle density substitution."|>];
  If[!FreeQ[rr,Derivative] && !FreeQ[rr,rho],
   solutions=DSolve[zeroeq,rho,rad];emit["DENSITY_ODE_SOLUTIONS_ON_OPEN_COMPONENT",solutions],
   emit["DENSITY_ZERO_SET_REPRESENTATION",zeroeq]];
  emit["HL_CONTINUITY_COMPARISON",<|"oracle"->"UNMATCHED: Hamilton-Lisle states no material continuity/exchange operand.","S9b"->ex,"residual"->"UNMATCHED"|>];
  check[{tag,"exchange_evaluated",caseIndex},FreeQ[ex,Inactive[Integrate]|Inactive[Limit]]],{cc,cases}];
 check[{tag,"squared_evaluated"},FreeQ[{general,squared},Inactive[Integrate]]];
];

opticalComparisons[eng_] := Module[{rad,ass,sourceDisp,disp,subsets,binding,sourceCone,metric,rawCone,normCone,factor,rawFactor,speed,rate,variants,variant,coneOp,transformed,cmp,deg},
 progress[{eng["name"],"optical characteristics"}];
 rad=eng["radius"];ass=rad>0 && eng["c0"]>0 && soundO>0 && clockRate>0 && norm53!=0 && norm55!=0 && densityO!=0;
 emit["OPTICAL_INPUT_STATUS","SUPPLIED leading-eikonal input; not an S9b result"];
 emit["S9B_DISPERSION_AS_PRINTED",eng["dispersion"]];
 emit["FLOW_ONLY_RESTRICTION",eng["restriction"]];
 sourceDisp=eng["dispersion"]/.eng["restriction"];
 emit["S9B_DISPERSION_AFTER_PROFILE_RESTRICTION",sourceDisp];
 emit["FIXED_DISPERSION_BINDINGS",eng["fixedBindings"]];
 emit["WAVE_COVECTOR_MAP",<|"map"->eng["waveMap"],"convention"->"phase covector (-omega,k); spherical angular momentum converted with the printed spatial metric. WL is compared on its printed two-momentum slice."|>];
 subsets=Subsets[eng["openBindings"]];
 Do[
  quoteMany[model["quotes"]];
  metric=model["metric"]/.orad->rad;
  factor=model["factor"]/.orad->rad;
  emit["ORACLE_OVERALL_COVARIANT_FACTOR",factor];
  emit["ORACLE_METRIC_OPERAND",metric];
  rawCone=cone[metric];rawFactor=Coefficient[Expand[rawCone],frequencyQ,2];
  emit["ORACLE_INVERSE_METRIC_CONE",rawCone];
  emit["ORACLE_CONE_NORMALIZATION_DIVISOR",rawFactor];
  normCone=normalize[rawCone];
  emit["ORACLE_NORMALIZED_CONE",normCone];
  emit["ORACLE_CONFORMAL_REMOVAL_OPERATION",<|"covariant_factor"->factor,"inverse_factor"->Power[factor,-1],"frequency_coefficient_divided_out"->rawFactor|>];
  speed=metricSpeedSquared[metric];rate=Simplify[Sqrt[eng["speedSquared"]/speed],ass];
  variants={{"open_clock",clockRate},{"conditional_speed_match",rate}};
  Do[
   transformed=clockPullback[metric,variant[[2]]];
   coneOp=normalize[cone[transformed]];
   If[eng["name"]=="WL",coneOp=coneOp/.otherQ->0];
   emit["OPTICAL_CLOCK_IDENTIFICATION",variant];
   emit["ORACLE_CONE_OPERAND",coneOp];
   Do[
    binding=Join[eng["fixedBindings"],conditional];
    disp=Fold[ReplaceAll,sourceDisp,binding];
    disp=disp/.eng["waveMap"];
    sourceCone=normalize[disp];
    emit["OPTICAL_IDENTIFICATION_MAP",<|"conditional_bindings"->conditional,"unbound_entries"->Complement[eng["openBindings"],conditional],"fixed_bindings"->eng["fixedBindings"],"oracle"->model["name"],"clock"->variant,"frame"->"UNMATCHED: conditional common coordinates only","mass"->"UNMATCHED: independent oracle parameters retained"|>];
    emit["OPTICAL_S9B_OPERAND",sourceCone];
    emit["OPTICAL_ORACLE_OPERAND",coneOp];
    cmp=residual[sourceCone,coneOp,ass];emit["OPTICAL_CHARACTERISTIC_RESIDUAL",cmp];
    (* A second explicitly conditional map identifies the drift components
       directly, without calibrating them from a profile residual. *)
    driftRules=If[model["name"]=="V4",Join[First[Solve[drift[transformed]==eng["velocity"],vxO]],{vyO->0,vzO->0}],{}];
    If[Length[driftRules]>0,
     emit["CONDITIONAL_RADIAL_VELOCITY_IDENTIFICATION",driftRules];
     emit["OPTICAL_S9B_OPERAND",sourceCone];
     emit["OPTICAL_ORACLE_OPERAND",coneOp/.driftRules];
     emit["OPTICAL_CHARACTERISTIC_RESIDUAL",residual[sourceCone,coneOp/.driftRules,ass]]];
    check[{eng["name"],model["name"],variant[[1]],"cone_polynomial"},PolynomialQ[cmp,momenta]],
    {conditional,subsets}],{variant,variants}];
  emit["LINE_ELEMENT_COMPARISON",<|"oracle"->metric,"S9b"->"UNMATCHED: no line element is printed by S9b.","residual"->"UNMATCHED"|>],{model,models}];
];

Do[Do[profileComparisons[eng,route],{route,{"DEFLECTION","RADAR"}}];opticalComparisons[eng],{eng,engineData}];
emit["UNAPPLIED_ORACLES",<|"Unruh1981"->"Full text unverified in the authorized scope; remains open.","Cahill2003"->"Not opened in scoping section 5; survey flags its status; remains open.","Sbitnev"->"Named in survey but not opened in scoping section 5; remains open.","second_route_history"->"Out of scope; remains open with the register."|>];
emit["COMPUTED_STRUCTURAL_GUARDS",checks];
If[!And@@checks[[All,2]],Exit[2]];
emit["INSTRUMENT_COMPLETION",Length[checks]];
Exit[0];
