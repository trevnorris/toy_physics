SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
(* Independent S15 construction from the supplied equations. No external input.
   Boolean-lattice finite differences partition polynomial terms by label support.
   Integrals are retained as exact inactive integrals of general live functions.
   No numerical or symbolic value is an acceptance criterion. *)
SetOptions[FullSimplify, TimeConstraint -> Infinity];
SetOptions[Simplify, TimeConstraint -> Infinity];
$HistoryLength = 0;
$MaxExtraPrecision = 0;
(* Execution selection contains only object names and file paths. Each kernel
   reconstructs the supplied fields and the selected closed object in this file;
   no expression, condition or result is imported from another kernel. *)
s15RunMode = Environment["S15_WL_MODE"];
If[!MemberQ[{"base", "item4"}, s15RunMode], s15RunMode = "single"];
s15SelectedObject = Environment["S15_WL_OBJECT"];
s15OutputPath = Environment["S15_WL_TRANSCRIPT"];
output = If[StringQ[s15OutputPath] && StringLength[s15OutputPath]>0,
 s15OutputPath, FileNameJoin[{DirectoryName[$InputFileName], "out", "S15_separation_departure.out"}]];
stream = If[s15RunMode === "item4", OpenAppend[output], OpenWrite[output]];
emit[name_, value_] := If[s15RunMode =!= "item4" ||
 StringStartsQ[name, "reading/"<>s15SelectedObject<>"/item4/"] ||
 MemberQ[{"WORKING_ON", "ENGINE_MESSAGES", "EXACT_OBJECT_GAP", "END_OF_ITEM4_PASS"}, name],
 WriteString[stream, name, ": ", ToString[value/.eulerAppell->Inactive[AppellF1], InputForm, PageWidth -> Infinity], "\n"]];
currentObject="construction";
stepProgress[step_] := (Print["WORKING_ON: ",currentObject," / ",step];
 emit["WORKING_ON",currentObject<>" / "<>step];Close[stream];stream=OpenAppend[output]);
progress[name_] := (currentObject=name;stepProgress["object"]);
SetAttributes[casStep,HoldRest];
casStep[step_,expression_] := Module[{result},stepProgress[step];result=expression;
 Print["FINISHED_STEP: ",currentObject," / ",step];result];
SetAttributes[trackedIntegrate,HoldAll]; trackedIntegrate[args___] := casStep["Integrate "<>ToString[HoldForm[{args}],InputForm],Integrate[args]];
SetAttributes[trackedSolve,HoldAll]; trackedSolve[args___] := casStep["Solve "<>ToString[HoldForm[{args}],InputForm],Solve[args]];
SetAttributes[trackedSolveAlways,HoldAll]; trackedSolveAlways[args___] := casStep["SolveAlways "<>ToString[HoldForm[{args}],InputForm],SolveAlways[args]];
SetAttributes[trackedResolve,HoldAll]; trackedResolve[args___] := casStep["Resolve "<>ToString[HoldForm[{args}],InputForm],Resolve[args]];
SetAttributes[trackedReduce,HoldAll]; trackedReduce[args___] := casStep["Reduce "<>ToString[HoldForm[{args}],InputForm],Reduce[args]];
SetAttributes[trackedPolynomialReduce,HoldAll]; trackedPolynomialReduce[args___] := casStep["PolynomialReduce "<>ToString[HoldForm[{args}],InputForm],PolynomialReduce[args]];
SetAttributes[trackedFullSimplify,HoldAll]; trackedFullSimplify[args___] := casStep["FullSimplify "<>ToString[HoldForm[{args}],InputForm],FullSimplify[args]];
qObjects=<||>;
qName[name_] := StringMatchQ[name, ("inner"|"outer")~~"/"~~("carried"|"pressure"|"viscous"|"other"|"transport_centre"|"transport_radius")~~"/"~~__] ||
 Module[{parts=StringSplit[name,"/"]}, Length[parts]>=3 && First[parts]==="region" &&
  (StringStartsQ[parts[[3]],"momentum"] || StringStartsQ[parts[[3]],"rate_labels_"] || parts[[3]]==="integrated_force")];
solo[name_, value_] := (If[qName[name],AssociateTo[qObjects,name->value]];emit[name, value]; emit[name <> "/comparison", "There is nothing to compare it with in this engine."]);
pair[name_, a_, b_] := (emit[name <> "/operand_a", a]; emit[name <> "/operand_b", b]; emit[name <> "/residual", Expand[a-b]]);

Check[
progress["domain and supplied fields"];
xyz = {x, y, z}; labels = {drain1, drain2, remainder}; marks = {a1, a2, ar};
sets = Rest[Subsets[Range[3]]];
centres = Table[xx[k, j][t], {k, 2}, {j, 3}];
radii = Table[Sqrt[Total[(xyz-centres[[k]])^2]], {k, 2}];
velocities = Table[uu[k][radii[[k]], t] (xyz-centres[[k]])/radii[[k]], {k, 2}];
densities = Table[rr[k][radii[[k]], t], {k, 2}];
AppendTo[velocities, Table[ww[j][x,y,z,t], {j,3}]];
AppendTo[densities, vr[x,y,z,t]];
velocity = Total[MapThread[Times, {marks, velocities}]];
density = marks.densities;
etaField = eta[x,y,z,t]; zetaField = zeta[x,y,z,t];
other = Table[oo[i,j][x,y,z,t], {i,3},{j,3}];
force = Table[ff[i][x,y,z,t], {i,3}];
untag[e_] := e /. Thread[marks -> 1];
grad = Table[D[velocity[[j]], xyz[[i]]], {i,3},{j,3}];
stress = -pb[density] IdentityMatrix[3] + etaField (grad+Transpose[grad] - (2/3) IdentityMatrix[3] Tr[grad]) + zetaField IdentityMatrix[3] Tr[grad];
momentum = density velocity;
current = density Outer[Times,velocity,velocity]-stress+other;
localLeft = D[momentum,t] + Table[Sum[D[current[[i,j]],xyz[[j]]],{j,3}],{i,3}];
sep = Sqrt[Total[(centres[[1]]-centres[[2]])^2]];
domain = sep > rout[t] > rin[t] > 0;
emit["DOMAIN_FOR_EVERY_OBJECT", {domain, "both cores excluded; drain 1 outside the outer sphere; projected slab domain of S14a part 1"}];
emit["PRESSURE_DOMAIN", "Whole pressure: pb defined at total density. Expansion and subtraction: also pb defined and differentiable at rr[2]+vr. Integration requires these pointwise on its surface."];
emit["REGULARITY", "Derivatives, surface/volume integrals and transport theorem are on their classical domains of existence; no steady or other limiting regime is imposed."];
emit["PREMISES", "Supplied P5 and P6 applied to projected fields; steady-origin P6 used as written for moving drains; departures reside in oo; ff remains undivided."];
emit["CHOICES", "Radial profiles about moving centres; labelled split with unrestricted remainders; concentric moving spheres; first order pressure expansion where differentiable."];
emit["BOOKKEEPING", "a1,a2,ar are dimensionless label markers. Their evaluations below project formal monomials, never assign physical fields. All printed physical objects have markers removed."];
emit["LIVE_PARAMETER_NOT_IN_SUPPLIED_BALANCE", ccomp[x,y,z,t]];
solo["radial_distances",radii]; solo["drain_and_remainder_velocities",velocities]; solo["drain_and_remainder_densities",densities];
solo["total_velocity",untag[velocity]]; solo["total_density",untag[density]];
solo["P6_stress",untag[stress]]; solo["P5_momentum_density",untag[momentum]]; solo["momentum_current",untag[current]];
pair["supplied_local_balance",untag[localLeft],force];
emit["supplied_local_balance/equation",Equal[untag[localLeft],force]];

(* Extract current channels by field occurrence, rather than writing their fluxes. *)
channel[term_] := Which[!FreeQ[term,pb], "pressure", !FreeQ[term,eta|zeta], "viscous", !FreeQ[term,oo], "other", True, "carried"];
channelNames = {"carried","pressure","viscous","other"};
terms[e_] := With[{v=Expand[e]}, If[Head[v]===Plus,List@@v,{v}]];
channels = Association@Table[name -> Table[Total[Select[terms[current[[i,j]]],channel[#]===name&]],{i,3},{j,3}],{name,channelNames}];
pair["current_channel_partition",current,Total[Values[channels]]];

(* Group by exact support using the Boolean subset lattice. *)
activate[e_,s_] := e /. Thread[marks -> Table[Boole[MemberQ[s,j]],{j,3}]];
supportPart[e_,s_] := Expand[Total[Table[(-1)^(Length[s]-Length[q]) activate[e,q],{q,Subsets[s]}]]];
partition[e_] := Association@Table[ToString[labels[[s]],InputForm] -> supportPart[e,s],{s,sets}];
pressurePoly = Normal[Series[pb[density],{a1,0,1}]];
(* Series is used only for its first two Taylor coefficients. Its O(a1^2)
   bookkeeping is NOT asserted under merely differentiable pb. *)
pressureTerms = terms[pressurePoly];
pressureGroups = Association@Table[ToString[labels[[s]],InputForm] -> untag[Total[Select[pressureTerms,Function[q,Select[Range[3],!FreeQ[q,marks[[#]]]&]===s]]]],{s,sets}];
solo["pressure_whole",untag[pb[density]]];
solo["pressure_expansion",untag[pressurePoly]];
emit["pressure_expansion/order_and_condition", "Labelled choice: degree one in drain-1 density about drain-2 density plus remainder. pb defined and differentiable there; local remainder o(drain-1 density), not an asserted quadratic bound."];
solo["pressure_whole_minus_expansion",untag[pb[density]-pressurePoly]];
pair["pressure_group_partition",untag[pressurePoly],Total[Values[pressureGroups]]];

(* Spherical charts. Metric determinants supply measures. *)
angles = {theta,phi};
direction = {Sin[theta] Cos[phi],Sin[theta] Sin[phi],Cos[theta]};
chart[b_] := centres[[2]]+b direction;
normal[orientation_] := orientation direction;
surfaceVelocity[b_] := D[chart[b],t];
metricMeasure[map_,coordinates_,assumptions_] := trackedFullSimplify[Sqrt[Det[Table[D[map,q].D[map,r],{q,coordinates},{r,coordinates}]]],Assumptions->assumptions];
assumptions = 0 < theta < Pi && Element[phi,Reals] && rin[t]>0 && rout[t]>0 && s>0;
volumeMap = chart[s];
volumeMeasure = metricMeasure[volumeMap,{s,theta,phi},assumptions];
pull[e_,b_] := e /. Thread[xyz -> chart[b]];
integralOrigins=<||>;rateOrigins=<||>;
sint[e_,b_,measure_] := Module[{obj=Inactive[Integrate][pull[e,b] measure,{theta,0,Pi},{phi,0,2 Pi}]},
 AssociateTo[integralOrigins,obj->{e,b,Cancel[measure/Sin[theta]],False}];obj];
vint[e_] := Module[{obj=Inactive[Integrate][pull[e,s] volumeMeasure,{s,rin[t],rout[t]},{theta,0,Pi},{phi,0,2 Pi}]},
 AssociateTo[integralOrigins,obj->{e,s,Cancel[volumeMeasure/Sin[theta]],True}];obj];
(* Generic Leibniz rule on the radial member of the chart integral. *)
rate[e_] := Module[{integrand=pull[e,s] volumeMeasure,limits={rin[t],rout[t]},sgn={-1,1},obj},
 obj=Inactive[Integrate][D[integrand,t],{s,limits[[1]],limits[[2]]},{theta,0,Pi},{phi,0,2 Pi}] +
 Total[Table[sgn[[j]] D[limits[[j]],t] Inactive[Integrate][integrand /. s->limits[[j]],{theta,0,Pi},{phi,0,2 Pi}],{j,2}]];
 AssociateTo[rateOrigins,obj->vint[e]];obj];

(* Dimensional algebra. Dens is an abstract density unit, with no mass assignment.
   Exponents are {Dens,Length,Time}. Field signatures include their arguments. *)
zeroDim={0,0,0}; lenDim={0,1,0}; timeDim={0,0,1}; densDim={1,0,0};
functionDims=<||>; argumentDims=<||>;
register[h_,d_,args_] := (AssociateTo[functionDims,ToString[h,InputForm]->d];AssociateTo[argumentDims,ToString[h,InputForm]->args]);
Do[register[xx[k,j],lenDim,{timeDim}],{k,2},{j,3}];
Do[register[uu[k],lenDim-timeDim,{lenDim,timeDim}];register[rr[k],densDim,{lenDim,timeDim}],{k,2}];
Do[register[ww[j],lenDim-timeDim,Join[ConstantArray[lenDim,3],{timeDim}]];register[ff[j],densDim+lenDim-2 timeDim,Join[ConstantArray[lenDim,3],{timeDim}]],{j,3}];
Do[register[oo[i,j],densDim+2 lenDim-2 timeDim,Join[ConstantArray[lenDim,3],{timeDim}]],{i,3},{j,3}];
Do[register[h,densDim+2 lenDim-timeDim,Join[ConstantArray[lenDim,3],{timeDim}]],{h,{eta,zeta}}];
register[vr,densDim,Join[ConstantArray[lenDim,3],{timeDim}]];
register[pb,densDim+2 lenDim-2 timeDim,{densDim}];
Do[register[h,lenDim,{timeDim}],{h,{rin,rout}}];
register[ccomp,lenDim-timeDim,Join[ConstantArray[lenDim,3],{timeDim}]];
(* A list of dimensions retains additive alternatives; {} represents literal 0. *)
dim[0] := {};
dim[e_?NumberQ] := {zeroDim};
dim[Pi] := {zeroDim};
dim[e_Symbol] := {Which[MemberQ[Join[xyz,{s}],e],lenDim,e===t,timeDim,MemberQ[Join[angles,marks],e],zeroDim,True,Throw[{"unregistered dimension symbol",e}]]};
dim[e_Plus] := DeleteDuplicates[Flatten[dim /@ (List@@e),1]];
dim[e_Times] := Module[{ds=dim /@ (List@@e)},If[MemberQ[ds,{}],{},DeleteDuplicates[Total /@ Tuples[ds]]]];
dim[Power[e_,n_?NumberQ]] := (# n& /@ dim[e]);
dim[Sin[_]|Cos[_]] := {zeroDim};
dim[Derivative[ns__][h_][args___]] := {functionDims[ToString[h,InputForm]] - Total[MapThread[Times,{{ns},argumentDims[ToString[h,InputForm]]}]]};
dim[Inactive[Integrate][e_,ranges__]] := Module[{base=dim[e],ds=Total[First[dim[First[#]]]& /@ {ranges}]}, (# + ds& /@ base)];
dim[e_] := Module[{key=ToString[Head[e],InputForm]},If[KeyExistsQ[functionDims,key],{functionDims[key]},Throw[{"unregistered dimension expression",e}]]];
emit["dimension_basis",{Dens,Length,Time}]; emit["field_dimension_signatures",functionDims]; emit["field_argument_dimension_signatures",argumentDims];
emit["dimension_convention", "Literal 0 has no preferred dimension (empty dimension list); sums retain all distinct term dimensions. Dens is left abstract, with no mass or normalization assigned."];
dimensionAudit[name_,e_] := Module[{ds=dim[e]},emit[name<>"/dimension_operands",ds];emit[name<>"/dimension_pair_residuals",Subtract@@#& /@ Subsets[ds,{2}]]];

progress["surface crossings"];
surfaceResults={};
Do[
 b=spec[[2]]; orientation=spec[[3]]; prefix=spec[[1]];
 progress[prefix];
 measure=metricMeasure[chart[b],angles,assumptions]; nn=normal[orientation]; sv=surfaceVelocity[b];
 (* Radius/centre velocity channels are derivatives of the allowed chart. *)
 cv=D[centres[[2]],t]; rv=sv-cv;
 solo[prefix<>"/chart",chart[b]];solo[prefix<>"/outward_normal",nn];solo[prefix<>"/area_measure",measure];
 solo[prefix<>"/surface_velocity",sv];solo[prefix<>"/centre_velocity",cv];solo[prefix<>"/radius_velocity",rv];
 pair[prefix<>"/velocity_partition",sv,cv+rv];
 crossing=Association@Table[name -> -(channels[name].nn),{name,channelNames}];
 AssociateTo[crossing,{"transport_centre" -> momentum (cv.nn),"transport_radius" -> momentum (rv.nn)}];
 Do[
 progress[prefix<>"/"<>name];
 Do[
  raw=crossing[name][[i]]; tag=prefix<>"/"<>name<>"/component_"<>ToString[i];
  obj=sint[untag[raw],b,measure];solo[tag,obj];dimensionAudit[tag,obj];
  If[!MemberQ[{"pressure","other"},name],
   groups=partition[raw];
   KeyValueMap[Function[{key,value},solo[tag<>"/labels_"<>key,sint[value,b,measure]];dimensionAudit[tag<>"/labels_"<>key,sint[value,b,measure]]],groups];
   pair[tag<>"/group_partition_integrand",untag[raw],Total[Values[groups]]]
  ];
  If[name==="pressure",
   expanded=raw /. pb[density]->pressurePoly;
   expansionIntegral=sint[untag[expanded],b,measure];
   solo[tag<>"/first_order_choice",expansionIntegral];dimensionAudit[tag<>"/first_order_choice",expansionIntegral];
   emit[tag<>"/first_order_domain", "pb defined and differentiable at rr[2]+vr pointwise; degree one in rr[1]; whole pb defined at total density for subtraction."];
   residualIntegral=sint[untag[raw-expanded],b,measure];solo[tag<>"/whole_minus_expansion",residualIntegral];dimensionAudit[tag<>"/whole_minus_expansion",residualIntegral];
   (* The pressure traction is linear in pb, so apply its computed coefficient. *)
   pressureCoefficient=Coefficient[raw,pb[density]];
   KeyValueMap[Function[{key,value},solo[tag<>"/expansion_labels_"<>key,sint[untag[pressureCoefficient] value,b,measure]];dimensionAudit[tag<>"/expansion_labels_"<>key,sint[untag[pressureCoefficient] value,b,measure]]],pressureGroups]
  ],{i,3}],{name,Keys[crossing]}];
 AppendTo[surfaceResults,Table[Association@Table[name->sint[untag[crossing[name][[i]]],b,measure],{name,Keys[crossing]}],{i,3}]],
 {spec,{{"inner",rin[t],-1},{"outer",rout[t],1}}}];

progress["region momentum and balance"];
solo["region/chart",volumeMap];solo["region/volume_measure",volumeMeasure];
Do[
 tag="region/component_"<>ToString[i];
 mom=vint[untag[momentum[[i]]]]; lhs=rate[untag[momentum[[i]]]]; forcing=vint[force[[i]]];
 solo[tag<>"/momentum",mom];dimensionAudit[tag<>"/momentum",mom];
 solo[tag<>"/momentum_rate",lhs];dimensionAudit[tag<>"/momentum_rate",lhs];
 groups=partition[momentum[[i]]];
 KeyValueMap[Function[{key,value},solo[tag<>"/momentum_labels_"<>key,vint[value]];solo[tag<>"/rate_labels_"<>key,rate[value]]],groups];
 pair[tag<>"/momentum_group_partition_integrand",untag[momentum[[i]]],Total[Values[groups]]];
 solo[tag<>"/integrated_force",forcing];dimensionAudit[tag<>"/integrated_force",forcing];
 rhs=Total[Flatten[Table[Values[surfaceResults[[j,i]]],{j,2}]]]+forcing;
 emit[tag<>"/balance_parts",{surfaceResults[[1,i]],surfaceResults[[2,i]],forcing}];
 pair[tag<>"/balance",lhs,rhs];emit[tag<>"/balance_equation",Equal[lhs,rhs]];
 emit[tag<>"/balance/reduction", "This balance residual has no reduction route in this engine using the supplied local balance."];
 dimensionAudit[tag<>"/balance_rhs",rhs];
 emit[tag<>"/balance_dimension_operand_a",dim[lhs]];emit[tag<>"/balance_dimension_operand_b",dim[rhs]];
 emit[tag<>"/balance_dimension_residuals",Flatten[Table[da-db,{da,dim[lhs]},{db,dim[rhs]}],1]],
 {i,3}];
emit["balance_construction", "Printed objects: the chart momentum integral, its Leibniz time derivative, each surface contribution, the undivided force integral, and the balance operands, equation and residual."];
emit["scope", "No force law, coefficient, normalization, mass, throat-force identification, or relation to J, Q, Theta_Q or N_infinity_3 is computed."];
(* Reading engine. This construction uses differential operators on independent
   field jets, rather than importing a perturbed expression or any other code. *)
progress["reading family"];
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
changingHeads=Join[Table[uu[k],{k,2}],Table[rr[k],{k,2}],Table[ww[j],{j,3}],{vr,eta,zeta},Flatten[Table[oo[i,j],{i,3},{j,3}]],Table[ff[j],{j,3}]];
Do[register[change[h],functionDims[ToString[h,InputForm]],argumentDims[ToString[h,InputForm]]],{h,changingHeads}];
Do[register[centreVelocityChange[j],lenDim-timeDim,{timeDim}],{j,3}];
centreVelocityAtoms=D[centres[[1]],t];
centreVelocityRules=Thread[centreVelocityAtoms->Table[centreVelocityChange[j][t],{j,3}]];
positionSymbols=Table[Unique["position"],{j,3}];
positionForward=Thread[centres[[1]]->positionSymbols];
positionBackward=Thread[positionSymbols->centres[[1]]];
separationVector=centres[[1]]-centres[[2]];
emit["family/held",{centres[[2]],D[centres[[2]],t],rin[t],rout[t],rin'[t],rout'[t],pb}];
emit["family/changing_heads",changingHeads];
emit["family/field_change_notation",Thread[changingHeads->(change/@changingHeads)]];
emit["family/velocity_changes",centreVelocityRules];
emit["family/dimension_signatures_including_changes",functionDims];
emit["family/separation_vector",separationVector];
emit["family/choice","Instantaneous comparison; X1 position and X1 velocity are independent jets. Held quantities do not acquire a response. Each change has the arguments and dimensions of its field."];

(* Differentiate under integrals with fixed family bounds; generic product and
   sum differentiation retains all boundary terms already produced by rate. *)
geometricDerivative[e:Inactive[Integrate][f_,ranges__]] := geometricDerivative[e] = Inactive[Integrate][geometricDerivative[f],ranges];
geometricDerivative[e_Plus] := geometricDerivative[e] = Total[geometricDerivative/@(List@@e)];
geometricDerivative[e_Times] /; !FreeQ[e,Inactive[Integrate]] := geometricDerivative[e] =
 Total[MapIndexed[geometricDerivative[#1] Times@@Delete[List@@e,First[#2]]&,List@@e]];
geometricDerivative[e_] := geometricDerivative[e] = Module[{lift=e/.positionForward},
 (Table[D[lift,positionSymbols[[j]]],{j,3}].separationVector)/.positionBackward];

responseJet[Derivative[orders__][h_][arguments___]] := Derivative[orders][change[h]][arguments];
responseJet[h_[arguments___]] := change[h][arguments];
responseDerivative[e:Inactive[Integrate][f_,ranges__]] := responseDerivative[e] = Inactive[Integrate][responseDerivative[f],ranges];
responseDerivative[e_Plus] := responseDerivative[e] = Total[responseDerivative/@(List@@e)];
responseDerivative[e_Times] /; !FreeQ[e,Inactive[Integrate]] := responseDerivative[e] =
 Total[MapIndexed[responseDerivative[#1] Times@@Delete[List@@e,First[#2]]&,List@@e]];
responseDerivative[e_] := responseDerivative[e] = Module[{atoms,slots,lift,changes},
 atoms=DeleteDuplicates[Join[
  Cases[e,f:h_[args___]/;MemberQ[changingHeads,h]:>f,{0,Infinity}],
  Cases[e,f:Derivative[ns__][h_][args___]/;MemberQ[changingHeads,h]:>f,{0,Infinity}],
  Select[centreVelocityAtoms,!FreeQ[e,#]&]]];
 slots=Table[Unique["jet"],{Length[atoms]}];
 lift=e/.Thread[atoms->slots];
 changes=Map[If[MemberQ[centreVelocityAtoms,#],#/.centreVelocityRules,responseJet[#]]&,atoms];
 Total[MapThread[(D[lift,#1]/.Thread[slots->atoms]) #2&,{slots,changes}]]];

(* Power profiles are definitions; every derivative of them is computed. *)
profileHeads=Join[Table[uu[k],{k,2}],Table[rr[k],{k,2}]];
profileDefinitions=Association[Join[Table[uu[k]->(AA[k][profileTime] profileRadius^aa[k]),{k,2}],Table[rr[k]->(BB[k][profileTime] profileRadius^bb[k]),{k,2}]]];
profileEvaluate[h_,orders_,args_] := (Fold[D[#1,#2]&,profileDefinitions[h],Flatten[MapThread[ConstantArray,{{profileRadius,profileTime},orders}]]]/.Thread[{profileRadius,profileTime}->args]);
profileRules={
 Derivative[ns__][h_][args___]/;MemberQ[profileHeads,h]:>profileEvaluate[h,{ns},{args}],
 h_[args___]/;MemberQ[profileHeads,h]:>profileEvaluate[h,{0,0},{args}]};
separationChart=Thread[centres[[1]]->(centres[[2]]+dist Table[directionSeparation[j],{j,3}])];
unitDirection=Total[Table[directionSeparation[j]^2,{j,3}]]==1;
familyObject[e_] := familyObject[e] = (e/.profileRules)/.separationChart;
emit["profile_family/definitions",profileDefinitions];
emit["profile_family/separation_chart",separationChart];
emit["profile_family/domain",{dist>rout[t],unitDirection,Element[Join[{dist},Table[directionSeparation[j],{j,3}],Table[aa[k],{k,2}],Table[bb[k],{k,2}]],Reals]}];
emit["profile_family/chart_note","Arbitrary unit direction in the original Cartesian basis. No coordinate alignment or velocity restriction is imposed."];
allowedHeads=Join[Flatten[Table[xx[k,j],{k,2},{j,3}]],{rin,rout},Table[AA[k],{k,2}],Table[BB[k],{k,2}]];
hasGeneral[e_] := !FreeQ[e,uu|rr|ww|vr|eta|zeta|oo|ff|pb|change|centreVelocityChange];

(* Item 3. Retain the pre-chart tensor integrands as provenance. The angular
   integration below is derived with a generating polynomial in an auxiliary
   source, independently of the Cartesian grouping and variation operators. *)
normalSymbols={normalX,normalY,normalZ};
sourceSymbols={sourceX,sourceY,sourceZ};
axisSymbols={axisX,axisY,axisZ};
axisAssumptions=Element[Join[axisSymbols,{mu}],Reals] && Total[axisSymbols^2]==1 && -1<mu<1;
emit["item3_kernel/generating_chart_condition",axisAssumptions];
angularGenerator[degree_] := angularGenerator[degree] = Module[{input,value},
 input=Expand[(axial mu+transverse Sqrt[1-mu^2] Cos[azimuth])^degree];
 value=trackedIntegrate[input,{azimuth,0,2 Pi},Assumptions->Element[{axial,transverse,mu},Reals] && -1<mu<1];
 solo["item3_kernel/generator_operand",Inactive[Integrate][input,{azimuth,0,2 Pi}]];
 solo["item3_kernel/generator_integral",value];
 value=Expand[value]/.transverse^n_Integer?EvenQ:>(sourceSymbols.sourceSymbols-(sourceSymbols.axisSymbols)^2)^(n/2);
 value=value/.axial->sourceSymbols.axisSymbols;
 Expand[value]];
angularMoment[powers_] := angularMoment[powers] = Module[{degree=Total[powers],generator,value},
 generator=angularGenerator[degree];
 value=Fold[D[#1,#2]&,generator,Flatten[MapThread[ConstantArray,{sourceSymbols,powers}]]]/Factorial[degree];
 value=Expand[value/.Thread[sourceSymbols->0]];
 solo["item3_kernel/angular_moment",{powers,value}];value];
(* Derive the distance substitution from dot products on the unit sphere. *)
(* The axis in this mathematical chart is the pole. The generating
   polynomial retains its arbitrary direction in the Cartesian output. *)
distanceChart={Sqrt[1-mu^2],0,mu};
squaredDistance=trackedFullSimplify[(radial distanceChart-distance {0,0,1}).(radial distanceChart-distance {0,0,1}),Assumptions->Element[{radial,distance,mu},Reals]];
muRule=First[trackedSolve[squared==squaredDistance,mu]];
distanceJacobian=-D[mu/.muRule,squared];
distanceBounds=Factor[squaredDistance/.mu->#]& /@ {1,-1};
solo["item3_kernel/distance_substitution",{squaredDistance,muRule,distanceJacobian,distanceBounds}];

(* Mellin primitive: Integrate supplies both ordinary and logarithmic branches. *)
powerPrimitive=trackedIntegrate[squared^power,squared,Assumptions->squared>0 && Element[power,Reals],GenerateConditions->True];
solo["item3_kernel/power_primitive_operand",squared^power];
solo["item3_kernel/power_primitive_return",powerPrimitive];
(* Compute the exceptional parameter from the primitive denominator and its
   exact integral there, rather than dropping Integrate's excluded set. *)
primitiveGeneric=trackedIntegrate[squared^power,squared,Assumptions->squared>0,GenerateConditions->False];
primitiveDenominator=Denominator[Together[primitiveGeneric]];
primitiveExceptions=trackedSolve[primitiveDenominator==0,power,Reals];
primitiveBranches=Join[Table[{trackedIntegrate[squared^power/.rule,squared,Assumptions->squared>0],And@@(Equal@@#& /@ (List@@@rule))},{rule,primitiveExceptions}],{{primitiveGeneric,True}}];
solo["item3_kernel/power_primitive_branches",primitiveBranches];
pair["item3_kernel/power_primitive_derivative",D[primitiveGeneric,squared],squared^power];

(* The Appell kernel uses fixed positive Euler parameters (1,2). Obtain its
   coefficient normalization by exact Beta integration. No series truncation
   is used to replace an object; continuation is along positive real bases. *)
eulerMoment=trackedIntegrate[unitVariable^(index1+index2),{unitVariable,0,1},Assumptions->Element[{index1,index2},Integers] && index1>=0 && index2>=0];
appellMoment=Pochhammer[1,index1+index2]/Pochhammer[2,index1+index2];
pair["item3_kernel/Euler_moment",eulerMoment,appellMoment];
solo["item3_kernel/Euler_moment_residual",trackedFullSimplify[eulerMoment-appellMoment,Assumptions->Element[{index1,index2},Integers] && index1>=0 && index2>=0]];
emit["item3_kernel/Appell_domain","Euler continuation between positive endpoints, distance > upper > lower > 0. AppellF1 has fixed first/lower parameters 1 and 2; all real radial and affine powers are included."];
(* Retain the Euler analytic function during differentiation. Its argument
   derivatives are computed from exact coefficient identities; no generic
   parameter reduction is allowed to introduce an untracked pole. *)
eulerCoefficient[av_,bv_,dv_,cv_,iv_,jv_] := Pochhammer[av,iv+jv] Pochhammer[bv,iv] Pochhammer[dv,jv]/(Pochhammer[cv,iv+jv] Factorial[iv] Factorial[jv]);
coefficientDomain=Element[{ei,ej},Integers] && ei>=0 && ej>=0;
argumentFactors=Table[
 originalCoefficient=If[coordinate==1,(ei+1) eulerCoefficient[ea,eb,ed,ec,ei+1,ej],(ej+1) eulerCoefficient[ea,eb,ed,ec,ei,ej+1]];
 shiftedCoefficient=If[coordinate==1,eulerCoefficient[ea+1,eb+1,ed,ec+1,ei,ej],eulerCoefficient[ea+1,eb,ed+1,ec+1,ei,ej]];
 pair["item3_kernel/Appell_argument_coefficient_operands",originalCoefficient,shiftedCoefficient];
 factor=trackedFullSimplify[FunctionExpand[originalCoefficient/shiftedCoefficient],Assumptions->coefficientDomain];
 solo["item3_kernel/Appell_argument_coefficient_ratio",factor];
 pair["item3_kernel/Appell_argument_coefficient_identity",originalCoefficient,factor shiftedCoefficient];factor,{coordinate,2}];
Derivative[0,0,0,0,1,0][eulerAppell][av_,bv_,dv_,cv_,xv_,yv_] := (argumentFactors[[1]]/.{ea->av,eb->bv,ed->dv,ec->cv}) eulerAppell[av+1,bv+1,dv,cv+1,xv,yv];
Derivative[0,0,0,0,0,1][eulerAppell][av_,bv_,dv_,cv_,xv_,yv_] := (argumentFactors[[2]]/.{ea->av,eb->bv,ed->dv,ec->cv}) eulerAppell[av+1,bv,dv+1,cv+1,xv,yv];
emit["item3_kernel/Appell_function_definition",HoldForm[eulerAppell[ea,eb,ed,ec,ex,ey]==Inactive[AppellF1][ea,eb,ed,ec,ex,ey]]];
emit["item3_kernel/Appell_notation","Inactive[AppellF1] denotes the standard analytic Appell function with native parameter simplification suppressed. The positive-base Euler continuation is used."];
emit["item3_kernel/Appell_function_domain","All uses have positive integral first parameter and lower parameter one larger. The Euler integral on positive affine bases defines the real analytic branch for every real pair of power parameters."];

affinePrimitiveIntegral[radialExponent_,affineExponent_,orientation_] := affinePrimitiveIntegral[radialExponent,affineExponent,orientation] = Module[{map,bases,constants,arguments},
 stepProgress["finite-interval affine power integration"];
 map=lower+(upper-lower) unitVariable;
 bases={map,distance+orientation map};
 constants=bases/.unitVariable->0;
 arguments=-D[bases,unitVariable]/constants;
 D[map,unitVariable] Times@@MapThread[Power,{constants,{radialExponent,affineExponent}}] eulerAppell[1,-radialExponent,-affineExponent,2,Sequence@@arguments]];
logPrimitive[orientation_] := logPrimitive[orientation] = trackedIntegrate[Log[1+orientation logVariable]/logVariable,logVariable,Assumptions->0<logVariable<1];
Do[solo["item3_kernel/log_primitive_operand",Log[1+orientation logVariable]/logVariable];solo["item3_kernel/log_primitive",logPrimitive[orientation]];
 pair["item3_kernel/log_primitive_derivative",D[logPrimitive[orientation],logVariable],Log[1+orientation logVariable]/logVariable],{orientation,{-1,1}}];
radialLogIntegral[exponent_,orientation_] := Module[{primitive,boundary,generic,special},
 stepProgress["radial logarithm integration by parts"];
 primitive=trackedIntegrate[radial^radialExponent,radial,Assumptions->radial>0,GenerateConditions->False]/.radialExponent->exponent;
 boundary=(primitive Log[distance+orientation radial]/.radial->upper)-(primitive Log[distance+orientation radial]/.radial->lower);
 generic=boundary-orientation affinePrimitiveIntegral[exponent+1,-1,orientation]/(exponent+1);
 special=Log[distance] Log[upper/lower]+(logPrimitive[orientation]/.logVariable->upper/distance)-(logPrimitive[orientation]/.logVariable->lower/distance);
 Piecewise[{{special,exponent==-1}},generic]];
monomialPower[expression_,variable_] := Total[Map[Which[#===variable,1,Head[#]===Power && First[#]===variable,Last[#],True,0]&,If[Head[expression]===Times,List@@expression,{expression}]]];
radialIntegral[expression_] := Module[{answer=0,rpower,rest,logarithms,factors,factor,base,exponent,orientation,coefficient},
 stepProgress["radial term integration"];
 Do[
  rpower=monomialPower[term,radial];
  (* Exponent works for arbitrary real symbolic monomials after expanding only sums. *)
  rest=Cancel[term/radial^rpower];
  logarithms=Cases[rest,_Log,{0,Infinity}];
  If[Length[logarithms]>0,
   factor=First[logarithms];orientation=D[First[factor],radial];coefficient=Cancel[rest/factor];
   emit["item3_kernel/radial_log_factors",{term,rpower,coefficient,factor,orientation}];
   pair["item3_kernel/radial_log_factorization",term,coefficient radial^rpower factor];
   If[!FreeQ[coefficient,radial] || !MemberQ[{-1,1},orientation],Throw[{"radial logarithm factorization",term}]];
   answer+=coefficient radialLogIntegral[rpower,orientation],
   factors=Select[If[Head[rest]===Times,List@@rest,{rest}],!FreeQ[#,radial]&];
   If[Length[factors]>0,
    factor=First[factors];{base,exponent}=If[Head[factor]===Power,List@@factor,{factor,1}];
    orientation=D[base,radial];coefficient=Cancel[rest/factor];
    emit["item3_kernel/radial_affine_factors",{term,rpower,coefficient,base,exponent,orientation}];
    pair["item3_kernel/radial_affine_factorization",term,coefficient radial^rpower base^exponent];
    If[!FreeQ[coefficient,radial] || !MemberQ[{-1,1},orientation] || Expand[base-orientation radial]=!=distance,Throw[{"radial affine factorization",term}]];
    answer+=coefficient affinePrimitiveIntegral[rpower,exponent,orientation],
    answer+=rest Piecewise[{{Log[upper/lower],rpower==-1}},(upper^(rpower+1)-lower^(rpower+1))/(rpower+1)]
   ]
  ],{term,terms[expression]}];answer];

angularKernel[powers_,profileExponent_,volume_] := angularKernel[powers,profileExponent,volume] = Module[{polynomial,answer=0,coefficient,degree,integral,branches,primitive,condition,endpoints,value},
 stepProgress["angular kernel integration "<>ToString[{powers,profileExponent,volume},InputForm]];
 If[NumberQ[profileExponent],
  value=trackedIntegrate[angularMoment[powers] squaredDistance^(profileExponent/2),{mu,-1,1},Assumptions->distance>radial>0];
  Return[If[volume,radialIntegral[Expand[value radial^radialPower]],value]]];
 polynomial=Expand[(angularMoment[powers]/.muRule) distanceJacobian];
 Do[
  degree=First[monomial][[1]];coefficient=Last[monomial];
  branches=Table[
   primitive=branch[[1]]/.power->profilePower/2+degree;
   condition=branch[[2]]/.power->profilePower/2+degree;
   endpoints=Table[PowerExpand[primitive/.squared->bound],{bound,distanceBounds}];
   value=Expand[coefficient (endpoints[[2]]-endpoints[[1]])];
   If[volume,value=radialIntegral[value radial^radialPower]];
   {value,condition},{branch,primitiveBranches}];
  answer+=Piecewise[Most[branches],Last[branches][[1]]],
  {monomial,CoefficientRules[polynomial,{squared}]}];
 If[NumberQ[profileExponent],trackedFullSimplify[answer/.profilePower->profileExponent,Assumptions->distance>radial>0 && Total[axisSymbols^2]==1],answer/.profilePower->profileExponent]];

trigNormalRules={Cos[theta]->normalZ,Sin[theta]->chartSine,Cos[phi]->normalX/chartSine,Sin[phi]->normalY/chartSine};
integratedOrigin[obj_] := integratedOrigin[obj] = Module[{origin=integralOrigins[obj],raw,b,measure,volume,lift,polynomial,answer=0,p,rpower,coefficient,kernel,rules},
 stepProgress["profile integral construction"];
 {raw,b,measure,volume}=origin;
 lift=raw/.Thread[radii^2->{profileDistance^2,radial^2}];
 lift=PowerExpand[lift/.profileRules];
 lift=lift/.Thread[xyz->(centres[[2]]+radial normalSymbols)];
 lift=Expand[PowerExpand[(lift/.separationChart)/.trigNormalRules]];
 measure=If[volume,measure/.s->radial,measure/.b->radial];
 polynomial=Expand[lift measure];
 Do[
  Do[
   p=monomialPower[term,profileDistance];rpower=monomialPower[term,radial];
   coefficient=Cancel[term/(profileDistance^p radial^rpower)];
   emit["item3_kernel/profile_factors",{term,p,rpower,coefficient}];
   pair["item3_kernel/profile_factorization",term,coefficient profileDistance^p radial^rpower];
   If[!FreeQ[coefficient,profileDistance|radial|theta|phi],Throw[{"profile factorization",coefficient}]];
   kernel=angularKernel[First[monomial],p,volume];
   rules=Join[{distance->dist},Thread[axisSymbols->Table[directionSeparation[j],{j,3}]],If[volume,{radialPower->rpower,lower->rin[t],upper->rout[t]},{radial->b}]];
   answer+=coefficient If[volume,1,b^rpower] (kernel/.rules),
   {term,terms[Last[monomial]]}],
  {monomial,CoefficientRules[polynomial,normalSymbols]}];answer];

(* Time chain rule computed from the live separation, with independent jets. *)
separationRate=trackedFullSimplify[D[sep,t]/.separationChart,Assumptions->dist>0 && unitDirection && Element[Table[directionSeparation[j],{j,3}],Reals]];
directionRate=trackedFullSimplify[D[separationVector/sep,t]/.separationChart,Assumptions->dist>0 && unitDirection && Element[Table[directionSeparation[j],{j,3}],Reals]];
solo["item3_kernel/separation_rate",separationRate];solo["item3_kernel/direction_rate",directionRate];
directionDifferentiationSymbols={directionSlot1,directionSlot2,directionSlot3};
directionDifferentiationRules=Thread[Table[directionSeparation[j],{j,3}]->directionDifferentiationSymbols];
directionDifferentiationBack=Reverse/@directionDifferentiationRules;
geometryDistance=trackedFullSimplify[familyObject[geometricDerivative[sep]],Assumptions->dist>0 && unitDirection && Element[Table[directionSeparation[j],{j,3}],Reals]];
geometryDirection=trackedFullSimplify[familyObject[geometricDerivative[#]]& /@ (separationVector/sep),Assumptions->dist>0 && unitDirection && Element[Table[directionSeparation[j],{j,3}],Reals]];
solo["item3_kernel/geometric_distance_derivative",geometryDistance];solo["item3_kernel/geometric_direction_derivative",geometryDirection];
closedGeometricDerivative[value_] := (stepProgress["closed geometric differentiation"];D[value,dist] geometryDistance+((D[value/.directionDifferentiationRules,#]& /@ directionDifferentiationSymbols)/.directionDifferentiationBack).geometryDirection);
closedFamilyQ[obj_] := closedFamilyQ[obj] = Module[{value},Which[
 KeyExistsQ[rateOrigins,obj],stepProgress["closed momentum time differentiation"];value=closedFamilyQ[rateOrigins[obj]];D[value,t]+D[value,dist] separationRate+((D[value/.directionDifferentiationRules,#]& /@ directionDifferentiationSymbols)/.directionDifferentiationBack).directionRate,
 MatchQ[obj,Inactive[Integrate][__]],integratedOrigin[obj],
 Head[obj]===Plus,Total[closedFamilyQ/@(List@@obj)],
 Head[obj]===Times,Times@@(closedFamilyQ/@(List@@obj)),True,familyObject[obj]]];

(* A mathematical chart for the direction identity, derived by intersecting
   a ray with the unit sphere. This does not specialize any physical field. *)
identityRay={planeU rayParameter,planeV rayParameter,1-rayParameter};
identityRayEquation=identityRay.identityRay==1;
solo["item4_kernel/ray_equation",identityRayEquation];
identityRaySolutions=trackedSolve[{identityRayEquation,rayParameter!=0},rayParameter,Reals];
solo["item4_kernel/ray_solutions",identityRaySolutions];
identityChart=identityRay/.First[identityRaySolutions];
solo["item4_kernel/direction_chart",identityChart];
pair["item4_kernel/chart_sphere_identity",Together[identityChart.identityChart],1];
identityChartDenominator=Denominator[Together[identityChart[[1]]]];
solo["item4_kernel/chart_denominator",identityChartDenominator];
identityChartExcluded=trackedReduce[identityChartDenominator==0,{planeU,planeV},Reals];
solo["item4_kernel/chart_excluded_real_parameters",identityChartExcluded];
emit["item4_kernel/chart_condition","planeU and planeV are real. The ray chart omits its pole; a polynomial identity on this dense chart extends to the pole by continuity. All chart parameters remain symbolic."];

(* Convert Piecewise values to exact indicator-weighted polynomials before
   eliminating directions. The indicators retain their full exponent strata. *)
directionPolynomial[e_] := Which[
 FreeQ[e,Alternatives@@directionDifferentiationSymbols],e,
 AtomQ[e],e,
 Head[e]===Piecewise,Module[{prior=False,rows=Append[e[[1]],{e[[2]],True}],weight},
  Total[Table[weight=Boole[row[[2]] && !prior];prior=prior || row[[2]];
   weight directionPolynomial[row[[1]]],{row,rows}]]],
 True,Map[directionPolynomial,e]];

identityDirectionDegree[e_,v_] := Which[FreeQ[e,v],0,e===v,1,
 Head[e]===Plus,Max[identityDirectionDegree[#,v]& /@ List@@e],
 Head[e]===Times,Total[identityDirectionDegree[#,v]& /@ List@@e],
 Head[e]===Power && IntegerQ[Last[e]] && Last[e]>=0,Last[e] identityDirectionDegree[First[e],v],
 Head[e]===Piecewise,Max[identityDirectionDegree[#,v]& /@ Append[e[[1,All,1]],e[[2]]]],
 True,Throw[{"direction degree not polynomial",e,v},identityFailure]];
identityDirectionTotalDegree[e_] := Which[
 FreeQ[e,Alternatives@@directionDifferentiationSymbols],0,
 MemberQ[directionDifferentiationSymbols,e],1,
 Head[e]===Plus,Max[identityDirectionTotalDegree /@ List@@e],
 Head[e]===Times,Total[identityDirectionTotalDegree /@ List@@e],
 Head[e]===Power && IntegerQ[Last[e]] && Last[e]>=0,Last[e] identityDirectionTotalDegree[First[e]],
 Head[e]===Piecewise,Max[identityDirectionTotalDegree /@ Append[e[[1,All,1]],e[[2]]]],
 True,Throw[{"direction total degree not polynomial",e},identityFailure]];
(* Independent item-4 construction: exact direction elimination, a native
   universal identity solve, and polynomial-monomial sufficient branches.
   No zero-multiplicity claim and no Euler-reference substitute for Q. *)
(* Exact parameter branching of a homogeneous linear system. *)
Clear[S15LinearParameterBranches];
S15LinearParameterBranches[inputMatrix_,inputVariables_,inputQ_,exponentVariables_,radiusVariables_,inputSupport_,writeFunction_,progressFunction_] := Module[
 {ratio=Unique["item4RadiusRatio",{Temporary}],scale=Unique["item4RadiusScale",{Temporary}],forward,inverse,domain,matrix,q,support,
  fractions,query,variables,answer,branches,branch,rules,pending,ready,closed,restricted,answers={},gaps={},serial=0,show,step,reconstruction},
 forward=Thread[radiusVariables->{ratio scale,scale}];
 inverse=Thread[{ratio,scale}->{First[radiusVariables]/Last[radiusVariables],Last[radiusVariables]}];
 domain=0<ratio<1 && scale>0 && Element[exponentVariables,Reals];
 show[key_,value_] := (writeFunction[key,value/.inverse];writeFunction[key<>"/comparison","There is nothing to compare it with in this engine."]);
 step[key_] := progressFunction[key];
 show["matrix_parameter_input",{inputMatrix,inputVariables,inputQ,inputSupport}];
 show["radius_coordinate_operands",{radiusVariables,forward,inverse}];
 reconstruction=Together/@(radiusVariables-(radiusVariables/.forward/.inverse));
 show["radius_coordinate_reconstruction_residual",reconstruction];
 If[!And@@(TrueQ[#===0]&/@reconstruction),Return[{{},{{"Radius coordinate reconstruction",reconstruction}}}]];
 step["transform polynomial matrix and Q to radius ratio and scale"];
 matrix=Map[Factor[PowerExpand[#/.forward,Assumptions->domain]]&,inputMatrix,{2}];
 q=PowerExpand[inputQ/.forward,Assumptions->domain];support=inputSupport/.forward;
 show["radius_coordinate_matrix",matrix];show["radius_coordinate_Q",q];
 fractions=NumeratorDenominator[Together[#]]&/@(matrix.inputVariables);
 query=And@@Thread[fractions[[All,1]]==0] && And@@Thread[fractions[[All,2]]!=0] && support && domain;
 variables=Join[Select[exponentVariables,!FreeQ[query,#]&],{ratio,scale},inputVariables];
 show["simultaneous_real_reduction_operands",{query,variables}];
 step["Reduce simultaneous exponents radii and product coordinates"];
 answer=Reduce[query,variables,Reals];show["simultaneous_real_reduction_return",answer];
 If[!FreeQ[answer,Reduce|Resolve],
  show["not_obtained",{"NOT OBTAINED EXACTLY: simultaneous real reduction including exponents returned unevaluated",answer}];
  Return[{{},{{"Simultaneous real reduction including exponents",answer}}}]];
 step["expand returned parameter alternatives"];
 answer=LogicalExpand[answer];branches=If[Head[answer]===Or,List@@answer,{answer}];
 show["simultaneous_real_reduction_alternatives",branches];
 Do[
  serial++;show["parameter_alternative_"<>ToString[serial]<>"/operand",branch];
  If[branch===False,Continue[]];
  rules=Cases[If[Head[branch]===And,List@@branch,{branch}],HoldPattern[Equal[left_,right_]]:>Which[
   Head[left]===Symbol && MemberQ[variables,left] && FreeQ[right,left],left->right,
   Head[right]===Symbol && MemberQ[variables,right] && FreeQ[left,right],right->left,
   True,Nothing],{1}];
  show["parameter_alternative_"<>ToString[serial]<>"/substitution_operand",rules];
  pending=rules;closed={};
  While[pending=!={},
   ready=Select[pending,FreeQ[Last[#],Alternatives@@(First/@pending)]&];
   show["parameter_alternative_"<>ToString[serial]<>"/dependency_ready",ready];
   If[ready==={},Break[]];
   closed=Join[closed,ready];pending=Select[pending,!MemberQ[First/@ready,First[#]]&];
   pending=(First[#]->(Last[#]/.ready))&/@pending;
  ];
  show["parameter_alternative_"<>ToString[serial]<>"/dependency_remainder",pending];
  If[pending=!={},
   show["not_obtained",{"NOT OBTAINED EXACTLY: returned parameter substitution dependencies did not close",branch,pending}];
   AppendTo[gaps,{"Parameter substitution dependencies",branch,pending}];Continue[]];
  show["parameter_alternative_"<>ToString[serial]<>"/substitution_return",closed];
  step["parameter alternative "<>ToString[serial]<>" substitute Q on returned condition"];
  restricted=q/.closed;
  show["parameter_alternative_"<>ToString[serial]<>"/Q_restriction_operands",{q,closed,restricted}];
  show["parameter_alternative_"<>ToString[serial]<>"/Q_restriction_residual",(q/.closed)-restricted];
  show["parameter_alternative_"<>ToString[serial]<>"/condition",branch];
  show["parameter_alternative_"<>ToString[serial]<>"/Q_under_condition",restricted];
  AppendTo[answers,{branch,restricted}],{branch,branches}];
 {answers,gaps}/.inverse
];

(* Item-4 construction using native Euler integration, differential matrices,
   polynomial elimination, and real parameter reduction. No external inputs. *)
Clear[S15NativeItem4];
S15NativeItem4[closedObject_, parameterList_, exponentList_, boundaryList_, distanceVariable_, directionList_, departureFunction_, writeFunction_, progressFunction_] := Module[
 {parameters=parameterList,exponents=exponentList,lo=boundaryList[[1]],hi=boundaryList[[2]],d=distanceVariable,normals=directionList,
 radial=Unique["item4Radial",{Temporary}],domain,baseDomain,sourceBounds=boundaryList,show,step,cas,normalize,directionReduce,factor,exactFactor,requiredConjuncts,reduceSupport,applyParameterRules,powerSubstitute,powerMultiplicity,liftEuler,directionCoefficients,endpointObject,differentialEliminate,primitiveRational,seriesNecessary,radialMoment,radialIntegrate,recordCondition,selectLaurentCoefficient,finiteParameterSystem,elementaryFactorEligible,
 factorAlternatives,parameterBranches,boundaryAlternatives,canonicalSupport,classify,laurentObject,removeTranscendentals,conditions={},unobtained={},counter=0,seen=<||>,integralCache=<||>,strata,partition,primitiveConditions,walkPartition,selectPrimitive,momentPower=Unique["item4MomentPower",{Temporary}],momentGeneric,momentPrimitive,momentPoles,momentSpecial,momentFormula,chartFormula,chartOrigin=Unique["item4ChartOrigin",{Temporary}],chartSlope=Unique["item4ChartSlope",{Temporary}],chartLower=Unique["item4ChartLower",{Temporary}],chartUpper=Unique["item4ChartUpper",{Temporary}],chartVariable=Unique["item4ChartVariable",{Temporary}],cleanupResult},
 domain=d>hi>lo>0 && Element[parameters,Reals];baseDomain=domain;
 factor[value_] := If[FreeQ[value,Root|AlgebraicNumber],Factor[Refine[value,domain]],Cancel[Together[Refine[value,domain]],Extension->Automatic]];
 exactFactor[value_] := If[FreeQ[value,Root|AlgebraicNumber],Factor[Refine[value,domain]],Factor[Refine[value,domain],Extension->Automatic]];
 show[key_,value_] := (writeFunction[key,value];writeFunction[key<>"/comparison","There is nothing to compare it with in this engine."]);
 step[key_] := progressFunction[key];
 SetAttributes[cas,HoldRest];cas[key_,value_] := Module[{result},step[key];result=value;show[key<>"/return",result];result];
 requiredConjuncts[value_] := If[Head[value]===And,List@@value,{value}];
 reduceSupport[value_] := Refine[value /.
  (relation:(_Equal|_Unequal))/;Length[relation]===2 :>
   With[{difference=factor[Together[PowerExpand[relation[[1]]-relation[[2]],Assumptions->domain]]]},
    Head[relation][difference,0]],domain];
 powerMultiplicity[power_,form_] := Module[{pivot,multiple,remainder},
  pivot=SelectFirst[exponents,Coefficient[form,#]=!=0&,Missing["NoAffineCoordinate"]];
  If[MissingQ[pivot],Return[Missing["NoAffineCoordinate"]]];
  multiple=Cancel[Coefficient[power,pivot]/Coefficient[form,pivot]];
  remainder=Expand[power-multiple form];
  If[(IntegerQ[multiple] || Head[multiple]===Rational) && multiple=!=0 && FreeQ[remainder,Alternatives@@exponents],
   {multiple,remainder},Missing["DifferentExponentForm"]]
 ];
 powerSubstitute[value_,target_,replacement_] := If[
  Head[target]===Power && !FreeQ[target[[2]],Alternatives@@exponents],
  value/.Power[base_,power_]/;SameQ[base,target[[1]]] && !MissingQ[powerMultiplicity[power,target[[2]]]]:>
   With[{parts=powerMultiplicity[power,target[[2]]]},replacement^First[parts] base^Last[parts]],
  value/.target->replacement];
 applyParameterRules[value_,substitutions_] := Module[{symbolRules,functionRules,prepared},
  symbolRules=Select[substitutions,MemberQ[parameters,First[#]]&];
  functionRules=Select[substitutions,!MemberQ[parameters,First[#]]&];
  prepared=Fold[ReplaceAll,value,Reverse[symbolRules]];
  prepared=PowerExpand[prepared/.Log[base_]:>Log[factor[base]],Assumptions->domain];
  functionRules=(First[#]->Fold[ReplaceAll,Last[#],Reverse[symbolRules]])& /@ functionRules;
  Fold[powerSubstitute[#1,First[#2],Last[#2]]&,prepared,functionRules]
 ];
 directionReduce[value_] := Module[{table,monomial,reduction,remainder=0,quotient=0,residual},
  show["normalization_direction_operand",value];
  table=CoefficientRules[Expand[value],normals];
  show["normalization_direction_table",table];
  Do[
   monomial=Times@@MapThread[Power,{normals,First[entry]}];
   reduction=cas["PolynomialReduce normalization direction monomial",PolynomialReduce[monomial,{normals.normals-1},normals]];
   show["normalization_direction_monomial_operands",{monomial,reduction}];
   residual=Expand[monomial-reduction[[1]].{normals.normals-1}-reduction[[2]]];
   show["normalization_direction_monomial_residual",residual];
   If[residual=!=0,Throw[{"Direction reduction residual",residual}]];
   remainder+=Last[entry] reduction[[2]];quotient+=Last[entry] First[reduction[[1]]],
  {entry,table}];
  show["normalization_direction_reconstruction_operands",{value,remainder,quotient,normals.normals-1}];
  step["rational direction reconstruction residual"];
  residual=exactFactor[Together[Expand[value-remainder-quotient (normals.normals-1)]]];
  show["normalization_direction_reconstruction_residual",residual];
  If[residual=!=0,Throw[{"Direction reconstruction residual",residual}]];
  Expand[remainder]
 ];
 normalize[value_,canonicalAlgebraicPowers_:False] := Module[{expanded,powers,coordinates,lift,rules,definitions,coefficient,constant,result,prepared,algebraicPowers,powerData,powerKeys,powerSlots,powerDefinitions,powerResiduals},
  prepared=directionReduce[value]/.Power[base_,pow_]/;!FreeQ[pow,Alternatives@@exponents]:>Module[{rat=Together[base],num,den,positive},
   num=Numerator[rat];den=Denominator[rat];
   positive=FullSimplify[{num>0,den>0},domain && lo<radial<hi];
   If[TrueQ[And@@positive],num^Expand[pow] den^(-Expand[pow]),
    positive=FullSimplify[{num<0,den<0},domain && lo<radial<hi];
    If[TrueQ[And@@positive],(-num)^Expand[pow] (-den)^(-Expand[pow]),base^Expand[pow]]]
  ];
  prepared=prepared/.Log[base_]:>Log[factor[base]];
  expanded=Expand[PowerExpand[prepared/.Power[base_,pow_]/;!FreeQ[pow,Alternatives@@exponents] || (!IntegerQ[pow] && !FreeQ[base,d|radial]):>factor[base]^pow,Assumptions->domain && lo<radial<hi]];
  powers=DeleteDuplicates[Cases[expanded,Power[base_,pow_]/;!FreeQ[pow,Alternatives@@exponents]:>base,Infinity]];
  coordinates=Table[Unique["item4Power",{Temporary}],{Length[powers]},{Length[exponents]}];
  definitions=Flatten[Table[coordinates[[i,j]]->powers[[i]]^exponents[[j]],{i,Length[powers]},{j,Length[exponents]}]];
  lift=expanded/.Power[base_,pow_]/;MemberQ[powers,base] && !FreeQ[pow,Alternatives@@exponents]:>Module[{index=First@FirstPosition[powers,base,Missing["absent"],{1}],coefs},
   coefs=Coefficient[pow,#]& /@ exponents;
   (Times@@MapThread[Power,{coordinates[[index]],coefs}]) * base^Expand[pow-coefs.exponents]];
  If[TrueQ[canonicalAlgebraicPowers],
  algebraicPowers=DeleteDuplicates[Cases[lift,
   p:Power[base_,pow_]/;!FreeQ[pow,Root|AlgebraicNumber] && FreeQ[pow,Alternatives@@exponents]:>p,Infinity]];
  step["canonical algebraic power exponents"];
  powerData=Table[With[{base=power[[1]],pow=power[[2]]},
   {power,base,Floor[pow],RootReduce[pow-Floor[pow]],Refine[base>0,domain && lo<radial<hi]}],
  {power,algebraicPowers}];
  show["algebraic_power_coordinate_operands",powerData];
  powerData=Select[powerData,IntegerQ[#[[3]]] && TrueQ[#[[5]]]&];
  powerResiduals=RootReduce[#[[1,2]]-#[[3]]-#[[4]]]&/@powerData;
  show["algebraic_power_exponent_residuals",powerResiduals];
  If[!And@@(TrueQ[#===0]&/@powerResiduals),Throw[{"Algebraic exponent reconstruction",powerResiduals}]];
  powerKeys=DeleteDuplicates[(#[[{2,4}]]&)/@powerData];
  powerSlots=Table[Unique["item4AlgebraicPower",{Temporary}],{Length[powerKeys]}];
  powerDefinitions=MapThread[Rule,{powerSlots,(Power@@#)&/@powerKeys}];
  show["algebraic_power_coordinate_definitions",powerDefinitions];
  lift=lift/.Table[entry[[1]]->entry[[2]]^entry[[3]] powerSlots[[First@FirstPosition[powerKeys,entry[[{2,4}]]]]],{entry,powerData}];
  definitions=Join[definitions,powerDefinitions];
  ];
  step["rational algebra in positive power coordinates"];
  result=If[FreeQ[lift,Root|AlgebraicNumber],factor[Together[lift]],Cancel[Together[Refine[lift,domain]],Extension->Automatic]]/.definitions;
  result
 ];
 liftEuler[value_] := Module[{kernels,slots,chart,integrands,polynomial,coefficients,body,boundary},
  kernels=DeleteDuplicates[Cases[value,HoldPattern[eulerAppell[__]],Infinity]];
  If[kernels==={},Return[{0,value}]];
  slots=Table[Unique["item4EulerIntegral",{Temporary}],{Length[kernels]}];
  If[!ValueQ[chartFormula],
   show["Euler_chart_equation_operands",{chartOrigin+chartSlope chartLower==0,chartOrigin+chartSlope chartUpper==1}];
   chartFormula=(chartOrigin+chartSlope chartVariable)/.First[cas["Solve Euler integration chart",Solve[{chartOrigin+chartSlope chartLower==0,chartOrigin+chartSlope chartUpper==1},{chartOrigin,chartSlope}]]];
   show["Euler_chart_template",chartFormula];
   show["Euler_chart_endpoint_residuals",Together/@{chartFormula/.chartVariable->chartLower,(chartFormula/.chartVariable->chartUpper)-1}]
  ];
  chart=factor[chartFormula/.{chartVariable->radial,chartLower->lo,chartUpper->hi}];show["Euler_chart_substitution",chart];
  integrands=Table[Module[{a,b,c,e,x,y,kernel=kernels[[index]],operand},
   {a,b,c,e,x,y}=List@@kernel;
   operand=a chart^(a-1)(1-x chart)^(-b)(1-y chart)^(-c) D[chart,radial];
   show["Euler_change_of_variable_operands",{kernel,chart,operand}];
   normalize[operand]],{index,Length[kernels]}];
  polynomial=Expand[value/.Thread[kernels->slots]];
  coefficients=Coefficient[polynomial,#]& /@ slots;
  boundary=polynomial/.Thread[slots->0];
  show["Euler_linear_reconstruction_operands",{polynomial,coefficients.slots+boundary}];
  show["Euler_linear_reconstruction_residual",Expand[polynomial-coefficients.slots-boundary]];
  body=normalize[coefficients.integrands];boundary=normalize[boundary];
  show["Euler_radial_form",{body,boundary}];{body,boundary}
 ];
 directionCoefficients[value_] := Module[{polynomial,table,remainder=0,monomial,reduction},
  polynomial=Expand[value];
  show["direction_polynomial_operand",polynomial];
  table=CoefficientRules[polynomial,normals];show["direction_coefficient_table",table];
  Do[
   monomial=Times@@MapThread[Power,{normals,First[entry]}];
   reduction=cas["PolynomialReduce direction monomial",PolynomialReduce[monomial,{normals.normals-1},normals]];
   show["direction_monomial_operands",{monomial,reduction}];
   show["direction_monomial_residual",Expand[monomial-reduction[[1]].{normals.normals-1}-reduction[[2]]]];
   remainder+=Last[entry] reduction[[2]],{entry,table}];
  Last/@CoefficientRules[Expand[remainder],normals]
 ];
 radialMoment[power_] := Module[{},
  If[!ValueQ[momentFormula],
   show["radial_moment_operand",radial^momentPower];
   momentPrimitive=cas["Integrate generic radial moment primitive",Integrate[radial^momentPower,radial,Assumptions->radial>0 && Element[momentPower,Reals],GenerateConditions->True]];
   show["radial_moment_primitive_residual",FullSimplify[D[momentPrimitive,radial]-radial^momentPower,radial>0 && Element[momentPower,Reals]]];
   momentPoles=cas["Solve radial moment primitive divisor",Solve[Denominator[Together[momentPrimitive]]==0,momentPower,Reals]];
   momentGeneric=(momentPrimitive/.radial->hi)-(momentPrimitive/.radial->lo);
   momentSpecial=Table[{cas["Integrate exceptional radial moment",Integrate[radial^momentPower/.rule,{radial,lo,hi},Assumptions->domain,GenerateConditions->True]],And@@(Equal@@#& /@ rule)},{rule,momentPoles}];
   momentFormula=Piecewise[momentSpecial,momentGeneric];show["computed_radial_moment_family",momentFormula]
  ];momentFormula/.momentPower->power
 ];
 radialIntegrate[value_] := Module[{terms,answer=0,power,coefficient},
  terms=If[Head[Expand[value]]===Plus,List@@Expand[value],{Expand[value]}];
  Do[
   step["extract radial monomial power"];
   power=Cancel[Together[PowerExpand[radial D[term,radial]/term,Assumptions->domain && lo<radial<hi]]];
   step["extract radial monomial coefficient"];
   coefficient=Cancel[Together[PowerExpand[term/radial^power,Assumptions->domain && lo<radial<hi]]];
   show["radial_power_decomposition",{term,power,coefficient}];
   show["radial_power_reconstruction_residual",FullSimplify[term-coefficient radial^power,domain && lo<radial<hi]];
   If[FreeQ[{power,coefficient},radial],answer+=coefficient radialMoment[power],
    answer+=cas["Integrate remaining radial coefficient",Integrate[term,{radial,lo,hi},Assumptions->domain,GenerateConditions->True]]],
  {term,DeleteCases[terms,0]}];answer
 ];
 laurentObject[body_,boundary_] := Module[{den,degree,poly,table,integrals,result},
  den=Denominator[Together[body]];
  If[!PolynomialQ[den,d],Return[Missing["NotLaurent"]]];
  degree=Exponent[den,d];
  If[!FreeQ[Cancel[den/d^degree],d] || !PolynomialQ[Together[body d^degree],d],Return[Missing["NotLaurent"]]];
  poly=Expand[Together[body d^degree]];
  table=CoefficientRules[poly,{d}];show["finite_Laurent_table",{degree,table,boundary}];
  show["finite_Laurent_reconstruction_residual",Together[body-Total[(Last[#] d^(First[#][[1]]-degree))& /@ table]]];
  result=boundary;
  Do[
   If[!KeyExistsQ[integralCache,Last[entry]],AssociateTo[integralCache,Last[entry]->radialIntegrate[Last[entry]]]];
   result+=d^(First[entry][[1]]-degree) integralCache[Last[entry]],{entry,table}];
  show["finite_Laurent_departure",result];result
 ];
 selectLaurentCoefficient[value_,assume_] := Module[{scalars,rows={},minimum,polynomial,table,reconstructed,tested,answer=0,usable=True},
  scalars=directionCoefficients[value];
  Do[
   minimum=Min[0,Exponent[Expand[scalar],d,Min]];
   polynomial=Expand[scalar d^(-minimum)];
   If[!PolynomialQ[polynomial,d],usable=False;Break[]];
   table=Reverse[CoefficientRules[polynomial,{d}]];
   table=SortBy[table,-First[#][[1]]&];
   reconstructed=Total[(Last[#] d^(First[#][[1]]+minimum))& /@ table];
   show["integrated_Laurent_coefficient_operands",{scalar,minimum,table}];
   show["integrated_Laurent_coefficient_residual",Together[scalar-reconstructed]];
   rows=Join[rows,Last/@table],{scalar,scalars}];
  If[!usable,Return[value]];
  Do[
   tested=FullSimplify[selectPrimitive[row,assume],domain];
   show["integrated_Laurent_selected_coefficient",{row,tested}];
   If[tested=!=0,answer=tested;Break[]],{row,rows}];
  show["integrated_Laurent_selection_return",answer];answer
 ];
 endpointObject[body_,boundary_] := Module[{h=body,b=boundary,terms,weights,weight,candidate},
  terms=If[Head[Expand[h]]===Plus,List@@Expand[h],{h}];
  weights=DeleteDuplicates[FullSimplify[(d D[#,d]+radial D[#,radial])/# +1,domain && lo<radial<hi]& /@ DeleteCases[terms,0]];
  show["Euler_radial_weights",weights];
  If[!FreeQ[weights,d|radial],Return[Missing["UnresolvedEulerWeight",weights]]];
  Do[step["eliminate Euler integral "<>ToString[weight,InputForm]];
   b=normalize[d D[b,d]-weight b-(radial h/.radial->hi)+(radial h/.radial->lo)];
   h=normalize[d D[h,d]+radial D[h,radial]+(1-weight)h];
   show["Euler_remaining_integrand",h],{weight,weights}];
  show["Euler_endpoint_object",b];
  If[h=!=0,Missing["UnresolvedEulerRemainder",{h,b}],b]
 ];
 removeTranscendentals[value_,assume_] := Module[{functions,next,derivatives,slots,equations,selected,reference,others,returns,coordinate,expanded,index,prepared},
  prepared=value/.Power[base_,pow_]/;!FreeQ[pow,Alternatives@@exponents]:>Module[{rat=Together[base],num,den,positive},
   num=Numerator[rat];den=Denominator[rat];
   positive=FullSimplify[{num>0,den>0},domain && lo<radial<hi];
   If[TrueQ[And@@positive],num^Expand[pow] den^(-Expand[pow]),
    positive=FullSimplify[{num<0,den<0},domain && lo<radial<hi];
    If[TrueQ[And@@positive],(-num)^Expand[pow] (-den)^(-Expand[pow]),base^Expand[pow]]]
  ];
  expanded=Expand[PowerExpand[(prepared/.Log[base_]:>Log[factor[base]])/.Power[base_,pow_]/;!IntegerQ[pow] && !FreeQ[base,d]:>factor[base]^pow,Assumptions->assume]];
  functions=DeleteDuplicates[Cases[expanded,(f:HoldPattern[Log[__]|PolyLog[__]])/;!FreeQ[f,d]:>f,Infinity]];
  If[functions==={},Return[expanded]];
  While[True,
   next=DeleteDuplicates@Join[functions,Flatten[Table[
    step["differentiate logarithmic closure kernel"];
    Cases[Expand[PowerExpand[FullSimplify[D[kernel,d],assume]/.Log[base_]:>Log[factor[base]],Assumptions->assume]],(f:HoldPattern[Log[__]|PolyLog[__]])/;!FreeQ[f,d]:>f,Infinity],{kernel,functions}]]];
   show["logarithmic_differential_closure",next];
   If[next===functions,Break[]];functions=next
  ];
  derivatives={expanded};
  Do[step["differentiate logarithmic identity "<>ToString[index]];
   AppendTo[derivatives,normalize[Expand[PowerExpand[D[Last[derivatives],d]/.Log[base_]:>Log[factor[base]],Assumptions->assume]]]],{index,Length[functions]}];
  slots=Table[Unique["item4LogarithmicCoordinate",{Temporary}],{Length[functions]}];
  equations=Numerator[Together[#]]& /@ (derivatives/.Thread[functions->slots]);
  show["logarithmic_coordinate_definitions",Thread[slots->functions]];
  show["differentiated_identity_operands",derivatives];show["logarithmic_polynomial_operands",equations];
  Do[
   selected=Select[equations,!FreeQ[#,coordinate]&];
   If[selected=!={},
    reference=First[selected];others=Delete[equations,FirstPosition[equations,reference]];
    returns=Table[Module[{rawResult,normalizedResult},
     show["logarithmic_resultant_operands",{reference,equation,coordinate}];
     rawResult=cas["Resultant logarithmic elimination",Resultant[reference,equation,coordinate]];
     step["normalize logarithmic resultant"];
     normalizedResult=If[FreeQ[rawResult,Root|AlgebraicNumber],factor[Together[rawResult]],Cancel[Together[Refine[rawResult,domain]],Extension->Automatic]];
     show["normalized_logarithmic_resultant",normalizedResult];normalizedResult],{equation,others}];
    equations=DeleteCases[returns,0];show["logarithmic_elimination_return",equations]
   ],{coordinate,slots}];
  If[equations==={},Return[Missing["UnresolvedLogarithmicElimination",{value,derivatives}]]];
  First[equations]
 ];
 primitiveRational[value_] := Module[{num,den,top,bottom,content,primitive},
  step["extract parameter content of rational coefficient"];
  {num,den}=NumeratorDenominator[Together[value]];
  top=FactorTermsList[num,d];bottom=FactorTermsList[den,d];
  content=First[top]/First[bottom];primitive=Last[top]/Last[bottom];
  show["rational_coefficient_content",{value,content,primitive}];
  show["rational_content_reconstruction_residual",Together[value-content primitive]];
  {content,primitive}
 ];
 differentialEliminate[coefficients_,basis_,assume_] := Module[{active,exceptions={},pivot,rest,ratio,next,polynomial,index=0,rows,contents,normal,variables,idealRows},
  active=Select[Transpose[{coefficients,basis}],First[#]=!=0&];
  While[Length[active]>1,
   index++;
   normal=Table[contents=primitiveRational[First[entry]];AppendTo[exceptions,Numerator[Together[First[contents]]]];{Last[contents],Last[entry]},{entry,active}];
   active=SortBy[normal,LeafCount[First[#]]&];pivot=First[active];rest=Rest[active];
   polynomial=Numerator[Together[First[pivot]]];
   show["differential_pivot_operand",{pivot,polynomial}];
   rows=CoefficientList[polynomial,d];AppendTo[exceptions,Last[rows]];
   ratio=factor[Together[D[First[pivot],d]/First[pivot]+D[Last[pivot],d]/Last[pivot]]];
   show["differential_pivot_logarithmic_derivative",ratio];
   next=Table[
    step["differential elimination "<>ToString[index]<>" coefficient "<>ToString[j]];
    {factor[Together[D[rest[[j,1]],d]+(D[rest[[j,2]],d]/rest[[j,2]]-ratio) rest[[j,1]]]],rest[[j,2]]},
   {j,Length[rest]}];
   show["differential_elimination_return",next];
   active=Select[next,First[#]=!=0&]
  ];
  If[active==={},Return[Missing["UnresolvedDifferentialRemainder",exceptions]]];
  show["differential_final_coefficient",active];
  polynomial=Numerator[Together[active[[1,1]]]];
  rows=CoefficientList[polynomial,d];show["differential_final_polynomial_coefficients",rows];
  variables=Select[parameters,!FreeQ[rows,#]&];
  show["differential_coefficient_ideal_operands",{rows,variables}];
  If[variables=!={} && And@@(PolynomialQ[#,variables]& /@ rows),
   idealRows=cas["GroebnerBasis differential coefficient ideal",GroebnerBasis[rows,variables]];
   show["differential_coefficient_ideal_selection",idealRows];
   If[ListQ[idealRows] && Length[idealRows]>0,rows=idealRows]
  ];
  Append[exceptions,First@SortBy[Select[rows,#=!=0&],LeafCount]]
 ];
 seriesNecessary[coefficientRows_,basisRows_,assume_] := Module[
  {z=Unique["item4SeriesCoordinate",{Temporary}],x=Unique["item4EndpointCoordinate",{Temporary}],lambda=Unique["item4EulerEigenvalue",{Temporary}],
   denominator,terms={},basis,coefficient,shift,power,table,expanded,groups,group,integerPart,form,offset,maxOffset,weight,
   normals,row,reduced,pivot,echelon,pivots,indices,order,factors={},spectra={},polynomial,resultants,reconstruction,seriesFactor},
  seriesFactor[value_] := Module[{result,parts,variables,collected},
   result=factor[value];
   If[FreeQ[result,Root|AlgebraicNumber],Return[result]];
   parts=NumeratorDenominator[Together[result]];variables=Variables/@parts;
   show["series_algebraic_coefficient_operands",{parts,variables}];
   step["collect exact series algebraic coefficients"];
   collected=MapThread[Collect[#1,#2,Function[c,If[RootReduce[c]===0,0,c]]]&,{parts,variables}];
   show["series_algebraic_coefficient_return",collected];First[collected]/Last[collected]
  ];
  denominator=PolynomialLCM@@(Denominator[Together[#]]& /@ coefficientRows);
  show["series_denominator_operand",{coefficientRows,denominator}];
  Do[
   basis=basisRows[[j]];coefficient=Cancel[coefficientRows[[j]] denominator];
   If[basis===1,shift=0;power=0,
    If[Head[basis]===Power,shift=Expand[basis[[1]]-d];power=Expand[basis[[2]]],shift=Expand[basis-d];power=1]];
   polynomial=Expand[coefficient/.d->x-shift];table=CoefficientRules[polynomial,{x}];
   expanded={shift,Expand[power+First[#][[1]]],Last[#]}& /@ table;
   show["endpoint_polynomial_terms",{basis,coefficient,expanded}];
   reconstruction=Total[#[[3]] (d+#[[1]])^#[[2]]& /@ expanded];
   show["endpoint_polynomial_reconstruction_residual",normalize[basis coefficient-reconstruction]];
   terms=Join[terms,expanded],{j,Length[basisRows]}];
  terms=({#[[1,1]],#[[1,2]],seriesFactor[Total[#[[All,3]]]]}& /@ GatherBy[terms,Take[#,2]&]);
  terms=Select[terms,Last[#]=!=0&];show["endpoint_constant_parameter_coefficients",terms];
  AppendTo[factors,Last[First[SortBy[terms,LeafCount[Last[#]]&]]]];
  terms=Table[integerPart=Floor[term[[2]]/.Thread[exponents->0]];
   If[!IntegerQ[integerPart],integerPart=0];{term[[1]],Expand[term[[2]]-integerPart],integerPart,term[[3]]},{term,terms}];
  groups=GatherBy[terms,#[[2]]&];
  Do[
   form=group[[1,2]];maxOffset=Max[group[[All,3]]];weight=Expand[form+maxOffset];
   normals=z^(maxOffset-#[[3]])(1+#[[1]] z)^(form+#[[3]])& /@ group;
   show["native_series_functions",{weight,normals}];
   echelon={};pivots={};indices={};order=0;
   While[Length[pivots]<Length[normals],
    row=Table[cas["SeriesCoefficient "<>ToString[order]<>" column "<>ToString[j],seriesFactor[FunctionExpand[SeriesCoefficient[normals[[j]],{z,0,order}]]]],{j,Length[normals]}];
    reduced=row;
    Do[step["eliminate series pivot "<>ToString[j]];reduced=seriesFactor[Together[#]]& /@ (reduced-reduced[[pivots[[j]]]] echelon[[j]]),{j,Length[pivots]}];
    show["native_series_row_operands",{order,row,reduced,pivots,echelon}];
    pivot=SelectFirst[Range[Length[reduced]],reduced[[#]]=!=0&,Missing["NoPivot"]];
    If[!MissingQ[pivot],
     AppendTo[factors,Numerator[Together[reduced[[pivot]]]]];
     AppendTo[pivots,pivot];AppendTo[indices,order];AppendTo[echelon,seriesFactor[Together[#/reduced[[pivot]]]]& /@ reduced];
     show["native_series_pivot",{pivot,reduced[[pivot]],Last[echelon]}]
    ];order++
   ];
   polynomial=CharacteristicPolynomial[DiagonalMatrix[weight-indices],lambda];
   show["series_euler_characteristic_polynomial",{weight,indices,polynomial}];AppendTo[spectra,polynomial],
  {group,groups}];
  Do[AppendTo[factors,cas["Resultant Euler spectra "<>ToString[i]<>" "<>ToString[j],seriesFactor[Resultant[spectra[[i]],spectra[[j]],lambda]]]],{i,Length[spectra]},{j,i+1,Length[spectra]}];
  show["native_series_candidate_factors",factors];Clear[seriesFactor];factors
 ];
 canonicalSupport[value_,scope_] := Module[{clauses,unique,slots,indices,lifted,canonical,residual},
  clauses=requiredConjuncts[value];unique=Sort[DeleteDuplicates[clauses]];
  slots=Table[Unique["item4BooleanClause",{Temporary}],{Length[unique]}];
  indices=Table[First@FirstPosition[Map[SameQ[clause,#]&,unique],True],{clause,clauses}];
  lifted=And@@slots[[indices]];canonical=And@@unique;
  show[scope<>"/support_canonicalization_operands",{value,canonical,Thread[slots->unique],lifted,And@@slots}];
  residual=cas[scope<>" BooleanConvert support canonicalization difference",BooleanConvert[Xor[lifted,And@@slots],"DNF"]];
  show[scope<>"/support_canonicalization_residual",residual];
  If[residual=!=False,Throw[{"Support canonicalization difference",residual}]];
  canonical
 ];
 boundaryAlternatives[value_,support_] := Module[{variables=sourceBounds,others,table,degrees,coordinate,relation,substitutions,transformed,scaleTable,polynomial,degree,residual,range,answer,solutions,newRules,newCondition,branches={}},
  others=Complement[parameters,variables];
  If[!PolynomialQ[value,variables] || AnyTrue[variables,FreeQ[value,#]&] || !FreeQ[value,Alternatives@@others],Return[Missing["NotBoundaryPolynomial"]]];
  table=CoefficientRules[value,variables];degrees=DeleteDuplicates[Total/@(First/@table)];
  show["boundary_polynomial_homogeneity_operands",{value,variables,table,degrees}];
  If[Length[degrees]=!=1,Return[Missing["NotHomogeneousBoundaryPolynomial"]]];
  coordinate=Unique["item4BoundaryCoordinate",{Temporary}];relation=coordinate==Last[variables]/First[variables];
  show["boundary_coordinate_relation",relation];
  substitutions=cas["Solve homogeneous boundary coordinate",Solve[relation,Last[variables]]];
  If[!MatchQ[substitutions,{{_Rule}}],Return[Missing["BoundaryCoordinateSubstitution",substitutions]]];
  transformed=Expand[value/.First[substitutions]];scaleTable=CoefficientRules[transformed,{First[variables]}];
  show["boundary_coordinate_polynomial_operands",{value,First[substitutions],transformed,scaleTable}];
  If[Length[scaleTable]=!=1,Return[Missing["BoundaryScaleMonomial",scaleTable]]];
  degree=First[First[First[scaleTable]]];polynomial=Last[First[scaleTable]];
  residual=Together[value-(First[variables]^degree polynomial/.coordinate->Last[variables]/First[variables])];
  show["boundary_coordinate_reconstruction_operands",{value,polynomial,degree,relation}];show["boundary_coordinate_reconstruction_residual",residual];
  If[residual=!=0,Throw[{"Boundary coordinate reconstruction residual",residual}]];
  range=cas["Reduce admissible boundary coordinate",With[{eliminated=DeleteDuplicates@Join[parameters,{d}],constraint=domain && relation},Reduce[Exists[eliminated,constraint],coordinate,Reals]]];
  If[!FreeQ[range,Reduce|Resolve|Exists] || !FreeQ[range,Alternatives@@Join[parameters,{d}]],Return[Missing["BoundaryCoordinateRange",range]]];
  show["boundary_coordinate_root_operands",{polynomial,range}];
  answer=cas["Reduce homogeneous boundary polynomial on coordinate range",Reduce[polynomial==0 && range,coordinate,Reals]];
  If[!FreeQ[answer,Reduce|Resolve],Return[Missing["BoundaryCoordinateRoots",answer]]];
  solutions=cas["Solve reduced boundary coordinate roots",Solve[answer,coordinate,Reals]];
  If[!ListQ[solutions] || !FreeQ[solutions,Solve],Return[Missing["BoundaryCoordinateRootRules",solutions]]];
  Do[
   newRules=First[substitutions]/.solution;
   show["boundary_parameter_substitution_operands",{First[substitutions],solution,newRules}];
   If[!FreeQ[newRules,coordinate],Return[Missing["BoundaryCoordinateRemaining",newRules]]];
   newCondition=Refine[support && And@@(Equal@@#&/@newRules),domain];
   show["boundary_parameter_branch",{newCondition,newRules}];
   If[newCondition=!=False,AppendTo[branches,<|"Condition"->newCondition,"Rules"->newRules|>]],{solution,solutions}];
  branches
 ];
 factorAlternatives[value_] := Module[{factored,numerator,factors,answer},
  show["parameter_factor_input",value];step["factor individual necessary parameter expression"];
  factored=exactFactor[value];numerator=Numerator[Together[factored]];
  factors=If[Head[exactFactor[numerator]]===Times,List@@exactFactor[numerator],{exactFactor[numerator]}];
  show["parameter_factor_product_operands",{numerator,factors}];
  show["parameter_factor_product_residual",Together[numerator-Times@@factors]];
  factors=Replace[#,Power[base_,power_Integer]/;power>0:>base,{0}]& /@ factors;
  show["parameter_factorization_operands",{value,factored,factors}];
  step["individual parameter factorization residual"];show["parameter_factorization_residual",Together[value-factored]];
  DeleteDuplicates[factors]
 ];
 parameterBranches[expressions_,support_] := Module[{factors,remaining=support,answer,branches={},query,known,linear,pivot,coefficient,constant,solutions,elementary,coordinate,lift,target,variables,solved,newCondition,newRules,rawBranches,summands,positive,negative,logarithmic,expandedFactor,boundaryAnswer},
  factors=SortBy[DeleteDuplicates[Flatten[factorAlternatives/@expressions]],{Boole[!MemberQ[Complement[parameters,exponents,sourceBounds],#]],LeafCount[#]}&];
  Do[
   show["parameter_factor_domain_operands",{factor,domain}];step["simplify individual parameter factor domain"];
   query=FullSimplify[factor==0,domain];show["parameter_factor_query",query];
   If[query===False,Continue[]];
   summands=If[Head[Expand[factor]]===Plus,List@@Expand[factor],{Expand[factor]}];
   If[Length[summands]==2,
    positive=Select[summands,TrueQ[FullSimplify[#>0,domain]]&];
    negative=Select[-summands,TrueQ[FullSimplify[#>0,domain]]&];
    If[Length[positive]==1 && Length[negative]==1,
     logarithmic=FullSimplify[PowerExpand[Log[First[positive]]-Log[First[negative]],Assumptions->domain],domain];
     show["positive_power_parameter_operands",{factor,First[positive],First[negative],logarithmic}];
     show["positive_power_parameter_ratio_residual",FullSimplify[Exp[logarithmic]-First[positive]/First[negative],domain]];
     If[logarithmic=!=factor && PolynomialQ[logarithmic,exponents],
      branches=Join[branches,parameterBranches[{logarithmic},remaining]];remaining=remaining && factor!=0;Continue[]]
    ]
   ];
   linear=Select[parameters,PolynomialQ[factor,#] && Exponent[factor,#]===1&];
   linear=SortBy[linear,{Boole[MemberQ[exponents,#]],Boole[MemberQ[sourceBounds,#]]}&];
   solved=False;
   If[linear=!={},
    pivot=First[linear];coefficient=Coefficient[factor,pivot];constant=factor/.pivot->0;
    show["linear_parameter_operands",{factor,pivot,coefficient,constant}];
    show["linear_parameter_reconstruction_residual",Expand[factor-coefficient pivot-constant]];
    solutions=cas["Solve linear named parameter",Solve[factor==0,pivot]];
    Do[newCondition=Refine[remaining && coefficient!=0 && And@@(Equal@@#& /@ solution),domain];If[newCondition=!=False,AppendTo[branches,<|"Condition"->newCondition,"Rules"->solution|>]],{solution,solutions}];
    branches=Join[branches,parameterBranches[{coefficient},remaining]];solved=True
   ];
   If[!solved,
    boundaryAnswer=boundaryAlternatives[factor,remaining];
    show["boundary_parameter_attempt_return",boundaryAnswer];
    If[ListQ[boundaryAnswer],branches=Join[branches,boundaryAnswer];solved=True]
   ];
   If[!solved,
    expandedFactor=PowerExpand[factor,Assumptions->domain];
    show["elementary_parameter_expansion_operands",{factor,expandedFactor}];
    elementary=DeleteDuplicates@Join[Cases[expandedFactor,HoldPattern[Log[arg_]]/;FreeQ[arg,d|radial]:>Log[arg],Infinity],Cases[expandedFactor,Power[base_,pow_]/;MemberQ[sourceBounds,base] && !FreeQ[pow,Alternatives@@exponents]:>base^Expand[pow-(pow/.Thread[exponents->0])],Infinity]];
    Do[
     coordinate=Unique["item4ElementaryValue",{Temporary}];lift=powerSubstitute[expandedFactor,target,coordinate];
     If[PolynomialQ[lift,coordinate] && IntegerQ[Exponent[lift,coordinate]] && Exponent[lift,coordinate]>=1,
      coefficient=Coefficient[lift,coordinate,Exponent[lift,coordinate]];constant=lift/.coordinate->0;
      show["elementary_parameter_operands",{factor,target,lift,coefficient,constant}];show["elementary_parameter_reconstruction_residual",exactFactor[Together[PowerExpand[factor-(lift/.coordinate->target),Assumptions->domain]]]];
      solutions=cas["Solve elementary named-parameter expression",Solve[lift==0,coordinate]];
      Do[newRules=solution/.coordinate->target;newCondition=Refine[remaining && coefficient!=0 && And@@(Equal@@#& /@ newRules),domain];If[newCondition=!=False,AppendTo[branches,<|"Condition"->newCondition,"Rules"->newRules|>]],{solution,solutions}];
      branches=Join[branches,parameterBranches[{coefficient},remaining]];solved=True;Break[]],
    {target,elementary}]
   ];
   If[!solved,
    variables=Select[parameters,!FreeQ[query,#]&];answer=cas["Reduce exact parameter factor",Quiet[Reduce[query,variables,Reals],Reduce::nsmet]];
    If[!FreeQ[answer,Reduce|Resolve],AppendTo[unobtained,{"Real parameter reduction returned unevaluated",query,answer}];show["not_obtained",{"NOT OBTAINED EXACTLY by native real parameter reduction including exponents",query,answer}],
     rawBranches=If[Head[LogicalExpand[answer]]===Or,List@@LogicalExpand[answer],{LogicalExpand[answer]}];
     Do[newCondition=Refine[remaining && raw,domain];
      newRules=Cases[requiredConjuncts[raw],HoldPattern[Equal[left_,right_]]:>Which[Head[left]===Symbol && MemberQ[parameters,left] && FreeQ[right,left],left->right,Head[right]===Symbol && MemberQ[parameters,right] && FreeQ[left,right],right->left,True,Nothing],{1}];
      If[newCondition=!=False,AppendTo[branches,<|"Condition"->newCondition,"Rules"->newRules|>]],{raw,rawBranches}]]
   ];remaining=remaining && factor!=0,
  {factor,factors}];DeleteCases[DeleteDuplicates[branches],False]
 ];
 recordCondition[original_,support_,rules_,body_,boundary_,scope_] := Module[{record,restricted=original,evaluated,symbolRules,consistency,decisions},
  record=support;
  symbolRules=Select[rules,MemberQ[parameters,First[#]]&];
  consistency=Table[Refine[normalize[applyParameterRules[First[rule]-Last[rule],symbolRules]],domain],{rule,Select[rules,!MemberQ[parameters,First[#]]&]}];
  show[scope<>"/parameter_relation_consistency_operands",{rules,symbolRules}];
  show[scope<>"/parameter_relation_consistency_residuals",consistency];
  decisions=Table[If[FreeQ[residual,Alternatives@@parameters],FullSimplify[residual==0,domain],residual==0],{residual,consistency}];
  show[scope<>"/parameter_relation_consistency_return",decisions];
  If[MemberQ[decisions,False],Return[]];
  If[AnyTrue[rules,!MemberQ[parameters,First[#]]&],
   show[scope<>"/Q_parameter_relation_operand",{original,rules}];
   evaluated=laurentObject[body,boundary];
   If[!MissingQ[evaluated],restricted=selectPrimitive[applyParameterRules[evaluated,rules],domain && support]];
   show[scope<>"/Q_parameter_relation_return",restricted]
  ];
  AppendTo[conditions,{record,restricted}];show[scope<>"/condition",record];show[scope<>"/Q_under_condition",restricted]
 ];
 elementaryFactorEligible[factor_] := Module[{query,terms,left,right,decisions,logarithmic,linear,eligible=False},
  query=FullSimplify[factor==0,domain];show["elementary_factor_domain_operand",factor];show["elementary_factor_domain_return",query];
  If[query===False,Return[True]];
  If[MemberQ[parameters,factor],Return[True]];
  linear=Select[exponents,PolynomialQ[factor,#] && Exponent[factor,#]===1 && FreeQ[Coefficient[factor,#],Alternatives@@parameters]&];
  If[linear=!={},Return[True]];
  terms=If[Head[Expand[factor]]===Plus,List@@Expand[factor],{Expand[factor]}];
  If[Length[terms]=!=2,Return[False]];
  Do[
   left=terms[[index]];right=-terms[[3-index]];
   decisions=FullSimplify[{left>0,right>0},domain];show["elementary_factor_log_operands",{left,right,decisions}];
   If[And@@(TrueQ/@decisions),
    logarithmic=FullSimplify[PowerExpand[Log[left]-Log[right],Assumptions->domain],domain];
    show["elementary_factor_log_return",logarithmic];
    show["elementary_factor_log_ratio_residual",FullSimplify[Exp[logarithmic]-left/right,domain]];
    If[PolynomialQ[logarithmic,exponents],eligible=True;Break[]]],{index,2}];eligible
 ];
 finiteParameterSystem[original_,departure_,body_,boundary_,support_,scope_] := Module[
  {qvalue,marker=Unique["item4JointObject",{Temporary}],joint,rows,entry,localCondition,exponentRules,localE,localQ,
   powerDefinitions,logDefinitions,powerSlots,logSlots,bases,logs,steps,positive,coordinateDomain,liftFunction,
   ampVars,ampTable,ampCoordinates,ampDefinitions,liftE,liftQ,jointTable,parts,coefficientList,fractions,equations,
   variables,query,result,branches,branch,substitution,qRestricted,namedCondition,namedQ,residual,expandedCondition,
   allDefinitions,base,index,form,coefficients,constant,amplitudeMarkers,ruleReturn,earlyScalars,earlyCoefficients,earlyCandidate,earlyFactors,earlyEligible,earlyBranches,finiteForms,finiteDecisions,productVariables,productPowers,productSlots,productDefinitions,productEquations,productQ,productMatrix,matrixAnswer},
  qvalue=laurentObject[body,boundary];
  show[scope<>"/joint_Q_integration_return",qvalue];
  If[MissingQ[qvalue] || !FreeQ[qvalue,Integrate],Return[False]];
  finiteForms=NumeratorDenominator[Together[#]]&/@{departure,qvalue};
  finiteDecisions=Table[PolynomialQ[parts[[1]],d] && PolynomialQ[parts[[2]],d] && FreeQ[Cancel[parts[[2]]/d^Exponent[parts[[2]],d]],d],{parts,finiteForms}];
  show[scope<>"/finite_Laurent_test_operands",{departure,qvalue,finiteForms}];show[scope<>"/finite_Laurent_test_return",finiteDecisions];
  If[!And@@finiteDecisions,Return[False]];
  joint=cas[scope<>" PiecewiseExpand Q and departure strata",PiecewiseExpand[PowerExpand[departure+marker qvalue,Assumptions->domain],domain]];
  rows=If[Head[joint]===Piecewise,Append[joint[[1]],{joint[[2]],Not[Or@@joint[[1,All,2]]]}],{{joint,True}}];
  Do[
   show[scope<>"/joint_stratum_operand",entry];
   localCondition=cas[scope<>" Reduce exponent stratum",Reduce[entry[[2]],exponents,Reals]];
   If[localCondition===False,Continue[]];
   exponentRules=cas[scope<>" Solve exponent stratum equalities",Quiet[Solve[And@@Cases[requiredConjuncts[localCondition],HoldPattern[Equal[_,_]],{1}],exponents,Reals],Solve::svars]];
   exponentRules=If[Length[exponentRules]===1,First[exponentRules],{}];
   localE=PowerExpand[(entry[[1]]/.marker->0)/.exponentRules,Assumptions->domain];
   localQ=PowerExpand[Coefficient[entry[[1]],marker]/.exponentRules,Assumptions->domain];
   earlyScalars=directionCoefficients[localE];
   earlyCoefficients=Flatten[Table[With[{minimum=Min[0,Exponent[Expand[value],d,Min]]},Last/@CoefficientRules[Expand[value d^(-minimum)],{d}]],{value,earlyScalars}]];
   earlyCandidate=Missing["NoElementaryParameterFactorization"];
   Do[
    earlyFactors=factorAlternatives[coefficient];
    earlyEligible=elementaryFactorEligible/@earlyFactors;
    show[scope<>"/elementary_coefficient_factorization",{coefficient,earlyFactors,earlyEligible}];
    If[coefficient=!=0 && And@@earlyEligible,earlyCandidate=coefficient;Break[]],{coefficient,earlyCoefficients}];
   If[!MissingQ[earlyCandidate],
    earlyBranches=parameterBranches[{earlyCandidate},support && localCondition];
    show[scope<>"/elementary_coefficient_parameter_return",earlyBranches];
    Do[classify[localQ,support && localCondition && candidate["Condition"],candidate["Rules"]],{candidate,earlyBranches}];Continue[]];
   bases=DeleteDuplicates[Cases[{localE,localQ},Power[base_,power_]/;!FreeQ[power,Alternatives@@exponents] && FreeQ[base,d|radial]:>base,Infinity]];
   powerSlots=Table[Unique["item4ParameterPower",{Temporary}],{Length[bases]},{Length[exponents]}];
   powerDefinitions=Flatten[Table[powerSlots[[i,j]]->bases[[i]]^exponents[[j]],{i,Length[bases]},{j,Length[exponents]}]];
   logs=DeleteDuplicates[Cases[{localE,localQ},HoldPattern[Log[arg_]]/;FreeQ[arg,d|radial]:>Log[arg],Infinity]];
   logSlots=Table[Unique["item4ParameterLog",{Temporary}],{Length[logs]}];logDefinitions=Thread[logSlots->logs];
   liftFunction[value_] := (value/.Power[base_,power_]/;MemberQ[bases,base] && !FreeQ[power,Alternatives@@exponents]:>Module[{i=First@FirstPosition[bases,base,Missing["absent"],{1}],vector,remainder},
    vector=Coefficient[power,#]&/@exponents;remainder=Expand[power-vector.exponents];
    If[FreeQ[remainder,Alternatives@@exponents],base^remainder Times@@MapThread[Power,{powerSlots[[i]],vector}],base^power]])/.Thread[logs->logSlots];
   liftE=liftFunction[localE];liftQ=liftFunction[localQ];
   powerDefinitions=Select[powerDefinitions,!FreeQ[{liftE,liftQ},First[#]]&];
   ampVars=Complement[parameters,exponents,sourceBounds];
   jointTable=CoefficientRules[Expand[liftE+marker liftQ],ampVars];
   ampCoordinates=Table[Unique["item4ParameterMonomial",{Temporary}],{Length[jointTable]}];
   ampDefinitions=MapThread[Rule,{ampCoordinates,(Times@@MapThread[Power,{ampVars,First[#]}])&/@jointTable}];
   parts=Total[MapThread[Times,{ampCoordinates,Last/@jointTable}]];
   liftE=parts/.marker->0;liftQ=Coefficient[parts,marker];
   allDefinitions=Join[ampDefinitions,powerDefinitions,logDefinitions];
   show[scope<>"/joint_coordinate_definitions",allDefinitions];
   show[scope<>"/joint_coordinate_operands",{localE,localQ,liftE,liftQ}];
   show[scope<>"/joint_coordinate_reconstruction_residuals",Together/@{localE-(liftE/.allDefinitions),localQ-(liftQ/.allDefinitions)}];
   coefficientList=directionCoefficients[liftE];
   coefficientList=Flatten[Table[With[{minimum=Min[0,Exponent[Expand[value],d,Min]]},Last/@CoefficientRules[Expand[value d^(-minimum)],{d}]],{value,coefficientList}]];
   fractions=NumeratorDenominator[Together[#]]&/@coefficientList;
   equations=And@@Join[Thread[fractions[[All,1]]==0],Thread[fractions[[All,2]]!=0]];
   variables=Join[ampCoordinates,Select[Join[exponents,Flatten[powerSlots],logSlots],!FreeQ[equations,#]&]];
   positive=Table[FullSimplify[Last[definition]>0,domain],{definition,powerDefinitions}];
   show[scope<>"/joint_coordinate_positivity_operands",powerDefinitions];show[scope<>"/joint_coordinate_positivity_return",positive];
   coordinateDomain=And@@MapThread[If[TrueQ[#1],First[#2]>0,True]&,{positive,powerDefinitions}];
   query=equations && localCondition && hi>lo>0 && coordinateDomain;
   show[scope<>"/joint_parameter_system_operands",{coefficientList,fractions,query,variables}];
   productVariables=Join[ampCoordinates,First/@powerDefinitions,logSlots];
   productPowers=DeleteDuplicates[Flatten[(First/@CoefficientRules[Expand[#],productVariables])&/@Join[fractions[[All,1]],{liftQ}],1]];
   productSlots=Table[Unique["item4ParameterProduct",{Temporary}],{Length[productPowers]}];
   productDefinitions=MapThread[Rule,{productSlots,(Times@@MapThread[Power,{productVariables,#}])&/@productPowers}];
   productEquations=Table[Total[(Last[#]*productSlots[[First@FirstPosition[productPowers,First[#]]]])&/@CoefficientRules[Expand[polynomial],productVariables]],{polynomial,fractions[[All,1]]}];
   productQ=Total[(Last[#]*productSlots[[First@FirstPosition[productPowers,First[#]]]])&/@CoefficientRules[Expand[liftQ],productVariables]];
   show[scope<>"/linear_product_definitions",productDefinitions];
   show[scope<>"/linear_product_operands",{fractions[[All,1]],liftQ,productEquations,productQ}];
   show[scope<>"/linear_product_reconstruction_residuals",Expand/@Join[fractions[[All,1]]-(productEquations/.productDefinitions),{liftQ-(productQ/.productDefinitions)}]];
   productMatrix=Table[Coefficient[productEquations[[i]],productSlots[[j]]],{i,Length[productEquations]},{j,Length[productSlots]}];
   show[scope<>"/linear_product_matrix",productMatrix];
   show[scope<>"/linear_product_matrix_residual",Expand[productEquations-productMatrix.productSlots]];
   matrixAnswer=S15LinearParameterBranches[productMatrix,productSlots,productQ,exponents,sourceBounds,localCondition && And@@Thread[fractions[[All,2]]!=0],Function[{k,v},show[scope<>"/"<>k,v]],step];
   show[scope<>"/linear_parameter_return",matrixAnswer];
   Do[
    namedCondition=support && ((answer[[1]]/.productDefinitions)/.allDefinitions);
    namedQ=(answer[[2]]/.productDefinitions)/.allDefinitions;
    AppendTo[conditions,{namedCondition,namedQ}];show[scope<>"/condition",namedCondition];show[scope<>"/Q_under_condition",namedQ],{answer,matrixAnswer[[1]]}];
   Do[AppendTo[unobtained,{support,localCondition,gap}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY by native exponent-sensitive matrix parameter reduction",gap}],{gap,matrixAnswer[[2]]}];
   Continue[];
,   {entry,rows}];Clear[liftFunction];True
 ];
 classify[value_,givenSupport_,rules_,inherited_:None,inheritedRules_:{}] := Module[{support=givenSupport,key,scope,assume,feasible,restricted,body,boundary,eh,eb,coefficients,endpoint,actual,parts,rows,kernels,slots,lift,bases,ratios,matrix,determinant,numerator,equations,candidateBranches,substitutions,branch,powerSpecifications,offsets,forms,record,primitiveRules,primitiveEqualities,necessaryFactors,candidateExpressions,activeRules,functionRule,rangeKinds,rangeQuery,rangeReturn,preflightRelations,preflightRelation,preflightValue,preflightRejected,branchBounds,branchDomain},
  counter++;scope="branch_"<>ToString[counter];
  support=canonicalSupport[support,scope];
  step[scope<>" release completed CAS cache entries"];ClearSystemCache[];
  assume=domain && support;
  show[scope<>"/parameter_support_operand",assume];
  If[TrueQ[assume===False],Return[]];
  preflightRelations=Select[requiredConjuncts[assume],
   MemberQ[{Unequal,Greater,Less,GreaterEqual,LessEqual},Head[#]] &&
   And@@(FreeQ[Denominator[Together[#]],Alternatives@@parameters]& /@ (List@@#))&];
  preflightRelations=SortBy[preflightRelations,LeafCount];preflightRejected=False;
  Do[
   show[scope<>"/parameter_support_preflight_operand",{preflightRelation,rules}];
   step["Refine parameter support before rule composition"];
   preflightValue=Refine[preflightRelation/.rules,domain];
   show[scope<>"/parameter_support_preflight_return",preflightValue];
   If[TrueQ[preflightValue===False],preflightRejected=True;Break[]],{preflightRelation,preflightRelations}];
  show[scope<>"/parameter_support_preflight_rejection",preflightRejected];
  If[preflightRejected,Return[]];
  activeRules=Join[(First[#]->(applyParameterRules[Last[#],rules]))& /@ Select[inheritedRules,!MemberQ[First/@rules,First[#]]&],rules];
  show[scope<>"/accumulated_parameter_rules",activeRules];
  primitiveEqualities=Select[Cases[requiredConjuncts[support],HoldPattern[Equal[_,_]],{1}],FreeQ[#,Alternatives@@Complement[parameters,exponents]]&];
  primitiveRules=If[primitiveEqualities==={},{},cas["Solve supplied branch equalities",Quiet[Solve[And@@primitiveEqualities,exponents,Reals],Solve::svars]]];
  primitiveRules=If[Length[primitiveRules]==1,First[primitiveRules],{}];
  preflightRelations=Select[requiredConjuncts[assume],
   MemberQ[{Unequal,Greater,Less,GreaterEqual,LessEqual},Head[#]] &&
   And@@(FreeQ[Denominator[Together[#]],Alternatives@@parameters]& /@ (List@@#))&];
  preflightRelations=SortBy[preflightRelations,LeafCount];preflightRejected=False;
  Do[
   show[scope<>"/primitive_parameter_support_preflight_operand",{preflightRelation,primitiveRules}];
   step["Refine parameter support before rule composition"];
   preflightValue=Refine[preflightRelation/.primitiveRules,domain];
   show[scope<>"/primitive_parameter_support_preflight_return",preflightValue];
   If[TrueQ[preflightValue===False],preflightRejected=True;Break[]],{preflightRelation,preflightRelations}];
  show[scope<>"/primitive_parameter_support_preflight_rejection",preflightRejected];
  If[preflightRejected,Return[]];
  activeRules=Join[(First[#]->applyParameterRules[Last[#],primitiveRules])& /@ Select[activeRules,!MemberQ[First/@primitiveRules,First[#]]&],primitiveRules];
  show[scope<>"/feasibility_operand",assume/.activeRules];
  step["local parameter-support substitution"];feasible=reduceSupport[assume/.activeRules];show[scope<>"/support",{support,rules,feasible}];
  If[TrueQ[feasible===False],Return[]];
  Do[
   show[scope<>"/functional_range_operands",{functionRule,domain,feasible}];
   step["FullSimplify functional expression range"];
   rangeKinds=FullSimplify[{First[functionRule]>0,First[functionRule]<0,First[functionRule]==0},domain];
   rangeQuery=And@@MapThread[If[TrueQ[#1],#2,True]&,
    {rangeKinds,{Last[functionRule]>0,Last[functionRule]<0,Last[functionRule]==0}}];
   show[scope<>"/functional_range_query",{rangeKinds,rangeQuery}];
   step["FullSimplify substituted functional range"];
   rangeReturn=FullSimplify[rangeQuery,feasible];show[scope<>"/functional_range_return",rangeReturn];
   feasible=reduceSupport[feasible && rangeReturn];
   show[scope<>"/functional_parameter_domain",feasible];
   If[TrueQ[feasible===False],Break[]],
  {functionRule,Select[activeRules,!MemberQ[parameters,First[#]]&]}];
  If[TrueQ[feasible===False],Return[]];
  restricted=PiecewiseExpand[applyParameterRules[value,Join[activeRules,primitiveRules]],domain];show[scope<>"/Q_substitution",restricted];
  key=Hash[{restricted,support,Sort[activeRules]}];show[scope<>"/iteration_operand",{restricted,support,Sort[rules]}];
  If[KeyExistsQ[seen,key],AppendTo[unobtained,{support,"Native parameter iteration repeated an unchanged input",seen[key]}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY: native exponent-sensitive parameter iteration repeated an unchanged input",seen[key]}];Return[]];AssociateTo[seen,key->scope];
  branchBounds=applyParameterRules[sourceBounds,Join[activeRules,primitiveRules]];
  branchDomain=baseDomain && applyParameterRules[baseDomain,Join[activeRules,primitiveRules]];
  show[scope<>"/radial_interval_substitution_operands",{sourceBounds,baseDomain,Join[activeRules,primitiveRules]}];
  show[scope<>"/radial_interval_substitution_return",{branchBounds,branchDomain}];
  Block[{lo=First[branchBounds],hi=Last[branchBounds],domain=branchDomain,integralCache=<||>,momentFormula,momentPrimitive,momentPoles,momentGeneric,momentSpecial},
  If[inherited===None,{body,boundary}=liftEuler[restricted],
   show[scope<>"/radial_substitution_operands",inherited];{body,boundary}=normalize[applyParameterRules[#,Join[activeRules,primitiveRules]]]& /@ inherited;show[scope<>"/radial_substitution_return",{body,boundary}]];
  eh=normalize[departureFunction[body]];eb=normalize[departureFunction[boundary]];
  show[scope<>"/departure_radial_form",{eh,eb}];
  coefficients=Join[directionCoefficients[eh],directionCoefficients[eb]];
  If[And@@(TrueQ[#===0]& /@ coefficients),
   recordCondition[restricted,support,Join[activeRules,primitiveRules],body,boundary,scope];Return[]];
  actual=laurentObject[eh,eb];
  If[MissingQ[actual],
   endpoint=endpointObject[eh,eb];
   If[MissingQ[endpoint],AppendTo[unobtained,{support,endpoint}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY by native Euler differential elimination",endpoint}];Return[]];
   coefficients=directionCoefficients[endpoint];
   If[And@@(TrueQ[#===0]& /@ coefficients),
    actual=cas["Integrate remaining radial departure",Integrate[eh,{radial,lo,hi},Assumptions->assume,GenerateConditions->True]]+eb
   ]
  ];
  If[!MissingQ[actual],
   show[scope<>"/integrated_departure",actual];
   If[finiteParameterSystem[restricted,actual,body,boundary,support,scope],Return[]];
   actual=selectLaurentCoefficient[applyParameterRules[actual,Join[activeRules,primitiveRules]],assume];
   actual=selectPrimitive[actual,assume];show[scope<>"/integrated_identity_operand",actual];
   If[!FreeQ[actual,Integrate],AppendTo[unobtained,{support,actual}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY by native radial integration",actual}];Return[]];
   If[!FreeQ[actual,Piecewise],
    parts=DeleteDuplicates[Cases[actual,HoldPattern[Piecewise[rows_,default_]]:>rows[[All,2]],Infinity]//Flatten];
    branch=First@Select[parts,!MemberQ[{True,False},FullSimplify[#,assume]]&];
    show[scope<>"/integral_parameter_partition",{branch,Not[branch]}];
    classify[restricted,support && branch,{}, {body,boundary},activeRules];
    classify[restricted,support && Not[branch],{}, {body,boundary},activeRules];Return[]];
   step["normalize integrated direction coefficients"];coefficients=normalize/@directionCoefficients[normalize[actual]];show[scope<>"/integrated_identity_reduction",coefficients];
   If[And@@(TrueQ[#===0]& /@ coefficients),recordCondition[restricted,support,Join[activeRules,primitiveRules],body,boundary,scope];Return[]]
  ];
  endpoint=First@SortBy[Select[coefficients,!TrueQ[#===0]&],LeafCount];show[scope<>"/identity_candidate_operand",endpoint];
  endpoint=removeTranscendentals[endpoint,assume];
  If[MissingQ[endpoint],AppendTo[unobtained,{support,endpoint}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY by native logarithmic resultant elimination",endpoint}];Return[]];
  endpoint=Expand[PowerExpand[endpoint,Assumptions->assume]];
  bases=DeleteDuplicates[Cases[endpoint,Power[base_,power_]/;!FreeQ[base,d] && !IntegerQ[power]:>{base,power},Infinity]];
  offsets=Floor[#[[2]]/.Thread[exponents->0]]& /@ bases;
  offsets=Replace[offsets,x_/;!IntegerQ[x]:>0,{1}];
  powerSpecifications=DeleteDuplicates[MapThread[{First[#1],Expand[Last[#1]-#2]}&,{bases,offsets}]];
  slots=Table[Unique["item4DifferentialKernel",{Temporary}],{Length[powerSpecifications]}];
  lift=endpoint/.Power[base_,power_]/;!FreeQ[base,d] && !IntegerQ[power]:>Module[{offset=Floor[power/.Thread[exponents->0]],index},
   If[!IntegerQ[offset],offset=0];index=First@FirstPosition[powerSpecifications,{base,Expand[power-offset]},Missing["absent"],{1}];slots[[index]] base^offset];
  step["normalize lifted kernel polynomial"];lift=factor[Together[lift]];
  If[slots==={},equations=CoefficientList[Numerator[Together[endpoint]],d];candidateExpressions=None,
   matrix=CoefficientRules[Expand[lift],slots];
   rows=factor[Last[#]]&/@matrix;
   bases=Power@@#&/@powerSpecifications;
   kernels=(Times@@MapThread[Power,{bases,First[#]}])&/@matrix;
   show[scope<>"/differential_kernel_operands",{endpoint,rows,kernels}];record=normalize[endpoint-rows.kernels,True];show[scope<>"/kernel_reconstruction_residual",record];If[record=!=0,AppendTo[unobtained,{support,record}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY: exponent-sensitive kernel reconstruction did not reduce to zero",record}];Return[]];
   forms=Table[If[kernel===1,True,Head[kernel]===Power && FreeQ[Expand[kernel[[1]]-d],d]],{kernel,kernels}];
   show[scope<>"/affine_kernel_operands",{kernels,forms}];
   necessaryFactors=If[And@@forms,seriesNecessary[rows,kernels,assume],differentialEliminate[rows,kernels,assume]];
   If[MissingQ[necessaryFactors],AppendTo[unobtained,{support,necessaryFactors}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY by native differential elimination including exponent parameters",necessaryFactors}];Return[]];
   equations={Times@@necessaryFactors};candidateExpressions=necessaryFactors

  ];
  show[scope<>"/necessary_parameter_polynomials",equations];
  If[And@@(TrueQ[#===0]& /@ equations),
   AppendTo[unobtained,{support,equations}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY: native differential determinant gave no nonzero parameter polynomial",equations}];Return[]];
  show[scope<>"/necessary_parameter_variable_check",FreeQ[equations,Alternatives@@Join[{d,radial},slots]]];
  If[!FreeQ[equations,Alternatives@@Join[{d,radial},slots]],
   AppendTo[unobtained,{support,equations}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY: exponent-sensitive elimination retained an independent variable or kernel coordinate",equations}];Return[]];
  If[candidateExpressions===None,candidateExpressions={First@SortBy[Select[equations,!TrueQ[#===0]&],LeafCount]}];
  show[scope<>"/necessary_parameter_factor_operands",candidateExpressions];
  candidateBranches=parameterBranches[candidateExpressions,support];
  Do[
   branch=candidate["Condition"];substitutions=candidate["Rules"];
   show[scope<>"/candidate_substitutions",{branch,substitutions}];
   If[substitutions==={},AppendTo[unobtained,{support,branch}];show[scope<>"/not_obtained",{"NOT OBTAINED EXACTLY: native parameter reduction gave no substitution",branch}],classify[restricted,support && branch,substitutions,{body,boundary},activeRules]],
  {candidate,candidateBranches}]
  ]
 ];
 selectPrimitive[value_,assume_] := FixedPoint[
  Function[expression,expression/.HoldPattern[Piecewise[rows_,default_]]:>Module[{tested,selected={}},
   tested={First[#],FullSimplify[Last[#],domain && assume]}& /@ rows;
   show["primitive_row_decisions",tested];
   Do[If[Last[row]=!=False,AppendTo[selected,row]];If[TrueQ[Last[row]],Break[]],{row,tested}];
   Piecewise[selected,default]
  ]],value];
 primitiveConditions=DeleteDuplicates[Cases[closedObject,HoldPattern[Piecewise[rows_,default_]]:>rows[[All,2]],Infinity]//Flatten];
 primitiveConditions=DeleteDuplicates[Cases[primitiveConditions,HoldPattern[Equal[left_,right_]|Unequal[left_,right_]]:>factor[left-right],Infinity]];
 show["primitive_parameter_boundaries",primitiveConditions];
 walkPartition[value_,support_,index_] := Module[{boundary,nonzero,zero,variables,solutions,rule,restricted},
  step["partition primitive boundary "<>ToString[index]];
  If[index>Length[primitiveConditions],
   restricted=selectPrimitive[value,support];show["selected_primitive_object",{support,restricted}];
   If[FreeQ[restricted,Piecewise],classify[restricted,support,{}],AppendTo[unobtained,{"Unresolved primitive selection",support,restricted}]];Return[]];
  boundary=primitiveConditions[[index]];
  nonzero=FullSimplify[support && boundary!=0,domain];zero=FullSimplify[support && boundary==0,domain];
  show["primitive_boundary_partition",{support,boundary,nonzero,zero}];
  If[nonzero=!=False,walkPartition[value,nonzero,index+1]];
  If[zero=!=False,walkPartition[value,zero,index+1]]
 ];
 walkPartition[closedObject,True,1];
 cleanupResult={conditions,unobtained};
 Clear[factor,exactFactor,show,step,requiredConjuncts,reduceSupport,powerMultiplicity,powerSubstitute,applyParameterRules,directionReduce,normalize,liftEuler,directionCoefficients,radialMoment,radialIntegrate,laurentObject,selectLaurentCoefficient,endpointObject,removeTranscendentals,primitiveRational,differentialEliminate,seriesNecessary,factorAlternatives,parameterBranches,boundaryAlternatives,canonicalSupport,recordCondition,elementaryFactorEligible,finiteParameterSystem,classify,selectPrimitive,walkPartition,cas];
 ClearSystemCache[];cleanupResult
];

parameterConditionAttempt[name_,object_,q_] := Module[
 {prefix="reading/"<>name<>"/item4",physical,slots,rules,inverse,lift,qlift,
 degrees,totalDegree,powers,values,coefficientSlots,polynomial,sphere,reduction,
 coefficients,query,answer,messages,monomials,polynomialSlots,rows,monomialQuery,
 monomialAnswer,branches,branch,substitutions,operands,residuals,restricted,
 valid,explicit,unresolved=True,conditionIndex=0,coefficientIndex,nativeAnswer,nativeIndex,nativeEmit,attemptFailure},
 physical=DeleteDuplicates@Join[
  Cases[{object,q},Derivative[__][h_][t]/;MemberQ[allowedHeads,h],{0,Infinity}],
  Cases[{object,q},h_[t]/;MemberQ[allowedHeads,h],{0,Infinity}],
  Table[aa[k],{k,2}],Table[bb[k],{k,2}],{rin[t],rout[t]}];
 slots=Table[Unique["item4Parameter"],{Length[physical]}];
 rules=Thread[physical->slots];inverse=Thread[slots->physical];
 lift=(object/.rules)/.directionDifferentiationRules;qlift=q/.rules;
 solo[prefix<>"/parameter_slots",inverse];
 solo[prefix<>"/condition_domain",Element[physical,Reals] && rout[t]>rin[t]>0];
 degrees=identityDirectionDegree[lift,#]& /@ directionDifferentiationSymbols;
 totalDegree=identityDirectionTotalDegree[lift];
 powers=Select[Tuples[Range[0,#]& /@ degrees],Total[#]<=totalDegree&];
 stepProgress["item4 exact direction coefficient differentiation"];
 values=Table[stepProgress["item4 direction coefficient "<>ToString[p,InputForm]];
  (Fold[D[#1,#2]&,lift,Flatten[MapThread[ConstantArray,{directionDifferentiationSymbols,p}]]]/(Times@@(Factorial/@p)))/.Thread[directionDifferentiationSymbols->0],{p,powers}];
 solo[prefix<>"/direction_coefficient_operands",Transpose[{powers,values}]/.inverse];
 coefficientSlots=Table[Unique["directionCoefficient"],{Length[values]}];
 polynomial=Total[MapThread[#1 Times@@MapThread[Power,{directionDifferentiationSymbols,#2}]&,{coefficientSlots,powers}]];
 sphere=Total[directionDifferentiationSymbols^2]-1;
 solo[prefix<>"/formal_polynomial_operand",polynomial];solo[prefix<>"/sphere_operand",sphere];
 reduction=trackedPolynomialReduce[polynomial,{sphere},directionDifferentiationSymbols];
 pair[prefix<>"/direction_reconstruction",polynomial,reduction[[1]].{sphere}+reduction[[2]]];
 coefficients=(Last/@CoefficientRules[reduction[[2]],directionDifferentiationSymbols])/.Thread[coefficientSlots->values];
 solo[prefix<>"/scalar_identity_operands",coefficients/.inverse];
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
  Do[
   conditionIndex++;
   substitutions=Cases[branch,Equal[a_Symbol,b_]/;MemberQ[slots,a] && FreeQ[b,a]:>(a->b),{0,Infinity}];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/parameter_branch",branch/.inverse];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/substitutions",substitutions/.inverse];
   operands=coefficients/.substitutions;
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/identity_operands",coefficients/.inverse];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/substituted_identities",operands/.inverse];
   residuals=trackedFullSimplify[operands,Assumptions->branch && (dist>rout[t]>rin[t]>0/.rules) && Element[slots,Reals]];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/identity_residuals",residuals/.inverse];
   valid=And@@Map[TrueQ[#==0]&,residuals];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/substitution_verification",valid];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/Q_operand",q];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/Q_substituted_operand",(qlift/.substitutions)/.inverse];
   restricted=trackedFullSimplify[qlift/.substitutions,Assumptions->branch && (dist>rout[t]>rin[t]>0/.rules) && unitDirection && Element[slots,Reals]];
   solo[prefix<>"/candidate_"<>ToString[conditionIndex]<>"/Q_restriction_return",restricted/.inverse];
   solo[prefix<>"/condition_"<>ToString[conditionIndex],If[valid,(branch && (rout[t]>rin[t]>0/.rules))/.inverse,Missing["NotObtained"]]];
   solo[prefix<>"/Q_under_condition_"<>ToString[conditionIndex],If[valid,restricted/.inverse,Missing["NotObtained"]]],
   {branch,branches}]
 ];
 solo[prefix<>"/Q_on_family",q];
 emit[prefix<>"/condition_scope","The printed verified parameter branches are sufficient branches obtained by Reduce of the physical monomial equations followed by substitution into every exact scalar identity. They are not asserted to exhaust the conditions."];
 attemptFailure=Environment["S15_WL_ATTEMPT_FAILURE"];
 If[StringQ[attemptFailure] && StringLength[attemptFailure]>0,
  solo[prefix<>"/failed_attempt",attemptFailure];
  emit[prefix<>"/label_kind",If[StringStartsQ[attemptFailure,"User decision 2026-11-06"],"user_decision_2026_11_06","resource_failure"]];
  emit[prefix<>"/status","NOT OBTAINED EXACTLY: "<>attemptFailure];
  AppendTo[exactGaps,{name,"item4",attemptFailure}];Return[Null]];
 nativeEmit[key_,value_] := solo[prefix<>"/exact_attempt/"<>key,value/.inverse/.directionDifferentiationBack];
 nativeAnswer=S15NativeItem4[
  qlift/.directionDifferentiationRules,slots,
  Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}]]/.rules,
  {rin[t],rout[t]}/.rules,dist,directionDifferentiationSymbols,
  Function[value,(closedGeometricDerivative[value]+2 value)/.directionDifferentiationRules],
  nativeEmit,Function[step,stepProgress["item4 exponent-sensitive / "<>step]]];
 solo[prefix<>"/exponent_sensitive_return",nativeAnswer/.inverse/.directionDifferentiationBack];
 Do[
  conditionIndex++;
  solo[prefix<>"/condition_"<>ToString[conditionIndex],
   (nativeAnswer[[1,nativeIndex,1]] && (rout[t]>rin[t]>0/.rules))/.inverse];
  solo[prefix<>"/Q_under_condition_"<>ToString[conditionIndex],
   nativeAnswer[[1,nativeIndex,2]]/.inverse/.directionDifferentiationBack],
 {nativeIndex,Length[nativeAnswer[[1]]]}];
 solo[prefix<>"/unresolved_parameter_branches",nativeAnswer[[2]]/.inverse];
 If[nativeAnswer[[2]]=!={},
  emit[prefix<>"/classification_return_issue", "Unresolved branches returned; no resource-failure label applies."];
  Close[stream];Print["ITEM4_CLASSIFICATION_RETURNED_UNRESOLVED: ",name];Quit[92],
  emit[prefix<>"/status","The exponent-sensitive traversal returned the printed conditions, restricted Q objects, and unresolved-branch list."]];
 Clear[nativeEmit];ClearSystemCache[];

];

(* Compute the exceptional loci of every displayed algebraic exponent divisor,
   with the actual disjoint Piecewise support carried along the traversal. *)
exponentDivisorSupports[e_,support_:True] := Which[
 AtomQ[e],{},
 Head[e]===Piecewise,Module[{prior=False,rows=Append[e[[1]],{e[[2]],True}],guard},
  Flatten[Table[guard=support && row[[2]] && !prior;prior=prior || row[[2]];
   exponentDivisorSupports[row[[1]],guard],{row,rows}],1]],
 Head[e]===Power && IntegerQ[Last[e]] && Last[e]<0 && !FreeQ[First[e],aa|bb],
  Join[{{First[e],support}},Flatten[exponentDivisorSupports[#,support]& /@ List@@e,1]],
 True,Flatten[exponentDivisorSupports[#,support]& /@ List@@e,1]];
exponentPoleLocus[den_,support_] := exponentPoleLocus[den,support] = Module[{rules,operand,answer},
 rules=Thread[Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}]]->{exponentSlot1,exponentSlot2,exponentSlot3,exponentSlot4}];
 operand=(support && den==0)/.rules;
 solo["item3_kernel/exponent_divisor_operand",operand];
 answer=trackedReduce[operand,{exponentSlot1,exponentSlot2,exponentSlot3,exponentSlot4},Reals];
 solo["item3_kernel/exponent_divisor_return",answer];answer/.(Reverse/@rules)];
item3ExponentAudit[name_,expr_] := Module[{pairs,regions},
 pairs=DeleteDuplicates[exponentDivisorSupports[expr]];
 solo["reading/"<>name<>"/item3/exponent_divisor_supports",pairs];
 regions=Table[{pair,exponentPoleLocus@@pair},{pair,pairs}];
 solo["reading/"<>name<>"/item3/exponent_divisor_loci",regions];
 If[AnyTrue[regions,Last[#]=!=False&],
  emit["reading/"<>name<>"/item3/exceptional_status","NOT OBTAINED EXACTLY on the displayed nonempty or unresolved exponent-divisor loci. Exact branch-support Reduce was attempted on every printed algebraic exponent divisor."];
  AppendTo[exactGaps,{name,"item3","Branch-support Reduce of displayed exponent divisors left the printed loci."}]]];

conditionAttempt[name_,object_,q_] := Module[{prefix="reading/"<>name<>"/item4",parameters,slots,lift,rules,inverse,query,answer,messages,branches,condition,substitutions,restricted,directionQuery,directionAnswer,polynomial,reduction,remainder,coefficientRows,coefficientSlots,formalPolynomial,formalAnswer,spherePolynomial,identityReturn,identityConditions,chartPullback,chartNumerator,chartCoefficients},
 If[!FreeQ[object,dist],
  If[s15RunMode === "base",
   emit["ITEM4_WORKER_MANIFEST_ENTRY",name];Return[Null]];
  Return[parameterConditionAttempt[name,object,q]]];
 parameters=DeleteDuplicates[Join[Cases[object,Derivative[__][h_][t]/;MemberQ[allowedHeads,h],{0,Infinity}],Cases[object,h_[t]/;MemberQ[allowedHeads,h],{0,Infinity}],Table[aa[k],{k,2}],Table[bb[k],{k,2}],{rin[t],rout[t]}]];
 slots=Table[Unique["parameter"],{Length[parameters]}];rules=Thread[parameters->slots];inverse=Thread[slots->parameters];
 lift=(object/.rules)/.directionDifferentiationRules;
 stepProgress["direction polynomial construction"];
 polynomial=Expand[directionPolynomial[lift]];
 spherePolynomial=(Subtract@@(List@@unitDirection))/.directionDifferentiationRules;
 solo[prefix<>"/direction_polynomial_operand",polynomial];
 coefficientRows=CoefficientRules[polynomial,directionDifferentiationSymbols];
 coefficientSlots=Table[Unique["directionCoefficient"],{Length[coefficientRows]}];
 formalPolynomial=Total[MapThread[#1 Times@@MapThread[Power,{directionDifferentiationSymbols,First[#2]}]&,{coefficientSlots,coefficientRows}]];
 solo[prefix<>"/formal_direction_polynomial",formalPolynomial];
 solo[prefix<>"/sphere_polynomial_operand",spherePolynomial];
 Block[{$MessageList={}},reduction=Quiet[trackedPolynomialReduce[formalPolynomial,{spherePolynomial},directionDifferentiationSymbols]];messages=$MessageList];
 solo[prefix<>"/polynomial_reduction_return",reduction];emit[prefix<>"/polynomial_reduction_messages",messages];
 If[!MatchQ[reduction,{_List,_}],
  emit[prefix<>"/status","NOT OBTAINED EXACTLY: indicator expansion and formal-coefficient PolynomialReduce were attempted; polynomial reduction did not return a quotient and remainder."];
  AppendTo[exactGaps,{name,"item4","Indicator expansion and formal-coefficient PolynomialReduce did not return a quotient and remainder."}];
  solo[prefix<>"/Q_on_family",q];Return[Null]];
 remainder=reduction[[2]];
 pair[prefix<>"/polynomial_reconstruction",formalPolynomial,reduction[[1]].{spherePolynomial}+remainder];
 chartPullback=Together[remainder/.Thread[directionDifferentiationSymbols->identityChart]];
 chartNumerator=Numerator[chartPullback];
 chartCoefficients=Flatten[CoefficientList[chartNumerator,{planeU,planeV}]];
 solo[prefix<>"/chart_pullback",chartPullback];
 solo[prefix<>"/chart_numerator",chartNumerator];
 solo[prefix<>"/chart_coefficients",chartCoefficients];
 directionQuery=And@@Map[#==0&,chartCoefficients];
 solo[prefix<>"/direction_identity_operand",directionQuery];emit[prefix<>"/parameter_slots",inverse];
 emit[prefix<>"/direction_coefficient_substitution",Thread[coefficientSlots->(Last/@coefficientRows)]];
 Block[{$MessageList={}},formalAnswer=Quiet[trackedReduce[directionQuery,coefficientSlots,Reals]];messages=$MessageList];
 solo[prefix<>"/direction_solver_return",formalAnswer];emit[prefix<>"/direction_solver_messages",messages];
 directionAnswer=formalAnswer/.Thread[coefficientSlots->(Last/@coefficientRows)];
 solo[prefix<>"/direction_elimination_return",directionAnswer];
 If[!FreeQ[directionAnswer,Reduce|Resolve|ForAll|Exists] || !FreeQ[directionAnswer,Alternatives@@directionDifferentiationSymbols],
  emit[prefix<>"/status","NOT OBTAINED EXACTLY: exact indicator expansion, polynomial reduction and Reduce of the rational-chart coefficient equations were attempted; direction elimination remains unresolved."];
  AppendTo[exactGaps,{name,"item4","Indicator expansion, PolynomialReduce and Reduce of the rational-chart coefficient equations returned unresolved."}];
  solo[prefix<>"/Q_on_family",q];Return[Null]];
 query=Implies[dist>rout[t]/.rules,directionAnswer];
 solo[prefix<>"/identity_operand",query];
 If[FreeQ[directionAnswer,dist],
  Block[{$MessageList={}},answer=Quiet[trackedReduce[(rout[t]>rin[t]>0/.rules) && directionAnswer,slots,Reals]];messages=$MessageList],
  (* SolveAlways is asked for the exact identity, including the separation
     domain as an implication. Its unsupported return is kept, not promoted
     to conditions. No numerical separation or exponent is substituted. *)
  Block[{$MessageList={}},identityReturn=Quiet[trackedSolveAlways[query,dist]];messages=$MessageList];
  solo[prefix<>"/identity_solver_return",identityReturn];emit[prefix<>"/identity_solver_messages",messages];
  If[MatchQ[identityReturn,{{(_Rule|_RuleDelayed)...}...}] && messages==={},
   identityConditions=Or@@Map[And@@(Equal@@#& /@ (List@@@#))&,identityReturn];
   solo[prefix<>"/identity_candidate_conditions",identityConditions];
   (* SolveAlways supplies generic solutions; a returned list is not a proof
      of completeness on all exceptional parameter strata. Keep that gap. *)
   answer=identityReturn,
   answer=identityReturn]
 ];
 solo[prefix<>"/solver_return",answer];emit[prefix<>"/solver_messages",messages];
 If[!FreeQ[answer,Reduce|SolveAlways|ForAll|Exists] || !FreeQ[answer,dist] || ListQ[answer],
  emit[prefix<>"/status","NOT OBTAINED EXACTLY: indicator-polynomial direction elimination was computed; SolveAlways on the separation-domain implication, or Reduce for a separation-independent condition, did not return a complete exact parameter classification; SolveAlways lists, if present, are only candidate conditions without a completeness proof."];
  AppendTo[exactGaps,{name,"item4","Indicator-polynomial direction elimination; SolveAlways on the separation-domain implication (Reduce for separation-independent conditions) returned unresolved."}],
  branches=If[Head[LogicalExpand[answer]]===Or,List@@LogicalExpand[answer],{LogicalExpand[answer]}];
  Do[
   condition=branch/.inverse;
   substitutions=Cases[branch,Equal[a_Symbol,b_]/;MemberQ[slots,a] && FreeQ[b,a]:>(a->b),{0,Infinity}]/.inverse;
   restricted=trackedFullSimplify[q/.substitutions,Assumptions->condition && dist>rout[t]>rin[t]>0 && unitDirection];
   solo[prefix<>"/condition",condition];emit[prefix<>"/restriction_substitutions",substitutions];
   pair[prefix<>"/restriction_operands",q,q/.substitutions];solo[prefix<>"/Q_under_condition",restricted],
   {branch,branches}];
  emit[prefix<>"/status","Reduce conditions and Q simplified under each computed condition are printed."]
 ];solo[prefix<>"/Q_on_family",q]];

existenceCondition="SUFFICIENT condition under which the departure exists, not a necessary or complete condition: fields and their responses are C^4 in their displayed arguments on open neighborhoods of the compact surfaces/shell avoiding both cores; centres and radii are C^2 in time; mixed family derivatives exist uniformly with integrable majorants, permitting differentiation under the integrals; pb is C^3 on open intervals containing all used total and base densities. The S15 domain and pointwise pressure conditions hold as well.";
exactGaps={};
emit["Q_manifest",Keys[qObjects]];
Do[
 readingTag="reading/"<>objectName;
 progress[readingTag]; qValue=qObjects[objectName];
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
 familyGeo=familyObject[geoDeparture]; familyQ=familyObject[qValue];
 solo[readingTag<>"/item3/family_integral",familyGeo];
 emit[readingTag<>"/item3/existence_condition",existenceCondition];
 If[hasGeneral[familyGeo],
  familyResult=familyGeo;closedQ=familyQ;
  item3Status="Exact family integral with general fields retained, as permitted by item 3.",
  closedQ=closedFamilyQ[qValue];
  solo[readingTag<>"/item3/Q_closed_operand",closedQ];
  familyResult=closedGeometricDerivative[closedQ]+2 closedQ;
  item3Status="Closed-form return from the angular, distance and finite-interval integration construction. The exponent-branch and retained-function inventories are printed separately."
 ];
 solo[readingTag<>"/item3/object",familyResult];
 If[!hasGeneral[familyGeo] && !FreeQ[familyResult,Inactive[Integrate]|Integrate|ConditionalExpression],
  item3Status="NOT OBTAINED IN CLOSED FORM ON ALL BRANCHES: generating-polynomial, distance-substitution and Appell integration was attempted; its return contains an integral or ConditionalExpression.";
  AppendTo[exactGaps,{objectName,"item3",item3Status}]];
 emit[readingTag<>"/item3/exponent_branch_inventory",DeleteDuplicates[Flatten[Cases[familyResult,HoldPattern[Piecewise[rows_,default_]]:>rows[[All,2]],{0,Infinity}]]]];
 emit[readingTag<>"/item3/retained_Appell_parameters",DeleteDuplicates[Cases[familyResult,eulerAppell[a_,b_,c_,e_,x_,y_]:>{a,b,c,e},{0,Infinity}]]];
 If[!hasGeneral[familyGeo],item3ExponentAudit[objectName,familyResult];item3ExponentAudit[objectName<>"/Q_operand",closedQ]];
 emit[readingTag<>"/item3/status",item3Status];
 If[hasGeneral[familyResult] || !FreeQ[familyResult,Inactive[Integrate]|Integrate|ConditionalExpression],
  emit[readingTag<>"/item4/status","Not applicable to the general-field family integral."];
  solo[readingTag<>"/item4/Q_on_family",familyQ],
  conditionAttempt[objectName,familyResult,closedQ]
 ],{objectName,If[s15RunMode === "item4",{s15SelectedObject},Keys[qObjects]]}];
Scan[emit["EXACT_OBJECT_GAP",#]&,exactGaps];
Switch[s15RunMode,
 "base",emit["END_OF_BASE_PASS","S15 separation departure Wolfram objects"],
 "item4",emit["END_OF_ITEM4_PASS",s15SelectedObject],
 _,emit["END_OF_TRANSCRIPT","S15 separation departure Wolfram objects"]];

, emit["ENGINE_MESSAGES", $MessageList]; Close[stream]; Exit[1]];
Close[stream];
