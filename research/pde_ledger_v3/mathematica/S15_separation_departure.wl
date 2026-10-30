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
output = FileNameJoin[{DirectoryName[$InputFileName], "out", "S15_separation_departure.out"}];
stream = OpenWrite[output];
emit[name_, value_] := (WriteString[stream, name, ": ", ToString[value/.eulerAppell->Inactive[AppellF1], InputForm, PageWidth -> Infinity], "\n"]);
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

(* Exact parameter equations from a regular polynomial differential system.
   This construction uses rational integrand identities and native D closure.
   No other engine's source or result is an input. *)
identityWeight=iv^(ia-1) (1-iv)^(ic-ia-1) (1-ix iv)^(-ib) (1-iy iv)^(-id);
identityLogDerivative=Together[D[identityWeight,iv]/identityWeight];
identityRatios={1,Together[D[identityWeight,ix]/identityWeight],Together[D[identityWeight,iy]/identityWeight]};
identitySecondRatios={Together[D[identityWeight,{ix,2}]/identityWeight],Together[D[identityWeight,{iy,2}]/identityWeight],Together[D[identityWeight,ix,iy]/identityWeight]};
identityClosures=Table[
 identityBoundaryMultiplier=If[coordinate==1,iv (1-iv)/(1-ix iv),iv (1-iv)/(1-iy iv)];
 identityBoundaryRatio=Together[D[identityBoundaryMultiplier,iv]+identityBoundaryMultiplier identityLogDerivative];
 (* The unknowns are bookkeeping coefficients in an algebraic identity. *)
 identityClosureOperand=Together[identityBoundaryRatio-{identityU,identityV,identityW}.identityRatios-identityLast identitySecondRatios[[coordinate]]];
 solo["item4_kernel/IBP_operand",identityClosureOperand];
 identityClosureRules=First[trackedSolve[Thread[CoefficientList[Numerator[identityClosureOperand],iv]==0],{identityU,identityV,identityW,identityLast}]];
 solo["item4_kernel/IBP_solution",identityClosureRules];
 pair["item4_kernel/IBP_identity",identityBoundaryRatio,({identityU,identityV,identityW}.identityRatios+identityLast identitySecondRatios[[coordinate]])/.identityClosureRules];
 Together[-{identityU,identityV,identityW}/identityLast/.identityClosureRules],{coordinate,2}];
solo["item4_kernel/IBP_derivative_closures",identityClosures];
identityCrossOperand=Together[identitySecondRatios[[3]]-{identityU,identityV,identityW}.identityRatios];
identityCrossRules=First[trackedSolve[Thread[CoefficientList[Numerator[identityCrossOperand],iv]==0],{identityU,identityV,identityW}]];
identityCross=Together[{identityU,identityV,identityW}/.identityCrossRules];
solo["item4_kernel/mixed_derivative_identity",identityCross];
pair["item4_kernel/mixed_integrand_identity",identitySecondRatios[[3]],identityCross.identityRatios];
(* A boundary derivative integrates to zero for a>0,c>a and the positive
   affine bases used here. The resulting rational identities are regular in
   both real power parameters. *)
solo["item4_kernel/reference_choice_operand",identityReference-identityBoundary==identityBoundary];
identityReferenceChoice=identityReference/.First[trackedSolve[identityReference-identityBoundary==identityBoundary,identityReference]];
solo["item4_kernel/reference_choice_return",identityReferenceChoice];
identityReferencePoint=identityReferenceChoice/.identityBoundary->rout[t];
solo["item4_kernel/reference_on_domain",identityReferencePoint];
identityEulerOperand=(closedGeometricDerivative[identityQ[dist]]+2 identityQ[dist])/.dist->identityD;
solo["item4_kernel/Euler_equation_operand",identityEulerOperand==0];
identityScalarSolution=casStep["DSolve scalar Euler equation",DSolve[{identityEulerOperand==0,identityQ[identityReference]==identityValue},identityQ[identityD],identityD]];
solo["item4_kernel/Euler_equation_solution",identityScalarSolution];
identityRestrictionTemplate=identityQ[identityD]/.First[identityScalarSolution];
pair["item4_kernel/Euler_equation_residual",identityD D[identityRestrictionTemplate,identityD]+2 identityRestrictionTemplate,0];
pair["item4_kernel/Euler_reference_residual",identityRestrictionTemplate/.identityD->identityReference,identityValue];

Clear[identityPrimitive];
identityPrimitive[atom_,dv_,rv_] := Module[{bv,pv,xv,yv,av,b1v,b2v,cv,sub,xx,yy,xy,jet,mat,zv},
 Which[
 Head[atom]===Power && FreeQ[Last[atom],dv],
  {bv,pv}=List@@atom;{{{Together[pv D[bv,dv]/bv]}},{atom/.dv->rv}},
 Head[atom]===Log,
  bv=First[atom];{{{0,Together[D[bv,dv]/bv]},{0,0}},{atom/.dv->rv,1}},
 Head[atom]===PolyLog && First[atom]==1,
  zv=Last[atom];{{{0,Together[D[zv,dv]/(1-zv)]},{0,0}},{atom/.dv->rv,1}},
 Head[atom]===PolyLog && First[atom]==2,
  zv=Last[atom];{{{0,Together[D[zv,dv]/zv],0},{0,0,Together[D[zv,dv]/(1-zv)]},{0,0,0}},{atom/.dv->rv,PolyLog[1,zv/.dv->rv],1}},
 Head[atom]===eulerAppell || Head[atom]===AppellF1,
  {av,b1v,b2v,cv,xv,yv}=List@@atom;
  sub={ia->av,ib->b1v,id->b2v,ic->cv,ix->xv,iy->yv};
  {xx,yy}=identityClosures/.sub;xy=identityCross/.sub;
  mat=Together[{{0,D[xv,dv],D[yv,dv]},D[xv,dv] xx+D[yv,dv] xy,D[xv,dv] xy+D[yv,dv] yy}];
  jet={eulerAppell[av,b1v,b2v,cv,ix,iy],D[eulerAppell[av,b1v,b2v,cv,ix,iy],ix],D[eulerAppell[av,b1v,b2v,cv,ix,iy],iy]}/.{ix->(xv/.dv->rv),iy->(yv/.dv->rv)};
  {mat,jet},
 True,Throw[{"identity primitive outside computed differential closure",atom},identityFailure]]];

(* Expand only algebraic products and exponent branches. A branch weight is
   carried on a coefficient, never multiplied into a singular generic value. *)
identityTerms[e_,dv_] := Which[
 FreeQ[e,dv],{{{},True,e}},
 Head[e]===Piecewise,Module[{prior=False,rows=Append[e[[1]],{e[[2]],True}],guard},
  Flatten[Table[guard=row[[2]] && !prior;prior=prior || row[[2]];
   ({#[[1]],guard && #[[2]],#[[3]]}& /@ identityTerms[row[[1]],dv]),{row,rows}],1]],
 Head[e]===Plus,Join@@(identityTerms[#,dv]& /@ List@@e),
 Head[e]===Times,Fold[Flatten[Table[{Sort[Join[a[[1]],b[[1]]]],a[[2]] && b[[2]],a[[3]] b[[3]]},{a,#1},{b,#2}],1]&,{{{},True,1}},identityTerms[#,dv]& /@ List@@e],
 Head[e]===Power && IntegerQ[Last[e]],If[PolynomialQ[Numerator[Together[e]],dv] && PolynomialQ[Denominator[Together[e]],dv],{{{},True,e}},{{{e},True,1}}],
 PolynomialQ[e,dv],{{{},True,e}},
 True,{{{e},True,1}}];

finiteIdentity[expression_,dv_,rv_,prefix_] := Module[
 {rows,groups,systems={},starts={},weights={},weightParts={},blocks={},atoms,coefficient,mat,jet,primitive,dims,offsets,size,a,y,p,h=1,den,m,pp,hd,md,pd,g,bound,z,hc,mc,pc,depth,fullSize,idx,ord,entries={},projection={},start={},value,condition,productName,projectionName,startName,poles,poleRules,regularity,shiftedM,shiftedP,coefficientAt},
 stepProgress["finite identity: branch and product decomposition"];
 rows=identityTerms[expression,dv];solo[prefix<>"/decomposition",rows];
 groups=GatherBy[rows,First];
 Do[
  stepProgress["finite identity: differential closure block "<>ToString[block]];
  atoms=groups[[block,1,1]];mat={{0}};jet={1};
  Do[primitive=identityPrimitive[atom,dv,rv];
   mat=KroneckerProduct[mat,IdentityMatrix[Length[primitive[[1]]]]]+KroneckerProduct[IdentityMatrix[Length[mat]],primitive[[1]]];
   jet=Flatten[Outer[Times,jet,primitive[[2]]]],{atom,atoms}];
  coefficient=Table[With[{r=Together[row[[3]]]},
   den=Denominator[r];
   With[{lc=Coefficient[den,dv,Exponent[den,dv]]},
    {Total[Map[(Piecewise[{{Last[#]/lc,row[[2]]}},0] dv^First[#][[1]])&,CoefficientRules[Numerator[r],{dv}]]],Cancel[den/lc]}]],{row,groups[[block]]}];
  AppendTo[weightParts,coefficient];coefficient=Total[First[#]/Last[#]& /@ coefficient];
  AppendTo[systems,mat];AppendTo[starts,jet];AppendTo[weights,coefficient];AppendTo[blocks,{atoms,Length[mat]}],
  {block,Length[groups]}];
 dims=Length/@systems;offsets=Most[FoldList[Plus,0,dims]];size=Total[dims];
 a=SparseArray[Flatten[Table[Map[({First[#][[1]]+offsets[[k]],First[#][[2]]+offsets[[k]]}->Last[#])&,Most[ArrayRules[SparseArray[systems[[k]]]]] ],{k,Length[systems]}],1],{size,size}];
 y=Flatten[starts];p=SparseArray[Table[{1,1+offsets[[k]]}->weights[[k]],{k,Length[weights]}],{1,size}];
 solo[prefix<>"/basis_blocks",blocks];solo[prefix<>"/differential_matrix",a];solo[prefix<>"/output_row",p];solo[prefix<>"/initial_values",y];
 stepProgress["finite identity: polynomial differential system"];
 (* Together treats parameter-only Piecewise expressions as coefficients. *)
 Do[den=Denominator[Together[value]];den=Cancel[den/Coefficient[den,dv,Exponent[den,dv]]];h=PolynomialLCM[h,den],{value,Last/@Most[ArrayRules[a]]}];
 Do[h=PolynomialLCM[h,pair[[2]]],{parts,weightParts},{pair,parts}];
 m=Map[Cancel[h #]&,a,{2}];
 pp=SparseArray[Table[{1,1+offsets[[k]]}->Total[First[#] PolynomialQuotient[h,Last[#],dv]& /@ weightParts[[k]]],{k,Length[weightParts]}],{1,size}];
 hd=Exponent[h,dv];md=Max[0,Exponent[#,dv]& /@ (Last/@Most[ArrayRules[m]])];pd=Max[0,Exponent[#,dv]& /@ (Last/@Most[ArrayRules[pp]])];
 g=Max[0,hd-1,md];bound=size pd+g size (size-1)/2+size-1;
 solo[prefix<>"/polynomial_denominator",h];
 poleRules=trackedSolve[h==0,dv];poles=(dv/.#)& /@ poleRules;
 solo[prefix<>"/denominator_solver_return",poleRules];
 poles=trackedFullSimplify[poles,Assumptions->rout[t]>rin[t]>0];
 regularity=trackedFullSimplify[And@@(#<=rv/2& /@ poles),Assumptions->rout[t]>rin[t]>0];
 solo[prefix<>"/denominator_roots",poles];
 solo[prefix<>"/regular_interval_condition",regularity];solo[prefix<>"/polynomial_matrix",m];solo[prefix<>"/polynomial_output",pp];
 solo[prefix<>"/degree_data",{size,hd,md,pd,g}];solo[prefix<>"/identity_derivative_bound",bound];
 stepProgress["finite identity: polynomial translation at symbolic reference"];
 hc=CoefficientList[h/.dv->rv+z,z];
 shiftedM=Table[If[Mod[k,50]==1,stepProgress["finite identity: matrix polynomial translation "<>ToString[k]]];With[{entry=Most[ArrayRules[m]][[k]]},{First[entry],CoefficientList[Last[entry]/.dv->rv+z,z]}],{k,Length[Most[ArrayRules[m]]]}];
 shiftedP=Table[stepProgress["finite identity: output polynomial translation "<>ToString[k]];With[{entry=Most[ArrayRules[pp]][[k]]},{First[entry],CoefficientList[Last[entry]/.dv->rv+z,z]}],{k,Length[Most[ArrayRules[pp]]]}];
 coefficientAt[list_,n_]:=If[n<Length[list],list[[n+1]],0];
 (* Sparse polynomial entries are shifted separately; no physical Taylor
    approximation is formed. Only finite polynomials are coefficient-extracted. *)
 depth=Max[hd,md+1,pd+1,1];fullSize=size depth;
 Do[
  stepProgress["finite identity: transfer block "<>ToString[lag]];
  Do[value=coefficientAt[Last[entry],lag]/((idx+1) First[hc]);
   AppendTo[entries,{First[entry][[1]],First[entry][[2]]+lag size}->value],{entry,shiftedM}];
  Do[AppendTo[entries,{i,i+lag size}->(-(idx-lag) coefficientAt[hc,lag+1]/((idx+1)First[hc]))],{i,size}],
  {lag,0,depth-1}];
 (* Sum duplicate sparse entries before construction. *)
 entries=Normal[Merge[Association /@ (List /@ entries),Total]];
 Do[AppendTo[entries,{i,i-size}->1],{i,size+1,fullSize}];
 Do[Do[AppendTo[projection,{1,First[entry][[2]]+lag size}->coefficientAt[Last[entry],lag]],{lag,0,depth-1}],{entry,shiftedP}];
 start=Table[{i,1}->y[[i]],{i,Length[y]}];
 solo[prefix<>"/reference_separation",rv];solo[prefix<>"/reference_denominator",First[hc]];
 solo[prefix<>"/transfer_index",idx];solo[prefix<>"/transfer_definition",SparseArray[entries,{fullSize,fullSize}]];
 solo[prefix<>"/projection_definition",SparseArray[projection,{1,fullSize}]];solo[prefix<>"/initial_state_definition",SparseArray[start,{fullSize,1}]];
 solo[prefix<>"/separation_occurrences_in_condition_definitions",Cases[{entries,projection,start},dv,Infinity]];
 productName=Unique["identityTransfer"];projectionName=Unique["identityProjection"];startName=Unique["identityInitial"];
 solo[prefix<>"/matrix_aliases",{productName->{prefix<>"/transfer_definition",{fullSize,fullSize},"index"->idx,"default"->0},projectionName->{prefix<>"/projection_definition",{1,fullSize},"default"->0},startName->{prefix<>"/initial_state_definition",{fullSize,1},"default"->0}}];
 emit[prefix<>"/ordered_product_definition","orderedMatrixProduct[T,j,0,n-1] is T(n-1)...T(0), and the identity matrix when n=0. matrixEntry[M,i,j] denotes entry (i,j) of M. All products and sums are finite; the complete sparse entries are printed above."];
 (* Part extraction is retained explicitly for the scalar of a 1 by 1 product.
    Positive reference powers give every squared summand the same dimension. *)
 condition=Inactive[Sum][Abs[rv^ord matrixEntry[Inactive[Dot][projectionName,orderedMatrixProduct[productName,idx,0,ord-1],startName],1,1]]^2,{ord,0,bound}]==0;
 solo[prefix<>"/condition",condition];And[regularity,condition]];

(* Multiplicity certificate. For H Y'=M Y and f=P Y put
 P_0=P, P_(j+1)=H P_j'+P_j M-j H' P_j. Then H^j f^(j)=P_j Y.
 With m=deg P and g=max(0,deg H-1,deg M), deg P_j<=m+j g.
 If their rational rank is r, the first r rows span all subsequent rows.
 Every minor of their wedge has degree <=r m+g r(r-1)/2. Multiplication
 by an invertible analytic fundamental matrix preserves the minimum order
 of vanishing of that wedge. Its nonzero Wronskians differ by constants.
 A scalar output with zero order k therefore has k<=r m+g r(r-1)/2+r-1.
 The printed bound uses size>=r and remains valid when any parameter causes
 a rank drop. Thus the finite coefficient equations are equivalent to the
 identity, with no division by a parameter-dependent rank minor. *)
emit["item4_kernel/finite_identity_certificate","For the printed polynomial system H Y'=M Y and scalar numerator P Y, the finite multiplicity bound is N deg(P)+max(0,deg(H)-1,deg(M)) N(N-1)/2+N-1. The source gives the wedge/Wronskian proof, including rank drops. Transfer entries are the exact coefficient recurrence of this system at the regular symbolic reference 2 rout[t]. Only finite polynomials are coefficient-extracted. No truncated object or numerical sample replaces a closed form."];
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
finiteConditionAttempt[name_,object_,q_] := Module[{prefix="reading/"<>name<>"/item4",lift,degrees,totalDegree,powers,values,slots,polynomial,sphere,reduction,coefficients,conditions={},result,referenceValue,restricted},
 lift=object/.directionDifferentiationRules;
 degrees=identityDirectionDegree[lift,#]& /@ directionDifferentiationSymbols;
 totalDegree=identityDirectionTotalDegree[lift];
 powers=Select[Tuples[Range[0,#]& /@ degrees],Total[#]<=totalDegree&];
 solo[prefix<>"/direction_total_degree",totalDegree];
 solo[prefix<>"/direction_degree",degrees];
 stepProgress["exact direction coefficient differentiation"];
 values=Table[stepProgress["direction coefficient "<>ToString[p,InputForm]];(Fold[D[#1,#2]&,lift,Flatten[MapThread[ConstantArray,{directionDifferentiationSymbols,p}]]]/(Times@@(Factorial/@p)))/.Thread[directionDifferentiationSymbols->0],{p,powers}];
 solo[prefix<>"/direction_coefficient_operands",Transpose[{powers,values}]];
 slots=Table[Unique["finiteDirectionCoefficient"],{Length[values]}];
 polynomial=Total[MapThread[#1 Times@@MapThread[Power,{directionDifferentiationSymbols,#2}]&,{slots,powers}]];
 sphere=Total[directionDifferentiationSymbols^2]-1;
 solo[prefix<>"/formal_polynomial_operand",polynomial];solo[prefix<>"/sphere_operand",sphere];
 reduction=trackedPolynomialReduce[polynomial,{sphere},directionDifferentiationSymbols];
 pair[prefix<>"/direction_reconstruction",polynomial,reduction[[1]].{sphere}+reduction[[2]]];
 coefficients=(Last/@CoefficientRules[reduction[[2]],directionDifferentiationSymbols])/.Thread[slots->values];
 solo[prefix<>"/scalar_identity_operands",coefficients];
 Do[result=Block[{currentObject=prefix<>"/scalar_"<>ToString[k]},Catch[finiteIdentity[coefficients[[k]],dist,identityReferencePoint,prefix<>"/scalar_"<>ToString[k]],identityFailure]];
  If[Head[result]===List,solo[prefix<>"/attempt_return",result];
   emit[prefix<>"/status","NOT OBTAINED EXACTLY: differential closure and finite multiplicity elimination were attempted; the displayed unsupported primitive or algebraic step prevented completion."];
   AppendTo[exactGaps,{name,"item4",result}];Return[Null]];
  AppendTo[conditions,result],{k,Length[coefficients]}];
 solo[prefix<>"/condition",rout[t]>rin[t]>0 && And@@conditions];
 referenceValue=q/.dist->identityReferencePoint;
 solo[prefix<>"/restriction_reference_operand",referenceValue];
 restricted=identityRestrictionTemplate/.{identityD->dist,identityReference->identityReferencePoint,identityValue->referenceValue};
 solo[prefix<>"/Q_under_condition",restricted];solo[prefix<>"/Q_on_family",q];
 emit[prefix<>"/status","Finite parameter equations computed from the closed-form differential system and exact direction coefficient elimination. Q under that condition is computed from the Euler equation with its symbolic reference value. Branches are retained in the printed coefficient and initial-state definitions."]];

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
 If[!FreeQ[object,dist],Return[finiteConditionAttempt[name,object,q]]];
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
 ],{objectName,Keys[qObjects]}];
Scan[emit["EXACT_OBJECT_GAP",#]&,exactGaps];
emit["END_OF_TRANSCRIPT","S15 separation departure Wolfram objects"];

, emit["ENGINE_MESSAGES", $MessageList]; Close[stream]; Exit[1]];
Close[stream];
