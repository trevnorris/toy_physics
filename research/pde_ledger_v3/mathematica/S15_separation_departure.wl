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
emit[name_, value_] := (WriteString[stream, name, ": ", ToString[value, InputForm, PageWidth -> Infinity], "\n"]);
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
affinePrimitiveIntegral[radialExponent_,affineExponent_,orientation_] := affinePrimitiveIntegral[radialExponent,affineExponent,orientation] = Module[{map,bases,constants,arguments},
 stepProgress["finite-interval affine power integration"];
 map=lower+(upper-lower) unitVariable;
 bases={map,distance+orientation map};
 constants=bases/.unitVariable->0;
 arguments=-D[bases,unitVariable]/constants;
 D[map,unitVariable] Times@@MapThread[Power,{constants,{radialExponent,affineExponent}}] AppellF1[1,-radialExponent,-affineExponent,2,Sequence@@arguments]];
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

conditionAttempt[name_,object_,q_] := Module[{prefix="reading/"<>name<>"/item4",parameters,slots,lift,rules,inverse,query,answer,messages,branches,condition,substitutions,restricted,directionQuery,directionAnswer,polynomial,reduction,remainder,coefficientRows,coefficientSlots,formalPolynomial,formalAnswer,spherePolynomial,identityReturn,identityConditions,chartPullback,chartNumerator,chartCoefficients},
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
  item3Status="Closed form from generating-polynomial angular moments, distance substitution and finite-interval Appell integration, including logarithmic branches."
 ];
 solo[readingTag<>"/item3/object",familyResult];
 If[!hasGeneral[familyGeo] && !FreeQ[familyResult,Inactive[Integrate]|Integrate|ConditionalExpression],
  item3Status="NOT OBTAINED IN CLOSED FORM ON ALL BRANCHES: generating-polynomial, distance-substitution and Appell integration was attempted; its return contains an integral or ConditionalExpression.";
  AppendTo[exactGaps,{objectName,"item3",item3Status}]];
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
