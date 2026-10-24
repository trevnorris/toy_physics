SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
(* Independent S15 construction from the supplied equations. No external input.
   Boolean-lattice finite differences partition polynomial terms by label support.
   Integrals are retained as exact inactive integrals of general live functions.
   No numerical or symbolic value is an acceptance criterion. *)
$HistoryLength = 0;
$MaxExtraPrecision = 0;
output = FileNameJoin[{DirectoryName[$InputFileName], "out", "S15_two_drain_region.out"}];
stream = OpenWrite[output];
emit[name_, value_] := (WriteString[stream, name, ": ", ToString[value, InputForm, PageWidth -> Infinity], "\n"]);
progress[name_] := (Print["WORKING_ON: ", name]; emit["WORKING_ON", name]);
solo[name_, value_] := (emit[name, value]; emit[name <> "/comparison", "There is nothing to compare it with in this engine."]);
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
metricMeasure[map_,coordinates_,assumptions_] := FullSimplify[Sqrt[Det[Table[D[map,q].D[map,r],{q,coordinates},{r,coordinates}]]],Assumptions->assumptions];
assumptions = 0 < theta < Pi && Element[phi,Reals] && rin[t]>0 && rout[t]>0 && s>0;
volumeMap = chart[s];
volumeMeasure = metricMeasure[volumeMap,{s,theta,phi},assumptions];
pull[e_,b_] := e /. Thread[xyz -> chart[b]];
sint[e_,b_,measure_] := Inactive[Integrate][pull[e,b] measure,{theta,0,Pi},{phi,0,2 Pi}];
vint[e_] := Inactive[Integrate][pull[e,s] volumeMeasure,{s,rin[t],rout[t]},{theta,0,Pi},{phi,0,2 Pi}];
(* Generic Leibniz rule on the radial member of the chart integral. *)
rate[e_] := Module[{integrand=pull[e,s] volumeMeasure,limits={rin[t],rout[t]},sgn={-1,1}},
 Inactive[Integrate][D[integrand,t],{s,limits[[1]],limits[[2]]},{theta,0,Pi},{phi,0,2 Pi}] +
 Total[Table[sgn[[j]] D[limits[[j]],t] Inactive[Integrate][integrand /. s->limits[[j]],{theta,0,Pi},{phi,0,2 Pi}],{j,2}]]];

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
 dimensionAudit[tag<>"/balance_rhs",rhs];
 emit[tag<>"/balance_dimension_operand_a",dim[lhs]];emit[tag<>"/balance_dimension_operand_b",dim[rhs]];
 emit[tag<>"/balance_dimension_residuals",Flatten[Table[da-db,{da,dim[lhs]},{db,dim[rhs]}],1]],
 {i,3}];
emit["balance_construction", "Leibniz differentiation of the moving spherical-chart momentum integral supplies the rate. Integrating the supplied local balance and applying divergence and moving-domain transport gives the separately printed surface parts plus undivided force. The residual is printed as an equation operand, not imposed on arbitrary independent live fields."];
emit["scope", "No force law, coefficient, normalization, mass, throat-force identification, or relation to J, Q, Theta_Q or N_infinity_3 is computed."];
emit["END_OF_TRANSCRIPT", "S15 Wolfram objects"];
, emit["ENGINE_MESSAGES", $MessageList]; Close[stream]; Exit[1]];
Close[stream];
