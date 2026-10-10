(* S9b blind builder.  Sole physical input: S9b_SHARED_PHYSICS.md v10, Parts A-C only.
   No Get, Needs, Import, file access, exports, or runtime configuration.
   General radial functions; no profile-family restriction.  All expressions
   below are formal coefficients about the flat ray; no straight-ray
   substitution is made for the higher retained grades of the travel time.
   Units: V is the contravariant radial coordinate velocity of the material.
   Signed bending is the central scattering angle (positive toward the mass).
   Build demonstrations, not production transcripts, belong in _scratch/.
*)
ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};
emitted = <||>; localNames = {};
emit[name_String, object_] := Module[{tag = If[StringStartsQ[name, "LOCAL_"],
    "WL_LOCAL_S9B_" <> StringDrop[name, 6], "WL_S9B_" <> name], stream},
  If[KeyExistsQ[emitted, tag], Quit[91]];
  AssociateTo[emitted, tag -> object];
  If[StringStartsQ[name, "LOCAL_"], AppendTo[localNames, tag]];
  stream = First[$Output];
  WriteString[stream, tag <> ": " <>
    ToString[object, InputForm, PageWidth -> Infinity] <> "\n"];
  Flush[stream];
];

(* Retained-ring operations.  xV counts V, not V^2. *)
markers = {xD, xV, xW}; bounds = {1, 2, 1};
grades = Tuples[Range[0, #] & /@ bounds];
evenGrades = Select[grades, EvenQ[#[[2]]] &];
nonzeroEven = DeleteCases[evenGrades, {0, 0, 0}];
firstGrades = {{1, 0, 0}, {0, 1, 0}, {0, 2, 0}, {0, 0, 1}};
mon[g_] := Times @@ MapThread[Power, {markers, g}];
coeff[z_, g_] := Fold[Coefficient[#1, #2[[1]], #2[[2]]] &,
  Expand[z], Transpose[{markers, g}]];
box[z_] := Expand[Fold[Normal[Series[#1, {#2[[1]], 0, #2[[2]]}]] &,
  z, Reverse[Transpose[{markers, bounds}]]]];
project[z_] := Total[(coeff[z, #] mon[#]) & /@ grades];
label[g_] := StringJoin[ToString /@ g];
atOne[z_] := z /. Thread[markers -> 1];

$Assumptions = Element[{c0, ell, b, bFar, ZE, ZR, r, rr, impact,
  GM, n, K, m, s, rho0, Phi, rAnchor}, Reals] && c0 > 0 && ell > 0 &&
  b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 &&
  rho0 > 0 && m > 0 && rAnchor > bFar;
(* Supplied action and general radial ansatz. *)
delta = deltaProfile[r]; velocity = velocityProfile[r];
xiAnsatz = xiProfile[r];
xiSlope = Simplify[D[xiAnsatz, r], $Assumptions]; (* K4c *)
densityFraction = fProfile[r];
returnDirection = -1; (* K5 *)
speed = c0 (1 + xD delta);
radialMetric = 1 + xW xiSlope^2;
metric3 = DiagonalMatrix[{radialMetric, r^2, r^2 Sin[theta]^2}];
metric2 = DiagonalMatrix[{aMetric, r^2}]; (* K2 *)
kCov = {kr, kphi}; vVector = {vLocal, 0}; (* K11 *)
kineticSymbol = (omega - vVector.kCov)^2; (* K1 *)
dispersionSpeed = cSquared; (* K3 *)
dispersion = kineticSymbol - dispersionSpeed kCov.Inverse[metric2].kCov;
localSpeedIdentification = cSquared == muPerp[r]/rhoBr[r];
farSpeedIdentification = c0^2 == muPerpInfinity/rhoBrInfinity;
massBalanceInput = divRhoV == -jn[r];
referenceTheta = (1 + gamma) 2 GM/(b c0^2);
referenceRadar = 2 (1 + gamma) GM/c0^3 Log[4 rE rR/b^2];
bulkPressure = K rho^n;
bulkSoundSquared = D[bulkPressure, rho]/m;
localSoundSquared = bulkSoundSquared /. rho -> rho0 (1 + fLocal); (* K9a *)
localBulkDensity = rho0 (1 + fLocal); (* K9b *)
responseInputs = {c0, c0 Sqrt[localSoundSquared/(bulkSoundSquared /. rho -> rho0)],
  c0 (localBulkDensity/rho0)^s};
rayDerivative[z_, x_] := D[z, x]; (* K4a,K4b: freeze only at differentiation *)
rayDerivative[z_, x_, j_Integer] := Nest[rayDerivative[#, x] &, z, j];

emit["LOCAL_SUPPLIED_INPUTS", <|"Dispersion" -> dispersion == 0,
  "Speed" -> localSpeedIdentification, "FarSpeed" -> farSpeedIdentification,
  "Metric" -> metric3, "Balance" -> massBalanceInput,
  "References" -> {referenceTheta, referenceRadar, gamma == 1},
  "Pressure" -> bulkPressure, "Responses" -> responseInputs,
  "Embedding" -> {w == xiLive[r], xiLive[r] == ell h[r]},
  "Polarizations" -> {cGamma1[r] == cGamma[r], cGamma2[r] == cGamma[r]},
  "ReferenceSpeed" -> (c0 == Inactive[Limit][cGamma[r], r -> Infinity]),
  "UniformAnchor" -> {rhoBrInfinity omega^2 == muPerpInfinity k^2,
    vUniform == 0, xiUniform == 0, gradMuUniform == 0, gradRhoUniform == 0},
  "BackgroundAnchoring" -> {qBackgroundLab[x,t] == qBackground[x],
    qBackgroundMaterial[x,t] == qBackground[chi[x,t]], cGammaLab[x,t] == cGamma[x]},
  "SelectedAnchoring" -> "LAB_HELD"|>];
emit["LOCAL_PROFILE_ANSATZ", <|"Delta" -> delta, "V" -> velocity,
  "XiDerivative" -> xiSlope, "RhoBr" -> rhoBr[r],
  "FractionalBulkDensity" -> densityFraction, "Domain" -> $Assumptions,
  "Counting" -> {"delta O(epsilon)", "V/c0 O(Sqrt[epsilon])", "xiPrimeSquared O(epsilon)"},
  "Anchoring" -> "LAB_HELD"|>];
emit["LOCAL_RETAINED_GRADES", <|"Markers" -> markers, "Grades" -> grades,
  "EpsilonPowers" -> (#.{1, 1/2, 1} & /@ grades)|>];

(* Hamilton velocity and the Fermat functional are derived from the symbol.
   The positive comoving root is selected by its value at vLocal=0.
*)
localAssumptions = c0 > 0 && aMetric > 0 && cSquared > 0 && r > 0 && kr > 0 &&
  Element[{aMetric, cSquared, r, kr, kphi, vLocal}, Reals];
frequencyRoots = omega /. Solve[dispersion == 0, omega];
positiveRoot = Select[frequencyRoots,
  TrueQ[Simplify[(# /. vLocal -> 0) > 0, localAssumptions]] &][[1]];
groupVelocity = D[positiveRoot, #] & /@ kCov;
advectedVelocity = Coefficient[positiveRoot, vLocal] vLocal;
constructedDrift = D[advectedVelocity, #] & /@ kCov;
relativeNorm = Simplify[(groupVelocity - constructedDrift).metric2.
  (groupVelocity - constructedDrift), localAssumptions];
tangent = {dr, dphi};
timePolynomial = Expand[(tangent - tau constructedDrift).metric2.
  (tangent - tau constructedDrift) - tau^2 relativeNorm];
timeCoefficients = Coefficient[timePolynomial, tau, #] & /@ Range[0, 2];
timeDiscriminant = Simplify[timeCoefficients[[2]]^2 -
  4 timeCoefficients[[1]] timeCoefficients[[3]]];
timeRoot = (-timeCoefficients[[2]] - Sqrt[timeDiscriminant])/
  (2 timeCoefficients[[3]]);
subcritical = Simplify[-timeCoefficients[[3]] > 0];
fermatOdd = Simplify[(timeRoot - (timeRoot /. vLocal -> -vLocal))/2,
  localAssumptions && subcritical];
fermatEvenSquared = Simplify[((timeRoot + (timeRoot /. vLocal -> -vLocal))/2)^2,
  localAssumptions && subcritical];
opticalMetric = Simplify[Table[D[fermatEvenSquared, tangent[[i]], tangent[[j]]]/2,
  {i, 2}, {j, 2}], localAssumptions && subcritical];
radialSpeeds = Simplify[(D[#, kr] & /@ frequencyRoots) /. kphi -> 0,
  localAssumptions];
propagatingLocal = relativeNorm > 0 && Det[metric2] > 0;
branchGateLocal = propagatingLocal; (* K10 *)
traversalLocal = Simplify[And @@ {Min @@ radialSpeeds < 0,
  Max @@ radialSpeeds > 0}, localAssumptions];
localRules = {aMetric -> radialMetric, cSquared -> speed^2,
  vLocal -> xV velocity};
physicalPropagation = atOne[branchGateLocal /. localRules];
physicalTraversal = atOne[subcritical /. localRules];
emit["LOCAL_BRANCH_FREQUENCIES", frequencyRoots];
branchTypes = Piecewise[{
  {Map[Function[root, Piecewise[{{"growing", Im[root] > 0},
       {"decaying", Im[root] < 0}}, "absent"]], frequencyRoots], relativeNorm < 0},
  {"absent", relativeNorm == 0 || Det[metric2] <= 0 || rhoBr[r] == 0},
  {"unable to traverse in a required direction", branchGateLocal && Not[subcritical]}},
  "real propagating"];
emit["LOCAL_BRANCH_TYPES", branchTypes];
emit["LOCAL_HAMILTON_VELOCITY", <|"Velocity" -> groupVelocity,
  "RelativeMetricNorm" -> relativeNorm|>];
emit["LOCAL_FERMAT_OBJECT", <|"TimePolynomial" -> timePolynomial,
  "TimeRoot" -> timeRoot, "EvenMetric" -> opticalMetric, "Odd" -> fermatOdd|>];

(* The actual odd Fermat one-form, before any radial reduction. The compact
   K11 evaluator mechanically extracts through the following boundary. *)
oddOneFormLocal = D[(1 - returnDirection) fermatOdd/2, #] & /@ tangent;
oddOneForm = (atOne[box[# /. localRules]] &) /@ oddOneFormLocal;
oddExterior = Table[D[oddOneForm[[j]], {r,phi}[[i]]] -
  D[oddOneForm[[i]], {r,phi}[[j]]], {i,2},{j,2}];
emit["A_NONRECIPROCAL_PATH_DEPENDENCE", <|"OneForm" -> oddOneForm,
  "ExteriorDerivative" -> Simplify[oddExterior],
  "Closedness" -> Simplify[And @@ Thread[Flatten[oddExterior] == 0]],
  "Domain" -> (r > 0 && atOne[(subcritical && branchGateLocal) /. localRules]),
  "RadialPrimitive" -> Inactive[Integrate][oddOneForm[[1]], {r,rAnchor,rr}]|>];
(* COMPACT_K11_BOUNDARY *)

(* Functional calculus. Integrals remain general radial functionals; their
   parameter derivatives are computed by the Leibniz rule. All derivative
   construction passes through rayDerivative, including K4a/b. *)
lim[z_, x_, target_, opts___] := Block[
  {$Assumptions = And @@ Select[List @@ $Assumptions, FreeQ[#,x] &]},
  Limit[z, x -> target, opts]];
fi[z_, range_List] := If[SameQ[z, 0], 0, Inactive[Integrate][z, range]];
fd[z_, x_] := Module[{ints, slots, algebra, result, v, lo, hi, q},
  ints = DeleteDuplicates[Cases[z, HoldPattern[Inactive[Integrate][_, {_, _, _}]], {0, Infinity}]];
  slots = Table[Unique["integralSlot"], {Length[ints]}];
  algebra = z /. Thread[ints -> slots];
  result = rayDerivative[algebra, x];
  Do[{v, lo, hi} = ints[[i, 2]]; q = ints[[i, 1]];
    result += D[algebra, slots[[i]]] (fi[rayDerivative[q, x], {v, lo, hi}] +
      (q /. v -> hi) rayDerivative[hi, x] - (q /. v -> lo) rayDerivative[lo, x]),
    {i, Length[ints]}];
  result /. Thread[slots -> ints]];
fd[z_, x_, j_Integer] := Nest[fd[#, x] &, z, j];
linearRules[polynomial_, variable_] := {Thread[{variable} ->
  LinearSolve[{{Coefficient[Expand[polynomial], variable]}},
    {-(polynomial /. variable -> 0)}]]};

(* Central optical-coordinate construction from the derived Fermat metric. *)
optical = Simplify[(c0^2 opticalMetric) /. localRules, $Assumptions];
circumferenceRadius = Sqrt[optical[[2, 2]]];
radiusShift = box[r (Sqrt[box[optical[[2, 2]]/r^2]] - 1)];
radialLengthDensity = box[Sqrt[box[optical[[1, 1]]]]];
opticalDensity = radialLengthDensity;
Do[opticalDensity += (-1)^j/Factorial[j] rayDerivative[
  project[radialLengthDensity project[radiusShift^j]], r, j], {j, 1, 3}];
opticalDensity = project[opticalDensity];
hCoefficient[g_] := hCoefficient[g] = Simplify[coeff[opticalDensity, g]];
emit["LOCAL_OPTICAL_COORDINATE", <|"Radius" -> circumferenceRadius,
  "Shift" -> radiusShift, "PushedDensity" -> opticalDensity|>];

(* The physical radii traversed by a central ray are constructed as the
   inverse images of its monotone optical-radius intervals. Turning radii
   are roots of the computed optical radius, not an undefined domain head. *)
endpointRadii = {Sqrt[b^2 + ZE^2], Sqrt[b^2 + ZR^2]};
rayRadius = atOne[circumferenceRadius];
angularSeparation = Total[ArcCos[b/#] & /@ endpointRadii];
exactRadialActionIntegrand = Sqrt[optical[[1,1]]] Sqrt[1-impact^2/optical[[2,2]]];
radarAngleIntegrand = atOne[-D[exactRadialActionIntegrand,impact] /. impact -> impactRadar];
radarEndpointEquation = Total[fi[radarAngleIntegrand,{r,rTurnRadar,#}] & /@ endpointRadii] == angularSeparation;
turningCondition[j_, t_] := (rayRadius /. r -> t) == j && t > bFar &&
  (D[rayRadius, r] /. r -> t) > 0;
flybyRadii = r >= rTurnFlyby;
radarRadii = rTurnRadar <= r <= Max[endpointRadii];
rayConstruction = radarEndpointEquation && turningCondition[b, rTurnFlyby] &&
  turningCondition[impactRadar, rTurnRadar] && rTurnRadar <= Min[endpointRadii] &&
  Inactive[ForAll][r, r >= rTurnFlyby, D[rayRadius, r] > 0] &&
  Inactive[ForAll][r, radarRadii, D[rayRadius, r] > 0];
traversedRadii = flybyRadii || radarRadii;
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
  physicalPropagation && atOne[speed] > 0 && rhoBr[r] != 0 &&
  Element[{deltaProfile[r],velocityProfile[r],xiProfile[r],rhoBr[r]},Reals]];
traversalGate = rayConstruction && Inactive[ForAll][r, traversedRadii, physicalTraversal];
rayGate = propagationGate && traversalGate;
rayDomain = <|"FlybyRadii" -> flybyRadii, "RadarRadii" -> radarRadii,
  "TurningAndExterior" -> rayConstruction, "Endpoints" -> endpointRadii,
  "FlatConnectedBranch" -> Thread[markers -> 0], "FarZone" -> b > bFar|>;
emit["LOCAL_RAY_DOMAIN", rayDomain];
gated[z_, rules_:{}] := <|"OnDomain" -> ConditionalExpression[z, rayGate /. rules],
  "OutsideDomain" -> ConditionalExpression["NOT_ESTABLISHED", Not[rayGate /. rules]],
  "Supplied" -> {"DISPERSION", "ADVECTION", "OPTICAL_RATIO", "EMBEDDING", "LAB_HELD"}|>;
emitGates[prefix_, rules_, wrapper_] := (
  emit[prefix <> "BRANCH_EXISTENCE", wrapper[propagationGate /. rules]];
  emit[prefix <> "PATH_TRAVERSAL", wrapper[traversalGate /. rules]];
  emit[prefix <> "BRANCH_TYPE", wrapper[<|"LocalClassification" ->
    (branchTypes /. localRules /. Thread[markers -> 1] /. rules),
    "Radii" -> (traversedRadii /. rules),
    "TraversalClassification" -> Piecewise[{{"unable to traverse in a required direction",
      Not[traversalGate /. rules]}}, "traversable"]|>]]);
emitGates["", {}, Identity];

(* Universal geometry kernels, obtained by coordinate transformation. *)
angleMap = ArcCos[b/r];
thetaKernel = FullSimplify[2 D[angleMap, r], r > b > 0];
radialKernel = Sqrt[1 - impact^2/rr^2];
angularMap = impact Sec[psi];
angularActionKernel = FullSimplify[(radialKernel /. rr -> angularMap)
  D[angularMap, psi], impact > 0 && 0 < psi < Pi/2];
flatPrimitive = Integrate[radialKernel, rr, Assumptions -> rr > impact > 0,
  GenerateConditions -> False];
flatKernel = Simplify[flatPrimitive - lim[flatPrimitive, rr, impact, Direction -> "FromAbove"], impact > 0 && rr > impact];
bending = Association[];
Do[AssociateTo[bending, label[g] -> If[g == {0,0,0},
  Simplify[2 Integrate[1, {psi,0,Pi/2}] - Pi],
  fi[hCoefficient[g] thetaKernel, {r,b,Infinity}]]], {g,grades}];

(* Expand physical endpoints and solve the finite-endpoint variational
   problem in the complete retained ring, including path displacement. *)
endpointShift = radiusShift /. r -> rr;
endpointDensity = opticalDensity /. r -> rr;
upperCorrection = 0;
Do[upperCorrection += project[project[endpointShift^j] rayDerivative[
  endpointDensity radialKernel, rr, j-1]]/Factorial[j], {j,1,3}];
upperCorrection = project[upperCorrection];
wEnd[g_] := fi[(hCoefficient[g] /. r -> integrationRadius)
  (radialKernel /. rr -> integrationRadius), {integrationRadius,impact,rr}] +
  coeff[upperCorrection,g];
wEndAngular[g_] := fi[(hCoefficient[g] /. r -> angularMap) angularActionKernel,
  {psi,0,ArcCos[impact/rr]}] + coeff[upperCorrection,g];
flatAction = Total[(flatKernel /. rr -> #) & /@ endpointRadii];
angularSeparation = Total[ArcCos[b/#] & /@ endpointRadii];
flatPrincipal = flatAction + impact angularSeparation;
flatDerivatives = Table[FullSimplify[D[flatPrincipal,{impact,j}] /. impact -> b,
  $Assumptions], {j,0,3}];
wPoly[j_] := Total[(wJet[label[#],j] mon[#]) & /@ nonzeroEven];
shift = 0;
Do[gradient = project[sHessian shift + sThird shift^2/2 + wPoly[1] +
  wPoly[2] shift + wPoly[3] shift^2/2];
  shift = project[shift - gradient/sHessian], {iteration,1,2}];
shift = Total[(coeff[shift,#] mon[#]) & /@
  Select[evenGrades, #.{1,1/2,1} <= 2 &]];
stationaryAction = project[wPoly[0] + wPoly[1] shift + wPoly[2] shift^2/2 +
  sHessian shift^2/2 + sThird shift^3/6];
jetRules = Flatten[Table[wJet[label[g],j] -> Total[
  ((If[j == 0, wEnd[g], fd[wEndAngular[g],impact,j]]) /.
    {impact -> b, rr -> #}) & /@ endpointRadii], {g,nonzeroEven},{j,0,2}]];
flatJetRules = {sHessian -> flatDerivatives[[3]], sThird -> flatDerivatives[[4]]};
evenTimes = Association[Table[label[g] ->
  ((coeff[stationaryAction,g] /. flatJetRules) /. jetRules)/c0, {g,grades}]];
oddRadial = box[Coefficient[fermatOdd,dr] /. localRules];
oddTimes = Association[Table[label[g] -> fi[coeff[oddRadial,g] /. r -> integrationRadius,
  {integrationRadius,endpointRadii[[1]],endpointRadii[[2]]}], {g,grades}]];
oneWayER = AssociationMap[evenTimes[#] + oddTimes[#] &, Keys[evenTimes]];
oneWayRE = AssociationMap[evenTimes[#] + returnDirection oddTimes[#] &, Keys[evenTimes]];
roundTrip = AssociationMap[oneWayER[#] + oneWayRE[#] &, Keys[evenTimes]];
nonreciprocal = AssociationMap[(oneWayER[#] - oneWayRE[#])/2 &, Keys[evenTimes]];
emit["LOCAL_ENDPOINT_ACTION", <|"GradeActions" -> (wEnd /@ nonzeroEven),
  "EndpointCorrection" -> upperCorrection, "StationaryAction" -> stationaryAction,
  "JetRules" -> jetRules, "FlatJetRules" -> flatJetRules,
  "RadarImpact" -> (impactRadar == b + (shift /. flatJetRules /. jetRules /. Thread[markers -> 1]))|>];

(* Amendment 4: differentiate the computed finite-endpoint source with ZE,
   ZR fixed, before taking the far-endpoint limit. Integral limits become
   infinity; all remaining endpoint terms are separately computed. The
   limiting endpoint profile values below are independent bounded symbols,
   justified by the supplied optical O(epsilon) counting, not a family. *)
radarSource = roundTrip; (* K6 *)
endpointBounds = Flatten[Table[{
  deltaProfile[endpointRadii[[i]]] -> boundedDelta[i],
  velocityProfile[endpointRadii[[i]]] -> boundedVelocity[i],
  xiProfile'[endpointRadii[[i]]] -> boundedSlope[i]}, {i,2}]];
radarFinite = <||>; radarBoundary = <||>; radarBoundaryLimit = <||>; radarKernels = <||>;
Do[finiteSlope = Expand[-b fd[radarSource[label[g]], b]/2];
  integrals = DeleteDuplicates[Cases[finiteSlope,
    HoldPattern[Inactive[Integrate][_,{_,_,_}]], {0,Infinity}]];
  boundary = Simplify[finiteSlope /. Thread[integrals -> 0]];
  kernel = Simplify[(finiteSlope - boundary) /.
    HoldPattern[Inactive[Integrate][q_, {v_,lo_,hi_}]] :> (q /. v -> r)];
  boundaryLimit = FullSimplify[lim[lim[boundary /. endpointBounds, ZE, Infinity], ZR, Infinity]];
  AssociateTo[radarFinite,label[g] -> finiteSlope];
  AssociateTo[radarBoundary,label[g] -> boundary];
  AssociateTo[radarBoundaryLimit,label[g] -> boundaryLimit];
  AssociateTo[radarKernels,label[g] -> kernel], {g,firstGrades}];
radarSlope = AssociationMap[fi[radarKernels[#],{r,b,Infinity}] + radarBoundaryLimit[#] &,
  Keys[radarKernels]];
emit["LOCAL_RADAR_SLOPE_DERIVATION", <|"Source" -> radarSource,
  "FiniteEndpointDerivative" -> radarFinite, "EndpointRemainder" -> radarBoundary,
  "EndpointLimit" -> radarBoundaryLimit, "IntegralKernels" -> radarKernels,
  "BoundedEndpointSymbols" -> endpointBounds|>];

(* Abel inversion on the exterior half-line. The double-integral kernel
   and the reference inverse are computed here, rather than prescribed.
   Conditions are local differential conditions on general functions,
   with the constant GM fixed at an arbitrary exterior anchor radius. *)
abelKernel = Integrate[1/Sqrt[(v-t)(t-u)], {t,u,v},
  Assumptions -> 0 < u < v, GenerateConditions -> False];
referenceSlope = FullSimplify[lim[lim[-b D[referenceRadar /.
  {rE -> endpointRadii[[1]],rR -> endpointRadii[[2]]},b]/2,
  ZE, Infinity], ZR, Infinity], $Assumptions];
abelInverse[ref_] := Module[{transformed, primitive},
  transformed = (ref /. b -> Sqrt[t])/(2 Sqrt[t]);
  primitive = Integrate[transformed/Sqrt[t-u], {t,u,Infinity},
    Assumptions -> u > 0 && c0 > 0 && Element[{GM,gamma},Reals], GenerateConditions -> False];
  FullSimplify[(-2 u D[primitive,u]/abelKernel) /. u -> r^2, r > 0 && c0 > 0]];
quantifierMode = "EVERY_B"; (* K7 *)
comparisonThetaKernel = Total[hCoefficient[#] thetaKernel & /@ firstGrades];
comparisonRadarKernel = Total[radarKernels[label[#]] & /@ firstGrades];
comparisonTheta = fi[comparisonThetaKernel,{r,b,Infinity}];
comparisonRadar = fi[comparisonRadarKernel,{r,b,Infinity}] + Total[Values[radarBoundaryLimit]];
regularityDomain = <|"Radius" -> r > bFar, "Anchor" -> rAnchor > bFar,
  "Profiles" -> {deltaProfile,velocityProfile,xiProfile,rhoBr,fProfile},
  "FunctionalDomain" -> {"required derivatives exist", "displayed improper integrals converge",
    "Abel inverse exists on the exterior half-line"},
  "OpticalCounting" -> {"delta O(1/r)", "V/c0 O(1/Sqrt[r])", "xiPrimeSquared O(1/r)"}|>;
reduceProfile[kernel_, reference_, rules_, extra_] := Module[{q, power, weight, target, mass, equation, rule},
  q = Simplify[(kernel /. rules)/thetaKernel];
  If[quantifierMode === "FIXED_B", Return[<|"Domain" -> extra && bFixed > bFar,
    "BranchDomain" -> (rayGate /. rules),
    "Cases" -> {<|"Domain" -> extra && bFixed > bFar,
      "GM" -> linearRules[(fi[kernel /. rules,{r,b,Infinity}] - reference) /. b -> bFixed,GM]|>}|>]];
  power = Exponent[Simplify[kernel/thetaKernel],b];
  weight = Simplify[q/b^power]; target = abelInverse[reference/b^power];
  rule = linearRules[weight-target,GM]; mass = GM /. First[rule];
  equation = Simplify[D[mass,r] == 0];
  <|"Domain" -> regularityDomain, "BranchDomain" -> (rayGate /. rules),
    "Cases" -> {<|"Domain" -> extra && r > bFar && rAnchor > bFar,
      "GM" -> (rule /. r -> rAnchor), "ProfileEquation" -> equation,
      "ProfileEquationRadius" -> r, "InverseWeight" -> weight,
      "ReferenceInverse" -> target|>}|>];
emit["LOCAL_ABEL_REDUCTION", <|"DoubleIntegralKernel" -> abelKernel,
  "DeflectionKernel" -> comparisonThetaKernel, "RadarKernel" -> comparisonRadarKernel,
  "ReferenceDeflectionInverse" -> abelInverse[referenceTheta /. gamma -> 1],
  "ReferenceRadarSlope" -> referenceSlope, "Domain" -> regularityDomain|>];

(* Exact linear rank strata include GM=0. No branch of a Piecewise slope
   is discarded; general functionals occupy the observable parameter. *)
gammaTemplate = Reduce[aa gamma + bb == qq, gamma, Reals];
gammaDifferenceTemplate = Reduce[Exists[{gd,gr},
  ad gd + bd == qd && ar gr + br == qr && gammaDifference == gd-gr],
  gammaDifference, Reals];
gammaSolution[obs_, ref_] := gammaTemplate /.
  {aa -> Coefficient[ref,gamma], bb -> (ref /. gamma -> 0), qq -> obs};
gammaDifferenceSolution[obsD_,obsR_] := gammaDifferenceTemplate /.
  {ad -> Coefficient[referenceTheta,gamma], bd -> (referenceTheta /. gamma -> 0), qd -> obsD,
   ar -> Coefficient[referenceSlope,gamma], br -> (referenceSlope /. gamma -> 0), qr -> obsR};
emit["LOCAL_GAMMA_STRATA", <|"Single" -> gammaTemplate, "Difference" -> gammaDifferenceTemplate|>];

(* Coordinate mass law, with live density inside the component derivatives. *)
cartesian = {x1,x2,x3}; cartRadius = Sqrt[cartesian.cartesian];
massDensity = rhoBr[r]; (* K8 *)
massVector = (massDensity velocity cartesian/r) /. r -> cartRadius;
divergence = Simplify[Total[MapThread[D,{massVector,cartesian}]] /.
  {x1 -> r,x2 -> 0,x3 -> 0}, $Assumptions];
exchange = jn[r] /. First[Solve[massBalanceInput /. divRhoV -> divergence, jn[r]]];
emit["LOCAL_MASS_BALANCE", <|"Measure" -> "coordinate d3x", "Components" -> massVector,
  "Divergence" -> divergence, "Exchange" -> exchange|>];
impliedExchange[condition_, rules_] := <|"ReducedProfileCondition" -> condition,
  "ImpliedExchange" -> (jn[r] == (exchange /. rules)),
  "SignedVelocity" -> (velocity /. rules), "Density" -> (massDensity /. rules)|>;

(* Shared vocabulary and all substitutions, including their branch gates. *)
conditionPair[rules_, domain_] := {
  reduceProfile[comparisonThetaKernel,referenceTheta /. gamma -> 1,rules,domain],
  reduceProfile[comparisonRadarKernel,referenceSlope /. gamma -> 1,rules,domain]};
observations = {"DEFLECTION","RADAR"};
emitComparisons[prefix_, rules_, domain_, wrapper_, withJn_] := Module[{obs, refs, cond},
  obs = {comparisonTheta,comparisonRadar} /. rules;
  refs = {referenceTheta,referenceSlope}; cond = conditionPair[rules,domain];
  Do[emit[prefix <> "B_GAMMA_" <> observations[[i]], wrapper[gated[gammaSolution[obs[[i]],refs[[i]]],rules]]];
    emit[prefix <> "B_RESIDUAL_" <> observations[[i]], wrapper[gated[obs[[i]]-(refs[[i]] /. gamma -> 1),rules]]];
    emit[prefix <> "B_CONDITION_" <> observations[[i]], wrapper[cond[[i]]]];
    If[withJn,emit[prefix <> "B_IMPLIED_JN_" <> observations[[i]], wrapper[impliedExchange[cond[[i]],rules]]]], {i,2}];
  emit[prefix <> "B_GAMMA_DIFFERENCE", wrapper[gated[gammaDifferenceSolution @@ obs,rules]]]];
allA = <|"DEFLECTION" -> bending, "ROUND_TRIP" -> roundTrip,
  "ONE_WAY_ER" -> oneWayER, "ONE_WAY_RE" -> oneWayRE, "NONRECIPROCAL" -> nonreciprocal|>;
emitA[prefix_,rules_,wrapper_,which_] := (
  Do[emit[prefix <> "A_" <> name <> "_G" <> label[g],
    wrapper[gated[allA[name][label[g]] /. rules,rules]]], {name,which},{g,grades}];
  Do[emit[prefix <> "A_RADAR_LOG_G" <> label[g],
    wrapper[gated[radarSlope[label[g]] /. rules,rules]]], {g,firstGrades}]);
emitA["",{},Identity,Keys[allA]];
emitComparisons["",{},True,Identity,True];
responseDeltas = (Simplify[Normal[Series[#/c0 - 1,{fLocal,0,1}]],
  rho0 > 0 && c0 > 0 && Element[{n,s,fLocal},Reals]] &) /@ responseInputs;
responseRules = Table[{deltaProfile -> Function[{r}, Evaluate[responseDeltas[[i]] /. fLocal -> fProfile[r]]]}, {i,3}];
responseNames = {"CONSTANT","FIXED_RATIO","POWER"};
responseDomain[i_] := With[{resp = responseInputs[[i]],change = responseDeltas[[i]]},
  Element[resp,Reals] && resp > 0 && Abs[change] < 1 &&
    If[i == 2, localSoundSquared > 0 && (bulkSoundSquared /. rho -> rho0) > 0,
      If[i == 3,localBulkDensity/rho0 > 0,True]]] /. fLocal -> fProfile[r];
emit["LOCAL_C_RESPONSES", <|"Responses" -> responseInputs,"DeltaSeries" -> responseDeltas,
  "Domains" -> Table[responseDomain[i],{i,3}]|>];
emitC[prefix_, additional_, wrapper_, withJn_, indices_] := Do[
  cRules = Join[responseRules[[i]],additional];
  cConditions = conditionPair[cRules,responseDomain[i]];
  Do[emit[prefix <> "C_" <> responseNames[[i]] <> "_CONDITION_" <> observations[[j]],wrapper[cConditions[[j]]]];
    If[withJn, emit[prefix <> "C_" <> responseNames[[i]] <> "_IMPLIED_JN_" <> observations[[j]],
      wrapper[impliedExchange[cConditions[[j]],cRules]]]], {j,2}];
  If[prefix == "", emit["C_" <> responseNames[[i]] <> "_N_DEPENDENCE", <|
    "ResponseDerivative" -> D[responseDeltas[[i]],n],
    "KernelDerivatives" -> D[{comparisonThetaKernel,comparisonRadarKernel} /. cRules,n],
    "ConditionGMderivatives" -> (D[(GM /. First[#["Cases"][[1]]["GM"]]),n] & /@ cConditions)|>]], {i,indices}];
emitC["",{},Identity,True,Range[3]];
zeroDelta = {deltaProfile -> Function[{r},0]};
zeroVelocity = {velocityProfile -> Function[{r},0]};
zeroXi = {xiProfile -> Function[{r},0]};
restrictions = <|"FLOW_ONLY" -> Join[zeroDelta,zeroXi],
  "SPEED_ONLY" -> Join[zeroVelocity,zeroXi], "TILT_ONLY" -> Join[zeroDelta,zeroVelocity]|>;
KeyValueMap[Function[{name,rules},
  emit["LOCAL_RESTRICTION_" <> name, rules];
  emitA["R_" <> name <> "_",rules,Identity,{"DEFLECTION"}];
  emitComparisons["R_" <> name <> "_",rules,True,Identity,True]],restrictions];
emitC["R_BULK_ONLY_",Join[zeroVelocity,zeroXi],Identity,True,{2,3}];

(* The forward case solves the general mass law, without selecting rhoBr. *)
forwardMassVector = (massDensity vForward[r] cartesian/r) /. r -> cartRadius;
forwardDiv = Simplify[Total[MapThread[D,{forwardMassVector,cartesian}]] /.
  {x1 -> r,x2 -> 0,x3 -> 0},$Assumptions];
forwardSolution = DSolve[forwardDiv == 0,vForward,r];
fluxInput = Phi == Integrate[massDensity vForward[r] r^2 Sin[theta],
  {theta,0,Pi},{phi,0,2 Pi}];
forwardV = vForward[r] /. First[Solve[fluxInput,vForward[r]]];
forwardRules = {velocityProfile -> Function[{r},Evaluate[forwardV]]};
forwardPremise = <|"NormalExchange" -> (jn[r] == 0),"Flux" -> fluxInput,
  "Measure" -> "coordinate d3x"|>;
forwardWrap[z_] := <|"Premise" -> forwardPremise,"Object" -> z|>;
emit["F_MASS_SOLUTION",forwardWrap[<|"GeneralSolution" -> forwardSolution,
  "FluxVelocity" -> forwardV,"Residual" -> Simplify[forwardDiv /.
    vForward -> Function[{r},Evaluate[forwardV]]]|>]];
Do[stageRules = Join[forwardRules,If[stage == 1,Join[zeroDelta,zeroXi],{}]];
  prefix = If[stage == 1,"F_FLOW_","F_LIVE_"];
  emitGates[prefix,stageRules,forwardWrap];
  emitA[prefix,stageRules,forwardWrap,Keys[allA]];
  emitComparisons[prefix,stageRules,massDensity != 0,forwardWrap,False],{stage,2}];
emitC["F_LIVE_",forwardRules,forwardWrap,False,Range[3]];
emit["LOCAL_ENGINE_LOCAL_NAMES", Append[localNames, "WL_LOCAL_S9B_ENGINE_LOCAL_NAMES"]];
