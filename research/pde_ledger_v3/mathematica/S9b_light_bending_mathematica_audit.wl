(* S9b blind builder.  Sole physical input: S9b_SHARED_PHYSICS.md v10, Parts A-C only.
   No Get, Needs, Import, file access, exports, or runtime configuration.
   Power-law ansatz, with symbolic amplitudes AND exponents.  All expressions
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
    ampD, ampV, ampW, powD, powV, powW, GM, n, K, m, s, rho0,
    ampF}, Reals] && c0 > 0 && ell > 0 && b > 0 && bFar > 0 &&
    ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 &&
    powD >= 1 && powV >= 1/2 && powW >= 1/2 && rho0 > 0 && m > 0;

(* Supplied equations and chosen ansatz: the only hand-combined physics.
   powW is the exponent of xi', so the embedding grade has exponent 2 powW.
   rhoBr is an arbitrary live radial function, not the bulk number density.
*)
delta = ampD (ell/r)^powD;
velocity = c0 ampV (ell/r)^powV;
xiAnsatz = ampW ell^powW Piecewise[{{Log[r], powW == 1}}, r^(1 - powW)/(1 - powW)];
xiSlope = Simplify[D[xiAnsatz, r], $Assumptions]; (* K4c *)
densityFraction = ampF (ell/r)^powD;
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

emit["SUPPLIED_INPUTS", <|"Dispersion" -> dispersion == 0,
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
emit["PROFILE_ANSATZ", <|"Delta" -> delta, "V" -> velocity,
  "XiDerivative" -> xiSlope, "RhoBr" -> rhoBr[r],
  "FractionalBulkDensity" -> densityFraction, "Domain" -> $Assumptions,
  "OrderCountingDomain" -> (powD >= 1 && powV >= 1/2 && powW >= 1/2),
  "Anchoring" -> "LAB_HELD"|>];
emit["RETAINED_GRADES", <|"Markers" -> markers, "Grades" -> grades,
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
emit["BRANCH_FREQUENCIES", frequencyRoots];
emit["BRANCH_EXISTENCE", <|"Local" -> propagatingLocal,
  "AlongRay" -> Inactive[ForAll][r, rayRadiusDomain[r],
    physicalPropagation && atOne[speed] > 0 && rhoBr[r] != 0],
  "SpeedIdentification" -> localSpeedIdentification|>];
emit["PATH_TRAVERSAL", <|"RadialGroupVelocities" -> radialSpeeds,
  "DirectionalCondition" -> traversalLocal,
  "SubcriticalDomain" -> Inactive[ForAll][r, rayRadiusDomain[r], physicalTraversal]|>];
branchTypes = Piecewise[{
  {Map[Function[root, Piecewise[{{"growing", Im[root] > 0},
       {"decaying", Im[root] < 0}}, "absent"]], frequencyRoots], relativeNorm < 0},
  {"absent", relativeNorm == 0 || Det[metric2] <= 0},
  {"unable to traverse in a required direction", branchGateLocal && Not[subcritical]}},
  "real propagating"];
emit["BRANCH_TYPES", branchTypes];
emit["LOCAL_HAMILTON_VELOCITY", <|"Velocity" -> groupVelocity,
  "RelativeMetricNorm" -> relativeNorm|>];
emit["LOCAL_FERMAT_OBJECT", <|"TimePolynomial" -> timePolynomial,
  "TimeRoot" -> timeRoot, "EvenMetric" -> opticalMetric, "Odd" -> fermatOdd|>];

(* The actual odd Fermat one-form, before any radial reduction. The compact
   K11 evaluator mechanically extracts through the following boundary. *)
oddOneFormLocal = D[fermatOdd, #] & /@ tangent;
oddExteriorLocal = Table[D[oddOneFormLocal[[j]], {r, phi}[[i]]] -
  D[oddOneFormLocal[[i]], {r, phi}[[j]]], {i, 2}, {j, 2}];
(* Local symbols stand for radial fields: include their chain-rule jets. *)
oddExteriorLocal = oddExteriorLocal + Table[If[i == 1,
  Total[MapThread[D[oddOneFormLocal[[j]], #1] #2 &,
    {{aMetric, cSquared, vLocal}, {aJet, cJet, vJet}}]], 0] - If[j == 1,
  Total[MapThread[D[oddOneFormLocal[[i]], #1] #2 &,
    {{aMetric, cSquared, vLocal}, {aJet, cJet, vJet}}]], 0], {i, 2}, {j, 2}];
emit["NONRECIPROCAL_DEPENDENCE", <|"OneForm" -> oddOneFormLocal,
  "ExteriorDerivative" -> Simplify[oddExteriorLocal],
  "Closedness" -> Simplify[And @@ Thread[Flatten[oddExteriorLocal] == 0]],
  "Domain" -> (r > 0 && subcritical && branchGateLocal),
  "RadialPrimitive" -> Inactive[Integrate][oddOneFormLocal[[1]], r]|>];
(* COMPACT_K11_BOUNDARY *)

(* Normalize the even optical metric to far-field length units.  R(r) is
   its circumferential radius.  The push-forward of a(r) dr under R=r+u(r)
   has formal density Sum[(-D_R)^j (a u^j)/j!].  Three iterations suffice:
   every even perturbation contains xD, xV^2 or xW; no fourth product
   survives the requested box.  This includes the shift of the periapsis.
*)
optical = Simplify[(c0^2 opticalMetric) /. localRules, $Assumptions];
circumferenceRadius = Simplify[Sqrt[optical[[2, 2]]], $Assumptions];
(* Derive the same coordinate from the computed angular metric, keeping
   the analytic branch connected to positive r at zero perturbation. *)
radiusShift = box[r (Sqrt[box[optical[[2, 2]]/r^2]] - 1)];
radialLengthDensity = box[Sqrt[box[optical[[1, 1]]]]];
opticalDensity = radialLengthDensity;
Do[opticalDensity += (-1)^j/Factorial[j] rayDerivative[
    project[radialLengthDensity project[radiusShift^j]], r, j], {j, 1, 3}];
opticalDensity = project[Expand[opticalDensity]];
exponent[g_] := g.{powD, powV, 2 powW};
hCoefficient[g_] := hCoefficient[g] = Simplify[coeff[opticalDensity, g] /. r -> ell,
  $Assumptions];
emit["LOCAL_OPTICAL_COORDINATE", <|"Radius" -> circumferenceRadius,
  "Shift" -> radiusShift, "RadialDensity" -> radialLengthDensity,
  "PushedDensity" -> opticalDensity|>];

(* Conditions include a simple exterior turning point and the connected
   perturbative branch of the boundary-value solution.  rayRadiusDomain
   is the union of radii traversed by that branch, not the throat interior.
*)
rayDomain = <|"LocalPropagation" -> physicalPropagation,
  "PositiveSpeed" -> atOne[speed] > 0,
  "Traversal" -> physicalTraversal,
  "TurningPoint" -> {atOne[circumferenceRadius] == impact,
    D[atOne[circumferenceRadius], r] > 0},
  "Exterior" -> Inactive[ForAll][r, r > rTurning,
    atOne[circumferenceRadius] > impact],
  "PerturbativeBranch" -> AnalyticContinuationFrom[Thread[markers -> 0]],
  "Endpoints" -> {rE == Sqrt[b^2 + ZE^2], rR == Sqrt[b^2 + ZR^2]},
  "FarZone" -> b > bFar|>;
emit["RAY_DOMAIN", rayDomain];

(* Universal quadratures are computed once using nonphysical dummy
   variables.  Their parameters are subsequently replaced by computed
   grade exponents.  The finite-endpoint radial action is regular at its
   turning point. *)
angularKernel = Integrate[Cos[psi]^nu, {psi, 0, Pi/2},
  Assumptions -> nu > 0, GenerateConditions -> False];
(* Convert the elementary integrand to its incomplete-beta integral by
   solving for the two exponents.  This representation stays regular at
   all positive nu, including nu=2; generic antiderivatives can introduce
   removable poles there.  The defining derivative is emitted below. *)
betaIntegrand = y^(1/2) (1 - y)^((nu - 3)/2);
betaParameters = First[SolveAlways[Together[D[betaIntegrand, y]/betaIntegrand -
  ((alphaBeta - 1)/y - (betaBeta - 1)/(1 - y))] == 0, y]];
betaNormalization = Simplify[betaIntegrand/
  (y^(alphaBeta - 1) (1 - y)^(betaBeta - 1)) /. betaParameters,
  0 < y < 1 && Element[nu, Reals]];
betaPrimitive = betaNormalization (Beta[xx, alphaBeta, betaBeta] /. betaParameters);
betaDerivative = Simplify[D[betaPrimitive, xx] - (betaIntegrand /. y -> xx),
  0 < xx < 1 && Element[nu, Reals]];
(* R=impact/Sqrt[1-y] transforms R^-nu Sqrt[1-impact^2/R^2] dR. *)
radialMap = impact/Sqrt[1 - y];
radialJacobian = Simplify[radialMap^(-nu) Sqrt[1 - impact^2/radialMap^2]
  D[radialMap, y], impact > 0 && 0 < y < 1 && Element[nu, Reals]];
kernelFactor = Simplify[radialJacobian/(y^(1/2) (1 - y)^((nu - 3)/2)),
  impact > 0 && 0 < y < 1 && Element[nu, Reals]];
radialKernel = kernelFactor (betaPrimitive /. xx -> 1 - impact^2/rr^2);
flatPrimitive = Integrate[Sqrt[1 - impact^2/rr^2], rr,
  Assumptions -> rr > impact > 0, GenerateConditions -> False];
flatKernel = Simplify[flatPrimitive - Block[{$Assumptions = impact > 0},
  Limit[flatPrimitive, rr -> impact, Direction -> "FromAbove"]], rr > impact > 0];
emit["LOCAL_QUADRATURE_KERNELS", <|"Angular" -> angularKernel,
  "RadialChangeOfVariable" -> radialJacobian, "Radial" -> radialKernel,
  "BetaIntegrand" -> betaIntegrand, "BetaParameters" -> betaParameters,
  "BetaDerivativeResidual" -> betaDerivative,
  "FlatRadial" -> flatKernel|>];

bending = Association[];
Do[AssociateTo[bending, label[g] -> If[g == {0, 0, 0},
    Simplify[2 Integrate[1, {psi, 0, Pi/2}] - Pi],
    Simplify[2 hCoefficient[g] (ell/b)^exponent[g]
      (angularKernel /. nu -> exponent[g]), $Assumptions]]], {g, grades}];

(* Fixed physical endpoints.  Expand the upper limit R(r_i) as well as
   the density.  The action W(impact)+impact Phi is stationary at the
   angular momentum selected by the endpoints.  Universal stationarity
   is solved in the retained ring below; this keeps all path corrections.
*)
endpointShift = radiusShift /. r -> rr;
endpointDensity = Total[(hCoefficient[#] mon[#] (ell/rr)^exponent[#]) & /@ evenGrades];
upperCorrection = 0;
Do[upperCorrection += project[project[endpointShift^j] rayDerivative[
    endpointDensity Sqrt[1 - impact^2/rr^2], rr, j - 1]]/Factorial[j], {j, 1, 3}];
upperCorrection = project[upperCorrection];
wEnd[g_List] := wEnd[g] = hCoefficient[g] ell^exponent[g] *
    (radialKernel /. nu -> exponent[g]) + coeff[upperCorrection, g];
endpointRadii = {Sqrt[b^2 + ZE^2], Sqrt[b^2 + ZR^2]};
flatAction = Total[(flatKernel /. rr -> #) & /@ endpointRadii];
angularSeparation = Total[(ArcCos[b/#]) & /@ endpointRadii];
flatPrincipal = flatAction + impact angularSeparation;
flatDerivatives = Table[FullSimplify[D[flatPrincipal, {impact, j}] /. impact -> b,
  $Assumptions], {j, 0, 3}];

(* Symbols wJet[grade,j] stand for the displayed derivatives of W.
   Solve the formal stationarity equation with a symbolic flat Hessian.
   The equations, solution, and full substitution map are emitted. *)
wPoly[j_] := Total[(wJet[label[#], j] mon[#]) & /@ nonzeroEven];
shift = 0;
Do[gradient = project[sHessian shift + sThird shift^2/2 +
    wPoly[1] + wPoly[2] shift + wPoly[3] shift^2/2];
  shift = project[shift - gradient/sHessian], {iteration, 1, 2}];
shift = Total[(coeff[shift, #] mon[#]) & /@
  Select[evenGrades, #.{1, 1/2, 1} <= 2 &]];
stationaryAction = project[wPoly[0] + wPoly[1] shift + wPoly[2] shift^2/2 +
  sHessian shift^2/2 + sThird shift^3/6];
jetRules = Flatten[Table[wJet[label[g], j] ->
    Total[(Simplify[D[wEnd[g], {impact, j}] /. {impact -> b, rr -> #},
      $Assumptions]) & /@ endpointRadii], {g, nonzeroEven}, {j, 0, 2}]];
flatJetRules = {sHessian -> flatDerivatives[[3]], sThird -> flatDerivatives[[4]]};
emit["LOCAL_ENDPOINT_ACTION", <|"UpperCorrection" -> upperCorrection,
  "GradeActions" -> Table[{g, wEnd[g]}, {g, nonzeroEven}],
  "FlatAction" -> flatPrincipal, "FlatDerivatives" -> flatDerivatives|>];
emit["LOCAL_STATIONARY_ACTION", <|"ShiftThroughEvenDegreeTwo" -> shift,
  "Action" -> stationaryAction, "JetSubstitution" -> jetRules,
  "FlatJetSubstitution" -> flatJetRules|>];

evenTimes = Association[];
Do[AssociateTo[evenTimes, label[g] ->
  ((coeff[stationaryAction, g] /. flatJetRules) /. jetRules)/c0], {g, grades}];

(* The odd Fermat one-form is integrated on the actual endpoints.  Its
   exterior derivative tests path dependence without choosing a ray. *)
oddRadial = box[(Coefficient[fermatOdd, dr] /. localRules)];
oddComponents = {oddRadial, 0, 0};
coords = {r, theta, phi};
oddExterior = Table[D[oddComponents[[j]], coords[[i]]] -
  D[oddComponents[[i]], coords[[j]]], {i, 3}, {j, 3}];
primitiveGeneral = Integrate[rr^(-nu), rr, GenerateConditions -> False];
primitiveResonance = Integrate[rr^(-nu) /. nu -> 1, rr];
oddTimes = Association[];
Do[oddPower = exponent[g]; oddAmp = Simplify[coeff[oddRadial, g] /. r -> ell];
  oddPrimitive = oddAmp ell^oddPower Piecewise[{{primitiveResonance, oddPower == 1}},
    primitiveGeneral /. nu -> oddPower];
  AssociateTo[oddTimes, label[g] -> ((oddPrimitive /. rr -> endpointRadii[[2]]) -
    (oddPrimitive /. rr -> endpointRadii[[1]]))], {g, grades}];
(* All observables are expressions on the computed branch, rather than
   values evaluated on a branch which cannot propagate. *)
rayGate = Inactive[ForAll][r, rayRadiusDomain[r],
  physicalPropagation && physicalTraversal && atOne[speed] > 0 && rhoBr[r] != 0];
gated[z_, gate_:rayGate, rules_:{}] := <|"OnDomain" -> ConditionalExpression[z, gate],
  "OutsideDomain" -> ConditionalExpression["NOT_ESTABLISHED", Not[gate]],
  "BranchType" -> (branchTypes /. localRules /. Thread[markers -> 1] /. rules)|>;
returnDirection = -1; (* K5 *)
oneWayER = AssociationMap[evenTimes[#] + oddTimes[#] &, Keys[evenTimes]];
oneWayRE = AssociationMap[evenTimes[#] + returnDirection oddTimes[#] &, Keys[evenTimes]];
roundTrip = AssociationMap[oneWayER[#] + oneWayRE[#] &, Keys[evenTimes]];
nonreciprocal = AssociationMap[(oneWayER[#] - oneWayRE[#])/2 &, Keys[evenTimes]];
Do[emit["A_DEFLECTION_G" <> label[g], gated[bending[label[g]]]];
  emit["A_ROUND_TRIP_G" <> label[g], gated[roundTrip[label[g]]]];
  emit["A_ONE_WAY_ER_G" <> label[g], gated[oneWayER[label[g]]]];
  emit["A_ONE_WAY_RE_G" <> label[g], gated[oneWayRE[label[g]]]];
  emit["A_NONRECIPROCAL_G" <> label[g], gated[nonreciprocal[label[g]]]], {g, grades}];

(* Extract the logarithm from the computed time, not from a separately
   assigned density coefficient. Its incomplete-beta primitive is expanded
   on the exponent-one stratum. The endpoint corrections are algebraic.
   Positive exponents exclude further logarithmic resonances. *)
resonantMap = impact Cosh[hyperbolicParameter];
resonantIntegrand = FullSimplify[(Sqrt[1 - impact^2/resonantMap^2]/resonantMap)
  D[resonantMap, hyperbolicParameter], impact > 0 && hyperbolicParameter > 0];
resonantKernel = Integrate[resonantIntegrand,
  {hyperbolicParameter, 0, ArcCosh[rr/impact]},
  Assumptions -> rr > impact > 0, GenerateConditions -> False];
resonantAsymptotic = FullSimplify[Normal[Series[resonantKernel /. rr -> impact/zeta,
  {zeta, 0, 0}]], impact > 0 && 0 < zeta < 1];
logPerEndpoint = FullSimplify[zeta D[resonantAsymptotic, zeta]/(-2),
  impact > 0 && 0 < zeta < 1];
radarSource = roundTrip; (* K6 *)
(* wEnd uses kernelFactor times Beta. This quotient obtains that Beta's
   log coefficient from the very same computed primitive. *)
betaLogCoefficient = Simplify[logPerEndpoint/(kernelFactor /. nu -> 1), impact > 0];
logFromTime[z_] := Module[{betas, answer = 0, q, factor, rest, resonanceDomain},
  betas = DeleteDuplicates[Cases[z, _Beta, Infinity]];
  Do[q = 2 bt[[3]] + 1;
    resonanceDomain = Simplify[q == 1, $Assumptions];
    If[!SameQ[resonanceDomain, False],
      factor = Coefficient[Expand[z], bt];
      answer += Piecewise[{{Simplify[factor betaLogCoefficient,
        $Assumptions && resonanceDomain], resonanceDomain}}, 0]], {bt, betas}];
  Simplify[answer /. impact -> b, $Assumptions]];
radarLog = AssociationMap[logFromTime[radarSource[#]] &, Keys[radarSource]];
emit["LOCAL_RADAR_LOG_EXTRACTION", <|"Source" -> radarSource,
  "ResonantIntegral" -> resonantKernel, "LargeEndpointExpansion" -> resonantAsymptotic,
  "BetaLogCoefficient" -> betaLogCoefficient, "Coefficients" -> radarLog,
  "Domain" -> (ZE/b > 1 && ZR/b > 1)|>];
thetaFirst = Total[bending[label[#]] & /@ firstGrades];
radarFirst = Total[radarLog[label[#]] & /@ firstGrades];
referenceLog = Coefficient[Expand[referenceRadar /. Log[4 rE rR/b^2] -> logBasis], logBasis];
gammaTheta = gamma /. First[Solve[referenceTheta == thetaFirst, gamma]];
gammaRadar = gamma /. First[Solve[referenceLog == radarFirst, gamma]];
thetaResidual = thetaFirst - (referenceTheta /. gamma -> 1);
radarResidual = radarFirst - (referenceLog /. gamma -> 1);

(* Coordinate Cartesian component calculus: no induced-measure object. *)
cartesian = {x1, x2, x3}; cartRadius = Sqrt[cartesian.cartesian];
massDensity = rhoBr[r]; (* K8 *)
massVector = (massDensity velocity cartesian/r) /. r -> cartRadius;
divergenceCartesian = Total[MapThread[D, {massVector, cartesian}]];
divergence = Simplify[divergenceCartesian /. {x1 -> r, x2 -> 0, x3 -> 0}, $Assumptions];
exchange = jn[r] /. First[Solve[massBalanceInput /. divRhoV -> divergence, jn[r]]];
emit["MASS_BALANCE", <|"Measure" -> "coordinate d3x", "Components" -> massVector,
  "Divergence" -> divergence, "Exchange" -> exchange|>];

(* Exhaustive exponent partitions reduce equality of finite power sums on
   an open interval. The exponent-one group solves GM; every other group
   has vanishing total coefficient. No ForAll survives this reduction.
   Equalities among the amplitudes are kept as exact algebraic constraints,
   so zero amplitudes, coincident tails and cancellations are all included. *)
linearRules[polynomial_, variable_] := {Thread[{variable} ->
  LinearSolve[{{Coefficient[Expand[polynomial], variable]}},
    {-(polynomial /. variable -> 0)}]]};
partitions[{}] = {{}};
partitions[list_List] := partitions[list] = Module[{a = First[list], rest},
  rest = partitions[Rest[list]];
  Flatten[Map[Function[p, Join[{Prepend[p, {a}]},
    Table[ReplacePart[p, i -> Prepend[p[[i]], a]], {i, Length[p]}]]], rest], 1]];
comparisonPowers = exponent /@ firstGrades;
thetaCoefficients = Simplify[Table[bending[label[g]] b^exponent[g], {g, firstGrades}], $Assumptions];
referenceThetaCoefficient = Simplify[b (referenceTheta /. gamma -> 1)];
quantifierMode = "EVERY_B"; (* K7 *)
reducePowers[cs_, ps_, ref_, domain_] := Module[{allP, allC, cases, groups, reps,
    equalities, inequalities, equations, gmGroup, gmSolution, other},
  If[quantifierMode === "FIXED_B", Return[<|"Domain" -> domain && bFixed > bFar,
    "Cases" -> linearRules[Total[MapThread[#1 bFixed^(-#2) &, {cs, ps}]] - ref/bFixed, GM]|>]];
  allP = Append[ps, 1]; allC = Append[cs, -ref];
  cases = Table[groups = pp; reps = First /@ groups;
    equalities = And @@ Flatten[Table[Thread[allP[[gg]] == allP[[First[gg]]]], {gg, groups}]];
    inequalities = And @@ (Unequal @@ # & /@ Subsets[allP[[reps]], {2}]);
    gmGroup = First[Select[groups, MemberQ[#, Length[allC]] &]];
    gmSolution = linearRules[Total[allC[[gmGroup]]], GM];
    other = DeleteCases[groups, gmGroup];
    equations = And @@ (Total[allC[[#]]] == 0 & /@ other);
    <|"Domain" -> Simplify[domain && equalities && inequalities, $Assumptions],
      "GM" -> gmSolution, "AmplitudeConditions" -> Simplify[equations, $Assumptions]|>,
    {pp, partitions[Range[Length[allP]]]}];
  <|"Domain" -> domain, "Cases" -> cases|>];
reduceRadar[rad_, domain_] := Module[{pw, cases, pred, value},
  pw = DeleteDuplicates[Flatten[(#[[1, All, 2]] &) /@ Cases[rad, _Piecewise, {0, Infinity}]]];
  cases = Table[pred = And @@ MapThread[If[#2, #1, Not[#1]] &, {pw, mask}];
    value = Simplify[rad /. Thread[pw -> mask], $Assumptions && pred];
    <|"Domain" -> Simplify[domain && pred, $Assumptions],
      "GM" -> linearRules[value - (referenceLog /. gamma -> 1), GM]|>,
    {mask, Tuples[{False, True}, Length[pw]]}];
  <|"Domain" -> domain, "Cases" -> cases|>];
profileDomain = powD >= 1 && powV >= 1/2 && powW >= 1/2 && b > bFar;
conditionPair[rules_, dom_] := (Append[#, "BranchDomain" -> (rayGate /. rules)] & /@ {
  reducePowers[thetaCoefficients /. rules, comparisonPowers /. rules, referenceThetaCoefficient, dom],
  reduceRadar[radarFirst /. rules, dom]});
conditions = conditionPair[{}, profileDomain];
constrainedExchange[cond_, rules_:{}] := <|"Domain" -> cond["Domain"],
  "BranchDomain" -> cond["BranchDomain"],
  "Cases" -> Map[Append[#, "ImpliedExchange" -> <|
    "Relation" -> (jn[r] == (exchange /. rules)),
    "Density" -> (massDensity /. rules), "SignedVelocity" -> (velocity /. rules)|>] &,
    cond["Cases"]]|>;
emit["B_FIRST_ORDER_DEFLECTION", gated[Table[{g, bending[label[g]], g.{1, 1/2, 1}}, {g, firstGrades}]]];
emit["B_FIRST_ORDER_RADAR_LOG", gated[Table[{g, radarLog[label[g]], g.{1, 1/2, 1}}, {g, firstGrades}]]];
emit["B_GAMMA_DEFLECTION", gated[gammaTheta]];
emit["B_GAMMA_RADAR", gated[gammaRadar]];
emit["B_GAMMA_DIFFERENCE", gated[gammaTheta - gammaRadar]];
emit["B_DEFLECTION_RESIDUAL", gated[<|"Computed" -> thetaFirst,
  "Reference" -> (referenceTheta /. gamma -> 1), "Residual" -> thetaResidual|>]];
emit["B_RADAR_LOG_RESIDUAL", gated[<|"Computed" -> radarFirst,
  "Reference" -> (referenceLog /. gamma -> 1), "Residual" -> radarResidual|>]];
emit["B_DEFLECTION_CONDITION", conditions[[1]]];
emit["B_RADAR_CONDITION", conditions[[2]]];
emit["B_DEFLECTION_EXCHANGE", constrainedExchange[conditions[[1]]]];
emit["B_RADAR_EXCHANGE", constrainedExchange[conditions[[2]]]];

responseDeltas = (Simplify[Normal[Series[#/c0 - 1, {fLocal, 0, 1}]],
  rho0 > 0 && c0 > 0 && Element[{n, s, fLocal}, Reals]] &) /@ responseInputs;
responseMultipliers = Coefficient[#, fLocal] & /@ responseDeltas;
responseDomain[j_] := With[{resp = responseInputs[[j]], change = responseDeltas[[j]]},
  Element[resp, Reals] && resp > 0 && Abs[change] < 1 &&
  If[j == 2, localSoundSquared > 0 && (bulkSoundSquared /. rho -> rho0) > 0,
    If[j == 3, localBulkDensity/rho0 > 0, True]]];
emit["C_RESPONSES", <|"SoundSquared" -> bulkSoundSquared,
  "DeltaSeries" -> responseDeltas, "Multipliers" -> responseMultipliers,
  "Domains" -> Table[responseDomain[j], {j, 3}]|>];
Do[cRules = {ampD -> responseMultipliers[[j]] ampF};
  cDomain = profileDomain && (responseDomain[j] /. fLocal -> densityFraction);
  cConditions = conditionPair[cRules, cDomain];
  emit["C_RESPONSE_" <> ToString[j] <> "_DEFLECTION_CONDITION", cConditions[[1]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_RADAR_CONDITION", cConditions[[2]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_DEFLECTION_EXCHANGE", constrainedExchange[cConditions[[1]], cRules]];
  emit["C_RESPONSE_" <> ToString[j] <> "_RADAR_EXCHANGE", constrainedExchange[cConditions[[2]], cRules]];
  emit["C_RESPONSE_" <> ToString[j] <> "_N_DEPENDENCE", <|
    "DeltaDerivative" -> D[responseDeltas[[j]], n],
    "ResidualDerivatives" -> D[{thetaResidual, radarResidual} /. cRules, n]|>], {j, 3}];

restrictions = <|"FLOW_ONLY" -> {ampD -> 0, ampW -> 0},
  "SPEED_ONLY" -> {ampV -> 0, ampW -> 0}, "TILT_ONLY" -> {ampD -> 0, ampV -> 0}|>;
KeyValueMap[Function[{name, rules}, Module[{cc = conditionPair[rules, profileDomain], prefix},
  prefix = "RESTRICTION_" <> name <> "_";
  emit[prefix <> "PROFILES", rules];
  emit[prefix <> "DEFLECTION", gated[Values[bending] /. rules, rayGate /. rules]];
  emit[prefix <> "RADAR_LOG", gated[Values[radarLog] /. rules, rayGate /. rules]];
  emit[prefix <> "GAMMA_DEFLECTION", gated[gammaTheta /. rules, rayGate /. rules, rules]];
  emit[prefix <> "GAMMA_RADAR", gated[gammaRadar /. rules, rayGate /. rules, rules]];
  emit[prefix <> "GAMMA_DIFFERENCE", gated[(gammaTheta - gammaRadar) /. rules, rayGate /. rules, rules]];
  emit[prefix <> "DEFLECTION_RESIDUAL", gated[thetaResidual /. rules, rayGate /. rules, rules]];
  emit[prefix <> "RADAR_RESIDUAL", gated[radarResidual /. rules, rayGate /. rules, rules]];
  emit[prefix <> "DEFLECTION_CONDITION", cc[[1]]]; emit[prefix <> "RADAR_CONDITION", cc[[2]]];
  emit[prefix <> "DEFLECTION_EXCHANGE", constrainedExchange[cc[[1]], rules]];
  emit[prefix <> "RADAR_EXCHANGE", constrainedExchange[cc[[2]], rules]]]], restrictions];
Do[bulkRules = {ampD -> responseMultipliers[[j]] ampF, ampV -> 0, ampW -> 0};
  bulkDomain = profileDomain && (responseDomain[j] /. fLocal -> densityFraction);
  bulkConditions = conditionPair[bulkRules, bulkDomain];
  emit["C_RESPONSE_" <> ToString[j] <> "_BULK_ONLY_DEFLECTION_CONDITION", bulkConditions[[1]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_BULK_ONLY_RADAR_CONDITION", bulkConditions[[2]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_BULK_ONLY_DEFLECTION_EXCHANGE", constrainedExchange[bulkConditions[[1]], bulkRules]];
  emit["C_RESPONSE_" <> ToString[j] <> "_BULK_ONLY_RADAR_EXCHANGE", constrainedExchange[bulkConditions[[2]], bulkRules]], {j, {2, 3}}];

(* Forward premise, solved before choosing a profile representation. The
   density power law below is a declared family restriction for evaluation,
   with live normalization/exponent; no identification with bulk rho0. *)
forwardMassVector = (massDensity vForward[r] cartesian/r) /. r -> cartRadius;
forwardDiv = Simplify[Total[MapThread[D, {forwardMassVector, cartesian}]] /.
  {x1 -> r, x2 -> 0, x3 -> 0}, $Assumptions];
forwardSolution = DSolve[forwardDiv == 0, vForward, r];
fluxInput = Phi == Integrate[(massDensity vForward[r]) r^2 Sin[theta],
  {theta, 0, Pi}, {phi, 0, 2 Pi}];
forwardV = Simplify[vForward[r] /. First[Solve[fluxInput, vForward[r]]]];
forwardResidual = Simplify[forwardDiv /. {vForward -> Function[{r}, Evaluate[forwardV]]}];
rhoForwardAnsatz = rhoScale (ell/r)^powRho;
forwardFamilyV = Simplify[forwardV /. rhoBr[r] -> rhoForwardAnsatz];
forwardPower = Simplify[-r D[forwardFamilyV, r]/forwardFamilyV];
forwardAmplitude = Simplify[(forwardFamilyV /. r -> ell)/c0];
forwardRules = {ampV -> forwardAmplitude, powV -> forwardPower, rhoBr[r] -> rhoForwardAnsatz};
forwardDomain = (profileDomain /. forwardRules) && rhoScale != 0 &&
  Element[{Phi, powRho, rhoScale}, Reals];
forwardPremise = <|"NormalExchange" -> (jn[r] == 0), "FluxDefinition" -> fluxInput,
  "Measure" -> "coordinate d3x", "DensityFamily" -> (rhoBr[r] == rhoForwardAnsatz)|>;
forwardEmit[name_, z_] := emit["FORWARD_NO_FAR_ZONE_LOSS_" <> name,
  <|"Premise" -> forwardPremise, "Object" -> z|>];
forwardEmit["MASS_SOLUTION", <|"Solution" -> forwardSolution, "FluxVelocity" -> forwardV,
  "DifferentialResidual" -> forwardResidual, "FamilyVelocity" -> forwardFamilyV,
  "Substitution" -> forwardRules, "Domain" -> forwardDomain|>];
Do[stageRules = Join[forwardRules, If[stage == 1, {ampD -> 0, ampW -> 0}, {}]];
  stageGate = rayGate /. stageRules;
  prefix = If[stage == 1, "FLOW_", "LIVE_OPTICS_"];
  forwardEmit[prefix <> "BRANCH_EXISTENCE", Inactive[ForAll][r, rayRadiusDomain[r], physicalPropagation /. stageRules]];
  forwardEmit[prefix <> "PATH_TRAVERSAL", Inactive[ForAll][r, rayRadiusDomain[r], physicalTraversal /. stageRules]];
  forwardEmit[prefix <> "DEFLECTION", gated[Values[bending] /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "RADAR_LOG", gated[Values[radarLog] /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "GAMMA_DEFLECTION", gated[gammaTheta /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "GAMMA_RADAR", gated[gammaRadar /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "GAMMA_DIFFERENCE", gated[(gammaTheta - gammaRadar) /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "DEFLECTION_RESIDUAL", gated[thetaResidual /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "RADAR_RESIDUAL", gated[radarResidual /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "ROUND_TRIP", gated[Values[roundTrip] /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "ONE_WAY_ER", gated[Values[oneWayER] /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "ONE_WAY_RE", gated[Values[oneWayRE] /. stageRules, stageGate, stageRules]];
  forwardEmit[prefix <> "NONRECIPROCAL", gated[Values[nonreciprocal] /. stageRules, stageGate, stageRules]];
  stageConditions = conditionPair[stageRules, forwardDomain];
  forwardEmit[prefix <> "DEFLECTION_CONDITION", stageConditions[[1]]];
  forwardEmit[prefix <> "RADAR_CONDITION", stageConditions[[2]]], {stage, 2}];
Do[forwardCRules = Join[forwardRules, {ampD -> responseMultipliers[[j]] ampF}];
  forwardCDomain = forwardDomain && (responseDomain[j] /. fLocal -> densityFraction);
  forwardCConditions = conditionPair[forwardCRules, forwardCDomain];
  forwardEmit["C_RESPONSE_" <> ToString[j] <> "_DEFLECTION_CONDITION", forwardCConditions[[1]]];
  forwardEmit["C_RESPONSE_" <> ToString[j] <> "_RADAR_CONDITION", forwardCConditions[[2]]], {j, 3}];
emit["SUPPLIED_DEPENDENCIES", <|"A" -> {"DISPERSION", "ADVECTION", "EMBEDDING",
    "ISOTROPIC_SPEED", "LAB_HELD"}, "B" -> {"ORDER_COUNTING", "GM_REFERENCE", "PPN_REFERENCE"},
  "EXCHANGE" -> {"COORDINATE_MASS_BALANCE"}, "C" -> {"BULK_EOS", "DENSITY_RESPONSE"}|>];
emit["OBSERVABLE_DOMAIN", <|"Domain" -> rayDomain, "Gate" -> rayGate,
  "Outside" -> ConditionalExpression["NOT_ESTABLISHED", Not[rayGate]],
  "ProfileFamily" -> {delta, velocity, xiSlope}, "LogComparison" -> ZE/b > 1 && ZR/b > 1|>];
emit["ENGINE_LOCAL_NAMES", localNames];
