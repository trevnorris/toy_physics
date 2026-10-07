(* S9b blind builder.  Sole physical input: S9b_SHARED_PHYSICS.md v8.
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
  AssociateTo[emitted, tag -> True];
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
    powD > 0 && powV > 0 && powW > 0 && GM > 0 && rho0 > 0 && m > 0;

(* Supplied equations and chosen ansatz: the only hand-combined physics.
   powW is the exponent of xi', so the embedding grade has exponent 2 powW.
   rhoBr is an arbitrary live radial function, not the bulk number density.
*)
delta = ampD (ell/r)^powD;
velocity = c0 ampV (ell/r)^powV;
xiSlope = ampW (ell/r)^powW;
densityFraction = ampF (ell/r)^powD;
speed = c0 (1 + xD delta);
radialMetric = 1 + xW xiSlope^2;
metric3 = DiagonalMatrix[{radialMetric, r^2, r^2 Sin[theta]^2}];
metric2 = DiagonalMatrix[{aMetric, r^2}];
kCov = {kr, kphi}; vVector = {vLocal, 0};
dispersion = (omega - vVector.kCov)^2 - cSquared kCov.Inverse[metric2].kCov;
localSpeedIdentification = cSquared == muPerp[r]/rhoBr[r];
farSpeedIdentification = c0^2 == muPerpInfinity/rhoBrInfinity;
massBalanceInput = divRhoV == -jn[r];
referenceTheta = (1 + gamma) 2 GM/(b c0^2);
referenceRadar = 2 (1 + gamma) GM/c0^3 Log[4 rE rR/b^2];
bulkPressure = K rho^n;
bulkSoundSquared = D[bulkPressure, rho]/m;
responseInputs = {c0, c0 Sqrt[(bulkSoundSquared /. rho -> rho0 (1 + fLocal))/
    (bulkSoundSquared /. rho -> rho0)], c0 (1 + fLocal)^s};

emit["SUPPLIED_INPUTS", <|"Dispersion" -> dispersion == 0,
  "Speed" -> localSpeedIdentification, "FarSpeed" -> farSpeedIdentification,
  "Metric" -> metric3, "Balance" -> massBalanceInput,
  "References" -> {referenceTheta, referenceRadar, gamma == 1},
  "Pressure" -> bulkPressure, "Responses" -> responseInputs|>];
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
localAssumptions = aMetric > 0 && cSquared > 0 && r > 0 && kr > 0 &&
  Element[{aMetric, cSquared, r, kr, kphi, vLocal}, Reals];
frequencyRoots = omega /. Solve[dispersion == 0, omega];
positiveRoot = Select[frequencyRoots,
  TrueQ[Simplify[(# /. vLocal -> 0) > 0, localAssumptions]] &][[1]];
groupVelocity = D[positiveRoot, #] & /@ kCov;
relativeNorm = Simplify[(groupVelocity - vVector).metric2.
  (groupVelocity - vVector), localAssumptions];
tangent = {dr, dphi};
timePolynomial = Expand[(tangent - tau vVector).metric2.
  (tangent - tau vVector) - tau^2 relativeNorm];
timeCoefficients = Coefficient[timePolynomial, tau, #] & /@ Range[0, 2];
timeDiscriminant = Simplify[timeCoefficients[[2]]^2 -
  4 timeCoefficients[[1]] timeCoefficients[[3]]];
timeRoot = (-timeCoefficients[[2]] - Sqrt[timeDiscriminant])/
  (2 timeCoefficients[[3]]);
subcritical = cSquared - aMetric vLocal^2 > 0;
fermatOdd = Simplify[(timeRoot - (timeRoot /. vLocal -> -vLocal))/2,
  localAssumptions && subcritical];
fermatEvenSquared = Simplify[((timeRoot + (timeRoot /. vLocal -> -vLocal))/2)^2,
  localAssumptions && subcritical];
opticalMetric = Simplify[Table[D[fermatEvenSquared, tangent[[i]], tangent[[j]]]/2,
  {i, 2}, {j, 2}], localAssumptions && subcritical];
radialSpeeds = Simplify[(D[#, kr] & /@ frequencyRoots) /. kphi -> 0,
  localAssumptions];
propagatingLocal = cSquared > 0 && aMetric > 0;
traversalLocal = Simplify[And @@ {Min @@ radialSpeeds < 0,
  Max @@ radialSpeeds > 0}, localAssumptions];
localRules = {aMetric -> radialMetric, cSquared -> speed^2,
  vLocal -> xV velocity};
physicalPropagation = atOne[propagatingLocal /. localRules];
physicalTraversal = atOne[subcritical /. localRules];
emit["BRANCH_FREQUENCIES", frequencyRoots];
emit["BRANCH_EXISTENCE", <|"Local" -> propagatingLocal,
  "AlongRay" -> Inactive[ForAll][r, rayRadiusDomain[r],
    physicalPropagation && atOne[speed] > 0 && rhoBr[r] != 0],
  "SpeedIdentification" -> localSpeedIdentification|>];
emit["PATH_TRAVERSAL", <|"RadialGroupVelocities" -> radialSpeeds,
  "DirectionalCondition" -> traversalLocal,
  "SubcriticalDomain" -> Inactive[ForAll][r, rayRadiusDomain[r], physicalTraversal]|>];
emit["BRANCH_TYPES", <|
  "Complex" -> ConditionalExpression[frequencyRoots, aMetric > 0 && cSquared < 0],
  "ComplexTypes" -> (Sign[Im[#]] & /@ frequencyRoots),
  "Absent" -> (cSquared == 0 || aMetric <= 0 || rhoBr[r] == 0),
  "UnableToTraverse" -> (propagatingLocal && Not[subcritical]),
  "ObservablesOutsideDomain" -> <|"Status" -> "NOT_ESTABLISHED",
    "Missing" -> "BRANCH_OR_TRAVERSAL"|>|>];
emit["LOCAL_HAMILTON_VELOCITY", <|"Velocity" -> groupVelocity,
  "RelativeMetricNorm" -> relativeNorm|>];
emit["LOCAL_FERMAT_OBJECT", <|"TimePolynomial" -> timePolynomial,
  "TimeRoot" -> timeRoot, "EvenMetric" -> opticalMetric, "Odd" -> fermatOdd|>];

(* Normalize the even optical metric to far-field length units.  R(r) is
   its circumferential radius.  The push-forward of a(r) dr under R=r+u(r)
   has formal density Sum[(-D_R)^j (a u^j)/j!].  Three iterations suffice:
   every even perturbation contains xD, xV^2 or xW; no fourth product
   survives the requested box.  This includes the shift of the periapsis.
*)
optical = Simplify[(c0^2 opticalMetric) /. localRules, $Assumptions];
circumferenceRadius = c0 r/Sqrt[(cSquared - aMetric vLocal^2) /. localRules];
(* Derive the same coordinate from the computed angular metric, keeping
   the analytic branch connected to positive r at zero perturbation. *)
radiusShift = box[r (Sqrt[box[optical[[2, 2]]/r^2]] - 1)];
radialLengthDensity = box[Sqrt[box[optical[[1, 1]]]]];
opticalDensity = radialLengthDensity;
Do[opticalDensity += (-1)^j/Factorial[j] D[
    project[radialLengthDensity project[radiusShift^j]], {r, j}], {j, 1, 3}];
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
Do[upperCorrection += project[project[endpointShift^j] D[
    endpointDensity Sqrt[1 - impact^2/rr^2], {rr, j - 1}]]/Factorial[j], {j, 1, 3}];
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
emit["NONRECIPROCAL_DEPENDENCE", <|"OneForm" -> oddComponents,
  "ExteriorDerivative" -> oddExterior,
  "ExactnessTest" -> Simplify[And @@ Thread[Flatten[oddExterior] == 0]],
  "PrimitiveKernels" -> {primitiveGeneral, primitiveResonance},
  "Domain" -> r > 0|>];

(* Every grade has its own tag, including the unperturbed excess and
   grades which the calculation cancels.  Conditional validity is shared
   by the entire block through RAY_DOMAIN and BRANCH_TYPES. *)
Do[emit["A_DEFLECTION_G" <> label[g], bending[label[g]]];
  emit["A_ROUND_TRIP_G" <> label[g], 2 evenTimes[label[g]]];
  emit["A_ONE_WAY_ER_G" <> label[g], evenTimes[label[g]] + oddTimes[label[g]]];
  emit["A_ONE_WAY_RE_G" <> label[g], evenTimes[label[g]] - oddTimes[label[g]]];
  emit["A_NONRECIPROCAL_G" <> label[g], oddTimes[label[g]]], {g, grades}];

(* For positive power-law exponents the large-endpoint radial kernel
   has a logarithm only on its exponent-one stratum.  Compute its
   logarithmic coefficient from the resonant integral, not from a
   reference value.  Nonresonant tails remain in the full Part A times.
*)
resonantMap = impact Cosh[hyperbolicParameter];
resonantIntegrand = FullSimplify[(Sqrt[1 - impact^2/resonantMap^2]/resonantMap)
  D[resonantMap, hyperbolicParameter], impact > 0 && hyperbolicParameter > 0];
resonantKernel = Integrate[resonantIntegrand,
  {hyperbolicParameter, 0, ArcCosh[rr/impact]},
  Assumptions -> rr > impact > 0, GenerateConditions -> False];
resonantAsymptotic = FullSimplify[Normal[Series[resonantKernel /. rr -> impact/zeta,
  {zeta, 0, 0}]], impact > 0 && 0 < zeta < 1];
(* On the retained constant-plus-log asymptotic expression, zeta D_zeta
   extracts the log coefficient even when the CAS combines its argument.
   Log[zeta]=Log[impact]-Log[endpoint]; the reported basis is Log[1/b^2]. *)
logPerEndpoint = FullSimplify[zeta D[resonantAsymptotic, zeta]/(-2),
  impact > 0 && 0 < zeta < 1];
radarLog = Association[];
Do[AssociateTo[radarLog, label[g] ->
  Simplify[2 Length[endpointRadii] logPerEndpoint hCoefficient[g] ell/c0
    Piecewise[{{1, exponent[g] == 1}}, 0], $Assumptions]], {g, firstGrades}];
emit["LOCAL_RADAR_LOG_EXTRACTION", <|"ResonantIntegral" -> resonantKernel,
  "LargeEndpointExpansion" -> resonantAsymptotic,
  "LogCoefficientPerEndpoint" -> logPerEndpoint,
  "TailExponents" -> (exponent /@ firstGrades),
  "NonresonantPrimitive" -> primitiveGeneral|>];
thetaFirst = Total[bending[label[#]] & /@ firstGrades];
radarFirst = Total[radarLog[label[#]] & /@ firstGrades];
referenceLog = Coefficient[Expand[referenceRadar /. Log[4 rE rR/b^2] -> logBasis], logBasis];
gammaTheta = gamma /. First[Solve[referenceTheta == thetaFirst, gamma]];
gammaRadar = gamma /. First[Solve[referenceLog == radarFirst, gamma]];
thetaResidual = thetaFirst - (referenceTheta /. gamma -> 1);
radarResidual = radarFirst - (referenceLog /. gamma -> 1);
everyB[z_] := Inactive[ForAll][b, b > bFar, z == 0];
conditions = {everyB[thetaResidual], everyB[radarResidual]};
emit["B_FIRST_ORDER_DEFLECTION", Table[{g, bending[label[g]], g.{1, 1/2, 1}}, {g, firstGrades}]];
emit["B_FIRST_ORDER_RADAR_LOG", Table[{g, radarLog[label[g]], g.{1, 1/2, 1}}, {g, firstGrades}]];
emit["B_GAMMA_DEFLECTION", gammaTheta];
emit["B_GAMMA_RADAR", gammaRadar];
emit["B_GAMMA_DIFFERENCE", gammaTheta - gammaRadar];
emit["B_DEFLECTION_RESIDUAL", <|"Computed" -> thetaFirst,
  "Reference" -> (referenceTheta /. gamma -> 1), "Residual" -> thetaResidual|>];
emit["B_RADAR_LOG_RESIDUAL", <|"Computed" -> radarFirst,
  "Reference" -> (referenceLog /. gamma -> 1), "Residual" -> radarResidual|>];
emit["B_DEFLECTION_CONDITION", conditions[[1]]];
emit["B_RADAR_CONDITION", conditions[[2]]];

(* Induced-volume mass balance.  rhoBr stays inside the derivative.
   Conditions constrain ampV^2; the signed flow and density remain free.
   Hence the implied exchange is a constrained set of radial functions,
   not a uniquely selected jn.  This preserves inward/outward ambiguity. *)
volume = Simplify[Sqrt[Det[metric3]], $Assumptions && 0 < theta < Pi && xW >= 0];
divergence = Simplify[D[volume rhoBr[r] xV velocity, r]/volume,
  $Assumptions && 0 < theta < Pi && xW >= 0];
exchange = jn[r] /. First[Solve[massBalanceInput /. divRhoV -> divergence, jn[r]]];
exchangeBox = box[exchange];
exchangeFlat = exchangeBox /. xW -> 0;
emit["MASS_BALANCE", <|"Volume" -> volume, "Divergence" -> divergence,
  "Exchange" -> atOne[exchangeBox], "FlatExchange" -> atOne[exchangeFlat],
  "InducedCorrection" -> atOne[exchangeBox - exchangeFlat],
  "Grades" -> Table[{g, coeff[exchangeBox, g]}, {g, grades}]|>];
constrainedExchange[condition_] := <|"ProfileCondition" -> condition,
  "ExchangeRelation" -> (jn[r] == atOne[exchangeBox]),
  "DensityDomain" -> rhoBr[r] != 0|>;
emit["B_DEFLECTION_EXCHANGE", constrainedExchange[conditions[[1]]]];
emit["B_RADAR_EXCHANGE", constrainedExchange[conditions[[2]]]];

(* Part C: only the density response is first order in f.  Its exponent
   and coefficient remain free, as do the velocity and embedding. *)
responseDeltas = (Simplify[Normal[Series[#/c0 - 1, {fLocal, 0, 1}]],
  rho0 > 0 && c0 > 0 && Element[{n, s, fLocal}, Reals]] &) /@ responseInputs;
responseMultipliers = Coefficient[#, fLocal] & /@ responseDeltas;
emit["C_RESPONSES", <|"SoundSquared" -> bulkSoundSquared,
  "DeltaSeries" -> responseDeltas, "Multipliers" -> responseMultipliers|>];
Do[cRules = {ampD -> responseMultipliers[[j]] ampF};
  cResiduals = {thetaResidual, radarResidual} /. cRules;
  cConditions = everyB /@ cResiduals;
  emit["C_RESPONSE_" <> ToString[j] <> "_DEFLECTION_CONDITION", cConditions[[1]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_RADAR_CONDITION", cConditions[[2]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_DEFLECTION_EXCHANGE", constrainedExchange[cConditions[[1]]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_RADAR_EXCHANGE", constrainedExchange[cConditions[[2]]]];
  emit["C_RESPONSE_" <> ToString[j] <> "_N_DEPENDENCE", <|
    "DeltaDerivative" -> D[responseDeltas[[j]], n],
    "ResidualDerivatives" -> D[cResiduals, n]|>], {j, 1, Length[responseInputs]}];
emit["SUPPLIED_DEPENDENCIES", <|"A" -> {"DISPERSION", "ADVECTION", "EMBEDDING",
    "ISOTROPIC_SPEED", "LAB_HELD"}, "B" -> {"ORDER_COUNTING", "GM_REFERENCE", "PPN_REFERENCE"},
  "EXCHANGE" -> {"INDUCED_METRIC", "MASS_BALANCE"}, "C" -> {"BULK_EOS", "DENSITY_RESPONSE"}|>];
emit["OBSERVABLE_DOMAIN", <|"Domain" -> rayDomain,
  "Outside" -> <|"Status" -> "NOT_ESTABLISHED", "Missing" -> "BRANCH_OR_TRAVERSAL"|>,
  "ProfileFamily" -> {delta, velocity, xiSlope},
  "LogComparison" -> ZE/b > 1 && ZR/b > 1|>];
emit["ENGINE_LOCAL_NAMES", localNames];
Quit[0];
