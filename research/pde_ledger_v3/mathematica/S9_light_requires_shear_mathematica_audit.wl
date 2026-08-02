(* S9 independent Mathematica audit: all reported physics is derived below. *)

ClearAll["Global`*"];
$HistoryLength = 0;

coordinates = {x, y, z};
allVariables = {t, x, y, z};
fieldHeads = {u1, u2, u3};
uVector = Through[fieldHeads[t, x, y, z]];
kVector = {kx, ky, kz};
amplitudeVector = {ax, ay, az};
kSquared = Expand[kVector . kVector];
planePhase = Exp[I (kVector . coordinates - omega t)];

baseAssumptions =
  rhoBr > 0 && muR > 0 && kSquared > 0 &&
   Element[Join[kVector, amplitudeVector, {omega, k}], Reals];

ClearAll[deriveFromAction];
deriveFromAction[inertiaCoefficient_, stiffnessCoefficient_] := Module[
  {
   curlU, lagrangian, variationalResidual, eomResidual,
   waveFunctions, waveRules, planeWaveResidual, dynamicMatrix,
   dynamicMatrixOmegaSquared, determinant, omegaRules,
   transverseGenerator, transverseRules, omegaSquaredResult,
   localAssumptions, derivationOK
   },

  localAssumptions =
    baseAssumptions && inertiaCoefficient > 0 &&
     stiffnessCoefficient > 0;

  curlU = Curl[uVector, coordinates];
  lagrangian =
    (1/2) inertiaCoefficient (D[uVector, t] . D[uVector, t]) -
     (1/2) stiffnessCoefficient (curlU . curlU);

  (* Euler-Lagrange variation with respect to every field component. *)
  variationalResidual = FullSimplify[
    Table[
     D[lagrangian, uVector[[component]]] -
      Sum[
       D[
        D[lagrangian,
         D[uVector[[component]], allVariables[[variable]]]],
        allVariables[[variable]]
        ],
       {variable, Length[allVariables]}
       ],
     {component, Length[uVector]}
     ],
    localAssumptions
    ];

  (* Multiplication by -1 gives the conventional positive-inertia EOM. *)
  eomResidual = FullSimplify[-variationalResidual, localAssumptions];

  waveFunctions = MapThread[
    Function[{amplitude},
      Function[{tt, xx, yy, zz},
       Evaluate[
        amplitude Exp[I (kx xx + ky yy + kz zz - omega tt)]]
       ]
      ],
    {amplitudeVector}
    ];
  waveRules = Thread[fieldHeads -> waveFunctions];

  planeWaveResidual = FullSimplify[
    (eomResidual /. waveRules)/planePhase,
    localAssumptions
    ];

  dynamicMatrix = Table[
    Coefficient[planeWaveResidual[[row]], amplitudeVector[[column]]],
    {row, Length[amplitudeVector]},
    {column, Length[amplitudeVector]}
    ];
  dynamicMatrixOmegaSquared =
    Expand[dynamicMatrix] /. omega^2 -> omegaSquaredSymbol;
  determinant = Factor[Det[dynamicMatrixOmegaSquared]];
  omegaRules = DeleteDuplicates[
    Solve[determinant == 0, omegaSquaredSymbol]
    ];

  (* k^2 I-k k^T spans the transverse image without dividing by k^2. *)
  transverseGenerator =
    kSquared IdentityMatrix[Length[kVector]] -
     Outer[Times, kVector, kVector];
  transverseRules = Select[
    omegaRules,
    TrueQ[FullSimplify[
        (dynamicMatrixOmegaSquared /. #) . transverseGenerator ==
         ConstantArray[0, {Length[kVector], Length[kVector]}],
        localAssumptions
        ]] &
    ];

  derivationOK = Length[transverseRules] == 1;
  omegaSquaredResult = If[
    derivationOK,
    FullSimplify[
     omegaSquaredSymbol /. First[transverseRules],
     localAssumptions
     ],
    Missing["NoUniqueTransverseRoot"]
    ];

  <|
   "Lagrangian" -> lagrangian,
   "VariationalResidual" -> variationalResidual,
   "EOMResidual" -> eomResidual,
   "PlaneWaveResidual" -> planeWaveResidual,
   "DynamicMatrix" -> dynamicMatrixOmegaSquared,
   "Determinant" -> determinant,
   "AllOmegaRules" -> omegaRules,
   "TransverseRules" -> transverseRules,
   "OmegaSquared" -> omegaSquaredResult,
   "DerivationOK" -> derivationOK
   |>
  ];

ClearAll[inspectDispersion];
inspectDispersion[omegaSquaredExpression_] := Module[
  {omegaSquaredInK, speedSquaredCandidate, isLinear},
  omegaSquaredInK =
    FullSimplify[Factor[omegaSquaredExpression] /. kSquared -> k^2,
     baseAssumptions];
  speedSquaredCandidate =
    FullSimplify[omegaSquaredInK/k^2, baseAssumptions];
  isLinear = TrueQ[
    FullSimplify[
      omegaSquaredInK == speedSquaredCandidate k^2,
      baseAssumptions
      ] && FreeQ[speedSquaredCandidate, k]
    ];
  <|
   "OmegaSquaredInK" -> omegaSquaredInK,
   "IsLinear" -> isLinear,
   "SpeedSquared" ->
    If[isLinear, speedSquaredCandidate, Missing["NotApplicable"]]
   |>
  ];

mainDerivation = deriveFromAction[rhoBr, muR];
mainDispersion = inspectDispersion[mainDerivation["OmegaSquared"]];

(* Dimensional derivation in the requested {L,T,M} exponent convention. *)
dimLength = {1, 0, 0};
dimTime = {0, 1, 0};
dimMass = {0, 0, 1};
dimAcceleration = dimLength - 2 dimTime;
dimForce = dimMass + dimAcceleration;
dimEnergy = dimForce + dimLength;
dimEnergyDensity3D = dimEnergy - 3 dimLength;
dimDisplacement = dimLength;
dimTimeDerivativeU = dimDisplacement - dimTime;
dimCurlU = dimDisplacement - dimLength;

dimRhoUnknown = {rhoLengthExponent, rhoTimeExponent, rhoMassExponent};
dimMuUnknown = {muLengthExponent, muTimeExponent, muMassExponent};
dimensionUnknowns = Join[dimRhoUnknown, dimMuUnknown];
dimensionEquations = Join[
   Thread[dimRhoUnknown + 2 dimTimeDerivativeU == dimEnergyDensity3D],
   Thread[dimMuUnknown + 2 dimCurlU == dimEnergyDensity3D]
   ];
dimensionSolutions = Solve[dimensionEquations, dimensionUnknowns];
dimensionDerivationOK = Length[dimensionSolutions] == 1;
dimensionRule = If[
   dimensionDerivationOK,
   First[dimensionSolutions],
   {}
   ];
dimRhoBr = dimRhoUnknown /. dimensionRule;
dimMuR = dimMuUnknown /. dimensionRule;

ClearAll[symbolPower];
symbolPower[expression_, symbol_] := Module[{rationalExpression},
  rationalExpression = Together[expression];
  Exponent[Numerator[rationalExpression], symbol] -
   Exponent[Denominator[rationalExpression], symbol]
  ];

dimReportedSpeedSquared =
  symbolPower[mainDispersion["SpeedSquared"], rhoBr] dimRhoBr +
   symbolPower[mainDispersion["SpeedSquared"], muR] dimMuR;
dimExpectedSpeedSquared = 2 (dimLength - dimTime);
dimensionAdmissible =
  TrueQ[dimReportedSpeedSquared == dimExpectedSpeedSquared];

(* Able-to-fail controls: change each action coefficient and re-run everything. *)
inertiaControlDerivation = deriveFromAction[2 rhoBr, muR];
inertiaControlDispersion =
  inspectDispersion[inertiaControlDerivation["OmegaSquared"]];
stiffnessControlDerivation = deriveFromAction[rhoBr, 3 muR];
stiffnessControlDispersion =
  inspectDispersion[stiffnessControlDerivation["OmegaSquared"]];

inertiaControlFired = TrueQ[
   inertiaControlDerivation["DerivationOK"] &&
    inertiaControlDispersion["IsLinear"] &&
    FullSimplify[
     inertiaControlDispersion["SpeedSquared"] !=
      mainDispersion["SpeedSquared"],
     baseAssumptions
     ]
   ];
stiffnessControlFired = TrueQ[
   stiffnessControlDerivation["DerivationOK"] &&
    stiffnessControlDispersion["IsLinear"] &&
    FullSimplify[
     stiffnessControlDispersion["SpeedSquared"] !=
      mainDispersion["SpeedSquared"],
     baseAssumptions
     ]
   ];

ClearAll[printAuditLine];
printAuditLine[label_, value_] :=
  Print[label <>
    If[StringQ[value], value, ToString[value, InputForm]]];

(* Readable intermediate algebra. *)
printAuditLine["WL_S9_INTERMEDIATE_ACTION: ",
  mainDerivation["Lagrangian"]];
printAuditLine["WL_S9_INTERMEDIATE_VARIATIONAL_RESIDUAL: ",
  mainDerivation["VariationalResidual"]];
printAuditLine["WL_S9_INTERMEDIATE_PLANE_WAVE: ",
  amplitudeVector planePhase];
printAuditLine["WL_S9_INTERMEDIATE_PLANE_WAVE_RESIDUAL: ",
  mainDerivation["PlaneWaveResidual"]];
printAuditLine["WL_S9_INTERMEDIATE_DYNAMIC_MATRIX: ",
  mainDerivation["DynamicMatrix"]];
printAuditLine["WL_S9_INTERMEDIATE_DETERMINANT: ",
  mainDerivation["Determinant"]];
printAuditLine["WL_S9_INTERMEDIATE_ALL_OMEGA_SQUARED_ROOTS: ",
  mainDerivation["AllOmegaRules"]];
printAuditLine["WL_S9_INTERMEDIATE_TRANSVERSE_ROOT: ",
  mainDerivation["TransverseRules"]];
printAuditLine["WL_S9_INTERMEDIATE_DIM_ENERGY_DENSITY_3D: ",
  dimEnergyDensity3D];
printAuditLine["WL_S9_INTERMEDIATE_DIMENSION_EQUATIONS: ",
  dimensionEquations];
printAuditLine["WL_S9_INTERMEDIATE_DIM_SPEED_SQUARED_EXPECTED: ",
  dimExpectedSpeedSquared];

(* Required output tokens. *)
printAuditLine["WL_S9_EOM: ",
  Thread[mainDerivation["EOMResidual"] == ConstantArray[0, 3]]];
printAuditLine["WL_S9_OMEGA_SQUARED: ",
  mainDispersion["OmegaSquaredInK"]];
printAuditLine["WL_S9_DISPERSION_FORM: ",
  If[mainDispersion["IsLinear"],
   "LINEAR_IN_K",
   "OTHER " <> ToString[mainDispersion["OmegaSquaredInK"], InputForm]
   ]];
printAuditLine["WL_S9_SPEED_SQUARED: ",
  If[mainDispersion["IsLinear"],
   mainDispersion["SpeedSquared"],
   "NOT_APPLICABLE"
   ]];
printAuditLine["WL_S9_DIM_RHOBR: ", dimRhoBr];
printAuditLine["WL_S9_DIM_MUR: ", dimMuR];
printAuditLine["WL_S9_DIM_SPEED: ", dimReportedSpeedSquared];

(* Control results are themselves derived values, never expected-value literals. *)
printAuditLine["WL_S9_CONTROL_INERTIA_ACTION: ",
  inertiaControlDerivation["Lagrangian"]];
printAuditLine["WL_S9_CONTROL_INERTIA_SPEED_SQUARED: ",
  inertiaControlDispersion["SpeedSquared"]];
printAuditLine["WL_S9_CONTROL_INERTIA_FIRED: ",
  If[inertiaControlFired, "TRUE", "FALSE"]];
printAuditLine["WL_S9_CONTROL_STIFFNESS_ACTION: ",
  stiffnessControlDerivation["Lagrangian"]];
printAuditLine["WL_S9_CONTROL_STIFFNESS_SPEED_SQUARED: ",
  stiffnessControlDispersion["SpeedSquared"]];
printAuditLine["WL_S9_CONTROL_STIFFNESS_FIRED: ",
  If[stiffnessControlFired, "TRUE", "FALSE"]];
printAuditLine["WL_S9_DIMENSION_ADMISSIBLE_AS_SPEED_SQUARED: ",
  If[dimensionAdmissible, "TRUE", "FALSE"]];

checks = <|
   "main transverse solve" -> mainDerivation["DerivationOK"],
   "linear dispersion inspection" -> mainDispersion["IsLinear"],
   "unique dimensional solve" -> dimensionDerivationOK,
   "speed-squared dimensional admissibility" -> dimensionAdmissible,
   "inertia perturbation control" -> inertiaControlFired,
   "stiffness perturbation control" -> stiffnessControlFired
   |>;
failedChecks = Keys[Select[checks, Not[TrueQ[#]] &]];
verdict = If[failedChecks === {}, "PASS", "FAIL"];
printAuditLine["WL_S9_VERDICT: ", verdict];

If[failedChecks =!= {},
  Scan[printAuditLine["WL_S9_FAILURE: ", #] &, failedChecks];
  Exit[1]
  ];

Exit[0];
