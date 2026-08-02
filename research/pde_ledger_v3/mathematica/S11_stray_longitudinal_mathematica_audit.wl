(* S11 blind Mathematica audit: standalone, print-only, and imports no files. *)

ClearAll["Global`*"];

printTag[tag_String, value_] :=
  Print[tag <> ": " <> ToString[value, InputForm]];

printBooleanTag[tag_String, value_] :=
  Print[tag <> ": " <> If[TrueQ[value], "TRUE", "FALSE"]];

positiveAssumptions =
  rhoBr > 0 && muR > 0 && Bcomp > 0 && muBr > 0 && cs0 > 0 &&
  kSq > 0 && Element[{rhoBr, muR, Bcomp, muBr, cs0, kSq}, Reals];

(* Construct the infinitesimal SO(d) action A |-> R.A.Transpose[R].
   A symmetric quadratic form Q is invariant exactly when G^T.Q+Q.G=0
   for every rotation generator G.  Solving those equations also retains
   the low-dimensional coincidences among irreducible pieces. *)
invariantQuadraticCount[d_Integer] := Module[
  {n, matrixBasis, rotationGenerators, representationGenerators, pairs,
   qVariables, qMatrix, equations, coefficientMatrix},
  n = d^2;
  matrixBasis = Table[
    SparseArray[{p -> 1}, {n}] // Normal // Partition[#, d] &,
    {p, n}
  ];
  rotationGenerators = Table[
    SparseArray[{{p[[1]], p[[2]]} -> 1, {p[[2]], p[[1]]} -> -1}, {d, d}],
    {p, Subsets[Range[d], {2}]}
  ];
  representationGenerators = Table[
    Transpose@Table[
      Flatten@Normal[generator . basisMatrix - basisMatrix . generator],
      {basisMatrix, matrixBasis}
    ],
    {generator, rotationGenerators}
  ];
  pairs = Flatten[Table[{i, j}, {i, 1, n}, {j, i, n}], 1];
  qVariables = Array[qCoefficient, Length[pairs]];
  qMatrix = ConstantArray[0, {n, n}];
  Do[
    qMatrix[[Sequence @@ pairs[[p]]]] = qVariables[[p]];
    qMatrix[[Sequence @@ Reverse[pairs[[p]]]]] = qVariables[[p]],
    {p, Length[pairs]}
  ];
  equations = Flatten[
    Transpose[#] . qMatrix + qMatrix . # & /@
      representationGenerators
  ];
  coefficientMatrix = Last@CoefficientArrays[equations, qVariables];
  Length[qVariables] - MatrixRank[coefficientMatrix]
];

(* The three canonical matrix pieces used in the construction. *)
matrixPieces[a_?MatrixQ] := Module[{d = Length[a], trace, antisymmetric,
    symmetricTraceless},
  trace = Tr[a] IdentityMatrix[d]/d;
  antisymmetric = (a - Transpose[a])/2;
  symmetricTraceless = (a + Transpose[a])/2 - trace;
  {trace, antisymmetric, symmetricTraceless}
];

pieceConstructionChecks = Table[
  Module[{a, pieces},
    a = Array[aEntry, {d, d}];
    pieces = matrixPieces[a];
    TrueQ[Simplify[Total[pieces] == a]] &&
      And @@ Flatten@Table[
        TrueQ[Simplify[Tr[Transpose[pieces[[i]]] . pieces[[j]]] == 0]],
        {i, 1, 3}, {j, i + 1, 3}
      ]
  ],
  {d, 2, 5}
];

invariantCounts = Table[d -> invariantQuadraticCount[d], {d, 2, 5}];

(* Build the dynamical matrix as the Hessian of the plane-wave potential
   energy with respect to its amplitude.  Omitting the common factor I in
   the gradient is harmless because the real quadratic amplitude is used. *)
originalMatrixFromEnergy[d_Integer, kVector_List] := Module[
  {amplitude, gradient, curlSquared, divergence, potential},
  amplitude = Array[amplitudeComponent, d];
  gradient = Outer[Times, kVector, amplitude];
  curlSquared = Sum[
    (gradient[[i, j]] - gradient[[j, i]])^2,
    {i, 1, d}, {j, 1, d}
  ]/2;
  divergence = Tr[gradient];
  potential = muR curlSquared/2 + Bcomp divergence^2/2;
  Table[
    D[potential, amplitude[[i]], amplitude[[j]]],
    {i, 1, d}, {j, 1, d}
  ] // Simplify
];

alignedDirection[d_Integer] := UnitVector[d, 1];
alignedOriginalMatrix[d_Integer] :=
  originalMatrixFromEnergy[d, Sqrt[kSq] alignedDirection[d]] // Simplify;

originalClosedMatrix[d_Integer] :=
  muR kSq IdentityMatrix[d] +
    (Bcomp - muR) kSq Outer[Times, alignedDirection[d], alignedDirection[d]];

dynamicalMatrixDerivationChecks = Table[
  TrueQ@FullSimplify[
    alignedOriginalMatrix[d] == originalClosedMatrix[d],
    positiveAssumptions
  ],
  {d, 2, 5}
];

(* Component form of the general, unaligned matrix. *)
dynamicalMatrixOutput =
  M[i, j] -> muR KroneckerDelta[i, j] Sum[k[ell]^2, {ell, 1, D}] +
    (Bcomp - muR) k[i] k[j];

zeroVectorQ[vector_List, assumptions_] :=
  And @@ (TrueQ[FullSimplify[# == 0, assumptions]] & /@ vector);

spectrumFromMatrix[matrix_List, direction_List, assumptions_] := Module[
  {d, characteristic, roots, details, reducedMatrix, kernel,
   perpendicularQ, parallelQ, orientation},
  d = Length[matrix];
  characteristic = Det[matrix - rhoBr wSquared IdentityMatrix[d]] // Factor;
  roots = DeleteDuplicates[
    FullSimplify[wSquared /. Solve[characteristic == 0, wSquared], assumptions],
    TrueQ[FullSimplify[#1 == #2, assumptions]] &
  ];
  details = Table[
    reducedMatrix = FullSimplify[
      matrix - rhoBr root IdentityMatrix[d], assumptions
    ];
    kernel = NullSpace[reducedMatrix];
    perpendicularQ = And @@
      (TrueQ[FullSimplify[direction . # == 0, assumptions]] & /@ kernel);
    parallelQ = And @@
      (zeroVectorQ[# - (direction . #) direction, assumptions] & /@ kernel);
    orientation = Which[
      perpendicularQ && ! parallelQ, PERPENDICULAR,
      parallelQ && ! perpendicularQ, PARALLEL,
      True, UNRESOLVED
    ];
    {root, Length[kernel], orientation, kernel,
     TrueQ[FullSimplify[
       (characteristic /. wSquared -> root) == 0, assumptions
     ]]},
    {root, roots}
  ];
  {characteristic, details}
];

spectrumRecords[spectrum_] :=
  (#[[1]] -> {#[[2]], #[[3]]} & /@ spectrum[[2]]);

rootWithOrientation[spectrum_, requested_] :=
  (SelectFirst[spectrum[[2]], #[[3]] === requested &])[[1]];

originalSpectra = Table[
  spectrumFromMatrix[
    alignedOriginalMatrix[d], alignedDirection[d],
    positiveAssumptions && Bcomp != muR
  ],
  {d, 2, 5}
];

originalPerpendicularRoot = rootWithOrientation[originalSpectra[[2]], PERPENDICULAR];
originalParallelRoot = rootWithOrientation[originalSpectra[[2]], PARALLEL];

perpendicularBcompResidual = FullSimplify[
  D[originalPerpendicularRoot, Bcomp], positiveAssumptions
];
parallelMuRResidual = FullSimplify[
  D[originalParallelRoot, muR], positiveAssumptions
];
perpendicularDependsOnBcomp = ! TrueQ[perpendicularBcompResidual == 0];
parallelDependsOnMuR = ! TrueQ[parallelMuRResidual == 0];

degeneracyCondition = FullSimplify[
  Bcomp == (Bcomp /. First@Solve[
    originalParallelRoot == originalPerpendicularRoot, Bcomp
  ]),
  positiveAssumptions
];

(* Dimension vectors are ordered [L,T,M]. *)
energyDimension = {2, -2, 1};
energyDensityDimension = energyDimension - {D, 0, 0};
displacementDimension = {1, 0, 0};
timeDerivativeDimension = displacementDimension - {0, 1, 0};
spatialDerivativeDimension = displacementDimension - {1, 0, 0};
dimRhoBr = Simplify[
  energyDensityDimension - 2 timeDerivativeDimension
];
dimMuR = Simplify[
  energyDensityDimension - 2 spatialDerivativeDimension
];
dimBcomp = Simplify[
  energyDensityDimension - 2 spatialDerivativeDimension
];
dimLongitudinalSpeed = Simplify[(dimBcomp - dimRhoBr)/2];
dimensionOutput[dimension_] :=
  {generalD -> dimension, atD3 -> (dimension /. D -> 3)};

(* These are dimensional relations supplied by the physical meanings of the
   variables, rather than inversions of the definitions above.  In particular,
   u is a displacement in the deformation map x |-> x+u, and rhoBr is mass per
   D-dimensional brane volume.  The three quadratic terms must then all have
   the energy-density dimension. *)
spatialCoordinateDimension = UnitVector[3, 1];
timeCoordinateDimension = UnitVector[3, 2];
massDimension = UnitVector[3, 3];
braneMeasureDimension = D spatialCoordinateDimension;
braneMassDensityDimension = massDimension - braneMeasureDimension;
coordinateSpeedDimension =
  spatialCoordinateDimension - timeCoordinateDimension;
lagrangianTermDimensions = {
  braneMassDensityDimension + 2 timeDerivativeDimension,
  dimMuR + 2 spatialDerivativeDimension,
  dimBcomp + 2 spatialDerivativeDimension
};
dimensionChecks = {
  TrueQ[Simplify[
    energyDimension == massDimension + 2 coordinateSpeedDimension
  ]],
  TrueQ[Simplify[
    energyDensityDimension == energyDimension - braneMeasureDimension
  ]],
  TrueQ[Simplify[displacementDimension == spatialCoordinateDimension]],
  TrueQ[Simplify[dimRhoBr == braneMassDensityDimension]],
  And @@ (TrueQ[Simplify[# == energyDensityDimension]] & /@
    lagrangianTermDimensions),
  TrueQ[Simplify[dimLongitudinalSpeed == coordinateSpeedDimension]]
};

(* Bulk matching, first for the parallel brane mode. *)
longitudinalPhaseSpeedSquared = FullSimplify[
  originalParallelRoot/kSq, positiveAssumptions
];
kwSquaredSolution = FullSimplify[
  kwSquared /. First@Solve[
    longitudinalPhaseSpeedSquared kSq == cs0^2 (kSq + kwSquared),
    kwSquared
  ],
  positiveAssumptions
];
boundCondition = FullSimplify[
  Reduce[kwSquaredSolution < 0 && positiveAssumptions,
    Bcomp, Reals],
  positiveAssumptions
];
(* Strip premises already declared positive, retaining the requested threshold. *)
boundConditionOutput = Bcomp < rhoBr cs0^2;
bulkResidual = FullSimplify[
  longitudinalPhaseSpeedSquared kSq - cs0^2 (kSq + kwSquaredSolution),
  positiveAssumptions
];

transversePhaseSpeedSquared = FullSimplify[
  originalPerpendicularRoot/kSq, positiveAssumptions
];
transverseFormalKwSquared = FullSimplify[
  kwSquared /. First@Solve[
    transversePhaseSpeedSquared kSq == cs0^2 (kSq + kwSquared),
    kwSquared
  ],
  positiveAssumptions
];
transverseKernel = (SelectFirst[
  originalSpectra[[2, 2]], #[[3]] === PERPENDICULAR &
])[[4]];
transverseScalarCouplingResiduals = FullSimplify[
  alignedDirection[3] . # & /@ transverseKernel,
  positiveAssumptions
];
transverseMatchingOutput = {
  formalDispersionKwSquared -> transverseFormalKwSquared,
  selectedPerpendicularKernel -> transverseKernel,
  alignedDirectionProjectionResiduals -> transverseScalarCouplingResiduals,
  scope -> "The residuals only restate the perpendicularity used to select this brane eigenspace; the algebra does not determine whether a physical transverse bulk channel couples.",
  absentPhysicalInput -> "No brane-bulk interaction operator or boundary condition specifying coupling to bulk polarizations was supplied."
};

(* FORM control: keep Curl2 and replace Bcomp Div[u]^2/2 by
   muBr times the norm squared of the symmetric-traceless gradient. *)
formControlMatrixFromEnergy[d_Integer, kVector_List] := Module[
  {amplitude, gradient, curlSquared, symmetricTraceless, potential},
  amplitude = Array[formAmplitudeComponent, d];
  gradient = Outer[Times, kVector, amplitude];
  curlSquared = Sum[
    (gradient[[i, j]] - gradient[[j, i]])^2,
    {i, 1, d}, {j, 1, d}
  ]/2;
  symmetricTraceless = (gradient + Transpose[gradient])/2 -
    Tr[gradient] IdentityMatrix[d]/d;
  potential = muR curlSquared/2 +
    muBr Total[Flatten[symmetricTraceless]^2];
  Table[
    D[potential, amplitude[[i]], amplitude[[j]]],
    {i, 1, d}, {j, 1, d}
  ] // Simplify
];

formControlMatrixD3 = formControlMatrixFromEnergy[
  3, Sqrt[kSq] alignedDirection[3]
];
formControlSpectrumD3 = spectrumFromMatrix[
  formControlMatrixD3, alignedDirection[3],
  positiveAssumptions && 4 muBr/3 != muR + muBr
];
formPerpendicularRoot = rootWithOrientation[formControlSpectrumD3, PERPENDICULAR];
formParallelRoot = rootWithOrientation[formControlSpectrumD3, PARALLEL];
formPerpendicularDependsOnMuBr = ! TrueQ[
  FullSimplify[D[formPerpendicularRoot, muBr], positiveAssumptions] == 0
];
formParallelDependsOnMuR = ! TrueQ[
  FullSimplify[D[formParallelRoot, muR], positiveAssumptions] == 0
];
formPerpendicularMoved = ! TrueQ@FullSimplify[
  formPerpendicularRoot == originalPerpendicularRoot, positiveAssumptions
];
formParallelMoved = ! TrueQ@FullSimplify[
  formParallelRoot == originalParallelRoot, positiveAssumptions
];
formMovedSectors = Pick[
  {PERPENDICULAR, PARALLEL},
  {formPerpendicularMoved, formParallelMoved}
];
formControlOutput = {
  ROOTS -> spectrumRecords[formControlSpectrumD3],
  perpendicularRootDependsOnMuBr -> formPerpendicularDependsOnMuBr,
  parallelRootDependsOnMuR -> formParallelDependsOnMuR,
  movedSectors -> formMovedSectors
};

(* COEFFICIENT control. *)
coefficientControlMatrixD3 = alignedOriginalMatrix[3] /. Bcomp -> 2 Bcomp;
coefficientControlSpectrumD3 = spectrumFromMatrix[
  coefficientControlMatrixD3, alignedDirection[3],
  positiveAssumptions && 2 Bcomp != muR
];
coefficientPerpendicularRoot = rootWithOrientation[
  coefficientControlSpectrumD3, PERPENDICULAR
];
coefficientParallelRoot = rootWithOrientation[
  coefficientControlSpectrumD3, PARALLEL
];
coefficientPerpendicularMoved = ! TrueQ@FullSimplify[
  coefficientPerpendicularRoot == originalPerpendicularRoot,
  positiveAssumptions
];
coefficientParallelMoved = ! TrueQ@FullSimplify[
  coefficientParallelRoot == originalParallelRoot,
  positiveAssumptions
];
coefficientMovedSectors = Pick[
  {PERPENDICULAR, PARALLEL},
  {coefficientPerpendicularMoved, coefficientParallelMoved}
];
coefficientControlOutput = {
  ROOTS -> spectrumRecords[coefficientControlSpectrumD3],
  movedSectors -> coefficientMovedSectors,
  unmovedSectors -> Complement[
    {PERPENDICULAR, PARALLEL}, coefficientMovedSectors
  ]
};

(* Independent internal consistency tests. *)
spectrumChecksByDimension = MapThread[
  Function[{spectrum, d},
    {
      Total[#[[2]] & /@ spectrum[[2]]] == d,
      Sort[#[[3]] & /@ spectrum[[2]]] ===
        Sort[{PERPENDICULAR, PARALLEL}],
      And @@ (TrueQ[#[[5]]] & /@ spectrum[[2]])
    }
  ],
  {originalSpectra, {2, 3, 4, 5}}
];
spectrumChecks = And @@ Flatten[spectrumChecksByDimension];

degeneracyCheck = TrueQ@FullSimplify[
  (originalParallelRoot == originalPerpendicularRoot) /. Bcomp -> muR,
  positiveAssumptions
];

controlNullityChecks =
  Total[#[[2]] & /@ formControlSpectrumD3[[2]]] == 3 &&
  Total[#[[2]] & /@ coefficientControlSpectrumD3[[2]]] == 3;

namedChecks = {
  pieceDecomposition -> And @@ pieceConstructionChecks,
  dynamicalMatrixDerivation -> And @@ dynamicalMatrixDerivationChecks,
  spectrumNullitiesAndRoots -> spectrumChecks,
  crossSectorResiduals ->
    TrueQ[perpendicularBcompResidual == 0] && TrueQ[parallelMuRResidual == 0],
  degeneracySubstitution -> degeneracyCheck,
  dimensionalClosure -> And @@ dimensionChecks,
  bulkDispersionResidual -> TrueQ[bulkResidual == 0],
  boundThresholdReduction -> TrueQ[boundCondition === boundConditionOutput],
  controlNullities -> controlNullityChecks,
  formParallelIndependence -> ! formParallelDependsOnMuR,
  coefficientPerpendicularUnmoved -> ! coefficientPerpendicularMoved,
  coefficientParallelMovedCheck -> coefficientParallelMoved
};
allChecks = And @@ (TrueQ[Last[#]] & /@ namedChecks);
failedChecks = First /@ Select[namedChecks, ! TrueQ[Last[#]] &];

printTag["WL_S11_INVARIANT_COUNT", invariantCounts];
printTag["WL_S11_DYNAMICAL_MATRIX", dynamicalMatrixOutput];
Do[
  printTag[
    "WL_S11_SPECTRUM_D" <> ToString[d],
    spectrumRecords[originalSpectra[[index]]]
  ];
  printTag[
    "WL_S11_NULLITY_SUM_D" <> ToString[d],
    Total[#[[2]] & /@ originalSpectra[[index, 2]]]
  ],
  {index, 1, 4}, {d, {2, 3, 4, 5}[[{index}]]}
];
printBooleanTag[
  "WL_S11_PERP_ROOT_DEPENDS_ON_BCOMP", perpendicularDependsOnBcomp
];
printBooleanTag[
  "WL_S11_PARA_ROOT_DEPENDS_ON_MUR", parallelDependsOnMuR
];
printTag["WL_S11_DEGENERACY_CONDITION", degeneracyCondition];
printTag["WL_S11_DIM_RHOBR", dimensionOutput[dimRhoBr]];
printTag["WL_S11_DIM_MUR", dimensionOutput[dimMuR]];
printTag["WL_S11_DIM_BCOMP", dimensionOutput[dimBcomp]];
printTag[
  "WL_S11_DIM_SPEED_LONGITUDINAL",
  dimensionOutput[dimLongitudinalSpeed]
];
printTag["WL_S11_KW_SQUARED", kwSquaredSolution];
printTag["WL_S11_BOUND_CONDITION", boundConditionOutput];
printTag["WL_S11_TRANSVERSE_MATCHING", transverseMatchingOutput];
printTag["WL_S11_FORM_CONTROL_D3", formControlOutput];
printTag["WL_S11_COEFFICIENT_CONTROL_D3", coefficientControlOutput];
Print[
  "WL_S11_VERDICT: " <>
    If[TrueQ[allChecks], "PASS", "FAIL " <> ToString[failedChecks, InputForm]]
];

If[! TrueQ[allChecks], Exit[1]];
