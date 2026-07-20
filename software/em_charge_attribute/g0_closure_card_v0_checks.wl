(*
  Checks the scalar-block admissibility, localized transverse zero mode and
  factorization, and dimensional homogeneity of the active action/interface
  terms in g0_closure_card_v0.md.

  This is the independent Mathematica dual-engine check for
  g0_closure_card_v0.md.  Every operator, coefficient, field dimension, and
  target density below is derived from the card itself, not from the SymPy
  checker.
*)

ClearAll["Global`*"];

format[expression_] := ToString[expression, InputForm];
status[passed_] := If[TrueQ[passed], "PASS", "FAIL"];

(* ---------------------------------------------------------------------- *)
(* 1. Scalar block from Q_s(omega,k) in the card.                         *)
(* ---------------------------------------------------------------------- *)

aEff = 1;
bEff = 2;
mH = 1;
kH = 1;
cHU = 1/2;

inertiaMatrix = {{aEff, 0}, {0, mH}};
stiffnessMatrix = {{bEff, cHU}, {cHU, kH}};
stabilityMargin = FullSimplify[bEff kH - cHU^2];
stiffnessEigenvalues = FullSimplify[Eigenvalues[stiffnessMatrix]];
stiffnessEigenvalues = Sort[stiffnessEigenvalues];

scalarWellPosed = TrueQ[aEff > 0] && TrueQ[mH > 0] &&
  TrueQ[bEff > 0] && TrueQ[stabilityMargin > 0] &&
  And @@ (TrueQ[# > 0] & /@ stiffnessEigenvalues);

detPolynomial = Collect[
  Expand[Det[stiffnessMatrix - z inertiaMatrix]], z
];
waveSpeedsSquared = cSquared /. Solve[
  (detPolynomial /. z -> cSquared) == 0,
  cSquared
];
waveSpeedsSquared = Sort[FullSimplify[waveSpeedsSquared]];
expectedDetPolynomial = (4 z^2 - 12 z + 7)/4;
expectedWaveSpeedsSquared = Sort[{(3 - Sqrt[2])/2, (3 + Sqrt[2])/2}];

waveSpeedCheck = TrueQ[
    FullSimplify[detPolynomial == expectedDetPolynomial]
  ] && And @@ MapThread[
    TrueQ[FullSimplify[#1 == #2]] &,
    {waveSpeedsSquared, expectedWaveSpeedsSquared}
  ] && And @@ (TrueQ[# > 0] & /@ waveSpeedsSquared);

scalarSymbolic = <|
  "M" -> inertiaMatrix,
  "K" -> stiffnessMatrix,
  "M_h" -> mH,
  "stability margin" -> stabilityMargin,
  "eig(K)" -> stiffnessEigenvalues
|>;
scalarNumeric = N[scalarSymbolic, 12];

speedSymbolic = <|
  "det(K-z M)" -> detPolynomial,
  "c^2" -> waveSpeedsSquared
|>;
speedNumeric = N[speedSymbolic, 12];

(* ---------------------------------------------------------------------- *)
(* 2. Transverse Poschl-Teller operator from section 2.3 of the card.     *)
(* ---------------------------------------------------------------------- *)

transverseAssumptions = ell > 0 && Element[{ell, w}, Reals];
transversePotential = (4 - 6 Sech[w/ell]^2)/ell^2;
zeroMode = Sech[w/ell]^2/ell;
superpotential = 2 Tanh[w/ell]/ell;

transverseOperator[expression_] := -D[expression, {w, 2}] +
  transversePotential expression;

zeroModeResidual = FullSimplify[
  transverseOperator[zeroMode],
  Assumptions -> transverseAssumptions
];

(* For A = d/dw + W and A^dagger = -d/dw + W, the potential in
   A^dagger A is W^2-W'. *)
factorizedPotential = FullSimplify[
  superpotential^2 - D[superpotential, w],
  Assumptions -> transverseAssumptions
];
factorizationResidual = FullSimplify[
  transversePotential - factorizedPotential,
  Assumptions -> transverseAssumptions
];
firstOrderZeroModeResidual = FullSimplify[
  D[zeroMode, w] + superpotential zeroMode,
  Assumptions -> transverseAssumptions
];

zeroModeNorm = FullSimplify[
  Integrate[
    2 zeroMode^2,
    {w, -Infinity, Infinity},
    Assumptions -> ell > 0,
    GenerateConditions -> False
  ],
  Assumptions -> ell > 0
];

ellStar = 1/20;
transverseCheck = TrueQ[zeroModeResidual === 0] &&
  TrueQ[factorizationResidual === 0] &&
  TrueQ[firstOrderZeroModeResidual === 0] &&
  TrueQ[FullSimplify[zeroModeNorm == 8/(3 ell), ell > 0]] &&
  TrueQ[FullSimplify[zeroModeNorm > 0, ell > 0]];

transverseSymbolic = <|
  "O_perp f0" -> zeroModeResidual,
  "V-(W^2-W')" -> factorizationResidual,
  "A f0" -> firstOrderZeroModeResidual,
  "N0" -> zeroModeNorm
|>;
transverseNumeric = <|
  "O_perp f0" -> N[zeroModeResidual /. {ell -> ellStar, w -> 1/7}, 12],
  "V-(W^2-W')" -> N[factorizationResidual /. {ell -> ellStar, w -> 1/7}, 12],
  "A f0" -> N[firstOrderZeroModeResidual /. {ell -> ellStar, w -> 1/7}, 12],
  "N0 at ell*=1/20" -> N[zeroModeNorm /. ell -> ellStar, 12]
|>;

(* ---------------------------------------------------------------------- *)
(* 3. (L,T,M) dimensions, derived from the card's field/coefficient       *)
(*    tables.  Bulk spatial densities target E/L^4; brane spatial         *)
(*    densities target E/L^3.                                             *)
(* ---------------------------------------------------------------------- *)

dimZero = {0, 0, 0};
dimL = {1, 0, 0};
dimT = {0, 1, 0};
dimMass = {0, 0, 1};
dimEnergy = {2, -2, 1};
dimAction = dimEnergy + dimT;
dimVelocity = dimL - dimT;
dimForce = dimEnergy - dimL;
dimPower = dimEnergy - dimT;
dimMomentum = dimMass + dimVelocity;

bulkDensityTarget = dimEnergy - 4 dimL;
braneDensityTarget = dimEnergy - 3 dimL;
twoSurfaceDensityTarget = dimEnergy - 2 dimL;
massSourceTarget = -4 dimL - dimT;
momentumSourceTarget = dimMomentum - 4 dimL - dimT;
energySourceTarget = dimEnergy - 4 dimL - dimT;
energyCurrentTarget = dimEnergy - 3 dimL - dimT;

dimRho = -4 dimL;
dimTheta = dimZero;
dimRB = dimZero;
dimParentH = -dimL;
dimReducedH = dimZero;
dimUL = dimL;
dimHbar = dimAction;
dimParticleMass = dimMass;
dimBulkK = bulkDensityTarget - 5 dimRho;

dimZChi = {-2, 0, 1};
dimKappaChi = {0, -2, 1};
dimLambdaChi = {-2, -2, 1};
dimLambdaSigma = {-1, -2, 1};
dimM4 = {0, 0, 1};
dimK4 = {2, -2, 1};
dimAEff = {-3, 0, 1};
dimBEff = {-1, -2, 1};
dimMH = {-1, 0, 1};
dimKH = {1, -2, 1};
dimCHU = {0, -2, 1};
dimEta = -3 dimL;
dimOddKernel = -2 dimL;
dimKm = {4, -2, 1};
dimJm = {3, -2, 1};
dimReducedKm = dimEnergy;
dimReducedG = dimEnergy;
dimGamma0 = -dimT;
dimD = -4 dimL;
dimReturnKernel = -4 dimL;
dimN0 = -dimL;
dimNChi = -dimL;

dimDt[fieldDimension_] := fieldDimension - dimT;
dimGrad[fieldDimension_] := fieldDimension - dimL;
dimLaplacian[fieldDimension_] := fieldDimension - 2 dimL;
dimCheck[name_, actual_, expected_] := name -> FullSimplify[actual - expected];

dimensionChecks = {
  (* Four-spatial-dimensional bulk action density S_bulk. *)
  dimCheck["bulk: hbar rho dt(theta)",
    dimHbar + dimRho + dimDt[dimTheta], bulkDensityTarget],
  dimCheck["bulk: (hbar^2 rho/m)|grad(theta)|^2",
    2 dimHbar + dimRho - dimParticleMass + 2 dimGrad[dimTheta],
    bulkDensityTarget],
  dimCheck["bulk: (K/4) rho^5",
    dimBulkK + 5 dimRho, bulkDensityTarget],
  dimCheck["bulk: hbar^2 |grad(rho)|^2/(m rho)",
    2 dimHbar - dimParticleMass - dimRho + 2 dimGrad[dimRho],
    bulkDensityTarget],

  (* Four-spatial-dimensional wall action density S_chi. *)
  dimCheck["bulk wall: Z_chi (dt r_B)^2",
    dimZChi + 2 dimDt[dimRB], bulkDensityTarget],
  dimCheck["bulk wall: kappa_chi |grad r_B|^2",
    dimKappaChi + 2 dimGrad[dimRB], bulkDensityTarget],
  dimCheck["bulk wall: lambda_chi r_B^2(1-r_B)^2",
    dimLambdaChi, bulkDensityTarget],

  (* Four-spatial-dimensional parent-H action density S_H. *)
  dimCheck["bulk parent H: M4 (dt H)^2",
    dimM4 + 2 dimDt[dimParentH], bulkDensityTarget],
  dimCheck["bulk parent H: K4 |grad_x H|^2",
    dimK4 + 2 dimGrad[dimParentH], bulkDensityTarget],
  dimCheck["bulk parent H: K4 H O_perp H",
    dimK4 + dimParentH + dimLaplacian[dimParentH],
    bulkDensityTarget],
  dimCheck["zero-mode norm: integral dw 2 f0^2",
    dimL + 2 dimParentH, dimN0],
  dimCheck["zero-mode reduction: N0 M4=M_h",
    dimN0 + dimM4, dimMH],
  dimCheck["zero-mode reduction: N0 K4=K_h",
    dimN0 + dimK4, dimKH],

  (* Three-spatial-dimensional reduced brane action density S_Lh. *)
  dimCheck["brane: A_eff (dt u_L)^2",
    dimAEff + 2 dimDt[dimUL], braneDensityTarget],
  dimCheck["brane: B_eff |grad u_L|^2",
    dimBEff + 2 dimGrad[dimUL], braneDensityTarget],
  dimCheck["brane: M_h (dt h)^2",
    dimMH + 2 dimDt[dimReducedH], braneDensityTarget],
  dimCheck["brane: K_h |grad h|^2",
    dimKH + 2 dimGrad[dimReducedH], braneDensityTarget],
  dimCheck["brane: C_hu grad(u_L).grad(h)",
    dimCHU + dimGrad[dimUL] + dimGrad[dimReducedH],
    braneDensityTarget],

  (* Frozen-wall holding and mouth interface functionals. *)
  dimCheck["holding surface density: lambda_Sigma(r_B-1/2)",
    dimLambdaSigma, braneDensityTarget],
  dimCheck["holding functional after d^3A",
    3 dimL + dimLambdaSigma, dimEnergy],
  dimCheck["wall conormal: kappa_chi partial_n r_B=lambda_Sigma",
    dimKappaChi + dimGrad[dimRB], dimLambdaSigma],
  dimCheck["orientation N_chi integral",
    3 dimL + dimEta + dimL + dimOddKernel, dimNChi],
  dimCheck["orientation Q_chi numerator/N_chi",
    3 dimL + dimEta + dimL + dimOddKernel - dimNChi, dimZero],
  dimCheck["mouth Robin term after d^3x eta",
    3 dimL + dimEta + dimKm + 2 dimParentH, dimEnergy],
  dimCheck["mouth source term after d^3x eta",
    3 dimL + dimEta + dimJm + dimParentH, dimEnergy],
  dimCheck["parent mouth jump: K4 partial_w H",
    dimK4 + dimGrad[dimParentH], twoSurfaceDensityTarget],
  dimCheck["parent mouth jump: eta K_m H",
    dimEta + dimKm + dimParentH, twoSurfaceDensityTarget],
  dimCheck["parent mouth jump: eta J_m Q_chi",
    dimEta + dimJm, twoSurfaceDensityTarget],
  dimCheck["reduced mouth h source: eta k_m h",
    dimEta + dimReducedKm + dimReducedH, braneDensityTarget],
  dimCheck["reduced mouth fixed source: eta g_chih s_i",
    dimEta + dimReducedG, braneDensityTarget],

  (* Reduced Euler operators and their natural conormals. *)
  dimCheck["u equation: B_eff laplacian(u_L)",
    dimBEff + dimLaplacian[dimUL], bulkDensityTarget],
  dimCheck["u equation: C_hu laplacian(h)",
    dimCHU + dimLaplacian[dimReducedH], bulkDensityTarget],
  dimCheck["h equation: K_h laplacian(h)",
    dimKH + dimLaplacian[dimReducedH], braneDensityTarget],
  dimCheck["h equation: C_hu laplacian(u_L)",
    dimCHU + dimLaplacian[dimUL], braneDensityTarget],
  dimCheck["u conormal: B_eff partial_n u_L",
    dimBEff + dimGrad[dimUL], braneDensityTarget],
  dimCheck["u conormal: C_hu partial_n h",
    dimCHU + dimGrad[dimReducedH], braneDensityTarget],
  dimCheck["h conormal: K_h partial_n h",
    dimKH + dimGrad[dimReducedH], twoSurfaceDensityTarget],
  dimCheck["h conormal: C_hu partial_n u_L",
    dimCHU + dimGrad[dimUL], twoSurfaceDensityTarget],

  (* The explicit E2 surface-variation functional in section 8.2. *)
  dimCheck["E2 interface: m rho (v-V).n delta(Phi)",
    dimParticleMass + dimRho + dimVelocity + (2 dimL - dimT),
    braneDensityTarget],
  dimCheck["E2 interface: hbar^2 (n.grad rho) delta(rho)/(m rho)",
    2 dimHbar - dimParticleMass - dimRho + dimGrad[dimRho] + dimRho,
    braneDensityTarget],
  dimCheck["E2 interface: (Pi n).delta(xi)",
    bulkDensityTarget + dimL, braneDensityTarget],
  dimCheck["E2/E3 normal mass flux rho(v-V).n",
    dimRho + dimVelocity, -3 dimL - dimT],
  dimCheck["E2 mouth integrated mass flux",
    3 dimL + dimRho + dimVelocity, dimGamma0],

  (* Bulk stress/energy package and non-Hamiltonian source functionals. *)
  dimCheck["quantum stress Pi_Q",
    2 dimHbar - dimParticleMass + dimRho - 2 dimL,
    bulkDensityTarget],
  dimCheck["quantum energy density epsilon_Q",
    2 dimHbar - dimParticleMass - dimRho + 2 dimGrad[dimRho],
    bulkDensityTarget],
  dimCheck["quantum energy current: (rho Q_B-epsilon_Q) v",
    bulkDensityTarget + dimVelocity, energyCurrentTarget],
  dimCheck["quantum energy current: hbar^2 (dt rho/rho) grad rho/m",
    2 dimHbar - dimParticleMass - dimT + dimGrad[dimRho],
    energyCurrentTarget],
  dimCheck["normalized drain/return kernel D_i or R_0",
    dimD, -4 dimL],
  dimCheck["mass source Gamma_0 D_i",
    dimGamma0 + dimD, dimDt[dimRho]],
  dimCheck["comoving momentum source m S v",
    dimParticleMass + massSourceTarget + dimVelocity,
    momentumSourceTarget],
  dimCheck["return momentum correction R_0 I_ret",
    dimReturnKernel + dimForce, momentumSourceTarget],
  dimCheck["carrier enthalpy per particle e_c",
    bulkDensityTarget - dimRho, dimEnergy],
  dimCheck["comoving energy source e_c S",
    dimEnergy + massSourceTarget, dimDt[bulkDensityTarget]],
  dimCheck["return energy correction R_0 P_ret",
    dimReturnKernel + dimPower, dimDt[bulkDensityTarget]],
  dimCheck["continuity lhs partial_t rho",
    dimDt[dimRho], massSourceTarget],
  dimCheck["momentum balance lhs partial_t(m rho v)",
    dimParticleMass + dimRho + dimVelocity - dimT,
    momentumSourceTarget],
  dimCheck["energy balance lhs partial_t epsilon",
    bulkDensityTarget - dimT, energySourceTarget],

  (* Separation-independent geon term. *)
  dimCheck["geon constant time integrand E_g", dimEnergy, dimEnergy],
  dimCheck["geon constant action after dt", dimT + dimEnergy, dimAction]
};

dimensionResiduals = Association[dimensionChecks];
nonzeroDimensionResiduals = Select[
  dimensionResiduals,
  ! TrueQ[# === dimZero] &
];
dimensionCheck = Length[nonzeroDimensionResiduals] == 0;
dimensionNumeric = Map[N[#, 6] &, dimensionResiduals];

(* ---------------------------------------------------------------------- *)
(* 4. Card-level static prerequisites and mutation/ablation harness.      *)
(*                                                                        *)
(* The existing headline calculations above are intentionally retained.  *)
(* This block adds one predicate per atomic card-level claim.  Every      *)
(* predicate is recomputed from a primitive association so a mutation of  *)
(* a card primitive can exercise the actual derivation.                   *)
(* ---------------------------------------------------------------------- *)

basePrimitives = <|
  "AEff" -> 1,
  "BEff" -> 2,
  "MHDeclared" -> 1,
  "KH" -> 1,
  "CHU" -> 1/2,
  "ZeroModePower" -> 2,
  "TransverseSuperpotentialCoefficient" -> 2,
  "WallLambdaRatio" -> 2,
  "BulkUPower" -> 5,
  "PhasePinningCoefficient" -> 0,
  "PhaseQuotientWeightIntegral" -> 1,
  "BumpNormalizationMultiplier" -> 1,
  "Gamma0" -> 1/1000,
  "MassReturnMultiplier" -> 1,
  "MomentumControllerSign" -> -1,
  "EnergyControllerSign" -> -1,
  "M4Multiplier" -> 1,
  "ReducedHMassCoefficient" -> 0,
  "JmOverEll" -> 1,
  "GChiH" -> 1,
  "UEquationCrossCoefficient" -> 1/2,
  "HEquationCrossCoefficient" -> 1/2,
  "DimLambdaChi" -> dimLambdaChi
|>;

checkOrder = {
  "SCALAR_WELL_POSEDNESS",
  "GENERALIZED_WAVE_SPEEDS",
  "PARENT_TRANSVERSE_ZERO_EIGENVALUE",
  "TRANSVERSE_FACTORIZATION",
  "PLANAR_WALL_FACTORIZATION",
  "POSITIVE_BULK_HESSIAN",
  "PHASE_CONSTANT_MODE_QUOTIENT",
  "DRAIN_RETURN_KERNEL_NORMALIZATION",
  "POSITIVE_DRAIN_THROUGHPUT",
  "DRAIN_RETURN_MASS_CLOSURE",
  "DRAIN_RETURN_MOMENTUM_CLOSURE",
  "DRAIN_RETURN_ENERGY_CLOSURE",
  "LOCALIZED_H_POSITIVE_NORM",
  "REDUCED_H_MASSLESSNESS",
  "BARE_MOUTH_MONOPOLE_PREREQUISITE",
  "CONSERVATIVE_HESSIAN_SYMMETRY_PREREQUISITE",
  "DIMENSIONAL_HOMOGENEITY"
};

mutationDescriptions = <|
  "SCALAR_WELL_POSEDNESS" -> "C_hu: 1/2 -> 2",
  "GENERALIZED_WAVE_SPEEDS" -> "A_eff: 1 -> 2",
  "PARENT_TRANSVERSE_ZERO_EIGENVALUE" ->
    "f0 profile: Sech[w/ell]^2/ell -> Sech[w/ell]/ell",
  "TRANSVERSE_FACTORIZATION" ->
    "A superpotential coefficient: 2 -> 3",
  "PLANAR_WALL_FACTORIZATION" ->
    "lambda_chi ell^2/kappa_chi: 2 -> 3",
  "POSITIVE_BULK_HESSIAN" -> "power in U(rho): 5 -> 4",
  "PHASE_CONSTANT_MODE_QUOTIENT" ->
    "bulk phase pinning coefficient: 0 -> 1",
  "DRAIN_RETURN_KERNEL_NORMALIZATION" ->
    "B_n normalization multiplier: 1 -> 2",
  "POSITIVE_DRAIN_THROUGHPUT" -> "Gamma_0: 1/1000 -> -1/1000",
  "DRAIN_RETURN_MASS_CLOSURE" ->
    "S_leakage return multiplier: 1 -> 3/2",
  "DRAIN_RETURN_MOMENTUM_CLOSURE" ->
    "I_ret defining sign: -1 -> +1",
  "DRAIN_RETURN_ENERGY_CLOSURE" ->
    "P_ret defining sign: -1 -> +1",
  "LOCALIZED_H_POSITIVE_NORM" ->
    "M4 normalization multiplier: 1 -> 1/2",
  "REDUCED_H_MASSLESSNESS" -> "reduced h^2 coefficient: 0 -> 1",
  "BARE_MOUTH_MONOPOLE_PREREQUISITE" ->
    "parent J_m/ell: 1 -> 0 (reduced g_chih held at 1)",
  "CONSERVATIVE_HESSIAN_SYMMETRY_PREREQUISITE" ->
    "u_L-Euler cross coefficient: C_hu -> 3/4",
  "DIMENSIONAL_HOMOGENEITY" ->
    "dimension of lambda_chi: (-2,-2,1) -> (-1,-2,1)"
|>;

mutatePrimitives[checkID_String] := Module[{p = Association[basePrimitives]},
  Switch[checkID,
    "SCALAR_WELL_POSEDNESS", p["CHU"] = 2,
    "GENERALIZED_WAVE_SPEEDS", p["AEff"] = 2,
    "PARENT_TRANSVERSE_ZERO_EIGENVALUE", p["ZeroModePower"] = 1,
    "TRANSVERSE_FACTORIZATION",
      p["TransverseSuperpotentialCoefficient"] = 3,
    "PLANAR_WALL_FACTORIZATION", p["WallLambdaRatio"] = 3,
    "POSITIVE_BULK_HESSIAN", p["BulkUPower"] = 4,
    "PHASE_CONSTANT_MODE_QUOTIENT", p["PhasePinningCoefficient"] = 1,
    "DRAIN_RETURN_KERNEL_NORMALIZATION",
      p["BumpNormalizationMultiplier"] = 2,
    "POSITIVE_DRAIN_THROUGHPUT", p["Gamma0"] = -1/1000,
    "DRAIN_RETURN_MASS_CLOSURE", p["MassReturnMultiplier"] = 3/2,
    "DRAIN_RETURN_MOMENTUM_CLOSURE", p["MomentumControllerSign"] = 1,
    "DRAIN_RETURN_ENERGY_CLOSURE", p["EnergyControllerSign"] = 1,
    "LOCALIZED_H_POSITIVE_NORM", p["M4Multiplier"] = 1/2,
    "REDUCED_H_MASSLESSNESS", p["ReducedHMassCoefficient"] = 1,
    "BARE_MOUTH_MONOPOLE_PREREQUISITE", p["JmOverEll"] = 0,
    "CONSERVATIVE_HESSIAN_SYMMETRY_PREREQUISITE",
      p["UEquationCrossCoefficient"] = 3/4,
    "DIMENSIONAL_HOMOGENEITY", p["DimLambdaChi"] = {-1, -2, 1},
    _, Null
  ];
  p
];

evaluatePredicates[p_Association] := Module[
  {
    localM, localK, localMargin, localEigenvalues, localDet, localSpeeds,
    localExpectedSpeeds, f0Local, parentPotential, parentResidual,
    transverseW, transverseFactorPotential, transverseFactorResidual,
    r0, r0Prime, wallV, wallVSecond, wallPotentialPerKappa, wallW,
    wallFactorPotential, wallFactorResidual, wallFirstOrderResidual,
    wallZeroResidual, rho0, particleMass, soundSpeed, xiQ, bulkKDerived,
    hbarDerived, muDerived, bulkU, bulkUSecond, fourierHessian,
    fourierConstant, fourierK2, fourierPositive, phaseConstantEigenvalue,
    phaseQuotientResidual, bumpNorm, bump1Norm, bump3Norm, drainNorm,
    returnShellNorm, returnNorm, gamma, massResidual, vD1, vD2, vR,
    momentumRaw, iRet, momentumResidual, eD1, eD2, eR, energyRaw,
    pRet, energyResidual, n0Local, m4Local, reducedMass,
    qScalar, qAtOrigin, annulusNorm, jmLocal, gChiH, projectedSources,
    conservativeDensity, conservativeHessian, eulerStiffness,
    conservativeResidual, mutatedDimensionResiduals,
    mutatedNonzeroDimensions, results, values
  },

  (* Scalar matrices are independently formed from the Q_s coefficients. *)
  localM = {{p["AEff"], 0}, {0, p["MHDeclared"]}};
  localK = {{p["BEff"], p["CHU"]}, {p["CHU"], p["KH"]}};
  localMargin = FullSimplify[p["BEff"] p["KH"] - p["CHU"]^2];
  localEigenvalues = Sort[FullSimplify[Eigenvalues[localK]]];
  localDet = Collect[Expand[Det[localK - z localM]], z];
  localSpeeds = Sort[FullSimplify[
    cSquared /. Solve[(localDet /. z -> cSquared) == 0, cSquared]
  ]];
  localExpectedSpeeds = Sort[{(3 - Sqrt[2])/2, (3 + Sqrt[2])/2}];

  (* Parent O_perp and f0 are separately taken from section 2.3. *)
  f0Local = Sech[w/ell]^p["ZeroModePower"]/ell;
  parentPotential = (4 - 6 Sech[w/ell]^2)/ell^2;
  parentResidual = FullSimplify[
    -D[f0Local, {w, 2}] + parentPotential f0Local,
    Assumptions -> ell > 0 && Element[w, Reals]
  ];
  transverseW = p["TransverseSuperpotentialCoefficient"] Tanh[w/ell]/ell;
  transverseFactorPotential = FullSimplify[
    transverseW^2 - D[transverseW, w],
    Assumptions -> ell > 0 && Element[w, Reals]
  ];
  transverseFactorResidual = FullSimplify[
    parentPotential - transverseFactorPotential,
    Assumptions -> ell > 0 && Element[w, Reals]
  ];

  (* V_chi'' supplies L_chi^(2); the kink derivative supplies A_chi. *)
  r0 = 1/(1 + Exp[-x/ell]);
  r0Prime = D[r0, x];
  wallV = (p["WallLambdaRatio"]/(4 ell^2))
    rWall^2 (1 - rWall)^2;
  wallVSecond = D[wallV, {rWall, 2}];
  wallPotentialPerKappa = FullSimplify[
    wallVSecond /. rWall -> r0,
    Assumptions -> ell > 0 && Element[x, Reals]
  ];
  wallW = FullSimplify[
    -D[r0Prime, x]/r0Prime,
    Assumptions -> ell > 0 && Element[x, Reals]
  ];
  wallFactorPotential = FullSimplify[
    wallW^2 - D[wallW, x],
    Assumptions -> ell > 0 && Element[x, Reals]
  ];
  wallFactorResidual = FullSimplify[
    wallPotentialPerKappa - wallFactorPotential,
    Assumptions -> ell > 0 && Element[x, Reals]
  ];
  wallFirstOrderResidual = FullSimplify[
    D[r0Prime, x] + wallW r0Prime,
    Assumptions -> ell > 0 && Element[x, Reals]
  ];
  wallZeroResidual = FullSimplify[
    -D[r0Prime, {x, 2}] + wallPotentialPerKappa r0Prime,
    Assumptions -> ell > 0 && Element[x, Reals]
  ];

  (* Section 7 witness primitives determine K, hbar, and mu_infinity. *)
  rho0 = 1;
  particleMass = 1;
  soundSpeed = 2;
  xiQ = 1/20;
  bulkKDerived = particleMass soundSpeed^2/(5 rho0^4);
  hbarDerived = Sqrt[2] particleMass soundSpeed xiQ;
  muDerived = (5 bulkKDerived/4) rho0^4;
  bulkU = (bulkKDerived/4) rho^p["BulkUPower"];
  bulkUSecond = FullSimplify[D[bulkU, {rho, 2}] /. rho -> rho0];
  fourierHessian = FullSimplify[
    bulkUSecond + hbarDerived^2 k^2/(4 particleMass rho0)
  ];
  fourierConstant = Coefficient[fourierHessian, k, 0];
  fourierK2 = Coefficient[fourierHessian, k, 2];
  fourierPositive = TrueQ[Reduce[fourierHessian <= 0, k, Reals] === False];

  (* The bulk phase has only its gradient operator and a normalized quotient. *)
  phaseConstantEigenvalue = p["PhasePinningCoefficient"];
  phaseQuotientResidual = p["PhaseQuotientWeightIntegral"] - 1;

  (* Independently integrate the actual B_n and return-shell forms. *)
  bumpNorm[n_Integer] := FullSimplify[
    p["BumpNormalizationMultiplier"]
      Gamma[n/2 + 3]/(2 Pi^(n/2) sigma^n)
      (2 Pi^(n/2)/Gamma[n/2])
      Integrate[
        q^(n - 1) (1 - q^2/sigma^2)^2,
        {q, 0, sigma},
        Assumptions -> sigma > 0,
        GenerateConditions -> False
      ],
    Assumptions -> sigma > 0
  ];
  bump1Norm = bumpNorm[1];
  bump3Norm = bumpNorm[3];
  drainNorm = FullSimplify[bump3Norm bump1Norm];
  returnShellNorm = FullSimplify[
    4 Pi (3/(4 Pi ((rRet + a)^3 - rRet^3)))
      Integrate[q^2, {q, rRet, rRet + a},
        Assumptions -> a > 0 && rRet > 0],
    Assumptions -> a > 0 && rRet > 0
  ];
  returnNorm = FullSimplify[
    returnShellNorm (bump1Norm + bump1Norm)/2
  ];

  (* Independent moments represent arbitrary spatially-varying v(y),e_c(y). *)
  gamma = p["Gamma0"];
  massResidual = FullSimplify[
    -2 gamma drainNorm +
      2 gamma p["MassReturnMultiplier"] returnNorm
  ];
  vD1 = Array[vD1Moment, 4];
  vD2 = Array[vD2Moment, 4];
  vR = Array[vReturnMoment, 4];
  momentumRaw = FullSimplify[
    particleMass gamma (-drainNorm vD1 - drainNorm vD2 + 2 returnNorm vR)
  ];
  iRet = p["MomentumControllerSign"] momentumRaw;
  momentumResidual = FullSimplify[momentumRaw + returnNorm iRet];
  eD1 = eD1Moment;
  eD2 = eD2Moment;
  eR = eReturnMoment;
  energyRaw = FullSimplify[
    gamma (-drainNorm eD1 - drainNorm eD2 + 2 returnNorm eR)
  ];
  pRet = p["EnergyControllerSign"] energyRaw;
  energyResidual = FullSimplify[energyRaw + returnNorm pRet];

  (* N0 is integrated from f0; M_h=1 is an independent card declaration. *)
  n0Local = FullSimplify[
    Integrate[2 f0Local^2, {w, -Infinity, Infinity},
      Assumptions -> ell > 0, GenerateConditions -> False],
    Assumptions -> ell > 0
  ];
  m4Local = p["M4Multiplier"] 3 ell/8;
  reducedMass = FullSimplify[n0Local m4Local, Assumptions -> ell > 0];

  (* Q_s includes an independently declared reduced h^2 coefficient. *)
  qScalar = {
    {p["AEff"] omega^2 - p["BEff"] k^2, -p["CHU"] k^2},
    {-p["CHU"] k^2,
      p["MHDeclared"] omega^2 - p["KH"] k^2 -
        p["ReducedHMassCoefficient"]}
  };
  qAtOrigin = qScalar /. {omega -> 0, k -> 0};

  (* Independently integrate the parent J_m f0(0) source and compare it
     with the separately declared reduced g_chih=J_m/ell coefficient. *)
  annulusNorm = FullSimplify[
    4 Pi (3/(4 Pi ((a + ell)^3 - a^3)))
      Integrate[q^2, {q, a, a + ell}, Assumptions -> a > 0 && ell > 0],
    Assumptions -> a > 0 && ell > 0
  ];
  jmLocal = p["JmOverEll"] ell;
  gChiH = p["GChiH"];
  projectedSources =
    FullSimplify[annulusNorm jmLocal (Sech[0]^2/ell) #] & /@ {-1, 1};

  (* Compare the action Hessian with the separately displayed Euler block. *)
  conservativeDensity = (p["BEff"]/2) gradU^2 +
    (p["KH"]/2) gradH^2 + p["CHU"] gradU gradH;
  conservativeHessian = FullSimplify[
    D[conservativeDensity, {{gradU, gradH}, 2}]
  ];
  eulerStiffness = {
    {p["BEff"], p["UEquationCrossCoefficient"]},
    {p["HEquationCrossCoefficient"], p["KH"]}
  };
  conservativeResidual = FullSimplify[conservativeHessian - eulerStiffness];

  mutatedDimensionResiduals = ReplacePart[
    dimensionResiduals,
    "bulk wall: lambda_chi r_B^2(1-r_B)^2" ->
      FullSimplify[p["DimLambdaChi"] - bulkDensityTarget]
  ];
  mutatedNonzeroDimensions = Select[
    mutatedDimensionResiduals, ! TrueQ[# === dimZero] &
  ];

  results = <|
    "SCALAR_WELL_POSEDNESS" -> (
      TrueQ[p["AEff"] > 0] && TrueQ[p["MHDeclared"] > 0] &&
      TrueQ[p["BEff"] > 0] && TrueQ[p["KH"] > 0] &&
      TrueQ[localMargin > 0] &&
      And @@ (TrueQ[# > 0] & /@ localEigenvalues)
    ),
    "GENERALIZED_WAVE_SPEEDS" -> (
      TrueQ[FullSimplify[localDet == (4 z^2 - 12 z + 7)/4]] &&
      And @@ MapThread[
        TrueQ[FullSimplify[#1 == #2]] &,
        {localSpeeds, localExpectedSpeeds}
      ] && And @@ (TrueQ[# > 0] & /@ localSpeeds)
    ),
    "PARENT_TRANSVERSE_ZERO_EIGENVALUE" -> TrueQ[parentResidual === 0],
    "TRANSVERSE_FACTORIZATION" -> TrueQ[transverseFactorResidual === 0],
    "PLANAR_WALL_FACTORIZATION" -> (
      TrueQ[FullSimplify[wallW == Tanh[x/(2 ell)]/ell,
        ell > 0 && Element[x, Reals]]] &&
      TrueQ[wallFactorResidual === 0] &&
      TrueQ[wallFirstOrderResidual === 0] &&
      TrueQ[wallZeroResidual === 0]
    ),
    "POSITIVE_BULK_HESSIAN" -> (
      TrueQ[bulkKDerived == 4/5] &&
      TrueQ[hbarDerived == Sqrt[2]/10] && TrueQ[muDerived == 1] &&
      TrueQ[bulkUSecond == 4] && TrueQ[bulkUSecond > 0] &&
      TrueQ[fourierConstant > 0] && TrueQ[fourierK2 > 0] &&
      fourierPositive
    ),
    "PHASE_CONSTANT_MODE_QUOTIENT" -> (
      TrueQ[phaseConstantEigenvalue == 0] &&
      TrueQ[phaseQuotientResidual == 0]
    ),
    "DRAIN_RETURN_KERNEL_NORMALIZATION" -> (
      TrueQ[bump1Norm == 1] && TrueQ[bump3Norm == 1] &&
      TrueQ[drainNorm == 1] && TrueQ[returnShellNorm == 1] &&
      TrueQ[returnNorm == 1]
    ),
    "POSITIVE_DRAIN_THROUGHPUT" -> TrueQ[gamma > 0],
    "DRAIN_RETURN_MASS_CLOSURE" -> TrueQ[massResidual == 0],
    "DRAIN_RETURN_MOMENTUM_CLOSURE" ->
      And @@ (TrueQ[# == 0] & /@ momentumResidual),
    "DRAIN_RETURN_ENERGY_CLOSURE" -> TrueQ[energyResidual == 0],
    "LOCALIZED_H_POSITIVE_NORM" -> (
      TrueQ[FullSimplify[n0Local == 8/(3 ell), ell > 0]] &&
      TrueQ[FullSimplify[n0Local > 0, ell > 0]] &&
      TrueQ[reducedMass == p["MHDeclared"]] && TrueQ[reducedMass > 0]
    ),
    "REDUCED_H_MASSLESSNESS" -> (
      TrueQ[p["ReducedHMassCoefficient"] == 0] &&
      TrueQ[qAtOrigin == ConstantArray[0, {2, 2}]]
    ),
    "BARE_MOUTH_MONOPOLE_PREREQUISITE" -> (
      TrueQ[annulusNorm == 1] && TrueQ[gChiH == p["JmOverEll"]] &&
      TrueQ[gChiH != 0] &&
      TrueQ[projectedSources == gChiH {-1, 1}] &&
      And @@ (TrueQ[# != 0] & /@ projectedSources)
    ),
    "CONSERVATIVE_HESSIAN_SYMMETRY_PREREQUISITE" -> (
      TrueQ[conservativeHessian == Transpose[conservativeHessian]] &&
      TrueQ[eulerStiffness == Transpose[eulerStiffness]] &&
      TrueQ[conservativeResidual == ConstantArray[0, {2, 2}]]
    ),
    "DIMENSIONAL_HOMOGENEITY" ->
      TrueQ[Length[mutatedNonzeroDimensions] == 0]
  |>;

  values = <|
    "SCALAR_WELL_POSEDNESS" -> <|
      "M" -> localM, "K" -> localK, "margin" -> localMargin,
      "eig(K)" -> localEigenvalues
    |>,
    "GENERALIZED_WAVE_SPEEDS" -> <|
      "det(K-z M)" -> localDet, "c^2" -> localSpeeds
    |>,
    "PARENT_TRANSVERSE_ZERO_EIGENVALUE" -> <|
      "O_perp f0" -> parentResidual
    |>,
    "TRANSVERSE_FACTORIZATION" -> <|
      "O_perp-(A^dagger A)" -> transverseFactorResidual
    |>,
    "PLANAR_WALL_FACTORIZATION" -> <|
      "V_chi''(r0)" -> wallPotentialPerKappa,
      "W_from_r0prime" -> wallW,
      "L/kappa-A_chi^dagger A_chi" -> wallFactorResidual,
      "A_chi r0prime" -> wallFirstOrderResidual,
      "L r0prime/kappa" -> wallZeroResidual
    |>,
    "POSITIVE_BULK_HESSIAN" -> <|
      "K" -> bulkKDerived, "hbar" -> hbarDerived,
      "mu_infinity" -> muDerived, "U''(rho0)" -> bulkUSecond,
      "Fourier Hessian" -> fourierHessian,
      "coefficients" -> {fourierConstant, fourierK2}
    |>,
    "PHASE_CONSTANT_MODE_QUOTIENT" -> <|
      "constant eigenvalue" -> phaseConstantEigenvalue,
      "integral W_IR" -> p["PhaseQuotientWeightIntegral"]
    |>,
    "DRAIN_RETURN_KERNEL_NORMALIZATION" -> <|
      "integral B1" -> bump1Norm, "integral B3" -> bump3Norm,
      "integral D_i" -> drainNorm, "integral eta_ret" -> returnShellNorm,
      "integral R_0" -> returnNorm
    |>,
    "POSITIVE_DRAIN_THROUGHPUT" -> <|"Gamma_0" -> gamma|>,
    "DRAIN_RETURN_MASS_CLOSURE" -> <|
      "integrated residual" -> massResidual
    |>,
    "DRAIN_RETURN_MOMENTUM_CLOSURE" -> <|
      "arbitrary v moments" -> {vD1, vD2, vR},
      "integrated residual" -> momentumResidual
    |>,
    "DRAIN_RETURN_ENERGY_CLOSURE" -> <|
      "arbitrary e_c moments" -> {eD1, eD2, eR},
      "integrated residual" -> energyResidual
    |>,
    "LOCALIZED_H_POSITIVE_NORM" -> <|
      "N0" -> n0Local, "M4" -> m4Local,
      "N0 M4" -> reducedMass, "declared M_h" -> p["MHDeclared"]
    |>,
    "REDUCED_H_MASSLESSNESS" -> <|
      "h^2 k^0 coefficient" -> p["ReducedHMassCoefficient"],
      "Q_s(0,0)" -> qAtOrigin
    |>,
    "BARE_MOUTH_MONOPOLE_PREREQUISITE" -> <|
      "integral eta_i" -> annulusNorm, "g_chih=J_m/ell" -> gChiH,
      "parent-projected sources s={-1,+1}" -> projectedSources,
      "assembled reduced sources s={-1,+1}" -> gChiH {-1, 1}
    |>,
    "CONSERVATIVE_HESSIAN_SYMMETRY_PREREQUISITE" -> <|
      "action Hessian" -> conservativeHessian,
      "Euler stiffness" -> eulerStiffness,
      "difference" -> conservativeResidual
    |>,
    "DIMENSIONAL_HOMOGENEITY" -> <|
      "checked" -> Length[mutatedDimensionResiduals],
      "nonzero residuals" -> mutatedNonzeroDimensions,
      "mass source compared with" -> dimDt[dimRho],
      "energy sources compared with" -> dimDt[bulkDensityTarget]
    |>
  |>;

  <|"Results" -> results, "Values" -> values|>
];

dispositionRows = {
  {"10.2-1", "localized H: normalizable N0>0 and physical M4 N0>0",
    "(1)", "LOCALIZED_H_POSITIVE_NORM"},
  {"10.2-2a", "M_h=8 M4/(3 ell)>0 at the card witness",
    "(1)", "LOCALIZED_H_POSITIVE_NORM"},
  {"10.2-2b", "positive reduced scalar stiffness",
    "(1)", "SCALAR_WELL_POSEDNESS"},
  {"10.2-2c", "planar unconstrained wall nonnegative/factorized",
    "(1)", "PLANAR_WALL_FACTORIZATION"},
  {"10.2-2d", "curved held sleeve-slab positivity after constraints and BCs",
    "(2)", "missing lambda_min of assembled constrained wall Hessian L_chi,Sigma^(1)"},
  {"10.2-2e", "homogeneous bulk density Fourier Hessian is positive",
    "(1)", "POSITIVE_BULK_HESSIAN"},
  {"10.2-2f", "bulk-phase constant mode is explicitly quotiented",
    "(1)", "PHASE_CONSTANT_MODE_QUOTIENT"},
  {"10.2-2g", "no negative quantum mode in the inhomogeneous one-body geometry",
    "(2)", "missing lambda_min of assembled one-body Madelung Hessian H_rho,phi^(1) with BC/quotient"},
  {"10.2-3a", "bare projected mouth source is nonzero",
    "(1)", "BARE_MOUTH_MONOPOLE_PREREQUISITE"},
  {"10.2-3b", "dressed response monopole is nonzero",
    "(3)", "missing dressed 1/R coefficient q_h,dressed from the pair/far-field solve"},
  {"10.2-4a", "reduced h has no k^0 h^2 term and Q_s(0,0)=0",
    "(1)", "REDUCED_H_MASSLESSNESS"},
  {"10.2-4b", "parent transverse mode has zero eigenvalue",
    "(1)", "PARENT_TRANSVERSE_ZERO_EIGENVALUE"},
  {"10.2-5a", "Hadamard force cell (1,0) vanishes or is contamination",
    "(3)", "missing four-orientation channel-resolved pair force F_ch^(s1,s2) and H_10 projection"},
  {"10.2-5b", "Hadamard force cell (0,1) vanishes or is contamination",
    "(3)", "missing four-orientation channel-resolved pair force F_ch^(s1,s2) and H_01 projection"},
  {"10.2-6a", "D_i and R_0 independently normalize to one",
    "(1)", "DRAIN_RETURN_KERNEL_NORMALIZATION"},
  {"10.2-6b", "Gamma_0 is positive",
    "(1)", "POSITIVE_DRAIN_THROUGHPUT"},
  {"10.2-6c", "global controller mass identity",
    "(1)", "DRAIN_RETURN_MASS_CLOSURE"},
  {"10.2-6d", "global controller momentum identity",
    "(1)", "DRAIN_RETURN_MOMENTUM_CLOSURE"},
  {"10.2-6e", "global controller energy identity",
    "(1)", "DRAIN_RETURN_ENERGY_CLOSURE"},
  {"10.2-6f", "one-body finite-volume mass residual closes",
    "(2)", "missing R_rho^(1)[C_i] evaluated on assembled one-body solution Psi_i^(1)"},
  {"10.2-6g", "one-body finite-volume momentum residual closes",
    "(2)", "missing R_p^(1)[C_i] evaluated on assembled one-body solution Psi_i^(1)"},
  {"10.2-6h", "one-body finite-volume energy residual closes",
    "(2)", "missing R_epsilon^(1)[C_i] evaluated on assembled one-body solution Psi_i^(1)"},
  {"10.2-6i", "drain ablation nulls its attributed one-body flux piece",
    "(2)", "missing Delta F_flux,drain^(1) from the assembled one-body ablation solve"},
  {"10.2-6j", "pair finite-volume mass residual closes",
    "(3)", "missing R_rho^(2) on the controller-coupled two-center solution Psi^(2)"},
  {"10.2-6k", "pair finite-volume momentum residual closes",
    "(3)", "missing R_p^(2) on the controller-coupled two-center solution Psi^(2)"},
  {"10.2-6l", "pair finite-volume energy residual closes",
    "(3)", "missing R_epsilon^(2) on the controller-coupled two-center solution Psi^(2)"},
  {"10.2-6m", "outer-return ablation nulls its attributed pair flux piece",
    "(3)", "missing Delta F_flux,return^(2) from the shared-return two-center ablation solve"},
  {"10.2-7a", "conservative operator has symmetric Hessian prerequisite",
    "(1)", "CONSERVATIVE_HESSIAN_SYMMETRY_PREREQUISITE"},
  {"10.2-7b", "solved conservative pair force is integrable before U_11",
    "(3)", "missing configuration-space curl/cross-derivatives of solved F_var,11^(2)"},
  {"10.2-7c", "solved stationary total is F_var+F_flux",
    "(3)", "missing channel-resolved pair readouts F_var,F_flux,F_B,F_rad on Psi^(2)"},
  {"10.4-PASS_STABILITY-a", "M and K equal the declared witness matrices",
    "(1)", "SCALAR_WELL_POSEDNESS"},
  {"10.4-PASS_STABILITY-b", "M_h=N0 M4=1>0",
    "(1)", "LOCALIZED_H_POSITIVE_NORM"},
  {"10.4-PASS_STABILITY-c", "scalar stiffness is positive; margin and eig(K) are cross-method evidence",
    "(1)", "SCALAR_WELL_POSEDNESS"},
  {"10.4-PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS", "two positive c_pm^2 roots",
    "(1)", "GENERALIZED_WAVE_SPEEDS"},
  {"10.4-PASS_TRANSVERSE_FACTORIZATION-a", "O_perp f0=0",
    "(1)", "PARENT_TRANSVERSE_ZERO_EIGENVALUE"},
  {"10.4-PASS_TRANSVERSE_FACTORIZATION-b", "O_perp=A^dagger A>=0",
    "(1)", "TRANSVERSE_FACTORIZATION"},
  {"10.4-PASS_WALL_FACTORIZATION", "planar L_chi^(2)/kappa=A_chi^dagger A_chi",
    "(1)", "PLANAR_WALL_FACTORIZATION"},
  {"10.4-PASS_POSITIVE_BULK_HESSIAN", "witness U'' and Fourier coefficients positive",
    "(1)", "POSITIVE_BULK_HESSIAN"},
  {"10.4-PASS_DIMENSIONAL_HOMOGENEITY", "active term/source dimensions agree",
    "(1)", "DIMENSIONAL_HOMOGENEITY"},
  {"10.4-null-a", "bulk-phase constant is the only assembled bulk-phase null",
    "(2)", "missing kernel basis of assembled one-body bulk phase operator L_phi^(1)"},
  {"10.4-null-b", "body translations/wall-shape modes exhaust wall nulls",
    "(2)", "missing kernel basis of assembled constrained wall operator L_chi,Sigma^(1)"},
  {"10.4-null-c", "massless h is the preserved physical far-field Green mode",
    "(3)", "missing zero-frequency far-field Green function G_hh^(2) and pole residue"},
  {"10.4-null-d", "rigid u_L is the only assembled longitudinal null",
    "(2)", "missing kernel basis of assembled one-body longitudinal operator L_u^(1)"},
  {"10.4-fix-a", "phase quotient removes the assembled bulk constant null",
    "(2)", "missing restricted kernel of L_phi^(1) on integral W_IR phi=0"},
  {"10.4-fix-b", "S_hold and translation quotient remove assembled wall nulls",
    "(2)", "missing restricted kernel of L_chi,Sigma^(1) with hold/moduli constraints"},
  {"10.4-fix-c", "IR boundary removes assembled rigid u_L null",
    "(2)", "missing restricted kernel of L_u^(1) with u_L=0 outer BC"},
  {"10.4-fix-d", "null fixing preserves the physical h Green response",
    "(3)", "missing constrained far-field Green function G_hh^(2)"},
  {"10.4-controller", "finite-rank I_ret,P_ret constraints do not spoil positive principal symbol",
    "(2)", "missing spectrum/principal-symbol audit of controller-augmented assembled one-body operator"},
  {"10.4-Fredholm-a", "finite-R_out assembled problem is elliptic/Fredholm",
    "(2)", "missing assembled operator/domain L^(1), Fredholm index, and range conditions"},
  {"10.4-Fredholm-b", "assembled nullspace is fixed and compatibility conditions hold",
    "(2)", "missing kernel/cokernel bases and one-body compatibility residuals"},
  {"10.4-IR", "required IR extrapolation converges",
    "(3)", "missing extrapolated pair observable as R_out,R_ret->infinity and R/R_ret->0"}
};

normalEvaluation = evaluatePredicates[basePrimitives];
normalResults = normalEvaluation["Results"];
normalValues = normalEvaluation["Values"];

scriptArguments = $CommandLine;
mutationPosition = FirstPosition[scriptArguments, "--mutation", Missing["NotFound"]];
mutationMode = ! MissingQ[mutationPosition];

If[mutationMode,
  mutationID = If[
    mutationPosition[[1]] < Length[scriptArguments],
    scriptArguments[[mutationPosition[[1]] + 1]],
    ""
  ];
  If[! MemberQ[checkOrder, mutationID],
    Print["FAIL_UNKNOWN_MUTATION | id=" <> mutationID];
    Quit[2]
  ];
  mutationEvaluation = evaluatePredicates[mutatePrimitives[mutationID]];
  mutationResults = mutationEvaluation["Results"];
  mutationValues = mutationEvaluation["Values"];
  Do[
    Print[
      status[mutationResults[id]] <> "_" <> id <> " | value=" <>
        format[mutationValues[id]]
    ];
    If[! TrueQ[mutationResults[id]],
      Print["FIRST_FAILURE=" <> id];
      Quit[1]
    ],
    {id, checkOrder}
  ];
  Print["FAIL_MUTATION_DID_NOT_FIRE | target=" <> mutationID];
  Quit[2]
];

Print["DISPOSITION_TABLE_BEGIN"];
Do[
  Print[StringRiffle[row, " | "]],
  {row, dispositionRows}
];
Print["DISPOSITION_TABLE_END"];

guard[checkID_String] := Module[{passed = normalResults[checkID]},
  Print[
    status[passed] <> "_" <> checkID <> " | value=" <>
      format[normalValues[checkID]]
  ];
  If[! TrueQ[passed],
    Print["FIRST_FAILURE=" <> checkID];
    Quit[1]
  ]
];

Scan[guard, checkOrder];

(* Each row below is backed by a real child process whose mutation mode     *)
(* exits nonzero at the first failed guard.                                 *)
mathExecutable = First[$CommandLine];
scriptPath = If[StringLength[$InputFileName] > 0, $InputFileName,
  SelectFirst[$ScriptCommandLine, StringMatchQ[#, ___ ~~ ".wl"] &, ""]
];
ablationRows = Table[
  child = RunProcess[
    {mathExecutable, "-script", scriptPath, "--mutation", id},
    All
  ];
  firstFailureText = FirstCase[
    StringSplit[child["StandardOutput"], {"\r\n", "\n"}],
    line_ /; StringStartsQ[line, "FIRST_FAILURE="] :>
      StringDrop[line, StringLength["FIRST_FAILURE="]],
    "NONE"
  ];
  rowPass = TrueQ[child["ExitCode"] != 0] && TrueQ[firstFailureText == id];
  <|
    "Predicate" -> id,
    "Mutated primitive" -> mutationDescriptions[id],
    "First failing check ID" -> firstFailureText,
    "Exit code" -> child["ExitCode"],
    "Pass" -> rowPass
  |>,
  {id, checkOrder}
];
ablationCheck = And @@ Lookup[ablationRows, "Pass"];

Print["ABLATION_MATRIX_BEGIN"];
Do[
  Print[
    row["Predicate"] <> " | " <> row["Mutated primitive"] <>
      " | first=" <> row["First failing check ID"] <>
      " | exit=" <> ToString[row["Exit code"]] <>
      " | " <> status[row["Pass"]]
  ],
  {row, ablationRows}
];
Print["ABLATION_MATRIX_END"];
Print[status[ablationCheck] <> "_MUTATION_SELF_TEST"];
Print[
  "SCOPE_STATIC_ALGEBRAIC_PREREQUISITES_ONLY | " <>
  "class-(2) assembled one-body eigenproblem/BVP and finite-volume claims " <>
  "and class-(3) dressed-monopole, Hadamard, pair-residual, pair-force, " <>
  "far-field/IR claims remain unproved"
];

(* ---------------------------------------------------------------------- *)
(* Clear, machine-readable headline lines and process exit status.         *)
(* ---------------------------------------------------------------------- *)

Print[
  status[scalarWellPosed] <> "_SCALAR_WELL_POSEDNESS | symbolic=" <>
  format[scalarSymbolic] <> " | numeric=" <> format[scalarNumeric]
];
Print[
  status[waveSpeedCheck] <> "_GENERALIZED_WAVE_SPEEDS | symbolic=" <>
  format[speedSymbolic] <> " | numeric=" <> format[speedNumeric]
];
Print[
  status[transverseCheck] <> "_TRANSVERSE_ZERO_MODE_AND_FACTORIZATION | symbolic=" <>
  format[transverseSymbolic] <> " | numeric=" <> format[transverseNumeric]
];
Print[
  status[dimensionCheck] <> "_DIMENSIONAL_HOMOGENEITY | checked=" <>
  ToString[Length[dimensionResiduals]] <> " | symbolic residuals=" <>
  format[dimensionResiduals] <> " | numeric residuals=" <>
  format[dimensionNumeric] <> " | nonzero=" <>
  format[nonzeroDimensionResiduals]
];

allChecksPass = scalarWellPosed && waveSpeedCheck && transverseCheck &&
  dimensionCheck && And @@ Values[normalResults] && ablationCheck;
Print[status[allChecksPass] <> "_G0_CLOSURE_CARD_DUAL_ENGINE"];

Quit[If[TrueQ[allChecksPass], 0, 1]];
