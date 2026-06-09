ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]]
);

fmt[expr_] := ToString[InputForm[expr]];
pass[name_String] := Print["PASS: ", name];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1]
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]]
];

expectVectorZero[name_String, expr_List] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === ConstantArray[0, Length[res]]], pass[name], fail[name, res]]
];

expectMatrixZero[name_String, expr_List] := Module[{res, dims},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  dims = Dimensions[res];
  Print[name, " = ", fmt[res]];
  If[Length[dims] == 2 && TrueQ[res === ConstantArray[0, dims]], pass[name], fail[name, res]]
];

expectNoBranchZero[name_String, expr_] := Module[{zeroSet},
  zeroSet = FullSimplify[
    Reduce[branchAssumptions && expr == 0, branchParameters, Reals],
    Assumptions -> branchAssumptions
  ];
  Print[name, " zero-set = ", fmt[zeroSet]];
  If[TrueQ[zeroSet === False], pass[name], fail[name, zeroSet]]
];

rowOf[expr_] := FullSimplify[
  Coefficient[Expand[expr], #] & /@ rawSlippages,
  Assumptions -> $Assumptions
];

banner["STAGE 183 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT"];

Clear[
  chi0, epsW, epsEta, deltaU,
  sigmaZ, sigmaChi, sigmaEps, sigmaDel, sigmaEta,
  thetaObs, xiObs, rhoObs, coordTr, coordNt, coordEta
];

branchParameters = {chi0, epsW, epsEta, deltaU};
branchAssumptions = (
  Element[branchParameters, Reals] &&
  chi0 > 0 && deltaU > 0 && epsW > 0 && 0 < epsEta < 1
);

rawSlippages = {sigmaZ, sigmaChi, sigmaEps, sigmaDel, sigmaEta};
$Assumptions = (
  branchAssumptions &&
  Element[rawSlippages, Reals] &&
  Element[{thetaObs, xiObs, rhoObs, coordTr, coordNt, coordEta}, Reals]
);

epsBranch = FullSimplify[
  epsW*(1 - (2 deltaU)/(11 (1 + deltaU))),
  Assumptions -> $Assumptions
];

eShapeCoeff = FullSimplify[
  (2 epsW/(1 - epsBranch))*((11 + 9 deltaU)/(11 (1 + deltaU))),
  Assumptions -> $Assumptions
];

deltaOnlyCoeff = FullSimplify[
  (4 epsW deltaU)/(11 (1 - epsBranch) (1 + deltaU)^2),
  Assumptions -> $Assumptions
];

trackingRaw = FullSimplify[
  (1 + chi0) sigmaDel + (1 + deltaU) sigmaChi,
  Assumptions -> $Assumptions
];

thetaRaw = FullSimplify[
  -(chi0 deltaU/((1 + chi0) (1 + chi0 + deltaU))) sigmaChi
  -(chi0 deltaU/((1 + deltaU) (1 + chi0 + deltaU))) sigmaDel,
  Assumptions -> $Assumptions
];

xiRaw = FullSimplify[
  sigmaZ
  + (2 chi0/(1 + chi0)) sigmaChi
  + eShapeCoeff sigmaEps
  - deltaOnlyCoeff sigmaDel,
  Assumptions -> $Assumptions
];

rhoRaw = FullSimplify[
  -(epsEta/(1 - epsEta)) sigmaEta,
  Assumptions -> $Assumptions
];

banner["Raw-slope compiler"];

cFromThetaChi = FullSimplify[
  -Coefficient[thetaRaw, sigmaChi]/(1 + deltaU),
  Assumptions -> $Assumptions
];
cFromThetaDel = FullSimplify[
  -Coefficient[thetaRaw, sigmaDel]/(1 + chi0),
  Assumptions -> $Assumptions
];
aFromXiChi = FullSimplify[
  Coefficient[xiRaw, sigmaChi]/(1 + deltaU),
  Assumptions -> $Assumptions
];
etaDressing = FullSimplify[-Coefficient[rhoRaw, sigmaEta], Assumptions -> $Assumptions];

expectZero["Theta raw row has one tracking prefactor", cFromThetaChi - cFromThetaDel];

cTr = cFromThetaChi;
aTr = aFromXiChi;
etaPref = etaDressing;

sigmaNTDerived = FullSimplify[xiRaw - aTr trackingRaw, Assumptions -> $Assumptions];
sigmaNTCanonical = FullSimplify[
  sigmaZ
  + eShapeCoeff sigmaEps
  - (2 chi0/(1 + deltaU) + deltaOnlyCoeff) sigmaDel,
  Assumptions -> $Assumptions
];

Print["Sigma_tr(raw) = ", fmt[trackingRaw]];
Print["Sigma_nt(raw) = ", fmt[sigmaNTDerived]];
Print["Sigma_eta(raw) = ", fmt[sigmaEta]];
expectZero["Canonical branch-adapted Sigma_nt", sigmaNTDerived - sigmaNTCanonical];

branchCompiler = rowOf /@ {trackingRaw, sigmaNTDerived, sigmaEta};
rawObservableRows = rowOf /@ {thetaRaw, xiRaw, rhoRaw};
triangularMatrix = FullSimplify[
  {
    {-cTr, 0, 0},
    {aTr, 1, 0},
    {0, 0, -etaPref}
  },
  Assumptions -> $Assumptions
];

banner["Triangular ledger by matrix factorization"];
Print["compiler rows = ", fmt[branchCompiler]];
Print["observable rows = ", fmt[rawObservableRows]];
Print["normal-form matrix = ", fmt[triangularMatrix]];

expectMatrixZero[
  "raw observables - triangular matrix.compiler",
  rawObservableRows - triangularMatrix.branchCompiler
];

expectVectorZero[
  "Xi raw - (A_tr Sigma_tr + Sigma_nt) row",
  rowOf[xiRaw - (aTr trackingRaw + sigmaNTDerived)]
];

expectVectorZero[
  "Rho raw + eps_eta/(1-eps_eta) Sigma_eta row",
  rowOf[rhoRaw + (epsEta/(1 - epsEta)) sigmaEta]
];

banner["Symbolic inverse from observable variables"];

inverseRules = First @ Solve[
  Thread[
    {thetaObs, xiObs, rhoObs} ==
    triangularMatrix.{coordTr, coordNt, coordEta}
  ],
  {coordTr, coordNt, coordEta}
];

sigmaTrInverse = FullSimplify[coordTr /. inverseRules, Assumptions -> $Assumptions];
sigmaNtInverse = FullSimplify[coordNt /. inverseRules, Assumptions -> $Assumptions];
sigmaEtaInverse = FullSimplify[coordEta /. inverseRules, Assumptions -> $Assumptions];

ratio = FullSimplify[aTr/cTr, Assumptions -> $Assumptions];
expectedRatio = FullSimplify[2 (1 + chi0 + deltaU)/deltaU, Assumptions -> $Assumptions];

Print["solved Sigma_tr = ", fmt[sigmaTrInverse]];
Print["solved Sigma_nt = ", fmt[sigmaNtInverse]];
Print["solved Sigma_eta = ", fmt[sigmaEtaInverse]];
Print["A_tr/C_tr = ", fmt[ratio]];

expectZero[
  "inverse Sigma_tr formula",
  sigmaTrInverse +
  ((1 + chi0) (1 + deltaU) (1 + chi0 + deltaU)/(chi0 deltaU)) thetaObs
];
expectZero["A_tr/C_tr closed form", ratio - expectedRatio];
expectZero[
  "inverse Sigma_nt formula",
  sigmaNtInverse - (xiObs + expectedRatio thetaObs)
];
expectZero[
  "inverse Sigma_eta formula",
  sigmaEtaInverse + ((1 - epsEta)/epsEta) rhoObs
];

inverseMatrix = FullSimplify[
  Coefficient[#, #2] & @@@ Tuples[{{sigmaTrInverse, sigmaNtInverse, sigmaEtaInverse}, {thetaObs, xiObs, rhoObs}}],
  Assumptions -> $Assumptions
];
inverseMatrix = Partition[inverseMatrix, 3];

expectMatrixZero[
  "inverse matrix.normal-form matrix - identity",
  inverseMatrix.triangularMatrix - IdentityMatrix[3]
];

banner["Triple-rigidity on the branch"];

Print["C_tr = ", fmt[cTr]];
Print["A_tr = ", fmt[aTr]];
Print["eps_eta/(1-eps_eta) = ", fmt[etaPref]];

expectNoBranchZero["C_tr diagonal prefactor", cTr];
expectNoBranchZero["A_tr feed-through prefactor", aTr];
expectNoBranchZero["dressing diagonal prefactor", etaPref];

zeroObservableSet = FullSimplify[
  Reduce[
    branchAssumptions &&
    And @@ Thread[triangularMatrix.{coordTr, coordNt, coordEta} == {0, 0, 0}],
    {coordTr, coordNt, coordEta},
    Reals
  ],
  Assumptions -> branchAssumptions
];
Print["zero-observable branch solution = ", fmt[zeroObservableSet]];
If[
  TrueQ[
    FullSimplify[
      zeroObservableSet == (coordTr == 0 && coordNt == 0 && coordEta == 0),
      Assumptions -> branchAssumptions
    ]
  ],
  pass["zero observables imply zero normal-form defect"],
  fail["zero observables imply zero normal-form defect", zeroObservableSet]
];

Print[""];
Print["Carry-forward formulas:"];
Print["  Sigma_tr = (1+chi_0) Sigma_del + (1+delta_U) Sigma_chi"];
Print["  Sigma_nt = Sigma_Z"];
Print["             + [2 eps_W/(1-eps)] * [(11+9 delta_U)/(11(1+delta_U))] Sigma_eps"];
Print["             - [2 chi_0/(1+delta_U) + 4 eps_W delta_U/(11(1-eps)(1+delta_U)^2)] Sigma_del"];
Print["  Theta_1   = -C_tr Sigma_tr"];
Print["  Xi_1      = A_tr Sigma_tr + Sigma_nt"];
Print["  R_1+Xi_1  = -(eps_eta/(1-eps_eta)) Sigma_eta"];
Print["  Sigma_tr  = -((1+chi_0)(1+delta_U)(1+chi_0+delta_U)/(chi_0 delta_U)) Theta_1"];
Print["  Sigma_nt  = Xi_1 + 2(1+chi_0+delta_U)/delta_U * Theta_1"];
Print["  Sigma_eta = -((1-eps_eta)/eps_eta) (R_1 + Xi_1)"];

Print[""];
Print["Stage 183 Mathematica audit passed."];

Exit[0];
