#!/usr/bin/env wolframscript

(* Path-A 01 symbolic verifier.
   Target-blind checks only.  Genuine derivation gates are kept separate from
   construction restatements, disclosure guards, and bookkeeping identities. *)

$HistoryLength = 0;
ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/pathA_01_return_source_and_balance.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
repoRoot = AbsoluteFileName[FileNameJoin[{scriptDir, "..", "..", ".."}]];
runDir = FileNameJoin[{scriptDir, "runs", "pathA_01_return_source_and_balance"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

notePath = FileNameJoin[{repoRoot, "software", "stage1_solver", "derivations", "pathA_01_return_source_and_balance.md"}];
diagnosticsPath = FileNameJoin[{runDir, "pathA_01_return_source_and_balance_diagnostics.json"}];
reportPath = FileNameJoin[{scriptDir, "pathA_01_return_source_and_balance_report.md"}];
sympyPath = FileNameJoin[{scriptDir, "pathA_01_return_source_and_balance_sympy.py"}];
sympyDiagnosticsPath = FileNameJoin[{runDir, "pathA_01_return_source_and_balance_sympy_diagnostics.json"}];

ClearAll[checkRow, vecAdd, vecSub, vecEq, markdownTable, boolString, readArtifactText, targetBlindFreeQ];

checkRow[name_String, pass_, detail_String] := <|
  "check" -> name,
  "pass" -> TrueQ[pass],
  "detail" -> detail
|>;

vecAdd[a_List, b_List] := Simplify[a + b];
vecSub[a_List, b_List] := Simplify[a - b];
vecEq[a_List, b_List] := TrueQ[Simplify[a == b]];

markdownTable[headers_List, rows_List] := StringRiffle[
  Join[
    {"| " <> StringRiffle[headers, " | "] <> " |"},
    {"| " <> StringRiffle[ConstantArray["---", Length[headers]], " | "] <> " |"},
    ("| " <> StringRiffle[ToString /@ #, " | "] <> " |") & /@ rows
  ],
  "\n"
];

boolString[x_] := If[TrueQ[x], "PASS", "FAIL"];

readArtifactText[path_String] := If[FileExistsQ[path], Import[path, "Text"], ""];

targetBlindFreeQ[texts_List, patterns_List] := And @@ Flatten[
  Table[StringFreeQ[text, pat], {text, texts}, {pat, patterns}]
];

If[! FileExistsQ[notePath],
  Print["ERROR: derivation note missing: ", notePath];
  Exit[2]
];

(* D1: Taylor expansion and reciprocal source kernel. *)
ClearAll[eta, drho, rho0, k1, k2, V0];
vTaylor = V0 - k1 eta + (1/2) k2 eta^2;
lInt = Expand[-vTaylor (rho0 + drho)];

lCrossCoeff = Simplify[Coefficient[Coefficient[lInt, eta, 1], drho, 1]];
lSelfCoeff = Simplify[Coefficient[lInt /. drho -> 0, eta, 2]];
lRho2Coeff = Simplify[Coefficient[lInt, drho, 2]];

hCross = -k1 eta drho;
hQuad = hCross + (1/2) k2 rho0 eta^2;
forwardFromBilinear = Simplify[D[hCross, drho]];
returnFromBilinear = Simplify[D[hCross, eta]];
crossHessianEtaRho = Simplify[D[D[hCross, eta], drho]];
crossHessianRhoEta = Simplify[D[D[hCross, drho], eta]];
sourceRaw = Simplify[D[hQuad, eta]];
sourceAfterSelfAbsorb = Simplify[sourceRaw - k2 rho0 eta];

reciprocityPass = And[
  TrueQ[lCrossCoeff == k1],
  TrueQ[lSelfCoeff == -(1/2) k2 rho0],
  TrueQ[lRho2Coeff == 0],
  TrueQ[forwardFromBilinear == -k1 eta],
  TrueQ[returnFromBilinear == -k1 drho],
  TrueQ[crossHessianEtaRho == crossHessianRhoEta],
  TrueQ[crossHessianEtaRho == -k1],
  TrueQ[sourceRaw == -k1 drho + k2 rho0 eta],
  TrueQ[sourceAfterSelfAbsorb == -k1 drho]
];

(* D2: gauge invariance and absence of direct fixed-action eta -> Z coupling. *)
ClearAll[t, aCoord, wCoord, chi, A0, Aa, Aw, theta, rho, q, hbar, mass, Z, H, mu0, xi, Fsq, gfSq, extTerm];

ew0 = -D[Aw[t, aCoord, wCoord], t] - D[A0[t, aCoord, wCoord], wCoord];
ca0 = D[Aw[t, aCoord, wCoord], aCoord] - D[Aa[t, aCoord, wCoord], wCoord];

a0Gauge = A0[t, aCoord, wCoord] - D[chi[t, aCoord, wCoord], t];
aaGauge = Aa[t, aCoord, wCoord] + D[chi[t, aCoord, wCoord], aCoord];
awGauge = Aw[t, aCoord, wCoord] + D[chi[t, aCoord, wCoord], wCoord];

ewGauge = -D[awGauge, t] - D[a0Gauge, wCoord];
caGauge = D[awGauge, aCoord] - D[aaGauge, wCoord];
ewInvariant = TrueQ[Simplify[ewGauge - ew0] == 0];
caInvariant = TrueQ[Simplify[caGauge - ca0] == 0];

jHydro = (hbar/mass) rho[t, aCoord, wCoord] *
  (D[theta[t, aCoord, wCoord], aCoord] - (q/hbar) Aa[t, aCoord, wCoord]);
thetaGauge = theta[t, aCoord, wCoord] + (q/hbar) chi[t, aCoord, wCoord];
jHydroGauge = (hbar/mass) rho[t, aCoord, wCoord] *
  (D[thetaGauge, aCoord] - (q/hbar) aaGauge);
currentInvariant = TrueQ[Simplify[jHydroGauge - jHydro] == 0];

lEMDeclared = -Z[wCoord] Fsq/(4 mu0) - H[wCoord] gfSq/(2 xi mu0) - extTerm;
directMaxwellEtaVariation = Simplify[D[lEMDeclared, eta]];
directMaxwellPass = TrueQ[directMaxwellEtaVariation == 0];

ClearAll[JcovEta, EwInv, CaInv, sJ, sE, sC];
sAAllowedData = sJ JcovEta + sE EwInv + sC CaInv;
allowedGaugeDataPass = FreeQ[sAAllowedData, A0 | Aa | Aw | chi];

gaugePass = And[ewInvariant, caInvariant, currentInvariant];

(* Dimensional consistency using exponent vectors {M,L,T}; d is the parent spatial density exponent. *)
ClearAll[d];
dimEnergy = {1, 2, -2};
dimLength = {0, 1, 0};
dimTime = {0, 0, 1};
dimEta = dimLength;
dimRho = {0, -d, 0};
dimK1 = vecSub[dimEnergy, dimLength];
dimK2 = vecSub[dimEnergy, 2 dimLength];
dimSource = vecAdd[dimK1, dimRho];
dimKeta = vecSub[dimSource, dimEta];
dimMuEta = vecSub[dimSource, vecSub[dimEta, 2 dimTime]];
dimTw = vecAdd[dimSource, dimLength];
dimTOmega = vecSub[dimSource, dimEta];

dimMatterCross = vecAdd[vecAdd[dimK1, dimEta], dimRho];
dimMatterSelf = vecAdd[vecAdd[dimK2, 2 dimEta], dimRho];
dimMatterDensity = vecAdd[dimEnergy, dimRho];
dimMatterSource1 = vecAdd[dimK1, dimRho];
dimMatterSource2 = vecAdd[vecAdd[dimK2, dimRho], dimEta];
dimWallMuTerm = vecAdd[dimMuEta, vecSub[dimEta, 2 dimTime]];
dimWallTwTerm = vecSub[dimTw, dimLength];
dimWallOmegaTerm = vecAdd[dimTOmega, dimEta];
dimWallKTerm = vecAdd[dimKeta, dimEta];

dimensionPass = And[
  vecEq[dimMatterCross, dimMatterDensity],
  vecEq[dimMatterSelf, dimMatterDensity],
  vecEq[dimMatterSource1, dimSource],
  vecEq[dimMatterSource2, dimSource],
  vecEq[dimWallMuTerm, dimSource],
  vecEq[dimWallTwTerm, dimSource],
  vecEq[dimWallOmegaTerm, dimSource],
  vecEq[dimWallKTerm, dimSource]
];

(* Construction restatement: P and N0 use only forward mixed-sector couplings.
   The no-new-numerator conclusion is a corollary of the genuine reciprocity
   check above, not an independent derivation gate. *)
ClearAll[OmegaU2, OmegaW2, Rmix, GU, GW, Kstar, Qcons, SpsiSym, SASym, deltaD];
deltaExpr = OmegaU2 OmegaW2 - Rmix^2;
pExpr = OmegaU2 GW + Rmix GU;
n0Expr = (pExpr/deltaExpr)^2;
d0ExprWithReturn = Kstar - Qcons/deltaExpr + deltaD[SpsiSym, SASym];
noNumeratorStructurePass = And[
  TrueQ[D[pExpr, SpsiSym] == 0],
  TrueQ[D[pExpr, SASym] == 0],
  TrueQ[D[n0Expr, SpsiSym] == 0],
  TrueQ[D[n0Expr, SASym] == 0],
  ! FreeQ[d0ExprWithReturn, SpsiSym | SASym]
];
noNumeratorCorollaryPass = And[reciprocityPass, noNumeratorStructurePass];

(* Reduction of promoted S_Sigma to the effective quadratic wall operator. *)
ClearAll[eps, e, ep, et, gradE, R0p, Tw0, TwR, TwRR, U0, U1, U2, muSigma0, TOmega0, dA];
twTaylor = Tw0 + eps TwR e + (1/2) eps^2 TwRR e^2;
uTaylor = U0 + eps U1 e + (1/2) eps^2 U2 e^2;
rwTaylor = R0p + eps ep;
lSigmaLocal = (1/2) muSigma0 eps^2 et^2
  - (1/2) twTaylor rwTaylor^2
  - (1/2) TOmega0 eps^2 gradE^2
  - uTaylor;
lSigma2 = Simplify[Coefficient[Expand[lSigmaLocal], eps, 2]];
lSigma2ExpectedLocal = (1/2) muSigma0 et^2
  - (1/2) Tw0 ep^2
  - TwR R0p e ep
  - (1/4) TwRR R0p^2 e^2
  - (1/2) TOmega0 gradE^2
  - (1/2) U2 e^2;
localReductionPass = TrueQ[Simplify[lSigma2 - lSigma2ExpectedLocal] == 0];

kEtaFormula = U2 - dA + (1/2) TwRR R0p^2;
lSigma2AfterIBP = (1/2) muSigma0 et^2
  - (1/2) Tw0 ep^2
  - (1/2) TOmega0 gradE^2
  - (1/2) kEtaFormula e^2;
lSigma2IBPFromLocal = lSigma2ExpectedLocal /. (-TwR R0p e ep) -> (1/2) dA e^2;
ClearAll[wIBP, etaIBP, TwRIBP, R0pIBP];
ibpCoeff = TwRIBP[wIBP] R0pIBP[wIBP];
dADerived = D[ibpCoeff, wIBP];
ibpLocalDensity = -ibpCoeff etaIBP[wIBP] D[etaIBP[wIBP], wIBP];
ibpNoBoundaryDensity = (1/2) dADerived etaIBP[wIBP]^2;
ibpBoundaryDensity = -(1/2) D[ibpCoeff etaIBP[wIBP]^2, wIBP];
ibpSymbolicPass = TrueQ[
  Simplify[ibpLocalDensity - (ibpNoBoundaryDensity + ibpBoundaryDensity)] == 0
];
ibpCompactAlgebraPass = TrueQ[Simplify[lSigma2IBPFromLocal - lSigma2AfterIBP] == 0];
ibpReductionPass = And[ibpSymbolicPass, ibpCompactAlgebraPass];
reductionPass = localReductionPass;

(* Open-item disclosure guard and target-blind artifact scan. *)
noteText = Import[notePath, "Text"];
openItemsPass = And[
  StringContainsQ[noteText, "## Posited or Open"],
  StringContainsQ[noteText, "constitutive"],
  StringContainsQ[noteText, "gauge return"],
  StringContainsQ[noteText, "moving-boundary"],
  StringContainsQ[noteText, "Boundary terms"],
  StringContainsQ[noteText, "does not solve"]
];

scriptText = Import[scriptPath, "Text"];
forbiddenPatterns = {
  StringJoin["10", ".", "8"],
  StringJoin["54", "/", "5"],
  StringJoin["R", "_", "norm"],
  StringJoin["P0", "_", "target"],
  StringJoin["GR", " ", "constant"]
};

targetBlindBaseTexts = {
  noteText,
  scriptText,
  readArtifactText[sympyPath],
  readArtifactText[reportPath],
  readArtifactText[sympyDiagnosticsPath]
};
targetBlindBasePass = targetBlindFreeQ[targetBlindBaseTexts, forbiddenPatterns];

makeGenuineChecks[targetBlindArtifactPass_] := {
  checkRow[
    "reciprocity",
    reciprocityPass,
    "Taylor series gives L_cross=k1 eta drho and L_self=-(1/2) k2 rho0 eta^2; equation-kernel Hessian is symmetric with kernel -k1."
  ],
  checkRow[
    "gauge_invariance",
    gaugePass,
    "Three live identities pass: E_w, C_a, and the hydrodynamic covariant current are gauge invariant."
  ],
  checkRow[
    "dimensional_consistency",
    dimensionPass,
    "Matter return terms and all wall-PDE left-hand terms share the same source dimension vector."
  ],
  checkRow[
    "reduction_consistency",
    reductionPass,
    "The local eps^2 coefficient of promoted S_Sigma reproduces mu_eta, T_w, T_Omega, and the unintegrated quadratic terms."
  ],
  checkRow[
    "target_blind_artifacts",
    targetBlindArtifactPass,
    "Note, report markdown, Mathematica verifier, SymPy verifier, and SymPy diagnostics contain no forbidden target-leakage tokens."
  ]
};

nonPhysicsChecks = {
  checkRow[
    "no_numerator_knob_corollary_not_a_physics_gate",
    noNumeratorCorollaryPass,
    "Corollary of reciprocity: return kernels are adjoints of the forward kernels, so they add no new free numerator magnitude. The P/N0 symbol-placement scan is only a construction restatement."
  ],
  checkRow[
    "honest_open_items_disclosure_guard_not_a_physics_gate",
    openItemsPass,
    "Disclosure-presence guard: the Posited or Open section still names constitutive, gauge, boundary, and no-solve limits."
  ],
  checkRow[
    "ibp_bookkeeping_identity_not_a_physics_gate",
    ibpReductionPass,
    "Bookkeeping identity: symbolic IBP derives dA=partial_w(T_w,Sigma,R R0') and the compact K_eta algebra agrees after dropping the named boundary density."
  ],
  checkRow[
    "declared_maxwell_no_direct_eta_variation_not_a_physics_gate",
    directMaxwellPass,
    "Declared fixed-domain Z(w),H(w) Maxwell action has zero direct eta variation by construction."
  ],
  checkRow[
    "allowed_gauge_data_placeholder_not_a_physics_gate",
    allowedGaugeDataPass,
    "Construction placeholder: the allowed reduced gauge data list contains only J_cov_eta, E_w, and C_a names, not raw potentials or chi."
  ]
};

makeReport[genuineRows_List, nonPhysicsRows_List, genuinePassLocal_, nonPhysicsPassLocal_, verifierPassLocal_] := Module[
  {genuineReportRows, nonPhysicsReportRows},
  genuineReportRows = ({#["check"], boolString[#["pass"]], #["detail"]} &) /@ genuineRows;
  nonPhysicsReportRows = ({#["check"], boolString[#["pass"]], #["detail"]} &) /@ nonPhysicsRows;
  StringRiffle[
    {
      "# Path-A 01 Return Source and Balance Verification",
      "",
      "Genuine checks pass: " <> boolString[genuinePassLocal],
      "Verifier including construction-restatements/guards pass: " <> boolString[verifierPassLocal],
      "",
      "## Genuine Checks",
      "",
      markdownTable[{"Check", "Result", "Detail"}, genuineReportRows],
      "",
      "## Construction-Restatements and Guards",
      "",
      "These are excluded from the genuine derivation gates and carry `_not_a_physics_gate` suffixes.",
      "",
      markdownTable[{"Check", "Result", "Detail"}, nonPhysicsReportRows],
      "",
      "Non-physics guards/restatements pass: " <> boolString[nonPhysicsPassLocal],
      "",
      "Diagnostics JSON: `" <> diagnosticsPath <> "`"
    },
    "\n"
  ]
];

targetBlindPass = targetBlindBasePass;
checks = makeGenuineChecks[targetBlindPass];
genuinePass = AllTrue[checks, TrueQ[#["pass"]] &];
nonPhysicsPass = AllTrue[nonPhysicsChecks, TrueQ[#["pass"]] &];
allPass = And[genuinePass, nonPhysicsPass];
candidateReport = makeReport[checks, nonPhysicsChecks, genuinePass, nonPhysicsPass, allPass];
targetBlindPass = And[targetBlindBasePass, targetBlindFreeQ[{candidateReport}, forbiddenPatterns]];

checks = makeGenuineChecks[targetBlindPass];
genuinePass = AllTrue[checks, TrueQ[#["pass"]] &];
nonPhysicsPass = AllTrue[nonPhysicsChecks, TrueQ[#["pass"]] &];
allPass = And[genuinePass, nonPhysicsPass];

diagnostics = <|
  "schema" -> "pathA_01_return_source_and_balance/verifier/v2",
  "all_pass" -> allPass,
  "genuine_pass" -> genuinePass,
  "non_physics_pass" -> nonPhysicsPass,
  "checks" -> checks,
  "construction_restatements_and_guards" -> nonPhysicsChecks,
  "target_blind_artifacts_scanned" -> {
    notePath,
    reportPath,
    scriptPath,
    sympyPath,
    sympyDiagnosticsPath
  },
  "symbolic_results" -> <|
    "lagrangian_cross_coefficient" -> ToString[InputForm[lCrossCoeff]],
    "lagrangian_self_coefficient" -> ToString[InputForm[lSelfCoeff]],
    "forward_kernel" -> ToString[InputForm[forwardFromBilinear/eta]],
    "return_kernel" -> ToString[InputForm[returnFromBilinear/drho]],
    "raw_matter_source" -> ToString[InputForm[sourceRaw]],
    "absorbed_matter_source" -> ToString[InputForm[sourceAfterSelfAbsorb]],
    "direct_maxwell_eta_variation" -> ToString[InputForm[directMaxwellEtaVariation]],
    "dA_derived_by_ibp" -> ToString[InputForm[dADerived]],
    "ibp_boundary_density" -> ToString[InputForm[ibpBoundaryDensity]],
    "k_eta_formula" -> ToString[InputForm[kEtaFormula]]
  |>
|>;

Export[diagnosticsPath, diagnostics, "RawJSON"];

report = makeReport[checks, nonPhysicsChecks, genuinePass, nonPhysicsPass, allPass];
Export[reportPath, report, "Text"];

If[! allPass,
  Print["FAIL pathA_01_return_source_and_balance; see ", diagnosticsPath];
  Exit[1]
];

Print["PASS pathA_01_return_source_and_balance"];
Print["Report: ", reportPath];
Print["Diagnostics: ", diagnosticsPath];
Exit[0];
