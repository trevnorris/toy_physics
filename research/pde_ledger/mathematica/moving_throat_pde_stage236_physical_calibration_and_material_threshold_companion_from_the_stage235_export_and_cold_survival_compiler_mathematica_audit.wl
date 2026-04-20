ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

normalizeExpr[expr_] := Module[{res},
  res = If[
    MatrixQ[expr],
    Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr, {2}],
    If[
      VectorQ[expr],
      Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr],
      FullSimplify[expr, Assumptions -> $Assumptions]
    ]
  ];
  res
];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# == 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr == 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

expectNear[name_String, value_, target_, tol_] := Module[{res},
  res = N[Abs[value - target], 20];
  Print[name, " diff = ", fmt[res], " (tol = ", fmt[tol], ")"];
  If[res <= tol, pass[name], fail[name, res]]
];

banner["STAGE 236 — PHYSICAL CALIBRATION AND MATERIAL-THRESHOLD COMPANION"];

Clear[
  sc, s0, muEta, fLat, zetaEp, tStar, UpsilonLat, lambdaEpOmegaD,
  gammaLatticeLegacy, lambdaPhys, lambdaRef, rTurn, EStar, VprimeTurnAbs,
  aInt, KTurnSym, KCorr, T, kEff
];

$Assumptions =
  Element[
    {
      sc, s0, muEta, fLat, zetaEp, tStar, UpsilonLat, lambdaEpOmegaD,
      gammaLatticeLegacy, lambdaPhys, lambdaRef, rTurn, EStar,
      VprimeTurnAbs, aInt, KTurnSym, KCorr, T, kEff
    },
    Reals
  ] &&
  sc > 0 && s0 > 0 && muEta > 0 && fLat > 0 && zetaEp > 0 && tStar > 0 &&
  UpsilonLat > 0 && lambdaEpOmegaD > 0 && gammaLatticeLegacy > 0 &&
  lambdaPhys > 0 && lambdaRef > 0 && rTurn > 0 && EStar > 0 &&
  VprimeTurnAbs > 0 && aInt > 0 && KTurnSym > 0 && KCorr > 0 && T > 0 &&
  kEff > 0;

subbanner["I. Exact lattice-turnover compiler"];

gammaLatSafeEq = FullSimplify[fLat muEta (s0^2 - sc^2)/sc, Assumptions -> $Assumptions];
tCrossPhys = FullSimplify[tStar/sc, Assumptions -> $Assumptions];
gammaLatTurnPhys = FullSimplify[
  gammaLatSafeEq/(UpsilonLat tStar),
  Assumptions -> $Assumptions
];
thresholdLambda = FullSimplify[gammaLatTurnPhys/zetaEp, Assumptions -> $Assumptions];
thresholdProduct = FullSimplify[thresholdLambda tCrossPhys, Assumptions -> $Assumptions];
thresholdLambdaExpected = FullSimplify[
  fLat muEta (s0^2 - sc^2)/(UpsilonLat zetaEp sc tStar),
  Assumptions -> $Assumptions
];
thresholdProductExpected = FullSimplify[
  fLat muEta (s0^2 - sc^2)/(UpsilonLat zetaEp sc^2),
  Assumptions -> $Assumptions
];

Print["gamma_lat,safe^eq = ", fmt[gammaLatSafeEq]];
Print["t_cross^phys = ", fmt[tCrossPhys]];
Print["gamma_lat,turn^phys = ", fmt[gammaLatTurnPhys]];
Print["(lambda_ep omega_D)_min = ", fmt[thresholdLambda]];
Print["product threshold = ", fmt[thresholdProduct]];

expectZero["lambda-threshold compiler", thresholdLambda - thresholdLambdaExpected];
expectZero["product-threshold compiler", thresholdProduct - thresholdProductExpected];
expectZero[
  "product-threshold reduced form",
  thresholdProduct - gammaLatSafeEq/(UpsilonLat zetaEp sc)
];

subbanner["II. Legacy Session-V slice recovery"];

UpsilonLegacy = FullSimplify[
  gammaLatSafeEq/gammaLatticeLegacy,
  Assumptions -> $Assumptions
];
thresholdLambdaLegacy = FullSimplify[
  thresholdLambda /. UpsilonLat -> UpsilonLegacy,
  Assumptions -> $Assumptions
];
thresholdProductLegacy = FullSimplify[
  thresholdProduct /. UpsilonLat -> UpsilonLegacy,
  Assumptions -> $Assumptions
];

Print["Upsilon_lat^(sess) = ", fmt[UpsilonLegacy]];
Print["legacy lambda threshold = ", fmt[thresholdLambdaLegacy]];
Print["legacy product threshold = ", fmt[thresholdProductLegacy]];

expectZero[
  "legacy lambda-threshold recovery",
  thresholdLambdaLegacy - gammaLatticeLegacy/(zetaEp tStar)
];
expectZero[
  "legacy product-threshold recovery",
  thresholdProductLegacy - gammaLatticeLegacy/(zetaEp sc)
];

subbanner["III. Harmonic trigger and force-matched stiffness compiler"];

rTurnPhys = FullSimplify[lambdaPhys rTurn/lambdaRef, Assumptions -> $Assumptions];
chiLambdaLattice = FullSimplify[2 lambdaPhys/rTurnPhys, Assumptions -> $Assumptions];
kEffReq = FullSimplify[
  EStar lambdaRef VprimeTurnAbs/(lambdaPhys rTurnPhys),
  Assumptions -> $Assumptions
];
KTurn = FullSimplify[
  lambdaRef^2 VprimeTurnAbs/rTurn,
  Assumptions -> $Assumptions
];
kEffReqKTurn = FullSimplify[
  kEffReq /. VprimeTurnAbs -> KTurnSym rTurn/lambdaRef^2,
  Assumptions -> $Assumptions
];
kEffReqAInt = FullSimplify[
  kEffReqKTurn /. lambdaPhys -> aInt/2,
  Assumptions -> $Assumptions
];

Print["r_turn^phys = ", fmt[rTurnPhys]];
Print["chi_lambda,lattice = ", fmt[chiLambdaLattice]];
Print["k_eff,req = ", fmt[kEffReq]];
Print["K_turn = ", fmt[KTurn]];
Print["k_eff,req(K_turn) = ", fmt[kEffReqKTurn]];
Print["k_eff,req(a_int/2) = ", fmt[kEffReqAInt]];

expectZero["turning-point radius map", rTurnPhys - lambdaPhys rTurn/lambdaRef];
expectZero["chi_lambda geometry ratio", chiLambdaLattice - 2 lambdaRef/rTurn];
expectZero[
  "force-matched stiffness compiler",
  kEffReq - EStar lambdaRef^2 VprimeTurnAbs/(lambdaPhys^2 rTurn)
];
expectZero["K_turn definition", KTurn - lambdaRef^2 VprimeTurnAbs/rTurn];
expectZero["K_turn stiffness rewrite", kEffReqKTurn - KTurnSym EStar/lambdaPhys^2];
expectZero["a_int stiffness rewrite", kEffReqAInt - 4 KTurnSym EStar/aInt^2];

subbanner["IV. Korringa-limited temperature ceiling"];

TMax = FullSimplify[KCorr/tCrossPhys, Assumptions -> $Assumptions];
TMaxExpected = FullSimplify[sc KCorr/tStar, Assumptions -> $Assumptions];
PiT = FullSimplify[TMax/T, Assumptions -> $Assumptions];

Print["T_max = ", fmt[TMax]];
Print["Pi_T = ", fmt[PiT]];

expectZero["Korringa ceiling", TMax - TMaxExpected];
expectZero["thermal screening ratio", PiT - KCorr/(T tCrossPhys)];

subbanner["V. Exact screening ratios"];

PiEp = FullSimplify[
  UpsilonLat zetaEp lambdaEpOmegaD tStar/gammaLatSafeEq,
  Assumptions -> $Assumptions
];
PiChi = FullSimplify[chiLambdaLattice, Assumptions -> $Assumptions];
PiK = FullSimplify[
  kEff lambdaPhys^2/(KTurnSym EStar),
  Assumptions -> $Assumptions
];

Print["Pi_ep = ", fmt[PiEp]];
Print["Pi_chi = ", fmt[PiChi]];
Print["Pi_k = ", fmt[PiK]];
Print["Pi_T = ", fmt[PiT]];

expectZero["Pi_ep threshold ratio", PiEp - lambdaEpOmegaD/thresholdLambda];
expectZero["Pi_chi reduced ratio", PiChi - 2 lambdaRef/rTurn];
expectZero["Pi_k stiffness ratio", PiK - kEff/kEffReqKTurn];

subbanner["VI. Stage-235 / Session-V benchmark-only specialization"];

scNum = 0.5489386551062235;
s0Num = 6.94311167;
fLatNum = 0.75;
muEtaNum = 1.0;
gammaLatticeLegacyNum = 4.79562976;
lambdaRefNum = 0.42826825;
rTurnNum = 0.39096144;
KTurnNum = 2.73855812;

gammaLatSafeEqNum = N[fLatNum muEtaNum (s0Num^2 - scNum^2)/scNum, 16];
UpsilonLegacyNum = N[gammaLatSafeEqNum/gammaLatticeLegacyNum, 16];
productMicroNum = N[gammaLatSafeEqNum/scNum, 16];
productLegacyNum = N[gammaLatticeLegacyNum/scNum, 16];
rTurnPhysCoeffNum = N[rTurnNum/lambdaRefNum, 16];
chiLambdaNum = N[2.0 lambdaRefNum/rTurnNum, 16];
TMaxCoeffNum = N[scNum, 16];
aIntCoeffNum = N[4.0 KTurnNum, 16];
VprimeTurnAbsNum = N[KTurnNum rTurnNum/lambdaRefNum^2, 16];

Print["s_c benchmark = ", fmt[scNum]];
Print["gamma_lat,safe^eq = ", fmt[gammaLatSafeEqNum]];
Print["Upsilon_lat^(sess) = ", fmt[UpsilonLegacyNum]];
Print["micro product threshold = ", fmt[productMicroNum], " / zeta_ep"];
Print["legacy product threshold = ", fmt[productLegacyNum], " / zeta_ep"];
Print["r_turn^phys coefficient = ", fmt[rTurnPhysCoeffNum], " * lambda_phys"];
Print["chi_lambda benchmark = ", fmt[chiLambdaNum]];
Print["K_turn benchmark = ", fmt[KTurnNum]];
Print["a_int stiffness coefficient = ", fmt[aIntCoeffNum]];
Print["T_max coefficient = ", fmt[TMaxCoeffNum], " * K_corr / t_*"];
Print["|V'_red|(r_turn) = ", fmt[VprimeTurnAbsNum]];

expectNear["gamma_lat,safe^eq(session)", gammaLatSafeEqNum, 65.45193925961132, 10^-9];
expectNear["Upsilon_lat^(sess)", UpsilonLegacyNum, 13.64824695299483, 10^-9];
expectNear["micro product threshold", productMicroNum, 119.23361317476524, 10^-9];
expectNear["legacy product threshold", productLegacyNum, 8.736185210116078, 10^-9];
expectNear["r_turn^phys coefficient", rTurnPhysCoeffNum, 0.9128891530016525, 10^-12];
expectNear["chi_lambda benchmark", chiLambdaNum, 2.1908464937104797, 10^-12];
expectNear["a_int stiffness coefficient", aIntCoeffNum, 10.95423248, 10^-12];
expectNear["T_max coefficient", TMaxCoeffNum, 0.5489386551062235, 10^-15];
expectNear["|V'_red|(r_turn)", VprimeTurnAbsNum, 5.837462857946154, 10^-12];

Print[""];
Print["All Stage-236 symbolic and benchmark-specialization checks passed."];
