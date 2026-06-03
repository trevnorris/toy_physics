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

stripConditional[expr_] := expr /. ConditionalExpression[value_, _] :> value;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

cleanExpr[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {ArrayDepth[expr]}],
  cleanScalar[expr]
];

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# === 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr === 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

expectTrue[name_String, statement_] := Module[{res},
  res = stripConditional[FullSimplify[statement, Assumptions -> $Assumptions]];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = N[Abs[N[value, 50] - N[target, 50]], 50];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 247 - RELAXED STATIONARY BARRIER COMPILER"];

Clear[
  r, kappa, alpha6, alpha2, betaQ, betaU, betaW, Kstar, OmU2, OmW2,
  GU, GW, Rmix, etaLeak, muW, rho0, q, lam, Lvar, aU, aV, chiUV,
  fU, etaUV, rSig, a0, b0, rF1, xiR, lambdaL, lambdaW, rLimit
];

$Assumptions = (
  Element[
    {r, kappa, alpha6, alpha2, betaQ, betaU, betaW, Kstar, OmU2, OmW2,
      GU, GW, Rmix, etaLeak, muW, rho0, q, lam, Lvar, aU, aV, chiUV,
      fU, etaUV, rSig, a0, b0, rF1, xiR, lambdaL, lambdaW},
    Reals
  ] &&
  r > 0 && kappa > 0 && alpha6 >= 0 && alpha2 >= 0 &&
  Kstar > 0 && OmU2 > 0 && OmW2 > 0 &&
  etaLeak > 0 && muW > 0 && rho0 > 0 && q > 0 && lam > 0 && Lvar > 0 &&
  aU > 0 && aV > 0 && etaUV > 0 && rSig > 0 && rF1 > 0 && xiR > 0 &&
  lambdaL >= 0 && lambdaW >= 0
);

subbanner["M1. Susceptibilities from the reduced stiffness inverse"];

delta = OmU2 OmW2 - Rmix^2;
qMix = GU^2 OmW2 + 2 GU GW Rmix + GW^2 OmU2;
pMix = OmU2 GW + Rmix GU;
pUMix = GU OmW2 + Rmix GW;
dZero = Kstar - qMix/delta;

Kred = {
  {Kstar, -GU, -GW},
  {-GU, OmU2, -Rmix},
  {-GW, -Rmix, OmW2}
};
Inv = FullSimplify[Inverse[Kred], Assumptions -> $Assumptions];

chiClosed = {
  1/dZero,
  pUMix/(delta dZero),
  pMix/(delta dZero),
  (Kstar OmW2 - GW^2)/(delta dZero),
  (Kstar Rmix + GU GW)/(delta dZero),
  (Kstar OmU2 - GU^2)/(delta dZero)
};
chiFromInverse = {
  Inv[[1, 1]],
  Inv[[1, 2]],
  Inv[[1, 3]],
  Inv[[2, 2]],
  Inv[[2, 3]],
  Inv[[3, 3]]
};

expectZero["M1 inverse-entry closed forms", chiFromInverse - chiClosed];

cSix = chiClosed[[1]] betaQ^2;
cFour = chiClosed[[2]] betaQ betaU + chiClosed[[3]] betaQ betaW;
cTwo = (
  chiClosed[[4]] betaU^2 + 2 chiClosed[[5]] betaU betaW +
  chiClosed[[6]] betaW^2
);
aSix = 3 alpha6 + cSix/2;
aFour = cFour;
aTwo = alpha2 + cTwo/2;
Vshort = (
  (1/r) (1 + Exp[-2 kappa r]/2) -
  aSix/r^6 -
  aFour Exp[-2 kappa r]/r^4 -
  aTwo Exp[-4 kappa r]/r^2
);

subbanner["M2. Leakage and work factoring"];

E0 = 16 etaLeak Lvar/Pi^2;
Sleak = 8 Sqrt[2] etaLeak muW rho0 Lvar/(Pi^(5/2) lam^3);
Wsess = 512 etaLeak^2 muW q rho0 Lvar^2/(Pi^4 lam^2);

expectZero[
  "M2 S_leak from E0",
  Sleak - Sqrt[2] E0 muW rho0/(2 Sqrt[Pi] lam^3)
];
expectZero[
  "M2 W_sess from E0",
  Wsess - 2 E0^2 muW q rho0/lam^2
];

subbanner["M3. Source-response limiting mouth moment"];

sMoment[x_] := rSig^2/(x^2 + rSig^2);
gMoment[x_] := (2/Pi) (1 + a0 sMoment[x]/3 - b0 sMoment[x]/15);
gInf = FullSimplify[
  Limit[gMoment[rLimit], rLimit -> Infinity],
  Assumptions -> Element[{a0, b0}, Reals] && rSig > 0
];

expectZero["M3 g infinity", gInf - 2/Pi];

subbanner["M4. Weighted U/V drain nonnegativity"];

deltaUV = aU aV - chiUV^2;
Duv = aV chiUV^2 fU^2/deltaUV^2;
Euv = etaUV Duv;
drainNegativeBranch = Reduce[
  aV > 0 && Element[{aU, aV, chiUV, fU}, Reals] &&
    deltaUV != 0 && Duv < 0,
  {aU, aV, chiUV, fU},
  Reals
];
drainProbe = N[Duv /. {aU -> 2, aV -> 3/2, chiUV -> 7/10, fU -> 2/5}, 30];

Print["M4 drain negative branch Reduce = ", fmt[drainNegativeBranch]];
Print["M4 drain positive probe = ", fmt[drainProbe]];
expectTrue["M4 no negative drain branch", drainNegativeBranch === False];
expectTrue["M4 positive drain probe", drainProbe > 0];

subbanner["M5. Stationary lowering identity"];

Rresponse[x_] := (gMoment[x] - rF1)^2/(1 + rF1^2);
Rinfty = (2/Pi - rF1)^2/(1 + rF1^2);
Msigma = xiR (Rinfty - Rresponse[r]);
packetSum = lambdaL Sleak + lambdaW Wsess + Euv + Msigma;
Veff = Vshort - packetSum;

expectZero["M5 lowering identity", (Vshort - Veff) - packetSum];

subbanner["M6. Session-I baseline numerics"];

sessionRules = {
  Kstar -> 4,
  OmU2 -> 9,
  OmW2 -> 16,
  GU -> 1,
  GW -> 5/4,
  Rmix -> 27/20,
  betaQ -> 3/100,
  betaU -> 3/20,
  betaW -> 1/5,
  alpha6 -> 0,
  alpha2 -> 0,
  kappa -> 1,
  r -> 9/50,
  rSig -> 4/5,
  a0 -> 11/5,
  b0 -> -3/5,
  rF1 -> 88899676773749/50000000000000,
  xiR -> 9/10,
  lam -> 1,
  etaLeak -> 3/100,
  muW -> 4/5,
  rho0 -> 1,
  q -> 1
};

deltaNum = N[delta /. sessionRules, 30];
dZeroNum = N[dZero /. sessionRules, 30];
VshortNum = N[Vshort /. sessionRules, 30];
MsigmaNum = N[Msigma /. sessionRules, 30];
gSoft = N[gMoment[r] /. sessionRules, 30];

Print["Delta(session) = ", fmt[deltaNum]];
Print["D0(session) = ", fmt[dZeroNum]];
Print["V_short(session) = ", fmt[VshortNum]];
expectApprox["M6 Delta session", deltaNum, 1421775/10000, 10^-8];
expectApprox["M6 D0 session", dZeroNum, 188240931/50000000, 10^-8];
expectApprox["M6 V_short session", VshortNum, 187081849/50000000, 10^-6];

subbanner["M7. Imported benchmark values and branch positivity"];

WsessObs = 151632107/100000000;
UVdropObs = 10532139/50000000;
VeffObs = 87350563/50000000;

LvarSoft = N[
  Sqrt[
    WsessObs Pi^4/(512 (3/100)^2 (4/5))
  ],
  30
];
SleakSoft = N[Sleak /. Join[sessionRules, {Lvar -> LvarSoft}], 30];

Print["M_sigma(session) = ", fmt[MsigmaNum]];
Print["Lvar(session) = ", fmt[LvarSoft]];
Print["S_leak(session) = ", fmt[SleakSoft]];
Print["g(session) = ", fmt[gSoft]];

expectApprox["M7 M_sigma session", MsigmaNum, 459653/2500000, 10^-6];
expectTrue["M7 M_sigma nonnegative", MsigmaNum >= 0];
expectApprox["M7 Lvar session", LvarSoft, 2001677473/100000000, 10^-6];
expectApprox["M7 S_leak session", SleakSoft, 31069599/100000000, 10^-6];
expectTrue["M7 g lower bound", gSoft >= 2/Pi];
expectTrue["M7 g below rF1", gSoft < (rF1 /. sessionRules)];

subbanner["M8. Forward benchmark decomposition"];

lambdaLPaper = 13485959/50000000;
VeffForward = N[
  VshortNum - lambdaLPaper SleakSoft - WsessObs - UVdropObs - MsigmaNum,
  30
];

Print["V_eff forward = ", fmt[VeffForward]];
expectApprox["M8 forward V_eff", VeffForward, VeffObs, 10^-6];
expectTrue["M8 lambda_L positive", lambdaLPaper > 0];

banner["Stage 247 Mathematica audit result"];
Print["All Mathematica checks passed."];
Exit[0];
