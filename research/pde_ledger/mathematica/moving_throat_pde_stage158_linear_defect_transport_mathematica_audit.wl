ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 158 — LINEAR DEFECT TRANSPORT AROUND THE RENORMALIZED CANONICAL POINT"];

Clear[g, r, dg];
$Assumptions = Element[{g, r, dg}, Reals];

gStar = r - Sqrt[1 + r^2]/2;
rFun = (g - r)^2/(1 + r^2);
rBase = FullSimplify[rFun /. g -> gStar, Assumptions -> $Assumptions];
rSlope = FullSimplify[D[rFun, g] /. g -> gStar, Assumptions -> $Assumptions];
rLin = Expand[rBase + rSlope*dg];
rExpected = Expand[1/4 - dg/Sqrt[1 + r^2]];
expectZero["linear delta R law", rLin - rExpected];

Clear[sigma0, dSigma0, rStar, dR, sigmaVar, rVar];
$Assumptions = Element[{sigma0, dSigma0, rStar, dR}, Reals];

mQBase = -sigmaVar*rVar;
mQ0 = mQBase /. {sigmaVar -> sigma0, rVar -> rStar};
mQLin = Expand[
  mQ0
  + (D[mQBase, sigmaVar] /. {sigmaVar -> sigma0, rVar -> rStar})*dSigma0
  + (D[mQBase, rVar] /. {sigmaVar -> sigma0, rVar -> rStar})*dR
];
expectZero["delta Mq law", (mQLin - mQ0) - (-rStar*dSigma0 - sigma0*dR)];

Clear[sStar, dS, piSigmaVar, piRVar, piSVar];
$Assumptions = Element[{sigma0, dSigma0, rStar, dR, sStar, dS}, Reals];

piBase = piSigmaVar*(1 - piRVar*piSVar);
pi0 = piBase /. {piSigmaVar -> sigma0, piRVar -> rStar, piSVar -> sStar};
piLin = Expand[
  pi0
  + (D[piBase, piSigmaVar] /. {piSigmaVar -> sigma0, piRVar -> rStar, piSVar -> sStar})*dSigma0
  + (D[piBase, piRVar] /. {piSigmaVar -> sigma0, piRVar -> rStar, piSVar -> sStar})*dR
  + (D[piBase, piSVar] /. {piSigmaVar -> sigma0, piRVar -> rStar, piSVar -> sStar})*dS
];
dPiExpected = (1 - rStar*sStar)*dSigma0 - sigma0*(rStar*dS + sStar*dR);
expectZero["delta Pi law", (piLin - pi0) - dPiExpected];

Clear[dgSym, rSym, gSym];
$Assumptions = Element[{sigma0, dSigma0, sStar, dS, dgSym, rSym}, Reals] && rSym > 0;
dRFromDg = FullSimplify[
  (D[(gSym - rSym)^2/(1 + rSym^2), gSym] /. gSym -> rSym - Sqrt[1 + rSym^2]/2)*dgSym,
  Assumptions -> $Assumptions
];

dMqComposed = -(1/4)*dSigma0 - sigma0*dRFromDg;
dMqBoxed = -(1/4)*dSigma0 + (sigma0/Sqrt[1 + rSym^2])*dgSym;
expectZero["composed delta Mq law", Expand[dMqComposed - dMqBoxed]];

dPiComposed = (1 - sStar/4)*dSigma0 - sigma0*((1/4)*dS + sStar*dRFromDg);
dPiBoxed = (1 - sStar/4)*dSigma0 - (sigma0/4)*dS + (sigma0*sStar/Sqrt[1 + rSym^2])*dgSym;
expectZero["composed delta Pi law", Expand[dPiComposed - dPiBoxed]];

Clear[eps, s, b, a0, a5];
$Assumptions = Element[{eps, s, b, a0, a5}, Reals];

sNorm = 1 + eps*s;
beta = 1 + eps*b;
sigmaZero = eps*a0;
sigmaFive = eps*a5;
chi = FullSimplify[3*(sNorm*beta^5 + 9*sigmaFive)/(3*sNorm - sigmaZero)];
chiBase = FullSimplify[chi /. eps -> 0, Assumptions -> $Assumptions];
chiSlope = FullSimplify[D[chi, eps] /. eps -> 0, Assumptions -> $Assumptions];
chiLin = Expand[chiBase + chiSlope*eps];
chiExpected = 1 + eps*(5*b + a0/3 + 9*a5);
expectZero["linear Delta_Q law", chiLin - chiExpected];

banner["Numerical coefficients at the renormalized canonical point"];

rF1 = SetPrecision[1.77799353547498, 30];
sigma0Can = SetPrecision[4.651033550168876, 30];
sCan = SetPrecision[0.6703621156734617, 30];
tCan = SetPrecision[1.4467083664567624, 30];

Clear[gNum, rNum, sigmaNum, rTransportNum, sNum, tNum];
rNumFun = (gNum - rNum)^2/(1 + rNum^2);
gNumStar = rNum - Sqrt[1 + rNum^2]/2;
rTransportValue = 1/4;
mQNum = -sigmaNum*rTransportNum;
piNum = sigmaNum*(1 - rTransportNum*sNum);
sigmaThatNum = (20/9)*tNum^2;
dRdGCan = FullSimplify[D[rNumFun, gNum] /. gNum -> gNumStar, Assumptions -> Element[rNum, Reals]];

coefDRDG = N[dRdGCan /. rNum -> rF1, 20];
coefDMQDSigma = D[mQNum, sigmaNum] /. {sigmaNum -> sigma0Can, rTransportNum -> rTransportValue};
coefDMQDg = N[
  (D[mQNum, rTransportNum] /. {sigmaNum -> sigma0Can, rTransportNum -> rTransportValue})*
    (dRdGCan /. rNum -> rF1),
  20
];
coefDPiDSigma = N[
  D[piNum, sigmaNum] /. {sigmaNum -> sigma0Can, rTransportNum -> rTransportValue, sNum -> sCan},
  20
];
coefDPiDS = N[
  D[piNum, sNum] /. {sigmaNum -> sigma0Can, rTransportNum -> rTransportValue, sNum -> sCan},
  20
];
coefDPiDG = N[
  (D[piNum, rTransportNum] /. {sigmaNum -> sigma0Can, rTransportNum -> rTransportValue, sNum -> sCan})*
    (dRdGCan /. rNum -> rF1),
  20
];
coefDSigmaDT = N[D[sigmaThatNum, tNum] /. tNum -> tCan, 20];
coefDPiDT = N[coefDPiDSigma*coefDSigmaDT, 20];

Print["dR/dg        = ", fmt[coefDRDG]];
Print["dMq/dSigma0  = ", fmt[coefDMQDSigma]];
Print["dMq/dg       = ", fmt[coefDMQDg]];
Print["dPi/dSigma0  = ", fmt[coefDPiDSigma]];
Print["dPi/dS       = ", fmt[coefDPiDS]];
Print["dPi/dg       = ", fmt[coefDPiDG]];
Print["dSigma0/dThat= ", fmt[coefDSigmaDT]];
Print["dPi/dThat    = ", fmt[coefDPiDT]];

Print[""];
Print["Carry-forward summary:"];
Print["  delta R  = -(delta g)/sqrt(1+r_F1^2) + O(delta g^2)"];
Print["  delta Mq = -(1/4) delta Sigma0 + Sigma0_can/sqrt(1+r_F1^2) delta g + O(2)"];
Print["  delta Pi = (1-S_can/4) delta Sigma0 - (Sigma0_can/4) delta S + Sigma0_can*S_can/sqrt(1+r_F1^2) delta g + O(2)"];
Print["  Delta_Q  = 5 b + a0/3 + 9 a5 + O(2)"];

Print[""];
Print["Stage 158 Mathematica audit passed."];

Exit[0];
