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
rShift = Expand[rFun /. g -> gStar + dg];
rLin = Normal[Series[rShift, {dg, 0, 1}]];
rExpected = Expand[1/4 - dg/Sqrt[1 + r^2]];
expectZero["linear delta R law", rLin - rExpected];

Clear[sigma0, dSigma0, rStar, dR];
$Assumptions = Element[{sigma0, dSigma0, rStar, dR}, Reals];

mQ = -(sigma0 + dSigma0)*(rStar + dR);
mQLin = Expand[mQ] /. dSigma0*dR -> 0;
mQ0 = -sigma0*rStar;
expectZero["delta Mq law", (mQLin - mQ0) - (-rStar*dSigma0 - sigma0*dR)];

Clear[sStar, dS];
$Assumptions = Element[{sigma0, dSigma0, rStar, dR, sStar, dS}, Reals];

piExpr = (sigma0 + dSigma0)*(1 - (rStar + dR)*(sStar + dS));
piLin = Expand[piExpr] /. {dSigma0*dR -> 0, dSigma0*dS -> 0, dR*dS -> 0};
pi0 = sigma0*(1 - rStar*sStar);
dPiExpected = (1 - rStar*sStar)*dSigma0 - sigma0*(rStar*dS + sStar*dR);
expectZero["delta Pi law", (piLin - pi0) - dPiExpected];

Clear[dgSym, rSym];
$Assumptions = Element[{sigma0, dSigma0, sStar, dS, dgSym, rSym}, Reals] && rSym > 0;
dRFromDg = -dgSym/Sqrt[1 + rSym^2];

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
chiLin = Expand[Normal[Series[chi, {eps, 0, 1}]]];
chiExpected = 1 + eps*(5*b + a0/3 + 9*a5);
expectZero["linear Delta_Q law", chiLin - chiExpected];

banner["Numerical coefficients at the renormalized canonical point"];

rF1 = SetPrecision[1.77799353547498, 30];
sigma0Can = SetPrecision[4.651033550168876, 30];
sCan = SetPrecision[0.6703621156734617, 30];
tCan = SetPrecision[1.4467083664567624, 30];
sqrt1 = Sqrt[1 + rF1^2];

coefDRDG = N[-1/sqrt1, 20];
coefDMQDSigma = -1/4;
coefDMQDg = N[sigma0Can/sqrt1, 20];
coefDPiDSigma = N[1 - sCan/4, 20];
coefDPiDS = N[-sigma0Can/4, 20];
coefDPiDG = N[sigma0Can*sCan/sqrt1, 20];
coefDSigmaDT = N[(40/9)*tCan, 20];
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
