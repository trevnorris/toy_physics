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
  If[MatrixQ[expr],
    res = Map[FullSimplify[Together[Expand[#]], Assumptions -> $Assumptions] &, expr, {2}];
    Print[name, " = ", fmt[res]];
    If[!AllTrue[Flatten[res], # === 0 &], fail[name, res], pass[name]],
    res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
    Print[name, " = ", fmt[res]];
    If[TrueQ[res === 0], pass[name], fail[name, res]]
  ]
];

banner["STAGE 15.1 — EXACT FINITE-THROAT AXIAL INTEGRALS"];

Clear[l, s];
$Assumptions = Element[{l, s}, Reals] && l > 0;
u0 = 1/Sqrt[l];
u1 = Sqrt[2/l]*Cos[Pi*s/l];
f0 = Sqrt[2/l]*Sin[Pi*s/(2*l)];

kappa0 = FullSimplify[Integrate[u0*f0, {s, 0, l}], Assumptions -> $Assumptions];
kappa1 = FullSimplify[Integrate[u1*f0, {s, 0, l}], Assumptions -> $Assumptions];

Print["u0 . u0 = ", fmt[FullSimplify[Integrate[u0*u0, {s, 0, l}], Assumptions -> $Assumptions]]];
Print["u1 . u1 = ", fmt[FullSimplify[Integrate[u1*u1, {s, 0, l}], Assumptions -> $Assumptions]]];
Print["u0 . u1 = ", fmt[FullSimplify[Integrate[u0*u1, {s, 0, l}], Assumptions -> $Assumptions]]];
Print["f0 . f0 = ", fmt[FullSimplify[Integrate[f0*f0, {s, 0, l}], Assumptions -> $Assumptions]]];
Print["kappa0 = ", fmt[kappa0]];
Print["kappa1 = ", fmt[kappa1]];
Print["sigma = ", fmt[FullSimplify[kappa0^2 + kappa1^2, Assumptions -> $Assumptions]]];
Print["sigma/kappa0^2 = ", fmt[FullSimplify[(kappa0^2 + kappa1^2)/kappa0^2, Assumptions -> $Assumptions]]];

expectZero["u0.u0 - 1", Integrate[u0*u0, {s, 0, l}] - 1];
expectZero["u1.u1 - 1", Integrate[u1*u1, {s, 0, l}] - 1];
expectZero["u0.u1", Integrate[u0*u1, {s, 0, l}]];
expectZero["f0.f0 - 1", Integrate[f0*f0, {s, 0, l}] - 1];
expectZero["kappa0 - 2 sqrt(2)/pi", kappa0 - 2*Sqrt[2]/Pi];
expectZero["kappa1 + 4/(3 pi)", kappa1 + 4/(3*Pi)];
expectZero["sigma - 88/(9 pi^2)", kappa0^2 + kappa1^2 - 88/(9*Pi^2)];
expectZero["sigma/kappa0^2 - 11/9", (kappa0^2 + kappa1^2)/kappa0^2 - 11/9];

banner["STAGE 15.2 — LOCAL-KERNEL MODE REDUCTION"];
Clear[q0, q1, u0c, u1c, phi, wField, qStf, gU, gB, gW, gR, gQ];
$Assumptions = Element[{q0, q1, u0c, u1c, phi, wField, qStf, gU, gB, gW, gR, gQ}, Reals] && l > 0;

eta = q0*u0 + q1*u1;
uField = u0c*u0 + u1c*u1;
phiField = phi*f0;
wMode = wField*f0;

lEtaU = FullSimplify[gU*Integrate[eta*uField, {s, 0, l}], Assumptions -> $Assumptions];
lEtaPhi = FullSimplify[gB*Integrate[eta*phiField, {s, 0, l}], Assumptions -> $Assumptions];
lEtaW = FullSimplify[gW*Integrate[eta*wMode, {s, 0, l}], Assumptions -> $Assumptions];
lUW = FullSimplify[-gR*Integrate[uField*wMode, {s, 0, l}], Assumptions -> $Assumptions];
lSrc = FullSimplify[gQ*qStf*Integrate[eta*f0, {s, 0, l}], Assumptions -> $Assumptions];

Print["L_etaU = ", fmt[lEtaU]];
Print["L_etaphi = ", fmt[lEtaPhi]];
Print["L_etaW = ", fmt[lEtaW]];
Print["L_UW = ", fmt[lUW]];
Print["L_src = ", fmt[lSrc]];

expectZero["L_etaU - gU (q.U)", lEtaU - gU*(q0*u0c + q1*u1c)];
expectZero["L_etaphi - gB (v.q) phi", lEtaPhi - gB*(kappa0*q0 + kappa1*q1)*phi];
expectZero["L_etaW - gW (v.q) W", lEtaW - gW*(kappa0*q0 + kappa1*q1)*wField];
expectZero["L_UW + gR (v.U) W", lUW + gR*(kappa0*u0c + kappa1*u1c)*wField];
expectZero["L_src - gQ Q (v.q)", lSrc - gQ*qStf*(kappa0*q0 + kappa1*q1)];

banner["STAGE 15.3 — EXACT SCHUR-COMPLEMENT DECOMPOSITION"];
Clear[d0, d1, aPhi, aU, aW];
$Assumptions = Element[{d0, d1, aPhi, aU, aW, gU, gB, gW, gR}, Reals] && aU != 0 && aPhi != 0;

v = {{kappa0}, {kappa1}};
i2 = IdentityMatrix[2];
kInt = {
  {aU, 0, 0, -gR*kappa0},
  {0, aU, 0, -gR*kappa1},
  {0, 0, aPhi, 0},
  {-gR*kappa0, -gR*kappa1, 0, aW}
};
bMat = {
  {gU, 0, gB*kappa0, gW*kappa0},
  {0, gU, gB*kappa1, gW*kappa1}
};
sigmaMat = FullSimplify[bMat.Inverse[kInt].Transpose[bMat], Assumptions -> $Assumptions];
sigma = FullSimplify[(Transpose[v].v)[[1, 1]], Assumptions -> $Assumptions];
xi = FullSimplify[gU^2/aU, Assumptions -> $Assumptions];
alphaCoeff = FullSimplify[
  gB^2/aPhi + (aU*gW + gR*gU)^2/(aU*(aU*aW - gR^2*sigma)),
  Assumptions -> $Assumptions
];
sigmaTarget = FullSimplify[xi*i2 + alphaCoeff*(v.Transpose[v]), Assumptions -> $Assumptions];

Print["Sigma = ", fmt[sigmaMat]];
Print["Xi = ", fmt[xi]];
Print["alpha = ", fmt[alphaCoeff]];
expectZero["Sigma - [Xi I + alpha vv^T]", sigmaMat - sigmaTarget];

banner["STAGE 15.4 — NATURAL D/N SOURCE MAP"];
Clear[a, dK, alpha0, sigmaSym, deltaKappa, kappaProd, beta0, gConst, cs, radius, cSpeed, mhat];
$Assumptions = Element[{a, dK, alpha0, sigmaSym, deltaKappa, kappaProd, beta0, gConst, cs, radius, cSpeed, mhat}, Reals] &&
  a > 0 && dK > 0 && alpha0 >= 0 && sigmaSym > 0 && kappaProd > 0 && beta0 > 0 &&
  gConst > 0 && cs > 0 && radius > 0 && cSpeed > 0 && mhat > 0;

r = Sqrt[(dK + alpha0*deltaKappa)^2 + 4*alpha0^2*kappaProd];
lamMinus = FullSimplify[(2*a + dK - alpha0*sigmaSym - r)/2, Assumptions -> $Assumptions];
sMinus = FullSimplify[(sigmaSym + ((dK + alpha0*deltaKappa)*deltaKappa + 4*alpha0*kappaProd)/r)/2, Assumptions -> $Assumptions];
subsNat = {sigmaSym -> kappa0^2 + kappa1^2, deltaKappa -> kappa0^2 - kappa1^2, kappaProd -> kappa0^2*kappa1^2};
sMinusNat = FullSimplify[sMinus /. subsNat, Assumptions -> $Assumptions];
mhatSq = FullSimplify[sMinusNat/kappa0^2, Assumptions -> $Assumptions];

Print["mhat_-^2 = ", fmt[mhatSq]];
expectZero["mhat_-^2(alpha=0) - 1", (mhatSq /. alpha0 -> 0) - 1];
expectZero["limit_{alpha->oo} mhat_-^2 - 11/9", Limit[mhatSq, alpha0 -> Infinity] - 11/9];

banner["STAGE 15.5 — ELIMINATION OF THE ABSTRACT SOURCE-MAP FACTOR"];
p0Minus = FullSimplify[beta0*sMinus/lamMinus, Assumptions -> $Assumptions];
nProd = FullSimplify[(sMinus/kappa0^2)*p0Minus, Assumptions -> $Assumptions];
nProdTarget = FullSimplify[beta0*sMinus^2/(kappa0^2*lamMinus), Assumptions -> $Assumptions];
expectZero["mhat^2 P0_- - beta0 s^2/(kappa0^2 lambda_-)", nProd - nProdTarget];

Print["All Stage 15 checks passed."];

Exit[0];
