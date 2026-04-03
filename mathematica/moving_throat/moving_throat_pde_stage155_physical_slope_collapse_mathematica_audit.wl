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

banner["STAGE 155 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM"];

Clear[d0, n0, u2, dD0, dD2, dD4, dN0, eps, lam, sigma];
$Assumptions = Element[{d0, n0, u2, dD0, dD2, dD4, dN0, eps, lam, sigma}, Reals] &&
  d0 != 0 && n0 != 0 && sigma != 0;

d2 = -u2*d0;
p0 = n0/d0;

dA0 = d0 + eps*lam*dD0;
dA2 = d2 + eps*lam*dD2;
nA0 = n0 + eps*lam*dN0;

u2A = -dA2/dA0;
p0A = nA0/dA0;

deltaU2 = FullSimplify[(Normal[Series[u2A, {eps, 0, 1}]] - u2)/(eps*lam), Assumptions -> $Assumptions];
deltaP0 = FullSimplify[(Normal[Series[p0A, {eps, 0, 1}]] - p0)/(eps*lam), Assumptions -> $Assumptions];

Print["delta u_2^(A) = ", fmt[deltaU2]];
Print["delta P_0^(A) = ", fmt[deltaP0]];

kA = dD2 + dD0/9;
gA = dN0 - p0*dD0;

banner["Exact obstruction-pair rewrite"];
expectZero[
  "K_A + D0*delta_u2 - (1/9-u2)*dD0",
  kA + d0*deltaU2 - (1/9 - u2)*dD0
];
expectZero["G_A - D0*delta_P0", gA - d0*deltaP0];

banner["Canonical-even collapse at u2 = 1/9"];
expectZero[
  "K_A + D0*delta_u2 on canonical branch",
  (kA + d0*deltaU2) /. u2 -> 1/9
];

banner["Physical form of the hidden-even consistency relation"];
u2Star = 1/9;
u4Star = 4/81;
d2Star = -u2Star*d0;
d4Star = d0*(u2Star^2 - u4Star);

dA2Star = d2Star + eps*lam*dD2;
dA4Star = d4Star + eps*lam*dD4;

u2AStar = -dA2Star/dA0;
u4AStar = (dA2Star^2 - dA0*dA4Star)/dA0^2;

deltaU2Star = FullSimplify[(Normal[Series[u2AStar, {eps, 0, 1}]] - u2Star)/(eps*lam), Assumptions -> $Assumptions];
deltaU4Star = FullSimplify[(Normal[Series[u4AStar, {eps, 0, 1}]] - u4Star)/(eps*lam), Assumptions -> $Assumptions];

Print["delta u_2^(A) on canonical branch = ", fmt[deltaU2Star]];
Print["delta u_4^(A) on canonical branch = ", fmt[deltaU4Star]];

evenRelation = 2*dD2/3 + dD0/27;
expectZero[
  "delta u_4 - (8/9) delta u_2 under microscopic even-consistency relation",
  (deltaU4Star - 8*deltaU2Star/9) /. dD4 -> evenRelation
];

banner["Direct outlet formulas in physical grouped variables"];
deltaKappa = FullSimplify[3*(1 - sigma)*kA/(sigma*d0), Assumptions -> $Assumptions];
deltaGamma = FullSimplify[-(1 - sigma)*gA/(9*sigma*n0), Assumptions -> $Assumptions];

expectZero[
  "delta kappa_W^(A) + 3(1-sigma)/sigma * delta u_2^(A)",
  (deltaKappa + 3*(1 - sigma)*deltaU2/sigma) /. u2 -> 1/9
];
expectZero[
  "delta gamma_W^(A) + (1-sigma)/(9 sigma P0) * delta P_0^(A)",
  deltaGamma + (1 - sigma)*deltaP0/(9*sigma*p0)
];

deltaQ = FullSimplify[-9*sigma*deltaGamma/(1 - sigma), Assumptions -> $Assumptions];
expectZero["Delta_Q^(A) - delta P_0^(A)/P0", deltaQ - deltaP0/p0];

Print[""];
Print["Carry-forward formulas:"];
Print["  K_A = -D_0 delta u_2^(A) on the canonical even branch"];
Print["  G_A = D_0 delta P_0^(A) exactly"];
Print["  delta u_4^(A) = (8/9) delta u_2^(A) iff dD_4 = (2/3) dD_2 + (1/27) dD_0"];
Print["  delta kappa_W^(A) = -3(1-sigma_*)/sigma_* * delta u_2^(A)"];
Print["  delta gamma_W^(A) = -(1-sigma_*)/(9 sigma_*) * delta P_0^(A)/P_0"];
Print["  Delta_Q^(A) = delta P_0^(A)/P_0 on the even-preserving branch"];

Print[""];
Print["Stage 155 Mathematica audit passed."];

Exit[0];
