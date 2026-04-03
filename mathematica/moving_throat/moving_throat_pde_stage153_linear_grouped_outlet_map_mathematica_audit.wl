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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 153 — LINEAR GROUPED-P2 DIRECT OUTLET MAP"];

Clear[
  D0, dD0, dD2, dD4, N0, dN0, sigma, P0, dkappa, dgamma, eps,
  aD0, aD2, aN0, bD0, bD2, bN0
];

$Assumptions =
  Element[
    {
      D0, dD0, dD2, dD4, N0, dN0, sigma, P0, dkappa, dgamma, eps,
      aD0, aD2, aN0, bD0, bD2, bN0
    },
    Reals
  ] && D0 != 0 && N0 != 0 && sigma != 0 && sigma != 1;

u2 = 1/9;
u4 = 4/81;
D2 = -u2*D0;
D4 = -D0/27;

u2Full = -(D2 + eps*dD2)/(D0 + eps*dD0);
u4Full = ((D2 + eps*dD2)^2 - (D0 + eps*dD0)*(D4 + eps*dD4))/(D0 + eps*dD0)^2;
P0Full = (N0 + eps*dN0)/(D0 + eps*dD0);

du2 = FullSimplify[Coefficient[Normal[Series[u2Full, {eps, 0, 1}]], eps, 1], Assumptions -> $Assumptions];
du4 = FullSimplify[Coefficient[Normal[Series[u4Full, {eps, 0, 1}]], eps, 1], Assumptions -> $Assumptions];
dP0 = FullSimplify[
  Coefficient[Normal[Series[P0Full /. N0 -> P0*D0, {eps, 0, 1}]], eps, 1],
  Assumptions -> $Assumptions
];

banner["Linear grouped conservative/output transport"];
expectZero["delta u2 + (dD2 + dD0/9)/D0", du2 + (dD2 + dD0/9)/D0];
expectZero["delta u4 + (dD4 + 2 dD2/9 + 5 dD0/81)/D0", du4 + (dD4 + 2*dD2/9 + 5*dD0/81)/D0];
expectZero["delta P0 - (dN0 - P0 dD0)/D0", dP0 - (dN0 - P0*dD0)/D0];

du2Hyb = -sigma*dkappa/(3*(1 - sigma));
dP0OverP0Hyb = -9*sigma*dgamma/(1 - sigma);

dkappaFromdu2 = FullSimplify[
  dkappa /. First[Solve[du2sym == du2Hyb, dkappa, Reals]] /. du2sym -> du2,
  Assumptions -> $Assumptions
];
dgammaFromdP0 = FullSimplify[
  dgamma /. First[Solve[dP0sym/P0 == dP0OverP0Hyb, dgamma, Reals]] /. dP0sym -> dP0,
  Assumptions -> $Assumptions
];

banner["Direct outlet coefficients"];
expectZero[
  "delta kappa_W - 3(1-sigma)(dD2 + dD0/9)/(sigma D0)",
  dkappaFromdu2 - 3*(1 - sigma)*(dD2 + dD0/9)/(sigma*D0)
];
expectZero[
  "delta gamma_W + (1-sigma)(dN0 - P0 dD0)/(9 sigma N0)",
  (dgammaFromdP0 /. P0 -> N0/D0) + (1 - sigma)*(dN0 - (N0/D0)*dD0)/(9*sigma*N0)
];

du4FromHyb = -8*sigma*dkappa/(27*(1 - sigma));
du4FromKappa = FullSimplify[du4FromHyb /. dkappa -> dkappaFromdu2, Assumptions -> $Assumptions];

banner["Even one-parameter consistency"];
expectZero["delta u4 - (8/9) delta u2", du4FromKappa - (8/9)*du2];
relation = FullSimplify[
  dD4 /. First[Solve[du4 == (8/9)*du2, dD4, Reals]],
  Assumptions -> $Assumptions
];
expectZero["delta D4 - (2/3) delta D2 - dD0/27", relation - ((2/3)*dD2 + dD0/27)];

aKappa = FullSimplify[3*(1 - sigma)*(aD2 + aD0/9)/(sigma*D0), Assumptions -> $Assumptions];
bKappa = FullSimplify[3*(1 - sigma)*(bD2 + bD0/9)/(sigma*D0), Assumptions -> $Assumptions];
aGamma = FullSimplify[-(1 - sigma)*(aN0 - P0*aD0)/(9*sigma*N0), Assumptions -> $Assumptions];
bGamma = FullSimplify[-(1 - sigma)*(bN0 - P0*bD0)/(9*sigma*N0), Assumptions -> $Assumptions];
P0Ref = FullSimplify[N0/D0, Assumptions -> $Assumptions];

aKappaFromMap = FullSimplify[Expand[dkappaFromdu2 /. {dD2 -> aD2, dD0 -> aD0}], Assumptions -> $Assumptions];
bKappaFromMap = FullSimplify[Expand[dkappaFromdu2 /. {dD2 -> bD2, dD0 -> bD0}], Assumptions -> $Assumptions];
aGammaFromMap = FullSimplify[Expand[dgammaFromdP0 /. {dN0 -> aN0, dD0 -> aD0}], Assumptions -> $Assumptions];
bGammaFromMap = FullSimplify[Expand[dgammaFromdP0 /. {dN0 -> bN0, dD0 -> bD0}], Assumptions -> $Assumptions];

banner["Grouped trace/anomaly transport"];
Print["a_kappa = ", fmt[aKappa]];
Print["b_kappa = ", fmt[bKappa]];
Print["a_gamma = ", fmt[aGamma]];
Print["b_gamma = ", fmt[bGamma]];
Print["a_kappa from map = ", fmt[aKappaFromMap]];
Print["b_kappa from map = ", fmt[bKappaFromMap]];
Print["a_gamma from map = ", fmt[aGammaFromMap]];
Print["b_gamma from map = ", fmt[bGammaFromMap]];

expectZero["trace kappa coefficient", aKappaFromMap - aKappa];
expectZero["anomaly kappa coefficient", bKappaFromMap - bKappa];
expectZero[
  "trace gamma coefficient",
  FullSimplify[(aGammaFromMap /. P0 -> P0Ref) - (aGamma /. P0 -> P0Ref), Assumptions -> $Assumptions]
];
expectZero[
  "anomaly gamma coefficient",
  FullSimplify[(bGammaFromMap /. P0 -> P0Ref) - (bGamma /. P0 -> P0Ref), Assumptions -> $Assumptions]
];

Print[""];
Print["Carry-forward formulas:"];
Print["  delta kappa_W^(A) = 3(1-sigma_*) [delta D_(A,2) + delta D_(A,0)/9] / (sigma_* D0)"];
Print["  delta gamma_W^(A) = -(1-sigma_*) [delta N_(A,0) - P0 delta D_(A,0)] / (9 sigma_* N0)"];
Print["  direct even-pole consistency: delta D_(A,4) = (2/3) delta D_(A,2) + delta D_(A,0)/27"];
Print["  so the full linear grouped-anisotropy outlet problem collapses to the pair"];
Print["      K_A := delta D_(A,2) + delta D_(A,0)/9"];
Print["      G_A := delta N_(A,0) - P0 delta D_(A,0)"];

Exit[0];
