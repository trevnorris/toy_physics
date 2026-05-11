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

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 022 — SPLIT-U SECTOR"];

Clear[kEtaEff, kU, kWEff, muEta, muW, cEtaU, cEtaW, cUW, tw, ell, gNewton, cLight, cS, aScale, delta0, deltaU, epsEta, epsW, rho0, zW, lambda];
$Assumptions =
  Element[{kEtaEff, kU, kWEff, muEta, muW, cEtaU, cEtaW, cUW, tw, ell, gNewton, cLight, cS, aScale, delta0, deltaU, epsEta, epsW, rho0, zW, lambda}, Reals] &&
  kEtaEff > 0 && kU > 0 && kWEff > 0 && muEta > 0 && muW > 0 && tw > 0 &&
  ell > 0 && gNewton > 0 && cLight > 0 && cS > 0 && aScale > 0 &&
  delta0 > 0 && deltaU > 0 && epsEta > 0 && epsW > 0 && rho0 > 0 &&
  zW > 0 && lambda > 0;

kappa0 = 2 Sqrt[2]/Pi;
kappa1 = -4/(3 Pi);
sigma = FullSimplify[kappa0^2 + kappa1^2, Assumptions -> $Assumptions];
lambda0 = FullSimplify[kappa1^2/kappa0^2, Assumptions -> $Assumptions];

Print["kappa0 = ", fmt[kappa0]];
Print["kappa1 = ", fmt[kappa1]];
Print["sigma = ", fmt[sigma]];
Print["lambda0 = ", fmt[lambda0]];

subbanner["1. Exact U-mode split and direct wall softening"];

kU1 = FullSimplify[kU (1 + deltaU), Assumptions -> $Assumptions];
a0 = FullSimplify[(kEtaEff - cEtaU^2/kU)/muEta, Assumptions -> $Assumptions];
a1 = FullSimplify[(kEtaEff (1 + delta0) - cEtaU^2/kU1)/muEta, Assumptions -> $Assumptions];

deltaSplit = FullSimplify[(delta0 + epsEta deltaU/(1 + deltaU))/(1 - epsEta), Assumptions -> $Assumptions];
a0Expected = FullSimplify[kEtaEff (1 - epsEta)/muEta, Assumptions -> $Assumptions];
a1Expected = FullSimplify[a0Expected (1 + deltaSplit), Assumptions -> $Assumptions];

Print["A0 = ", fmt[a0Expected]];
Print["A1 = ", fmt[a1Expected]];
Print["delta_split = ", fmt[deltaSplit]];
expectZero["A0 direct - expected", (a0 /. cEtaU^2 -> epsEta kU kEtaEff) - a0Expected];
expectZero["A1 direct - expected", (a1 /. cEtaU^2 -> epsEta kU kEtaEff) - a1Expected];

subbanner["2. Exact mixed blocking ratio with split U sector"];

sU = FullSimplify[kappa0^2/kU + kappa1^2/kU1, Assumptions -> $Assumptions];
epsWDirect = FullSimplify[cUW^2 sU/kWEff, Assumptions -> $Assumptions];
epsWSplit = FullSimplify[epsW (1 - (2/11) deltaU/(1 + deltaU)), Assumptions -> $Assumptions];

Print["S_U = ", fmt[sU]];
Print["eps_W_split = ", fmt[epsWSplit]];
expectZero["eps_W direct - split formula", (epsWDirect /. cUW^2 -> epsW kU kWEff/sigma) - epsWSplit];

subbanner["3. Mixed loading vector and exact direction-splitting invariant"];

gW = FullSimplify[cEtaW/Sqrt[muEta*muW], Assumptions -> $Assumptions];
z0 = FullSimplify[kappa0 gW (1 + rho0), Assumptions -> $Assumptions];
z1 = FullSimplify[kappa1 gW (1 + rho0/(1 + deltaU)), Assumptions -> $Assumptions];
rU = FullSimplify[(1 + rho0/(1 + deltaU))/(1 + rho0), Assumptions -> $Assumptions];

Print["z0 = ", fmt[z0]];
Print["z1 = ", fmt[z1]];
Print["R_U = ", fmt[rU]];
expectZero["z1/z0 - (kappa1/kappa0) R_U", z1/z0 - (kappa1/kappa0) rU];

dDir = FullSimplify[kappa0 z1 - kappa1 z0, Assumptions -> $Assumptions];
dDirExpected = FullSimplify[-kappa0 kappa1 gW rho0 deltaU/(1 + deltaU), Assumptions -> $Assumptions];
Print["D_dir = ", fmt[dDir]];
expectZero["direction-splitting invariant", dDir - dDirExpected];

subbanner["4. Split-U continuum placement map"];

mMixSplit = FullSimplify[8 zW (1 + rho0)^2/(Pi^2 (1 - epsEta) (1 - epsWSplit)), Assumptions -> $Assumptions];
rTargetSplit = FullSimplify[lambda (1 - epsEta) (1 - epsWSplit)^2/(zW (1 + rho0)^2), Assumptions -> $Assumptions];
product = FullSimplify[mMixSplit rTargetSplit, Assumptions -> $Assumptions];

Print["M_mix^(split U) = ", fmt[mMixSplit]];
Print["R_target^(split U) = ", fmt[rTargetSplit]];
Print["product = ", fmt[product]];
expectZero["product law", product - 8 lambda (1 - epsWSplit)/Pi^2];

subbanner["5. Small-splitting expansions"];

deltaSplitSeries = FullSimplify[Normal[Series[deltaSplit, {deltaU, 0, 1}]], Assumptions -> $Assumptions];
epsWSeries = FullSimplify[Normal[Series[epsWSplit, {deltaU, 0, 1}]], Assumptions -> $Assumptions];
rUSeries = FullSimplify[Normal[Series[rU, {deltaU, 0, 1}]], Assumptions -> $Assumptions];
m0 = FullSimplify[mMixSplit /. deltaU -> 0, Assumptions -> $Assumptions];
r0 = FullSimplify[rTargetSplit /. deltaU -> 0, Assumptions -> $Assumptions];
mRatio = FullSimplify[Normal[Series[mMixSplit/m0, {deltaU, 0, 1}]], Assumptions -> $Assumptions];
rRatio = FullSimplify[Normal[Series[rTargetSplit/r0, {deltaU, 0, 1}]], Assumptions -> $Assumptions];

Print["delta_split = ", fmt[deltaSplitSeries], " + O(deltaU^2)"];
Print["eps_W_split = ", fmt[epsWSeries], " + O(deltaU^2)"];
Print["R_U = ", fmt[rUSeries], " + O(deltaU^2)"];
Print["M_mix_split / M_mix_flat = ", fmt[mRatio], " + O(deltaU^2)"];
Print["R_target_split / R_target_flat = ", fmt[rRatio], " + O(deltaU^2)"];

Print[""];
Print["Stage 039 Mathematica audit passed."];

Exit[0];
