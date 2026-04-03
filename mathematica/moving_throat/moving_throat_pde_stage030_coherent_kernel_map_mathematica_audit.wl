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

banner["STAGE 030 — COHERENT KERNEL MAP"];

Clear[kEtaEff, kU, kWEff, kPhiEff, lamW, lamPhi, gamma, cEtaU, tw, ell, tU, gConst, cs, a, c, muW, zeta];
$Assumptions =
  Element[{kEtaEff, kU, kWEff, kPhiEff, lamW, lamPhi, gamma, cEtaU, tw, ell, tU, gConst, cs, a, c, muW, zeta}, Reals] &&
  kEtaEff > 0 && kU > 0 && kWEff > 0 && kPhiEff > 0 && lamW > 0 && lamPhi > 0 &&
  gamma > 0 && cEtaU > 0 && tw > 0 && ell > 0 && tU > 0 && gConst > 0 &&
  cs > 0 && a > 0 && c > 0 && muW > 0 && zeta > 0;

sigma = 88/(9*Pi^2);
rho0 = FullSimplify[(gamma*lamW)*cEtaU/(kU*lamW), Assumptions -> $Assumptions];
sigma0 = FullSimplify[cEtaU*(gamma*lamPhi)/(kU*lamPhi), Assumptions -> $Assumptions];
chi0 = FullSimplify[gamma*cEtaU/kU, Assumptions -> $Assumptions];

Print["rho_0 = ", fmt[rho0]];
Print["sigma_0 = ", fmt[sigma0]];
Print["chi_0 = ", fmt[chi0]];
expectZero["rho_0 - chi_0", rho0 - chi0];
expectZero["sigma_0 - chi_0", sigma0 - chi0];

epsEta = FullSimplify[cEtaU^2/(kU*kEtaEff), Assumptions -> $Assumptions];
epsW = FullSimplify[gamma^2*lamW^2*sigma/(kU*kWEff), Assumptions -> $Assumptions];
zetaDef = FullSimplify[lamPhi^2*kWEff/(lamW^2*kPhiEff), Assumptions -> $Assumptions];
zW = FullSimplify[lamW^2/(kEtaEff*kWEff), Assumptions -> $Assumptions];
delta0 = FullSimplify[Pi^2*tw/(ell^2*kEtaEff), Assumptions -> $Assumptions];
deltaU = FullSimplify[Pi^2*tU/(ell^2*kU), Assumptions -> $Assumptions];
lambdaScale = FullSimplify[27*Pi^2*gConst*cs^5*kWEff/(20*a^5*c^5*muW), Assumptions -> $Assumptions];
epsPhi = FullSimplify[gamma^2*lamPhi^2*sigma/(kU*kPhiEff), Assumptions -> $Assumptions];
zPhi = FullSimplify[lamPhi^2/(kEtaEff*kPhiEff), Assumptions -> $Assumptions];

Print["zeta_def = ", fmt[zetaDef]];
Print["eps_W = ", fmt[epsW]];
Print["eps_phi = ", fmt[epsPhi]];
Print["Z_W = ", fmt[zW]];
Print["Z_phi = ", fmt[zPhi]];
expectZero["eps_phi - zeta_def eps_W", epsPhi - zetaDef*epsW];
expectZero["Z_phi - zeta_def Z_W", zPhi - zetaDef*zW];

eps = FullSimplify[epsW*(1 - (2/11)*deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
rTr = FullSimplify[(1 + chi0/(1 + deltaU))/(1 + chi0), Assumptions -> $Assumptions];
delta = FullSimplify[(delta0 + epsEta*deltaU/(1 + deltaU))/(1 - epsEta), Assumptions -> $Assumptions];

Print["R_tr = ", fmt[rTr]];
Print["eps = ", fmt[eps]];
Print["delta = ", fmt[delta]];

mMix = FullSimplify[8*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - eps)), Assumptions -> $Assumptions];
mSupp = FullSimplify[8*zeta*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps)), Assumptions -> $Assumptions];
sEnhance = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps), Assumptions -> $Assumptions];
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
rTarget = FullSimplify[lambdaScale*(1 - epsEta)*(1 - eps)^2/(zW*(1 + chi0)^2), Assumptions -> $Assumptions];

Print["M_mix = ", fmt[mMix]];
Print["M_supp = ", fmt[mSupp]];
Print["S(zeta;eps) = ", fmt[sEnhance]];
Print["M_tr = ", fmt[mTr]];
Print["R_target = ", fmt[rTarget]];
expectZero["M_tr - M_mix S", mTr - mMix*sEnhance];

productExpected = FullSimplify[8*lambdaScale*(1 - eps)/Pi^2*sEnhance, Assumptions -> $Assumptions];
productActual = FullSimplify[rTarget*mTr, Assumptions -> $Assumptions];
Print["R_target M_tr = ", fmt[productActual]];
expectZero["product law", productActual - productExpected];
expectZero["dR_target/dzeta", D[rTarget, zeta]];
expectZero["dS/dzeta - (1-eps)/(1-zeta eps)^2", D[sEnhance, zeta] - (1 - eps)/(1 - zeta*eps)^2];
expectZero["S(zeta=0)-1", FullSimplify[sEnhance /. zeta -> 0, Assumptions -> $Assumptions] - 1];

Print[""];
Print["Stage 030 Mathematica audit passed."];

Exit[0];
