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

expectMatrixZero[name_String, expr_] := Module[{res, target},
  res = FullSimplify[Map[Together[Expand[#]] &, expr, {2}], Assumptions -> $Assumptions];
  target = ConstantArray[0, Dimensions[res]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === target], pass[name], fail[name, res]];
];

banner["STAGE 020 — CONTINUUM-KERNEL EXTRACTION"];

subbanner["1. Exact N/N and D/N modes"];

Clear[s, ell];
$Assumptions = Element[{s, ell}, Reals] && ell > 0;

u0 = 1/Sqrt[ell];
u1 = Sqrt[2/ell] Cos[Pi s/ell];
f0 = Sqrt[2/ell] Sin[Pi s/(2 ell)];

expectZero["u0 normalization - 1", Integrate[u0*u0, {s, 0, ell}] - 1];
expectZero["u1 normalization - 1", Integrate[u1*u1, {s, 0, ell}] - 1];
expectZero["f0 normalization - 1", Integrate[f0*f0, {s, 0, ell}] - 1];
expectZero["u0-u1 orthogonality", Integrate[u0*u1, {s, 0, ell}]];
expectZero["u0 N/N boundary slope @0", D[u0, s] /. s -> 0];
expectZero["u1 N/N boundary slope @0", D[u1, s] /. s -> 0];
expectZero["u1 N/N boundary slope @ell", D[u1, s] /. s -> ell];
expectZero["f0 D/N boundary value @0", f0 /. s -> 0];
expectZero["f0 D/N boundary slope @ell", D[f0, s] /. s -> ell];

kappa0 = FullSimplify[Integrate[u0*f0, {s, 0, ell}], Assumptions -> $Assumptions];
kappa1 = FullSimplify[Integrate[u1*f0, {s, 0, ell}], Assumptions -> $Assumptions];
sigma = FullSimplify[kappa0^2 + kappa1^2, Assumptions -> $Assumptions];

Print["kappa0 = ", fmt[kappa0]];
Print["kappa1 = ", fmt[kappa1]];
Print["sigma = ", fmt[sigma]];
expectZero["kappa0 - 2 Sqrt[2]/Pi", kappa0 - 2 Sqrt[2]/Pi];
expectZero["kappa1 + 4/(3 Pi)", kappa1 + 4/(3 Pi)];
expectZero["sigma - 88/(9 Pi^2)", sigma - 88/(9 Pi^2)];

subbanner["2. Mass-normalized projected kernels"];

Clear[muEta, muU, muPhi, muW, kEta, tOmega, tw, kU, kPhi, tPhi, kW, tW];
$Assumptions =
  Element[{muEta, muU, muPhi, muW, kEta, tOmega, tw, kU, kPhi, tPhi, kW, tW, ell}, Reals] &&
  muEta > 0 && muU > 0 && muPhi > 0 && muW > 0 && kEta > 0 && tOmega > 0 &&
  tw > 0 && kU > 0 && kPhi > 0 && tPhi > 0 && kW > 0 && tW > 0 && ell > 0;

kEtaEff = FullSimplify[kEta + 6 tOmega, Assumptions -> $Assumptions];
kWEff = FullSimplify[kW + Pi^2 tW/(4 ell^2), Assumptions -> $Assumptions];

k0 = FullSimplify[kEtaEff/muEta, Assumptions -> $Assumptions];
deltaKAx = FullSimplify[Pi^2 tw/(muEta ell^2), Assumptions -> $Assumptions];
varpi2 = FullSimplify[(kPhi + Pi^2 tPhi/(4 ell^2))/muPhi, Assumptions -> $Assumptions];
omegaU2 = FullSimplify[kU/muU, Assumptions -> $Assumptions];
omegaW2 = FullSimplify[kWEff/muW, Assumptions -> $Assumptions];

Print["K0 = ", fmt[k0]];
Print["DeltaK_ax = ", fmt[deltaKAx]];
Print["varpi^2 = ", fmt[varpi2]];
Print["Omega_U^2 = ", fmt[omegaU2]];
Print["Omega_W^2 = ", fmt[omegaW2]];

subbanner["3. Schur-complement factorization"];

Clear[aU, aPhi, aW, gU, gB, gW, gR];
$Assumptions =
  Element[{aU, aPhi, aW, gU, gB, gW, gR}, Reals] && aU != 0 && aPhi != 0 &&
  aW != 0 && (aU*aW - gR^2 sigma) != 0;

bMat = {
  {aU, 0, 0, -gR*kappa0},
  {0, aU, 0, -gR*kappa1},
  {0, 0, aPhi, 0},
  {-gR*kappa0, -gR*kappa1, 0, aW}
};

cMat = {
  {gU, 0, gB*kappa0, gW*kappa0},
  {0, gU, gB*kappa1, gW*kappa1}
};

sigmaWall = FullSimplify[cMat . LinearSolve[bMat, Transpose[cMat]], Assumptions -> $Assumptions];
v = {{kappa0}, {kappa1}};
alphaSolved = FullSimplify[
  alpha /. First[Solve[sigmaWall[[1, 2]] == alpha*kappa0*kappa1, alpha]],
  Assumptions -> $Assumptions
];
xiSolved = FullSimplify[
  xi /. First[Solve[sigmaWall[[1, 1]] == xi + alphaSolved*kappa0^2, xi]],
  Assumptions -> $Assumptions
];

expectZero["Sigma_wall (2,2) consistency with ansatz", sigmaWall[[2, 2]] - xiSolved - alphaSolved*kappa1^2];
Print["Xi (recovered) = ", fmt[xiSolved]];
Print["alpha (recovered) = ", fmt[alphaSolved]];

subbanner["4. Continuum extraction of A, alpha_mix, beta0, M_mix, delta"];

Clear[cEtaU, cEtaPhi, cEtaW, cUW];
$Assumptions =
  Element[{muEta, muU, muW, kEta, tOmega, tw, kU, kW, tW, ell, cEtaU, cEtaPhi, cEtaW, cUW}, Reals] &&
  muEta > 0 && muU > 0 && muW > 0 && kEta > 0 && tOmega > 0 && tw > 0 &&
  kU > 0 && kW > 0 && tW > 0 && ell > 0;

gUCont = FullSimplify[cEtaU/Sqrt[muEta*muU], Assumptions -> $Assumptions];
gBCont = FullSimplify[cEtaPhi/Sqrt[muEta*muPhi], Assumptions -> $Assumptions];
gWCont = FullSimplify[cEtaW/Sqrt[muEta*muW], Assumptions -> $Assumptions];
gRCont = FullSimplify[cUW/Sqrt[muU*muW], Assumptions -> $Assumptions];

a = FullSimplify[k0 - gUCont^2/omegaU2, Assumptions -> $Assumptions];
delta0 = FullSimplify[omegaU2*omegaW2 - gRCont^2*sigma, Assumptions -> $Assumptions];
chi = FullSimplify[omegaU2*gWCont + gRCont*gUCont, Assumptions -> $Assumptions];
beta0 = FullSimplify[chi^2/delta0^2, Assumptions -> $Assumptions];
alphaMix = FullSimplify[chi^2/(omegaU2*delta0), Assumptions -> $Assumptions];
mMix = FullSimplify[8 alphaMix/(Pi^2 a), Assumptions -> $Assumptions];
delta = FullSimplify[deltaKAx/a, Assumptions -> $Assumptions];

aDerived = FullSimplify[Together[k0 - gUCont^2/omegaU2], Assumptions -> $Assumptions];
delta0Expected = FullSimplify[(kU*kWEff - cUW^2 sigma)/(muU*muW), Assumptions -> $Assumptions];
chiExpected = FullSimplify[(kU*cEtaW + cUW*cEtaU)/(muU*Sqrt[muEta*muW]), Assumptions -> $Assumptions];
beta0Expected = FullSimplify[
  (muW/muEta) (kU*cEtaW + cUW*cEtaU)^2/(kU*kWEff - cUW^2 sigma)^2,
  Assumptions -> $Assumptions
];
alphaMixExpected = FullSimplify[
  (kU*cEtaW + cUW*cEtaU)^2/(muEta*kU*(kU*kWEff - cUW^2 sigma)),
  Assumptions -> $Assumptions
];
mMixExpected = FullSimplify[
  8 (kU*cEtaW + cUW*cEtaU)^2/(Pi^2 (kU*kEtaEff - cEtaU^2) (kU*kWEff - cUW^2 sigma)),
  Assumptions -> $Assumptions
];
deltaDerived = FullSimplify[Together[deltaKAx/aDerived], Assumptions -> $Assumptions];

expectZero["A reduces to closed form", a - aDerived];
expectZero[
  "A numerator matches Schur form",
  FullSimplify[
    Numerator[Together[aDerived]]*(muEta*kU) - Denominator[Together[aDerived]]*(kU*kEtaEff - cEtaU^2),
    Assumptions -> $Assumptions
  ]
];
expectZero["Delta0 continuum formula", delta0 - delta0Expected];
expectZero["Chi continuum formula", chi - chiExpected];
expectZero["beta0 continuum formula", beta0 - beta0Expected];
expectZero["alpha_mix continuum formula", alphaMix - alphaMixExpected];
expectZero["M_mix continuum formula", mMix - mMixExpected];
expectZero[
  "delta numerator matches closed form",
  FullSimplify[
    Numerator[Together[deltaDerived]]*(ell^2*(kU*kEtaEff - cEtaU^2)) -
      Denominator[Together[deltaDerived]]*(kU*Pi^2*tw),
    Assumptions -> $Assumptions
  ]
];

Print["A = ", fmt[aDerived]];
Print["Delta0 = ", fmt[delta0Expected]];
Print["Chi = ", fmt[chiExpected]];
Print["beta0 = ", fmt[beta0Expected]];
Print["alpha_mix = ", fmt[alphaMixExpected]];
Print["M_mix = ", fmt[mMixExpected]];
Print["delta = ", fmt[deltaDerived]];

subbanner["5. Exact continuum stability surfaces"];

expectZero[
  "A - [K_U K_eta^(eff) - c_(etaU)^2]/(mu_eta K_U)",
  a - (kU*kEtaEff - cEtaU^2)/(muEta*kU)
];
expectZero[
  "Delta0 - [K_U K_W^(eff) - c_(UW)^2 sigma]/(mu_U mu_W)",
  delta0 - (kU*kWEff - cUW^2 sigma)/(muU*muW)
];

Print["Stable selected branch requires:"];
Print["  K_U K_eta^(eff) > c_(etaU)^2"];
Print["  K_U K_W^(eff)   > c_(UW)^2 sigma"];
Print["  sigma = ", fmt[sigma]];

Print[""];
Print["Stage 037 Mathematica audit passed."];

Exit[0];
