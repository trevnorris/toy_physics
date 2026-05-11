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

banner["STAGE 026 — CONTINUUM EXTRACTION OF THE ACTUAL SUPPORT DIRECTION"];

Clear[kU, kEtaEff, kPhiEff, deltaU, gU, gB, gS, gW, gR, kappa0, kappa1, muEta, muPhi, cEtaU, cB, cUphi, cW, cUW, epsEta, zPhi, sigma];
$Assumptions =
  Element[{kU, kEtaEff, kPhiEff, deltaU, gU, gB, gS, gW, gR, kappa0, kappa1, muEta, muPhi, cEtaU, cB, cUphi, cW, cUW, epsEta, zPhi, sigma}, Reals] &&
  kU > 0 && kEtaEff > 0 && kPhiEff > 0 && deltaU > 0 && muEta > 0 && muPhi > 0 &&
  sigma > 0 && gB != 0 && gW != 0 && kappa0 != 0 && kappa1 != 0;

subbanner["1. Exact effective support vector after eliminating the split U doublet"];

dU = DiagonalMatrix[{1/kU, 1/(kU (1 + deltaU))}];
v = {kappa0, kappa1};
y = FullSimplify[gB v + gU gS (dU.v), Assumptions -> $Assumptions];
y0 = FullSimplify[y[[1]], Assumptions -> $Assumptions];
y1 = FullSimplify[y[[2]], Assumptions -> $Assumptions];

sigma0 = FullSimplify[gU gS/(kU gB), Assumptions -> $Assumptions];
rPhi = FullSimplify[(y1/y0)/(kappa1/kappa0), Assumptions -> $Assumptions];
rPhiExpected = FullSimplify[(1 + sigma0/(1 + deltaU))/(1 + sigma0), Assumptions -> $Assumptions];
dPhi = FullSimplify[kappa0 y1 - kappa1 y0, Assumptions -> $Assumptions];
dPhiExpected = FullSimplify[-kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), Assumptions -> $Assumptions];

Print["y = ", fmt[y]];
Print["R_phi = ", fmt[rPhi]];
Print["D_phi = ", fmt[dPhi]];
expectZero["R_phi - expected", rPhi - rPhiExpected];
expectZero["D_phi - expected", dPhi - dPhiExpected];

subbanner["2. Exact support pole shift and split support-blocking ratio"];

sU = FullSimplify[v.dU.v, Assumptions -> $Assumptions];
sUSub = FullSimplify[
  sU /. {kappa1^2 -> (2/11) sigma, kappa0^2 -> sigma - (2/11) sigma},
  Assumptions -> $Assumptions
];
sUExpected = FullSimplify[(sigma/kU) (1 - (2/11) deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
epsPhi = Symbol["epsPhi"];
epsPhiSplit = FullSimplify[epsPhi (1 - (2/11) deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
aPhiEff = FullSimplify[kPhiEff - cUphi^2 sUExpected, Assumptions -> $Assumptions];
aPhiEffExpected = FullSimplify[
  (kPhiEff (1 - epsPhiSplit)) /. epsPhi -> cUphi^2 sigma/(kU kPhiEff),
  Assumptions -> $Assumptions
];

Print["v.D_U.v = ", fmt[sUSub]];
Print["A_phi^(eff) = ", fmt[aPhiEff]];
expectZero["support overlap contraction", sUSub - sUExpected];
expectZero["A_phi^(eff) - expected", aPhiEff - aPhiEffExpected];

subbanner["3. Exact physical support baseline"];

epsPhiSplitPhys = FullSimplify[(epsPhiSplit /. epsPhi -> cUphi^2 sigma/(kU kPhiEff)), Assumptions -> $Assumptions];
mSuppCont = FullSimplify[
  (kappa0^2 cB^2 (1 + cEtaU cUphi/(kU cB))^2/(muEta muPhi))/
    ((kEtaEff (1 - epsEta)/muEta) (kPhiEff (1 - epsPhiSplitPhys)/muPhi)),
  Assumptions -> $Assumptions
];
mSuppExpected = FullSimplify[
  (8/Pi^2) (cB^2/(kEtaEff kPhiEff)) (1 + cEtaU cUphi/(kU cB))^2/
    ((1 - epsEta) (1 - epsPhiSplitPhys)),
  Assumptions -> $Assumptions
];
mSuppContEval = FullSimplify[mSuppCont /. kappa0^2 -> 8/Pi^2, Assumptions -> $Assumptions];

Print["M_supp = ", fmt[mSuppContEval]];
expectZero["M_supp - expected", mSuppContEval - mSuppExpected];

subbanner["4. Exact tracking condition relative to the mixed vector"];

z0 = FullSimplify[kappa0 (gW + gU gR/kU), Assumptions -> $Assumptions];
z1 = FullSimplify[kappa1 (gW + gU gR/(kU (1 + deltaU))), Assumptions -> $Assumptions];
dPhiZ = FullSimplify[y0 z1 - y1 z0, Assumptions -> $Assumptions];
dPhiZExpected = FullSimplify[
  -deltaU gU kappa0 kappa1 (gB gR - gW gS)/(kU (1 + deltaU)),
  Assumptions -> $Assumptions
];
rho0 = FullSimplify[gU gR/(kU gW), Assumptions -> $Assumptions];
rU = FullSimplify[(1 + rho0/(1 + deltaU))/(1 + rho0), Assumptions -> $Assumptions];

Print["D_(phi z) = ", fmt[dPhiZ]];
expectZero["D_(phi z) - expected", dPhiZ - dPhiZExpected];
expectZero["tracking condition via g_B g_R = g_W g_S", (rPhiExpected - rU) /. gS -> gB gR/gW];

subbanner["5. Exact mismatch formula"];

mismatch = FullSimplify[rPhiExpected - rU, Assumptions -> $Assumptions];
mismatchExpected = FullSimplify[
  deltaU (rho0 - sigma0)/((1 + deltaU) (1 + rho0) (1 + sigma0)),
  Assumptions -> $Assumptions
];

Print["R_phi - R_U = ", fmt[mismatch]];
expectZero["mismatch formula", mismatch - mismatchExpected];

Print[""];
Print["Stage 043 Mathematica audit passed."];

Exit[0];
