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

banner["STAGE 043 — CONTINUUM EXTRACTION OF THE ACTUAL SUPPORT DIRECTION"];

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
rPhi = FullSimplify[(y[[2]]/kappa1) / (y[[1]]/kappa0), Assumptions -> $Assumptions];
rPhiExpected = FullSimplify[(1 + sigma0/(1 + deltaU))/(1 + sigma0), Assumptions -> $Assumptions];
dPhi = FullSimplify[Det[{{kappa0, kappa1}, {y0, y1}}], Assumptions -> $Assumptions];
dPhiExpected = FullSimplify[-kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), Assumptions -> $Assumptions];

Print["y = ", fmt[y]];
Print["R_phi = ", fmt[rPhi]];
Print["D_phi = ", fmt[dPhi]];
expectZero["R_phi - expected", rPhi - rPhiExpected];
expectZero["D_phi - expected", dPhi - dPhiExpected];

subbanner["2. Exact support pole shift and split support-blocking ratio"];

sU = FullSimplify[v.dU.v, Assumptions -> $Assumptions];
sUSub = FullSimplify[
  sU /. {kappa1^2 -> (2/11) sigma, kappa0^2 -> (9/11) sigma},
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
sUEndpointZero = FullSimplify[(sU /. {kappa1^2 -> (2/11) sigma, kappa0^2 -> (9/11) sigma}) /. deltaU -> 0, Assumptions -> $Assumptions];
sUEndpointZeroExpected = sigma/kU;
sUEndpointInf = FullSimplify[Limit[sU /. {kappa1^2 -> (2/11) sigma, kappa0^2 -> (9/11) sigma}, deltaU -> Infinity], Assumptions -> $Assumptions];
sUEndpointInfExpected = (9/11) sigma/kU;
Print["v.D_U.v at deltaU=0 = ", fmt[sUEndpointZero]];
Print["v.D_U.v as deltaU->Infinity = ", fmt[sUEndpointInf]];
expectZero["overlap endpoint deltaU=0", sUEndpointZero - sUEndpointZeroExpected];
expectZero["overlap endpoint deltaU->Infinity", sUEndpointInf - sUEndpointInfExpected];
expectZero["support overlap contraction", sUSub - sUExpected];
expectZero["A_phi^(eff) - expected", aPhiEff - aPhiEffExpected];

aPhiEffMin = FullSimplify[Limit[aPhiEff, deltaU -> 0], Assumptions -> $Assumptions];
aPhiEffMinExpected = FullSimplify[kPhiEff - cUphi^2 sigma/kU, Assumptions -> $Assumptions];
Print["A_phi^(eff) at deltaU=0 = ", fmt[aPhiEffMin]];
expectZero["A_phi^(eff) at deltaU=0 (minimal)", aPhiEffMin - aPhiEffMinExpected];

overlapRatio = FullSimplify[(kPhiEff - aPhiEff)/(kPhiEff - aPhiEffMin), Assumptions -> $Assumptions];
overlapRatioExpected = FullSimplify[1 - (2/11) deltaU/(1 + deltaU), Assumptions -> $Assumptions];
Print["split-vs-minimal overlap ratio = ", fmt[overlapRatio]];
expectZero["split-vs-minimal overlap ratio", overlapRatio - overlapRatioExpected];

subbanner["3. Exact physical support baseline"];

epsPhiSplitPhys = FullSimplify[(epsPhiSplit /. epsPhi -> cUphi^2 sigma/(kU kPhiEff)), Assumptions -> $Assumptions];
mSuppCont = FullSimplify[
  (kappa0^2 cB^2 (1 + cEtaU cUphi/(kU cB))^2/(muEta muPhi))/
    ((kEtaEff (1 - epsEta)/muEta) (kPhiEff (1 - epsPhiSplitPhys)/muPhi)),
  Assumptions -> $Assumptions
];

(* F3: M_supp must not depend on muEta or muPhi (they must cancel). *)
expectZero["M_supp independent of muEta", FullSimplify[D[mSuppCont, muEta], Assumptions -> $Assumptions]];
expectZero["M_supp independent of muPhi", FullSimplify[D[mSuppCont, muPhi], Assumptions -> $Assumptions]];

(* F3: structural-form check with a free baseline symbol bBaseline. *)
Clear[bBaseline];
$Assumptions = $Assumptions && bBaseline > 0;
mSuppContInB = FullSimplify[mSuppCont /. kappa0^2 -> bBaseline, Assumptions -> $Assumptions];
mSuppStructExpected = FullSimplify[
  bBaseline (cB^2/(kEtaEff kPhiEff)) (1 + cEtaU cUphi/(kU cB))^2/
    ((1 - epsEta) (1 - epsPhiSplitPhys)),
  Assumptions -> $Assumptions
];
Print["M_supp structural form (free baseline) = ", fmt[mSuppContInB]];
expectZero["M_supp structural form (free baseline)", mSuppContInB - mSuppStructExpected];

(* F3: baseline value identification, isolated from the structural check. *)
mSuppContEval = FullSimplify[mSuppContInB /. bBaseline -> 8/Pi^2, Assumptions -> $Assumptions];
mSuppExpected = FullSimplify[mSuppStructExpected /. bBaseline -> 8/Pi^2, Assumptions -> $Assumptions];
Print["M_supp at baseline B = 8/Pi^2 = ", fmt[mSuppContEval]];
expectZero["M_supp at baseline B = 8/Pi^2", mSuppContEval - mSuppExpected];

subbanner["4. Exact tracking condition relative to the mixed vector"];

z0 = FullSimplify[kappa0 (gW + gU gR/kU), Assumptions -> $Assumptions];
z1 = FullSimplify[kappa1 (gW + gU gR/(kU (1 + deltaU))), Assumptions -> $Assumptions];
dPhiZ = FullSimplify[Det[{{y0, y1}, {z0, z1}}], Assumptions -> $Assumptions];
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
mismatchLeading = FullSimplify[Series[mismatch, {deltaU, 0, 1}] // Normal, Assumptions -> $Assumptions];
mismatchLeadingExpected = FullSimplify[deltaU (rho0 - sigma0)/((1 + rho0) (1 + sigma0)), Assumptions -> $Assumptions];
Print["mismatch leading in deltaU = ", fmt[mismatchLeading]];
expectZero["mismatch leading-in-deltaU coefficient", mismatchLeading - mismatchLeadingExpected];

Print["R_phi - R_U = ", fmt[mismatch]];
expectZero["mismatch formula", mismatch - mismatchExpected];

(* F2 cross-check 2: derive the sign of (R_phi - R_U) at sigma0 < rho0 and sigma0 > rho0
   by direct numerical-symbolic limits, independent of the closed-form mismatch formula.
   This catches a hidden sign error in the closed-form expected. *)
(* At deltaU = 1, sigma0 = 0, rho0 = 1: rho0 > sigma0, mismatch should be positive. *)
mismatchAtTestPoint1 = FullSimplify[
  (rPhi - rU) /. {deltaU -> 1, gS -> 0, gW -> gU*gR/kU},
  Assumptions -> $Assumptions
];
(* expected: deltaU*(rho0-sigma0)/((1+deltaU)(1+sigma0)(1+rho0)) at (1, 0, 1) = 1*1/(2*1*2) = 1/4 *)
expectZero["mismatch sign at deltaU=1, sigma0=0, rho0=1", mismatchAtTestPoint1 - 1/4];
(* At deltaU = 1, sigma0 = 1, rho0 = 0: sigma0 > rho0, mismatch should be negative. *)
mismatchAtTestPoint2 = FullSimplify[
  (rPhi - rU) /. {deltaU -> 1, gU -> 1, kU -> 1, gB -> 1, gS -> 1, gR -> 0, gW -> 1},
  Assumptions -> $Assumptions
];
(* expected: 1*(-1)/(2*2*1) = -1/4 *)
expectZero["mismatch sign at deltaU=1, sigma0=1, rho0=0", mismatchAtTestPoint2 - (-1/4)];
(* At sigma0 = rho0 (tracking): mismatch must vanish for ANY deltaU. *)
mismatchAtTracking = FullSimplify[(rPhi - rU) /. {gS -> gB*gR/gW}, Assumptions -> $Assumptions];
expectZero["mismatch vanishes at tracking sigma0=rho0", mismatchAtTracking];

Print[""];
Print["Stage 043 Mathematica audit passed."];

Exit[0];
