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

banner["STAGE 046 — PARENT THRESHOLDS"];

Clear[rhoStar, csStarSq, gPhi, kX, tX, ell, nSS, nPP, oSP, c2, gFail, gSuff, peReq, delta0, deltaInf, kappa, m];
$Assumptions =
  Element[{rhoStar, csStarSq, gPhi, kX, tX, ell, nSS, nPP, oSP, c2, gFail, gSuff, peReq, delta0, deltaInf, kappa, m}, Reals] &&
  rhoStar > 0 && csStarSq > 0 && gPhi > 0 && kX > 0 && tX > 0 && ell > 0 &&
  nSS > 0 && nPP > 0 && oSP > 0 && c2 > 0 && gFail > 0 && gSuff > 0 &&
  peReq > 0 && delta0 > 0 && deltaInf > 0 && kappa > 0 && m > 0;

gMicro = FullSimplify[rhoStar*gPhi^2*oSP^2/(m*csStarSq*kX*nSS), Assumptions -> $Assumptions];
gFailSq = FullSimplify[m*csStarSq*kX*nSS*gFail/(rhoStar*oSP^2), Assumptions -> $Assumptions];
gSuffSq = FullSimplify[m*csStarSq*kX*nSS*gSuff/(rhoStar*oSP^2), Assumptions -> $Assumptions];

Print["g_(phi,fail)^2 = ", fmt[gFailSq]];
Print["g_(phi,suff)^2 = ", fmt[gSuffSq]];
c2Def = FullSimplify[oSP^2/(nSS*nPP), Assumptions -> $Assumptions];
expectZero[
  "coherence substitution in g_fail^2",
  gFailSq - FullSimplify[m*csStarSq*kX*gFail/(rhoStar*nPP*c2Def), Assumptions -> $Assumptions]
];
expectZero[
  "coherence substitution in g_suff^2",
  gSuffSq - FullSimplify[m*csStarSq*kX*gSuff/(rhoStar*nPP*c2Def), Assumptions -> $Assumptions]
];

cFailSq = FullSimplify[m*csStarSq*kX*gFail/(rhoStar*gPhi^2*nPP), Assumptions -> $Assumptions];
cSuffSq = FullSimplify[m*csStarSq*kX*gSuff/(rhoStar*gPhi^2*nPP), Assumptions -> $Assumptions];
Print["C_fail^2 = ", fmt[cFailSq]];
Print["C_suff^2 = ", fmt[cSuffSq]];
expectZero["C_suff^2/C_fail^2 - G_suff/G_fail", FullSimplify[cSuffSq/cFailSq, Assumptions -> $Assumptions] - gSuff/gFail];

gMax = FullSimplify[rhoStar*gPhi^2*nPP/(m*csStarSq*kX), Assumptions -> $Assumptions];
expectZero["G_micro - G_max C^2", FullSimplify[gMicro /. oSP^2 -> c2*nSS*nPP, Assumptions -> $Assumptions] - gMax*c2];

gFailSub = peReq/(kappa*deltaInf);
gSuffSub = peReq/(kappa*delta0);
gFailSubKappa = FullSimplify[gFailSub /. kappa -> kX*ell^2/tX, Assumptions -> $Assumptions];
gSuffSubKappa = FullSimplify[gSuffSub /. kappa -> kX*ell^2/tX, Assumptions -> $Assumptions];
expectZero[
  "KX*g_fail threshold with kappa inserted",
  FullSimplify[gFailSq /. gFail -> gFailSubKappa, Assumptions -> $Assumptions] -
    m*csStarSq*tX*nSS*peReq/(rhoStar*oSP^2*ell^2*deltaInf)
];
expectZero[
  "KX*g_suff threshold with kappa inserted",
  FullSimplify[gSuffSq /. gSuff -> gSuffSubKappa, Assumptions -> $Assumptions] -
    m*csStarSq*tX*nSS*peReq/(rhoStar*oSP^2*ell^2*delta0)
];
expectZero[
  "coherence-form fail threshold with kappa inserted",
  FullSimplify[m*csStarSq*kX*gFailSubKappa/(rhoStar*nPP*c2Def), Assumptions -> $Assumptions] -
    FullSimplify[m*csStarSq*tX*peReq/(rhoStar*nPP*c2Def*ell^2*deltaInf), Assumptions -> $Assumptions]
];
expectZero[
  "coherence-form suff threshold with kappa inserted",
  FullSimplify[m*csStarSq*kX*gSuffSubKappa/(rhoStar*nPP*c2Def), Assumptions -> $Assumptions] -
    FullSimplify[m*csStarSq*tX*peReq/(rhoStar*nPP*c2Def*ell^2*delta0), Assumptions -> $Assumptions]
];

Print[""];
Print["Stage 046 Mathematica audit passed."];

Exit[0];
