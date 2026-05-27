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

banner["STAGE 085 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS"];

Clear[aNorm, beta0, nQ, alphaReq, alphaMix, mhat, sMinus, lamMinus, epsBlk, rhoAlpha];
$Assumptions =
  Element[{aNorm, beta0, nQ, alphaReq, alphaMix, mhat, sMinus, lamMinus, epsBlk, rhoAlpha}, Reals] &&
  aNorm > 0 && beta0 > 0 && nQ > 0 && alphaReq > 0 && alphaMix > 0 && mhat > 0 && sMinus > 0 &&
  lamMinus > 0 && rhoAlpha > 0;

kappa0Sq = 8/Pi^2;
rTarget = FullSimplify[nQ*aNorm/(beta0*kappa0Sq), Assumptions -> $Assumptions];
gTr = FullSimplify[8*alphaReq/(Pi^2*aNorm), Assumptions -> $Assumptions];
mMix = FullSimplify[8*alphaMix/(Pi^2*aNorm), Assumptions -> $Assumptions];

piTr = FullSimplify[rTarget*gTr, Assumptions -> $Assumptions];
cMix = FullSimplify[rTarget*mMix, Assumptions -> $Assumptions];

Print["R_target = ", fmt[rTarget]];
Print["G_tr = ", fmt[gTr]];
Print["M_mix = ", fmt[mMix]];
Print["Pi_tr = ", fmt[piTr]];
Print["C_mix = ", fmt[cMix]];

expectZero["Pi_tr - N_Q alpha_req/beta0", piTr - nQ*alphaReq/beta0];
expectZero["C_mix - N_Q alpha_mix/beta0", cMix - nQ*alphaMix/beta0];
expectZero["Pi_tr/C_mix - alpha_req/alpha_mix", piTr/cMix - alphaReq/alphaMix];

nQSelected = FullSimplify[mhat^2*beta0*sMinus/lamMinus, Assumptions -> $Assumptions];
piSel = FullSimplify[piTr /. nQ -> nQSelected, Assumptions -> $Assumptions];
cSel = FullSimplify[cMix /. nQ -> nQSelected, Assumptions -> $Assumptions];

Print["N_Q^(target) = ", fmt[nQSelected]];
Print["Pi_sel = ", fmt[piSel]];
Print["C_sel = ", fmt[cSel]];
expectZero["Pi_sel - mhat^2 s_- alpha_req/lambda_-", piSel - mhat^2*sMinus*alphaReq/lamMinus];
expectZero["C_sel - mhat^2 s_- alpha_mix/lambda_-", cSel - mhat^2*sMinus*alphaMix/lamMinus];

zetaReq = FullSimplify[(piTr - cMix)/(cMix - epsBlk*(2*cMix - piTr)), Assumptions -> $Assumptions];
zetaExpected = FullSimplify[(alphaReq - alphaMix)/(alphaMix - epsBlk*(2*alphaMix - alphaReq)), Assumptions -> $Assumptions];
zetaRho = FullSimplify[(rhoAlpha - 1)/(1 - epsBlk*(2 - rhoAlpha)), Assumptions -> $Assumptions];

expectZero["zeta_req product form - loading form", zetaReq - zetaExpected];
expectZero["zeta_req loading form - rho_alpha form", (zetaExpected /. alphaReq -> rhoAlpha*alphaMix) - zetaRho];
expectZero["unblocked limit", (zetaExpected /. epsBlk -> 0) - (alphaReq/alphaMix - 1)];

Print[""];
Print["Stage 085 Mathematica audit passed."];

Exit[0];
