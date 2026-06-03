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
stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[stripCE[expr]]], Assumptions -> $Assumptions];
  res = FullSimplify[stripCE[res], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[stripCE[cond], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 240 - SELECTED-BRANCH LOADING RATIO FROM THE QUADRUPOLE PRECURSOR"];

Clear[
  omega, omegaQ, alphaReq, alphaMix, beta0, nqTarget, mhatMinus,
  sMinus, lambdaMinus, cMix, lambdaSupport, epsStar
];

$Assumptions =
  Element[
    {
      omegaQ, alphaReq, alphaMix, beta0, nqTarget, mhatMinus,
      sMinus, lambdaMinus, cMix, lambdaSupport, epsStar
    },
    Reals
  ] &&
  omegaQ > 0 && alphaReq > alphaMix > 0 && beta0 > 0 &&
  nqTarget > 0 && mhatMinus > 0 && sMinus > 0 && lambdaMinus > 0 &&
  cMix > 0 && lambdaSupport > 0 && 0 < epsStar < 1;

(* M1: product-ratio identity and spectral substitution. *)
piTrProduct = nqTarget*alphaReq/beta0;
cMixProduct = nqTarget*alphaMix/beta0;
rhoAlpha = FullSimplify[alphaReq/alphaMix, Assumptions -> $Assumptions];

expectZero[
  "M1 Pi_tr/C_mix = alpha_req/alpha_mix",
  piTrProduct/cMixProduct - rhoAlpha
];

nqSpectral = mhatMinus^2*beta0*sMinus/lambdaMinus;
expectZero[
  "M1 spectral substitution leaves Pi_tr/C_mix unchanged",
  (piTrProduct/cMixProduct /. nqTarget -> nqSpectral) - rhoAlpha
];

(* M2: contact-plus-pole compiler from an Apart/Limit extraction. *)
poleFactor = 1 - omega^2/omegaQ^2;
ySupport = Together[
  alphaMix/alphaReq + ((alphaReq - alphaMix)/alphaReq)/(1 - omega^2/omegaQ^2)
];
yApart = Apart[ySupport, omega];

poleProbe = poleFactor*yApart;
If[FreeQ[poleProbe, omegaQ], fail["M2 pole extraction path sees omegaQ", poleProbe]];
c1Extracted = FullSimplify[
  stripCE[Limit[poleProbe, omega -> omegaQ]],
  Assumptions -> $Assumptions
];

contactProbe = yApart - c1Extracted/poleFactor;
If[FreeQ[contactProbe, omegaQ], fail["M2 contact extraction path sees omegaQ", contactProbe]];
c0Extracted = FullSimplify[
  stripCE[Limit[contactProbe, omega -> 0]],
  Assumptions -> $Assumptions
];

Print["M2 ySupport = ", fmt[ySupport]];
Print["M2 yApart = ", fmt[yApart]];
Print["M2 c0Extracted = ", fmt[c0Extracted]];
Print["M2 c1Extracted = ", fmt[c1Extracted]];

expectZero["M2 extracted c0 = 1/rho_alpha", c0Extracted - 1/rhoAlpha];
expectZero["M2 extracted c1 = (rho_alpha - 1)/rho_alpha",
           c1Extracted - (rhoAlpha - 1)/rhoAlpha];
expectZero["M2 c0 + c1 = 1", c0Extracted + c1Extracted - 1];
expectZero["M2 inverse rho_alpha = 1/c0", 1/c0Extracted - rhoAlpha];
expectZero["M2 inverse rho_alpha = 1/(1-c1)", 1/(1 - c1Extracted) - rhoAlpha];

zetaReq = FullSimplify[c1Extracted/c0Extracted, Assumptions -> $Assumptions];
expectZero["M2 zeta_req = c1/c0 = rho_alpha - 1", zetaReq - (rhoAlpha - 1)];
expectZero["M2 Omega_Q independence of extracted c0", D[c0Extracted, omegaQ]];
expectZero["M2 Omega_Q independence of extracted c1", D[c1Extracted, omegaQ]];

(* M3: minimal isotropic specialization. *)
c0Minimal = 3/4;
c1Minimal = 1/4;
rhoMinimalFromC0 = FullSimplify[1/c0Minimal];
rhoMinimalFromC1 = FullSimplify[1/(1 - c1Minimal)];
zetaMinimal = FullSimplify[c1Minimal/c0Minimal];

expectZero["M3 rho_alpha from c0 = 3/4 is 4/3", rhoMinimalFromC0 - 4/3];
expectZero["M3 rho_alpha from c1 = 1/4 is 4/3", rhoMinimalFromC1 - 4/3];
expectZero["M3 zeta_req = 1/3", zetaMinimal - 1/3];

(* M4: selected demand product, derived from the M3 loading ratio. *)
rhoMinimal = rhoMinimalFromC0;
piTrMinimal = FullSimplify[rhoMinimal*cMix, Assumptions -> $Assumptions];
sReqMinimal = FullSimplify[piTrMinimal/cMix, Assumptions -> $Assumptions];

expectZero["M4 Pi_tr = (4/3) C_mix on the minimal-isotropic branch",
           piTrMinimal - (4/3)*cMix];
expectZero["M4 S_req = Pi_tr/C_mix = 4/3", sReqMinimal - 4/3];

(* M5: support-selector reduction. *)
cMixSelector = 8*lambdaSupport*(1 - epsStar)/Pi^2;
piTrSelector = FullSimplify[piTrMinimal /. cMix -> cMixSelector, Assumptions -> $Assumptions];
varrhoFromSelector = FullSimplify[
  Pi^2*piTrSelector/(16*lambdaSupport),
  Assumptions -> $Assumptions
];

expectZero["M5 varrho = 2(1-eps_star)/3",
           varrhoFromSelector - 2*(1 - epsStar)/3];

(* M6: regime classification. *)
expectTrue["M6 regime ratio 1 < 4/3 < 2",
           Resolve[1 < sReqMinimal < 2, Reals]];
expectTrue["M6 C_mix < Pi_tr < 2 C_mix",
           cMix < piTrMinimal < 2*cMix];

Print[""];
Print["Stage 240 Mathematica audit passed."];

Exit[0];
