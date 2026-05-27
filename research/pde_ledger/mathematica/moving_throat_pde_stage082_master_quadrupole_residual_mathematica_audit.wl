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

banner["STAGE 082 — MASTER QUADRUPOLE RESIDUAL"];

Clear[PiTr, cMix, epsBlk, zeta, zetaMinus, zetaPlus, zetaPhys, thetaW, upsilonW];
$Assumptions =
  cMix > 0 && thetaW > 0 && upsilonW > 0 &&
  Element[{PiTr, epsBlk, zeta, zetaMinus, zetaPlus, zetaPhys}, Reals];

(* Independent re-derivation: solve PiTr == cMix * qMap(zetaSym) for zetaSym,
   then rename to zeta. This forces Mathematica to find the inverse of qMap
   rather than restating the SymPy-side closed form. *)
qMap = FullSimplify[(1 + (1 - 2*epsBlk)*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
zetaSym = Unique["zetaSym"];
zetaReqSolList = Solve[PiTr == cMix*((1 + (1 - 2*epsBlk)*zetaSym)/(1 - epsBlk*zetaSym)), zetaSym];
zetaReq = FullSimplify[(zetaSym /. First[zetaReqSolList]) /. ConditionalExpression[x_, _] :> x,
                       Assumptions -> $Assumptions];

Print["zeta_req(Pi_tr,C_mix,eps_blk) = ", fmt[zetaReq]];
Print["Q(zeta;eps_blk) = ", fmt[qMap]];

expectZero[
  "inverse map zeta_req(C_mix*Q(zeta)) - zeta",
  FullSimplify[(zetaReq /. PiTr -> cMix*qMap) - zeta, Assumptions -> $Assumptions]
];

piSuff = FullSimplify[cMix*(qMap /. zeta -> zetaMinus), Assumptions -> $Assumptions];
piFail = FullSimplify[cMix*(qMap /. zeta -> zetaPlus), Assumptions -> $Assumptions];

Print["Pi_suff = ", fmt[piSuff]];
Print["Pi_fail = ", fmt[piFail]];

expectZero["zeta_req(Pi_suff) - zeta_-", FullSimplify[(zetaReq /. PiTr -> piSuff) - zetaMinus, Assumptions -> $Assumptions]];
expectZero["zeta_req(Pi_fail) - zeta_+", FullSimplify[(zetaReq /. PiTr -> piFail) - zetaPlus, Assumptions -> $Assumptions]];

rQuad = FullSimplify[zetaReq - zetaPhys, Assumptions -> $Assumptions];
Print["R_quad = ", fmt[rQuad]];

expectZero[
  "R_quad(Pi_suff, zeta_phys=zeta_-)",
  FullSimplify[rQuad /. {PiTr -> piSuff, zetaPhys -> zetaMinus}, Assumptions -> $Assumptions]
];
expectZero[
  "R_quad(Pi_fail, zeta_phys=zeta_+)",
  FullSimplify[rQuad /. {PiTr -> piFail, zetaPhys -> zetaPlus}, Assumptions -> $Assumptions]
];

(* F3 (v2): directional content of zeta_req. The inverse-map theorem
   (notes section 4) relies on zeta_req being strictly increasing in Pi_tr
   on the physical branch (where PiTr, cMix, epsBlk are positive). Verify
   by factoring d zeta_req / d Pi_tr into a sign-controlled numerator /
   denominator pair under those assumptions. *)
dZetaReqDPi = FullSimplify[D[zetaReq, PiTr], Assumptions -> $Assumptions];
Print["dzeta_req/dPi_tr = ", fmt[dZetaReqDPi]];
{numD, denD} = {Numerator[Together[dZetaReqDPi]], Denominator[Together[dZetaReqDPi]]};
expectZero[
  "numerator(d zeta_req/d Pi_tr) - C_mix*(1 - eps_blk)",
  FullSimplify[numD - cMix*(1 - epsBlk), Assumptions -> $Assumptions]
];
expectZero[
  "denominator(d zeta_req/d Pi_tr) - (C_mix - eps_blk*(2*C_mix - Pi_tr))^2",
  FullSimplify[denD - (cMix - epsBlk*(2*cMix - PiTr))^2, Assumptions -> $Assumptions]
];

(* F1 (v2 paper-alignment Q3 direction (a)): closed-form pin for zeta_phys.
   Paper eq. app-stage082-zeta-phys:
     zeta_phys(Pe, eta; kappa) = Omega_Pe(Pe)^2 * (kappa + Pi^2/4) / (kappa + y(eta)^2)
   with Omega_Pe(Pe) = Pi*Pe*(2*Pe*Exp[Pe] + Pi) / ((4*Pe^2 + Pi^2)*(Exp[Pe] - 1))
   and y(eta) the smallest positive root of y*Tan[y] = eta.
   Verify by reproducing the Pe->Infinity limit at Family-1 (eta, kappa) = (37, 12321/5),
   matching stage 084's zeta_max^(F1) constant. *)
Module[{peSym, kappaSym, ySym, OmegaPe, OmegaPeLimit,
         yF1, kappaF1, zetaPhysF1Limit, zetaMaxF1Reference, diffToReference},
  ClearAll[peSym, kappaSym, ySym];
  OmegaPe = Pi*peSym*(2*peSym*Exp[peSym] + Pi) / ((4*peSym^2 + Pi^2)*(Exp[peSym] - 1));
  OmegaPeLimit = Limit[OmegaPe, peSym -> Infinity];
  Print["Omega_Pe -> ", fmt[OmegaPeLimit], " as Pe -> oo"];
  expectZero["Omega_Pe -> Pi/2 as Pe -> oo", OmegaPeLimit - Pi/2];
  yF1 = ySym /. FindRoot[ySym*Tan[ySym] - 37 == 0, {ySym, 1.527}, WorkingPrecision -> 30];
  Print["y_F1 (root of y tan y = 37, smallest positive) = ", fmt[yF1]];
  kappaF1 = 12321/5;
  zetaPhysF1Limit = FullSimplify[(Pi^2/4) * (kappaF1 + Pi^2/4) / (kappaF1 + yF1^2)];
  Print["zeta_phys(Pe->oo, kappa_F1, y_F1) = ", fmt[N[zetaPhysF1Limit, 20]]];
  zetaMaxF1Reference = ToExpression["2.467529229455835`30"];
  diffToReference = Abs[N[zetaPhysF1Limit - zetaMaxF1Reference, 30]];
  Print["|zeta_phys(F1, Pe->oo) - zeta_max^(F1)| = ", fmt[diffToReference]];
  If[TrueQ[diffToReference < 10^-10],
    pass["Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10"],
    fail["Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10", diffToReference]];
];

(* Carry-forward constants (lambdaEll = 37 from stages 056/073; upsilonW = 100 thetaW
   from stage 075 with alphaR = 10). After v2 paper-alignment Q2 direction (a):
   paper/stages/stage_075.tex Inputs line and notes/.../stage075...md were updated to
   state Upsilon_w = 100 Theta_w, fully consistent with the script's value here. *)
lambdaEll = 37;
xiF1FromUpsilon = FullSimplify[upsilonW*lambdaEll^2, Assumptions -> $Assumptions];
xiF1FromTheta = FullSimplify[100*thetaW*lambdaEll^2, Assumptions -> $Assumptions];

Print["Xi_F1 from Upsilon_w = ", fmt[xiF1FromUpsilon]];
Print["Xi_F1 from Theta_w = ", fmt[xiF1FromTheta]];
(* Note: the two equalities below are arithmetic on the hand-supplied integers
   37, 100, 1369, 136900. They are not independent verifications of the
   Family-1 strength identity — the upstream stage that derives lambdaEll = 37
   and the convention upsilonW = 100 * thetaW is responsible for those facts.
   Here we only display the resulting Xi_F1 forms for downstream readers. *)
Print["Xi_F1(Theta_w) - 136900 Theta_w = ", fmt[FullSimplify[xiF1FromTheta - 136900*thetaW]], "  (display only)"];
Print["Xi_F1(Upsilon_w=100 Theta_w) - Xi_F1(Theta_w) = ",
      fmt[FullSimplify[(xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta, Assumptions -> thetaW > 0]],
      "  (display only)"];

Print[""];
Print["Theorem ledger:"];
Print["  Pi_tr <= C_mix Q(zeta_-,eps_blk)  -> guaranteed success if zeta_phys >= zeta_-"];
Print["  Pi_tr >= C_mix Q(zeta_+,eps_blk)  -> guaranteed failure if zeta_phys <= zeta_+"];
Print["  Xi_F1 = 1369 Upsilon_w = 136900 Theta_w"];

Exit[0];
