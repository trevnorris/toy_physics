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

banner["STAGE 065 — MASTER QUADRUPOLE RESIDUAL"];

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

(* Directional content of R_quad: verify the partial derivatives that
   underwrite the "guaranteed success / guaranteed failure" theorems. *)
dRDzetaPhys = FullSimplify[D[rQuad, zetaPhys], Assumptions -> $Assumptions];
expectZero["dR_quad/dzeta_phys + 1", dRDzetaPhys + 1];

dZetaReqDPi = FullSimplify[D[zetaReq, PiTr], Assumptions -> $Assumptions];
dRDPiAtZetaMinus = FullSimplify[D[rQuad /. zetaPhys -> zetaMinus, PiTr], Assumptions -> $Assumptions];
expectZero[
  "dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)",
  FullSimplify[dRDPiAtZetaMinus - dZetaReqDPi, Assumptions -> $Assumptions]
];

(* TODO(provenance): lambdaEll = 37 and the convention upsilonW = 100 * thetaW
   are carry-forward constants. Cite the upstream stage's script (likely an
   earlier moving_throat_pde_stage*_mathematica_audit.wl) that derives them.
   Until then, this stage treats them as inputs and only displays their
   consequences. *)
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
