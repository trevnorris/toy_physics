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

banner["STAGE 137 — EXPLICIT CORE-TO-MOUTH GAIN MAP"];

Clear[lM, thetaSigma, kS, kQ, lam, gS, gQ, dSch, kappa0, gamma0, mCore, vCoup,
  deltaLambdaSchur, rhoCSchur, sigmaCSchur, rC, kappaC, gammaC, zVar, piVar,
  sqVar, piMap, mixedContribution, mQFromSchur, vacuityResidual];
$Assumptions =
  Element[{lM, thetaSigma, kS, kQ, lam, gS, gQ, dSch, kappa0, gamma0, kappaC, gammaC,
    zVar, piVar, sqVar}, Reals] &&
  lM > 0 && thetaSigma > 0 && kS > 0 && kQ > 0 && gS > 0 && gQ > 0 &&
  dSch > 0 && kappa0 > 0 && gamma0 > 0 && kappaC > 0 && gammaC > 0;

rhoC = FullSimplify[gS^2/kS, Assumptions -> $Assumptions];
sigmaC = FullSimplify[(kS*gQ - lam*gS)^2/(kS*(kS*kQ + lam^2)), Assumptions -> $Assumptions];

(* --- F1 (R3): independent matrix-Schur reconstruction of rhoC, sigmaC via Inverse. --- *)
(* Physical core stiffness matrix (notes stage097, owner script stage114); rhoC, sigmaC *)
(* are NOT inputs to it. Inverting it is the independent derivation primitive. *)
mCore = {{kS, lam}, {lam, -kQ*dSch}};
vCoup = {gS, gQ};
deltaLambdaSchur = Apart[(vCoup . Inverse[mCore] . vCoup), dSch];
rhoCSchur = FullSimplify[Limit[deltaLambdaSchur, dSch -> Infinity], Assumptions -> $Assumptions];
sigmaCSchur = FullSimplify[rhoCSchur - (deltaLambdaSchur /. dSch -> 1), Assumptions -> $Assumptions];
expectZero["rho_c equals M_core Schur residue (D -> Infinity)", rhoC - rhoCSchur];
expectZero["sigma_c equals M_core Schur residue (static D = 1)", sigmaC - sigmaCSchur];

mS = FullSimplify[lM*rhoC/thetaSigma, Assumptions -> $Assumptions];
mQ = FullSimplify[-lM*sigmaC/thetaSigma, Assumptions -> $Assumptions];

Print["rho_c = ", fmt[rhoC]];
Print["sigma_c = ", fmt[sigmaC]];
Print["M_s = ", fmt[mS]];
Print["M_q = ", fmt[mQ]];

(* Anchor M_s, M_q against direct paper-card closed forms (independent route). *)
mSPaper = lM * gS^2 / (kS * thetaSigma);
mQPaper = -lM * (kS*gQ - lam*gS)^2 / (kS * (kS*kQ + lam^2) * thetaSigma);
expectZero["M_s matches paper card", mS - mSPaper];
expectZero["M_q matches paper card", mQ - mQPaper];

(* F2 (R1): full-susceptibility anchor against the matrix-Schur source (NOT X - X). *)
(* Reduced envelope must equal Inverse[mCore] source on the bare denominator, using *)
(* the Stage 97/114 coefficient relations (notes stage097, stage114). *)
rC = lam^2/(kS*kQ);
kappaC = kappa0/(1 + rC);
gammaC = gamma0/(1 + rC);
dWbare = 1 - kappa0*zVar^2 - I*gamma0*zVar^5;
deltaLambdaMatrix = FullSimplify[deltaLambdaSchur /. dSch -> dWbare, Assumptions -> $Assumptions];
deltaLambdaReduced = rhoC - sigmaC/(1 - kappaC*zVar^2 - I*gammaC*zVar^5);
expectZero["reduced core susceptibility equals matrix-Schur source", deltaLambdaMatrix - deltaLambdaReduced];
(* Static specialization via Series (independent route from SymPy's Limit). *)
staticLimit = Normal[Series[deltaLambdaMatrix, {zVar, 0, 0}]];
expectZero["static core residue equals rho_c_schur - sigma_c_schur", staticLimit - (rhoCSchur - sigmaCSchur)];

(* F3 (R2): outlet consistency with a NONZERO mixed channel (paper Checks item 1). *)
(* Family-1 fixed-point law (notes stage137): Pi = mS + mQ*sqVar; the mixed term *)
(* mQ*sqVar must equal the matrix-Schur reconstruction -lM*sigmaCSchur*sqVar/thetaSigma. *)
piMap = mS + mQ*sqVar;
mixedContribution = FullSimplify[piMap - mS, Assumptions -> $Assumptions];
mQFromSchur = -lM*sigmaCSchur/thetaSigma;
expectZero["outlet mixed channel equals matrix-Schur M_q (sq nonzero)", mixedContribution - mQFromSchur*sqVar];
(* Non-vacuity: +mQ and -mQ must NOT both pass; the flipped-sign residual is nonzero. *)
vacuityResidual = FullSimplify[mixedContribution - (-mQFromSchur)*sqVar, Assumptions -> $Assumptions];
If[TrueQ[vacuityResidual === 0],
  fail["outlet check vacuous: +mQ and -mQ both pass", vacuityResidual],
  pass["outlet consistency non-vacuous (sign of M_q is exercised)"]
];

sigmaCRc = FullSimplify[(kS*gQ - lam*gS)^2/(kS^2*kQ*(1 + lam^2/(kS*kQ))), Assumptions -> $Assumptions];
Print["sigma_c (r_c form) = ", fmt[sigmaCRc]];
expectZero["sigma_c equivalence with r_c form", sigmaC - sigmaCRc];

Print[""];
Print["Final explicit gain map verified."];

Exit[0];
