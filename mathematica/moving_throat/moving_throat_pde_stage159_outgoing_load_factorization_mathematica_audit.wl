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

banner["STAGE 159 — OUTGOING LOAD-FACTOR FACTORIZATION"];

Clear[k, ou2, ow2, r, gu, gw];
$Assumptions = Element[{k, ou2, ow2, r, gu, gw}, Reals] &&
  k > 0 && ou2 > 0 && ow2 > 0 && r > 0 && gu > 0 && gw > 0;

lambda = (ou2*gw + r*gu)/(ou2*ow2 - r^2);
mCal = gw/(Sqrt[k]*ow2);
iCal = r*gu/(ou2*gw);
hCal = r^2/(ou2*ow2);

expectZero[
  "exact factorization of Lambda^2/K",
  lambda^2/k - mCal^2*(1 + iCal)^2/(1 - hCal)^2
];

Print[""];
Print["Carry-forward exact identity:"];
Print["  Lambda^2/K = M^2 (1+I)^2 / (1-H)^2"];

banner["First-order logarithmic transport"];

Clear[eps, dK, dOU, dOW, dR, dGU, dGW];
$Assumptions = Element[{k, ou2, ow2, r, gu, gw, eps, dK, dOU, dOW, dR, dGU, dGW}, Reals] &&
  k > 0 && ou2 > 0 && ow2 > 0 && r > 0 && gu > 0 && gw > 0;

kP = k*Exp[eps*dK];
ou2P = ou2*Exp[eps*dOU];
ow2P = ow2*Exp[eps*dOW];
rP = r*Exp[eps*dR];
guP = gu*Exp[eps*dGU];
gwP = gw*Exp[eps*dGW];

lambdaP = (ou2P*gwP + rP*guP)/(ou2P*ow2P - rP^2);
sigmaExact = FullSimplify[
  D[Log[((lambdaP^2/kP)/(lambda^2/k))], eps] /. eps -> 0,
  Assumptions -> $Assumptions
];

Clear[dlnM, dlnI, dlnH];
dlnMExpr = dGW - dOW - dK/2;
dlnIExpr = dR + dGU - dOU - dGW;
dlnHExpr = 2*dR - dOU - dOW;

sigmaFactoredForm = 2*dlnM + 2*iCal*dlnI/(1 + iCal) + 2*hCal*dlnH/(1 - hCal);
sigmaFactored = sigmaFactoredForm /. {dlnM -> dlnMExpr, dlnI -> dlnIExpr, dlnH -> dlnHExpr};
expectZero["factored first-order defect formula", sigmaExact - sigmaFactored];

sigmaExpanded = FullSimplify[
  -dK +
  2*dGW/(1 + iCal) +
  2*iCal*dGU/(1 + iCal) +
  2*(iCal/(1 + iCal) + 2*hCal/(1 - hCal))*dR -
  2*(iCal/(1 + iCal) + hCal/(1 - hCal))*dOU -
  2*dOW/(1 - hCal),
  Assumptions -> $Assumptions
];
expectZero["expanded primitive-variable transport", sigmaExact - sigmaExpanded];

Print[""];
Print["First-order formulas:"];
Print["  d ln M = d ln G_W - d ln Omega_W^2 - 1/2 dK"];
Print["  d ln I = d ln R + d ln G_U - d ln Omega_U^2 - d ln G_W"];
Print["  d ln H = 2 d ln R - d ln Omega_U^2 - d ln Omega_W^2"];
Print["  Sigma^(N) = 2 d ln M + 2 I/(1+I) d ln I + 2 H/(1-H) d ln H"];

banner["Rigidity corollary"];
sigmaRigid = FullSimplify[sigmaFactoredForm /. {dlnI -> 0, dlnH -> 0}, Assumptions -> $Assumptions];
expectZero["rigidity reduction to 2 d ln M", sigmaRigid - 2*dlnM];

Print[""];
Print["Conclusion:"];
Print["  If the interference and hybridization ratios are rigid, then"];
Print["      Sigma^(N) = 2 d ln[ G_W / (Omega_W^2 sqrt(K)) ]."];
Print["  So the remaining linear grouped defect vanishes when the raw mixed leg"];
Print["  obeys the square-root wall-loading law G_W / Omega_W^2 ∝ sqrt(K)."];

Print[""];
Print["Stage 159 Mathematica audit passed."];

Exit[0];
