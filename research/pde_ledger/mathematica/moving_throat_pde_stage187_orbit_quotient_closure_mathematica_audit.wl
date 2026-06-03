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

banner["STAGE 187 — EXACT ORBIT-QUOTIENT CLOSURE"];

Clear[dl, dc, dg, du, deta, dw, dm, dt, chiStar, deltaStar, eStar, fStar];

$Assumptions =
  Element[{dl, dc, dg, du, deta, dw, dm, dt, chiStar, deltaStar, eStar, fStar}, Reals] &&
  chiStar > 0 && deltaStar > 0 && eStar > 0 && fStar > 0;

rowTr = (1 + deltaStar)*(dg + dc - du) + (1 + chiStar)*(dt - du);
rowNt = 2*(1 + eStar)*dl + 2*eStar*dg + (fStar - eStar)*du - deta - (2 + eStar)*dw + dm - fStar*dt;
rowEta = 2*dc - du - deta;

(* Positive primitive ratios (xtilde/x); declare positivity for Log expansion. *)
$Assumptions =
  $Assumptions &&
  Element[{rL, rC, rG, rU, rEta, rW, rM, rT}, Reals] &&
  rL > 0 && rC > 0 && rG > 0 && rU > 0 &&
  rEta > 0 && rW > 0 && rM > 0 && rT > 0;
logSubs = {
  dl -> Log[rL], dc -> Log[rC], dg -> Log[rG], du -> Log[rU],
  deta -> Log[rEta], dw -> Log[rW], dm -> Log[rM], dt -> Log[rT]
};
ctrRatio = (rG*rC/rU)^(1 + deltaStar) * (rT/rU)^(1 + chiStar);
cntRatio = (rL^2*rM/(rEta*rW^2)) * (rG^2*rL^2/(rU*rW))^eStar * (rT/rU)^(-fStar);
epsEtaRatio = rC^2/(rU*rEta);

Print["Exact finite log-ratio equations:"];
Print["row_tr  = ", fmt[rowTr]];
Print["row_nt  = ", fmt[rowNt]];
Print["row_eta = ", fmt[rowEta]];

expectZero["log C_tr ratio - row_tr", PowerExpand[Log[ctrRatio]] - (rowTr /. logSubs)];
expectZero["log C_nt ratio - row_nt", PowerExpand[Log[cntRatio]] - (rowNt /. logSubs)];
expectZero["log epsilon_eta ratio - row_eta", PowerExpand[Log[epsEtaRatio]] - (rowEta /. logSubs)];

m = {
  {0, 1 + deltaStar, 1 + deltaStar, -(2 + chiStar + deltaStar), 0, 0, 0, 1 + chiStar},
  {2*(1 + eStar), 0, 2*eStar, fStar - eStar, -1, -(2 + eStar), 1, -fStar},
  {0, 2, 0, -1, -1, 0, 0, 0}
};
dx = {dl, dc, dg, du, deta, dw, dm, dt};
mx = Expand[m.dx];

expectZero["matrix row 1 - exact row_tr", mx[[1]] - rowTr];
expectZero["matrix row 2 - exact row_nt", mx[[2]] - rowNt];
expectZero["matrix row 3 - exact row_eta", mx[[3]] - rowEta];

minor = {
  {0, 0, 1 + chiStar},
  {-1, 1, -fStar},
  {-1, 0, 0}
};
Print["det selected minor (Delta_eta, Delta_mu, Delta_T) = ", fmt[Det[minor]]];
expectZero["selected minor determinant", Det[minor] - (1 + chiStar)];

sol = First[Solve[{rowTr == 0, rowNt == 0, rowEta == 0}, {deta, dt, dm}, Reals]];

Print[""];
Print["Exact finite fibre solution:"];
Print["Delta_eta = ", fmt[FullSimplify[deta /. sol, Assumptions -> $Assumptions]]];
Print["Delta_T = ", fmt[FullSimplify[dt /. sol, Assumptions -> $Assumptions]]];
Print["Delta_mu = ", fmt[FullSimplify[dm /. sol, Assumptions -> $Assumptions]]];

detaExpected = 2*dc - du;
dtExpected = du - (1 + deltaStar)*(dg + dc - du)/(1 + chiStar);
dmExpected = FullSimplify[
  2*dc - du + 2*dw - 2*dl
  - eStar*(2*dg + 2*dl - du - dw)
  - fStar*((1 + deltaStar)/(1 + chiStar))*(dg + dc - du),
  Assumptions -> $Assumptions
];

expectZero["Delta_eta finite law", (deta /. sol) - detaExpected];
expectZero["Delta_T finite law", (dt /. sol) - dtExpected];
expectZero["Delta_mu finite law", (dm /. sol) - dmExpected];

expectZero["row_tr after solve", rowTr /. sol];
expectZero["row_nt after solve", rowNt /. sol];
expectZero["row_eta after solve", rowEta /. sol];

banner["Finite orbit interpretation"];
Print["The three monomial equalities reduce exactly to the same rank-3 matrix condition"];
Print["M_* Delta x = 0, but now Delta x is a finite log-ratio vector rather than an"];
Print["infinitesimal drift. Therefore the Stage 186 similarity orbit integrates exactly:"];
Print["its fibres are the full level sets of (C_tr, C_nt, epsilon_eta)."];

banner["Carry-forward formulas"];
Print["  Delta_eta = 2 Delta_c - Delta_U"];
Print["  Delta_T   = Delta_U - ((1+deltaU_*)/(1+chi_*))(Delta_gamma + Delta_c - Delta_U)"];
Print["  Delta_mu  = 2 Delta_c - Delta_U + 2 Delta_W - 2 Delta_lambda"];
Print["              - E_*(2 Delta_gamma + 2 Delta_lambda - Delta_U - Delta_W)"];
Print["              - F_*((1+deltaU_*)/(1+chi_*))(Delta_gamma + Delta_c - Delta_U)"];
Print[""];
Print["Conclusion:"];
Print["  Equal values of (C_tr, C_nt, epsilon_eta) are exactly equivalent to lying on"];
Print["  the same finite similarity orbit G_*. The weak-axisymmetric defect therefore"];
Print["  lives in the exact three-dimensional orbit quotient."];

Exit[0];
