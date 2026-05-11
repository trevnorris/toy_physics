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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 070 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE"];

Clear[rhoAlpha, epsBlk];
$Assumptions = Element[{rhoAlpha, epsBlk}, Reals] && rhoAlpha > 0 && epsBlk >= 0;

zetaReq = FullSimplify[(rhoAlpha - 1)/(1 - epsBlk*(2 - rhoAlpha)), Assumptions -> $Assumptions];
dZeta = FullSimplify[D[zetaReq, rhoAlpha], Assumptions -> $Assumptions];
dZetaExpected = FullSimplify[(1 - epsBlk)/(1 - epsBlk*(2 - rhoAlpha))^2, Assumptions -> $Assumptions];

Print["zeta_req(rho_alpha,eps_blk) = ", fmt[zetaReq]];
Print["d zeta_req / d rho_alpha = ", fmt[dZeta]];
expectZero["d zeta_req exact formula", dZeta - dZetaExpected];
expectZero["unblocked zeta_req", (zetaReq /. epsBlk -> 0) - (rhoAlpha - 1)];

rhoSuff = ToExpression["3.46622291347846`20"];
rhoFail = ToExpression["3.46752913273870`20"];
rhoMax = ToExpression["3.46752922945601`20"];

zetaSuff = N[zetaReq /. {rhoAlpha -> rhoSuff, epsBlk -> 0}, 30];
zetaFail = N[zetaReq /. {rhoAlpha -> rhoFail, epsBlk -> 0}, 30];
zetaMax = N[zetaReq /. {rhoAlpha -> rhoMax, epsBlk -> 0}, 30];

Print["zeta at success ratio = ", fmt[zetaSuff]];
Print["zeta at failure ratio = ", fmt[zetaFail]];
Print["zeta at max ratio = ", fmt[zetaMax]];

expectApprox["zeta_suff numeric check", zetaSuff, 2.466222913478460121439184, 10^-14];
expectApprox["zeta_fail numeric check", zetaFail, 2.467529132738699892968270, 10^-14];
expectApprox["zeta_max numeric check", zetaMax, 2.467529229456010053667114, 10^-14];

Print[""];
Print["Stage 087 Mathematica audit passed."];

Exit[0];
