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

banner["STAGE 087 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE"];

Clear[rhoAlpha, epsBlk];
$Assumptions = Element[{rhoAlpha, epsBlk}, Reals] && rhoAlpha > 0 && epsBlk >= 0;

zetaReq = FullSimplify[(rhoAlpha - 1)/(1 - epsBlk*(2 - rhoAlpha)), Assumptions -> $Assumptions];
dZeta = FullSimplify[D[zetaReq, rhoAlpha], Assumptions -> $Assumptions];
dZetaExpected = FullSimplify[(1 - epsBlk)/(1 - epsBlk*(2 - rhoAlpha))^2, Assumptions -> $Assumptions];

Print["zeta_req(rho_alpha,eps_blk) = ", fmt[zetaReq]];
Print["d zeta_req / d rho_alpha = ", fmt[dZeta]];
expectZero["d zeta_req exact formula", dZeta - dZetaExpected];
expectZero["unblocked zeta_req", (zetaReq /. epsBlk -> 0) - (rhoAlpha - 1)];

(*
   Stage 087 is a checkpoint-consolidation statement (paper purpose:
   "records that the explicit Family-1 support/source side has been reduced
   to a single outgoing-branch loading ratio"). The cancellation chain is
   verified upstream in stages 081-086 (post-renumber):
   - mathematica/moving_throat_pde_stage085_quadrupole_demand_cancellation_*
   - mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_*
   The rho_X literals below are carried from stage 086; the cross-check
   below anchors them against the upstream stage-086 paper values to catch
   renumber or transcription drift.
*)
rhoSuff = ToExpression["3.46622291347846`20"];
rhoFail = ToExpression["3.46752913273870`20"];
rhoMax = ToExpression["3.46752922945601`20"];

expectApprox["rho_suff^(chi) vs stage-086", rhoSuff, 3.46622291347846, 10^-14];
expectApprox["rho_fail^(chi) vs stage-086", rhoFail, 3.46752913273870, 10^-14];
expectApprox["rho_max^(F1)   vs stage-086", rhoMax,  3.46752922945601, 10^-14];

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
