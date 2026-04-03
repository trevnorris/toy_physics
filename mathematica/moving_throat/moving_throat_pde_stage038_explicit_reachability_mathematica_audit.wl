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

banner["STAGE 038 — EXPLICIT LOWEST-LANE REACHABILITY"];

Clear[alpha, x, y, zetaReq];
$Assumptions =
  Element[{alpha, x, y, zetaReq}, Reals] && alpha > 0 && 0 < x < 4 && y > 0 && zetaReq > 0;

omegaExp = FullSimplify[
  Pi alpha (2 alpha Exp[alpha] + Pi)/((4 alpha^2 + Pi^2) (Exp[alpha] - 1)),
  Assumptions -> $Assumptions
];
aK = FullSimplify[1/(1 - x/4 + x y^2/Pi^2), Assumptions -> $Assumptions];
zetaFamily = FullSimplify[omegaExp^2 aK, Assumptions -> $Assumptions];

Print["Omega_exp(alpha) = ", fmt[omegaExp]];
Print["A_K(y,x) = ", fmt[aK]];
Print["zeta_family(alpha,y;x) = ", fmt[zetaFamily]];

zetaTwin = FullSimplify[(Limit[zetaFamily, alpha -> 0]) /. y -> Pi/2, Assumptions -> 0 < x < 4];
zetaMax = FullSimplify[Limit[Limit[zetaFamily, alpha -> Infinity], y -> 0, Direction -> "FromAbove"], Assumptions -> 0 < x < 4];
xFloor = FullSimplify[4 - Pi^2/zetaReq, Assumptions -> zetaReq > 0];

Print["symmetric twin point = ", fmt[zetaTwin]];
Print["closure maximum = ", fmt[zetaMax]];
Print["x floor from zeta_max = zeta_req = ", fmt[xFloor]];
expectZero["twin value", zetaTwin - 1];
expectZero["closure maximum", zetaMax - Pi^2/(4 - x)];
expectZero["x floor = 4 - Pi^2/zeta_req", xFloor - (4 - Pi^2/zetaReq)];
expectZero["zeta_max(x_floor) - zeta_req", (zetaMax /. x -> xFloor) - zetaReq];
expectZero["KX/KW equivalence", ((1 - x/4) /. x -> xFloor) - Pi^2/(4 zetaReq)];

Print[""];
Print["Stage 038 Mathematica audit passed."];

Exit[0];
