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

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 066 — DIMENSIONLESS WALL FIGURE OF MERIT"];

Clear[v0, ell, a, len, j1, tx, peReq, delta0, deltaInf, hw, ifMom, v0Sq, kx, kappa];
$Assumptions =
  Element[{v0, ell, a, len, j1, tx, peReq, delta0, deltaInf, hw, ifMom, v0Sq, kx, kappa}, Reals] &&
  v0 > 0 && ell > 0 && a > 0 && len > 0 && j1 > 0 && tx > 0 && peReq > 0 &&
  delta0 > 0 && deltaInf > 0 && hw > 0 && ifMom > 0 && v0Sq > 0 && kx > 0 && kappa > 0;

wWall = FullSimplify[4*Pi*a^2*len^2*j1*v0^2/(tx*ell), Assumptions -> $Assumptions];
wFail = FullSimplify[peReq/deltaInf, Assumptions -> $Assumptions];
wSuff = FullSimplify[peReq/delta0, Assumptions -> $Assumptions];
gEqTw = FullSimplify[4*Pi*a^2*j1*v0^2/(kx*ell), Assumptions -> $Assumptions];

Print["W_wall = ", fmt[wWall]];
Print["W_fail = ", fmt[wFail]];
Print["W_suff = ", fmt[wSuff]];
expectZero["W_wall - kappa G_eq^(tw)", wWall - (gEqTw /. kx -> tx*kappa/len^2)*kappa];

v0FailSq = FullSimplify[tx*ell*peReq/(4*Pi*a^2*len^2*j1*deltaInf), Assumptions -> $Assumptions];
v0SuffSq = FullSimplify[tx*ell*peReq/(4*Pi*a^2*len^2*j1*delta0), Assumptions -> $Assumptions];

expectZero["W_wall(V0_fail)^2 - W_fail", (wWall /. v0^2 -> v0FailSq) - wFail];
expectZero["W_wall(V0_suff)^2 - W_suff", (wWall /. v0^2 -> v0SuffSq) - wSuff];

banner["MONOTONICITY"];

wV0Sq = FullSimplify[wWall /. v0^2 -> v0Sq, Assumptions -> $Assumptions];
dV0Sq = FullSimplify[D[wV0Sq, v0Sq], Assumptions -> $Assumptions];
dA = FullSimplify[D[wWall, a], Assumptions -> $Assumptions];
dLen = FullSimplify[D[wWall, len], Assumptions -> $Assumptions];
dEll = FullSimplify[D[wWall, ell], Assumptions -> $Assumptions];
dJ1 = FullSimplify[D[wWall, j1], Assumptions -> $Assumptions];
dTx = FullSimplify[D[wWall, tx], Assumptions -> $Assumptions];

Print["dW/d(V0^2) = ", fmt[dV0Sq]];
Print["dW/da = ", fmt[dA]];
Print["dW/dL = ", fmt[dLen]];
Print["dW/dell = ", fmt[dEll]];
Print["dW/dJ1 = ", fmt[dJ1]];
Print["dW/dT_X = ", fmt[dTx]];

expectTrue["dW/d(V0^2) > 0", dV0Sq > 0];
expectTrue["dW/da > 0", dA > 0];
expectTrue["dW/dL > 0", dLen > 0];
expectTrue["dW/dell < 0", dEll < 0];
expectTrue["dW/dJ1 > 0", dJ1 > 0];
expectTrue["dW/dT_X < 0", dTx < 0];

banner["CONSTANT-COMPRESSIBILITY FORM"];

wH = FullSimplify[4*Pi*a^2*len^2*ifMom*v0^2/(hw*tx*ell), Assumptions -> $Assumptions];
Print["W_H = ", fmt[wH]];

expectZero["J1 = I_f/H_w reduction", (wWall /. j1 -> ifMom/hw) - wH];
expectZero[
  "W_H(V0_fail)^2 - W_fail",
  (wH /. v0^2 -> hw*tx*ell*peReq/(4*Pi*a^2*len^2*ifMom*deltaInf)) - wFail
];
expectZero[
  "W_H(V0_suff)^2 - W_suff",
  (wH /. v0^2 -> hw*tx*ell*peReq/(4*Pi*a^2*len^2*ifMom*delta0)) - wSuff
];

Print[""];
Print["Stage 066 Mathematica audit passed."];

Exit[0];
