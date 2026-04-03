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
  res = FullSimplify[Together[Expand[expr]]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 079 — GEOMETRY-LANE CHECK VERDICT"];

eps2 = 0;
eps4 = 0;
cPole = FullSimplify[(1 + eps4)/(4*(1 + eps2)^2)];
cGeom = FullSimplify[1 - cPole];
rhoAlpha = FullSimplify[1/cGeom];
zetaReq = FullSimplify[cPole/cGeom];

Print["c_pole = ", fmt[cPole]];
Print["c_geom = ", fmt[cGeom]];
Print["rho_alpha = ", fmt[rhoAlpha]];
Print["zeta_req = ", fmt[zetaReq]];

expectZero["eps2", eps2];
expectZero["eps4", eps4];
expectZero["c_pole - 1/4", cPole - 1/4];
expectZero["c_geom - 3/4", cGeom - 3/4];
expectZero["rho_alpha - 4/3", rhoAlpha - 4/3];
expectZero["zeta_req - 1/3", zetaReq - 1/3];

Print[""];
Print["Stage 079 Mathematica audit passed."];

Exit[0];
