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

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 073 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION"];

rhoAlpha = 4/3;
zetaReq = rhoAlpha - 1;
rhoSuffChi = SetPrecision[3.46622291347846, 20];
zetaMaxF1 = SetPrecision[2.46752922945601, 20];
aF1 = SetPrecision[1.00005192880220, 20];

Print["rho_alpha = ", fmt[N[rhoAlpha, 25]]];
Print["zeta_req = ", fmt[N[zetaReq, 25]]];
Print["rho_suff^(chi) = ", fmt[rhoSuffChi]];
Print["zeta_max^(F1) = ", fmt[zetaMaxF1]];
Print["A_F1 = ", fmt[aF1]];

expectZero["zeta_req - (rho_alpha - 1)", zetaReq - (rhoAlpha - 1)];
expectTrue["rho_alpha lies inside the exact Family-1 success region", rhoAlpha < rhoSuffChi];
expectTrue["zeta_req lies below the hard Family-1 ceiling", zetaReq < zetaMaxF1];
expectTrue["zeta_req lies below the zero-bias Family-1 baseline", zetaReq < aF1];

Print[""];
Print["Stage 090 Mathematica audit passed."];

Exit[0];
