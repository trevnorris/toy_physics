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

banner["STAGE 088 — EXACT FIXING OF chi_Q"];

Clear[aThroat, cSound, omega, chiQ];
$Assumptions =
  Element[{aThroat, cSound, omega, chiQ}, Reals] &&
  aThroat > 0 && cSound > 0 && chiQ > 0;

omegaQ = FullSimplify[3*cSound/(2*aThroat), Assumptions -> $Assumptions];
sigmaCan = FullSimplify[(9/8)/omegaQ^5, Assumptions -> $Assumptions];

Print["Omega_Q = ", fmt[omegaQ]];
Print["sigma_Q^can = ", fmt[sigmaCan]];
expectZero["sigma_Q^can - 4 a^5/(27 c_s^5)", sigmaCan - 4*aThroat^5/(27*cSound^5)];

yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5);
ySeries = Expand[Normal[Series[yRet, {omega, 0, 5}]]];

Print["Y_Q^ret(omega) = ", fmt[ySeries]];
expectZero["omega^2 coefficient", Coefficient[ySeries, omega, 2] - aThroat^2/(9*cSound^2)];
expectZero["omega^4 coefficient", Coefficient[ySeries, omega, 4] - 4*aThroat^4/(81*cSound^4)];
expectZero["imag omega^5 coefficient", Coefficient[ySeries, omega, 5]/I - chiQ*aThroat^5/(27*cSound^5)];

chiSol = First[Solve[(Coefficient[ySeries, omega, 5]/I) == aThroat^5/(27*cSound^5), chiQ]];
Print["chi_Q from exact outgoing match = ", fmt[chiQ /. chiSol]];
expectZero["chi_Q - 1", (chiQ /. chiSol) - 1];

Clear[z, xiQ];
$Assumptions = Element[{z, xiQ}, Reals];

lamDef = -3 + z^2/3 + z^4/9 + I*xiQ*z^5/9;
yDef = Expand[Normal[Series[-3/lamDef, {z, 0, 5}]]];

Print["Deformed DtN normalized branch = ", fmt[yDef]];
expectZero["deformed z^2 coefficient", Coefficient[yDef, z, 2] - 1/9];
expectZero["deformed z^4 coefficient", Coefficient[yDef, z, 4] - 4/81];
expectZero["deformed imag z^5 coefficient", Coefficient[yDef, z, 5]/I - xiQ/27];

Print[""];
Print["Stage 088 Mathematica audit passed."];

Exit[0];
