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

banner["STAGE 087 — EXACT OUTGOING l=2 DTN FINGERPRINT"];

Clear[z];
$Assumptions = Element[z, Reals];

h2 = FullSimplify[SphericalHankelH1[2, z], Assumptions -> $Assumptions];
lam = FullSimplify[z*D[h2, z]/h2, Assumptions -> $Assumptions];
yOut = FullSimplify[-3/lam, Assumptions -> $Assumptions];

lamSeries = Expand[Normal[Series[lam, {z, 0, 7}]]];
ySeries = Expand[Normal[Series[yOut, {z, 0, 7}]]];

Print["h_2^(1)(z) = ", fmt[h2]];
Print["Lambda_2^out(z) = ", fmt[lamSeries]];
Print["Y_2^out(z) = ", fmt[ySeries]];

expectZero["static DtN slot", Coefficient[lamSeries, z, 0] + 3];
expectZero["Y z^2 coefficient", Coefficient[ySeries, z, 2] - 1/9];
expectZero["Y z^4 coefficient", Coefficient[ySeries, z, 4] - 4/81];
expectZero["Y imag z^5 coefficient", Coefficient[ySeries, z, 5]/I - 1/27];
expectZero["Y z^6 coefficient", Coefficient[ySeries, z, 6] + 11/729];
expectZero["Y imag z^7 coefficient", Coefficient[ySeries, z, 7]/I + 1/243];

Clear[aThroat, omega, cSound];
$Assumptions =
  Element[{aThroat, omega, cSound}, Reals] &&
  aThroat > 0 && cSound > 0;

zSub = aThroat*omega/cSound;
yOmega = Expand[Normal[Series[yOut /. z -> zSub, {omega, 0, 5}]]];
Print["Y_2^out(omega) through O(omega^5) = ", fmt[yOmega]];
expectZero["omega^2 coefficient", Coefficient[yOmega, omega, 2] - aThroat^2/(9*cSound^2)];
expectZero["omega^4 coefficient", Coefficient[yOmega, omega, 4] - 4*aThroat^4/(81*cSound^4)];
expectZero["imag omega^5 coefficient", Coefficient[yOmega, omega, 5]/I - aThroat^5/(27*cSound^5)];

Print[""];
Print["Stage 087 Mathematica audit passed."];

Exit[0];
