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

banner["STAGE 083 — OUTGOING NORMALIZATION FACTORIZATION"];

Clear[gNewton, cLight, omegaQ, k0, mHat0, chiQ, omega, nQSym];
$Assumptions =
  Element[{gNewton, cLight, omegaQ, k0, mHat0, chiQ, omega, nQSym}, Reals] &&
  gNewton > 0 && cLight > 0 && omegaQ > 0 && k0 > 0 && mHat0 > 0 && chiQ > 0 && nQSym > 0;

sigmaCan = FullSimplify[(9/8)/omegaQ^5, Assumptions -> $Assumptions];
yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5);
ySeries = Expand[Normal[Series[yRet, {omega, 0, 5}]]];

k2 = FullSimplify[k0*Coefficient[ySeries, omega, 2], Assumptions -> $Assumptions];
k4 = FullSimplify[k0*Coefficient[ySeries, omega, 4], Assumptions -> $Assumptions];
gamma5 = FullSimplify[(Coefficient[ySeries, omega, 5]/I)*k0, Assumptions -> $Assumptions];

k0Target = FullSimplify[64*gNewton*omegaQ^5/(45*cLight^5), Assumptions -> $Assumptions];
k2Target = FullSimplify[k0Target/(4*omegaQ^2), Assumptions -> $Assumptions];
k4Target = FullSimplify[k0Target/(4*omegaQ^4), Assumptions -> $Assumptions];
gamma5Target = FullSimplify[2*gNewton/(5*cLight^5), Assumptions -> $Assumptions];
nQ = FullSimplify[k0/k0Target, Assumptions -> $Assumptions];

Print["Yhat_Q^ret series = ", fmt[ySeries]];
Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];
Print["Gamma5 = ", fmt[gamma5]];
Print["NQ = ", fmt[nQ]];

expectZero["K2/K2_target - NQ", k2/k2Target - nQ];
expectZero["K4/K4_target - NQ", k4/k4Target - nQ];
expectZero["Gamma5/Gamma5_target - chiQ*NQ", gamma5/gamma5Target - chiQ*nQ];
expectZero[
  "mhat0^2*Gamma5/Gamma5_target - mhat0^2*chiQ*NQ",
  mHat0^2*gamma5/gamma5Target - mHat0^2*chiQ*nQ
];

nQSol = First[Solve[mHat0^2*chiQ*nQSym == 1, nQSym]];
expectZero["NQ - 1/(mhat0^2*chiQ) after odd normalization", (nQSym /. nQSol) - 1/(mHat0^2*chiQ)];

Print[""];
Print["Stage 083 Mathematica audit passed."];

Exit[0];
