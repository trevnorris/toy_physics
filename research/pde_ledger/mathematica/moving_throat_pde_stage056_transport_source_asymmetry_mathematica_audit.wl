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

banner["STAGE 039 — TRANSPORT ORIGIN OF THE SOURCE-SHAPE ASYMMETRY"];

Clear[s, ell, dSigma, vSigma, pe];
$Assumptions =
  Element[{s, ell, dSigma, vSigma, pe}, Reals] && ell > 0 && dSigma > 0 && vSigma > 0 && pe > 0;

sigma = Exp[vSigma s/dSigma];
jFlux = FullSimplify[-dSigma D[sigma, s] + vSigma sigma, Assumptions -> $Assumptions];
expectZero["zero-flux transport residual", jFlux];

sigmaPe = FullSimplify[pe Exp[pe s/ell]/(Exp[pe] - 1), Assumptions -> $Assumptions];
norm = FullSimplify[Integrate[sigmaPe, {s, 0, ell}], Assumptions -> $Assumptions];
chi0 = Sqrt[2/ell] Sin[Pi s/(2 ell)];
iW = FullSimplify[Integrate[chi0, {s, 0, ell}], Assumptions -> $Assumptions];
iPe = FullSimplify[Integrate[sigmaPe chi0, {s, 0, ell}], Assumptions -> $Assumptions];
omegaPe = FullSimplify[iPe/iW, Assumptions -> $Assumptions];
omegaExpected = FullSimplify[
  Pi pe (2 pe Exp[pe] + Pi)/((4 pe^2 + Pi^2) (Exp[pe] - 1)),
  Assumptions -> $Assumptions
];

Print["I_W = ", fmt[iW]];
Print["I_Pe = ", fmt[iPe]];
Print["Omega_Pe = ", fmt[omegaPe]];
expectZero["normalization int sigma_Pe ds - ell", norm - ell];
expectZero["Omega_Pe - expected formula", omegaPe - omegaExpected];

omega0 = FullSimplify[Limit[omegaPe, pe -> 0], Assumptions -> ell > 0];
omegaInf = FullSimplify[Limit[omegaPe, pe -> Infinity], Assumptions -> ell > 0];
Print["Omega_Pe(0) = ", fmt[omega0]];
Print["lim Pe->+infty Omega_Pe = ", fmt[omegaInf]];
expectZero["twin baseline", omega0 - 1];
expectZero["upper finite-throat overlap limit", omegaInf - Pi/2];

pPe = FullSimplify[sigmaPe/ell, Assumptions -> $Assumptions];
eChi = FullSimplify[Integrate[pPe chi0, {s, 0, ell}], Assumptions -> $Assumptions];
eS = FullSimplify[Integrate[pPe s, {s, 0, ell}], Assumptions -> $Assumptions];
eChiS = FullSimplify[Integrate[pPe chi0 s, {s, 0, ell}], Assumptions -> $Assumptions];
cov = FullSimplify[eChiS - eChi eS, Assumptions -> $Assumptions];
expectZero["dOmega/dPe - Cov/I_W", D[omegaPe, pe] - cov/iW];

omegaSmall = FullSimplify[Normal[Series[omegaPe, {pe, 0, 2}]], Assumptions -> ell > 0];
smallExpected = Expand[1 + ((4 - Pi)/(2 Pi)) pe];
omegaLarge = FullSimplify[Normal[Series[omegaPe, {pe, Infinity, 2}]], Assumptions -> ell > 0];
largeCoeff = FullSimplify[Limit[pe^2 (omegaPe - Pi/2), pe -> Infinity], Assumptions -> ell > 0];

Print["Omega_Pe small-Pe series = ", fmt[Expand[omegaSmall]]];
Print["Omega_Pe large-Pe series = ", fmt[Expand[omegaLarge]]];
Print["large-Pe second-order coefficient = ", fmt[largeCoeff]];
expectZero[
  "small-Pe linear coefficient",
  Coefficient[Expand[omegaSmall], pe, 1] - Coefficient[smallExpected, pe, 1]
];
expectZero[
  "large-Pe second-order coefficient + Pi^3/8",
  largeCoeff + Pi^3/8
];

Print[""];
Print["Stage 056 Mathematica audit passed."];

Exit[0];
