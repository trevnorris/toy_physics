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

banner["STAGE 075 — DYNAMIC-GEOMETRY OBSTRUCTION"];

Clear[omega, kg0, kg2, kg4, kPole, omegaQ, eps2, eps4];
$Assumptions =
  Element[{omega, kg0, kg2, kg4, kPole, omegaQ, eps2, eps4}, Reals] &&
  kPole > 0 && omegaQ > 0;

kCons = FullSimplify[kg0 + kg2*omega^2 + kg4*omega^4 + kPole/(1 - omega^2/omegaQ^2), Assumptions -> $Assumptions];
series = Expand[Normal[Series[kCons, {omega, 0, 4}]]];
k0 = FullSimplify[Coefficient[series, omega, 0], Assumptions -> $Assumptions];
k2 = FullSimplify[Coefficient[series, omega, 2], Assumptions -> $Assumptions];
k4 = FullSimplify[Coefficient[series, omega, 4], Assumptions -> $Assumptions];

Print["K_Q^cons(omega) = ", fmt[kCons]];
Print["Series = ", fmt[series]];
Print["K0 = ", fmt[k0]];
Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];

branch = FullSimplify[k0*k4 - 4*k2^2, Assumptions -> $Assumptions];
Print["Branch obstruction = ", fmt[branch]];

kg0Sol = FullSimplify[First[Solve[branch == 0, kg0]], Assumptions -> $Assumptions];
kg0Branch = FullSimplify[kg0 /. kg0Sol, Assumptions -> $Assumptions];
Print["K_g0 on branch = ", fmt[Factor[kg0Branch]]];
expectZero["static-geometry limit gives 3 K_pole", (kg0Branch /. {kg2 -> 0, kg4 -> 0}) - 3*kPole];

cPole = FullSimplify[kPole/(kg0Branch + kPole), Assumptions -> $Assumptions];
cPoleDimless = FullSimplify[cPole /. {kg2 -> eps2*kPole/omegaQ^2, kg4 -> eps4*kPole/omegaQ^4}, Assumptions -> $Assumptions];
cPoleExpected = FullSimplify[(1 + eps4)/(4*(1 + eps2)^2), Assumptions -> $Assumptions];

Print["c_pole = ", fmt[Factor[cPole]]];
Print["c_pole in (eps2,eps4) variables = ", fmt[Factor[cPoleDimless]]];
expectZero["c_pole - (1+eps4)/(4(1+eps2)^2)", cPoleDimless - cPoleExpected];

cGeomDimless = FullSimplify[1 - cPoleDimless, Assumptions -> $Assumptions];
Print["c_geom in (eps2,eps4) variables = ", fmt[Factor[cGeomDimless]]];

smallSeries = Expand[Normal[Series[Normal[Series[cPoleExpected, {eps2, 0, 1}]], {eps4, 0, 1}]]];
linearPart = Expand[(1/4)*(1 + eps4 - 2*eps2)];
Print["First-order expansion of c_pole = ", fmt[smallSeries]];
Print["Linear part = ", fmt[linearPart]];
Print["Dropped higher-order tail = ", fmt[Expand[smallSeries - linearPart]]];

Print[""];
Print["Stage 092 Mathematica audit passed."];

Exit[0];
