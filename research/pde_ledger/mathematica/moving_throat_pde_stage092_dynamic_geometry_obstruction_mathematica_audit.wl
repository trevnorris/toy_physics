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

banner["STAGE 092 — DYNAMIC-GEOMETRY OBSTRUCTION"];

Clear[omega, kg0, kg2, kg4, kPole, omegaQ, eps2, eps4];
$Assumptions =
  Element[{omega, kg0, kg2, kg4, kPole, omegaQ, eps2, eps4}, Reals] &&
  kPole > 0 && omegaQ > 0;

(* Work in dimensionless (eps2, eps4) variables from the outset.
   Notes Section 3 gives K_0 = 4 K_pole (1+eps_2)^2 / (1+eps_4) on the branch.
   We verify this is consistent with the original symbolic K_0, K_2, K_4
   coefficients of the conservative-carrier expansion, then read off c_pole. *)

(* Step A: define K_0, K_2, K_4 from the conservative-carrier expansion,
   then substitute eps2, eps4 directly without first solving for kg0. *)
k0Full = kg0 + kPole;
k2Full = kg2 + kPole/omegaQ^2;
k4Full = kg4 + kPole/omegaQ^4;

(* Dimensionless rewrite: kg2 = eps2 kPole/omegaQ^2, kg4 = eps4 kPole/omegaQ^4. *)
k2Eps = FullSimplify[k2Full /. {kg2 -> eps2*kPole/omegaQ^2}, Assumptions -> $Assumptions];
k4Eps = FullSimplify[k4Full /. {kg4 -> eps4*kPole/omegaQ^4}, Assumptions -> $Assumptions];
Print["K_2 in eps variables = ", fmt[k2Eps]];
Print["K_4 in eps variables = ", fmt[k4Eps]];

(* The branch identity K_0 K_4 = 4 K_2^2 gives K_0 = 4 K_2^2/K_4.
   Compute this directly without Solve[]. *)
k0FromBranch = FullSimplify[4*k2Eps^2/k4Eps, Assumptions -> $Assumptions];
Print["K_0 from branch (eps form) = ", fmt[k0FromBranch]];

(* Notes Section 3 prediction: K_0 = 4 K_pole (1+eps2)^2 / (1+eps4). *)
k0Predicted = FullSimplify[4*kPole*(1 + eps2)^2/(1 + eps4), Assumptions -> $Assumptions];
expectZero["K_0 closed form matches 4 K_pole (1+eps2)^2/(1+eps4)", k0FromBranch - k0Predicted];

(* c_pole = K_pole / K_0. *)
cPole = FullSimplify[kPole/k0FromBranch, Assumptions -> $Assumptions];
cPoleExpected = FullSimplify[(1 + eps4)/(4*(1 + eps2)^2), Assumptions -> $Assumptions];
Print["c_pole = ", fmt[Factor[cPole]]];
expectZero["c_pole - (1+eps4)/(4(1+eps2)^2)", cPole - cPoleExpected];

(* Static-geometry limit: eps2 = eps4 = 0 should give c_pole = 1/4. *)
expectZero["static-geometry limit c_pole = 1/4", (cPole /. {eps2 -> 0, eps4 -> 0}) - 1/4];

cGeom = FullSimplify[1 - cPole, Assumptions -> $Assumptions];
Print["c_geom in (eps2,eps4) variables = ", fmt[Factor[cGeom]]];
expectZero["static-geometry limit c_geom = 3/4", (cGeom /. {eps2 -> 0, eps4 -> 0}) - 3/4];

(* Small-contamination first-order expansion (notes Section 4). *)
smallSeries = Expand[Normal[Series[Normal[Series[cPoleExpected, {eps2, 0, 1}]], {eps4, 0, 1}]]];
linearPart = Expand[(1/4)*(1 + eps4 - 2*eps2)];
Print["First-order expansion of c_pole = ", fmt[smallSeries]];
Print["Linear part = ", fmt[linearPart]];
Print["Dropped higher-order tail = ", fmt[Expand[smallSeries - linearPart]]];

expectZero["first-order eps^0 coefficient", Coefficient[Coefficient[smallSeries, eps2, 0], eps4, 0] - 1/4];
expectZero["first-order eps2 coefficient", Coefficient[Coefficient[smallSeries, eps2, 1], eps4, 0] - (-1/2)];
expectZero["first-order eps4 coefficient", Coefficient[Coefficient[smallSeries, eps4, 1], eps2, 0] - 1/4];

Print[""];
Print["Stage 092 Mathematica audit passed."];

Exit[0];
