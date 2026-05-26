(* Unit 018 Mathematica audit.
   Mathematica counterpart of the SymPy audit for the parent throat action bundle master.
   Claims:
   M1 Gaussian wall inertia and stiffness integrals.
*)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Print["STAGE 018 PARENT THROAT ACTION BUNDLE MASTER MATHEMATICA AUDIT"];

Clear[w];
$Assumptions = Element[w, Reals];
beta = Exp[-w^2/2];
massIntegral = Integrate[beta^2, {w, -Infinity, Infinity}];
stiffnessIntegral = Integrate[D[beta, w]^2 + beta^2, {w, -Infinity, Infinity}];
Print["M1 Gaussian inertia integral residual = ",
  fmt[FullSimplify[massIntegral - Sqrt[Pi]]]];
If[FullSimplify[massIntegral - Sqrt[Pi]] =!= 0,
  (Print["FAIL: M1 Gaussian inertia integral"]; Exit[1])];

Print["M1 Gaussian stiffness integral residual = ",
  fmt[FullSimplify[stiffnessIntegral - 3*Sqrt[Pi]/2]]];
If[FullSimplify[stiffnessIntegral - 3*Sqrt[Pi]/2] =!= 0,
  (Print["FAIL: M1 Gaussian stiffness integral"]; Exit[1])];

Print["STAGE 018 MATHEMATICA AUDIT PASS"];
Exit[0];
