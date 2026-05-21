(*
  Stage 005 Mathematica audit.
  This script independently verifies the projected inhomogeneous Maxwell law
  and projected charge continuity for a Gaussian transverse projection. The
  concrete profiles are chosen separately from the SymPy witness, with explicit
  boundary-decay checks and sign-mutant guards.
*)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Clear[t, x, y, z, w];
projectionAssumptions = {t, x, y, z} \[Element] Reals;

Wg[w_] := Exp[-w^2]/Sqrt[Pi];
Wgp[w_] := D[Wg[w], w];

projW[expr_] := FullSimplify[
  Integrate[Wg[w] expr, {w, -Infinity, Infinity}, Assumptions -> projectionAssumptions],
  Assumptions -> projectionAssumptions
];

projWPrime[expr_] := FullSimplify[
  Integrate[Wgp[w] expr, {w, -Infinity, Infinity}, Assumptions -> projectionAssumptions],
  Assumptions -> projectionAssumptions
];

boundaryValue[expr_] := FullSimplify[
  Limit[expr, w -> Infinity] - Limit[expr, w -> -Infinity],
  Assumptions -> projectionAssumptions
];

Print["STAGE 005 -- projected Maxwell covariance audit"];

(* M1: brane derivative commutation. *)
Phi = (Sin[t] + x^2)*(w^2 + 3);
lhsM1 = projW[D[Phi, t]];
rhsM1 = D[projW[Phi], t];
Print["M1 residual = ", fmt[Simplify[lhsM1 - rhsM1]]];
If[Simplify[lhsM1 - rhsM1] =!= 0,
  Print["FAIL: M1 commutation with brane derivative"]; Exit[1],
  Print["PASS: M1 commutation with brane derivative"]
];
m1MutantResidual = lhsM1 + rhsM1;
Print["M1 mutant residual = ", fmt[Simplify[m1MutantResidual]]];
If[Simplify[m1MutantResidual] === 0,
  Print["MUTANT FAIL: M1 commutation sign guard"]; Exit[1],
  Print["PASS: M1 mutant guard"]
];

(* M2: projected integration by parts with decaying boundary. *)
Q = w*Exp[-w^2/4];
m2Boundary = boundaryValue[Wg[w]*Q];
Print["M2 boundary value = ", fmt[m2Boundary]];
If[Simplify[m2Boundary] =!= 0,
  Print["FAIL: M2 boundary decay"]; Exit[1],
  Print["PASS: M2 boundary decay"]
];
lhsM2 = projW[D[Q, w]] + projWPrime[Q];
rhsM2 = 0;
Print["M2 residual = ", fmt[Simplify[lhsM2 - rhsM2]]];
If[Simplify[lhsM2 - rhsM2] =!= 0,
  Print["FAIL: M2 projected IBP identity"]; Exit[1],
  Print["PASS: M2 projected IBP identity"]
];
m2MutantResidual = projW[D[Q, w]] - projWPrime[Q];
Print["M2 mutant residual = ", fmt[Simplify[m2MutantResidual]]];
If[Simplify[m2MutantResidual] === 0,
  Print["MUTANT FAIL: M2 IBP sign guard"]; Exit[1],
  Print["PASS: M2 mutant guard"]
];

(* M3: projected inhomogeneous Maxwell law with transverse leakage. *)
G0 = Cos[t]*(w^2 + 2);
Gx = x^2*w^2;
Gy = y*(w^4 + 1);
Gz = z*w^2;
Gw = w + w^3/3;
gammaTerm = (Sin[t] + x)*(w^2 + 1);

m3Boundary = boundaryValue[Wg[w]*Gw];
Print["M3 boundary value = ", fmt[m3Boundary]];
If[Simplify[m3Boundary] =!= 0,
  Print["FAIL: M3 boundary decay"]; Exit[1],
  Print["PASS: M3 boundary decay"]
];
lhsM3 = projW[D[G0, t] + D[Gx, x] + D[Gy, y] + D[Gz, z] + D[Gw, w] + gammaTerm];
rhsM3 = D[projW[G0], t] + D[projW[Gx], x] + D[projW[Gy], y] +
  D[projW[Gz], z] - projWPrime[Gw] + projW[gammaTerm];
Print["M3 residual = ", fmt[Simplify[lhsM3 - rhsM3]]];
If[Simplify[lhsM3 - rhsM3] =!= 0,
  Print["FAIL: M3 projected inhomogeneous Maxwell law"]; Exit[1],
  Print["PASS: M3 projected inhomogeneous Maxwell law"]
];
m3MutantResidual = lhsM3 - (
  D[projW[G0], t] + D[projW[Gx], x] + D[projW[Gy], y] +
  D[projW[Gz], z] + projWPrime[Gw] + projW[gammaTerm]
);
Print["M3 mutant residual = ", fmt[Simplify[m3MutantResidual]]];
If[Simplify[m3MutantResidual] === 0,
  Print["MUTANT FAIL: M3 leakage sign guard"]; Exit[1],
  Print["PASS: M3 mutant guard"]
];

(* M4: exact Gaussian leakage value. *)
leakageGw = w;
leakage = -projWPrime[leakageGw];
Print["M4 leakage cross-check = ", fmt[leakage]];
lhsM4 = leakage;
rhsM4 = 1;
Print["M4 residual = ", fmt[Simplify[lhsM4 - rhsM4]]];
If[Simplify[lhsM4 - rhsM4] =!= 0,
  Print["FAIL: M4 analytic transverse leakage value"]; Exit[1],
  Print["PASS: M4 analytic transverse leakage value"]
];
m4MutantResidual = projWPrime[leakageGw] - 1;
Print["M4 mutant residual = ", fmt[Simplify[m4MutantResidual]]];
If[Simplify[m4MutantResidual] === 0,
  Print["MUTANT FAIL: M4 leakage sign guard"]; Exit[1],
  Print["PASS: M4 mutant guard"]
];
If[Simplify[leakage] === 0,
  Print["MUTANT FAIL: M4 nonzero leakage guard"]; Exit[1],
  Print["PASS: M4 nonzero leakage guard"]
];

(* M5: projected charge continuity in open-system form. *)
J0 = Cos[t]*(w^2 + 2);
Jx = x^2*w^2;
Jy = y*(w^4 + 1);
Jz = z*w^2;
Jw = w + w^3/3;

m5Boundary = boundaryValue[Wg[w]*Jw];
Print["M5 boundary value = ", fmt[m5Boundary]];
If[Simplify[m5Boundary] =!= 0,
  Print["FAIL: M5 boundary decay"]; Exit[1],
  Print["PASS: M5 boundary decay"]
];
lhsM5 = projW[D[J0, t] + D[Jx, x] + D[Jy, y] + D[Jz, z] + D[Jw, w]];
rhsM5 = D[projW[J0], t] + D[projW[Jx], x] + D[projW[Jy], y] +
  D[projW[Jz], z] - projWPrime[Jw];
Print["M5 residual = ", fmt[Simplify[lhsM5 - rhsM5]]];
If[Simplify[lhsM5 - rhsM5] =!= 0,
  Print["FAIL: M5 projected charge continuity"]; Exit[1],
  Print["PASS: M5 projected charge continuity"]
];
m5MutantResidual = lhsM5 - (
  D[projW[J0], t] + D[projW[Jx], x] + D[projW[Jy], y] +
  D[projW[Jz], z] + projWPrime[Jw]
);
Print["M5 mutant residual = ", fmt[Simplify[m5MutantResidual]]];
If[Simplify[m5MutantResidual] === 0,
  Print["MUTANT FAIL: M5 leakage sign guard"]; Exit[1],
  Print["PASS: M5 mutant guard"]
];

Print["STATUS: PASS"];
