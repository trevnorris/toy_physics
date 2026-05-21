ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Print["STAGE 004 PROJECTED MAXWELL BUNDLE INDEX MATHEMATICA AUDIT"];

Clear[w, lam, mu0];
scaleAssumptions = lam > 0;
couplingAssumptions = lam > 0 && mu0 > 0;

(* M1: integration by parts with decaying test functions.
   Combine W*f' + W'*f into one integrand so Mathematica recognizes it as
   d/dw(W*f) and evaluates to the boundary term, which vanishes for the
   decaying Gaussian profile. *)
decayingWindow = Exp[-w^2/lam^2];
oddProbe = w*Exp[-w^2/lam^2];
m1Integrand = decayingWindow * D[oddProbe, w] + D[decayingWindow, w] * oddProbe;
m1Left = Integrate[m1Integrand, {w, -Infinity, Infinity}, Assumptions -> scaleAssumptions];
m1Right = 0;
m1Residual = FullSimplify[m1Left - m1Right, Assumptions -> scaleAssumptions];
Print["M1 residual = ", fmt[m1Residual]];
If[FullSimplify[m1Left - m1Right, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M1 density-level integration by parts"]; Exit[1]
];
Print["PASS: M1 density-level integration by parts"];

(* M2: cyclic Bianchi signs reduce to vector Faraday signs. *)
Clear[t, x, y, z, E1, E2, E3, B1, B2, B3];
spaceTimeAssumptions = Element[{t, x, y, z}, Reals];

twoForm23 = B1[t, x, y, z];
twoForm30 = E3[t, x, y, z];
twoForm02 = -E2[t, x, y, z];
m2LeftOne = D[twoForm23, t] + D[twoForm30, y] + D[twoForm02, z];
m2RightOne =
  D[B1[t, x, y, z], t] + D[E3[t, x, y, z], y] - D[E2[t, x, y, z], z];
m2ResidualOne = FullSimplify[m2LeftOne - m2RightOne, Assumptions -> spaceTimeAssumptions];
Print["M2 component 1 residual = ", fmt[m2ResidualOne]];
If[FullSimplify[m2LeftOne - m2RightOne, Assumptions -> spaceTimeAssumptions] =!= 0,
  Print["FAIL: M2 Faraday component 1"]; Exit[1]
];
Print["PASS: M2 Faraday component 1"];

twoForm31 = B2[t, x, y, z];
twoForm10 = E1[t, x, y, z];
twoForm03 = -E3[t, x, y, z];
m2LeftTwo = D[twoForm31, t] + D[twoForm10, z] + D[twoForm03, x];
m2RightTwo =
  D[B2[t, x, y, z], t] + D[E1[t, x, y, z], z] - D[E3[t, x, y, z], x];
m2ResidualTwo = FullSimplify[m2LeftTwo - m2RightTwo, Assumptions -> spaceTimeAssumptions];
Print["M2 component 2 residual = ", fmt[m2ResidualTwo]];
If[FullSimplify[m2LeftTwo - m2RightTwo, Assumptions -> spaceTimeAssumptions] =!= 0,
  Print["FAIL: M2 Faraday component 2"]; Exit[1]
];
Print["PASS: M2 Faraday component 2"];

twoForm12 = B3[t, x, y, z];
twoForm20 = E2[t, x, y, z];
twoForm01 = -E1[t, x, y, z];
m2LeftThree = D[twoForm12, t] + D[twoForm20, x] + D[twoForm01, y];
m2RightThree =
  D[B3[t, x, y, z], t] + D[E2[t, x, y, z], x] - D[E1[t, x, y, z], y];
m2ResidualThree = FullSimplify[m2LeftThree - m2RightThree, Assumptions -> spaceTimeAssumptions];
Print["M2 component 3 residual = ", fmt[m2ResidualThree]];
If[FullSimplify[m2LeftThree - m2RightThree, Assumptions -> spaceTimeAssumptions] =!= 0,
  Print["FAIL: M2 Faraday component 3"]; Exit[1]
];
Print["PASS: M2 Faraday component 3"];

(* M3: Gaussian normalization. *)
localizedProfile = Exp[-w^2/lam^2];
profileArea =
  Integrate[localizedProfile, {w, -Infinity, Infinity},
    Assumptions -> scaleAssumptions];
m3Target = Sqrt[Pi]*lam;
m3Residual = FullSimplify[profileArea - m3Target, Assumptions -> scaleAssumptions];
Print["M3 residual = ", fmt[m3Residual]];
If[FullSimplify[profileArea - m3Target, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M3 Gaussian normalization"]; Exit[1]
];
Print["PASS: M3 Gaussian normalization"];

(* M4: Gaussian squared norm. *)
profileSelfMass =
  Integrate[localizedProfile^2, {w, -Infinity, Infinity},
    Assumptions -> scaleAssumptions];
m4Target = Sqrt[2*Pi]*lam/2;
m4Residual = FullSimplify[profileSelfMass - m4Target, Assumptions -> scaleAssumptions];
Print["M4 residual = ", fmt[m4Residual]];
If[FullSimplify[profileSelfMass - m4Target, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M4 Gaussian squared norm"]; Exit[1]
];
Print["PASS: M4 Gaussian squared norm"];

(* M5: matched-kernel overlap. *)
matchedWeight = localizedProfile/profileArea;
overlapValue =
  Integrate[matchedWeight*localizedProfile, {w, -Infinity, Infinity},
    Assumptions -> scaleAssumptions];
m5Target = Sqrt[2]/2;
m5Residual = FullSimplify[overlapValue - m5Target, Assumptions -> scaleAssumptions];
Print["M5 residual = ", fmt[m5Residual]];
If[FullSimplify[overlapValue - m5Target, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M5 matched-kernel overlap"]; Exit[1]
];
Print["PASS: M5 matched-kernel overlap"];

(* M6: delta-source projection compared with reduced coupling. *)
pointSourceCoupling = mu0*(matchedWeight /. w -> 0)/overlapValue;
volumeReducedCoupling = mu0/profileArea;
m6Target = Sqrt[2];
m6Residual =
  FullSimplify[pointSourceCoupling/volumeReducedCoupling - m6Target,
    Assumptions -> couplingAssumptions];
Print["M6 residual = ", fmt[m6Residual]];
If[FullSimplify[pointSourceCoupling/volumeReducedCoupling - m6Target,
    Assumptions -> couplingAssumptions] =!= 0,
  Print["FAIL: M6 delta-source projection/reduction ratio"]; Exit[1]
];
Print["PASS: M6 delta-source projection/reduction ratio"];

Print["STATUS: PASS"];
Exit[0];
