(* 4d_2pn_axisymmetric_robin_wall_pde_prelim.wl *)

ClearAll["Global`*"];

section[s_String] := Print["\n=== ", s, " ==="];
pass[name_String, cond_] := Print[If[TrueQ[cond], "PASS: ", "FAIL: "], name];
show[name_String, expr_] := Print[name, " = ", expr];

approxEq[a_, b_, tol_: 10^-20] := Abs[N[a - b, 80]] < tol;

doubleFactorialSafe[n_Integer] := If[n <= 0, 1, n!!];

sphereMoment[a_Integer, b_Integer, c_Integer] :=
  doubleFactorialSafe[2 a - 1]*
  doubleFactorialSafe[2 b - 1]*
  doubleFactorialSafe[2 c - 1]/
  doubleFactorialSafe[2 (a + b + c) + 1];

polyAvg[poly_] := Module[{rules},
  rules = CoefficientRules[Expand[poly], {X, Y, Z}];
  FullSimplify @ Total[(Last[#] * Apply[sphereMoment, First[#]]) & /@ rules]
];

section["Axisymmetric P2 and P2^2 expectation values on the raw mouth basis"];

P2 = (3 Z - 1)/2;

weightRules = {
  "1perp" -> X,
  "10"    -> Z,
  "20"    -> ((3 Z - 1)/2)^2,
  "21"    -> 3 X Z,
  "22"    -> 3 X Y
};

cVals = Association @ Map[
   First[#] -> FullSimplify[polyAvg[Last[#] * P2]/polyAvg[Last[#]]] &,
   weightRules
];
dVals = Association @ Map[
   First[#] -> FullSimplify[polyAvg[Last[#] * P2^2]/polyAvg[Last[#]]] &,
   weightRules
];

show["cVals", cVals];
show["dVals", dVals];

pass["<P2> dipole and quadrupole expectations",
  And[
    cVals["1perp"] === -1/5,
    cVals["10"] === 2/5,
    cVals["20"] === 2/7,
    cVals["21"] === 1/7,
    cVals["22"] === -2/7
  ]
];

pass["<P2^2> expectations",
  And[
    dVals["1perp"] === 1/7,
    dVals["10"] === 11/35,
    dVals["20"] === 3/7,
    dVals["21"] === 1/7,
    dVals["22"] === 1/7
  ]
];

section["Static 2PN raw support targets"];

yTarget = <|
  "0" -> 45/4,
  "1perp" -> 7/2,
  "10" -> 4,
  "20" -> 9/4,
  "21" -> 3/2,
  "22" -> 3/8
|>;

zTarget = Association @ KeyValueMap[#1 -> FullSimplify[1/#2] &, yTarget];

show["zTarget", zTarget];

section["Odd sector: exact first-order P2 wall fit"];

solOdd = First @ Solve[
   {
    b1 + a1 cVals["1perp"] == zTarget["1perp"],
    b1 + a1 cVals["10"] == zTarget["10"]
   },
   {b1, a1}
];

show["solOdd", solOdd];

pass["Odd first-order P2 fit exact",
  And[
    FullSimplify[(b1 + a1 cVals["1perp"]) /. solOdd] === zTarget["1perp"],
    FullSimplify[(b1 + a1 cVals["10"]) /. solOdd] === zTarget["10"],
    (b1 /. solOdd) === 23/84,
    (a1 /. solOdd) === -5/84
  ]
];

section["Even sector: first-order P2 no-go"];

solEven1Full = Solve[
   {
    b2 + a2 cVals["20"] == zTarget["20"],
    b2 + a2 cVals["21"] == zTarget["21"],
    b2 + a2 cVals["22"] == zTarget["22"]
   },
   {b2, a2}
];

solEven1Partial = First @ Solve[
   {
    b2 + a2 cVals["20"] == zTarget["20"],
    b2 + a2 cVals["21"] == zTarget["21"]
   },
   {b2, a2}
];

predZ22First = FullSimplify[(b2 + a2 cVals["22"]) /. solEven1Partial];

show["solEven1Partial", solEven1Partial];
show["predZ22First", predZ22First];

pass["No exact first-order quadrupole fit exists", solEven1Full === {}];
pass["First-order quadrupole mismatch is 4/3 vs 8/3",
  predZ22First === 4/3 && zTarget["22"] === 8/3
];

section["Even sector: exact second-order P2^2 wall fit"];

solEven2 = First @ Solve[
   {
    b2 + a2 cVals["20"] + c2w dVals["20"] == zTarget["20"],
    b2 + a2 cVals["21"] + c2w dVals["21"] == zTarget["21"],
    b2 + a2 cVals["22"] + c2w dVals["22"] == zTarget["22"]
   },
   {b2, a2, c2w}
];

show["solEven2", solEven2];

pass["Second-order quadrupole fit exact",
  And[
    FullSimplify[(b2 + a2 cVals["20"] + c2w dVals["20"]) /. solEven2] === zTarget["20"],
    FullSimplify[(b2 + a2 cVals["21"] + c2w dVals["21"]) /. solEven2] === zTarget["21"],
    FullSimplify[(b2 + a2 cVals["22"] + c2w dVals["22"]) /. solEven2] === zTarget["22"],
    (b2 /. solEven2) === 10/9,
    (a2 /. solEven2) === -14/3,
    (c2w /. solEven2) === 14/9
  ]
];

section["Monopole stiffness and axisymmetric source"];

b0Val = zTarget["0"];

sourceSol = First @ Solve[
   {
    (2/Sqrt[5]) (pIso + qAx/3) == 4/Sqrt[5],
    (2/3) qAx == 5/4
   },
   {pIso, qAx}
];

show["b0Val", b0Val];
show["sourceSol", sourceSol];

pass["Finite monopole stiffness exact", b0Val === 4/45];
pass["Axisymmetric source profile exact",
  And[
    (pIso /. sourceSol) === 11/8,
    (qAx /. sourceSol) === 15/8
  ]
];

section["4D-ball DtN pole equations"];

Lambda4Ball[l_Integer, z_?NumericQ] := Module[{nu = l + 1},
  -1 + z (BesselJ[nu - 1, z] - BesselJ[nu + 1, z])/(2 BesselJ[nu, z])
];

LambdaPrime[l_Integer, z_?NumericQ] := N[
  D[-1 + zz (BesselJ[l, zz] - BesselJ[l + 2, zz])/(2 BesselJ[l + 1, zz]), zz] /. zz -> z,
  80
];

b0Num = N[b0Val, 80];
b1Base = N[b1 /. solOdd, 80];
a1Num = N[a1 /. solOdd, 80];
b2Base = N[b2 /. solEven2, 80];
a2Num = N[a2 /. solEven2, 80];
c2Num = N[c2w /. solEven2, 80];

z0Exact = zz /. FindRoot[Lambda4Ball[0, zz] + b0Num == 0, {zz, 0.6},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];

z1Base = zz /. FindRoot[Lambda4Ball[1, zz] + b1Base == 0, {zz, 2.55},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
z1PerpExact = zz /. FindRoot[Lambda4Ball[1, zz] + N[zTarget["1perp"], 80] == 0, {zz, 2.56},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
z10Exact = zz /. FindRoot[Lambda4Ball[1, zz] + N[zTarget["10"], 80] == 0, {zz, 2.53},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];

z2Base = zz /. FindRoot[Lambda4Ball[2, zz] + b2Base == 0, {zz, 4.25},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
z20Exact = zz /. FindRoot[Lambda4Ball[2, zz] + N[zTarget["20"], 80] == 0, {zz, 3.90},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
z21Exact = zz /. FindRoot[Lambda4Ball[2, zz] + N[zTarget["21"], 80] == 0, {zz, 4.03},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
z22Exact = zz /. FindRoot[Lambda4Ball[2, zz] + N[zTarget["22"], 80] == 0, {zz, 4.82},
   WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];

z1PerpPert = N[z1Base - (a1Num N[cVals["1perp"], 80])/LambdaPrime[1, z1Base], 50];
z10Pert = N[z1Base - (a1Num N[cVals["10"], 80])/LambdaPrime[1, z1Base], 50];

z20Pert = N[z2Base - (a2Num N[cVals["20"], 80] + c2Num N[dVals["20"], 80])/LambdaPrime[2, z2Base], 50];
z21Pert = N[z2Base - (a2Num N[cVals["21"], 80] + c2Num N[dVals["21"], 80])/LambdaPrime[2, z2Base], 50];
z22Pert = N[z2Base - (a2Num N[cVals["22"], 80] + c2Num N[dVals["22"], 80])/LambdaPrime[2, z2Base], 50];

z0Small = N[4/Sqrt[45], 50];

show["z0Exact", z0Exact];
show["z0Small", z0Small];
show["z1Base", z1Base];
show["z1PerpExact", z1PerpExact];
show["z1PerpPert", z1PerpPert];
show["z10Exact", z10Exact];
show["z10Pert", z10Pert];
show["z2Base", z2Base];
show["z20Exact", z20Exact];
show["z20Pert", z20Pert];
show["z21Exact", z21Exact];
show["z21Pert", z21Pert];
show["z22Exact", z22Exact];
show["z22Pert", z22Pert];

pass["Finite monopole pole exists and small-z estimate is close",
  z0Exact > 0 && Abs[N[z0Exact - z0Small, 50]] < 0.01
];

pass["Dipole pole splitting exists",
  z1PerpExact > z10Exact > 0
];

pass["Odd perturbative root shifts are accurate",
  approxEq[z1PerpPert, z1PerpExact, 5*10^-4] &&
  approxEq[z10Pert, z10Exact, 5*10^-4]
];

pass["Quadrupole pole ordering exists",
  z20Exact < z21Exact < z22Exact
];

pass["Quadrupole perturbative ordering is at least qualitatively correct",
  z20Pert < z21Pert < z22Pert
];

section["Summary"];

Print[
  "Minimal static throat-wall PDE scaffold:\n",
  "  Z0(0) = 4/45,\n",
  "  Z1m(0) = 23/84 - (5/84) <P2>_(1m),\n",
  "  Z2m(0) = 10/9 - (14/3) <P2>_(2m) + (14/9) <P2^2>_(2m).\n",
  "This reproduces the full raw static support data, gives a finite monopole pole, ",
  "and yields explicit axisymmetric dipole/quadrupole pole splittings in the 4D-ball DtN model."
];
