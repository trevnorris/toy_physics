(* ::Package:: *)

ClearAll["Global`*"];

section[s_] := Print["\n=== ", s, " ==="];

(* Support-channel data from the passed Robin-wall step *)
supportData = {
  <|"Name" -> "1perp", "l" -> 1, "Lambda" -> 2, "P2" -> -1/5, "P22" -> 1/7,   "Z" -> 2/7|>,
  <|"Name" -> "10",    "l" -> 1, "Lambda" -> 2, "P2" ->  2/5, "P22" -> 11/35, "Z" -> 1/4|>,
  <|"Name" -> "20",    "l" -> 2, "Lambda" -> 6, "P2" ->  2/7, "P22" -> 3/7,   "Z" -> 4/9|>,
  <|"Name" -> "21",    "l" -> 2, "Lambda" -> 6, "P2" ->  1/7, "P22" -> 1/7,   "Z" -> 2/3|>,
  <|"Name" -> "22",    "l" -> 2, "Lambda" -> 6, "P2" -> -2/7, "P22" -> 1/7,   "Z" -> 8/3|>
};

section["Minimal generalized Robin / wall-stress support operator"];

eqns = Table[
   With[{Lm2 = d["Lambda"] - 2},
    A0 + A1*Lm2 + (B0 + B1*Lm2)*d["P2"] + C*Lm2*d["P22"] == d["Z"]
   ],
   {d, supportData}
];

sol = First @ Solve[eqns, {A0, A1, B0, B1, C}] // FullSimplify;

Print["Ansatz:"];
Print["  Z_lm = A0 + A1 (lambda_l-2) + (B0 + B1 (lambda_l-2)) <P2>_lm + C (lambda_l-2) <P2^2>_lm"];
Print[""];
Print["Solved coefficients:"];
Print["  A0 = ", A0 /. sol];
Print["  A1 = ", A1 /. sol];
Print["  B0 = ", B0 /. sol];
Print["  B1 = ", B1 /. sol];
Print["  C  = ", C  /. sol];

section["Exact verification on passed support data"];

verifyRows = Table[
   With[
    {
      pred = FullSimplify[
        (A0 + A1*(d["Lambda"] - 2) + (B0 + B1*(d["Lambda"] - 2))*d["P2"] + C*(d["Lambda"] - 2)*d["P22"]) /. sol
      ]
    },
    <|
      "Name" -> d["Name"],
      "PredictedZ" -> pred,
      "TargetZ" -> d["Z"],
      "Residual" -> FullSimplify[pred - d["Z"]]
    |>
   ],
   {d, supportData}
];

Do[
  Print[row["Name"], ": predicted Z = ", row["PredictedZ"], ", target Z = ", row["TargetZ"], ", residual = ", row["Residual"]],
  {row, verifyRows}
];

section["Equivalent profile decomposition in mu = Cos[theta]"];

mu = Symbol["mu"];
P2mu = (3 mu^2 - 1)/2;

zBase = FullSimplify[(A0 + B0 P2mu) /. sol] // Expand;
zCurv = FullSimplify[(A1 + B1 P2mu + C P2mu^2) /. sol] // Expand;

Print["Equivalent diagonal form on l=1,2:"];
Print["  Z_lm = < zBase(mu) >_lm + (lambda_l-2) < zCurv(mu) >_lm"];
Print[""];
Print["Base profile:"];
Print["  zBase(mu) = ", zBase];
Print[""];
Print["Curvature profile:"];
Print["  zCurv(mu) = ", zCurv];

section["Family-1 flare / soft-wall mouth basis"];

h1 = -mu^2;
h2 = mu^4;

baseFit = First @ SolveAlways[c0 + c2 h1 == zBase, mu] // FullSimplify;
curvFit = First @ SolveAlways[d0 + d2 h1 + d4 h2 == zCurv, mu] // FullSimplify;

Print["Use h1(mu) = -mu^2 and h2(mu) = mu^4."];
Print[""];
Print["Base profile in {1, h1}:"];
Print["  zBase(mu) = ", FullSimplify[c0 /. baseFit], " + ", FullSimplify[c2 /. baseFit], " h1(mu)"];
Print[""];
Print["Curvature profile in {1, h1, h2}:"];
Print["  zCurv(mu) = ", FullSimplify[d0 /. curvFit], " + ", FullSimplify[d2 /. curvFit], " h1(mu) + ", FullSimplify[d4 /. curvFit], " h2(mu)"];

section["Family-1 flare expansion near the mouth"];

(* a(z)/a_m = 1 - q mu^2 + r mu^4 + O(mu^6) *)
flareExpr = -q mu^2 + r mu^4;

flareFit = First @ SolveAlways[a + b P2mu + c P2mu^2 == flareExpr, mu] // FullSimplify;

Print["Assume a(z)/a_m = 1 - q mu^2 + r mu^4 + O(mu^6), with mu = Cos[theta]."];
Print[""];
Print["Then"];
Print["  -q mu^2 + r mu^4 = ", FullSimplify[a /. flareFit], " + ", FullSimplify[b /. flareFit], " P2(mu) + ", FullSimplify[c /. flareFit], " P2(mu)^2"];

section["Axisymmetric source profile in the same mu^2 basis"];

sourceP2 = 11/8 + (15/8) P2mu // FullSimplify // Expand;
sourceFit = First @ SolveAlways[s0 + s2 mu^2 == sourceP2, mu] // FullSimplify;

Print["S(mu) = 11/8 + (15/8) P2(mu) = ", sourceP2];
Print["     = ", FullSimplify[s0 /. sourceFit], " + ", FullSimplify[s2 /. sourceFit], " mu^2"];

section["Suggested minimal support-sector boundary operator"];

Print["Suggested diagonal-matrix-element form on the passed l=1,2 support sector:"];
Print["  B_eff = d_n + zBase(mu) + 1/2 {(-Delta_S - 2), zCurv(mu)}"];
Print[""];
Print["with"];
Print["  zBase(mu) = ", zBase];
Print["  zCurv(mu) = ", zCurv];
Print[""];
Print["This exactly reproduces the passed support impedances for the l=1,2 channels."];
Print[""];
Print["Note: the previously reported z^8 correction in the isotropic 4D-ball unit-test series does not affect this static operator reconstruction, which uses the exact support-channel data rather than the old truncated z^8 quotient."];
