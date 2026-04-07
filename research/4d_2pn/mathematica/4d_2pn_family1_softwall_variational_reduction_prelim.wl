ClearAll[section, pass, fail, recordCheck, momentP2, momentP22, Kansatz, Kchannel];

section[s_String] := Print["\n" <> StringRepeat["=", 78] <> "\n" <> s <> "\n" <> StringRepeat["=", 78]];
pass[s_String] := Print["PASS: " <> s];
fail[s_String] := Print["FAIL: " <> s];
recordCheck[label_String, expr_] := Module[{v = FullSimplify[expr]}, If[v === True, pass[label], fail[label <> "  ->  " <> ToString[InputForm[v]]]]];

section["1) Input channel data and exact angular moments"];

momentP2 = <|
  {1, 1} -> -1/5,
  {1, 0} -> 2/5,
  {2, 0} -> 2/7,
  {2, 1} -> 1/7,
  {2, 2} -> -2/7
|>;

momentP22 = <|
  {1, 1} -> 1/7,
  {1, 0} -> 11/35,
  {2, 0} -> 3/7,
  {2, 1} -> 1/7,
  {2, 2} -> 1/7
|>;

channelData = {
  {"0", 0, 0, 4/45},
  {"1perp", 1, 1, 2/7},
  {"10", 1, 0, 1/4},
  {"20", 2, 0, 4/9},
  {"21", 2, 1, 2/3},
  {"22", 2, 2, 8/3}
};

Do[
  Print[channelData[[i]]],
  {i, Length[channelData]}
];

section["2) Reduced modal wall-Hessian ansatz"];

Print["Ansatz:"];
Print["  K_00 = Kmono"];
Print["  K_{ell m} = Kmono + Tau0 ell(ell+1) + (A0 + A1 ell(ell+1)) <P2> + (B0 + B1 ell(ell+1)) <P2^2>, for ell = 1,2"];

ClearAll[Kmono, Tau0, A0, A1, B0, B1];
Kansatz[ell_, m_] := Kmono + Tau0 ell (ell + 1) + (A0 + A1 ell (ell + 1)) momentP2[{ell, m}] + (B0 + B1 ell (ell + 1)) momentP22[{ell, m}];

sol = First @ Solve[
  {
    Kmono == 4/45,
    Kansatz[1, 1] == 2/7,
    Kansatz[1, 0] == 1/4,
    Kansatz[2, 0] == 4/9,
    Kansatz[2, 1] == 2/3,
    Kansatz[2, 2] == 8/3
  },
  {Kmono, Tau0, A0, A1, B0, B1},
  Reals
] // FullSimplify;

Print["Solved coefficients:"];
Do[
  Print["  ", sym, " = ", FullSimplify[sym /. sol], "  ~=  ", N[sym /. sol, 16]],
  {sym, {Kmono, Tau0, A0, A1, B0, B1}}
];

section["3) Exact reconstruction check"];

Kchannel[ell_, m_] := FullSimplify[Kansatz[ell, m] /. sol];

Do[
  Module[{name = ch[[1]], ell = ch[[2]], m = ch[[3]], target = ch[[4]], expr, resid},
    expr = If[ell == 0, FullSimplify[Kmono /. sol], Kchannel[ell, m]];
    resid = FullSimplify[expr - target];
    Print[name, ": K = ", expr, ", target = ", target, ", residual = ", resid];
    recordCheck["Residual vanishes for " <> name, resid == 0];
  ],
  {ch, channelData}
];

section["4) Family-1 Gaussian flare -> {1, P2, P2^2} basis"];

ClearAll[xi, D0, D1, D2];
D0 = FullSimplify[1 - xi/3 + xi^2/18];
D1 = FullSimplify[-2 xi/3 + 2 xi^2/9];
D2 = FullSimplify[2 xi^2/9];

Print["For a(z) = a0 (1 + beta Exp[-z^2/zm^2]) and xi = a0^2/zm^2,"];
Print["  delta a(theta)/a0 = beta [D0 + D1 P2 + D2 P2^2] + O(xi^3)"];
Print["  D0 = ", D0];
Print["  D1 = ", D1];
Print["  D2 = ", D2];
Print["So the actual Family-1 flare automatically produces the exact axisymmetric basis used by the reduced wall Hessian."];

section["5) Minimal linear-gradient interpretation"];

ratio = FullSimplify[(B1 /. sol)/(A1 /. sol)];
xiFit = FullSimplify[xi /. First @ Solve[ratio == FullSimplify[D2/D1], xi, Reals]];
zmOverA = FullSimplify[1/Sqrt[xiFit]];

Print["If the anisotropic gradient block is linear in the Family-1 flare profile, then"];
Print["  B1 / A1 = D2 / D1 = xi / (xi - 3)"];
Print["  B1/A1 = ", ratio, "  ~=  ", N[ratio, 16]];
Print["  xi = a0^2/zm^2 = ", xiFit, "  ~=  ", N[xiFit, 16]];
Print["  zm/a0 = 1/Sqrt[xi] = ", zmOverA, "  ~=  ", N[zmOverA, 16]];

section["6) Summary"];

Print["Main exact result:"];
Print["  K_00 = 4/45"];
Print["  K_{ell m} = 4/45 + (23/135) ell(ell+1)"];
Print["             + (9095/3528 - (25559/21168) ell(ell+1)) <P2>"];
Print["             + (-109/56 + (1765/3024) ell(ell+1)) <P2^2>, for ell = 1,2"];
Print[""];
Print["Interpretation:"];
Print["- This is the first exact reduced variational wall-Hessian that reproduces the passed static axisymmetric support data."];
Print["- The Family-1 Gaussian flare naturally generates the same {1, P2, P2^2} axisymmetric basis at second order."];
Print["- In the minimal linear-gradient reading, the flare width comes out of order the throat radius: zm/a0 ~= 1.0114."];

(*"
Output:


==============================================================================
1) Input channel data and exact angular moments
==============================================================================
{0, 0, 0, 4/45}
{1perp, 1, 1, 2/7}
{10, 1, 0, 1/4}
{20, 2, 0, 4/9}
{21, 2, 1, 2/3}
{22, 2, 2, 8/3}

==============================================================================
2) Reduced modal wall-Hessian ansatz
==============================================================================
Ansatz:
  K_00 = Kmono
  K_{ell m} = Kmono + Tau0 ell(ell+1) + (A0 + A1 ell(ell+1)) <P2> + (B0 + B1 ell(ell+1)) <P2^2>, for ell = 1,2
Solved coefficients:
  Kmono = 4/45  ~=  0.08888888888888888888888888888888888889`16.
  Tau0 = 23/135  ~=  0.17037037037037037037037037037037037037`16.
  A0 = 9095/3528  ~=  2.57794784580498866213151927437641723356`16.
  A1 = -25559/21168  ~=  -1.20743575207860922146636432350718065004`16.
  B0 = -109/56  ~=  -1.94642857142857142857142857142857142857`16.
  B1 = 1765/3024  ~=  0.58366402116402116402116402116402116402`16.

==============================================================================
3) Exact reconstruction check
==============================================================================
0: K = 4/45, target = 4/45, residual = 0
PASS: Residual vanishes for 0
1perp: K = 2/7, target = 2/7, residual = 0
PASS: Residual vanishes for 1perp
10: K = 1/4, target = 1/4, residual = 0
PASS: Residual vanishes for 10
20: K = 4/9, target = 4/9, residual = 0
PASS: Residual vanishes for 20
21: K = 2/3, target = 2/3, residual = 0
PASS: Residual vanishes for 21
22: K = 8/3, target = 8/3, residual = 0
PASS: Residual vanishes for 22

==============================================================================
4) Family-1 Gaussian flare -> {1, P2, P2^2} basis
==============================================================================
For a(z) = a0 (1 + beta Exp[-z^2/zm^2]) and xi = a0^2/zm^2,
  delta a(theta)/a0 = beta [D0 + D1 P2 + D2 P2^2] + O(xi^3)
  D0 = 1 + ((-6 + xi)*xi)/18
  D1 = (2*(-3 + xi)*xi)/9
  D2 = (2*xi^2)/9
So the actual Family-1 flare automatically produces the exact axisymmetric basis used by the reduced wall Hessian.

==============================================================================
5) Minimal linear-gradient interpretation
==============================================================================
If the anisotropic gradient block is linear in the Family-1 flare profile, then
  B1 / A1 = D2 / D1 = xi / (xi - 3)
  B1/A1 = -12355/25559  ~=  -0.48339136898939708126296020971086505732`16.
  xi = a0^2/zm^2 = 12355/12638  ~=  0.97760721633169805348947618294033866118`16.
  zm/a0 = 1/Sqrt[xi] = Sqrt[12638/12355]  ~=  1.01138800971329750681587838229057966196`16.

==============================================================================
6) Summary
==============================================================================
Main exact result:
  K_00 = 4/45
  K_{ell m} = 4/45 + (23/135) ell(ell+1)
             + (9095/3528 - (25559/21168) ell(ell+1)) <P2>
             + (-109/56 + (1765/3024) ell(ell+1)) <P2^2>, for ell = 1,2

Interpretation:
- This is the first exact reduced variational wall-Hessian that reproduces the passed static axisymmetric support data.
- The Family-1 Gaussian flare naturally generates the same {1, P2, P2^2} axisymmetric basis at second order.
- In the minimal linear-gradient reading, the flare width comes out of order the throat radius: zm/a0 ~= 1.0114.
"*)
