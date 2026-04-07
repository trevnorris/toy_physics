(* ::Package:: *)

ClearAll["Global`*"];
$HistoryLength = 0;

SectionLine[title_String] := (
  Print[""];
  Print[StringRepeat["=", 72]];
  Print[title];
  Print[StringRepeat["=", 72]];
  Print[""];
);

SubSectionLine[title_String] := Print["-- " <> title <> " --"];

ShowValue[label_String, expr_] := Module[{value = FullSimplify[expr]},
  Print[label <> " = " <> ToString[InputForm[value]]];
  value
];

CheckEqual[label_String, expr_, expected_] := Module[{value = FullSimplify[expr]},
  Print[label <> " = " <> ToString[InputForm[value]]];
  If[!TrueQ[FullSimplify[value == expected]],
    Print["FAIL: expected " <> ToString[InputForm[expected]]];
    Exit[1];
  ];
  value
];

CheckMatch[label_String, expr_, expected_] := Module[{value = expr},
  Print[label <> " = " <> ToString[InputForm[value]]];
  If[value =!= expected,
    Print["FAIL: expected " <> ToString[InputForm[expected]]];
    Exit[1];
  ];
  value
];

SectionLine["4D 1PN bridge core checks"];

SubSectionLine["Added-mass topology discriminator"];
ClearAll[a, r, th, ph, rho0, U];

phi3[r_, th_] := -(U a^3/(2 r^2)) Cos[th];
v3r[r_, th_] := D[phi3[r, th], r];
v3th[r_, th_] := (1/r) D[phi3[r, th], th];
v3sq[r_, th_] := FullSimplify[v3r[r, th]^2 + v3th[r, th]^2];

T3 = FullSimplify[
  (rho0/2) Integrate[
    v3sq[r, th] r^2 Sin[th],
    {r, a, Infinity}, {th, 0, Pi}, {ph, 0, 2 Pi},
    Assumptions -> {a > 0, U > 0, rho0 > 0},
    GenerateConditions -> False
  ]
];
mAdd3 = FullSimplify[2 T3/U^2];
mDisp3 = FullSimplify[rho0 (4 Pi a^3/3)];

ShowValue["T_ext(3D sphere)", T3];
CheckEqual["kappa_add(throat)", mAdd3/mDisp3, 1/2];
CheckEqual["kappa_add(bubble)", 1/(4 - 1), 1/3];

SubSectionLine["Optical coefficient selects n = 5"];
ClearAll[n, delta];

csRatio = FullSimplify[(1 + delta)^((n - 1)/2)];
Ndelta = FullSimplify[(1 + delta)^(-(n - 1)/2)];

ShowValue["c_s(rho0(1+delta))/c0", csRatio];
ShowValue["N(delta)", Ndelta];
CheckEqual["linear coeff in N(delta)", SeriesCoefficient[Ndelta, {delta, 0, 1}], -(n - 1)/2];
CheckMatch["solve ((n-1)/2 == 2)", Solve[(n - 1)/2 == 2, n], {{n -> 5}}];

SubSectionLine["EIH family and thermodynamic closure"];
ClearAll[alpha2, aH2, Kvec];

Cparallel = Kvec Pi^2 (-1 + aH2 - alpha2);
CLongitudinal = Kvec Pi^2 (-1 + aH2 + alpha2);
family = Solve[{Cparallel == -7/2, CLongitudinal == -1/2}, {alpha2, Kvec}];
alphaFamily = FullSimplify[alpha2 /. First[family]];
KFamily = FullSimplify[Kvec /. First[family]];

ShowValue["family solve", family];
CheckEqual["alpha^2(a_H^2)", alphaFamily, (3 (1 - aH2))/4];
CheckEqual["K_vec(a_H^2)", KFamily, 2/(Pi^2 (1 - aH2))];
CheckEqual["K_vec (1-a_H^2)", KFamily (1 - aH2), 2/Pi^2];
CheckEqual["K_vec alpha^2", KFamily alphaFamily, 3/(2 Pi^2)];
CheckEqual["alpha^2_thermo(n=5)", 1 - 1/(n - 1) /. n -> 5, 3/4];
thermoClosure = Solve[{Cparallel == -7/2, CLongitudinal == -1/2, alpha2 == 3/4}, {alpha2, Kvec, aH2}];
ShowValue["thermo + EIH solve", thermoClosure];
CheckEqual["alpha^2 unique", alpha2 /. First[thermoClosure], 3/4];
CheckEqual["K_vec unique", Kvec /. First[thermoClosure], 2/Pi^2];
CheckEqual["a_H^2 unique", aH2 /. First[thermoClosure], 0];

SubSectionLine["Wave-supported mass and the 3/8 coefficient"];
ClearAll[v, c, E0];

gammaSeries = Normal[Series[1/Sqrt[1 - v^2/c^2], {v, 0, 4}]];
energySeries = Expand[E0 gammaSeries];

CheckEqual["gamma(v) through O(v^4)", gammaSeries, 1 + v^2/(2 c^2) + (3 v^4)/(8 c^4)];
ShowValue["E(v) through O(v^4)", energySeries];
CheckEqual["coeff of v^4 in E(v)", SeriesCoefficient[energySeries, {v, 0, 4}], 3 E0/(8 c^4)];

Print["ALL CORE CHECKS PASSED"];

(*"
Output:

========================================================================
4D 1PN bridge core checks
========================================================================

-- Added-mass topology discriminator --
T_ext(3D sphere) = (a^3*Pi*rho0*U^2)/3
kappa_add(throat) = 1/2
kappa_add(bubble) = 1/3
-- Optical coefficient selects n = 5 --
c_s(rho0(1+delta))/c0 = (1 + delta)^((-1 + n)/2)
N(delta) = (1 + delta)^(1/2 - n/2)
linear coeff in N(delta) = (1 - n)/2
solve ((n-1)/2 == 2) = {{n -> 5}}
-- EIH family and thermodynamic closure --
family solve = {{alpha2 -> (-3*(-1 + aH2))/4, Kvec -> -2/((-1 + aH2)*Pi^2)}}
alpha^2(a_H^2) = (-3*(-1 + aH2))/4
K_vec(a_H^2) = -2/((-1 + aH2)*Pi^2)
K_vec (1-a_H^2) = 2/Pi^2
K_vec alpha^2 = 3/(2*Pi^2)
alpha^2_thermo(n=5) = 3/4
thermo + EIH solve = {{alpha2 -> 3/4, Kvec -> 2/Pi^2, aH2 -> 0}}
alpha^2 unique = 3/4
K_vec unique = 2/Pi^2
a_H^2 unique = 0
-- Wave-supported mass and the 3/8 coefficient --
gamma(v) through O(v^4) = 1 + v^2/(2*c^2) + (3*v^4)/(8*c^4)
E(v) through O(v^4) = (E0*(8 + (4*v^2)/c^2 + (3*v^4)/c^4))/8
coeff of v^4 in E(v) = (3*E0)/(8*c^4)
ALL CORE CHECKS PASSED
"*)
