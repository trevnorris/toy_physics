(* Projected Maxwell -> P2 compatibility bridge. *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

check[label_, expr_] := If[
  FullSimplify[expr, Assumptions -> $Assumptions] =!= 0,
  Print["FAIL: ", label, " residual = ", fmt[FullSimplify[expr, Assumptions -> $Assumptions]]]; Exit[1],
  Print[label, " residual = ", fmt[FullSimplify[expr, Assumptions -> $Assumptions]]]
];

Print["STAGE 011 PROJECTED MAXWELL P2 BRIDGE MATHEMATICA AUDIT"];

Clear[K, B0, Z0slot, N0, Ptarget, S, T, z0, z2, z4, n0, eps];

$Assumptions =
  Element[{K, B0, Z0slot, N0, Ptarget, S, T, z0, z2, z4, n0, eps}, Reals] &&
    T != 0 && Ptarget != 0;

firstShift[expr_] :=
  Coefficient[Normal[Series[expr - (expr /. eps -> 0), {eps, 0, 1}]], eps, 1];

kPole =
  K /. First[
    Solve[(K - B0 - (Z0slot + eps z0)) (T + eps z4) ==
      3 (S + eps z2)^2, K]
  ];
kNorm =
  K /. First[
    Solve[(N0 + eps n0)/(K - B0 - (Z0slot + eps z0)) == Ptarget, K]
  ];

compatFixed = FullSimplify[kNorm - kPole, Assumptions -> $Assumptions];
compatFixedClosed =
  (N0 + eps n0)/Ptarget - 3 (S + eps z2)^2/(T + eps z4);
check[
  "M1 fixed-target K-eliminated surface",
  compatFixed - compatFixedClosed
];

deltaKPole = firstShift[kPole];
check[
  "M2 one-pole K shift",
  deltaKPole - (z0 + 6 S z2/T - 3 S^2 z4/T^2)
];

deltaKNorm = firstShift[kNorm];
check[
  "M3 normalization K shift",
  deltaKNorm - (z0 + n0/Ptarget)
];

compatFixedShift = firstShift[compatFixed];
check[
  "M4 fixed-target compatibility first variation",
  compatFixedShift - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2)
];
check["M5 fixed-target z0 independence", D[compatFixedShift, z0]];

Print["STATUS: PASS"];
Exit[0];
