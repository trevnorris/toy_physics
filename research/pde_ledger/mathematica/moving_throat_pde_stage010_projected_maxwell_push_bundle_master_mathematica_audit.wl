(* moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Print["STAGE 010 PROJECTED MAXWELL PUSH BUNDLE MASTER MATHEMATICA AUDIT"];

Clear[
  eps, D0, D2, D4, N0, N2, N4, z0, z2, z4, n0, n2, n4,
  K, B0, Z0slot, Ptarget, S, T, D0target, x0, x1,
  Q, Delta, P, q1, d1, p1, D0sym, N0sym
];
$Assumptions =
  Element[
    {eps, D0, D2, D4, N0, N2, N4, z0, z2, z4, n0, n2, n4,
      K, B0, Z0slot, Ptarget, S, T, D0target, x0, x1,
      Q, Delta, P, q1, d1, p1, D0sym, N0sym},
    Reals
  ] && D0 != 0 && T != 0 && Ptarget != 0 && D0target != 0 &&
    Delta != 0 && P != 0 && D0sym != 0 && N0sym != 0;

den0[e_] := D0 - e z0;
den2[e_] := D2 - e z2;
den4[e_] := D4 - e z4;
num0[e_] := N0 + e n0;
num2[e_] := N2 + e n2;
num4[e_] := N4 + e n4;

u2slot[e_] := -den2[e]/den0[e];
u4slot[e_] := den2[e]^2/den0[e]^2 - den4[e]/den0[e];
slot0[e_] := num0[e]/den0[e];
slot2[e_] :=
  (den0[e] num2[e] - 2 den2[e] num0[e])/den0[e]^2;
slot4[e_] :=
  (den0[e]^2 num4[e] -
      2 den0[e] (den2[e] num2[e] + den4[e] num0[e]) +
      3 den2[e]^2 num0[e])/den0[e]^3;

u2Linear = Coefficient[Normal[Series[u2slot[eps], {eps, 0, 1}]], eps, 1];
u4Linear = Coefficient[Normal[Series[u4slot[eps], {eps, 0, 1}]], eps, 1];
slot0Linear = Coefficient[Normal[Series[slot0[eps], {eps, 0, 1}]], eps, 1];
slot2Linear = Coefficient[Normal[Series[slot2[eps], {eps, 0, 1}]], eps, 1];
slot4Linear = Coefficient[Normal[Series[slot4[eps], {eps, 0, 1}]], eps, 1];

m0aResidual = FullSimplify[u2Linear - (D0 z2 - D2 z0)/D0^2];
Print["M0a residual = ", fmt[m0aResidual]];
If[FullSimplify[m0aResidual] =!= 0, Print["FAIL: M0a"]; Exit[1]];

m0bResidual =
  FullSimplify[
    u4Linear -
      (D0^2 z4 - D0 (2 D2 z2 + D4 z0) + 2 D2^2 z0)/D0^3
  ];
Print["M0b residual = ", fmt[m0bResidual]];
If[FullSimplify[m0bResidual] =!= 0, Print["FAIL: M0b"]; Exit[1]];

m1Residual =
  FullSimplify[slot0Linear - (n0/D0 + N0 z0/D0^2)];
Print["M1 residual = ", fmt[m1Residual]];
If[FullSimplify[m1Residual] =!= 0, Print["FAIL: M1"]; Exit[1]];

m2Target =
  n2/D0 + N2 z0/D0^2 + 2 N0 z2/D0^2 -
    2 D2 n0/D0^2 - 4 D2 N0 z0/D0^3;
m2Residual = FullSimplify[slot2Linear - m2Target];
Print["M2 residual = ", fmt[m2Residual]];
If[FullSimplify[m2Residual] =!= 0, Print["FAIL: M2"]; Exit[1]];

m3Target =
  n4/D0 + N4 z0/D0^2 + 2 N2 z2/D0^2 -
    2 D2 n2/D0^2 + 2 N0 z4/D0^2 - 2 D4 n0/D0^2 -
    4 (D2 N2 + D4 N0) z0/D0^3 -
    6 D2 N0 z2/D0^3 + 3 D2^2 n0/D0^3 +
    9 D2^2 N0 z0/D0^4;
m3Residual = FullSimplify[slot4Linear - m3Target];
Print["M3 residual = ", fmt[m3Residual]];
If[FullSimplify[m3Residual] =!= 0, Print["FAIL: M3"]; Exit[1]];

poleEquation =
  (K - B0 - Z0slot - eps z0) (T + eps z4) ==
    3 (S + eps z2)^2;
poleSolutions = Solve[poleEquation, K];
If[FullSimplify[Length[poleSolutions] - 1] =!= 0,
  Print["FAIL: M4"]; Exit[1]
];
onePoleSurface = K /. First[poleSolutions];
onePoleClosed =
  B0 + Z0slot + eps z0 + 3 (S + eps z2)^2/(T + eps z4);
m4Residual = FullSimplify[onePoleSurface - onePoleClosed];
Print["M4 residual = ", fmt[m4Residual]];
If[FullSimplify[m4Residual] =!= 0, Print["FAIL: M4"]; Exit[1]];

onePoleShift =
  Coefficient[Normal[Series[onePoleSurface, {eps, 0, 1}]], eps, 1];
m5Residual =
  FullSimplify[onePoleShift - (z0 + 6 S z2/T - 3 S^2 z4/T^2)];
Print["M5 residual = ", fmt[m5Residual]];
If[FullSimplify[m5Residual] =!= 0, Print["FAIL: M5"]; Exit[1]];

normEquation =
  (N0 + eps n0)/(K - B0 - Z0slot - eps z0) == Ptarget;
normSolutions = Solve[normEquation, K];
If[FullSimplify[Length[normSolutions] - 1] =!= 0,
  Print["FAIL: M6"]; Exit[1]
];
normSurface = K /. First[normSolutions];
normClosed = B0 + Z0slot + eps z0 + (N0 + eps n0)/Ptarget;
m6Residual = FullSimplify[normSurface - normClosed];
Print["M6 residual = ", fmt[m6Residual]];
If[FullSimplify[m6Residual] =!= 0, Print["FAIL: M6"]; Exit[1]];

normShift =
  Coefficient[Normal[Series[normSurface, {eps, 0, 1}]], eps, 1];
m7Residual = FullSimplify[normShift - (z0 + n0/Ptarget)];
Print["M7 residual = ", fmt[m7Residual]];
If[FullSimplify[m7Residual] =!= 0, Print["FAIL: M7"]; Exit[1]];

compatSurface = FullSimplify[normSurface - onePoleSurface];
compatDirect =
  (N0 + eps n0)/Ptarget - 3 (S + eps z2)^2/(T + eps z4);
m8Residual = FullSimplify[compatSurface - compatDirect];
Print["M8 residual = ", fmt[m8Residual]];
If[FullSimplify[m8Residual] =!= 0, Print["FAIL: M8"]; Exit[1]];

compatShift =
  Coefficient[Normal[Series[compatSurface, {eps, 0, 1}]], eps, 1];
m9Residual =
  FullSimplify[compatShift - (n0/Ptarget - 6 S z2/T + 3 S^2 z4/T^2)];
Print["M9 residual = ", fmt[m9Residual]];
If[FullSimplify[m9Residual] =!= 0, Print["FAIL: M9"]; Exit[1]];

transportTarget[e_] := (N0 + e n0)/D0target;
transportEquation =
  (N0 + eps n0)/(K - B0 - Z0slot - eps z0) ==
    transportTarget[eps];
transportSolutions = Solve[transportEquation, K];
transportSurface = K /. First[transportSolutions];
m10Residual =
  FullSimplify[transportSurface - (B0 + Z0slot + eps z0 + D0target)];
Print["M10 residual = ", fmt[m10Residual]];
If[FullSimplify[m10Residual] =!= 0, Print["FAIL: M10"]; Exit[1]];

transportCompat = FullSimplify[transportSurface - onePoleSurface];
m11Residual =
  FullSimplify[
    transportCompat -
      (D0target - 3 (S + eps z2)^2/(T + eps z4))
  ];
Print["M11 residual = ", fmt[m11Residual]];
If[FullSimplify[m11Residual] =!= 0, Print["FAIL: M11"]; Exit[1]];

transportShift =
  Coefficient[Normal[Series[transportCompat, {eps, 0, 1}]], eps, 1];
m12Residual =
  FullSimplify[transportShift - (-6 S z2/T + 3 S^2 z4/T^2)];
Print["M12 residual = ", fmt[m12Residual]];
If[FullSimplify[m12Residual] =!= 0, Print["FAIL: M12"]; Exit[1]];

gauntByThreeJ[l1_, l2_, l3_, m1_, m2_, m3_] :=
  If[m1 + m2 + m3 != 0,
    0,
    Sqrt[(2 l1 + 1) (2 l2 + 1) (2 l3 + 1)/(4 Pi)] *
      ThreeJSymbol[{l1, 0}, {l2, 0}, {l3, 0}] *
      ThreeJSymbol[{l1, m1}, {l2, m2}, {l3, m3}]
  ];
lambda2[m_] :=
  FullSimplify[
    (-1)^m gauntByThreeJ[2, 2, 2, 0, m, -m]/
      gauntByThreeJ[2, 2, 2, 0, 0, 0]
  ];
m13Residuals = FullSimplify[
  {
    lambda2[0] - 1,
    lambda2[1] - 1/2,
    lambda2[2] + 1,
    gauntByThreeJ[2, 2, 2, 0, 1, 1],
    gauntByThreeJ[2, 2, 2, 0, 2, 2]
  }
];
Print["M13 residuals = ", fmt[m13Residuals]];
If[FullSimplify[m13Residuals] =!= {0, 0, 0, 0, 0},
  Print["FAIL: M13"]; Exit[1]
];

lane0 = x0 + eps lambda2[0] x1;
lane1 = x0 + eps lambda2[1] x1;
lane2 = x0 + eps lambda2[2] x1;
meanTrace = FullSimplify[(lane0 + 2 lane1 + 2 lane2)/5];
axisTrace = FullSimplify[(2 lane0 - lane1 - lane2)/10];
branchTrace = FullSimplify[(lane1 - lane2)/2];
m14Residuals =
  FullSimplify[
    {
      meanTrace - x0,
      axisTrace - eps x1/4,
      branchTrace - 3 eps x1/4,
      branchTrace - 3 axisTrace
    }
  ];
Print["M14 residuals = ", fmt[m14Residuals]];
If[FullSimplify[m14Residuals] =!= {0, 0, 0, 0},
  Print["FAIL: M14"]; Exit[1]
];

primZ = (Delta q1 - Q d1)/Delta^2;
primN = 2 P (Delta p1 - P d1)/Delta^3;
m15Residual =
  FullSimplify[
    (primN/N0sym + primZ/D0sym /. N0sym -> P^2/Delta^2) -
      (2 p1/P - 2 d1/Delta + (Delta q1 - Q d1)/(D0sym Delta^2))
  ];
Print["M15 residual = ", fmt[m15Residual]];
If[FullSimplify[m15Residual] =!= 0, Print["FAIL: M15"]; Exit[1]];

compatDirectShift =
  Coefficient[Normal[Series[compatDirect, {eps, 0, 1}]], eps, 1];
m16Mutation =
  FullSimplify[
    compatDirectShift -
      (n0/Ptarget - 6 S z2/T - 3 S^2 z4/T^2)
  ];
Print["M16 mutation residual = ", fmt[m16Mutation]];
If[FullSimplify[m16Mutation] === 0,
  Print["FAIL: M16 mutation passed unexpectedly"]; Exit[1]
];

m17Mutation =
  FullSimplify[
    transportShift - (-6 S z2/T - 3 S^2 z4/T^2)
  ];
Print["M17 mutation residual = ", fmt[m17Mutation]];
If[FullSimplify[m17Mutation] === 0,
  Print["FAIL: M17 mutation passed unexpectedly"]; Exit[1]
];

Print["STATUS: PASS"];
Exit[0];
