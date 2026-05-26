(* moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

assertZero[label_String, expr_] := Module[{res},
  res = FullSimplify[Together[expr], Assumptions -> $Assumptions];
  If[TrueQ[res == 0],
    Print["OK ", label, " residual = ", fmt[res]],
    Print["FAIL ", label, ": ", fmt[res]]; Exit[1]
  ];
];

assertNonzero[label_String, expr_] := Module[{res},
  res = FullSimplify[Together[expr], Assumptions -> $Assumptions];
  If[TrueQ[res == 0],
    Print["FAIL ", label, ": unexpectedly zero"]; Exit[1],
    Print["OK ", label, " residual = ", fmt[res]]
  ];
];

Print["STAGE 013 PROJECTED MAXWELL MOUTH-TAYLOR MASTER MATHEMATICA AUDIT"];

Clear[
  u, ell, t, X0, X1, X2, Q, S2, Hport, Delta, P, Gw, q1, s1, h1,
  d1, p1, g1, D0, mu1, Qx, Sx, Hx, Dx, Px, Gx,
  qFun, sFun, hFun, dFun, pFun, gFun
];

$Assumptions =
  Element[
    {u, ell, t, X0, X1, X2, Q, S2, Hport, Delta, P, Gw, q1, s1,
      h1, d1, p1, g1, D0, mu1, Qx, Sx, Hx, Dx, Px, Gx},
    Reals
  ] && Delta != 0 && P != 0 && D0 != 0 && mu1 != 0;

X = X0 + ell u X1 + ell^2 u^2 X2/2;
Xproj = Integrate[Exp[-u] X, {u, 0, Infinity}];
assertZero[
  "M1 one-sided Taylor first moment",
  Normal[Series[Xproj, {ell, 0, 1}]] - (X0 + ell X1)
];

assertZero[
  "M2 W2 second moment",
  Integrate[u^2 Exp[-u], {u, 0, Infinity}] - 2
];
Xproj2 = Integrate[u Exp[-u] X, {u, 0, Infinity}];
assertZero[
  "M2 W2 first moment recovery",
  Normal[Series[Xproj2, {ell, 0, 1}]] - (X0 + 2 ell X1)
];

valueRules = {
  qFun[0] -> Q, sFun[0] -> S2, hFun[0] -> Hport, dFun[0] -> Delta,
  pFun[0] -> P, gFun[0] -> Gw,
  qFun'[0] -> q1, sFun'[0] -> s1, hFun'[0] -> h1, dFun'[0] -> d1,
  pFun'[0] -> p1, gFun'[0] -> g1
};

timeDerivativeAtZero[expr_] :=
  FullSimplify[(D[expr, t] /. t -> 0) /. valueRules, Assumptions -> $Assumptions];

Zsource[x_] := (qFun[t] - hFun[t] x^2)/(dFun[t] - sFun[t] x^2 + x^4);
Zexpansion = Normal[Series[Zsource[ell], {ell, 0, 4}]];
z0 = timeDerivativeAtZero[Coefficient[Zexpansion, ell, 0]];
z2 = timeDerivativeAtZero[Coefficient[Zexpansion, ell, 2]];
z4 = timeDerivativeAtZero[Coefficient[Zexpansion, ell, 4]];

z0Expected = (Delta q1 - Q d1)/Delta^2;
z2Expected =
  (-Delta^2 h1 + Delta (Hport d1 + Q s1 + S2 q1) -
      2 Q S2 d1)/Delta^3;
z4Expected =
  (-Delta^2 Hport s1 - Delta^2 S2 h1 - Delta^2 q1 +
      2 Delta Hport S2 d1 + 2 Delta Q S2 s1 + 2 Delta Q d1 +
      Delta S2^2 q1 - 3 Q S2^2 d1)/Delta^4;

assertZero["M3 z0 chain-rule coefficient", z0 - z0Expected];
assertZero["M3 z2 chain-rule coefficient", z2 - z2Expected];
assertZero["M3 z4 chain-rule coefficient", z4 - z4Expected];

Nsource[x_] := ((pFun[t] - gFun[t] x^2)^2)/(dFun[t] - sFun[t] x^2 + x^4)^2;
Nexpansion = Normal[Series[Nsource[ell], {ell, 0, 4}]];
n0 = timeDerivativeAtZero[Coefficient[Nexpansion, ell, 0]];
n2 = timeDerivativeAtZero[Coefficient[Nexpansion, ell, 2]];
n4 = timeDerivativeAtZero[Coefficient[Nexpansion, ell, 4]];

n0Expected = 2 P (Delta p1 - P d1)/Delta^3;
n2Expected =
  -(2 Delta^2 (Gw p1 + P g1) -
      2 Delta P (2 Gw d1 + P s1 + 2 S2 p1) +
      6 P^2 S2 d1)/Delta^4;
n4Expected =
  2 (Delta^3 Gw g1 - Delta^2 Gw^2 d1 -
      2 Delta^2 Gw P s1 - 2 Delta^2 Gw S2 p1 -
      2 Delta^2 P S2 g1 - 2 Delta^2 P p1 +
      6 Delta Gw P S2 d1 + 3 Delta P^2 S2 s1 +
      3 Delta P^2 d1 + 3 Delta P S2^2 p1 -
      6 P^2 S2^2 d1)/Delta^5;

assertZero["M4 n0 chain-rule coefficient", n0 - n0Expected];
assertZero["M4 n2 chain-rule coefficient", n2 - n2Expected];
assertZero["M4 n4 chain-rule coefficient", n4 - n4Expected];

subsDer = {
  q1 -> mu1 Qx, s1 -> mu1 Sx, h1 -> mu1 Hx, d1 -> mu1 Dx,
  p1 -> mu1 Px, g1 -> mu1 Gx
};

Xi =
  FullSimplify[
    ((2 p1/P - 2 d1/Delta + q1/(D0 Delta) -
        Q d1/(D0 Delta^2)) /. subsDer)/mu1,
    Assumptions -> $Assumptions
  ];
(* Paper round-trip: verify Xi matches Xi_load = n0/N0 + z0/D0 with N0 = P^2/Delta^2. *)
z0Form = (Delta q1 - Q d1)/Delta^2;
n0Form = 2 P (Delta p1 - P d1)/Delta^3;
N0Form = P^2/Delta^2;
XiPaper =
  FullSimplify[
    ((n0Form/N0Form + z0Form/D0) /. subsDer)/mu1,
    Assumptions -> $Assumptions
  ];
assertZero["M5 Xi matches paper closed form n0/N0 + z0/D0", Xi - XiPaper];
assertZero["M5 dXi/dPprime", D[Xi, Px] - 2/P];

Print["STAGE 013 MATHEMATICA AUDIT: PASS"];
Exit[0];
