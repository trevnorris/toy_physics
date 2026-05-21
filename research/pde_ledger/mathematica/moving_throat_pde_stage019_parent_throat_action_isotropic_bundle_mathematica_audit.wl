(* Unit 019 Mathematica audit: parent throat action isotropic bundle. *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];
reduce[expr_] := FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
reducePower[expr_] := FullSimplify[PowerExpand[Together[Expand[expr]]], Assumptions -> $Assumptions];

expectZero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res === 0],
    Null,
    Print["FAIL: ", label, " residual = ", fmt[res]];
    Exit[1]
  ];
];

expectZeroPower[label_String, expr_] := Module[{res},
  res = reducePower[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res === 0],
    Null,
    Print["FAIL: ", label, " residual = ", fmt[res]];
    Exit[1]
  ];
];

expectNonzero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res === 0],
    Print["FAIL: ", label, " residual = ", fmt[res]];
    Exit[1]
  ];
];

Print["STAGE 019 PARENT THROAT ACTION ISOTROPIC BUNDLE MATHEMATICA AUDIT"];

Clear[
  KSigma, MSigma, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4,
  mhat0, G, cs, a, c, eps, x, w
];

D0 = KSigma - B0 - Z0;
D2 = -(MSigma + B2 + Z2);
D4 = -(B4 + Z4);

P0target = 54 G cs^5/(5 a^5 c^5 mhat0^2);
Mgap = Sqrt[D0 (B4 + Z4)/3];
Mplus = Mgap - (B2 + Z2);
Mminus = -Mgap - (B2 + Z2);
N2closed = 2 N0 (B2 + MSigma + Z2)/(B0 - KSigma + Z0);
N4closed =
  N0 (
    2 B0 B4 + 2 B0 Z4 + B2^2 + 2 B2 MSigma + 2 B2 Z2 -
    2 B4 KSigma + 2 B4 Z0 - 2 KSigma Z4 + MSigma^2 +
    2 MSigma Z2 + 2 Z0 Z4 + Z2^2
  )/(B0 - KSigma + Z0)^2;

$Assumptions =
  Element[
    {KSigma, MSigma, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4,
      mhat0, G, cs, a, c, eps},
    Reals
  ] && D0 != 0 && B4 + Z4 != 0 && N0 != 0 && G != 0 &&
    cs != 0 && a != 0 && c != 0 && mhat0 != 0 && P0target != 0 &&
    eps != 0;

den = D0 + D2*x + D4*x^2;
poleSeries = Normal[Series[1/den, {x, 0, 2}]];
normalizedPoleSeries = FullSimplify[D0*poleSeries, Assumptions -> $Assumptions];
u2 = FullSimplify[Coefficient[normalizedPoleSeries, x, 1], Assumptions -> $Assumptions];
u4 = FullSimplify[Coefficient[normalizedPoleSeries, x, 2], Assumptions -> $Assumptions];

bundleSeries =
  Normal[Series[D0 (N0 + N2*x + N4*x^2)/den^2, {x, 0, 2}]];
P0 = N0/D0;
P0series = FullSimplify[Coefficient[bundleSeries, x, 0], Assumptions -> $Assumptions];
P2 = FullSimplify[Coefficient[bundleSeries, x, 1], Assumptions -> $Assumptions];
P4 = FullSimplify[Coefficient[bundleSeries, x, 2], Assumptions -> $Assumptions];

expectZero["P0 series normalization", P0series - P0];

expectZero[
  "M1 one-pole numerator",
  (u4 - 4 u2^2) -
    (D0 (B4 + Z4) - 3 (MSigma + B2 + Z2)^2)/D0^2
];
Print["M1 OK"];

KOnePoleSolve = KSigma /. First[Solve[u4 - 4 u2^2 == 0, KSigma]];
expectZero[
  "M2 one-pole KSigma",
  KOnePoleSolve -
    (B0 + Z0 + 3 (MSigma + B2 + Z2)^2/(B4 + Z4))
];
Print["M2 OK"];

KNormalizationSolve = KSigma /. First[Solve[P0 == P0target, KSigma]];
expectZero[
  "M3 normalization KSigma",
  KNormalizationSolve - (B0 + Z0 + N0/P0target)
];
Print["M3 OK"];

N2solve = N2 /. First[Solve[P2 == 0, N2]];
expectZero[
  "M4 N2 constant-prefactor closure",
  N2solve - N2closed
];
Print["M4 OK"];

N4solve = N4 /. First[Solve[(P4 /. N2 -> N2solve) == 0, N4]];
expectZero[
  "M5 N4 constant-prefactor closure",
  N4solve - N4closed
];
Print["M5 OK"];

P2zeroEq = Expand[D0^2 P2];
P4zeroEq = Expand[D0^3 P4];
constantPrefactorJacobian = {
  {D[P2zeroEq, N2], D[P2zeroEq, N4]},
  {D[P4zeroEq, N2], D[P4zeroEq, N4]}
};
expectZero[
  "M6 constant-prefactor Jacobian determinant",
  Det[constantPrefactorJacobian] - D0^3
];
Print["M6 OK"];

expectZero[
  "M7 P2 factorization",
  P2 - (N2 - N2closed)/D0
];
Print["M7 OK"];

expectZero[
  "M8 P4 factorization on N2 closure",
  (P4 /. N2 -> N2closed) - (N4 - N4closed)/D0
];
Print["M8 OK"];

expectNonzero[
  "M9 mutated N2 closure guard",
  P2 - (N2 - (N2closed + eps))/D0
];
expectNonzero[
  "M9 mutated N4 closure guard",
  (P4 /. N2 -> N2closed) - (N4 - (N4closed + eps))/D0
];
Print["M9 OK"];

onePoleNumerator = D0 (B4 + Z4) - 3 (MSigma + B2 + Z2)^2;
expectZeroPower[
  "M10 M-root factorization",
  onePoleNumerator + 3 (MSigma - Mplus) (MSigma - Mminus)
];
expectZeroPower[
  "M10 M-root Vieta sum",
  Mplus + Mminus + 2 (B2 + Z2)
];
expectZeroPower[
  "M10 M-root Vieta product",
  Mplus Mminus - ((B2 + Z2)^2 - D0 (B4 + Z4)/3)
];
Print["M10 OK"];

expectZeroPower[
  "M11 u2 on positive M root",
  (u2 /. MSigma -> Mplus) - Sqrt[D0 (B4 + Z4)/3]/D0
];
expectZeroPower[
  "M11 u2 on negative M root",
  (u2 /. MSigma -> Mminus) + Sqrt[D0 (B4 + Z4)/3]/D0
];
Print["M11 OK"];

Clear[beta, muEta, Tw, Tomega, Keta];
beta[w_] := Exp[-w^2/2];
muEta = 1;
Tw = 1;
Tomega = 1/6;
Keta = 0;
MSigmaExample =
  Integrate[muEta beta[w]^2, {w, -Infinity, Infinity},
    Assumptions -> Element[w, Reals]];
KSigmaExample =
  Integrate[
    Tw D[beta[w], w]^2 + (Keta + 6 Tomega) beta[w]^2,
    {w, -Infinity, Infinity},
    Assumptions -> Element[w, Reals]
  ];
expectZero["M12 Gaussian wall inertia integral", MSigmaExample - Sqrt[Pi]];
expectZero["M12 Gaussian wall stiffness integral", KSigmaExample - 3 Sqrt[Pi]/2];
Print["M12 OK"];

Print["STATUS: PASS"];
Exit[0];
