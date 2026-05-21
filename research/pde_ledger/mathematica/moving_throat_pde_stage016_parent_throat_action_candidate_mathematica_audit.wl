(* moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl *)
(* Header: paths, claim summary *)
(* Verifies M1-M11 for the parent-throat action candidate. *)

ClearAll["Global`*"];
$HistoryLength = 0;
On[Assert];

fmt[expr_] := ToString[InputForm[expr]];
reduce[expr_] := FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];

expectZero[tag_String, expr_] := Module[{ok},
  ok = reduce[expr];
  Print[tag, " residual = ", fmt[ok]];
  If[ok =!= 0,
    Print[tag, " FAIL: ", fmt[ok]];
    Exit[1]
  ];
  Print[tag, " PASS"];
];

expectNonzero[tag_String, expr_] := Module[{muto},
  muto = reduce[expr];
  Print[tag, " residual = ", fmt[muto]];
  If[muto === 0,
    Print[tag, " FAIL: residual collapsed"];
    Exit[1]
  ];
  Print[tag, " PASS"];
];

Print["STAGE 016 PARENT THROAT ACTION CANDIDATE MATHEMATICA AUDIT"];
Print["Script: mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl"];
Print["Claims: exact EL, fluctuation K_eta, boundary discharges, Y20 angular sector, modal EL."];

(* M1: exact nonlinear Euler-Lagrange equation. *)
Clear[t, w, u, v, R, mu, Tw, TO, USig];
$Assumptions = Element[{t, w, u, v}, Reals];
r = R[t, w, u, v];
muField = mu[r, w];
TwField = Tw[r, w];
TOField = TO[r, w];
USigField = USig[r, w];
rt = D[r, t];
rw = D[r, w];
ru = D[r, u];
rv = D[r, v];
lagExact =
  muField*rt^2/2 - TwField*rw^2/2 - TOField*(ru^2 + rv^2)/2 -
    USigField;
elRaw =
  D[lagExact, r] - D[D[lagExact, D[r, t]], t] -
    D[D[lagExact, D[r, w]], w] - D[D[lagExact, D[r, u]], u] -
    D[D[lagExact, D[r, v]], v];
elViaD = -elRaw;
handExact =
  D[muField*rt, t] - D[TwField*rw, w] - D[TOField*ru, u] -
    D[TOField*rv, v] - D[muField, r]*rt^2/2 +
    D[TwField, r]*rw^2/2 + D[TOField, r]*(ru^2 + rv^2)/2 +
    D[USigField, r];
expectZero["[M1]", elViaD - handExact];
expectNonzero["[M1 mutation]", elViaD - (handExact - 2*D[USigField, r])];

(* M2-M5: epsilon expansion and cross-term integration by parts. *)
Clear[
  eps, eta, etaT, etaW, grad2, R0p, mu0, Tw0, TO0, U0, TwR0,
  TwRR0, UR0, URR0, dTw0p, dTwR0p
];
$Assumptions =
  Element[
    {eps, eta, etaT, etaW, grad2, R0p, mu0, Tw0, TO0, U0, TwR0,
      TwRR0, UR0, URR0, dTw0p, dTwR0p},
    Reals
  ];
muExp = mu0;
TwExp = Tw0 + eps*TwR0*eta + eps^2*TwRR0*eta^2/2;
TOExp = TO0;
USigExp = U0 + eps*UR0*eta + eps^2*URR0*eta^2/2;
rtExp = eps*etaT;
rwExp = R0p + eps*etaW;
gradExp = eps^2*grad2;
lagExp =
  muExp*rtExp^2/2 - TwExp*rwExp^2/2 - TOExp*gradExp/2 -
    USigExp;
linDensity = Expand[D[lagExp, eps] /. eps -> 0];
quadRaw = Expand[(D[lagExp, {eps, 2}] /. eps -> 0)/2];

expectZero[
  "[M2]",
  linDensity - (-TwR0*R0p^2*eta/2 - Tw0*R0p*etaW - UR0*eta)
];
expectNonzero[
  "[M2 mutation]",
  linDensity - (-TwR0*R0p^2*eta/2 - Tw0*R0p*etaW + UR0*eta)
];

eBg = dTw0p - TwR0*R0p^2/2 - UR0;
expectZero[
  "[M3]",
  (linDensity + Tw0*R0p*etaW) + dTw0p*eta - eBg*eta
];
expectNonzero[
  "[M3 mutation]",
  (linDensity + Tw0*R0p*etaW) - dTw0p*eta - eBg*eta
];
Print["[M3 boundary] carried term = [-Tw(R0,w) R0' eta]"];

expectZero[
  "[M4]",
  quadRaw -
    (mu0*etaT^2/2 - Tw0*etaW^2/2 - TO0*grad2/2 -
      TwR0*R0p*eta*etaW - TwRR0*R0p^2*eta^2/4 -
      URR0*eta^2/2)
];
expectZero["[M4 cross]", D[D[quadRaw, eta], etaW] + TwR0*R0p];
expectNonzero["[M4 mutation]", D[D[quadRaw, eta], etaW] - TwR0*R0p];

Clear[aCheck, etaCheck, x];
$Assumptions = Element[x, Reals];
aEtaRule =
  -aCheck[x]*etaCheck[x]*D[etaCheck[x], x] -
    (D[-aCheck[x]*etaCheck[x]^2/2, x] +
      D[aCheck[x], x]*etaCheck[x]^2/2);
expectZero["[M5 product-rule]", aEtaRule];

$Assumptions =
  Element[
    {eta, etaT, etaW, grad2, R0p, mu0, Tw0, TO0, TwR0, TwRR0, URR0,
      dTwR0p},
    Reals
  ];
kEta = URR0 - dTwR0p + TwRR0*R0p^2/2;
quadAfterIBP = Expand[quadRaw - (-TwR0*R0p*eta*etaW) + dTwR0p*eta^2/2];
canonicalQuad =
  mu0*etaT^2/2 - Tw0*etaW^2/2 - TO0*grad2/2 - kEta*eta^2/2;
expectZero["[M5]", quadAfterIBP - canonicalQuad];
Print["L2 after IBP = -R0p^2*TwRR0*eta^2/4 - TO0*grad2/2 - Tw0*eta_w^2/2 - URR0*eta^2/2 + d_TwR_R0p*eta^2/2 + eta_t^2*mu0/2"];
Print["K_eta(w)  = R0p^2*TwRR0/2 + URR0 - d_TwR_R0p"];
expectNonzero[
  "[M5 mutation]",
  quadAfterIBP -
    (mu0*etaT^2/2 - Tw0*etaW^2/2 - TO0*grad2/2 -
      (URR0 + dTwR0p + TwRR0*R0p^2/2)*eta^2/2)
];

(* M6: concrete boundary discharges and IBP integrals. *)
Clear[w];
$Assumptions = Element[w, Reals];
limitDischarge[expr_] :=
  Quiet[
    Limit[expr, w -> Infinity] - Limit[expr, w -> -Infinity],
    Limit::alimv
  ];
Bcon = (1 + w^2)*Exp[-w^2];
Acon = (1 + w^2/2)*Exp[-w^2];
etaG = Exp[-w^2/2];
etaL = 1/(1 + w^2);
Blor = Exp[-w^2];
linGaussian = limitDischarge[-Bcon*etaG];
linLorentz = limitDischarge[-Blor*etaL];
quadGaussian = limitDischarge[-Acon*etaG^2/2];
quadLorentz = limitDischarge[-Acon*etaL^2/2];
linIBP =
  Integrate[-Bcon*D[etaG, w], {w, -Infinity, Infinity}] -
    Integrate[D[Bcon, w]*etaG, {w, -Infinity, Infinity}];
quadIBP =
  Integrate[-Acon*etaG*D[etaG, w], {w, -Infinity, Infinity}] -
    Integrate[D[Acon, w]*etaG^2/2, {w, -Infinity, Infinity}];
expectZero[
  "[M6]",
  linGaussian^2 + linLorentz^2 + quadGaussian^2 + quadLorentz^2
];
expectZero["[M6 linear IBP]", linIBP];
expectZero["[M6 quadratic IBP]", quadIBP];
expectNonzero["[M6 mutation]", linGaussian - 1];

(* M7: Lorentzian finite-endpoint discharge. *)
Bend = w*Sqrt[1 + w^2];
finiteDischarge = limitDischarge[-Bend*etaL];
expectZero["[M7]", finiteDischarge + 2];
expectNonzero["[M7 mutation]", finiteDischarge];

(* M8: nonzero boundary probe and nontrivial Lorentzian denominator. *)
atanDischarge = limitDischarge[ArcTan[w]];
lorentzTogether = Together[-Exp[-w^2]/(1 + w^2)];
expectZero["[M8]", atanDischarge - Pi];
Print["[M8 denominator] Together form = ", fmt[lorentzTogether]];
expectZero["[M8 denominator]", Together[(1 + w^2)*lorentzTogether + Exp[-w^2]]];
expectNonzero["[M8 mutation]", atanDischarge];
expectNonzero["[M8 denominator mutation]", Together[lorentzTogether + Exp[-w^2]]];

(* M9-M10: Y20 angular sector. *)
Clear[th, ph];
$Assumptions = Element[{th, ph}, Reals] && 0 < th < Pi;
y20Raw = SphericalHarmonicY[2, 0, th, ph];
y20 = FullSimplify[ExpToTrig[FunctionExpand[y20Raw]], Assumptions -> $Assumptions];
lapY20 =
  FullSimplify[
    D[Sin[th]*D[y20, th], th]/Sin[th] +
      D[y20, {ph, 2}]/Sin[th]^2,
    Assumptions -> $Assumptions
  ];
expectZero["[M9]", lapY20 + 6*y20];
expectNonzero["[M9 mutation]", lapY20 + 5*y20];

$Assumptions = Element[{th, ph}, Reals];
angularNorm =
  FullSimplify[
    Integrate[Sin[th]*y20^2, {ph, 0, 2*Pi}, {th, 0, Pi}]
  ];
angularStiff =
  FullSimplify[
    Integrate[
      Sin[th]*(D[y20, th]^2 + D[y20, ph]^2/Sin[th]^2),
      {ph, 0, 2*Pi}, {th, 0, Pi}
    ]
  ];
expectZero["[M10]", (angularNorm - 1)^2 + (angularStiff - 6)^2];
Print["Y20 angular norm = ", fmt[angularNorm]];
Print["Y20 angular stiffness = ", fmt[angularStiff]];
Print["Y20 angular sector contributes +6 T_Omega"];
expectNonzero["[M10 mutation]", angularStiff - 5];

(* M11: closed modal Euler-Lagrange equation. *)
Clear[tm, wm, q, mumode, Twmode, TOmode, Kmode];
$Assumptions = Element[{tm, wm}, Reals];
qMode = q[tm, wm];
muMode = mumode[wm];
TwMode = Twmode[wm];
TOMode = TOmode[wm];
KMode = Kmode[wm];
modalDensity =
  (muMode*D[qMode, tm]^2 - TwMode*D[qMode, wm]^2 -
    (KMode + 6*TOMode)*qMode^2)/2;
modalEL =
  D[modalDensity, qMode] - D[D[modalDensity, D[qMode, tm]], tm] -
    D[D[modalDensity, D[qMode, wm]], wm];
modalHand =
  D[muMode*D[qMode, tm], tm] - D[TwMode*D[qMode, wm], wm] +
    (KMode + 6*TOMode)*qMode;
expectZero["[M11]", modalEL + modalHand];
expectNonzero[
  "[M11 mutation]",
  modalEL +
    (D[muMode*D[qMode, tm], tm] - D[TwMode*D[qMode, wm], wm] +
      (KMode + 5*TOMode)*qMode)
];

Print["STATUS: PASS"];
Exit[0];
