(* moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];
reduce[expr_] := FullSimplify[Together[expr], Assumptions -> $Assumptions];

expectZero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res == 0],
    Print["PASS: ", label],
    Print["FAIL: ", label]; Exit[1]
  ];
];

expectNonzero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res == 0],
    Print["FAIL: ", label]; Exit[1],
    Print["PASS: ", label]
  ];
];

expectEqual[label_String, got_, expected_] := Module[{ok},
  ok = FullSimplify[got == expected, Assumptions -> $Assumptions];
  Print[label, " residual = ", fmt[ok]];
  If[TrueQ[ok],
    Print["PASS: ", label],
    Print["FAIL: ", label];
    Print["  got      = ", fmt[got]];
    Print["  expected = ", fmt[expected]];
    Exit[1]
  ];
];

Print["STAGE 015 PARENT THROAT ACTION MASTER MATHEMATICA AUDIT"];

Clear[
  w, eps, gateCoeff, mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0,
  R0p, dTwRR0p, eta, etat, etaw, grad2, dK, dM, b01, b21, b41,
  z01, z21, z41, x0, eps1
];

$Assumptions =
  Element[
    {w, eps, gateCoeff, mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0,
      R0p, dTwRR0p, eta, etat, etaw, grad2, dK, dM, b01, b21, b41,
      z01, z21, z41, x0, eps1},
    Reals
  ];

aFun[w_] := A[w];
etaFun[w_] := eta[w];
ibpProductRule =
  -aFun[w]*etaFun[w]*etaFun'[w] -
    (D[-aFun[w]*etaFun[w]^2/2, w] + aFun'[w]*etaFun[w]^2/2);
expectZero["M1 generic IBP product-rule identity", ibpProductRule];

ibpMutated =
  -aFun[w]*etaFun[w]*etaFun'[w] -
    (D[-(-aFun[w]*etaFun[w]^2/2), w] + aFun'[w]*etaFun[w]^2/2);
expectNonzero["M1 mutated IBP boundary sign", ibpMutated];

aConcrete = Exp[-w^2];
etaConcrete = Exp[-w^2/2];
boundaryGaussian =
  Quiet[
    Limit[-aConcrete*etaConcrete^2/2, w -> Infinity] -
      Limit[-aConcrete*etaConcrete^2/2, w -> -Infinity],
    Limit::alimv
  ];
expectZero["M2 Gaussian IBP boundary discharge", boundaryGaussian];

crossGaussian =
  Integrate[-aConcrete*etaConcrete*D[etaConcrete, w], {w, -Infinity, Infinity}];
bulkGaussian =
  Integrate[D[aConcrete, w]*etaConcrete^2/2, {w, -Infinity, Infinity}];
expectZero["M2 Gaussian IBP cross equals bulk", crossGaussian - bulkGaussian];

(* Second concrete probe with asymmetric A profile so cross and bulk are
   individually nonzero. Asserts a real IBP cancellation, not 0 = 0 + 0. *)
aConcreteAsym = w*Exp[-w^2];
etaConcreteAsym = Exp[-w^2/2];
boundaryGaussianAsym =
  Quiet[
    Limit[-aConcreteAsym*etaConcreteAsym^2/2, w -> Infinity] -
      Limit[-aConcreteAsym*etaConcreteAsym^2/2, w -> -Infinity],
    Limit::alimv
  ];
crossGaussianAsym =
  Integrate[-aConcreteAsym*etaConcreteAsym*D[etaConcreteAsym, w], {w, -Infinity, Infinity}];
bulkGaussianAsym =
  Integrate[D[aConcreteAsym, w]*etaConcreteAsym^2/2, {w, -Infinity, Infinity}];
expectNonzero["M2 asymmetric IBP cross nontrivial", crossGaussianAsym];
expectNonzero["M2 asymmetric IBP bulk nontrivial", bulkGaussianAsym];
expectZero["M2 asymmetric IBP boundary discharge", boundaryGaussianAsym];
expectZero["M2 asymmetric IBP cross equals bulk", crossGaussianAsym - bulkGaussianAsym];

(* M3: K_eta via direct Euler-Lagrange linearization on the parent
   Lagrangian density. A temporary profile twR[w] carries the w-dependence of
   Tw_R0, so the product derivative d/dw(Tw_R0 R0') is computed explicitly
   without relying on Dt's Constants bookkeeping. *)

ClearAll[R0, twR, TwSig, USig, KetaFromEL, rSlot, rtSlot, rwSlot, gOSlot, R0pp, t];
TwSig[r_, w_] := Tw0 + (r - R0[w])*twR[w] + (r - R0[w])^2*TwRR0/2;
USig[r_, w_] := U0 + (r - R0[w])*UR0 + (r - R0[w])^2*URR0/2;

LDensity[R_, Rt_, Rw_, gO_, w_] :=
  mu0*Rt^2/2 - TwSig[R, w]*Rw^2/2 - TO0*gO/2 - USig[R, w];

slotDensity = LDensity[rSlot, rtSlot, rwSlot, gOSlot, w];
dLdRSlot = D[slotDensity, rSlot];
dLdRtSlot = D[slotDensity, rtSlot];
dLdRwSlot = D[slotDensity, rwSlot];
fieldRules = {
  rSlot -> R0[w] + eps*eta,
  rtSlot -> eps*etat,
  rwSlot -> R0'[w] + eps*etaw,
  gOSlot -> eps^2*grad2
};

ELLinearized =
  Expand[
    (dLdRSlot /. fieldRules) -
      D[dLdRtSlot /. fieldRules, t] -
      D[dLdRwSlot /. fieldRules, w]
  ] /. {Derivative[1][R0][w] -> R0p, Derivative[2][R0][w] -> R0pp,
        twR[w] -> TwR0};

ELOrderEps = Coefficient[Expand[ELLinearized], eps, 1];
ELMassCoeff = Coefficient[ELOrderEps, eta];
KetaFromEL =
  FullSimplify[
    -ELMassCoeff /. {R0p*Derivative[1][twR][w] -> dTwRR0p - TwR0*R0pp},
    Assumptions -> $Assumptions
  ];
expectZero["M3 K_eta via EL linearization matches IBP form",
  KetaFromEL - (URR0 - dTwRR0p + TwRR0*R0p^2/2)];

(* Sign-mutation guard: a + instead of - on dTwRR0p must fail. *)
expectNonzero["M3 K_eta via EL dTwRR0p sign mutation",
  KetaFromEL - (URR0 + dTwRR0p + TwRR0*R0p^2/2)];

Print["STATUS: PASS"];
Exit[0];
