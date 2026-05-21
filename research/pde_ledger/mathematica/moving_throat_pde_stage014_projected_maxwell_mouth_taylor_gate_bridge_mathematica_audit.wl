(* moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];
reduce[expr_] := FullSimplify[Together[expr], Assumptions -> $Assumptions];

expectZero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[res =!= 0,
    Print["FAIL: ", label]; Exit[1],
    Print["PASS: ", label]
  ];
];

expectNonzero[label_String, expr_] := Module[{res},
  res = reduce[expr];
  Print[label, " residual = ", fmt[res]];
  If[res === 0,
    Print["FAIL: ", label]; Exit[1],
    Print["PASS: ", label]
  ];
];

expectEqual[label_String, got_, expected_] := Module[{res},
  res = FullSimplify[got == expected, Assumptions -> $Assumptions];
  Print[label, " residual = ", fmt[res]];
  If[TrueQ[res],
    Print["PASS: ", label],
    Print["FAIL: ", label]; Print["  got      = ", fmt[got]];
    Print["  expected = ", fmt[expected]]; Exit[1]
  ];
];

Print["STAGE 014 PROJECTED MAXWELL MOUTH-TAYLOR GATE BRIDGE MATHEMATICA AUDIT"];

Clear[
  Q, S2, Hport, Delta, P, Gw, D0, D2, D4, Nbase, mu1, ell,
  Qx, Sx, Hx, Deltax, Px, Gx
];

$Assumptions =
  Element[
    {Q, S2, Hport, Delta, P, Gw, D0, D2, D4, Nbase, mu1, ell,
      Qx, Sx, Hx, Deltax, Px, Gx},
    Reals
  ] && Delta != 0 && P != 0 && D0 != 0 && mu1 != 0;

(* Primitive identities from the Stage 012 one-port physical source forms. *)
Z0[Q_, Delta_] := Q/Delta;
Z2[Q_, S2_, Hport_, Delta_] := (Q*S2 - Hport*Delta)/Delta^2;
Z4[Q_, S2_, Hport_, Delta_] := (Q*(S2^2 - Delta) - S2*Hport*Delta)/Delta^3;
N0[P_, Delta_] := P^2/Delta^2;
N2[P_, Gw_, S2_, Delta_] := 2*P*(P*S2 - Delta*Gw)/Delta^3;
N4[P_, Gw_, S2_, Delta_] :=
  (Delta^2*Gw^2 - 2*Delta*P^2 - 4*Delta*P*S2*Gw + 3*P^2*S2^2)/Delta^4;

z0d =
  reduce[(D[Z0[Q + mu1*Qx*ell, Delta + mu1*Deltax*ell], ell] /. ell -> 0)/mu1];
z2d =
  reduce[
    (D[
        Z2[
          Q + mu1*Qx*ell,
          S2 + mu1*Sx*ell,
          Hport + mu1*Hx*ell,
          Delta + mu1*Deltax*ell
        ],
        ell
      ] /. ell -> 0)/mu1
  ];
z4d =
  reduce[
    (D[
        Z4[
          Q + mu1*Qx*ell,
          S2 + mu1*Sx*ell,
          Hport + mu1*Hx*ell,
          Delta + mu1*Deltax*ell
        ],
        ell
      ] /. ell -> 0)/mu1
  ];
n0d =
  reduce[(D[N0[P + mu1*Px*ell, Delta + mu1*Deltax*ell], ell] /. ell -> 0)/mu1];
n2d =
  reduce[
    (D[
        N2[
          P + mu1*Px*ell,
          Gw + mu1*Gx*ell,
          S2 + mu1*Sx*ell,
          Delta + mu1*Deltax*ell
        ],
        ell
      ] /. ell -> 0)/mu1
  ];
n4d =
  reduce[
    (D[
        N4[
          P + mu1*Px*ell,
          Gw + mu1*Gx*ell,
          S2 + mu1*Sx*ell,
          Delta + mu1*Deltax*ell
        ],
        ell
      ] /. ell -> 0)/mu1
  ];

XiLoad = reduce[n0d/N0[P, Delta] + z0d/D0];
K1 = reduce[-(z2d + z0d/9)];
He = reduce[-z4d + (2/3)*z2d - z0d/27];

deltaP2 =
  reduce[
    (D0^2*n2d - 2*D0*D2*n0d + 2*D0*Nbase*z2d - 2*D2*Nbase*z0d)/D0^3
  ];
deltaP4 =
  reduce[
    (D0^3*n4d - 2*D0^2*D2*n2d - 2*D0^2*D4*n0d +
        2*D0^2*Nbase*z4d + 3*D0*D2^2*n0d - 2*D0*D2*Nbase*z2d -
        2*D0*D4*Nbase*z0d + 2*D2^2*Nbase*z0d)/D0^4
  ];

Do[
  expectZero["M1 XiLoad independence from " <> ToString[sym], D[XiLoad, sym]],
  {sym, {Sx, Hx, Gx}}
];
Do[
  expectZero["M1 K1 independence from " <> ToString[sym], D[K1, sym]];
  expectZero["M1 He independence from " <> ToString[sym], D[He, sym]],
  {sym, {Px, Gx}}
];

expectZero["M2 dXiLoad/dPx", D[XiLoad, Px] - 2/P];
expectZero["M3 d(deltaP2)/dGx", D[deltaP2, Gx] - (-2*P/(D0*Delta^2))];
expectNonzero["M4 deltaP4 Gx dependence", D[deltaP4, Gx]];

qdSolve =
  FullSimplify[
    Solve[{K1 == 0, He == 0} /. {Sx -> 0, Hx -> 0}, {Qx, Deltax}],
    Assumptions -> $Assumptions
  ];
expectEqual["M5 source/denominator solve", qdSolve, {{Qx -> 0, Deltax -> 0}}];

shSolve =
  FullSimplify[
    Solve[{K1 == 0, He == 0} /. {Qx -> 0, Deltax -> 0}, {Sx, Hx}],
    Assumptions -> $Assumptions
  ];
expectEqual["M6 spectral solve", shSolve, {{Sx -> 0, Hx -> 0}}];

qdExprs = {K1, He} /. {Sx -> 0, Hx -> 0};
shExprs = {K1, He} /. {Qx -> 0, Deltax -> 0};
jacobianQD = Table[D[qdExprs[[i]], var], {i, 1, 2}, {var, {Qx, Deltax}}];
jacobianSH = Table[D[shExprs[[i]], var], {i, 1, 2}, {var, {Sx, Hx}}];

expectNonzero["M7 source/denominator Jacobian determinant", Det[jacobianQD]];
expectNonzero["M8 spectral Jacobian determinant", Det[jacobianSH]];

compSurface =
  FullSimplify[Solve[{K1 == 0, He == 0}, {Hx, Sx}], Assumptions -> $Assumptions];
hxDen = Factor[Denominator[Together[Hx /. compSurface[[1]]]]];
sxDen = Factor[Denominator[Together[Sx /. compSurface[[1]]]]];
expectZero[
  "M9 Hx compensation denominator",
  hxDen - 9*Delta^2*(Delta*Hport - Q*S2)
];
expectZero[
  "M9 Sx compensation denominator",
  sxDen - 9*Delta*(Delta*Hport - Q*S2)
];
expectNonzero[
  "M10 sign-flip denominator mutation",
  Denominator[Together[Hx /. compSurface[[1]]]] - 9*Delta^2*(Delta*Hport + Q*S2)
];

Print["STATUS: PASS"];
Exit[0];
