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
  Q, S2, Hport, Delta, D0, mu1, ell,
  Qx, Sx, Hx, Deltax, Px, Gx
];

$Assumptions =
  Element[
    {Q, S2, Hport, Delta, D0, mu1, ell,
      Qx, Sx, Hx, Deltax, Px, Gx},
    Reals
  ] && Delta != 0 && D0 != 0 && mu1 != 0;

(* Primitive identities from the Stage 012 one-port physical source forms. *)
Z0[Q_, Delta_] := Q/Delta;
Z2[Q_, S2_, Hport_, Delta_] := (Q*S2 - Hport*Delta)/Delta^2;
Z4[Q_, S2_, Hport_, Delta_] := (Q*(S2^2 - Delta) - S2*Hport*Delta)/Delta^3;

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
K1 = reduce[-(z2d + z0d/9)];
He = reduce[-z4d + (2/3)*z2d - z0d/27];

Do[
  expectZero["M1 K1 independence from " <> ToString[sym], D[K1, sym]];
  expectZero["M1 He independence from " <> ToString[sym], D[He, sym]],
  {sym, {Px, Gx}}
];

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
