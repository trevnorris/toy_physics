ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 154 — MICROSCOPIC GROUPED OUTLET OBSTRUCTIONS"];

Clear[p0, dK, dM, dB0, dB2, dZ0, dZ2, dN0];
$Assumptions = Element[{p0, dK, dM, dB0, dB2, dZ0, dZ2, dN0}, Reals] && p0 != 0;

kExact = Expand[(-dM - dB2 - dZ2) + (dK - dB0 - dZ0)/9];
kSplit = Expand[(dK/9 - dM) - (dB2 + dB0/9) - (dZ2 + dZ0/9)];
expectZero["K_A split", kExact - kSplit];

gExact = Expand[dN0 - p0*(dK - dB0 - dZ0)];
gSplit = Expand[-p0*dK + p0*dB0 + (dN0 + p0*dZ0)];
expectZero["G_A split", gExact - gSplit];

Clear[c, w, dc, dw];
$Assumptions = Element[{c, w, dc, dw}, Reals] && c != 0 && w != 0;

b0 = c^2/w^2;
b2 = c^2/w^4;

dB0Exact = D[b0, c]*dc + D[b0, w]*dw;
dB2Exact = D[b2, c]*dc + D[b2, w]*dw;

dB0Formula = 2*c*dc/w^2 - 2*c^2*dw/w^3;
dB2Formula = 2*c*dc/w^4 - 4*c^2*dw/w^5;
expectZero["delta B0 formula", dB0Exact - dB0Formula];
expectZero["delta B2 formula", dB2Exact - dB2Formula];

bCombExact = Expand[dB2Exact + dB0Exact/9];
bCombFormula = Expand[
  2*c*(1/w^4 + 1/(9*w^2))*dc
  - 2*c^2*(2/w^5 + 1/(9*w^3))*dw
];
expectZero["BdG obstruction bundle", bCombExact - bCombFormula];

Clear[u, wVar, r, gu, gw, dU, dW, dR, dgu, dgw];
$Assumptions = Element[{u, wVar, r, gu, gw, dU, dW, dR, dgu, dgw}, Reals] &&
  u != 0 && wVar != 0 && r != 0 && gu != 0 && gw != 0;

deltaExpr = u*wVar - r^2;
sExpr = u + wVar;
qExpr = gu^2*wVar + 2*gu*gw*r + gw^2*u;
gExpr = gu^2 + gw^2;
pExpr = u*gw + r*gu;

dDeltaExpr = D[deltaExpr, u]*dU + D[deltaExpr, wVar]*dW + D[deltaExpr, r]*dR;
dSExpr = D[sExpr, u]*dU + D[sExpr, wVar]*dW;
dQExpr =
  D[qExpr, u]*dU + D[qExpr, wVar]*dW + D[qExpr, r]*dR +
  D[qExpr, gu]*dgu + D[qExpr, gw]*dgw;
dGExpr = D[gExpr, gu]*dgu + D[gExpr, gw]*dgw;
dPExpr = D[pExpr, u]*dU + D[pExpr, r]*dR + D[pExpr, gu]*dgu + D[pExpr, gw]*dgw;

expectZero["delta Delta formula", dDeltaExpr - (wVar*dU + u*dW - 2*r*dR)];
expectZero["delta S formula", dSExpr - (dU + dW)];
expectZero["delta G formula", dGExpr - (2*gu*dgu + 2*gw*dgw)];
expectZero["delta P formula", dPExpr - (gw*dU + gu*dR + r*dgu + u*dgw)];
expectZero[
  "delta Q formula",
  dQExpr - (
    gw^2*dU + gu^2*dW + 2*gu*gw*dR +
    2*(gu*wVar + gw*r)*dgu +
    2*(gw*u + gu*r)*dgw
  )
];

Clear[delta, s, q, gSym, p, dDelta, dS, dQ, dG, dP];
$Assumptions = Element[{delta, s, q, gSym, p, dDelta, dS, dQ, dG, dP, p0}, Reals] &&
  delta != 0 && p != 0 && p0 != 0;

z0 = q/delta;
z2 = (q*s - gSym*delta)/delta^2;
n0Expr = p^2/delta^2;

dZ0Exact = D[z0, q]*dQ + D[z0, delta]*dDelta;
dZ2Exact = D[z2, q]*dQ + D[z2, s]*dS + D[z2, gSym]*dG + D[z2, delta]*dDelta;
dN0Exact = D[n0Expr, p]*dP + D[n0Expr, delta]*dDelta;

expectZero["delta Z0 formula", dZ0Exact - (delta*dQ - q*dDelta)/delta^2];
expectZero[
  "delta Z2 formula",
  dZ2Exact - (
    s*dQ/delta^2 + q*dS/delta^2 - dG/delta +
    (gSym/delta^2 - 2*q*s/delta^3)*dDelta
  )
];
expectZero["delta N0 formula", dN0Exact - (2*p*dP/delta^2 - 2*p^2*dDelta/delta^3)];

zCombExact = Expand[dZ2Exact + dZ0Exact/9];
zCombFormula = Expand[
  (s/delta^2 + 1/(9*delta))*dQ +
  (q/delta^2)*dS -
  dG/delta +
  (gSym/delta^2 - q/(9*delta^2) - 2*q*s/delta^3)*dDelta
];
expectZero["Z obstruction bundle", zCombExact - zCombFormula];

nCombExact = Expand[dN0Exact + p0*dZ0Exact];
nCombFormula = Expand[
  (p0/delta)*dQ + (2*p/delta^2)*dP -
  (p0*q/delta^2 + 2*p^2/delta^3)*dDelta
];
expectZero["N obstruction bundle", nCombExact - nCombFormula];

Clear[eps, k1, m1, b01, b21, z01, z21, n01];
$Assumptions = Element[{eps, k1, m1, b01, b21, z01, z21, n01, p0}, Reals] && p0 != 0;

Do[
  lamVal = lam;
  kMicro = eps*lamVal*(k1/9 - m1 - b21 - b01/9 - z21 - z01/9);
  gMicro = eps*lamVal*(n01 - p0*k1 + p0*b01 + p0*z01);
  kRebuilt = eps*lamVal*((-m1 - b21 - z21) + (k1 - b01 - z01)/9);
  gRebuilt = eps*lamVal*(n01 - p0*(k1 - b01 - z01));
  expectZero[
    "weak-axisymmetric K obstruction lambda=" <> ToString[InputForm[lamVal]],
    kMicro - kRebuilt
  ];
  expectZero[
    "weak-axisymmetric G obstruction lambda=" <> ToString[InputForm[lamVal]],
    gMicro - gRebuilt
  ],
  {lam, {1, 1/2, -1}}
];

Print[""];
Print["Carry-forward formulas:"];
Print["  K_A = (1/9 delta K_A - delta M_A) - (delta B_A2 + 1/9 delta B_A0) - (delta Z_A2 + 1/9 delta Z_A0)"];
Print["  G_A = -P0 delta K_A + P0 delta B_A0 + (delta N_A0 + P0 delta Z_A0)"];
Print["  The weak grouped branch collapses to the scalar pair (K_1, G_1)."];

Print[""];
Print["Stage 154 Mathematica audit passed."];

Exit[0];
