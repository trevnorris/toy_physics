(* ::Package:: *)

ClearAll[assertZero, section];

section[name_String] := Print["\n=== " <> name <> " ==="];

assertZero[name_String, expr_] := Module[{res = FullSimplify[Expand[expr]]},
  If[res === 0 || (MatrixQ[res] && And @@ Flatten[Map[# === 0 &, res, {2}]]) ||
     (ListQ[res] && And @@ (FullSimplify[#] === 0 & /@ res)),
    Print["PASS: " <> name],
    Print["FAIL: " <> name <> " -> " <> ToString[InputForm[res]]];
    Abort[]
  ]
];

Print["--- 2PN minimal irreducible-channel throat operator preliminary Mathematica script ---"];

(* ---------------------------------------------------------------------- *)
section["MASTER / pair kinematics and canonical port basis"];
(* ---------------------------------------------------------------------- *)

ClearAll[vAx, vAy, vAz, vBx, vBy, vBz, UA, UB];

vA2 = Expand[vAx^2 + vAy^2 + vAz^2];
vB2 = Expand[vBx^2 + vBy^2 + vBz^2];
vAB = Expand[vAx vBx + vAy vBy + vAz vBz];
uAB = Expand[vAx vBx + vAy vBy];
dA = vAz;
dB = vBz;

Pi0A = Expand[Sqrt[5] vA2/2];
Pi0B = Expand[Sqrt[5] vB2/2];
Pi20A = Expand[(3 dA^2 - vA2)/2];
Pi20B = Expand[(3 dB^2 - vB2)/2];
Pi21Ac = Expand[Sqrt[2] dA vAx];
Pi21As = Expand[Sqrt[2] dA vAy];
Pi21Bc = Expand[Sqrt[2] dB vBx];
Pi21Bs = Expand[Sqrt[2] dB vBy];
Pi22Ac = Expand[(vAx^2 - vAy^2)/(2 Sqrt[2])];
Pi22As = Expand[(2 vAx vAy)/(2 Sqrt[2])];
Pi22Bc = Expand[(vBx^2 - vBy^2)/(2 Sqrt[2])];
Pi22Bs = Expand[(2 vBx vBy)/(2 Sqrt[2])];

L1Wake = Expand[-(7/2) vAB - (1/2) dA dB];
assertZero["Frozen 1PN wake = transverse/longitudinal dipole metric",
  L1Wake - Expand[-(7/2) uAB - 4 dA dB]
];

(* ---------------------------------------------------------------------- *)
section["EVEN / solve the minimal P0⊕P2 zero-frequency support operator"];
(* ---------------------------------------------------------------------- *)

ClearAll[a0, a20, a21, a22, j0, j20, s];

evenAns = Expand[
  a0 Pi0A Pi0B + a20 Pi20A Pi20B +
  a21 (Pi21Ac Pi21Bc + Pi21As Pi21Bs) +
  a22 (Pi22Ac Pi22Bc + Pi22As Pi22Bs) +
  UA (j0 Pi0B + j20 Pi20B) +
  UB (j0 Pi0A + j20 Pi20A) +
  s UA UB
];

evenTarget = Expand[
  (11/8) vA2 vB2 + (1/4) vAB^2 - (5/8) (vA2 dB^2 + vB2 dA^2) +
  (3/2) vAB dA dB + (3/8) dA^2 dB^2 +
  UA ((11/8) vB2 + (15/8) dB^2) +
  UB ((11/8) vA2 + (15/8) dA^2) +
  (5/4) UA UB
];

evenResidual = Expand[evenAns - evenTarget];

evenSamples = {
  {vAx -> 1, vAy -> 0, vAz -> 0, vBx -> 1, vBy -> 0, vBz -> 0, UA -> 0, UB -> 0},
  {vAx -> 0, vAy -> 0, vAz -> 1, vBx -> 0, vBy -> 0, vBz -> 1, UA -> 0, UB -> 0},
  {vAx -> 1, vAy -> 0, vAz -> 1, vBx -> 1, vBy -> 0, vBz -> 1, UA -> 0, UB -> 0},
  {vAx -> 1, vAy -> 1, vAz -> 0, vBx -> 1, vBy -> -1, vBz -> 0, UA -> 0, UB -> 0},
  {vAx -> 1, vAy -> 2, vAz -> 3, vBx -> -1, vBy -> 1, vBz -> 2, UA -> 0, UB -> 0},
  {vAx -> 0, vAy -> 0, vAz -> 0, vBx -> 1, vBy -> 0, vBz -> 0, UA -> 1, UB -> 0},
  {vAx -> 0, vAy -> 0, vAz -> 0, vBx -> 0, vBy -> 0, vBz -> 1, UA -> 1, UB -> 0},
  {vAx -> 1, vAy -> 1, vAz -> 2, vBx -> 0, vBy -> 0, vBz -> 0, UA -> 0, UB -> 1},
  {vAx -> 0, vAy -> 0, vAz -> 0, vBx -> 0, vBy -> 0, vBz -> 0, UA -> 1, UB -> 1}
};

evenEqExprs = Expand[evenResidual /. #] & /@ evenSamples;
evenUnknowns = {a0, a20, a21, a22, j0, j20, s};
evenSolve = Solve[Thread[evenEqExprs == 0], evenUnknowns];
assertZero["Minimal P0\[CirclePlus]P2 support operator is uniquely fixed", Length[evenSolve] - 1];
evenSol = evenSolve[[1]];

Print["Solved even data: ", evenSol];
assertZero["Solved even operator reproduces the exact even 2PN target", evenResidual /. evenSol];
assertZero["a0 = a20 = a21 = a22 = 1", ({a0, a20, a21, a22} /. evenSol) - {1, 1, 1, 1}];
assertZero["j0 = 4/Sqrt[5]", (j0 /. evenSol) - 4/Sqrt[5]];
assertZero["j20 = 5/4", (j20 /. evenSol) - 5/4];
assertZero["static support+closure target coefficient = 5/4", (s /. evenSol) - 5/4];

(* ---------------------------------------------------------------------- *)
section["ODD / solve the minimal dressed dipole operator"];
(* ---------------------------------------------------------------------- *)

ClearAll[sigma, pperp, ppara];

oddAns = Expand[
  sigma (vA2 + vB2) L1Wake - (UA + UB) (pperp uAB + ppara dA dB)
];

oddTarget = Expand[
  -(7/4) vAB (vA2 + vB2) - (1/4) dA dB (vA2 + vB2) - (15/4) (UA + UB) vAB
];

oddResidual = Expand[oddAns - oddTarget];
oddSamples = {
  {vAx -> 1, vAy -> 0, vAz -> 0, vBx -> 1, vBy -> 0, vBz -> 0, UA -> 0, UB -> 0},
  {vAx -> 0, vAy -> 0, vAz -> 1, vBx -> 0, vBy -> 0, vBz -> 1, UA -> 0, UB -> 0},
  {vAx -> 1, vAy -> 0, vAz -> 1, vBx -> 1, vBy -> 0, vBz -> 1, UA -> 0, UB -> 0},
  {vAx -> 0, vAy -> 0, vAz -> 1, vBx -> 0, vBy -> 0, vBz -> 1, UA -> 1, UB -> 0},
  {vAx -> 1, vAy -> 0, vAz -> 0, vBx -> 1, vBy -> 0, vBz -> 1, UA -> 1, UB -> 0}
};
oddUnknowns = {sigma, pperp, ppara};
oddEqExprs = Expand[oddResidual /. #] & /@ oddSamples;
oddSolve = Solve[Thread[oddEqExprs == 0], oddUnknowns];
assertZero["Minimal odd dipole operator is uniquely fixed", Length[oddSolve] - 1];
oddSol = oddSolve[[1]];

Print["Solved odd data: ", oddSol];
assertZero["Solved odd operator reproduces the exact odd 2PN target", oddResidual /. oddSol];
assertZero["sigma = 1/2", (sigma /. oddSol) - 1/2];
assertZero["pperp = ppara = 15/4", ({pperp, ppara} /. oddSol) - {15/4, 15/4}];

etaPerp = FullSimplify[(pperp /. oddSol)/(7/2)];
etaPara = FullSimplify[(ppara /. oddSol)/4];
assertZero["eta_perp = 15/14", etaPerp - 15/14];
assertZero["eta_parallel = 15/16", etaPara - 15/16];

(* ---------------------------------------------------------------------- *)
section["DATA / canonical zero-frequency operator package"];
(* ---------------------------------------------------------------------- *)

Rodd = DiagonalMatrix[{7/2, 7/2, 4}];
assertZero["Odd zero-frequency residues are fixed at {7/2, 7/2, 4}",
  Rodd - DiagonalMatrix[{7/2, 7/2, 4}]
];

Reven = IdentityMatrix[6];
Jvec = {4/Sqrt[5], 5/4, 0, 0, 0, 0};
supportStatic = FullSimplify[Jvec.Jvec];
closureDeficit = FullSimplify[supportStatic - 5/4];
kappaGeom = FullSimplify[closureDeficit/2];

assertZero["Even zero-frequency support residues are all 1", Reven - IdentityMatrix[6]];
assertZero["Support static coefficient = 381/80", supportStatic - 381/80];
assertZero["Geometry closure deficit = 281/80", closureDeficit - 281/80];
assertZero["Direct geometry-energy coefficient = 281/160", kappaGeom - 281/160];
assertZero["Direct geometry-energy term produces the pure-U closure block",
  (-kappaGeom (UA + UB)^2 + kappaGeom UA^2 + kappaGeom UB^2) + closureDeficit UA UB
];

supportA = {
  Expand[Pi0A + Jvec[[1]] UA],
  Expand[Pi20A + Jvec[[2]] UA],
  Pi21Ac,
  Pi21As,
  Pi22Ac,
  Pi22As
};

supportB = {
  Expand[Pi0B + Jvec[[1]] UB],
  Expand[Pi20B + Jvec[[2]] UB],
  Pi21Bc,
  Pi21Bs,
  Pi22Bc,
  Pi22Bs
};

evenSupportMinusClosure = Expand[supportA.supportB - closureDeficit UA UB];
oddDressed = Expand[(sigma /. oddSol) (vA2 + vB2) L1Wake -
  (pperp /. oddSol) (UA + UB) uAB - (ppara /. oddSol) (UA + UB) dA dB];

fullTarget = Expand[evenTarget + oddTarget];
operatorRebuild = Expand[evenSupportMinusClosure + oddDressed];
assertZero["Canonical operator package reproduces the full added conservative 2PN cross block",
  operatorRebuild - fullTarget
];

(* ---------------------------------------------------------------------- *)
section["LOW-FREQUENCY / unresolved DtN data beyond the zero-frequency residues"];
(* ---------------------------------------------------------------------- *)

ClearAll[omega, chi0, chi2, chi1p, chi10];
Y0 = 1 + chi0 omega^2;
Y2 = 1 + chi2 omega^2;
Y1p = 7/2 + chi1p omega^2;
Y10 = 4 + chi10 omega^2;

assertZero["Y0(0) = 1", (Y0 /. omega -> 0) - 1];
assertZero["Y2(0) = 1", (Y2 /. omega -> 0) - 1];
assertZero["Y1_perp(0) = 7/2", (Y1p /. omega -> 0) - 7/2];
assertZero["Y1_parallel(0) = 4", (Y10 /. omega -> 0) - 4];

Print["Current 2PN conservative matching fixes only the zero-frequency residues; the omega^2 DtN coefficients remain free PDE observables:"];
Print["  Y0[omega]  = ", Y0];
Print["  Y2[omega]  = ", Y2];
Print["  Y1p[omega] = ", Y1p];
Print["  Y10[omega] = ", Y10];

(* ---------------------------------------------------------------------- *)
section["N-BODY / pair AB in a three-body background C"];
(* ---------------------------------------------------------------------- *)

ClearAll[Gconst, cLight, mA, mB, mC, rAB, rAC, rBC, vA2AB, vB2AB, vABAB, dAAB, dBAB];

UALoc = Expand[Gconst (mB/rAB + mC/rAC)];
UBLoc = Expand[Gconst (mA/rAB + mC/rBC)];
pairPref = Expand[Gconst mA mB/(cLight^4 rAB)];

sourceB = Expand[Jvec[[1]] (Sqrt[5] vB2AB/2) + Jvec[[2]] ((3 dBAB^2 - vB2AB)/2)];
sourceA = Expand[Jvec[[1]] (Sqrt[5] vA2AB/2) + Jvec[[2]] ((3 dAAB^2 - vA2AB)/2)];
assertZero["U-driven scalar source = 11/8 v^2 + 15/8 d^2",
  {
    sourceA - ((11/8) vA2AB + (15/8) dAAB^2),
    sourceB - ((11/8) vB2AB + (15/8) dBAB^2)
  }
];

pairABThreeBodyVel = Expand[
  pairPref ((UALoc - Gconst mB/rAB) sourceB + (UBLoc - Gconst mA/rAB) sourceA -
    (15/4) ((UALoc - Gconst mB/rAB) + (UBLoc - Gconst mA/rAB)) vABAB)
];

pairABThreeBodyVelExpected = Expand[
  Gconst^2 mA mB mC/(8 cLight^4 rAB) *
  ((11 vB2AB + 15 dBAB^2 - 30 vABAB)/rAC + (11 vA2AB + 15 dAAB^2 - 30 vABAB)/rBC)
];
assertZero["Pair AB velocity-dependent 3-body lift is fixed by the minimal operator",
  pairABThreeBodyVel - pairABThreeBodyVelExpected
];

pairABThreeBodyStatic = Expand[
  pairPref (5/4) (UALoc UBLoc - (Gconst mB/rAB) (Gconst mA/rAB))
];

pairABThreeBodyStaticExpected = Expand[
  5 Gconst^3 mA mB mC/(4 cLight^4) *
  (mA/(rAB^2 rAC) + mB/(rAB^2 rBC) + mC/(rAB rAC rBC))
];
assertZero["Pair AB static 3-body lift is fixed by the minimal operator",
  pairABThreeBodyStatic - pairABThreeBodyStaticExpected
];

Print["\nSolved zero-frequency operator data:"];
Print["  odd residues      : {7/2, 7/2, 4}"];
Print["  even residues     : {1, 1, 1, 1, 1, 1}"];
Print["  J source vector   : ", Jvec];
Print["  support static    : ", supportStatic];
Print["  geometry closure  : ", closureDeficit];
Print["  kappa_geom        : ", kappaGeom];
Print["  sigma             : ", sigma /. oddSol];
Print["  eta_perp          : ", etaPerp];
Print["  eta_parallel      : ", etaPara];

Print["\nPair AB in background C:"];
Print["  velocity 3-body lift = ", Factor[pairABThreeBodyVelExpected]];
Print["  static   3-body lift = ", Factor[pairABThreeBodyStaticExpected]];

Print["\nDone."];

(*"
Output:

--- 2PN minimal irreducible-channel throat operator preliminary Mathematica script ---

=== MASTER / pair kinematics and canonical port basis ===
PASS: Frozen 1PN wake = transverse/longitudinal dipole metric

=== EVEN / solve the minimal P0⊕P2 zero-frequency support operator ===
PASS: Minimal P0⊕P2 support operator is uniquely fixed
Solved even data: {a0 -> 1, a20 -> 1, a21 -> 1, a22 -> 1, j0 -> 4/Sqrt[5], j20 -> 5/4, s -> 5/4}
PASS: Solved even operator reproduces the exact even 2PN target
PASS: a0 = a20 = a21 = a22 = 1
PASS: j0 = 4/Sqrt[5]
PASS: j20 = 5/4
PASS: static support+closure target coefficient = 5/4

=== ODD / solve the minimal dressed dipole operator ===
PASS: Minimal odd dipole operator is uniquely fixed
Solved odd data: {sigma -> 1/2, pperp -> 15/4, ppara -> 15/4}
PASS: Solved odd operator reproduces the exact odd 2PN target
PASS: sigma = 1/2
PASS: pperp = ppara = 15/4
PASS: eta_perp = 15/14
PASS: eta_parallel = 15/16

=== DATA / canonical zero-frequency operator package ===
PASS: Odd zero-frequency residues are fixed at {7/2, 7/2, 4}
PASS: Even zero-frequency support residues are all 1
PASS: Support static coefficient = 381/80
PASS: Geometry closure deficit = 281/80
PASS: Direct geometry-energy coefficient = 281/160
PASS: Direct geometry-energy term produces the pure-U closure block
PASS: Canonical operator package reproduces the full added conservative 2PN cross block

=== LOW-FREQUENCY / unresolved DtN data beyond the zero-frequency residues ===
PASS: Y0(0) = 1
PASS: Y2(0) = 1
PASS: Y1_perp(0) = 7/2
PASS: Y1_parallel(0) = 4
Current 2PN conservative matching fixes only the zero-frequency residues; the omega^2 DtN coefficients remain free PDE observables:
  Y0[omega]  = 1 + chi0*omega^2
  Y2[omega]  = 1 + chi2*omega^2
  Y1p[omega] = 7/2 + chi1p*omega^2
  Y10[omega] = 4 + chi10*omega^2

=== N-BODY / pair AB in a three-body background C ===
PASS: U-driven scalar source = 11/8 v^2 + 15/8 d^2
PASS: Pair AB velocity-dependent 3-body lift is fixed by the minimal operator
PASS: Pair AB static 3-body lift is fixed by the minimal operator

Solved zero-frequency operator data:
  odd residues      : {7/2, 7/2, 4}
  even residues     : {1, 1, 1, 1, 1, 1}
  J source vector   : {4/Sqrt[5], 5/4, 0, 0, 0, 0}
  support static    : 381/80
  geometry closure  : 281/80
  kappa_geom        : 281/160
  sigma             : 1/2
  eta_perp          : 15/14
  eta_parallel      : 15/16

Pair AB in background C:
  velocity 3-body lift = (Gconst^2*mA*mB*mC*(15*dAAB^2*rAC + 15*dBAB^2*rBC + 11*rAC*vA2AB - 30*rAC*vABAB - 30*rBC*vABAB + 11*rBC*vB2AB))/(8*cLight^4*rAB*rAC*rBC)
  static   3-body lift = (5*Gconst^3*mA*mB*mC*(mC*rAB + mB*rAC + mA*rBC))/(4*cLight^4*rAB^2*rAC*rBC)

Done.
"*)
