(***
  4D -> preliminary 2PN P0/P2 mouth-port operator rebuild
  --------------------------------------------------------
  Purpose:
    Continue the constructive 2PN program past the channel-kernel solve by
    showing that the solved added conservative 2PN comparable-mass cross block
    closes exactly in a small body-local mouth-port basis.

  Main result proved here:
    - the frozen 1PN wake is already a body-local dipole kernel;
    - after removing the universal 1PN leg dressing, the genuinely new quartic
      2PN tensor sector is exactly a positive overlap of a monopole-like P0 port
      plus the full real P2 multiplet {m=0, ±1, ±2};
    - the scalar-potential dressing excites only the scalar P0/P20 ports;
    - the static cross term is a pure monopole-potential overlap.

  This is the first explicit algebraic bridge from the solved 2PN comparable-
  mass cross block to the inner-throat / mouth-port language.
***)

ClearAll["Global`*"];
Print["--- 4D preliminary 2PN P0/P2 mouth-port operator rebuild ---"];

passCount = 0;
failCount = 0;
section[name_String] := Print["\n=== ", name, " ==="];
pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);
info[msg_String] := Print["INFO: ", msg];

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet @ FullSimplify[res == 0, assum]]
];
checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  If[zeroScalarQ[expr, assum], pass[name], fail[name, res]]
];
checkEqScalar[name_String, lhs_, rhs_, assum_: True] := checkScalar[name, lhs - rhs, assum];
checkMatrix[name_String, lhs_, rhs_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[lhs - rhs, assum];
  If[And @@ Flatten[Map[zeroScalarQ[#, assum] &, res, {2}]], pass[name], fail[name, res]]
];

positiveAssumptions = Gconst > 0 && cLight > 0 && mA > 0 && mB > 0 && rAB > 0;

(* ---------------------------------------------------------------------- *)
section["MASTER / component kinematics and solved targets"];
(* ---------------------------------------------------------------------- *)

vA2 = Expand[vAx^2 + vAy^2 + vAz^2];
vB2 = Expand[vBx^2 + vBy^2 + vBz^2];
vAB = Expand[vAx vBx + vAy vBy + vAz vBz];

dA = vAz;
dB = vBz;

uAdotuB = Expand[vAx vBx + vAy vBy];

UA = Expand[Gconst mB/rAB];
UB = Expand[Gconst mA/rAB];

L1WakeTarget = Expand[-(7/2) vAB - (1/2) dA dB];

QuarticTarget = Expand[
  -(7/4) vAB (vA2 + vB2)
  - (1/4) dA dB (vA2 + vB2)
  + (11/8) vA2 vB2
  + (1/4) vAB^2
  - (5/8) (vA2 dB^2 + vB2 dA^2)
  + (3/2) vAB dA dB
  + (3/8) dA^2 dB^2
];

QuadraticTarget = Expand[
  (11/8) (mA vA2 + mB vB2)
  - (15/4) (mA + mB) vAB
  + (15/8) (mA dA^2 + mB dB^2)
];

AddedCrossTarget = Expand[
  (Gconst mA mB)/(cLight^4 rAB) QuarticTarget
  + (Gconst^2 mA mB)/(cLight^4 rAB^2) QuadraticTarget
  + 5 Gconst^3 mA^2 mB^2/(4 cLight^4 rAB^3)
];

(* ---------------------------------------------------------------------- *)
section["1PN / frozen wake as a body-local dipole kernel"];
(* ---------------------------------------------------------------------- *)

Pi1pmA = {
  Sqrt[7/2] vAx,
  Sqrt[7/2] vAy
};
Pi1pmB = {
  Sqrt[7/2] vBx,
  Sqrt[7/2] vBy
};
Pi10A = 2 dA;
Pi10B = 2 dB;

L1WakeDipole = Expand[-(Pi1pmA.Pi1pmB + Pi10A Pi10B)];

checkEqScalar[
  "The frozen 1PN wake is exactly a body-local dipole kernel",
  L1WakeDipole,
  L1WakeTarget,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["2PN / universal leg dressing plus quartic residual"];
(* ---------------------------------------------------------------------- *)

LegDressed1PN = Expand[(1/2) (vA2 + vB2) L1WakeDipole];
QuarticResidual = Expand[QuarticTarget - LegDressed1PN];

info["QuarticResidual is the genuinely new quartic tensor sector once the universal 1PN leg dressing is removed."];

(* ---------------------------------------------------------------------- *)
section["2PN / exact P0 plus full P2 mouth-port factorization"];
(* ---------------------------------------------------------------------- *)

Pi0A = Expand[Sqrt[5] vA2/2];
Pi0B = Expand[Sqrt[5] vB2/2];

P20A = Expand[3 dA^2 - vA2];
P20B = Expand[3 dB^2 - vB2];
Pi20A = Expand[P20A/2];
Pi20B = Expand[P20B/2];

Pi21A = {
  Sqrt[2] dA vAx,
  Sqrt[2] dA vAy
};
Pi21B = {
  Sqrt[2] dB vBx,
  Sqrt[2] dB vBy
};

Pi22A = {
  (vAx^2 - vAy^2)/(2 Sqrt[2]),
  (2 vAx vAy)/(2 Sqrt[2])
};
Pi22B = {
  (vBx^2 - vBy^2)/(2 Sqrt[2]),
  (2 vBx vBy)/(2 Sqrt[2])
};

QuarticPorts = Expand[
  Pi0A Pi0B
  + Pi20A Pi20B
  + Pi21A.Pi21B
  + Pi22A.Pi22B
];

checkEqScalar[
  "The new quartic residual is exactly the positive overlap of P0 plus the full real P2 multiplet",
  QuarticPorts,
  QuarticResidual,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["DIAGONALIZATION / old scalar TL sector becomes P0 and P20"];
(* ---------------------------------------------------------------------- *)

KTL = {{3/2, 3/4}, {3/4, 9/4}};
Mport = {{1, 1}, {-1, 2}};  (* {P0, P20} = Mport . {T,L} *)
KPort = FullSimplify[Transpose[Inverse[Mport]].KTL.Inverse[Mport]];

checkMatrix[
  "The TL scalar kernel diagonalizes to diag(5/4, 1/4) in the {P0,P20} basis",
  KPort,
  {{5/4, 0}, {0, 1/4}},
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["2PN / scalar-potential dressing excites only P0 and P20"];
(* ---------------------------------------------------------------------- *)

QuadraticPorts = Expand[
  -(15/4) (mA + mB) vAB
  + mA (2 vA2 + (5/8) P20A)
  + mB (2 vB2 + (5/8) P20B)
];

checkEqScalar[
  "The quadratic G^2/r^2 velocity block equals an isotropic dipole term plus scalar P0/P20 driving",
  QuadraticPorts,
  QuadraticTarget,
  positiveAssumptions
];

QuadraticPortsNormalizedBlock = Expand[
  (Gconst mA mB)/(cLight^4 rAB)
  (
    -(15/4) (UA + UB) vAB
    + UA ((4/Sqrt[5]) Pi0B + (5/4) Pi20B)
    + UB ((4/Sqrt[5]) Pi0A + (5/4) Pi20A)
  )
];

QuadraticTargetBlock = Expand[(Gconst^2 mA mB)/(cLight^4 rAB^2) QuadraticTarget];

checkEqScalar[
  "The normalized scalar-potential block is exact in the P0/P20 port basis",
  QuadraticPortsNormalizedBlock,
  QuadraticTargetBlock,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["STATIC / cross term is a pure monopole-potential overlap"];
(* ---------------------------------------------------------------------- *)

StaticCrossPortsBlock = Expand[(Gconst mA mB)/(cLight^4 rAB) (5/4) UA UB];
StaticCrossTargetBlock = Expand[5 Gconst^3 mA^2 mB^2/(4 cLight^4 rAB^3)];

checkEqScalar[
  "The static cross term is exactly +(5/4) U_A U_B in the overall cross normalization",
  StaticCrossPortsBlock,
  StaticCrossTargetBlock,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["ASSEMBLY / full added 2PN cross block in dipole plus P0/P2 language"];
(* ---------------------------------------------------------------------- *)

AddedCrossPorts = Expand[
  (Gconst mA mB)/(cLight^4 rAB)
  (
    LegDressed1PN
    + QuarticPorts
    - (15/4) (UA + UB) vAB
    + UA ((4/Sqrt[5]) Pi0B + (5/4) Pi20B)
    + UB ((4/Sqrt[5]) Pi0A + (5/4) Pi20A)
    + (5/4) UA UB
  )
];

checkEqScalar[
  "The full added 2PN cross block matches exactly in dipole plus P0/P2 mouth-port language",
  AddedCrossPorts,
  AddedCrossTarget,
  positiveAssumptions
];

(* Export useful symbols for later import / assembly. *)
Port2PNDipoleWake = L1WakeDipole;
Port2PNLegDressedWake = LegDressed1PN;
Port2PNP0A = Pi0A;
Port2PNP0B = Pi0B;
Port2PNP20A = Pi20A;
Port2PNP20B = Pi20B;
Port2PNP21A = Pi21A;
Port2PNP21B = Pi21B;
Port2PNP22A = Pi22A;
Port2PNP22B = Pi22B;
Port2PNQuarticResidual = QuarticPorts;
Port2PNQuadraticPotentialBlock = QuadraticPortsNormalizedBlock;
Port2PNStaticCrossBlock = StaticCrossPortsBlock;
Port2PNAddedCrossBlock = AddedCrossPorts;

(* ---------------------------------------------------------------------- *)
section["SUMMARY"];
(* ---------------------------------------------------------------------- *)

info["The frozen 1PN wake is a body-local dipole kernel with transverse m=±1 and longitudinal m=0 pieces."];
info["After removing the universal leg dressing, the genuinely new quartic 2PN tensor sector is exactly a positive overlap of one P0 scalar port plus the full real P2 multiplet {m=0,±1,±2}."];
info["The old TL scalar kernel diagonalizes exactly to diag(5/4,1/4) in the {P0,P20} basis."];
info["The scalar-potential dressing excites only the scalar ports P0 and P20, while the static cross term is +(5/4) U_A U_B."];
info["Therefore the solved added conservative 2PN cross block already has a minimal constructive mouth-port interpretation: carried-forward dipole wake plus a new P0⊕P2 response layer."];

Print["\n--- FINAL COUNTS ---"];
Print["Passes: ", passCount];
Print["Fails:  ", failCount];
