(* ::Package:: *)

ClearAll[assertZero, section, ballDtN4DSeries];

section[name_String] := Print["\n=== " <> name <> " ==="];

assertZero[name_String, expr_] := Module[{res = FullSimplify[Expand[expr]]},
  If[res === 0 || (MatrixQ[res] && And @@ Flatten[Map[# === 0 &, res, {2}]]) ||
     (ListQ[res] && And @@ (FullSimplify[#] === 0 & /@ res)),
    Print["PASS: " <> name],
    Print["FAIL: " <> name <> " -> " <> ToString[InputForm[res]]];
    Abort[]
  ]
];

ballDtN4DSeries[l_Integer] := Expand @ Normal @ Series[
  -1/a + z D[BesselJ[l + 1, z], z]/(a BesselJ[l + 1, z]),
  {z, 0, 8}
];

Print["--- 2PN inner-throat modal DtN from PDE-side scaffolding preliminary Mathematica script ---"];

(* ---------------------------------------------------------------------- *)
section["PDE unit test / isotropic 4D-ball DtN branch"];
(* ---------------------------------------------------------------------- *)

Lam0 = ballDtN4DSeries[0];
Lam1 = ballDtN4DSeries[1];
Lam2 = ballDtN4DSeries[2];

Lam0Expected = Expand[-z^2/(4 a) - z^4/(96 a) - z^6/(1536 a) - z^8/(23040 a)];
Lam1Expected = Expand[1/a - z^2/(6 a) - z^4/(288 a) - z^6/(8640 a) - 7 z^8/(1658880 a)];
Lam2Expected = Expand[2/a - z^2/(8 a) - z^4/(640 a) - z^6/(30720 a) - 13 z^8/(17203200 a)];

assertZero["4D-ball monopole DtN series", Lam0 - Lam0Expected];
assertZero["4D-ball dipole DtN series", Lam1 - Lam1Expected];
assertZero["4D-ball quadrupole DtN series", Lam2 - Lam2Expected];
assertZero["Laplace limit: Lambda_0(0)=0", Lam0 /. z -> 0];
assertZero["Laplace limit: Lambda_1(0)=1/a", (Lam1 /. z -> 0) - 1/a];
assertZero["Laplace limit: Lambda_2(0)=2/a", (Lam2 /. z -> 0) - 2/a];

Y1Ball = Expand @ Normal @ Series[1/Lam1, {z, 0, 4}];
Y2Ball = Expand @ Normal @ Series[1/Lam2, {z, 0, 4}];
Y1BallExpected = Expand[a + a z^2/6 + a z^4/32];
Y2BallExpected = Expand[a/2 + a z^2/32 + 3 a z^4/1280];
assertZero["4D-ball dipole admittance series", Y1Ball - Y1BallExpected];
assertZero["4D-ball quadrupole admittance series", Y2Ball - Y2BallExpected];

Print["4D-ball low-k series:"];
Print["  Lambda0[z] = ", Lam0];
Print["  Lambda1[z] = ", Lam1];
Print["  Lambda2[z] = ", Lam2];
Print["  Y1Ball[z]  = ", Y1Ball];
Print["  Y2Ball[z]  = ", Y2Ball];
Print["Implication: isotropic spherical support is a valid PDE unit test, but it cannot furnish finite static monopole stiffness or dipole splitting."];

(* ---------------------------------------------------------------------- *)
section["DTN completion / minimal axisymmetric one-pole-per-channel model"];
(* ---------------------------------------------------------------------- *)

ClearAll[omega, Omega1p, Omega10, Omega0, Omega20, Omega21, Omega22, OmegaGeom];
J0 = 4/Sqrt[5];
J20 = 5/4;
DeltaGeom = 281/80;

Y1p = (7/2)/(1 - omega^2/Omega1p^2);
Y10 = 4/(1 - omega^2/Omega10^2);
Y0 = 1/(1 - omega^2/Omega0^2);
Y20 = 1/(1 - omega^2/Omega20^2);
Y21 = 1/(1 - omega^2/Omega21^2);
Y22 = 1/(1 - omega^2/Omega22^2);
YGeom = 1/(1 - omega^2/OmegaGeom^2);

Z1p = FullSimplify[1/Y1p];
Z10 = FullSimplify[1/Y10];
Z0 = FullSimplify[1/Y0];
Z20 = FullSimplify[1/Y20];
Z21 = FullSimplify[1/Y21];
Z22 = FullSimplify[1/Y22];
ZGeom = FullSimplify[1/YGeom];

assertZero["Y1p(0)=7/2", (Y1p /. omega -> 0) - 7/2];
assertZero["Y10(0)=4", (Y10 /. omega -> 0) - 4];
assertZero["Y0(0)=Y20(0)=Y21(0)=Y22(0)=1", {(Y0 /. omega -> 0) - 1, (Y20 /. omega -> 0) - 1, (Y21 /. omega -> 0) - 1, (Y22 /. omega -> 0) - 1}];
assertZero["YGeom(0)=1", (YGeom /. omega -> 0) - 1];
assertZero["Z1p(0)=2/7", (Z1p /. omega -> 0) - 2/7];
assertZero["Z10(0)=1/4", (Z10 /. omega -> 0) - 1/4];
assertZero["Z0(0)=Z20(0)=Z21(0)=Z22(0)=1", {(Z0 /. omega -> 0) - 1, (Z20 /. omega -> 0) - 1, (Z21 /. omega -> 0) - 1, (Z22 /. omega -> 0) - 1}];
assertZero["ZGeom(0)=1", (ZGeom /. omega -> 0) - 1];

Meven = {
  {Y0, 0, 0, 0, 0, 0, J0 Y0},
  {0, Y20, 0, 0, 0, 0, J20 Y20},
  {0, 0, Y21, 0, 0, 0, 0},
  {0, 0, 0, Y21, 0, 0, 0},
  {0, 0, 0, 0, Y22, 0, 0},
  {0, 0, 0, 0, 0, Y22, 0},
  {J0 Y0, J20 Y20, 0, 0, 0, 0, J0^2 Y0 + J20^2 Y20 - DeltaGeom YGeom}
};

Rsupport = {
  {Sqrt[Y0], 0, 0, 0, 0, 0, J0 Sqrt[Y0]},
  {0, Sqrt[Y20], 0, 0, 0, 0, J20 Sqrt[Y20]},
  {0, 0, Sqrt[Y21], 0, 0, 0, 0},
  {0, 0, 0, Sqrt[Y21], 0, 0, 0},
  {0, 0, 0, 0, Sqrt[Y22], 0, 0},
  {0, 0, 0, 0, 0, Sqrt[Y22], 0}
};

Rgeom = {{0, 0, 0, 0, 0, 0, Sqrt[DeltaGeom YGeom]}};
assertZero["Meven[omega] = support PSD block minus one pure-U closure block", Transpose[Rsupport].Rsupport - Transpose[Rgeom].Rgeom - Meven];

MevenStaticExpected = {
  {1, 0, 0, 0, 0, 0, J0},
  {0, 1, 0, 0, 0, 0, J20},
  {0, 0, 1, 0, 0, 0, 0},
  {0, 0, 0, 1, 0, 0, 0},
  {0, 0, 0, 0, 1, 0, 0},
  {0, 0, 0, 0, 0, 1, 0},
  {J0, J20, 0, 0, 0, 0, 5/4}
};
assertZero["Static even matrix reproduces the solved 2PN support/closure data", (Meven /. omega -> 0) - MevenStaticExpected];

Jeff = {J0 Y0, J20 Y20, 0, 0, 0, 0};
Seff = FullSimplify[J0^2 Y0 + J20^2 Y20 - DeltaGeom YGeom];
assertZero["Static U-U coefficient is 5/4", (Seff /. omega -> 0) - 5/4];

Print["Bare channel DtN kernels:"];
Print["  Z1p[omega]   = ", Z1p];
Print["  Z10[omega]   = ", Z10];
Print["  Z0[omega]    = ", Z0];
Print["  Z20[omega]   = ", Z20];
Print["  Z21[omega]   = ", Z21];
Print["  Z22[omega]   = ", Z22];
Print["  ZGeom[omega] = ", ZGeom];
Print["Jeff[omega] = ", Jeff];
Print["Seff[omega] = ", Seff];

(* ---------------------------------------------------------------------- *)
section["STATIC limit / recover the solved full conservative 2PN cross operator"];
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

sigma = 1/2;
etaPerp = 15/14;
etaPara = 15/16;
L1Wake = Expand[-(7/2) vAB - (1/2) dA dB];
LoddAdded = Expand[
  sigma (vA2 + vB2) L1Wake - (UA + UB) (etaPerp (7/2) uAB + etaPara 4 dA dB)
];

LevenDynamic = Expand[
  Y0 (Pi0A + J0 UA) (Pi0B + J0 UB) +
  Y20 (Pi20A + J20 UA) (Pi20B + J20 UB) +
  Y21 (Pi21Ac Pi21Bc + Pi21As Pi21Bs) +
  Y22 (Pi22Ac Pi22Bc + Pi22As Pi22Bs) -
  DeltaGeom YGeom UA UB
];

LfullStaticFromDtN = Expand[-Y1p uAB - Y10 dA dB + LoddAdded + LevenDynamic] /. omega -> 0;
LfullStaticTarget = Expand[L1Wake + LoddAdded + LevenDynamic /. omega -> 0];
assertZero["Static limit of the dynamic model reproduces the solved full cross operator", LfullStaticFromDtN - LfullStaticTarget];

(* ---------------------------------------------------------------------- *)
section["LOW-FREQUENCY / what 2PN fixed and what remains PDE data"];
(* ---------------------------------------------------------------------- *)

Y1pSeries = Expand @ Normal @ Series[Y1p, {omega, 0, 2}];
Y10Series = Expand @ Normal @ Series[Y10, {omega, 0, 2}];
Y0Series = Expand @ Normal @ Series[Y0, {omega, 0, 2}];
Y20Series = Expand @ Normal @ Series[Y20, {omega, 0, 2}];
Y21Series = Expand @ Normal @ Series[Y21, {omega, 0, 2}];
Y22Series = Expand @ Normal @ Series[Y22, {omega, 0, 2}];
YGeomSeries = Expand @ Normal @ Series[YGeom, {omega, 0, 2}];

Z1pSeries = Expand @ Normal @ Series[Z1p, {omega, 0, 2}];
Z10Series = Expand @ Normal @ Series[Z10, {omega, 0, 2}];
Z0Series = Expand @ Normal @ Series[Z0, {omega, 0, 2}];
Z20Series = Expand @ Normal @ Series[Z20, {omega, 0, 2}];
Z21Series = Expand @ Normal @ Series[Z21, {omega, 0, 2}];
Z22Series = Expand @ Normal @ Series[Z22, {omega, 0, 2}];
ZGeomSeries = Expand @ Normal @ Series[ZGeom, {omega, 0, 2}];

JeffSeries = Expand /@ (Normal @ Series[#, {omega, 0, 2}] & /@ Jeff);
SeffSeries = Expand @ Normal @ Series[Seff, {omega, 0, 2}];
SeffExpected = 5/4 + omega^2 (16/(5 Omega0^2) + 25/(16 Omega20^2) - 281/(80 OmegaGeom^2));
assertZero["Seff low-frequency coefficient is fixed channel-by-channel", SeffSeries - SeffExpected];

Print["Low-frequency admittance series:"];
Print["  Y1p[omega]   = ", Y1pSeries];
Print["  Y10[omega]   = ", Y10Series];
Print["  Y0[omega]    = ", Y0Series];
Print["  Y20[omega]   = ", Y20Series];
Print["  Y21[omega]   = ", Y21Series];
Print["  Y22[omega]   = ", Y22Series];
Print["  YGeom[omega] = ", YGeomSeries];
Print["Low-frequency DtN series:"];
Print["  Z1p[omega]   = ", Z1pSeries];
Print["  Z10[omega]   = ", Z10Series];
Print["  Z0[omega]    = ", Z0Series];
Print["  Z20[omega]   = ", Z20Series];
Print["  Z21[omega]   = ", Z21Series];
Print["  Z22[omega]   = ", Z22Series];
Print["  ZGeom[omega] = ", ZGeomSeries];
Print["JeffSeries = ", JeffSeries];
Print["SeffSeries = ", SeffSeries];
Print["Remaining PDE observables: {Omega1p, Omega10, Omega0, Omega20, Omega21, Omega22, OmegaGeom}.\nOptional near-spherical reduction: set Omega20 = Omega21 = Omega22."];

(* ---------------------------------------------------------------------- *)
section["SUMMARY"];
(* ---------------------------------------------------------------------- *)

Print["1. The isotropic 4D-ball DtN branch is a valid PDE unit test but fails two key checks:"];
Print["   - no finite static monopole support, and"];
Print["   - no dipole splitting."];
Print["2. The minimal throat completion is therefore axisymmetric and must include:"];
Print["   - odd dipole channels {1perp, 10},"];
Print["   - even support channels {0, 20, 21, 22},"];
Print["   - one pure-U geometry closure channel."];
Print["3. The solved 2PN derivation fixes all zero-frequency residues and the scalar source vector, but not the seven pole scales. Those are now clean PDE/DtN observables."];

Print["\nDone."];

(*"
Output:

--- 2PN inner-throat modal DtN from PDE-side scaffolding preliminary Mathematica script ---

=== PDE unit test / isotropic 4D-ball DtN branch ===
PASS: 4D-ball monopole DtN series
PASS: 4D-ball dipole DtN series
PASS: 4D-ball quadrupole DtN series
PASS: Laplace limit: Lambda_0(0)=0
PASS: Laplace limit: Lambda_1(0)=1/a
PASS: Laplace limit: Lambda_2(0)=2/a
PASS: 4D-ball dipole admittance series
PASS: 4D-ball quadrupole admittance series
4D-ball low-k series:
  Lambda0[z] = -1/4*z^2/a - z^4/(96*a) - z^6/(1536*a) - z^8/(23040*a)
  Lambda1[z] = a^(-1) - z^2/(6*a) - z^4/(288*a) - z^6/(8640*a) - (7*z^8)/(1658880*a)
  Lambda2[z] = 2/a - z^2/(8*a) - z^4/(640*a) - z^6/(30720*a) - (13*z^8)/(17203200*a)
  Y1Ball[z]  = a + (a*z^2)/6 + (a*z^4)/32
  Y2Ball[z]  = a/2 + (a*z^2)/32 + (3*a*z^4)/1280
Implication: isotropic spherical support is a valid PDE unit test, but it cannot furnish finite static monopole stiffness or dipole splitting.

=== DTN completion / minimal axisymmetric one-pole-per-channel model ===
PASS: Y1p(0)=7/2
PASS: Y10(0)=4
PASS: Y0(0)=Y20(0)=Y21(0)=Y22(0)=1
PASS: YGeom(0)=1
PASS: Z1p(0)=2/7
PASS: Z10(0)=1/4
PASS: Z0(0)=Z20(0)=Z21(0)=Z22(0)=1
PASS: ZGeom(0)=1
PASS: Meven[omega] = support PSD block minus one pure-U closure block
PASS: Static even matrix reproduces the solved 2PN support/closure data
PASS: Static U-U coefficient is 5/4
Bare channel DtN kernels:
  Z1p[omega]   = (2*(1 - omega^2/Omega1p^2))/7
  Z10[omega]   = (1 - omega^2/Omega10^2)/4
  Z0[omega]    = 1 - omega^2/Omega0^2
  Z20[omega]   = 1 - omega^2/Omega20^2
  Z21[omega]   = 1 - omega^2/Omega21^2
  Z22[omega]   = 1 - omega^2/Omega22^2
  ZGeom[omega] = 1 - omega^2/OmegaGeom^2
Jeff[omega] = {4/(Sqrt[5]*(1 - omega^2/Omega0^2)), 5/(4*(1 - omega^2/Omega20^2)), 0, 0, 0, 0}
Seff[omega] = (256/(1 - omega^2/Omega0^2) + 125/(1 - omega^2/Omega20^2) + (281*OmegaGeom^2)/(omega^2 - OmegaGeom^2))/80

=== STATIC limit / recover the solved full conservative 2PN cross operator ===
PASS: Static limit of the dynamic model reproduces the solved full cross operator

=== LOW-FREQUENCY / what 2PN fixed and what remains PDE data ===
PASS: Seff low-frequency coefficient is fixed channel-by-channel
Low-frequency admittance series:
  Y1p[omega]   = 7/2 + (7*omega^2)/(2*Omega1p^2)
  Y10[omega]   = 4 + (4*omega^2)/Omega10^2
  Y0[omega]    = 1 + omega^2/Omega0^2
  Y20[omega]   = 1 + omega^2/Omega20^2
  Y21[omega]   = 1 + omega^2/Omega21^2
  Y22[omega]   = 1 + omega^2/Omega22^2
  YGeom[omega] = 1 + omega^2/OmegaGeom^2
Low-frequency DtN series:
  Z1p[omega]   = 2/7 - (2*omega^2)/(7*Omega1p^2)
  Z10[omega]   = 1/4 - omega^2/(4*Omega10^2)
  Z0[omega]    = 1 - omega^2/Omega0^2
  Z20[omega]   = 1 - omega^2/Omega20^2
  Z21[omega]   = 1 - omega^2/Omega21^2
  Z22[omega]   = 1 - omega^2/Omega22^2
  ZGeom[omega] = 1 - omega^2/OmegaGeom^2
JeffSeries = {4/Sqrt[5] + (4*omega^2)/(Sqrt[5]*Omega0^2), 5/4 + (5*omega^2)/(4*Omega20^2), 0, 0, 0, 0}
SeffSeries = 5/4 + (16*omega^2)/(5*Omega0^2) + (25*omega^2)/(16*Omega20^2) - (281*omega^2)/(80*OmegaGeom^2)
Remaining PDE observables: {Omega1p, Omega10, Omega0, Omega20, Omega21, Omega22, OmegaGeom}.
Optional near-spherical reduction: set Omega20 = Omega21 = Omega22.

=== SUMMARY ===
1. The isotropic 4D-ball DtN branch is a valid PDE unit test but fails two key checks:
   - no finite static monopole support, and
   - no dipole splitting.
2. The minimal throat completion is therefore axisymmetric and must include:
   - odd dipole channels {1perp, 10},
   - even support channels {0, 20, 21, 22},
   - one pure-U geometry closure channel.
3. The solved 2PN derivation fixes all zero-frequency residues and the scalar source vector, but not the seven pole scales. Those are now clean PDE/DtN observables.

Done.
"*)
