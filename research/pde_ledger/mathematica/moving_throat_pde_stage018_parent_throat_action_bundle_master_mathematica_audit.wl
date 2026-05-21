(* Unit 018 Mathematica audit.
   Mathematica counterpart of the SymPy audit for the parent throat action bundle master.
   Claims:
   M1 one-pole numerator identity.
   M2 one-pole KSigma closure.
   M3 normalization closure.
   M3-mut normalization mutation.
   M4 compatibility cross-closure.
   M5 even-gate determinant.
   M5-mut determinant mutation.
   M6 wall-stiffness and wall-inertia closed-form slopes.
   M7 residual amplitude from closed-form slopes.
   M7-mut residual sign mutation.
   M8 Gaussian wall inertia and stiffness integrals.
*)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Print["STAGE 018 PARENT THROAT ACTION BUNDLE MASTER MATHEMATICA AUDIT"];

Clear[KSigma, MSigma, B0, B2, B4, Z0, Z2, Z4, N0, Ptarget];
$Assumptions =
  Element[{KSigma, MSigma, B0, B2, B4, Z0, Z2, Z4, N0, Ptarget}, Reals] &&
    B4 + Z4 != 0 && Ptarget != 0 && KSigma - B0 - Z0 != 0;

D0 = KSigma - B0 - Z0;
D2 = -(MSigma + B2 + Z2);
D4 = -(B4 + Z4);
poleSeries = 1/(D0 + D2*x^2 + D4*x^4);
seriesToQuartic = Normal[Series[poleSeries, {x, 0, 4}]];
u2 = Coefficient[seriesToQuartic*D0, x, 2];
u4 = Coefficient[seriesToQuartic*D0, x, 4];
onePoleNumerator =
  (D0*(B4 + Z4) - 3*(MSigma + B2 + Z2)^2)/D0^2;
Print["M1 one-pole numerator identity residual = ",
  fmt[FullSimplify[(u4 - 4*u2^2) - onePoleNumerator]]];
If[FullSimplify[(u4 - 4*u2^2) - onePoleNumerator] =!= 0,
  (Print["FAIL: M1 one-pole numerator identity"]; Exit[1])];

KFromOnePole = B0 + Z0 + 3*(MSigma + B2 + Z2)^2/(B4 + Z4);
Print["M2 one-pole KSigma closure residual = ",
  fmt[FullSimplify[(u4 - 4*u2^2) /. KSigma -> KFromOnePole]]];
If[FullSimplify[(u4 - 4*u2^2) /. KSigma -> KFromOnePole] =!= 0,
  (Print["FAIL: M2 one-pole KSigma closure"]; Exit[1])];

KFromNorm = B0 + Z0 + N0/Ptarget;
Print["M3 normalization closure residual = ",
  fmt[FullSimplify[((N0/D0) /. KSigma -> KFromNorm) - Ptarget]]];
If[FullSimplify[((N0/D0) /. KSigma -> KFromNorm) - Ptarget] =!= 0,
  (Print["FAIL: M3 normalization closure"]; Exit[1])];

Print["M3-mut normalization mutation residual = ",
  fmt[FullSimplify[((N0/D0) /. KSigma -> KFromNorm) - 2*Ptarget]]];
If[FullSimplify[((N0/D0) /. KSigma -> KFromNorm) - 2*Ptarget] === 0,
  (Print["FAIL: M3-mut normalization mutation unexpectedly vanished"]; Exit[1])];

compatibility = N0/Ptarget - 3*(MSigma + B2 + Z2)^2/(B4 + Z4);
N0FromCompatibility = N0 /. First[Solve[compatibility == 0, N0]];
N0FromEquality = N0 /. First[Solve[KFromNorm - KFromOnePole == 0, N0]];
Print["M4 compatibility cross-closure residual = ",
  fmt[FullSimplify[N0FromCompatibility - N0FromEquality]]];
If[FullSimplify[N0FromCompatibility - N0FromEquality] =!= 0,
  (Print["FAIL: M4 compatibility cross-closure"]; Exit[1])];

Clear[dKSigma, dMSigma, B01, B21, B41, Z01, Z21, Z41, N01];
$Assumptions =
  Element[
    {KSigma, B0, Z0, N0, N01, dKSigma, dMSigma, B01, B21, B41, Z01, Z21, Z41},
    Reals
  ] && KSigma - B0 - Z0 != 0;

D01 = dKSigma - B01 - Z01;
D21 = -(dMSigma + B21 + Z21);
D41 = -(B41 + Z41);
gateVector = {
  D21 + D01/9,
  D41 - (2/3)*D21 - D01/27
};
K1 = gateVector[[1]];
HEven = gateVector[[2]];
gateJacobian = D[gateVector, {{dKSigma, dMSigma}}];
detGate = Det[gateJacobian];
Print["M5 even-gate determinant residual = ",
  fmt[FullSimplify[detGate - 1/27]]];
If[FullSimplify[detGate - 1/27] =!= 0,
  (Print["FAIL: M5 even-gate determinant"]; Exit[1])];

Print["M5-mut determinant mutation residual = ",
  fmt[FullSimplify[detGate + 1/27]]];
If[FullSimplify[detGate + 1/27] === 0,
  (Print["FAIL: M5-mut determinant mutation unexpectedly vanished"]; Exit[1])];

closedSlopeRules = {
  dKSigma -> B01 + Z01 + 27*(B41 + Z41),
  dMSigma -> -(B21 + Z21) + 3*(B41 + Z41)
};
Print["M6 closed slopes K1 residual = ",
  fmt[FullSimplify[K1 /. closedSlopeRules]]];
If[FullSimplify[K1 /. closedSlopeRules] =!= 0,
  (Print["FAIL: M6 closed slopes K1"]; Exit[1])];

Print["M6 closed slopes H_even residual = ",
  fmt[FullSimplify[HEven /. closedSlopeRules]]];
If[FullSimplify[HEven /. closedSlopeRules] =!= 0,
  (Print["FAIL: M6 closed slopes H_even"]; Exit[1])];

Xi1 = N01/N0 - D01/D0;
expectedXi1 = N01/N0 - 27*(B41 + Z41)/(KSigma - B0 - Z0);
Print["M7 residual amplitude residual = ",
  fmt[FullSimplify[(Xi1 /. closedSlopeRules) - expectedXi1]]];
If[FullSimplify[(Xi1 /. closedSlopeRules) - expectedXi1] =!= 0,
  (Print["FAIL: M7 residual amplitude"]; Exit[1])];

Print["M7-mut residual sign mutation residual = ",
  fmt[FullSimplify[(Xi1 /. closedSlopeRules) -
    (N01/N0 + 27*(B41 + Z41)/(KSigma - B0 - Z0))]]];
If[FullSimplify[(Xi1 /. closedSlopeRules) -
    (N01/N0 + 27*(B41 + Z41)/(KSigma - B0 - Z0))] === 0,
  (Print["FAIL: M7-mut residual sign mutation unexpectedly vanished"]; Exit[1])];

Clear[w];
$Assumptions = Element[w, Reals];
beta = Exp[-w^2/2];
massIntegral = Integrate[beta^2, {w, -Infinity, Infinity}];
stiffnessIntegral = Integrate[D[beta, w]^2 + beta^2, {w, -Infinity, Infinity}];
Print["M8 Gaussian inertia integral residual = ",
  fmt[FullSimplify[massIntegral - Sqrt[Pi]]]];
If[FullSimplify[massIntegral - Sqrt[Pi]] =!= 0,
  (Print["FAIL: M8 Gaussian inertia integral"]; Exit[1])];

Print["M8 Gaussian stiffness integral residual = ",
  fmt[FullSimplify[stiffnessIntegral - 3*Sqrt[Pi]/2]]];
If[FullSimplify[stiffnessIntegral - 3*Sqrt[Pi]/2] =!= 0,
  (Print["FAIL: M8 Gaussian stiffness integral"]; Exit[1])];

Print["STAGE 018 MATHEMATICA AUDIT PASS"];
Exit[0];
