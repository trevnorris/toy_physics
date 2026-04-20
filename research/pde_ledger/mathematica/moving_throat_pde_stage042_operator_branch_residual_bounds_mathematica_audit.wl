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

expectPositive[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[FullSimplify[res > 0, Assumptions -> $Assumptions]], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 50] - N[target, 50]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 042 — OPERATOR BRANCH RESIDUAL BOUNDS"];

Clear[capitalXi, delta0, deltaInf, kappa, y, zetaReq, peReq, pe];
$Assumptions =
  Element[{capitalXi, delta0, deltaInf, kappa, y, zetaReq, peReq, pe}, Reals] &&
  capitalXi > 0 && delta0 > 0 && deltaInf > 0 && kappa > 0 && y > 0 &&
  zetaReq > 0 && peReq > 0 && pe > 0;

omegaFn = Symbol["Omega"];
aK = FullSimplify[(kappa + Pi^2/4)/(kappa + y^2), Assumptions -> $Assumptions];

zetaLo = FullSimplify[aK*omegaFn[capitalXi*delta0]^2, Assumptions -> $Assumptions];
zetaHi = FullSimplify[aK*omegaFn[capitalXi*deltaInf]^2, Assumptions -> $Assumptions];
rLo = FullSimplify[zetaReq - zetaHi, Assumptions -> $Assumptions];
rHi = FullSimplify[zetaReq - zetaLo, Assumptions -> $Assumptions];
xiFail = FullSimplify[peReq/deltaInf, Assumptions -> $Assumptions];
xiSuff = FullSimplify[peReq/delta0, Assumptions -> $Assumptions];

Print["zeta_- = ", fmt[zetaLo]];
Print["zeta_+ = ", fmt[zetaHi]];
Print["R_- = ", fmt[rLo]];
Print["R_+ = ", fmt[rHi]];
Print["Xi_fail = ", fmt[xiFail]];
Print["Xi_suff = ", fmt[xiSuff]];
expectZero["residual bracket center identity", rHi - rLo - (zetaHi - zetaLo)];
Clear[deltaGap];
$Assumptions = $Assumptions && deltaGap > 0;
deltaInfOrdered = FullSimplify[delta0 + deltaGap, Assumptions -> $Assumptions];
xiFailOrdered = FullSimplify[peReq/deltaInfOrdered, Assumptions -> $Assumptions];
xiSuffOrdered = FullSimplify[peReq/delta0, Assumptions -> $Assumptions];
expectPositive["Xi_suff - Xi_fail on the ordered branch", xiSuffOrdered - xiFailOrdered];
expectPositive["Pe_req - Xi_fail Delta_0", FullSimplify[peReq - xiFailOrdered*delta0, Assumptions -> $Assumptions]];
expectPositive["Xi_suff Delta_inf - Pe_req", FullSimplify[xiSuffOrdered*deltaInfOrdered - peReq, Assumptions -> $Assumptions]];

omegaPe = FullSimplify[
  Pi*pe*(2*pe*Exp[pe] + Pi)/((4*pe^2 + Pi^2)*(Exp[pe] - 1)),
  Assumptions -> pe > 0
];
omegaSqSeries = FullSimplify[Normal[Series[omegaPe^2, {pe, 0, 1}]], Assumptions -> pe > 0];
Clear[xiNum];
kappaProbe = 2;
yProbe = 1;
delta0Probe = 3/5;
deltaGapProbe = 2/5;
deltaInfProbe = delta0Probe + deltaGapProbe;
peReqProbe = 7/10;
aKProbe = N[(kappaProbe + Pi^2/4)/(kappaProbe + yProbe^2), 80];
zetaReqProbe = N[aKProbe*(omegaPe /. pe -> peReqProbe)^2, 80];
xiFailExpected = N[peReqProbe/deltaInfProbe, 80];
xiSuffExpected = N[peReqProbe/delta0Probe, 80];
xiFailSolved = xiNum /. FindRoot[
  N[aKProbe*(omegaPe /. pe -> xiNum*deltaInfProbe)^2 - zetaReqProbe, 80],
  {xiNum, xiFailExpected},
  WorkingPrecision -> 70
];
xiSuffSolved = xiNum /. FindRoot[
  N[aKProbe*(omegaPe /. pe -> xiNum*delta0Probe)^2 - zetaReqProbe, 80],
  {xiNum, xiSuffExpected},
  WorkingPrecision -> 70
];
expectPositive["nsolve-style Xi_fail root stayed positive", xiFailSolved];
expectPositive["nsolve-style Xi_suff root stayed positive", xiSuffSolved];
expectApprox["FindRoot Xi_fail saturation solver", xiFailSolved, xiFailExpected, 10^-40];
expectApprox["FindRoot Xi_suff saturation solver", xiSuffSolved, xiSuffExpected, 10^-40];
weakZeta = Expand[aK*(1 + ((4 - Pi)/Pi)*capitalXi*delta0)];

Print["Omega_Pe^2 small-Pe series = ", fmt[omegaSqSeries]];
Print["weak-coupling zeta_phys = ", fmt[weakZeta]];
expectZero["Omega^2 linear coefficient", Coefficient[omegaSqSeries, pe, 1] - (4 - Pi)/Pi];

Print[""];
Print["Stage 042 Mathematica audit passed."];

Exit[0];
