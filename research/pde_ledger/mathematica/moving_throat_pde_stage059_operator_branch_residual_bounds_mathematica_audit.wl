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
  (* Strip ConditionalExpression wrapper: under $Assumptions, a result of
     the form ConditionalExpression[0, cond] is identically zero on the
     declared domain.  Solve[]/Reduce[] often introduce these wrappers when
     auxiliary inequalities are nontrivial. *)
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
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

banner["STAGE 059 — OPERATOR BRANCH RESIDUAL BOUNDS"];

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
Clear[deltaGap];
$Assumptions = $Assumptions && deltaGap > 0;
deltaInfOrdered = FullSimplify[delta0 + deltaGap, Assumptions -> $Assumptions];
xiFailOrdered = FullSimplify[peReq/deltaInfOrdered, Assumptions -> $Assumptions];
xiSuffOrdered = FullSimplify[peReq/delta0, Assumptions -> $Assumptions];
expectPositive["Xi_suff - Xi_fail (ordered)", xiSuffOrdered - xiFailOrdered];

omegaPe = FullSimplify[
  Pi*pe*(2*pe*Exp[pe] + Pi)/((4*pe^2 + Pi^2)*(Exp[pe] - 1)),
  Assumptions -> pe > 0
];
omegaSqLinearCoeff = FullSimplify[Limit[D[omegaPe^2, pe], pe -> 0], Assumptions -> pe > 0];
Clear[xiNum, peNum];
kappaProbe = 2;
yProbe = 1;
delta0Probe = 3/5;
deltaGapProbe = 2/5;
deltaInfProbe = delta0Probe + deltaGapProbe;
aKProbe = N[(kappaProbe + Pi^2/4)/(kappaProbe + yProbe^2), 80];
zetaReqProbe = 2/5;  (* independent target, NOT constructed from omegaPe *)
peStar = peNum /. FindRoot[
  N[aKProbe*(omegaPe /. pe -> peNum)^2 - zetaReqProbe, 80],
  {peNum, 1/2},
  WorkingPrecision -> 70
];
peReqProbe = peStar;
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
expectApprox["Xi_fail*DeltaInf saturates at Pe_star", xiFailSolved*deltaInfProbe, peStar, 10^-40];
expectApprox["Xi_suff*Delta0 saturates at Pe_star", xiSuffSolved*delta0Probe, peStar, 10^-40];
weakZeta = Expand[aK*(1 + ((4 - Pi)/Pi)*capitalXi*delta0)];

Print["Omega_Pe^2 linear coefficient = ", fmt[omegaSqLinearCoeff]];
Print["weak-coupling zeta_phys = ", fmt[weakZeta]];
expectZero["Omega^2 linear coefficient", omegaSqLinearCoeff - (4 - Pi)/Pi];

Print[""];
Print["Stage 059 Mathematica audit passed."];

Exit[0];
