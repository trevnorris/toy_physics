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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"];

rF1 = 1.77799353547498;
gStar = 0.758035078944663;
piStar = 1.50882951349316;
sigma0Star = 1.80594111095636;

n = 2401;
xGrid = Subdivide[0.0, 1.0, n - 1];
dx = xGrid[[2]] - xGrid[[1]];
weights = ConstantArray[dx, n];
weights[[1]] = dx/2;
weights[[-1]] = dx/2;

kappa = N[Pi/2];
cGrid = Cos[(Pi*xGrid)/2];
kqGrid = Cosh[kappa*(1 - xGrid)]/Cosh[kappa];

normalize[sig_List] := sig/Total[sig*weights];

tsOperator[sig_List] := Module[{cumSig, cumY},
  cumSig = Accumulate[sig*weights];
  cumY = Accumulate[xGrid*sig*weights];
  cumY + xGrid*(Last[cumSig] - cumSig)
];

tqOperator[sig_List] := Module[{aTerm, bTerm},
  aTerm = Accumulate[Sinh[kappa*xGrid]*sig*weights];
  bTerm = Reverse[Accumulate[Reverse[Cosh[kappa*(1 - xGrid)]*sig*weights]]];
  (Cosh[kappa*(1 - xGrid)]*aTerm + Sinh[kappa*xGrid]*bTerm)/(kappa*Cosh[kappa])
];

gMoment[sig_List] := Total[cGrid*sig*weights];
sMoment[sig_List] := Total[kqGrid*sig*weights];
rMoment[sig_List] := ((gMoment[sig] - rF1)^2)/(1 + rF1^2);
phi[sig_List, sigma0_?NumericQ] := sigma0*(tsOperator[sig] - rMoment[sig]*tqOperator[sig]);
nextSigma[sig_List, sigma0_?NumericQ] := Module[{ph, phShift},
  ph = phi[sig, sigma0];
  phShift = ph - Min[ph];
  normalize[Exp[-phShift]]
];

solveFixedPoint[sigma0_?NumericQ, maxIt_: 400, tol_: 1.*^-13] := Module[
  {sig, sigNew, err = Infinity, iter},
  sig = normalize[piStar*Exp[-piStar*xGrid]];
  For[iter = 1, iter <= maxIt, iter++,
    sigNew = nextSigma[sig, sigma0];
    err = Max[Abs[sigNew - sig]];
    sig = sigNew;
    If[err < tol,
      Return[<|"sigma" -> sig, "iterations" -> iter, "error" -> err|>]
    ];
  ];
  <|"sigma" -> sig, "iterations" -> maxIt, "error" -> err|>
];

fp = solveFixedPoint[sigma0Star];
sigFp = fp["sigma"];
gFp = gMoment[sigFp];
sFp = sMoment[sigFp];
rFp = rMoment[sigFp];
piFp = sigma0Star*(1 - rFp*sFp);
deltaRDirect = rFp - 1/4;
deltaG = gFp - gStar;
deltaRPred = -deltaG/Sqrt[1 + rF1^2] + deltaG^2/(1 + rF1^2);

Print["iterations = ", fp["iterations"]];
Print["max residual = ", fp["error"]];
Print["g_fp = ", N[gFp, 20]];
Print["S_fp = ", N[sFp, 20]];
Print["R_fp = ", N[rFp, 20]];
Print["Pi_fp = ", N[piFp, 20]];
Print["delta R from direct solve = ", N[deltaRDirect, 30]];
Print["delta R from exact transport law = ", N[deltaRPred, 30]];

expectApprox["g_fp numeric check", gFp, 0.693352419668063, 10^-12];
expectApprox["S_fp numeric check", sFp, 0.6216013167514007, 10^-12];
expectApprox["R_fp numeric check", rFp, 0.2827139049082381, 10^-12];
expectApprox["Pi_fp numeric check", piFp, 1.4885734438300713, 10^-12];
expectApprox["transport law consistency", deltaRPred, deltaRDirect, 10^-10];

Print[""];
Print["Conclusion:"];
Print["  At frozen canonical traction, the co-evolving fixed point stays close in Pi"];
Print["  but drifts to R > 1/4, so exact compensation is lost."];

Exit[0];
