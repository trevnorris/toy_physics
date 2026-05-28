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

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 156 — RENORMALIZED CANONICAL BRANCH"];

rF1 = 1.77799353547498;
gStar = 0.758035078944663;
piStar = 1.50882951349316;
sigma0Star = 1.80594111095636;
tHatStar = 0.901484054174204;

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

solveFixedPoint[sigma0_?NumericQ, maxIt_: 500, tol_: 1.*^-13] := Module[
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

fixedPointAt[sigma0_?NumericQ] := fixedPointAt[sigma0] = solveFixedPoint[sigma0];
gFp[sigma0_?NumericQ] := gFp[sigma0] = gMoment[fixedPointAt[sigma0]["sigma"]];

bracketRoot[lo_?NumericQ, hi_?NumericQ, target_?NumericQ] := Module[{flo, fhi},
  flo = gFp[lo] - target;
  fhi = gFp[hi] - target;
  If[flo == 0, Return[{lo, lo, flo, flo}]];
  If[fhi == 0, Return[{hi, hi, fhi, fhi}]];
  If[flo*fhi > 0, fail["root bracket", {lo, hi, flo, fhi}]];
  {lo, hi, flo, fhi}
];

bisect[lo_?NumericQ, hi_?NumericQ, target_?NumericQ, iters_: 55] := Module[
  {left, right, fLeft, fRight, mid, fMid, data},
  data = bracketRoot[lo, hi, target];
  {left, right, fLeft, fRight} = data;
  Do[
    mid = (left + right)/2;
    fMid = gFp[mid] - target;
    If[fLeft*fMid <= 0,
      right = mid;
      fRight = fMid;,
      left = mid;
      fLeft = fMid;
    ],
    {iters}
  ];
  (left + right)/2
];

scanPoints = N[Range[3.0, 6.0, 0.5], 20];
scanValues = gFp /@ scanPoints;
scanDiffs = Differences[scanValues];

Print["g_fp scan on [3,6] = ", N[Transpose[{scanPoints, scanValues}], 16]];
expectTrue["monotone increase on scan grid", Min[scanDiffs] > 0];

sigma0Can = bisect[3.0, 6.0, gStar, 55];
fpCan = fixedPointAt[sigma0Can];
sigCan = fpCan["sigma"];
gCan = gMoment[sigCan];
sCan = sMoment[sigCan];
rCan = rMoment[sigCan];
piCan = sigma0Can*(1 - rCan*sCan);
tHatCan = Sqrt[9*sigma0Can/20];

Print["Sigma0_can = ", N[sigma0Can, 20]];
Print["g_can      = ", N[gCan, 20]];
Print["S_can      = ", N[sCan, 20]];
Print["R_can      = ", N[rCan, 20]];
Print["Pi_can     = ", N[piCan, 20]];
Print["T_hat_can  = ", N[tHatCan, 20]];
Print[""];
Print["Relative shifts from original canonical point:"];
Print["Sigma0 ratio - 1 = ", N[sigma0Can/sigma0Star - 1, 20]];
Print["Pi ratio - 1     = ", N[piCan/piStar - 1, 20]];
Print["T_hat ratio - 1  = ", N[tHatCan/tHatStar - 1, 20]];

expectApprox["Sigma0_can numeric check", sigma0Can, 4.651033550168867, 10^-11];
expectApprox["g_can numeric check", gCan, 0.758035078944663, 10^-12];
expectApprox["S_can numeric check", sCan, 0.6703621156734617, 10^-12];
expectApprox["R_can numeric check", rCan, 0.25, 10^-10];
expectApprox["Pi_can numeric check", piCan, 3.871564377479002, 10^-11];
expectApprox["T_hat_can numeric check", tHatCan, 1.4467083664567613, 10^-11];

Print[""];
Print["Conclusion:"];
Print["  Full co-evolution preserves the lower compensated branch,"];
Print["  but only after a unique upward traction renormalization."];

Exit[0];
