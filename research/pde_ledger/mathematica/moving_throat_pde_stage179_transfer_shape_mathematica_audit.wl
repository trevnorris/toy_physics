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
  res = FullSimplify[Together[Cancel[expr]], Assumptions -> $Assumptions];
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

logSlope[expr_, pairs_List] := Module[{raw},
  raw = Total[(#[[2]]*#[[1]]*D[expr, #[[1]]]/expr) & /@ pairs];
  FullSimplify[Together[Cancel[raw]], Assumptions -> $Assumptions]
];

banner["STAGE 179 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM"];

Clear[k, omegaU, omegaW, gWall, gMixed, coupling];
$Assumptions = Element[{k, omegaU, omegaW, gWall, gMixed, coupling}, Reals] &&
  k > 0 && omegaU > 0 && omegaW > 0 && gWall > 0 && gMixed > 0;

portNumerator = omegaU^2*gWall + coupling*gMixed;
portDetuning = omegaU^2*omegaW^2 - coupling^2;
staticTransfer = portNumerator^2/portDetuning^2;

wallLeg = gWall/(omegaW^2*Sqrt[k]);
mixedLeg = gMixed/(omegaU*omegaW*Sqrt[k]);
interference = coupling/(omegaU*omegaW);
shapeFromHats = FullSimplify[
  (wallLeg + interference*mixedLeg)/(1 - interference^2),
  Assumptions -> $Assumptions
];

Print["P   = ", fmt[portNumerator]];
Print["Delta = ", fmt[portDetuning]];
Print["T   = ", fmt[shapeFromHats]];
expectZero["N0/K - T^2", staticTransfer/k - shapeFromHats^2];

banner["Weak-axisymmetric slope identity"];

Clear[kk, xWall, xMixed, xCoupling, kappa1, wSlope, uSlope, cSlope];
$Assumptions = Element[{kk, xWall, xMixed, xCoupling, kappa1, wSlope, uSlope, cSlope}, Reals] &&
  kk > 0 && xWall > 0 && xMixed > 0;

shape = (xWall + xCoupling*xMixed)/(1 - xCoupling^2);
tauFromShape = logSlope[
  shape,
  {{xWall, wSlope}, {xMixed, uSlope}, {xCoupling, cSlope}}
];

alphaHat = xWall/(xWall + xCoupling*xMixed);
betaHat = xCoupling*xMixed/(xWall + xCoupling*xMixed);
tauFormula = FullSimplify[
  alphaHat*wSlope + betaHat*(uSlope + cSlope) +
    2*xCoupling^2*cSlope/(1 - xCoupling^2),
  Assumptions -> $Assumptions
];
nuFromFactor = logSlope[
  kk*shape^2,
  {{kk, kappa1}, {xWall, wSlope}, {xMixed, uSlope}, {xCoupling, cSlope}}
];

Print["tau_shape = ", fmt[tauFromShape]];
Print["nu_from_factor = ", fmt[nuFromFactor]];
expectZero["tau directional - alpha/beta form", tauFromShape - tauFormula];
expectZero["nu_factor - (kappa1 + 2 tau)", nuFromFactor - (kappa1 + 2*tauFromShape)];

banner["Exact equivalence to the Stage 176/160/161 slippage formulas"];

iHat = xCoupling*xMixed/xWall;
hHat = xCoupling^2;
mSlip = wSlope;
iSlip = (uSlope + cSlope) - wSlope;
hSlip = 2*cSlope;
tauSlippage = FullSimplify[
  mSlip + iHat*iSlip/(1 + iHat) + hHat*hSlip/(1 - hHat),
  Assumptions -> $Assumptions
];
expectZero["tau - slippage form", tauFromShape - tauSlippage];
expectZero["(nu-kappa1) - 2*tau_slippage", (nuFromFactor - kappa1) - 2*tauSlippage];

banner["Weighted defect identity"];

Clear[
  rho1, rho2,
  xWall1, xWall2, xWall3, xMixed1, xMixed2, xMixed3,
  xCoupling1, xCoupling2, xCoupling3,
  wSlope1, wSlope2, wSlope3, uSlope1, uSlope2, uSlope3,
  cSlope1, cSlope2, cSlope3
];
$Assumptions = Element[
    {
      rho1, rho2, kappa1,
      xWall1, xWall2, xWall3, xMixed1, xMixed2, xMixed3,
      xCoupling1, xCoupling2, xCoupling3,
      wSlope1, wSlope2, wSlope3, uSlope1, uSlope2, uSlope3,
      cSlope1, cSlope2, cSlope3
    },
    Reals
  ] && xWall1 > 0 && xWall2 > 0 && xWall3 > 0 &&
    xMixed1 > 0 && xMixed2 > 0 && xMixed3 > 0 && kk > 0;

rho3 = 1 - rho1 - rho2;
portRules = {
  {xWall -> xWall1, xMixed -> xMixed1, xCoupling -> xCoupling1,
    wSlope -> wSlope1, uSlope -> uSlope1, cSlope -> cSlope1},
  {xWall -> xWall2, xMixed -> xMixed2, xCoupling -> xCoupling2,
    wSlope -> wSlope2, uSlope -> uSlope2, cSlope -> cSlope2},
  {xWall -> xWall3, xMixed -> xMixed3, xCoupling -> xCoupling3,
    wSlope -> wSlope3, uSlope -> uSlope3, cSlope -> cSlope3}
};

tauPorts = FullSimplify[tauFromShape /. portRules, Assumptions -> $Assumptions];
nuPorts = FullSimplify[nuFromFactor /. portRules, Assumptions -> $Assumptions];

Do[
  expectZero[
    "nu_" <> ToString[j] <> " - (kappa1 + 2 tau_" <> ToString[j] <> ")",
    nuPorts[[j]] - (kappa1 + 2*tauPorts[[j]])
  ],
  {j, 1, 3}
];

xi = rho1*(nuPorts[[1]] - kappa1) +
  rho2*(nuPorts[[2]] - kappa1) +
  rho3*(nuPorts[[3]] - kappa1);
xiExpected = 2*(rho1*tauPorts[[1]] + rho2*tauPorts[[2]] + rho3*tauPorts[[3]]);
expectZero["Xi_1 - 2 weighted tau", xi - xiExpected];

Print[""];
Print["Carry-forward formulas:"];
Print["  N0^(r) = K * T_r^2"];
Print["  T_r = (Ghat_W + Rhat Ghat_U)/(1-Rhat^2)"];
Print["  nu_r = kappa1 + 2 tau_r"];
Print["  tau_r = m_r + I_r/(1+I_r) i_r + H_r/(1-H_r) h_r"];
Print["  Xi_1 = 2 sum_r rho_r^(N) tau_r"];

Print[""];
Print["Stage 179 Mathematica audit passed."];

Exit[0];
