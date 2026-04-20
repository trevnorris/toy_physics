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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION"];

Clear[g, r, s, xi, csw, cs, a, lW, vw0, tM, epsL, epsv, epsT, sigma0, sCan, dSigma0, dS, sigma, dKappaW, dGammaW];
$Assumptions =
  Element[
    {g, r, s, xi, csw, cs, a, lW, vw0, tM, epsL, epsv, epsT, sigma0, sCan, dSigma0, dS, sigma, dKappaW, dGammaW},
    Reals
  ] && g > 0 && r > 0 && s == Sqrt[1 + r^2];

aCoeff = g + 1/(4*s);
bCoeff = 1/(2*s);
cCoeff = 2*g + 3/(4*s);

deltaPerp = Expand[
  aCoeff*xi
  + 3*aCoeff*csw
  + bCoeff*cs
  - g*tM
  - (g + bCoeff)*vw0
  - 2*aCoeff*a
  - cCoeff*lW
];

Print["delta_perp = ", fmt[deltaPerp]];

banner["Stage-148 lower-branch transport laws"];
tMBr = xi/2 + 3*csw/2 - cs - 3*a/2;
vw0Br = xi/2 + 3*csw/2 + cs - 5*a/2;
lWBr = a;
Print["Tm_br  = ", fmt[tMBr]];
Print["vw0_br = ", fmt[vw0Br]];
Print["LW_br  = ", fmt[lWBr]];
expectZero[
  "bundle tangency (delta_perp on exact lower branch)",
  deltaPerp /. {tM -> tMBr, vw0 -> vw0Br, lW -> lWBr}
];

banner["Introduce the three off-bundle slippages"];
deltaPerpSlip = Expand[
  deltaPerp /. {
    tM -> tMBr + epsT,
    vw0 -> vw0Br + epsv,
    lW -> lWBr + epsL
  }
];
epsPerp = g*epsT + (g + bCoeff)*epsv + cCoeff*epsL;
Print["delta_perp with slippages = ", fmt[deltaPerpSlip]];
expectZero["delta_perp + eps_perp", deltaPerpSlip + epsPerp];

banner["Mouth-bias transport"];
deltaPi = (1 - sCan/4)*dSigma0 - sigma0*dS/4 + sigma0*sCan*deltaPerpSlip/s;
deltaPiExpected = (1 - sCan/4)*dSigma0 - sigma0*dS/4 - sigma0*sCan*epsPerp/s;
expectZero["deltaPi transport identity", deltaPi - deltaPiExpected];

banner["Outlet defects"];
dC = 16*sigma*deltaPerpSlip/s;
dE2 = sigma*(16*deltaPerpSlip/s - 9*dKappaW)/(27*(1 - sigma));
dE4 = sigma*(80*deltaPerpSlip/s - 72*dKappaW)/(243*(1 - sigma));
deltaQ = sigma*(16*deltaPerpSlip/s - 27*dGammaW)/(3*(1 - sigma));
expectZero["dC + 16 sigma eps_perp / s", dC + 16*sigma*epsPerp/s];
expectZero[
  "dE2 - sigma(-16 eps_perp/s - 9 dKappaW)/(27(1-sigma))",
  dE2 - sigma*(-16*epsPerp/s - 9*dKappaW)/(27*(1 - sigma))
];
expectZero[
  "dE4 - sigma(-80 eps_perp/s - 72 dKappaW)/(243(1-sigma))",
  dE4 - sigma*(-80*epsPerp/s - 72*dKappaW)/(243*(1 - sigma))
];
expectZero[
  "DeltaQ - sigma(-16 eps_perp/s - 27 dGammaW)/(3(1-sigma))",
  deltaQ - sigma*(-16*epsPerp/s - 27*dGammaW)/(3*(1 - sigma))
];

banner["Numeric Family-1 coefficients"];
rNum = SetPrecision[1.77799353547498, 30];
gNum = SetPrecision[0.758035078944663, 30];
sNum = Sqrt[1 + rNum^2];
bNum = 1/(2*sNum);
cNum = 2*gNum + 3/(4*sNum);
sigma0Num = SetPrecision[4.651033550168876, 30];
sCanNum = SetPrecision[0.6703621156734617, 30];

Print["coeff epsT in delta_perp = ", fmt[N[-gNum, 20]]];
Print["coeff epsv in delta_perp = ", fmt[N[-(gNum + bNum), 20]]];
Print["coeff epsL in delta_perp = ", fmt[N[-cNum, 20]]];
Print["coeff epsT in deltaPi = ", fmt[N[-sigma0Num*sCanNum*gNum/sNum, 20]]];
Print["coeff epsv in deltaPi = ", fmt[N[-sigma0Num*sCanNum*(gNum + bNum)/sNum, 20]]];
Print["coeff epsL in deltaPi = ", fmt[N[-sigma0Num*sCanNum*cNum/sNum, 20]]];

Print[""];
Print["Carry-forward formulas:"];
Print["  eps_L := d ln L_W - d ln a"];
Print["  eps_v := d ln v_w0 - [1/2 d ln(Z_q/rho_w) + 3/2 d ln c_sw + d ln c_s - 5/2 d ln a]"];
Print["  eps_T := d ln T_m - [1/2 d ln(Z_q/rho_w) + 3/2 d ln c_sw - d ln c_s - 3/2 d ln a]"];
Print["  delta_perp = -[ g_* eps_T + (g_* + 1/(2 sqrt(1+r_*^2))) eps_v + (2 g_* + 3/(4 sqrt(1+r_*^2))) eps_L ]"];
Print["  All first-order off-bundle conservative-even defects depend only on eps_perp."];

Print[""];
Print["Stage 151 Mathematica audit passed."];

Exit[0];
