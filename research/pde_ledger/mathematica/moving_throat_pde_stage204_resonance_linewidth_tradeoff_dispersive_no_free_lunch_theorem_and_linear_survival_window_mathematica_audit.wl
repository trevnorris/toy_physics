ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

normalizeExpr[expr_] := FullSimplify[expr, Assumptions -> $Assumptions];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res == 0], pass[name], fail[name, res]];
];

banner["STAGE 204 — RESONANCE / LINEWIDTH TRADEOFF"];

Clear[
  delta, GammaOut, Zstar, Fprime, Numstar, portPi,
  K, M, Ccoef, varpi, OmegaU, OmegaW, R, GU, GW, omega,
  D0prime, Nstar, Aabs, gamma, r, eta, Sfam, omegaDrive,
  omegaStar, Qfac, DeltaVreq
];

$Assumptions =
  Element[
    {
      delta, GammaOut, Zstar, Fprime, Numstar,
      K, M, Ccoef, varpi, OmegaU, OmegaW, R, GU, GW, omega,
      D0prime, Nstar, Aabs, gamma, r, eta, Sfam, omegaDrive,
      omegaStar, Qfac, DeltaVreq
    },
    Reals
  ] &&
  GammaOut > 0 && Zstar > 0 && Fprime > 0 &&
  varpi != 0 && OmegaU != 0 && OmegaW != 0 && R != 0 &&
  GU != 0 && GW != 0 &&
  D0prime > 0 && Nstar > 0 &&
  Aabs > 0 && gamma > 0 && r > 0 && eta > 0 && eta <= 1 &&
  Sfam > 0 && omegaDrive > 0 && omegaStar > 0 && Qfac > 0 &&
  DeltaVreq > 0;

subbanner["I. Generic simple-pole normal form"];

Flin = Fprime delta - portPi Zstar;
chiLin = Together[Numstar/Flin];

Astar = normalizeExpr[Numstar/Fprime];
gammaStar = normalizeExpr[GammaOut Zstar/Fprime];
chiBW = Together[Astar/(delta - portPi Zstar/Fprime)];

expectZero["chi_lin - chi_BW", chiLin - chiBW];

chiPassive = Together[chiLin /. portPi -> I GammaOut];
chiPassiveExpected = Together[Astar/(delta - I gammaStar)];

expectZero["passive simple-pole normal form", chiPassive - chiPassiveExpected];

Print["chi(omega) = ", fmt[chiPassiveExpected]];
Print["A_* = ", fmt[Astar]];
Print["gamma_* = ", fmt[gammaStar]];

subbanner["II. Exact Stage-203 derivative identity"];

KB = K - M omega^2 - Ccoef^2/(varpi^2 - omega^2);
Afun = OmegaU^2 - omega^2;
Wfun = OmegaW^2 - omega^2 - portPi;

DeltaPi = Afun Wfun - R^2;
QPi = GU^2 Wfun + 2 GU GW R + GW^2 Afun;
DPi = Together[KB - QPi/DeltaPi];
Nfun = Together[(Afun GW + R GU)^2/DeltaPi^2];

expectZero["dD_Pi/dPi + N(omega)", D[DPi, portPi] + Nfun];

subbanner["III. Wall-like pole specialization"];

Dlin = D0prime delta - portPi Nstar;
chiqqLin = Together[1/Dlin];
gammaWall = normalizeExpr[GammaOut Nstar/D0prime];
chiqqExpected = Together[(1/D0prime)/(delta - I gammaWall)];

expectZero[
  "wall-like pole specialization",
  (chiqqLin /. portPi -> I GammaOut) - chiqqExpected
];

Print["chi_qq(omega) = ", fmt[chiqqExpected]];
Print["gamma_wall = ", fmt[gammaWall]];

subbanner["IV. Exact line-shape algebra"];

chiR = Together[Aabs/(r gamma - I gamma)];
reR = normalizeExpr[ComplexExpand[Re[chiR]]];
imR = normalizeExpr[ComplexExpand[Im[chiR]]];

reExpected = normalizeExpr[Aabs r/(gamma (1 + r^2))];
imExpected = normalizeExpr[Aabs/(gamma (1 + r^2))];

expectZero["Re chi - expected", reR - reExpected];
expectZero["Im chi - expected", imR - imExpected];
expectZero["|Re|/|Im| - r", reExpected/imExpected - r];

f = normalizeExpr[r/(1 + r^2)];
expectZero[
  "maximum identity",
  (f - 1/2) + (r - 1)^2/(2 (1 + r^2))
];
expectZero[
  "low-loss identity",
  (f - eta/(1 + eta^2)) + (r - eta) (eta r - 1)/((1 + r^2) (1 + eta^2))
];

subbanner["V. Barrier / absorbed-power ratio"];

Udisp = normalizeExpr[-(1/2) reExpected Sfam^2];
Pabs = normalizeExpr[(1/2) omegaDrive imExpected Sfam^2];
ratioBarrier = normalizeExpr[Pabs/(omegaDrive (-Udisp))];

expectZero["P_abs / (omega |U_disp|) - 1/r", ratioBarrier - 1/r];

Print["U_disp = ", fmt[Udisp]];
Print["P_abs = ", fmt[Pabs]];
Print["P_abs / (omega |U_disp|) = ", fmt[ratioBarrier]];

subbanner["VI. Quality-factor detuning and linear survival window"];

detuningFraction = normalizeExpr[(r gamma)/omegaStar];
detuningQForm = normalizeExpr[detuningFraction /. gamma -> omegaStar/(2 Qfac)];
lowLossDetuning = normalizeExpr[detuningQForm /. r -> 1/eta];

expectZero["detuning fraction - r/(2 Q_*)", detuningQForm - r/(2 Qfac)];
expectZero["low-loss detuning boundary", lowLossDetuning - 1/(2 Qfac eta)];

UdispLowLossMax = normalizeExpr[(1/2) Aabs/gamma * eta/(1 + eta^2) * Sfam^2];
survivalLeft = normalizeExpr[Aabs/gamma * eta/(1 + eta^2) * Sfam^2];
residueRequirement = normalizeExpr[
  2 DeltaVreq/Sfam^2 * (1 + eta^2)/eta
];

expectZero["survival left side - 2 |U_disp|_max", survivalLeft - 2 UdispLowLossMax];
expectZero[
  "residue requirement saturates the survival window",
  residueRequirement * eta/(1 + eta^2) * Sfam^2 - 2 DeltaVreq
];

Print["|U_disp|_max in the low-loss window = ", fmt[UdispLowLossMax]];
Print["Necessary residue/linewidth ratio >= ", fmt[residueRequirement]];

subbanner["VII. Probe-only numerical slice"];

probeValues = {
  Aabs -> 7,
  gamma -> 2,
  r -> 3,
  eta -> 1/4,
  Sfam -> 5,
  omegaDrive -> 11,
  Qfac -> 40,
  omegaStar -> 80,
  DeltaVreq -> 3/5
};

Print["Re chi = ", N[reExpected /. probeValues]];
Print["Im chi = ", N[imExpected /. probeValues]];
Print["|Re|/|Im| = ", N[(reExpected/imExpected) /. probeValues]];
Print["|U_disp| = ", N[(-Udisp) /. probeValues]];
Print["P_abs = ", N[Pabs /. probeValues]];
Print["P_abs/(omega|U|) = ", N[ratioBarrier /. probeValues]];
Print["low-loss |U|_max = ", N[UdispLowLossMax /. probeValues]];
Print["required |A|/gamma = ", N[residueRequirement /. probeValues]];

Print[""];
Print["All Stage-204 identities, derivative formulas, tradeoff laws, and survival-window relations verified."];
