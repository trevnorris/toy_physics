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

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[Expand[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanScalar[expr];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = N[Abs[value - target], 20];
  Print[name, " diff = ", fmt[diff], " (tol = ", fmt[tol], ")"];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 252 -- VACUUM/LATTICE HEAT PARTITION"];

Clear[
  Gamma3vac, Gamma3lat, Gamma5vac, Gamma5lat, I2, I3, rV, s, T, Vin,
  sc, s0, mueta, t, phi, G3T, G5T
];

$Assumptions =
  Element[
    {
      Gamma3vac, Gamma3lat, Gamma5vac, Gamma5lat, I2, I3, rV, s, T,
      Vin, sc, s0, mueta
    },
    Reals
  ] &&
  Gamma3vac >= 0 && Gamma3lat >= 0 && Gamma5vac >= 0 && Gamma5lat >= 0 &&
  I2 > 0 && rV > 0 && s > 0 && T > 0 && Vin > 0 &&
  sc > 0 && s0 > 0 && mueta > 0;

limitAssumptions =
  Element[{Gamma3vac, Gamma3lat, Gamma5vac, Gamma5lat}, Reals] &&
  Gamma3vac >= 0 && Gamma3lat >= 0 && Gamma5vac >= 0 && Gamma5lat >= 0;

Gamma3tot = Gamma3vac + Gamma3lat;
Gamma5tot = Gamma5vac + Gamma5lat;

subbanner["M1. Energy ledger and partition fractions"];

Evac = Gamma3vac I2 + Gamma5vac I3;
Elat = Gamma3lat I2 + Gamma5lat I3;
Eexp = Evac + Elat;
shapeRule = I3 -> rV I2;

fvac = FullSimplify[Cancel[Together[(Evac/Eexp) /. shapeRule]], Assumptions -> $Assumptions];
flat = FullSimplify[Cancel[Together[(Elat/Eexp) /. shapeRule]], Assumptions -> $Assumptions];
fvacExpected = (Gamma3vac + Gamma5vac rV)/(Gamma3tot + Gamma5tot rV);
flatExpected = (Gamma3lat + Gamma5lat rV)/(Gamma3tot + Gamma5tot rV);

Print["E_vac = ", fmt[Evac]];
Print["E_lat = ", fmt[Elat]];
Print["E_exp = ", fmt[Eexp]];
Print["f_vac(rV) = ", fmt[fvac]];
Print["f_lat(rV) = ", fmt[flat]];

expectZero["M1 vacuum fraction", fvac - fvacExpected];
expectZero["M1 lattice fraction", flat - flatExpected];
expectZero["M1 partition sum", fvac + flat - 1];

subbanner["M2. Drift law"];

flatDrift = FullSimplify[Together[D[flat, rV]], Assumptions -> $Assumptions];
flatDriftExpected =
  (Gamma5lat Gamma3vac - Gamma3lat Gamma5vac)/(Gamma3tot + Gamma5tot rV)^2;

Print["D[f_lat, rV] = ", fmt[flatDrift]];
expectZero["M2 lattice drift law", flatDrift - flatDriftExpected];

subbanner["M3. Endpoint limits"];

flatAtZero = stripConditional[Limit[flat, rV -> 0, Assumptions -> limitAssumptions]];
flatAtInfinity =
  stripConditional[Limit[flat, rV -> Infinity, Assumptions -> limitAssumptions]];

Print["f_lat(0) = ", fmt[flatAtZero]];
Print["f_lat(infinity) = ", fmt[flatAtInfinity]];

expectZero["M3 lattice slow endpoint", flatAtZero - Gamma3lat/Gamma3tot];
expectZero["M3 lattice fast endpoint", flatAtInfinity - Gamma5lat/Gamma5tot];

subbanner["M4. Exponential shape integrals"];

Vevent = Vin Exp[s t];
I1exp = FullSimplify[
  Integrate[D[Vevent, {t, 1}]^2, {t, 0, T}, Assumptions -> {s > 0, T > 0, Vin > 0}],
  Assumptions -> $Assumptions
];
I2exp = FullSimplify[
  Integrate[D[Vevent, {t, 2}]^2, {t, 0, T}, Assumptions -> {s > 0, T > 0, Vin > 0}],
  Assumptions -> $Assumptions
];
I3exp = FullSimplify[
  Integrate[D[Vevent, {t, 3}]^2, {t, 0, T}, Assumptions -> {s > 0, T > 0, Vin > 0}],
  Assumptions -> $Assumptions
];

Print["I1exp = ", fmt[I1exp]];
Print["I2exp = ", fmt[I2exp]];
Print["I3exp = ", fmt[I3exp]];

expectZero["M4 I2 exponential integral", I2exp - Vin^2 s^3 (Exp[2 s T] - 1)/2];
expectZero["M4 I3 = s^2 I2", I3exp - s^2 I2exp];
expectZero["M4 shape quotient", I3exp/I2exp - s^2];

subbanner["M5. Event-equivalent rates"];

gammaVacEq = Gamma3vac s^2 + Gamma5vac s^4;
gammaLatEq = Gamma3lat s^2 + Gamma5lat s^4;
EvacExp = FullSimplify[Evac /. {I2 -> I2exp, I3 -> I3exp}, Assumptions -> $Assumptions];
ElatExp = FullSimplify[Elat /. {I2 -> I2exp, I3 -> I3exp}, Assumptions -> $Assumptions];

Print["gamma_vac^eq(s) = ", fmt[gammaVacEq]];
Print["gamma_lat^eq(s) = ", fmt[gammaLatEq]];

expectZero["M5 vacuum event-equivalent rate", EvacExp - gammaVacEq I1exp];
expectZero["M5 lattice event-equivalent rate", ElatExp - gammaLatEq I1exp];

subbanner["M6. Safe-edge exported-energy theorem"];

gammaEffEq = Gamma3tot s^2 + Gamma5tot s^4;
EexpEvent = FullSimplify[Eexp /. {I2 -> I2exp, I3 -> I3exp}, Assumptions -> $Assumptions];
EexpSafe = FullSimplify[EexpEvent /. {s -> sc, T -> 1/sc}, Assumptions -> $Assumptions];
safeCombo = Gamma3tot sc^3 + Gamma5tot sc^5;
safeRule = safeCombo -> mueta (s0^2 - sc^2);
safeEnergyRaw = Vin^2 (E^2 - 1) safeCombo/2;
safeEnergyReduced = FullSimplify[safeEnergyRaw /. safeRule, Assumptions -> $Assumptions];

Print["EexpSafe = ", fmt[EexpSafe]];
Print["safeCombo = ", fmt[safeCombo]];
Print["safeEnergyReduced = ", fmt[safeEnergyReduced]];

expectZero["M6 safe-edge energy before rule", EexpSafe - safeEnergyRaw];
expectZero[
  "M6 safe-edge energy after safe rule",
  safeEnergyReduced - Vin^2 (E^2 - 1) mueta (s0^2 - sc^2)/2
];

subbanner["M7. Safe-edge rate identity"];

safeRateRaw = safeCombo/sc;
safeRateReduced = FullSimplify[safeRateRaw /. safeRule, Assumptions -> $Assumptions];

Print["safeRateRaw = ", fmt[safeRateRaw]];
Print["safeRateReduced = ", fmt[safeRateReduced]];

expectZero["M7 safe rate raw quotient", safeRateRaw - (Gamma3tot sc^2 + Gamma5tot sc^4)];
expectZero[
  "M7 safe rate after safe rule",
  safeRateReduced - mueta (s0^2 - sc^2)/sc
];

subbanner["M8. Three-to-one surface and phi-family"];

flatSc = FullSimplify[flat /. rV -> sc^2, Assumptions -> $Assumptions];
threeOneResid = Together[flatSc - 3/4];
threeOneNumerator = Numerator[Cancel[threeOneResid]];
threeOneSurface = Expand[(Gamma3lat + Gamma5lat sc^2) - 3 (Gamma3vac + Gamma5vac sc^2)];
threeOneFactor = 1;

flatPhi = FullSimplify[
  Cancel[
    Together[
      flat /. {
        Gamma3lat -> phi G3T,
        Gamma3vac -> (1 - phi) G3T,
        Gamma5lat -> phi G5T,
        Gamma5vac -> (1 - phi) G5T
      }
    ]
  ]
];
fvacPhi = FullSimplify[
  Cancel[
    Together[
      fvac /. {
        Gamma3lat -> phi G3T,
        Gamma3vac -> (1 - phi) G3T,
        Gamma5lat -> phi G5T,
        Gamma5vac -> (1 - phi) G5T
      }
    ]
  ]
];

Print["f_lat(sc) - 3/4 numerator = ", fmt[threeOneNumerator]];
Print["3:1 surface = ", fmt[threeOneSurface]];
Print["f_lat(phi family) = ", fmt[flatPhi]];
Print["f_vac(phi family) = ", fmt[fvacPhi]];

expectZero[
  "M8 three-to-one numerator is exact positive factor times surface",
  threeOneNumerator - threeOneFactor threeOneSurface
];
expectZero["M8 phi-family lattice fraction", flatPhi - phi];
expectZero["M8 phi-family vacuum fraction", fvacPhi - (1 - phi)];

subbanner["M9. Session-IV benchmark reproduction"];

phiVal = 3/4;
fracVacNum = N[fvacPhi /. phi -> phiVal, 20];
fracLatNum = N[flatPhi /. phi -> phiVal, 20];
tCrossNum = 1.82169718;
s0Num = 6.94311167;
muetaNum = 1.0;
eDissNum = 0.01033460;

scNum = N[1.0/tCrossNum, 20];
sc2Num = N[scNum^2, 20];
gammaEffSafeNum = N[muetaNum (s0Num^2 - sc2Num)/scNum, 20];
gammaVacSafeNum = N[fracVacNum gammaEffSafeNum, 20];
gammaLatSafeNum = N[fracLatNum gammaEffSafeNum, 20];
safePrefactorNum = N[1/2 (E^2 - 1.0) muetaNum (s0Num^2 - sc2Num), 20];
VinMatchNum = N[Sqrt[eDissNum/safePrefactorNum], 20];
EvacBenchNum = N[fracVacNum eDissNum, 20];
ElatBenchNum = N[fracLatNum eDissNum, 20];

Print["phiVal = ", fmt[phiVal]];
Print["fracVacNum = ", fmt[fracVacNum]];
Print["fracLatNum = ", fmt[fracLatNum]];
Print["scNum = ", fmt[scNum]];
Print["gammaEffSafeNum = ", fmt[gammaEffSafeNum]];
Print["safePrefactorNum = ", fmt[safePrefactorNum]];
Print["VinMatchNum = ", fmt[VinMatchNum]];
Print["EvacBenchNum = ", fmt[EvacBenchNum]];
Print["ElatBenchNum = ", fmt[ElatBenchNum]];

expectApprox["M9 sc benchmark", scNum, 0.5489386551062235, 10^-12];
expectApprox["M9 gamma_eff safe benchmark", gammaEffSafeNum, 87.26925234614843, 10^-9];
expectApprox["M9 gamma_vac safe benchmark", gammaVacSafeNum, 21.81731308653711, 10^-9];
expectApprox["M9 gamma_lat safe benchmark", gammaLatSafeNum, 65.45193925961132, 10^-9];
expectApprox["M9 safe prefactor benchmark", safePrefactorNum, 153.03535490769042, 10^-9];
expectApprox["M9 Vin match benchmark", VinMatchNum, 0.008217712598897912, 10^-15];
expectApprox["M9 vacuum energy benchmark", EvacBenchNum, 0.00258365, 10^-12];
expectApprox["M9 lattice energy benchmark", ElatBenchNum, 0.00775095, 10^-12];
expectApprox["M9 total energy benchmark", EvacBenchNum + ElatBenchNum, eDissNum, 10^-12];

banner["Stage 252 Mathematica audit result"];
Print["All Mathematica checks passed for M1-M9."];
