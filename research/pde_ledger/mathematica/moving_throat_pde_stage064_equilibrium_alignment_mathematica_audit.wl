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

banner["STAGE 064 — EQUILIBRIUM ALIGNMENT"];

Needs["VariationalMethods`"];

Clear[gPhi, kX, npp, i1, i2, hw, L, y];
$Assumptions =
  Element[{gPhi, kX, npp, i1, i2, hw, L}, Reals] &&
  gPhi > 0 && kX > 0 && npp > 0 && i1 > 0 && i2 > 0 && hw > 0 && L > 0;

(* --- Local static linear-response closure by variational derivative --- *)
Clear[yLoc, chiPhiLoc, hLoc, sigmaFun];
energyDensity = (1/2) hLoc[yLoc] sigmaFun[yLoc]^2 - gPhi chiPhiLoc[yLoc] sigmaFun[yLoc];
eulerLagrange = VariationalD[energyDensity, sigmaFun[yLoc], yLoc];
closureSolutions = Solve[eulerLagrange == 0, sigmaFun[yLoc]];
If[Length[closureSolutions] =!= 1, fail["linear-response minimiser must be unique"]];
chiSigmaClosure = sigmaFun[yLoc] /. First[closureSolutions];
Print["closure chi_sigma = ", fmt[chiSigmaClosure]];
expectZero["closure law chi_sigma = g_phi chi_phi/H",
  chiSigmaClosure - gPhi*chiPhiLoc[yLoc]/hLoc[yLoc]];

(* --- Lorentzian integral overlap identities --- *)
Clear[y];
chiPhiL[z_] := 1/(1 + z^2/L^2);
hLayer[z_] := hw;
chiSigmaLayer[z_] := gPhi*chiPhiL[z]/hLayer[z];

nppInt = FullSimplify[
  Integrate[chiPhiL[y]^2, {y, -Infinity, Infinity},
    Assumptions -> L > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
i1Int = FullSimplify[
  Integrate[chiPhiL[y]^2/hLayer[y], {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
i2Int = FullSimplify[
  Integrate[chiPhiL[y]^2/hLayer[y]^2, {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
ospInt = FullSimplify[
  Integrate[chiPhiL[y]*chiSigmaLayer[y], {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
nssInt = FullSimplify[
  Integrate[chiSigmaLayer[y]^2, {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];

Print["Npp_int = ", fmt[nppInt]];
Print["I1_int  = ", fmt[i1Int]];
Print["I2_int  = ", fmt[i2Int]];
expectZero["overlap O = g_phi * I1 (Lorentzian integral form)", ospInt - gPhi*i1Int];
expectZero["overlap N_ss = g_phi^2 * I2 (Lorentzian integral form)", nssInt - gPhi^2*i2Int];

osp = FullSimplify[gPhi*i1, Assumptions -> $Assumptions];
nss = FullSimplify[gPhi^2*i2, Assumptions -> $Assumptions];
c2 = FullSimplify[osp^2/(npp*nss), Assumptions -> $Assumptions];
gEq = FullSimplify[gPhi^2*i1/kX, Assumptions -> $Assumptions];

Print["O_(sigma phi) = ", fmt[osp]];
Print["N_(sigma sigma) = ", fmt[nss]];
Print["C_(sigma phi)^2 = ", fmt[c2]];
Print["G_eq = ", fmt[gEq]];
expectZero["coherence formula C^2 = I1^2/(Npp I2)", c2 - i1^2/(npp*i2)];

banner["MATCHED-LAYER INTEGRAL REDUCTION (Lorentzian profile)"];

expectZero["matched-layer I1 reduction", i1Int - nppInt/hw];
expectZero["matched-layer I2 reduction", i2Int - nppInt/hw^2];
c2Const = FullSimplify[(gPhi*i1Int)^2/(nppInt*gPhi^2*i2Int), Assumptions -> $Assumptions];
gEqConst = FullSimplify[gPhi^2*i1Int/kX, Assumptions -> $Assumptions];
Print["C^2 | H=const = ", fmt[c2Const]];
Print["G_eq | H=const = ", fmt[gEqConst]];
expectZero["matched-layer coherence", c2Const - 1];
expectZero["matched-layer gain vs Stage-062 best-alignment formula", gEqConst - gPhi^2*nppInt/(kX*hw)];

banner["CONTINUOUS CAUCHY BOUND CHECK"];

hVariable[z_] := hw*(1 + z^2/L^2);
fCauchy[z_] := chiPhiL[z]/Sqrt[hVariable[z]];
gCauchy[z_] := chiPhiL[z]/hVariable[z];
i1Var = FullSimplify[
  Integrate[chiPhiL[y]^2/hVariable[y], {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
i2Var = FullSimplify[
  Integrate[chiPhiL[y]^2/hVariable[y]^2, {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
fgIntegral = FullSimplify[
  Integrate[fCauchy[y]*gCauchy[y], {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
ffIntegral = FullSimplify[
  Integrate[fCauchy[y]^2, {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
ggIntegral = FullSimplify[
  Integrate[gCauchy[y]^2, {y, -Infinity, Infinity},
    Assumptions -> L > 0 && hw > 0, GenerateConditions -> False],
  Assumptions -> $Assumptions
];
pairGap = FullSimplify[ffIntegral*ggIntegral - fgIntegral^2, Assumptions -> $Assumptions];
pairExpected = FullSimplify[
  L^2*(15*Pi^2/128 - 256/225)/hw^3,
  Assumptions -> $Assumptions
];
cauchyGap = FullSimplify[nppInt*i2Var - i1Var^2, Assumptions -> $Assumptions];
cauchyExpected = FullSimplify[Pi^2*L^2/(64*hw^2), Assumptions -> $Assumptions];
Print["continuous pair Cauchy gap = ", fmt[pairGap]];
expectZero["continuous f-g Cauchy-Schwarz residual", pairGap - pairExpected];
If[!TrueQ[FullSimplify[pairGap >= 0, Assumptions -> $Assumptions]],
  fail["continuous f-g Cauchy gap is not nonnegative", pairGap]
];
Print["N_pp I2 - I1^2 = ", fmt[cauchyGap]];
Print["expected continuous gap = ", fmt[cauchyExpected]];
expectZero["continuous Cauchy bound C^2 <= 1", cauchyGap - cauchyExpected];
If[!TrueQ[FullSimplify[cauchyGap >= 0, Assumptions -> $Assumptions]],
  fail["continuous Cauchy gap is not nonnegative", cauchyGap]
];

banner["GENERAL-H EQUILIBRIUM SOFTENING CHECK"];

(* General-H equilibrium softening derived from the variational source-energy.    *)
(* The aligned closure is chi_sigma[y] = gPhi chiPhi[y] / H[y]; the parent source *)
(* self-energy is (1/2) Integrate[H[y] chi_sigma[y]^2, y].                        *)
Clear[chiPhi, hFun, ySoft];
chiSigmaFun[z_] := gPhi*chiPhi[z]/hFun[z];
$Assumptions = Element[{ySoft, gPhi}, Reals] && gPhi > 0;
(* Mathematica's Integrate does not pull constant gPhi^2 outside Integrate[...]   *)
(* when the integrand contains unspecified functions chiPhi[y] and hFun[y].       *)
(* Verify the integrand equality first; the integral equality then follows.       *)
thetaIntegrand = FullSimplify[hFun[ySoft]*chiSigmaFun[ySoft]^2];
lambdaIntegrand = FullSimplify[gPhi*chiPhi[ySoft]*chiSigmaFun[ySoft]];
i1Integrand = chiPhi[ySoft]^2/hFun[ySoft];
expectZero["general Theta integrand equals gPhi^2 chiPhi^2/hFun",
  FullSimplify[thetaIntegrand - gPhi^2*i1Integrand]];
expectZero["general Lambda integrand equals gPhi^2 chiPhi^2/hFun",
  FullSimplify[lambdaIntegrand - gPhi^2*i1Integrand]];
(* By the integrand equality, the integrals are equal too.                        *)
i1Integral = Integrate[i1Integrand, {ySoft, -Infinity, Infinity}];
thetaGeneral = gPhi^2*i1Integral;
lambdaGeneral = gPhi^2*i1Integral;
(* Equilibrium softening = Lambda^2/Theta = (gPhi^2 I_1)^2 / (gPhi^2 I_1) = gPhi^2 I_1 *)
softGeneral = FullSimplify[lambdaGeneral^2/thetaGeneral];
expectZero["general equilibrium softening equals gPhi^2 I_1",
  softGeneral - gPhi^2*i1Integral];

banner["ELIMINATED-SOURCE SOFTENING CHECK"];

Clear[sourceCoeff, mixCoeff, supportAmp, sourceAmp];
$Assumptions =
  Element[{sourceCoeff, mixCoeff, supportAmp, sourceAmp, kX}, Reals] &&
  sourceCoeff > 0 && mixCoeff > 0 && kX > 0;

sourceEnergy = 1/2*sourceCoeff*sourceAmp^2 - mixCoeff*sourceAmp*supportAmp + 1/2*kX*supportAmp^2;
sourceStat = FullSimplify[
  sourceAmp /. First[Solve[D[sourceEnergy, sourceAmp] == 0, sourceAmp]],
  Assumptions -> $Assumptions
];
effEnergy = FullSimplify[sourceEnergy /. sourceAmp -> sourceStat, Assumptions -> $Assumptions];

Print["source_stat = ", fmt[sourceStat]];
Print["F_eff = ", fmt[effEnergy]];
expectZero[
  "effective support softening",
  effEnergy - 1/2*(kX - mixCoeff^2/sourceCoeff)*supportAmp^2
];

Print[""];
Print["Stage 064 Mathematica audit passed."];

Exit[0];
