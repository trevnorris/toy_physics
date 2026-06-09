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

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

canonical[expr_, asm_] := Module[{res},
  res = FullSimplify[Together[Cancel[expr]], Assumptions -> asm];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> asm]
];

expectZeroUnder[name_String, expr_, asm_] := Module[{res},
  res = canonical[expr, asm];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectZero[name_String, expr_] := expectZeroUnder[name, expr, $Assumptions];

scaledFirstVariation[expr_, increments_List] := Module[{delta},
  delta = Total[(D[expr, #[[1]]] * #[[2]]) & /@ increments];
  FullSimplify[Together[Cancel[delta/expr]], Assumptions -> $Assumptions]
];

banner["STAGE 180 — EFFECTIVE TRANSFER-SHAPE COLLAPSE"];

Clear[kA, t1, t2, tau1, tau2];
$Assumptions = Element[{kA, t1, t2, tau1, tau2}, Reals] && kA > 0 && t1 > 0 && t2 > 0;
teff2Base = t1^2 + t2^2;
nPort1 = kA*t1^2;
nPort2 = kA*t2^2;
rho1 = FullSimplify[nPort1/(nPort1 + nPort2), Assumptions -> $Assumptions];
rho2 = FullSimplify[nPort2/(nPort1 + nPort2), Assumptions -> $Assumptions];
xiEff = scaledFirstVariation[teff2Base, {{t1, t1*tau1}, {t2, t2*tau2}}];
xiExpected = FullSimplify[2*(rho1*tau1 + rho2*tau2), Assumptions -> $Assumptions];
expectZero["multi-port effective-shape identity", xiEff - xiExpected];

banner["One-port continuum transfer shape"];
Clear[muW, muEta, kEta, kW, zW, rho, epsW, omegaW2, k0Sym, betaSym, transferSq];
$Assumptions = Element[{muW, muEta, kEta, kW, zW, rho, epsW, omegaW2}, Reals] &&
  muW > 0 && muEta > 0 && kEta > 0 && kW > 0 && zW > 0 &&
  omegaW2 > 0 && 1 + rho != 0 && 1 - epsW != 0;

continuumEquations = {
  k0Sym*muEta == kEta,
  betaSym*muEta*kW*(1 - epsW)^2 == muW*kEta*zW*(1 + rho)^2,
  transferSq*k0Sym == betaSym
};
continuumRules = First[Solve[continuumEquations, {k0Sym, betaSym, transferSq}]];
t2FromContinuum = FullSimplify[transferSq /. continuumRules, Assumptions -> $Assumptions];
t2MuKForm = FullSimplify[(muW/kW)*zW*(1 + rho)^2/(1 - epsW)^2, Assumptions -> $Assumptions];
t2OmegaForm = FullSimplify[zW*(1 + rho)^2/(omegaW2*(1 - epsW)^2), Assumptions -> $Assumptions];
expectZero["T^2 = beta0/K0 -> muW/KW form", t2FromContinuum - t2MuKForm];
expectZeroUnder[
  "T^2 = ZW(1+rho)^2 / [OmegaW^2 (1-epsW)^2]",
  t2FromContinuum - t2OmegaForm,
  $Assumptions && omegaW2 == kW/muW
];

banner["Selected-branch reformulation"];
Clear[gConst, cs, a, cSpeed, epsEta, rTarget, lambdaBranch, selectedTransferSq];
$Assumptions = Element[{gConst, cs, a, cSpeed, epsEta, rTarget, muW, kW, zW, rho, epsW}, Reals] &&
  gConst > 0 && cs > 0 && a > 0 && cSpeed > 0 && muW > 0 && kW > 0 &&
  zW > 0 && omegaW2 > 0 && rTarget != 0 && 1 + rho != 0 && 1 - epsW != 0;

frontFactor = FullSimplify[27*Pi^2*gConst*cs^5/(20*a^5*cSpeed^5), Assumptions -> $Assumptions];
selectedEquations = {
  lambdaBranch == frontFactor*omegaW2,
  rTarget*zW*(1 + rho)^2 == lambdaBranch*(1 - epsEta)*(1 - epsW)^2,
  selectedTransferSq == zW*(1 + rho)^2/(omegaW2*(1 - epsW)^2)
};
selectedRules = First[Solve[selectedEquations, {lambdaBranch, omegaW2, selectedTransferSq}]];
t2FromSelectedBranch = FullSimplify[selectedTransferSq /. selectedRules, Assumptions -> $Assumptions];
t2SelectedExpected = FullSimplify[frontFactor*(1 - epsEta)/rTarget, Assumptions -> $Assumptions];
expectZero["selected-branch T^2 identity", t2FromSelectedBranch - t2SelectedExpected];

banner["Weak-axisymmetric slope laws"];
Clear[zetaW, omegaW, rho1s, epsW1, eta1, r1];
$Assumptions = Element[
    {zetaW, omegaW, rho1s, epsW1, eta1, r1, zW, rho, epsW, omegaW2,
     gConst, cs, a, cSpeed, epsEta, rTarget},
    Reals
  ] &&
  zW > 0 && omegaW2 > 0 && gConst > 0 && cs > 0 && a > 0 && cSpeed > 0 &&
  rTarget != 0 && 1 + rho != 0 && 1 - epsW != 0;

directTransferShape = zW*(1 + rho)^2/(omegaW2*(1 - epsW)^2);
xiDirect = scaledFirstVariation[
  directTransferShape,
  {{zW, zW*zetaW}, {omegaW2, omegaW2*omegaW}, {rho, rho1s}, {epsW, epsW1}}
];
xiDirectExpected = FullSimplify[zetaW - omegaW + 2*rho1s/(1 + rho) + 2*epsW1/(1 - epsW), Assumptions -> $Assumptions];
expectZero["direct slope law", xiDirect - xiDirectExpected];

selectedTransferShape = frontFactor*(1 - epsEta)/rTarget;
xiSel = scaledFirstVariation[
  selectedTransferShape,
  {{epsEta, eta1}, {rTarget, rTarget*r1}}
];
xiSelExpected = FullSimplify[-eta1/(1 - epsEta) - r1, Assumptions -> $Assumptions];
expectZero["selected-branch slope law", xiSel - xiSelExpected];

Print[""];
Print["Carry-forward formulas:"];
Print["  T_eff^2 = sum_r T_r^2 = N_0/K"];
Print["  Xi_1    = delta ln(T_eff^2)/(eps lambda_A)"];
Print["  On the one-port continuum branch,"];
Print["    T^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-eps_W)^2]"];
Print["        = (27 pi^2 G c_s^5 / (20 a^5 c^5)) * (1-eps_eta)/R_target"];
Print["  Hence"];
Print["    Xi_1 = zeta_W - omega_W + 2 rho_1/(1+rho) + 2 epsW_1/(1-eps_W)"];
Print["        = - eta_1/(1-eps_eta) - R_1"];

Print[""];
Print["Stage 180 Mathematica audit passed."];

Exit[0];
