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

banner["STAGE 054 — CANONICAL TANH-WALL BRANCH"];

Clear[xi, t, a, L, ell, rhoW, cSw, V0, m, hbar, chiS, lambdaEll, upsilonW];
$Assumptions =
  Element[{xi, t}, Reals] &&
  Element[{a, L, ell, rhoW, cSw, V0, m, hbar, chiS, lambdaEll, upsilonW}, Reals] &&
  And @@ Thread[{a, L, ell, rhoW, cSw, V0, m, hbar, chiS, lambdaEll, upsilonW} > 0];

f = (1 + Tanh[xi])/2;
fp = FullSimplify[D[f, xi], Assumptions -> $Assumptions];
fpp = FullSimplify[D[fp, xi], Assumptions -> $Assumptions];

ifDirect = FullSimplify[Integrate[fp^2, {xi, -Infinity, Infinity}], Assumptions -> $Assumptions];
ifSub = FullSimplify[Integrate[(1 - t^2)/4, {t, -1, 1}], Assumptions -> $Assumptions];
igDirect = FullSimplify[Integrate[fpp^2, {xi, -Infinity, Infinity}], Assumptions -> $Assumptions];
igSub = FullSimplify[Integrate[t^2*(1 - t^2), {t, -1, 1}], Assumptions -> $Assumptions];

Print["f'(xi) = ", fmt[fp]];
Print["f''(xi) = ", fmt[fpp]];
Print["I_f = ", fmt[ifDirect]];
Print["I_f_sub = ", fmt[ifSub]];
Print["I_g = ", fmt[igDirect]];
Print["I_g_sub = ", fmt[igSub]];

expectZero["I_f - 1/3", ifDirect - 1/3];
expectZero["I_f direct - substitution", ifDirect - ifSub];
expectZero["I_g - 4/15", igDirect - 4/15];
expectZero["I_g direct - substitution", igDirect - igSub];
expectZero["I_g/I_f - 4/5", igDirect/ifDirect - 4/5];

hw = FullSimplify[m*cSw^2/rhoW, Assumptions -> $Assumptions];
tx = FullSimplify[Pi*a^2*ell*ifDirect*hbar^2/(m*rhoW), Assumptions -> $Assumptions];
kx = FullSimplify[
  4*Pi*a^2*ell*ifDirect*hw + Pi*a^2*igDirect*hbar^2/(m*rhoW*ell),
  Assumptions -> $Assumptions
];
j1 = FullSimplify[ifDirect/hw, Assumptions -> $Assumptions];
wwall = FullSimplify[4*rhoW^2*V0^2*L^2/(hbar^2*cSw^2*ell^2), Assumptions -> $Assumptions];

Print["T_X = ", fmt[tx]];
Print["K_X = ", fmt[kx]];
Print["J_1 = ", fmt[j1]];
Print["W_wall = ", fmt[wwall]];

expectZero["T_X exact formula", tx - Pi*a^2*ell*hbar^2/(3*m*rhoW)];
expectZero[
  "K_X exact formula",
  kx - 4*Pi*a^2*(5*m^2*cSw^2*ell^2 + hbar^2)/(15*ell*m*rhoW)
];
expectZero["J_1 exact formula", j1 - rhoW/(3*m*cSw^2)];

km = FullSimplify[tx/ell, Assumptions -> $Assumptions];
eta = FullSimplify[km*L/tx, Assumptions -> $Assumptions];
Print["K_m = ", fmt[km]];
Print["eta = ", fmt[eta]];
kmExpected = Pi*a^2*hbar^2/(3*m*rhoW);
expectZero["K_m - pi a^2 hbar^2 / (3 m rhoW)", km - kmExpected];
expectZero["eta - L/ell (from closed-form K_m)", (kmExpected*L/tx) - L/ell];

kappa = FullSimplify[kx*L^2/tx, Assumptions -> $Assumptions];
kappaExpected = FullSimplify[4*(m*cSw*L/hbar)^2 + (4/5)*(L/ell)^2, Assumptions -> $Assumptions];
wExpected = FullSimplify[(4*rhoW^2*V0^2/(hbar^2*cSw^2))*(L/ell)^2, Assumptions -> $Assumptions];

Print["kappa = ", fmt[kappa]];
Print["kappa reduced = ", fmt[4*chiS^2 + (4/5)*lambdaEll^2]];
expectZero["kappa reduced law", kappa - kappaExpected];

Print["W_wall reduced = ", fmt[upsilonW*lambdaEll^2]];
expectZero["W_wall reduced law", wwall - wExpected];

banner["STAGE 054 THEOREM LEDGER"];
Print["I_f = 1/3"];
Print["I_g = 4/15"];
Print["I_g/I_f = 4/5"];
Print["eta = L/ell under K_m = T_X/ell"];
Print["kappa = 4 chi_s^2 + (4/5) Lambda_ell^2"];
Print["W_wall = Upsilon_w Lambda_ell^2"];

Exit[0];
