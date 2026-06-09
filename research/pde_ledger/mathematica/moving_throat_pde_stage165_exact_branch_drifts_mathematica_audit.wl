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

banner["STAGE 165 — EXACT LOWER-BRANCH DRIFT LAWS"];

Clear[dZ, drho, dcsw, dcs, dT, dv, da, dLW, r, a];
$Assumptions = Element[{dZ, drho, dcsw, dcs, dT, dv, da, dLW, r, a}, Reals] && r > 0 && a > 0;

lwLaw = Pi*a*Sqrt[1 + r^2]/(2*Sqrt[3]);
Print["D/N law: L_W = ", fmt[lwLaw]];
expectZero["d ln L_W - d ln a at fixed r_*", a*D[Log[lwLaw], a] - 1];

dEll = -dcsw;
dLam = 2*da + dEll + dLW/2 + dv;
dKq = dZ + 2*dcs - 2*dLW;
dKs = 2*da - drho - dEll;
dJs = 2*da + dEll;
dGs = dT + dJs;
dGq = dZ - (3*dLW)/2;
chanR = FullSimplify[dKs + dKq - 2*dLam, Assumptions -> $Assumptions];
chanG = FullSimplify[dGq + dKs - dGs - dLam, Assumptions -> $Assumptions];
expectZero["eqR matches Stage-164 fixed-r channel", chanR - (dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW)];
expectZero["eqG matches Stage-164 fixed-g channel", chanG - (dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW)];
eqR = (chanR == 0);
eqG = (chanG == 0);
sol = Solve[{eqR, eqG}, {dv, dT}, Reals][[1]];
dvSol = FullSimplify[dv /. sol, Assumptions -> $Assumptions];
dTSol = FullSimplify[dT /. sol, Assumptions -> $Assumptions];

banner["Exact fixed-r / fixed-g drift laws"];
Print["d ln v_w0 = ", fmt[dvSol]];
Print["d ln T_m  = ", fmt[dTSol]];

dvDN = FullSimplify[dvSol /. dLW -> da, Assumptions -> $Assumptions];
dTDN = FullSimplify[dTSol /. dLW -> da, Assumptions -> $Assumptions];

banner["After using d ln L_W = d ln a"];
Print["d ln v_w0 = ", fmt[dvDN]];
Print["d ln T_m  = ", fmt[dTDN]];

ratioLaw = FullSimplify[dvDN - dTDN, Assumptions -> $Assumptions];
prodLaw = FullSimplify[dvDN + dTDN, Assumptions -> $Assumptions];
expectZero["ratio law - (2 d ln c_s - d ln a)", ratioLaw - (2*dcs - da)];
expectZero["product law - (dZ + 3 dcsw - drho - 4 da)", prodLaw - (dZ + 3*dcsw - drho - 4*da)];

banner["n=5 wall EOS"];
dvN5 = FullSimplify[dvDN /. dcsw -> 2*drho, Assumptions -> $Assumptions];
dTN5 = FullSimplify[dTDN /. dcsw -> 2*drho, Assumptions -> $Assumptions];
Print["d ln v_w0 = ", fmt[dvN5]];
Print["d ln T_m  = ", fmt[dTN5]];
expectZero[
  "n=5 product law - (dZ + 5 drho - 4 da)",
  (dvN5 + dTN5) - (dZ + 5*drho - 4*da)
];

banner["Stage 164 channel closure"];
(* channelG/channelR are the LHS of eqG/eqR evaluated at the solved (dv,dT); *)
(* they vanish by construction, so this is a solver-consistency print, not an *)
(* independent verification of delta_perp = 0. Reported, not asserted. *)
channelG = FullSimplify[(dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW) /. sol, Assumptions -> $Assumptions];
channelR = FullSimplify[(dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW) /. sol, Assumptions -> $Assumptions];
Print["fixed-g channel (solver consistency) = ", fmt[channelG]];
Print["fixed-r channel (solver consistency) = ", fmt[channelR]];

Print[""];
Print["Carry-forward formulas:"];
Print["  d ln L_W = d ln a"];
Print["  d ln v_w0 = 1/2 d ln(Z_q/rho_w) + 3/2 d ln c_s,w + d ln c_s - d ln a - 3/2 d ln L_W"];
Print["  d ln T_m  = 1/2 d ln(Z_q/rho_w) + 3/2 d ln c_s,w - d ln c_s - d ln a - 1/2 d ln L_W"];
Print["  d ln(v_w0/T_m) = 2 d ln c_s - d ln a"];
Print["  d ln(v_w0 T_m) = d ln Z_q + 3 d ln c_s,w - d ln rho_w - 4 d ln a"];
Print["  n=5 wall EOS: d ln c_s,w = 2 d ln rho_w"];

banner["Lower-branch numeric prefactors"];
rstar = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gstar = ToExpression["0.758035078944663`30"];
tmPref = N[3*Sqrt[10]*3^(3/4)/(5*Pi*gstar*(1 + rstar^2)^(1/4)), 30];
vPref = N[9*Sqrt[10]*3^(1/4)*rstar/(20*(1 + rstar^2)^(3/4)), 30];
ratioPref = N[Sqrt[3]*Pi*gstar*rstar/(4*Sqrt[1 + rstar^2]), 30];
prodPref = N[81*rstar/(10*Pi*gstar*(1 + rstar^2)), 30];
checkPref[name_String, value_, target_] := Module[{diff},
  diff = N[value - target, 40];
  Print[name, " = ", fmt[N[value, 30]]];
  If[TrueQ[Abs[diff] < 10^-12],
    pass[name <> " numeric check"],
    fail[name <> " numeric check", diff]
  ];
];
checkPref["Tm_pref", tmPref, ToExpression["1.2715890393387603`30"]];
checkPref["v_pref", vPref, ToExpression["1.1428896163056477`30"]];
checkPref["ratio_pref", ratioPref, ToExpression["0.8987885086678338`30"]];
checkPref["prod_pref", prodPref, ToExpression["1.4532859092683434`30"]];

Print[""];
Print["Stage 165 Mathematica audit passed."];

Exit[0];
