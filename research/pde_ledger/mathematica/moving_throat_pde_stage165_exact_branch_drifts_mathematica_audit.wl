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

eqR = dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW == 0;
eqG = dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW == 0;
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

Print[""];
Print["Stage 165 Mathematica audit passed."];

Exit[0];
