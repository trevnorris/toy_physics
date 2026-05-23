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

banner["STAGE 053 — EXPLICIT GNLS WALL-SHELL REDUCTION"];

Clear[a, L, ell, rhoW, cSw, V0, m, hbar, IfMoment, Ig, Hw];
$Assumptions =
  Element[{a, L, ell, rhoW, cSw, V0, m, hbar, IfMoment, Ig, Hw}, Reals] &&
  a > 0 && L > 0 && ell > 0 && rhoW > 0 && cSw > 0 && V0 > 0 && m > 0 && hbar > 0 && IfMoment > 0 && Ig > 0;

kappaClosed = FullSimplify[
  4*(m*cSw*L/hbar)^2 + (Ig/IfMoment)*(L/ell)^2,
  Assumptions -> $Assumptions
];
WwallClosed = FullSimplify[
  4*rhoW^2*V0^2*L^2/(hbar^2*cSw^2*ell^2),
  Assumptions -> $Assumptions
];

Hw = FullSimplify[m*cSw^2/rhoW, Assumptions -> $Assumptions];
Print["H_w = ", fmt[Hw]];

Nphiphi = FullSimplify[4*Pi*a^2*ell*IfMoment, Assumptions -> $Assumptions];
Gphiphi = FullSimplify[4*Pi*a^2*Ig/ell, Assumptions -> $Assumptions];
Tx = FullSimplify[hbar^2*Nphiphi/(4*m*rhoW), Assumptions -> $Assumptions];
Kx = FullSimplify[Hw*Nphiphi + hbar^2*Gphiphi/(4*m*rhoW), Assumptions -> $Assumptions];
kappaAssembled = FullSimplify[Kx*L^2/Tx, Assumptions -> $Assumptions];

Print["N_(phi phi) = ", fmt[Nphiphi]];
Print["G_(phi phi) = ", fmt[Gphiphi]];
Print["T_X = ", fmt[Tx]];
Print["K_X = ", fmt[Kx]];
Print["kappa = ", fmt[kappaAssembled]];

expectZero["kappa - expected", kappaAssembled - kappaClosed];

WwallAssembled = FullSimplify[
  4*Pi*a^2*L^2*(IfMoment*rhoW/(m*cSw^2))*V0^2/(Tx*ell),
  Assumptions -> $Assumptions
];

Print["W_wall = ", fmt[WwallAssembled]];
expectZero["W_wall - expected", WwallAssembled - WwallClosed];

Xi = FullSimplify[
  (V0/ell)^2 * (4*Pi*a^2*ell*IfMoment*rhoW/(m*cSw^2)) * L^2 / Tx,
  Assumptions -> $Assumptions
];

Print["Xi = ", fmt[Xi]];
expectZero["Xi - W_wall", Xi - WwallAssembled];

banner["STAGE 053 THEOREM LEDGER"];
Print["T_X = pi a^2 ell I_f hbar^2 / (m rho_w)"];
Print["K_X = 4 pi a^2 ell I_f H_w + pi a^2 I_g hbar^2 / (m rho_w ell)"];
Print["kappa = 4 (m c_(s,w) L / hbar)^2 + (I_g/I_f)(L/ell)^2"];
Print["W_wall = Xi = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2)"];

Exit[0];
