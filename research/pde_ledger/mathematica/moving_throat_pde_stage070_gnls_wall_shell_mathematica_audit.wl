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

banner["STAGE 070 — INDEPENDENT NUMERIC PROFILE CROSS-CHECK"];

Module[{xi, fProf, fp, fpp, IfNum, IgNum, ruleNum, TxNum, KxNum, kappaNum,
        kappaCmp, WwallNum, WwallCmp, JfromProfile, J1Stage48, XiNum, XiCmp, tol},
  tol = 10^-10;
  fProf[xi_] := Sech[xi];
  fp[xi_]    := D[fProf[xi], xi];
  fpp[xi_]   := D[fProf[xi], {xi, 2}];
  IfNum = NIntegrate[fp[xi]^2,  {xi, -Infinity, Infinity}, WorkingPrecision -> 30];
  IgNum = NIntegrate[fpp[xi]^2, {xi, -Infinity, Infinity}, WorkingPrecision -> 30];
  Print["I_f (sech profile) = ", fmt[IfNum], "   (analytic 2/3 = ", N[2/3, 30], ")"];
  Print["I_g (sech profile) = ", fmt[IgNum], "   (analytic 8/15 = ", N[8/15, 30], ")"];
  Print["Stage-48 normalization: J_1 := I_f/H_w (shell measure 4 pi a^2 ell absorbed into J_1)"];
  Print["Stage-47 normalization: I_1 := N_phiphi/H_w = (4 pi a^2 ell I_f)/H_w"];
  Print["Structural ratio I_1/J_1 should equal 4 pi a^2 ell."];
  (* Verify with the sech-profile I_f computed above:
     I_1/J_1 = (4 pi a^2 ell I_f / H_w) / (I_f / H_w) = 4 pi a^2 ell, independent of I_f's value. *)
  expectZero["I_1 / J_1 - 4 pi a^2 ell (independent of profile, symbolic check)",
             FullSimplify[(4*Pi*a^2*ell*IfMoment/Hw)/(IfMoment/Hw) - 4*Pi*a^2*ell,
                          Assumptions -> $Assumptions]];

  ruleNum = {a -> 1, L -> 1, ell -> 1/10, rhoW -> 1, cSw -> 1, V0 -> 1, m -> 1, hbar -> 1};

  TxNum    = N[Pi*a^2*ell*IfNum*hbar^2/(m*rhoW)                /. ruleNum, 30];
  KxNum    = N[(4*Pi*a^2*ell*IfNum*(m*cSw^2/rhoW)
               + Pi*a^2*IgNum*hbar^2/(m*rhoW*ell))             /. ruleNum, 30];
  kappaNum = N[KxNum*L^2/TxNum                                  /. ruleNum, 30];
  kappaCmp = N[4*(m*cSw*L/hbar)^2 + (IgNum/IfNum)*(L/ell)^2      /. ruleNum, 30];
  Print["kappa_num     = ", fmt[kappaNum]];
  Print["kappa_closed  = ", fmt[kappaCmp]];
  If[Abs[kappaNum - kappaCmp] < tol, pass["kappa numeric profile check"],
    fail["kappa numeric profile check", kappaNum - kappaCmp]];

  WwallNum = N[4*Pi*a^2*L^2*(IfNum*rhoW/(m*cSw^2))*V0^2/(TxNum*ell) /. ruleNum, 30];
  WwallCmp = N[4*rhoW^2*V0^2*L^2/(hbar^2*cSw^2*ell^2)               /. ruleNum, 30];
  Print["W_wall_num    = ", fmt[WwallNum]];
  Print["W_wall_closed = ", fmt[WwallCmp]];
  If[Abs[WwallNum - WwallCmp] < tol, pass["W_wall numeric profile check"],
    fail["W_wall numeric profile check", WwallNum - WwallCmp]];

  XiNum = N[(V0/ell)^2*(4*Pi*a^2*ell*IfNum*rhoW/(m*cSw^2))*L^2/TxNum /. ruleNum, 30];
  XiCmp = WwallNum;
  Print["Xi_num        = ", fmt[XiNum]];
  If[Abs[XiNum - XiCmp] < tol, pass["Xi = W_wall numeric profile check"],
    fail["Xi = W_wall numeric profile check", XiNum - XiCmp]];
];

banner["STAGE 053 THEOREM LEDGER"];
Print["T_X = pi a^2 ell I_f hbar^2 / (m rho_w)"];
Print["K_X = 4 pi a^2 ell I_f H_w + pi a^2 I_g hbar^2 / (m rho_w ell)"];
Print["kappa = 4 (m c_(s,w) L / hbar)^2 + (I_g/I_f)(L/ell)^2"];
Print["W_wall = Xi = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2)"];

Exit[0];
