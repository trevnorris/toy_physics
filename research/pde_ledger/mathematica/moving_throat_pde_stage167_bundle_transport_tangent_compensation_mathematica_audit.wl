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

banner["STAGE 150 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION"];

Clear[dTheta, dKs, dKq, dP, drho, da, dcs, dZ, dcsw, dell, dLW, dv, dT, dgq, dgs, dIsq, dlam, drc, dr, dg, gstar, rstar];
$Assumptions = Element[{dTheta, dKs, dKq, dP, drho, da, dcs, dZ, dcsw, dell, dLW, dv, dT, dgq, dgs, dIsq, dlam, drc, dr, dg, gstar, rstar}, Reals] && gstar > 0 && rstar > 0;

drho = dTheta/2;
da = dKs/2 - dTheta/4;
dcs = dKs/2 - dTheta/4 + dP/5;
dZ = dKq - 2*dP/5;
dcsw = 2*drho;
dell = -dcsw;
dLW = da;

Print["delta ln rho_w  = ", fmt[drho]];
Print["delta ln a      = ", fmt[da]];
Print["delta ln c_s    = ", fmt[dcs]];
Print["delta ln Z_q    = ", fmt[dZ]];
Print["delta ln c_s,w  = ", fmt[dcsw]];
Print["delta ln ell    = ", fmt[dell]];
Print["delta ln L_W    = ", fmt[dLW]];

dv = FullSimplify[(dZ - drho)/2 + 3*dcsw/2 + dcs - 5*da/2, Assumptions -> $Assumptions];
dT = FullSimplify[(dZ - drho)/2 + 3*dcsw/2 - dcs - 3*da/2, Assumptions -> $Assumptions];

banner["Branch transport"];
Print["delta ln v_w0 = ", fmt[dv]];
Print["delta ln T_m  = ", fmt[dT]];
expectZero["delta ln v_w0 carry-forward law", dv - (-3*dKs/4 + dKq/2 + 13*dTheta/8)];
expectZero["delta ln T_m carry-forward law", dT - (-5*dKs/4 + dKq/2 + 15*dTheta/8 - 2*dP/5)];
expectZero["delta ln(v/T) - (2 d ln c_s - d ln a)", (dv - dT) - (2*dcs - da)];
expectZero["delta ln(v*T) - (d ln Z_q + 5 d ln rho_w - 4 d ln a)", (dv + dT) - (dZ + 5*drho - 4*da)];

dgq = FullSimplify[dZ - 3*dLW/2, Assumptions -> $Assumptions];
dgs = FullSimplify[dT + 2*da + dell, Assumptions -> $Assumptions];
dIsq = FullSimplify[2*da + dell + dLW/2, Assumptions -> $Assumptions];
dlam = FullSimplify[dv + dIsq, Assumptions -> $Assumptions];

banner["Parent couplings"];
Print["delta ln g_q    = ", fmt[dgq]];
Print["delta ln g_s    = ", fmt[dgs]];
Print["delta ln I_sq   = ", fmt[dIsq]];
Print["delta ln lambda = ", fmt[dlam]];
expectZero["delta ln g_s carry-forward law", dgs - (-dKs/4 + dKq/2 + 3*dTheta/8 - 2*dP/5)];
expectZero["delta ln g_q carry-forward law", dgq - (-3*dKs/4 + dKq + 3*dTheta/8 - 2*dP/5)];
expectZero["delta ln lambda carry-forward law", dlam - (dKs + dKq)/2];

banner["Parent invariants"];
drc = FullSimplify[2*dlam - dKs - dKq, Assumptions -> $Assumptions];
dr = FullSimplify[dlam - (dKs + dKq)/2, Assumptions -> $Assumptions];
dg = FullSimplify[dgq + dKs/2 - dgs - dKq/2, Assumptions -> $Assumptions];
expectZero["delta ln r_c", drc];
expectZero["delta ln frak r", dr];
expectZero["delta ln frak g", dg];

banner["Stage 163 off-family channels"];
chan1 = FullSimplify[dgq + dKs - dgs - dlam, Assumptions -> $Assumptions];
chan2 = FullSimplify[dKs + dKq - 2*dlam, Assumptions -> $Assumptions];
expectZero["delta ln(g_q K_s/(g_s lambda))", chan1];
expectZero["delta ln(K_s K_q/lambda^2)", chan2];

deltaPerp = FullSimplify[gstar*chan1 + chan2/(4*Sqrt[1 + rstar^2]), Assumptions -> $Assumptions];
expectZero["delta_perp", deltaPerp];

Print[""];
Print["Carry-forward formulas:"];
Print["  delta ln v_w0   = -3/4 delta ln K_s + 1/2 delta ln K_q + 13/8 delta ln Theta_w"];
Print["  delta ln T_m    = -5/4 delta ln K_s + 1/2 delta ln K_q + 15/8 delta ln Theta_w - 2/5 delta ln P_0"];
Print["  delta ln g_s    = -1/4 delta ln K_s + 1/2 delta ln K_q + 3/8 delta ln Theta_w - 2/5 delta ln P_0"];
Print["  delta ln g_q    = -3/4 delta ln K_s + delta ln K_q + 3/8 delta ln Theta_w - 2/5 delta ln P_0"];
Print["  delta ln lambda = 1/2 (delta ln K_s + delta ln K_q)"];
Print["  delta ln r_c = delta ln frak r = delta ln frak g = 0"];
Print["  delta_perp = 0"];

Print[""];
Print["Stage 167 Mathematica audit passed."];

Exit[0];
