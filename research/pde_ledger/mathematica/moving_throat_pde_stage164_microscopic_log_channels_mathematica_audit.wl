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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

banner["STAGE 147 — MICROSCOPIC LOG-IMBALANCE CHANNELS"];

Clear[
  ks, kq, lam, gs, gq, tm, js, zq, mu0, lw, qstar, vw0, isq, cs,
  a, ell, mpsi, hbar, rhoW, csw, dlnZ, dlnrho, dlncsw, dlncs,
  dlnTm, dlnv, dlna, dlnLw, rstar, gstar
];

$Assumptions =
  Element[
    {
      ks, kq, lam, gs, gq, tm, js, zq, mu0, lw, qstar, vw0, isq, cs,
      a, ell, mpsi, hbar, rhoW, csw, dlnZ, dlnrho, dlncsw, dlncs,
      dlnTm, dlnv, dlna, dlnLw, rstar, gstar
    },
    Reals
  ] && ks > 0 && kq > 0 && lam > 0 && gs > 0 && gq > 0 && tm > 0 && js > 0 &&
  zq > 0 && mu0 > 0 && lw > 0 && qstar > 0 && vw0 > 0 && isq > 0 && cs > 0 &&
  a > 0 && ell > 0 && mpsi > 0 && hbar > 0 && rhoW > 0 && csw > 0 &&
  rstar > 0 && gstar > 0;

r = lam/Sqrt[ks*kq];
g = gq*Sqrt[ks]/(gs*Sqrt[kq]);
rc = lam^2/(ks*kq);

expectZero["g/r - (g_q K_s)/(g_s lam)", g/r - gq*ks/(gs*lam)];
expectZero["1/r_c - (K_s K_q)/lam^2", 1/rc - ks*kq/lam^2];

banner["General Stage 118 formulas"];

gsGen = tm*js;
gqGen = zq*Pi/(Sqrt[2]*mu0*lw^(3/2));
kqGen = zq*Pi^2*cs^2/(4*mu0*lw^2);
lamGen = -qstar*vw0*isq;

firstProdGeneral = FullSimplify[gqGen*ks/(gsGen*lamGen), Assumptions -> $Assumptions];
secondProdGeneral = FullSimplify[ks*kqGen/lamGen^2, Assumptions -> $Assumptions];

Print["general first product = ", fmt[firstProdGeneral]];
Print["general second product = ", fmt[secondProdGeneral]];

banner["Uniform-overlap + D/N simplification"];

jsClosure = 4*Pi*a^2*ell/3;
iq = 2*Sqrt[2*lw]/Pi;
isqClosure = jsClosure*iq;

firstProdUniform = FullSimplify[
  firstProdGeneral /. {js -> jsClosure, isq -> isqClosure},
  Assumptions -> $Assumptions
];
secondProdUniform = FullSimplify[
  secondProdGeneral /. {isq -> isqClosure},
  Assumptions -> $Assumptions
];

firstUniformExpected = -9*ks*zq/(64*lw^2*tm*a^4*ell^2*mu0*qstar*vw0);
secondUniformExpected = 9*Pi^2*ks*zq*cs^2/(512*lw^3*a^4*ell^2*mu0*qstar^2*vw0^2);

Print["uniform-overlap first product = ", fmt[firstProdUniform]];
Print["uniform-overlap second product = ", fmt[secondProdUniform]];
expectZero["uniform first product", firstProdUniform - firstUniformExpected];
expectZero["uniform second product", secondProdUniform - secondUniformExpected];

banner["Healing-locked shell branch"];

ellLock = hbar/(2*mpsi*csw);
ksLock = FullSimplify[3*Pi*a^2*hbar^2/(5*mpsi*rhoW*ellLock), Assumptions -> $Assumptions];

firstProdHeal = FullSimplify[
  firstProdUniform /. {ks -> ksLock, ell -> ellLock},
  Assumptions -> $Assumptions
];
secondProdHeal = FullSimplify[
  secondProdUniform /. {ks -> ksLock, ell -> ellLock},
  Assumptions -> $Assumptions
];

firstHealExpected = -(27/40)*Pi*mpsi^2*zq*csw^3/(hbar*mu0*qstar*rhoW*tm*vw0*a^2*lw^2);
secondHealExpected = (27/320)*Pi^3*mpsi^2*zq*cs^2*csw^3/(hbar*mu0*qstar^2*rhoW*vw0^2*a^2*lw^3);

Print["healing-locked first product = ", fmt[firstProdHeal]];
Print["healing-locked second product = ", fmt[secondProdHeal]];
expectZero["healing first product exact formula", firstProdHeal - firstHealExpected];
expectZero["healing second product exact formula", secondProdHeal - secondHealExpected];

banner["delta_perp on the healing-locked branch"];

b = 1/(4*Sqrt[1 + rstar^2]);
firstHeal = dlnZ + 3*dlncsw - dlnrho - dlnTm - dlnv - 2*dlna - 2*dlnLw;
secondHeal = dlnZ + 2*dlncs + 3*dlncsw - dlnrho - 2*dlnv - 2*dlna - 3*dlnLw;
deltaPerp = Expand[gstar*firstHeal + b*secondHeal];
deltaPerpExpected = Expand[
  (gstar + b)*(dlnZ - dlnrho)
  + 3*(gstar + b)*dlncsw
  + 2*b*dlncs
  - gstar*dlnTm
  - (gstar + 2*b)*dlnv
  - 2*(gstar + b)*dlna
  - (2*gstar + 3*b)*dlnLw
];
expectZero["delta_perp compressed form", deltaPerp - deltaPerpExpected];
Print["delta_perp = ", fmt[deltaPerpExpected]];

dlnTmSol = FullSimplify[
  dlnTm /. First[Solve[deltaPerpExpected == 0, dlnTm, Reals]],
  Assumptions -> $Assumptions
];
Print["tangency-law dln T_m = ", fmt[dlnTmSol]];

banner["Family-1 numerical coefficients"];

gNum = 0.758035078944663;
rNum = 1.77799353547498;
bNum = N[1/(4*Sqrt[1 + rNum^2]), 20];
aNum = N[gNum + bNum, 20];
bbNum = N[2*bNum, 20];

Print["g_* = ", N[gNum, 20]];
Print["1/(4 sqrt(1+r_*^2)) = ", bNum];
Print["A_* = g_* + 1/(4 sqrt(...)) = ", aNum];
Print["B_* = 1/(2 sqrt(1+r_*^2)) = ", bbNum];
Print["coeff[ln(Z_q/rho_w)] = ", aNum];
Print["coeff[ln(c_sw)]      = ", N[3*aNum, 20]];
Print["coeff[ln(c_s)]       = ", bbNum];
Print["coeff[ln(T_m)]       = ", N[-gNum, 20]];
Print["coeff[ln(v_w0)]      = ", N[-(gNum + bbNum), 20]];
Print["coeff[ln(a)]         = ", N[-2*aNum, 20]];
Print["coeff[ln(L_W)]       = ", N[-(2*gNum + 3*bNum), 20]];

expectApprox["A_* numeric check", aNum, 0.88058909004156951300, 10^-14];
expectApprox["B_* numeric check", bbNum, 0.24510802219381302601, 10^-14];
expectApprox["c_sw coefficient numeric check", 3*aNum, 2.6417672701247085390, 10^-14];
expectApprox["v_w0 coefficient numeric check", -(gNum + bbNum), -1.0031431011384760260, 10^-14];
expectApprox["L_W coefficient numeric check", -(2*gNum + 3*bNum), -1.8837321911800455390, 10^-14];

Print[""];
Print["Carry-forward formulas:"];
Print["  dln[(g_q K_s)/(g_s lam)] = dln(g/r)"];
Print["  dln[(K_s K_q)/lam^2] = - dln(r_c)"];
Print["  On the healing-locked branch the two channels reduce to direct microscopic drifts."];
Print["  delta_perp is the weighted sum of those two channels."];

Exit[0];
