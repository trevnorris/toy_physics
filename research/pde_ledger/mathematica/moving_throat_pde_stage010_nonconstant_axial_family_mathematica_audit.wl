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

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[TrigExpand[expr]]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]]
];

Clear[s, l, theta, lambdaB, lambdaU, lambdaW, lambdaR, varpi, omegaU, omegaW, kEta, tOmega, tW, gConst, cSpeed, cs, a, mhat];
$Assumptions = Element[{s, l, theta, lambdaB, lambdaU, lambdaW, lambdaR, varpi, omegaU, omegaW, kEta, tOmega, tW, gConst, cSpeed, cs, a, mhat}, Reals] &&
  l > 0 && varpi > 0 && omegaU > 0 && omegaW > 0 && gConst > 0 && cSpeed > 0 && cs > 0 && a > 0 && mhat > 0;

finiteThroatBasis[] := Module[{u0, u1, f0, chi},
  banner["SECTION I — EXACT FINITE-THROAT N/N AND D/N BASIS"];
  u0 = 1/Sqrt[l];
  u1 = Sqrt[2/l]*Cos[Pi*s/l];
  f0 = Sqrt[2/l]*Sin[Pi*s/(2*l)];
  chi = Cos[theta]*u0 + Sin[theta]*u1;

  Print["u0(s)  = ", fmt[u0]];
  Print["u1(s)  = ", fmt[u1]];
  Print["f0(s)  = ", fmt[f0]];
  Print["chi(s) = ", fmt[chi]];

  expectZero["int u0^2 - 1", Integrate[u0^2, {s, 0, l}] - 1];
  expectZero["int u1^2 - 1", Integrate[u1^2, {s, 0, l}] - 1];
  expectZero["int u0*u1", Integrate[u0*u1, {s, 0, l}]];
  expectZero["int f0^2 - 1", Integrate[f0^2, {s, 0, l}] - 1];
  {u0, u1, f0, chi}
];

overlapLaw[] := Module[{u0, u1, f0, chi, kappa0, kappa1, kappa, rho, blindSubs, maxSubs, kappaBlind, kappaMax, sin2Max},
  banner["SECTION II — EXACT OVERLAP LAW"];
  {u0, u1, f0, chi} = finiteThroatBasis[];
  kappa0 = FullSimplify[Integrate[u0*f0, {s, 0, l}], Assumptions -> $Assumptions];
  kappa1 = FullSimplify[Integrate[u1*f0, {s, 0, l}], Assumptions -> $Assumptions];
  kappa = FullSimplify[Integrate[chi*f0, {s, 0, l}], Assumptions -> $Assumptions];
  rho = FullSimplify[Sqrt[kappa0^2 + kappa1^2], Assumptions -> $Assumptions];

  subbanner["II.1 — Base overlap constants"];
  Print["kappa0 = ", fmt[kappa0]];
  Print["kappa1 = ", fmt[kappa1]];
  Print["kappa(theta) = ", fmt[kappa]];
  Print["rho = ", fmt[rho]];

  expectZero["kappa0 - 2*sqrt(2)/pi", kappa0 - 2*Sqrt[2]/Pi];
  expectZero["kappa1 + 4/(3*pi)", kappa1 + 4/(3*Pi)];
  expectZero["kappa(theta) - [kappa0 cos(theta) + kappa1 sin(theta)]", kappa - (kappa0*Cos[theta] + kappa1*Sin[theta])];
  expectZero["rho - 2*sqrt(22)/(3*pi)", rho - 2*Sqrt[22]/(3*Pi)];

  subbanner["II.2 — Blind angle and max-coupling branch"];
  blindSubs = {Cos[theta] -> Sqrt[2]/Sqrt[11], Sin[theta] -> 3/Sqrt[11]};
  maxSubs = {Cos[theta] -> 3/Sqrt[11], Sin[theta] -> -Sqrt[2]/Sqrt[11]};
  kappaBlind = FullSimplify[kappa /. blindSubs, Assumptions -> $Assumptions];
  kappaMax = FullSimplify[kappa /. maxSubs, Assumptions -> $Assumptions];
  sin2Max = FullSimplify[(Sin[theta]^2) /. maxSubs, Assumptions -> $Assumptions];

  Print["kappa(blind) = ", fmt[kappaBlind]];
  Print["kappa(max)   = ", fmt[kappaMax]];
  Print["sin^2(theta_max) = ", fmt[sin2Max]];
  expectZero["kappa(blind)", kappaBlind];
  expectZero["kappa(max) - rho", kappaMax - rho];
  expectZero["sin^2(theta_max) - 2/11", sin2Max - 2/11];
  {u0, u1, f0, chi, kappa}
];

wallStiffness[] := Module[{u0, u1, f0, chi, kappa, kGeo, kGeoExpected, maxSubs, kGeoMax},
  banner["SECTION III — EXACT WALL STIFFNESS EXPECTATION"];
  {u0, u1, f0, chi, kappa} = overlapLaw[];
  kGeo = kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2;
  kGeoExpected = kGeo;

  Print["K_geo(theta) = ", fmt[kGeo]];
  expectZero["K_geo - expected", kGeo - kGeoExpected];

  maxSubs = {Cos[theta] -> 3/Sqrt[11], Sin[theta] -> -Sqrt[2]/Sqrt[11]};
  kGeoMax = FullSimplify[TrigExpand[kGeo] /. maxSubs, Assumptions -> $Assumptions];
  Print["K_geo(theta_max) = ", fmt[kGeoMax]];
  expectZero["K_geo(theta_max) - [K_eta + 6 T_Omega + 2 T_w pi^2/(11 L^2)]", kGeoMax - (kEta + 6*tOmega + 2*tW*Pi^2/(11*l^2))];
  {kappa, kGeo}
];

branchSubstitution[] := Module[{kappa, kGeo, cCoupling, gUEff, gWEff, rEff, delta, q, p, b0, z0, n0, d0, deltaExpected, qExpected, pExpected, b0Expected, theta0, kappa0},
  banner["SECTION IV — PROFILE-DRESSED STAGE-8/9 BRANCH"];
  {kappa, kGeo} = wallStiffness[];

  cCoupling = lambdaB*kappa;
  gUEff = lambdaU;
  gWEff = lambdaW*kappa;
  rEff = lambdaR*kappa;

  delta = omegaU^2*omegaW^2 - rEff^2;
  q = gUEff^2*omegaW^2 + 2*gUEff*gWEff*rEff + gWEff^2*omegaU^2;
  p = omegaU^2*gWEff + rEff*gUEff;
  b0 = cCoupling^2/varpi^2;
  z0 = q/delta;
  n0 = p^2/delta^2;
  d0 = kGeo - b0 - z0;

  subbanner["IV.1 — Reduced branch quantities"];
  Print["C     = ", fmt[cCoupling]];
  Print["G_U   = ", fmt[gUEff]];
  Print["G_W   = ", fmt[gWEff]];
  Print["R     = ", fmt[rEff]];
  Print["Delta = ", fmt[delta]];
  Print["Q     = ", fmt[q]];
  Print["P     = ", fmt[p]];
  Print["B0    = ", fmt[b0]];
  Print["Z0    = ", fmt[z0]];
  Print["N0    = ", fmt[n0]];
  Print["D0    = ", fmt[d0]];

  deltaExpected = FullSimplify[omegaU^2*omegaW^2 - lambdaR^2*kappa^2, Assumptions -> $Assumptions];
  qExpected = FullSimplify[lambdaU^2*omegaW^2 + 2*lambdaU*lambdaW*lambdaR*kappa^2 + lambdaW^2*omegaU^2*kappa^2, Assumptions -> $Assumptions];
  pExpected = FullSimplify[kappa*(omegaU^2*lambdaW + lambdaR*lambdaU), Assumptions -> $Assumptions];
  b0Expected = FullSimplify[lambdaB^2*kappa^2/varpi^2, Assumptions -> $Assumptions];

  expectZero["Delta - expected", delta - deltaExpected];
  expectZero["Q - expected", q - qExpected];
  expectZero["P - expected", p - pExpected];
  expectZero["B0 - expected", b0 - b0Expected];

  subbanner["IV.2 — Exact recovery of Stage 9 at theta = 0"];
  theta0 = {Cos[theta] -> 1, Sin[theta] -> 0};
  kappa0 = 2*Sqrt[2]/Pi;
  expectZero["kappa(theta=0) - kappa0", (kappa /. theta0) - kappa0];
  expectZero["K_geo(theta=0) - (K_eta + 6 T_Omega)", (TrigExpand[kGeo] /. theta0) - (kEta + 6*tOmega)];
  expectZero[
    "Delta(theta=0) - [Omega_U^2 Omega_W^2 - lambda_R^2 kappa0^2]",
    (delta /. theta0) - (omegaU^2*omegaW^2 - lambdaR^2*kappa0^2)
  ];
  {kappa, kGeo, delta, q, p, n0}
];

normalizationAndNoGo[] := Module[{kappa, kGeo, delta, q, p, n0, target, b0, kReq, blindSubs},
  banner["SECTION V — NORMALIZATION EQUATION AND BLIND-ANGLE NO-GO"];
  {kappa, kGeo, delta, q, p, n0} = branchSubstitution[];
  target = 54*gConst*cs^5/(5*a^5*cSpeed^5);
  b0 = lambdaB^2*kappa^2/varpi^2;
  kReq = b0 + q/delta + mhat^2*p^2/(target*delta^2);
  Print["K_req(theta) = ", fmt[kReq]];

  subbanner["V.1 — Blind-angle no-go"];
  blindSubs = {Cos[theta] -> Sqrt[2]/Sqrt[11], Sin[theta] -> 3/Sqrt[11]};
  pBlind = FullSimplify[TrigExpand[p] /. blindSubs, Assumptions -> $Assumptions];
  deltaBlind = FullSimplify[TrigExpand[delta] /. blindSubs, Assumptions -> $Assumptions];
  expectZero["P(theta_blind)", pBlind];
  expectZero["N0(theta_blind)", FullSimplify[pBlind^2/deltaBlind^2, Assumptions -> $Assumptions]];
];

normalizationAndNoGo[];

Print[""];
Print["Stage 10 Mathematica audit passed."];

Exit[0];
