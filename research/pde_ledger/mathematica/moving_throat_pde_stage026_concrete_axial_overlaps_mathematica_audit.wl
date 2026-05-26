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
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

Clear[s, l, n, lambdaB, lambdaU, lambdaW, lambdaR, varpi, omegaU, omegaW, k, gConst, cSpeed, cs, a, mhat, kEta, tOmega];
$Assumptions = Element[{s, l, lambdaB, lambdaU, lambdaW, lambdaR, varpi, omegaU, omegaW, k, gConst, cSpeed, cs, a, mhat, kEta, tOmega}, Reals] &&
  Element[n, Integers] &&
  l > 0 && n >= 0 && varpi > 0 && omegaU > 0 && omegaW > 0 && gConst > 0 && cSpeed > 0 && cs > 0 && a > 0 && mhat > 0;

concreteModes[] := Module[{u0, fN, f0},
  banner["SECTION I — CONCRETE FINITE-THROAT MODES"];
  u0 = 1/Sqrt[l];
  fN = Sqrt[2/l]*Sin[(n + 1/2)*Pi*s/l];
  f0 = FullSimplify[fN /. n -> 0, Assumptions -> $Assumptions];

  Print["u0(s) = ", fmt[u0]];
  Print["f_n(s) = ", fmt[fN]];
  Print["f0(s) = ", fmt[f0]];

  expectZero["int u0^2 - 1", Integrate[u0^2, {s, 0, l}] - 1];
  expectZero["int f0^2 - 1", Integrate[f0^2, {s, 0, l}] - 1];
  {u0, fN, f0}
];

overlapLaw[] := Module[{u0, fN, f0, indef, kappaN, kappaNViaFundamentalThm, kappa, kappaExpected},
  banner["SECTION II — EXACT OVERLAP LAW"];
  {u0, fN, f0} = concreteModes[];

  (* Independent path: compute the indefinite integral of u0 * fN with respect
     to s, then evaluate at the boundary s = l and s = 0 by the fundamental
     theorem of calculus. This avoids the single Integrate[..., {s,0,l}] call
     used by the SymPy script and exercises a different code path in Mathematica
     (Integrate without limits + boundary evaluation, vs Integrate with limits). *)
  indef = FullSimplify[Integrate[u0*fN, s], Assumptions -> $Assumptions];
  kappaNViaFundamentalThm = FullSimplify[
    (indef /. s -> l) - (indef /. s -> 0),
    Assumptions -> $Assumptions
  ];

  (* Independent path 2: derive kappa_n algebraically from the trigonometric
     antiderivative of Sin[(n+1/2) Pi s/l], which is -Cos[(n+1/2) Pi s/l] *
     l/((n+1/2) Pi). At s = l the cosine is Cos[(n+1/2) Pi] = 0 for integer n;
     at s = 0 it is 1. So the boundary contribution is l/((n+1/2) Pi), and
     multiplying by the mode normalisation Sqrt[2/l] * (1/Sqrt[l]) = Sqrt[2]/l
     yields Sqrt[2]/((n+1/2) Pi). This is the analytic short form. *)
  kappaN = FullSimplify[
    Sqrt[2]/((n + 1/2)*Pi),
    Assumptions -> $Assumptions
  ];

  subbanner["II.1 — General D/N overlap with the constant zero mode"];
  Print["kappa_n (analytic short form) = ", fmt[kappaN]];
  Print["kappa_n (via fundamental thm) = ", fmt[kappaNViaFundamentalThm]];
  expectZero["kappa_n (analytic) - (fundamental thm)", kappaN - kappaNViaFundamentalThm];

  subbanner["II.2 — Lowest-branch overlap"];
  kappa = FullSimplify[kappaN /. n -> 0, Assumptions -> $Assumptions];
  kappaExpected = FullSimplify[2*Sqrt[2]/Pi, Assumptions -> $Assumptions];
  Print["kappa = ", fmt[kappa]];
  expectZero["kappa - 2*sqrt(2)/pi", kappa - kappaExpected];
  Print["kappa numeric = ", fmt[N[kappa, 15]]];
  {u0, f0, kappa}
];

branchSubstitution[] := Module[{u0, f0, kappa, iEtaPhi, iEtaU, iEtaW, iUW, cCoupling, gUeff, gWeff, rEff, delta, q, p, b0, z0, n0, d0, p0},
  banner["SECTION III — CONCRETE BRANCH SUBSTITUTION INTO STAGE-8 QUANTITIES"];
  {u0, f0, kappa} = overlapLaw[];

  overlapU0F0 = FullSimplify[Integrate[u0*f0, {s, 0, l}], Assumptions -> $Assumptions];
  overlapU0U0 = FullSimplify[Integrate[u0*u0, {s, 0, l}], Assumptions -> $Assumptions];

  iEtaPhi = overlapU0F0;
  iEtaU = overlapU0U0;
  iEtaW = overlapU0F0;
  iUW = overlapU0F0;

  subbanner["III.1 — Explicit overlap integrals on the branch"];
  Print["overlap_u0_f0 = ", fmt[overlapU0F0]];
  Print["overlap_u0_u0 = ", fmt[overlapU0U0]];
  Print["I_(eta,phi)   = ", fmt[iEtaPhi]];
  Print["I_(eta,u)     = ", fmt[iEtaU]];
  Print["I_(eta,w)     = ", fmt[iEtaW]];
  Print["I_(u,w)       = ", fmt[iUW]];

  expectZero["overlap_u0_f0 - kappa", overlapU0F0 - kappa];
  expectZero["overlap_u0_u0 - 1", overlapU0U0 - 1];

  cCoupling = FullSimplify[lambdaB*iEtaPhi, Assumptions -> $Assumptions];
  gUeff = FullSimplify[lambdaU*iEtaU, Assumptions -> $Assumptions];
  gWeff = FullSimplify[lambdaW*iEtaW, Assumptions -> $Assumptions];
  rEff = FullSimplify[lambdaR*iUW, Assumptions -> $Assumptions];

  delta = FullSimplify[omegaU^2*omegaW^2 - rEff^2, Assumptions -> $Assumptions];
  q = FullSimplify[gUeff^2*omegaW^2 + 2*gUeff*gWeff*rEff + gWeff^2*omegaU^2, Assumptions -> $Assumptions];
  p = FullSimplify[omegaU^2*gWeff + rEff*gUeff, Assumptions -> $Assumptions];

  b0 = FullSimplify[cCoupling^2/varpi^2, Assumptions -> $Assumptions];
  z0 = FullSimplify[q/delta, Assumptions -> $Assumptions];
  n0 = FullSimplify[p^2/delta^2, Assumptions -> $Assumptions];
  d0 = FullSimplify[k - b0 - z0, Assumptions -> $Assumptions];
  p0 = FullSimplify[n0/d0, Assumptions -> $Assumptions];

  subbanner["III.2 — Reduced coefficients"];
  Print["C     = ", fmt[cCoupling]];
  Print["G_U   = ", fmt[gUeff]];
  Print["G_W   = ", fmt[gWeff]];
  Print["R     = ", fmt[rEff]];
  Print["Delta = ", fmt[delta]];
  Print["Q     = ", fmt[q]];
  Print["P     = ", fmt[p]];
  Print["B0    = ", fmt[b0]];
  Print["Z0    = ", fmt[z0]];
  Print["N0    = ", fmt[n0]];
  Print["D0    = ", fmt[d0]];
  Print["P0    = ", fmt[p0]];

  {kappa, delta, q, p, b0, z0, n0, d0}
];

normalizationTest[] := Module[{kappa, delta, q, p, b0, z0, n0, d0, target, residual, kReq, kGeom},
  banner["SECTION IV — BRANCH-LEVEL NORMALIZATION TEST"];
  {kappa, delta, q, p, b0, z0, n0, d0} = branchSubstitution[];
  target = FullSimplify[54*gConst*cs^5/(5*a^5*cSpeed^5), Assumptions -> $Assumptions];
  residual = FullSimplify[mhat^2*n0/d0 - target, Assumptions -> $Assumptions];

  subbanner["IV.1 — Target residual"];
  Print["Target residual = ", fmt[residual]];

  kReq = k /. First[Solve[residual == 0, k]];
  kReq = FullSimplify[kReq, Assumptions -> $Assumptions];

  subbanner["IV.2 — Exact required wall stiffness"];
  Print["K_req = ", fmt[kReq]];

  (* Independent structural check: the paper's eq:app-stage026-Kreq states the
     three-term decomposition K_req = B0 + Q/Delta + mhat^2 * kappa^2 * (...)^2
     / (target * Delta^2). Verify the solver's output matches that form. *)
  kReqPaper = FullSimplify[
    b0
    + q/delta
    + mhat^2 * kappa^2 * (omegaU^2*lambdaW + lambdaR*lambdaU)^2
      / (target * delta^2),
    Assumptions -> $Assumptions
  ];
  expectZero["K_req - K_req_paper", kReq - kReqPaper];
  expectZero["residual @ K_req", residual /. k -> kReq];

  kGeom = FullSimplify[kEta + 6*tOmega, Assumptions -> $Assumptions];
  Print["On the constant wall branch, bare quadrupole wall stiffness is K = ", fmt[kGeom]];
];

normalizationTest[];

Print[""];
Print["Stage 9 Mathematica audit passed."];

Exit[0];
