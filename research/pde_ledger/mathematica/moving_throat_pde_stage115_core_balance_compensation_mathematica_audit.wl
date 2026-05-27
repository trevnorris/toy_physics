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

banner["STAGE 115 — EXACT CORE-BALANCE COMPENSATION THEOREM"];

Clear[kS, kQ, lam, gS, gQ, kappa0, gamma0, z];
$Assumptions =
  Element[{kS, kQ, lam, gS, gQ, kappa0, gamma0, z}, Reals] &&
  kS > 0 && kQ > 0 && kappa0 > 0 && gamma0 > 0;

rC = FullSimplify[lam^2/(kS*kQ), Assumptions -> $Assumptions];
rhoC = FullSimplify[gS^2/kS, Assumptions -> $Assumptions];
sigmaC = FullSimplify[(kS*gQ - lam*gS)^2/(kS^2*kQ*(1 + rC)), Assumptions -> $Assumptions];

balanceEq = Expand[rhoC - 4*sigmaC];
Print["rho_c - 4 sigma_c = ", fmt[Factor[balanceEq]]];

gQSolutions = Solve[balanceEq == 0, gQ, Reals];
If[Length[gQSolutions] =!= 2, fail["expected two coupling-balance branches", gQSolutions]];
Print["Exact coupling-balance solutions for g_q = ", fmt[gQ /. gQSolutions]];

gQBranch = FullSimplify[gQ /. First[gQSolutions], Assumptions -> $Assumptions];
sigmaStar = FullSimplify[gS^2/(4*kS), Assumptions -> $Assumptions];
expectZero["sigma_c on balance surface", (sigmaC /. gQ -> gQBranch) - sigmaStar];

(* Independent derivation: parent-overlap reparametrization (Part-IV appendix
   eqs eq:app-part04-r-g-parent-ratios and eq:app-part04-parent-compensation-family). *)
frakR = lam/Sqrt[kS*kQ];
frakG = gQ*Sqrt[kS]/(gS*Sqrt[kQ]);
parentFamilyResidual = 1 + frakR^2 - 4*(frakG - frakR)^2;
expectZero[
  "independent: parent family identical balance equation",
  FullSimplify[
    parentFamilyResidual - balanceEq*(kS*kQ + lam^2)/(gS^2*kQ),
    Assumptions -> $Assumptions
  ]
];
(* Solve the parent family for a fresh variable, then translate back. *)
Clear[gVar];
parentResidualGen = 1 + frakR^2 - 4*(gVar - frakR)^2;
frakGRoots = gVar /. Solve[parentResidualGen == 0, gVar];
If[Length[frakGRoots] =!= 2,
  fail["expected two parent-family roots", frakGRoots]];
Print["Parent-family roots for frakG = ", fmt[frakGRoots]];
frakGMinus = FullSimplify[frakR - Sqrt[1 + frakR^2]/2,
  Assumptions -> $Assumptions];
expectZero[
  "independent: frakG_- root matches Solve output",
  FullSimplify[(frakGRoots[[1]]) - frakGMinus,
    Assumptions -> $Assumptions] *
  FullSimplify[(frakGRoots[[2]]) - frakGMinus,
    Assumptions -> $Assumptions]
];
gQFromFrakMinus = FullSimplify[
  frakGMinus*gS*Sqrt[kQ]/Sqrt[kS],
  Assumptions -> $Assumptions
];
expectZero[
  "independent: sigma_c = sigma_* via parent reparametrization",
  FullSimplify[
    (sigmaC /. gQ -> gQFromFrakMinus) - sigmaStar,
    Assumptions -> $Assumptions
  ]
];

kappa0Can = FullSimplify[(1 + rC)/3, Assumptions -> $Assumptions];
gamma0Can = FullSimplify[(1 + rC)/9, Assumptions -> $Assumptions];

deltaCore = FullSimplify[
  (rhoC - sigmaC/(1 - (kappa0/(1 + rC))*z^2 - I*(gamma0/(1 + rC))*z^5)) /.
    {gQ -> gQBranch, kappa0 -> kappa0Can, gamma0 -> gamma0Can},
  Assumptions -> $Assumptions
];
targetDelta = FullSimplify[4*sigmaStar - sigmaStar/(1 - z^2/3 - I*z^5/9), Assumptions -> $Assumptions];
expectZero["exact collapse to compensated branch", deltaCore - targetDelta];

lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9;
lambdaEff = Expand[Normal[Series[lambdaOut + targetDelta, {z, 0, 5}]]];
yEff = FullSimplify[(lambdaEff /. z -> 0)/lambdaEff, Assumptions -> $Assumptions];
yTarget = 1 + z^2/9 + 4*z^4/81 + I*z^5/27;
expectZero["normalized outgoing fingerprint preserved", Expand[Normal[Series[yEff, {z, 0, 5}]]] - yTarget];

Print[""];
Print["Derived branch data:"];
Print["sigma_* = ", fmt[sigmaStar]];
Print["kappa0  = ", fmt[kappa0Can]];
Print["gamma0  = ", fmt[gamma0Can]];

Print[""];
Print["Stage 115 Mathematica audit passed."];

Exit[0];
