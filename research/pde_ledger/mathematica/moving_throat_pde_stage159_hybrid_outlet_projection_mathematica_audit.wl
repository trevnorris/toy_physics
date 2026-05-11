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

banner["STAGE 142 — HYBRID OUTLET LINEARIZATION"];

Clear[rho, sigma, kappa, gamma, sigma0, drho, dsigma, dkappa, dgamma];
$Assumptions = Element[{rho, sigma, kappa, gamma, sigma0, drho, dsigma, dkappa, dgamma}, Reals];

l0 = -3 + rho - sigma;
l2 = 1/3 - sigma*kappa;
l4 = 1/9 - sigma*kappa^2;
chi = FullSimplify[3*(1 - 9*sigma*gamma)/(3 - rho + sigma)];
e2 = FullSimplify[-l2/l0 - 1/9];
e4 = FullSimplify[l2^2/l0^2 - l4/l0 - 4/81];

subs = {rho -> 4*sigma0 + drho, sigma -> sigma0 + dsigma, kappa -> 1/3 + dkappa, gamma -> 1/9 + dgamma};
vars = {drho, dsigma, dkappa, dgamma};
linearize[expr_] := Module[{f = expr /. subs},
  FullSimplify[
    (f /. Thread[vars -> 0]) +
      Sum[(D[f, vars[[i]]] /. Thread[vars -> 0]) * vars[[i]], {i, Length[vars]}]
  ]
];

chiLin = linearize[chi];
e2Lin = linearize[e2];
e4Lin = linearize[e4];

deltaChi = Expand[chiLin - 1];
deltaChiExpected = Expand[(drho - 4*dsigma - 27*sigma0*dgamma)/(3*(1 - sigma0))];
e2Expected = Expand[(drho - 4*dsigma - 9*sigma0*dkappa)/(27*(1 - sigma0))];
e4Expected = Expand[(5*drho - 20*dsigma - 72*sigma0*dkappa)/(243*(1 - sigma0))];

expectZero["delta chi formula", deltaChi - deltaChiExpected];
expectZero["delta E2 formula", e2Lin - e2Expected];
expectZero["delta E4 formula", e4Lin - e4Expected];

banner["MOUTH-GAIN -> HYBRID-LOADING TRANSPORT"];

Clear[xi, sigma0Can, dSigma0, dR];
$Assumptions = Element[{xi, sigma0Can, dSigma0, dR}, Reals];

dMs = dSigma0;
dMq = -1/4*dSigma0 - sigma0Can*dR;
drhoExpr = xi*dMs;
dsigmaExpr = -xi*dMq;
deltaCExpr = Expand[drhoExpr - 4*dsigmaExpr];
deltaCExpected = Expand[-4*xi*sigma0Can*dR];
expectZero["deltaC mouth transport", deltaCExpr - deltaCExpected];

Clear[sigmaStar];
$Assumptions = Element[{xi, sigma0Can, dSigma0, dR, sigmaStar}, Reals];
expectZero[
  "sigma_star substitution",
  (deltaCExpected /. xi -> 4*sigmaStar/sigma0Can) + 16*sigmaStar*dR
];

banner["CANONICAL-EVEN PRESERVATION"];

Clear[deltaC, dk];
$Assumptions = Element[{deltaC, dk, sigma0}, Reals];
eq1 = deltaC - 9*sigma0*dk;
eq2 = 5*deltaC - 72*sigma0*dk;
sol = Solve[{eq1 == 0, eq2 == 0}, {deltaC, dk}, Reals];
Print["solution = ", fmt[sol]];
If[sol =!= {{deltaC -> 0, dk -> 0}}, fail["canonical-even preservation solution", sol]];
det = Det[{{1, -9*sigma0}, {5, -72*sigma0}}];
Print["determinant = ", fmt[Factor[det]]];

banner["FINAL REDUCED DEFECT"];
finalDefect = FullSimplify[((deltaC - 27*sigma0*dgamma)/(3*(1 - sigma0))) /. deltaC -> 0];
expectZero["final Delta_Q + 9 sigma* dgamma /(1-sigma*)", finalDefect + 9*sigma0*dgamma/(1 - sigma0)];

Print[""];
Print["Carry-forward formulas:"];
Print["  delta C := delta rho_R - 4 delta sigma_W"];
Print["  delta E2 = (delta C - 9 sigma_* delta kappa_W)/(27(1-sigma_*))"];
Print["  delta E4 = (5 delta C - 72 sigma_* delta kappa_W)/(243(1-sigma_*))"];
Print["  Delta_Q  = (delta C - 27 sigma_* delta gamma_W)/(3(1-sigma_*))"];
Print["  Canonical-even preservation => delta C = delta kappa_W = 0"];
Print["  Hence Delta_Q = -9 sigma_* delta gamma_W/(1-sigma_*)"];

Print[""];
Print["Stage 159 Mathematica audit passed."];

Exit[0];
