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

banner["STAGE 107 — GENERAL ISOTROPIC DTN DEFORMATION ALGEBRA"];

Clear[z, sNorm, beta, sigma0, sigma2, sigma4, sigma5, y2, y4, y5, ell0, ell2, ell4, ell5];
$Assumptions =
  Element[{z, sNorm, beta, sigma0, sigma2, sigma4, sigma5}, Reals] &&
  sNorm != 0 && beta != 0 && 3*sNorm - sigma0 != 0;

lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9;
lambdaDef = Expand[sNorm*(lambdaOut /. z -> beta*z) + sigma0 + sigma2*z^2 + sigma4*z^4 + I*sigma5*z^5];

l0 = FullSimplify[lambdaDef /. z -> 0, Assumptions -> $Assumptions];
l2 = FullSimplify[Coefficient[lambdaDef, z, 2], Assumptions -> $Assumptions];
l4 = FullSimplify[Coefficient[lambdaDef, z, 4], Assumptions -> $Assumptions];
l5 = FullSimplify[Coefficient[lambdaDef, z, 5]/I, Assumptions -> $Assumptions];

branchJet = ell0 + ell2*z^2 + ell4*z^4 + I*ell5*z^5;
branchAnsatz = 1 + y2*z^2 + y4*z^4 + I*y5*z^5;
branchResidual = Expand[branchJet*branchAnsatz - ell0];
branchEquations =
  Thread[(CoefficientList[branchResidual, z][[# + 1]] & /@ {2, 4, 5}) == 0];
branchSol = Solve[branchEquations, {y2, y4, y5}];
If[Length[branchSol] =!= 1, fail["unique branch-coefficient solution", branchSol]];
branchSol = First[branchSol];
branchSubstitution = {ell0 -> l0, ell2 -> l2, ell4 -> l4, ell5 -> l5};
branchY2 = y2 /. branchSol;
branchY4 = branchY2^2 - ell4/ell0;
branchY5 = y5 /. branchSol;
yDirect = Expand[branchAnsatz /. branchSol /. branchSubstitution];
yFormula = 1 - (l2/l0)*z^2 + (l2^2/l0^2 - l4/l0)*z^4 - I*(l5/l0)*z^5;

Print["Lambda_out(z) = ", fmt[lambdaOut]];
Print["Lambda_def(z) = ", fmt[lambdaDef]];
Print["L0 = ", fmt[l0]];
Print["L2 = ", fmt[l2]];
Print["L4 = ", fmt[l4]];
Print["L5 = ", fmt[l5]];
Print["Y_def(z) = ", fmt[Expand[yFormula]]];

expectZero["normalized expansion direct-formula", yDirect - yFormula];

m2 = FullSimplify[(branchY2 /. branchSubstitution), Assumptions -> $Assumptions];
m4 = FullSimplify[(branchY4 /. branchSubstitution), Assumptions -> $Assumptions];
chiQ = FullSimplify[((branchY5 /. branchSubstitution)/(1/27)), Assumptions -> $Assumptions];

Print["m2 = ", fmt[m2]];
Print["m4 = ", fmt[m4]];
Print["chi_Q = ", fmt[chiQ]];

sol = Solve[{m2 == 1/9, m4 == 4/81}, {sigma2, sigma4}, Reals];
If[Length[sol] =!= 1, fail["unique (sigma2, sigma4) solution", sol]];
sol = First[sol];

Print["Sigma2_evenmatch = ", fmt[sigma2 /. sol]];
Print["Sigma4_evenmatch = ", fmt[sigma4 /. sol]];
expectZero["Sigma2 exact formula", (sigma2 /. sol) + (3*sNorm*beta^2 - 3*sNorm + sigma0)/9];
expectZero["Sigma4 exact formula", (sigma4 /. sol) + (3*sNorm*beta^4 - 3*sNorm + sigma0)/27];

chiEven = FullSimplify[chiQ /. sol, Assumptions -> $Assumptions];
Print["chi_Q under canonical-even matching = ", fmt[Factor[chiEven]]];
expectZero[
  "chi_Q - 3(sNorm beta^5 + 9 sigma5)/(3 sNorm - sigma0)",
  chiEven - 3*(sNorm*beta^5 + 9*sigma5)/(3*sNorm - sigma0)
];

Print[""];
Print["Stage 107 Mathematica audit passed."];

Exit[0];
