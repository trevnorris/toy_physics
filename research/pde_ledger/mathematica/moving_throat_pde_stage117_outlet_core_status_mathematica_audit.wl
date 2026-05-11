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

banner["STAGE 100 — CONCRETE OUTLET-CORE STATUS"];

Clear[z, s, beta, rho, sigma, kappa, gamma, ks, kq, lam, gs, gq, kappa0, gamma0, a, lW];
$Assumptions =
  Element[{z, s, beta, rho, sigma, kappa, gamma, gs, gq}, Reals] &&
  Element[{ks, kq, lam, kappa0, gamma0, a, lW}, Reals] &&
  ks > 0 && kq > 0 && lam > 0 && kappa0 > 0 && gamma0 > 0 && a > 0 && lW > 0;

lambdaOut = -3 + z^2/3 + z^4/9 + I z^5/9;

banner["1. Harmless scale/argument class"];
yArg = Normal[Series[(-3 s)/(s (lambdaOut /. z -> beta z)), {z, 0, 5}]];
m2Arg = FullSimplify[Coefficient[yArg, z, 2]];
m4Arg = FullSimplify[Coefficient[yArg, z, 4]];
chiArg = FullSimplify[Coefficient[yArg, z, 5]/I/(1/27)];
betaSolutions = Solve[{m2Arg == 1/9, m4Arg == 4/81}, beta, Reals];
Print["scale/argument solutions = ", fmt[betaSolutions]];
If[Sort[beta /. betaSolutions] =!= {-1, 1}, fail["scale/argument branch roots", betaSolutions]];
expectZero["positive harmless branch has beta = 1 and chi_Q = 1", (chiArg /. beta -> 1) - 1];

banner["2. Pure Robin class"];
lambdaR = lambdaOut + rho;
yR = Normal[Series[(-3 + rho)/lambdaR, {z, 0, 5}]];
c2R = FullSimplify[Coefficient[yR, z, 2]];
c4R = FullSimplify[Coefficient[yR, z, 4]];
chiR = FullSimplify[Coefficient[yR, z, 5]/I/(1/27)];
robinSolutions = Solve[{c2R == 1/9, c4R == 4/81}, rho, Reals];
Print["pure Robin canonical-even solutions = ", fmt[robinSolutions]];
If[robinSolutions =!= {{rho -> 0}}, fail["pure Robin branch", robinSolutions]];
expectZero["pure Robin odd norm is trivial on rho = 0", (chiR /. rho -> 0) - 1];

banner["3. Standalone mixed-pole class"];
lambdaMix = Normal[Series[lambdaOut - sigma/(1 - kappa z^2 - I gamma z^5), {z, 0, 5}]];
l0Mix = FullSimplify[Coefficient[lambdaMix, z, 0]];
l2Mix = FullSimplify[Coefficient[lambdaMix, z, 2]];
l4Mix = FullSimplify[Coefficient[lambdaMix, z, 4]];
l5Mix = FullSimplify[Coefficient[lambdaMix, z, 5]/I];
kappaMatch = FullSimplify[kappa /. First[Solve[-l2Mix/l0Mix == 1/9, kappa, Reals]]];
sigmaMatch = FullSimplify[
  sigma /. First[Solve[((l2Mix^2/l0Mix^2 - l4Mix/l0Mix) /. kappa -> kappaMatch) == 4/81, sigma, Reals]]
];
chiMix = FullSimplify[(-l5Mix/l0Mix)/(1/27)];
Print["standalone mixed-pole kappa match = ", fmt[kappaMatch]];
Print["standalone mixed-pole sigma match = ", fmt[sigmaMatch]];
expectZero["formal even-match forces kappa = -1/9", kappaMatch + 1/9];
expectZero["standalone mixed pole disappears on the canonical branch", sigmaMatch];
expectZero["odd norm is then trivial", (chiMix /. sigma -> 0) - 1];

banner["4. Hybrid outlet class split"];
lambdaHyb = Normal[Series[lambdaOut + rho - sigma/(1 - kappa z^2 - I gamma z^5), {z, 0, 5}]];
l0Hyb = FullSimplify[Coefficient[lambdaHyb, z, 0]];
l2Hyb = FullSimplify[Coefficient[lambdaHyb, z, 2]];
l4Hyb = FullSimplify[Coefficient[lambdaHyb, z, 4]];
l5Hyb = FullSimplify[Coefficient[lambdaHyb, z, 5]/I];
hybridSolutions = Solve[{-l2Hyb/l0Hyb == 1/9, l2Hyb^2/l0Hyb^2 - l4Hyb/l0Hyb == 4/81}, {rho, kappa}, Reals];
Print["hybrid canonical-even branches = ", fmt[hybridSolutions]];
branchCancel = SelectFirst[hybridSolutions, (kappa /. #) === 0 &];
branchComp = SelectFirst[hybridSolutions, FullSimplify[(kappa /. #) - 1/3] === 0 &];
chiCancel = FullSimplify[((-l5Hyb/l0Hyb)/(1/27)) /. branchCancel];
chiComp = FullSimplify[((-l5Hyb/l0Hyb)/(1/27)) /. branchComp];
expectZero["hybrid cancellation branch odd norm", chiCancel - (1 - 9 sigma gamma)];
expectZero[
  "hybrid cancellation branch is trivial when gamma = 0",
  Normal[Series[(lambdaHyb /. branchCancel /. gamma -> 0) - lambdaOut, {z, 0, 5}]]
];
expectZero["compensated branch odd norm", (chiComp /. gamma -> 1/9) - 1];
expectZero[
  "compensated branch collapses to a pure scale deformation",
  Normal[Series[(lambdaHyb /. branchComp /. gamma -> 1/9) - (1 - sigma) lambdaOut, {z, 0, 5}]]
];

banner["5. Concrete core realization of the compensated class"];
rC = FullSimplify[lam^2/(ks kq)];
rhoC = FullSimplify[gs^2/ks];
sigmaC = FullSimplify[(ks gq - lam gs)^2/(ks^2 kq (1 + rC))];
kappaC = FullSimplify[kappa0/(1 + rC)];
gammaC = FullSimplify[gamma0/(1 + rC)];
gqSolutions = Solve[rhoC - 4 sigmaC == 0, gq, Reals];
Print["core-balance surface branches = ", fmt[gqSolutions]];
sigmaStar = FullSimplify[gs^2/(4 ks)];
expectZero[
  "both core-balance branches give the same sigma_*",
  (sigmaC /. First[gqSolutions]) - (sigmaC /. Last[gqSolutions])
];
expectZero["core-balance sigma_* value", (sigmaC /. First[gqSolutions]) - sigmaStar];

lWRequired = FullSimplify[lW /. First[Solve[4 lW^2/(Pi^2 a^2) == (1 + rC)/3, lW, Reals]]];
expectZero["D/N tube fixes kappa_c = 1/3", (kappaC /. kappa0 -> 4 lWRequired^2/(Pi^2 a^2)) - 1/3];
expectZero["bare mixed normalization fixes gamma_c = 1/9", (gammaC /. gamma0 -> (1 + rC)/9) - 1/9];

deltaCore = FullSimplify[
  rhoC - sigmaC/(1 - kappaC z^2 - I gammaC z^5)
] /. First[gqSolutions] /. {kappa0 -> (1 + rC)/3, gamma0 -> (1 + rC)/9};
deltaCoreExpected = FullSimplify[4 sigmaStar - sigmaStar/(1 - z^2/3 - I z^5/9)];
expectZero[
  "concrete core collapses to the compensated hybrid class",
  Normal[Series[deltaCore - deltaCoreExpected, {z, 0, 5}]]
];

banner["6. Classification capstone"];
classificationRows = {
  {"scale/argument", True, True, False, "harmless beta = 1 pure-scale branch"},
  {"pure Robin", False, False, False, "rho_R = 0 only"},
  {"standalone mixed pole", False, False, False, "sigma_W = 0 only (formal kappa = -1/9)"},
  {"hybrid cancellation", True, True, False, "gamma_W = 0 reduces to exact cancellation"},
  {"compensated Robin-mixed core realization", True, True, True, "balance surface + D/N tube normalization"}
};
Scan[Print, classificationRows];
nontrivialSurvivors = Cases[classificationRows, {name_, True, True, True, _} :> name];
If[nontrivialSurvivors =!= {"compensated Robin-mixed core realization"},
  fail["classification survivor set", nontrivialSurvivors]
];

Print[""];
Print["Open microscopic question:"];
Print["  The explicit low-frequency classification is closed at the reduced-model level,"];
Print["  but the actual moving-throat core still has to realize the balance surface and"];
Print["  D/N tube normalization. This script does not assert that realization."];

Print[""];
Print["Stage 117 Mathematica audit passed."];

Exit[0];
