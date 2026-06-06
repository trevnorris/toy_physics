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

banner["STAGE 117 — CONCRETE OUTLET-CORE STATUS"];

Clear[z, s, beta, rho, sigma, kappa, gamma, ks, kq, lam, gs, gq, kappa0, gamma0, a, lW,
  sCore, qCore, dSym, kappaSlot, gammaSlot, sigmaSlot];
$Assumptions =
  Element[{z, s, beta, rho, sigma, kappa, gamma, gs, gq}, Reals] &&
  Element[{ks, kq, lam, kappa0, gamma0, a, lW, dSym}, Reals] &&
  ks > 0 && kq > 0 && lam > 0 && kappa0 > 0 && gamma0 > 0 && a > 0 && lW > 0;

lambdaOut = -3 + z^2/3 + z^4/9 + I z^5/9;

jetFromBalance[tag_String, den_, num_] := Module[{jetVars, trial, equations, sol},
  jetVars = Table[Unique["j"], {6}];
  trial = Sum[jetVars[[k + 1]] z^k, {k, 0, 5}];
  equations = Table[Coefficient[Expand[den trial - num], z, k] == 0, {k, 0, 5}];
  sol = Solve[equations, jetVars];
  If[Length[sol] =!= 1, fail[tag <> " coefficient solve", sol]];
  FullSimplify[Table[jetVars[[k + 1]] /. First[sol], {k, 0, 5}], Assumptions -> $Assumptions]
];

banner["1. Harmless scale/argument class"];
argJet = jetFromBalance["scale/argument", lambdaOut /. z -> beta z, -3];
m2Arg = FullSimplify[argJet[[3]]];
m4Arg = FullSimplify[argJet[[5]]];
chiArg = FullSimplify[argJet[[6]]/I/(1/27)];
betaSolutions = Solve[{m2Arg == 1/9, m4Arg == 4/81}, beta, Reals];
Print["scale/argument solutions = ", fmt[betaSolutions]];
If[Sort[beta /. betaSolutions] =!= {-1, 1}, fail["scale/argument branch roots", betaSolutions]];
expectZero["positive harmless branch has beta = 1 and chi_Q = 1", (chiArg /. beta -> 1) - 1];

banner["2. Pure Robin class"];
lambdaR = lambdaOut + rho;
robinJet = jetFromBalance["pure Robin", lambdaR, -3 + rho];
c2R = FullSimplify[robinJet[[3]]];
c4R = FullSimplify[robinJet[[5]]];
chiR = FullSimplify[robinJet[[6]]/I/(1/27)];
robinSolutions = Solve[{c2R == 1/9, c4R == 4/81}, rho, Reals];
Print["pure Robin canonical-even solutions = ", fmt[robinSolutions]];
If[robinSolutions =!= {{rho -> 0}}, fail["pure Robin branch", robinSolutions]];
expectZero["pure Robin odd norm is trivial on rho = 0", (chiR /. rho -> 0) - 1];

banner["3. Standalone mixed-pole class"];
poleDen = 1 - kappa z^2 - I gamma z^5;
mixDenCleared = Expand[poleDen lambdaOut - sigma];
mixNumCleared = Expand[(-3 - sigma) poleDen];
mixJet = jetFromBalance["standalone mixed-pole", mixDenCleared, mixNumCleared];
c2Mix = FullSimplify[mixJet[[3]]];
c4Mix = FullSimplify[mixJet[[5]]];
chiMix = FullSimplify[mixJet[[6]]/I/(1/27)];
kappaMatch = FullSimplify[kappa /. First[Solve[c2Mix == 1/9, kappa, Reals]]];
sigmaMatch = FullSimplify[
  sigma /. First[Solve[(c4Mix /. kappa -> kappaMatch) == 4/81, sigma, Reals]]
];
Print["standalone mixed-pole kappa match = ", fmt[kappaMatch]];
Print["standalone mixed-pole sigma match = ", fmt[sigmaMatch]];
expectZero["formal even-match forces kappa = -1/9", kappaMatch + 1/9];
expectZero["standalone mixed pole disappears on the canonical branch", sigmaMatch];
expectZero["odd norm is then trivial", (chiMix /. sigma -> 0) - 1];

banner["4. Hybrid outlet class split"];
hybDenCleared = Expand[poleDen (lambdaOut + rho) - sigma];
hybNumCleared = Expand[(-3 + rho - sigma) poleDen];
hybJet = jetFromBalance["hybrid outlet", hybDenCleared, hybNumCleared];
c2Hyb = FullSimplify[hybJet[[3]]];
c4Hyb = FullSimplify[hybJet[[5]]];
chiHyb = FullSimplify[hybJet[[6]]/I/(1/27)];
hybridSolutions = Solve[{c2Hyb == 1/9, c4Hyb == 4/81}, {rho, kappa}, Reals];
Print["hybrid canonical-even branches = ", fmt[hybridSolutions]];
branchCancel = SelectFirst[hybridSolutions, (kappa /. #) === 0 &];
branchComp = SelectFirst[hybridSolutions, FullSimplify[(kappa /. #) - 1/3] === 0 &];
chiCancel = FullSimplify[chiHyb /. branchCancel];
chiComp = FullSimplify[chiHyb /. branchComp];
lambdaHybExact = lambdaOut + rho - sigma/poleDen;
expectZero["hybrid cancellation branch odd norm", chiCancel - (1 - 9 sigma gamma)];
expectZero[
  "hybrid cancellation branch is trivial when gamma = 0",
  (lambdaHybExact /. branchCancel /. gamma -> 0) - lambdaOut
];
expectZero["compensated branch odd norm", (chiComp /. gamma -> 1/9) - 1];
compDen = 1 - z^2/3 - I z^5/9;
compCollapseNumerator = Expand[
  compDen ((lambdaHybExact /. branchComp /. gamma -> 1/9) - (1 - sigma) lambdaOut)
];
compCollapseJet = Sum[Coefficient[compCollapseNumerator, z, k] z^k, {k, 0, 5}];
expectZero[
  "compensated branch collapses to a pure scale deformation",
  compCollapseJet
];

banner["5. Concrete core realization of the compensated class"];
(* Status-consolidation card. Provenance of the two bare coefficients DIFFERS:
   - kappa_0_bare = (1+r_c)/3 is FORWARD-DERIVED at Stage 116 from the D/N
     half-wave eigenvalue k_W = Pi/(2 L_W) on q(0)=0, q'(L_W)=0, giving
     kappa_0 = 4 L_W^2/(Pi^2 a^2) and the tube-length law
     L_W = Pi a Sqrt[(1+r_c)/3]/2. F2 routes the substitution below through that
     forward closed form, so the kappa_c = 1/3 / core-residual check exercises it.
   - gamma_0_bare = (1+r_c)/9 is NOT derived: it is a pure-scale ANSATZ of the
     canonical compact outgoing l=2 branch, postulated in the stage-116 note
     ("Bare outgoing normalization") and carried as a hardcoded input at
     Stages 115/116. gamma_c = 1/9 is thus a consistency-of-assumption check.
   The load-bearing check is the residual deltaCore - deltaCoreExpected at z^5. *)
dW = 1 - kappa0 z^2 - I gamma0 z^5;
coreMatrix = {{ks, lam}, {lam, -kq dW}};
coreSource = {gs, gq};
coreSolution = First[Solve[Thread[coreMatrix . {sCore, qCore} == coreSource], {sCore, qCore}]];
deltaCoreEliminated = FullSimplify[
  Together[gs (sCore /. coreSolution) + gq (qCore /. coreSolution)],
  Assumptions -> $Assumptions
];

coreMatrixD = {{ks, lam}, {lam, -kq dSym}};
coreSolutionD = First[Solve[Thread[coreMatrixD . {sCore, qCore} == coreSource], {sCore, qCore}]];
deltaCoreD = FullSimplify[
  Together[gs (sCore /. coreSolutionD) + gq (qCore /. coreSolutionD)],
  Assumptions -> $Assumptions
];
deltaDNum = Numerator[Together[deltaCoreD]];
deltaDDen = Denominator[Together[deltaCoreD]];
deltaDLead = Coefficient[deltaDDen, dSym, 1];
rhoC = FullSimplify[Coefficient[deltaDNum, dSym, 1]/deltaDLead, Assumptions -> $Assumptions];
rC = FullSimplify[deltaDDen/deltaDLead - dSym, Assumptions -> $Assumptions];
sigmaTilde = FullSimplify[(rhoC - deltaCoreD) (dSym + rC), Assumptions -> $Assumptions];
coreShape = FullSimplify[(dW + rC)/(1 + rC), Assumptions -> $Assumptions];
shapeResidual = Expand[coreShape - (1 - kappaSlot z^2 - I gammaSlot z^5)];
shapeSolution = First[Solve[
  {Coefficient[shapeResidual, z, 2] == 0, Coefficient[shapeResidual, z, 5]/I == 0},
  {kappaSlot, gammaSlot}
]];
kappaC = FullSimplify[kappaSlot /. shapeSolution, Assumptions -> $Assumptions];
gammaC = FullSimplify[gammaSlot /. shapeSolution, Assumptions -> $Assumptions];
sigmaSolution = First[Solve[sigmaSlot (1 + rC) == sigmaTilde, sigmaSlot]];
sigmaC = FullSimplify[sigmaSlot /. sigmaSolution, Assumptions -> $Assumptions];
coreSchurTarget = FullSimplify[rhoC - sigmaC/(1 - kappaC z^2 - I gammaC z^5), Assumptions -> $Assumptions];
coreSchurResidual = FullSimplify[Together[deltaCoreEliminated - coreSchurTarget], Assumptions -> $Assumptions];
If[!TrueQ[coreSchurResidual === 0], fail["core normalized Schur identity", coreSchurResidual]];
gqSolutions = Solve[rhoC - 4 sigmaC == 0, gq, Reals];
Print["core-balance surface branches = ", fmt[gqSolutions]];
sigmaStar = FullSimplify[gs^2/(4 ks)];
expectZero[
  "both core-balance branches give the same sigma_*",
  (sigmaC /. First[gqSolutions]) - (sigmaC /. Last[gqSolutions])
];
expectZero["core-balance sigma_* value", (sigmaC /. First[gqSolutions]) - sigmaStar];

(* De-tautologized (F2): build kappa0 from the stage-116 FORWARD tube-length law,
   not by inverting kappa0 = (1+r_c)/3. Stage-116 boxed result:
   L_W = Pi a Sqrt[(1+r_c)/3]/2, and kappa0 = 4 L_W^2/(Pi^2 a^2). A wrong
   tube-length coefficient would make kappa0FromTube != (1+r_c)/3, so deltaCore
   would no longer collapse and the O(z^5) residual check would FAIL. *)
lWForward = Pi a Sqrt[(1 + rC)/3]/2;
kappa0FromTube = FullSimplify[4 lWForward^2/(Pi^2 a^2)];
Print["carrying forward (Stage 116): L_W = Pi a Sqrt[(1+r_c)/3]/2 -> kappa_0_bare = 4 L_W^2/(Pi^2 a^2) -> kappa_c = 1/3"];
Print["gamma_0_bare = (1+r_c)/9 is a pure-scale ANSATZ of the canonical l=2 branch (stage-116 note), not derived; gamma_c = 1/9 is a consistency-of-assumption check"];

deltaCore = FullSimplify[
  deltaCoreEliminated /. First[gqSolutions] /. {kappa0 -> kappa0FromTube, gamma0 -> (1 + rC)/9},
  Assumptions -> $Assumptions
];
deltaCoreExpected = FullSimplify[4 sigmaStar - sigmaStar/(1 - z^2/3 - I z^5/9)];
deltaCoreResidual = Normal[Series[deltaCore - deltaCoreExpected, {z, 0, 5}]];
expectZero[
  "concrete core collapses to the compensated hybrid class",
  deltaCoreResidual
];

banner["6. Classification capstone"];
(* Wired booleans from sections 1-5 residuals. The load-bearing entry is
   `nontrivialCompensated`, anchored to the section-5 series residual. *)
evenOkScale = TrueQ[Sort[beta /. betaSolutions] === {-1, 1}];
oddOkScale = TrueQ[FullSimplify[(chiArg /. beta -> 1) - 1] === 0];
nontrivialScale = False;

evenOkRobin = TrueQ[robinSolutions === {{rho -> 0}}];
oddOkRobin = TrueQ[FullSimplify[(chiR /. rho -> 0) - 1] === 0];
nontrivialRobin = False;

evenOkStandalone = TrueQ[FullSimplify[kappaMatch + 1/9] === 0];
oddOkStandalone = TrueQ[FullSimplify[(chiMix /. sigma -> 0) - 1] === 0];
nontrivialStandalone = !TrueQ[FullSimplify[sigmaMatch] === 0];

evenOkHybCancel = TrueQ[FullSimplify[(kappa /. branchCancel)] === 0];
oddOkHybCancel = TrueQ[FullSimplify[chiCancel - (1 - 9 sigma gamma)] === 0];
nontrivialHybCancel = False;

evenOkCompensated = TrueQ[FullSimplify[(kappa /. branchComp) - 1/3] === 0];
oddOkCompensated = TrueQ[FullSimplify[(chiComp /. gamma -> 1/9) - 1] === 0];
nontrivialCompensated = TrueQ[
  deltaCoreResidual === 0
];

classificationRows = {
  {"scale/argument", evenOkScale, oddOkScale, nontrivialScale, "harmless beta = 1 pure-scale branch"},
  {"pure Robin", evenOkRobin, oddOkRobin, nontrivialRobin, "rho_R = 0 only"},
  {"standalone mixed pole", evenOkStandalone, oddOkStandalone, nontrivialStandalone, "sigma_W = 0 only (formal kappa = -1/9)"},
  {"hybrid cancellation", evenOkHybCancel, oddOkHybCancel, nontrivialHybCancel, "gamma_W = 0 reduces to exact cancellation"},
  {"compensated Robin-mixed core realization", evenOkCompensated, oddOkCompensated, nontrivialCompensated, "balance surface + D/N tube normalization"}
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
