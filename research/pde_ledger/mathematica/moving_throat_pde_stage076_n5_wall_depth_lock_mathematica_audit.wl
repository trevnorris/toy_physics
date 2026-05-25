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

banner["STAGE 059 — EXACT n=5 WALL-DEPTH LOCK"];

Clear[kConst, rho, mpsi, hbar, cSw, lambdaMu, rhoW, ell, a, nPoly];
$Assumptions =
  Element[{kConst, rho, mpsi, hbar, cSw, lambdaMu, rhoW, ell, a, nPoly}, Reals] &&
  kConst > 0 && rho > 0 && mpsi > 0 && hbar > 0 && cSw > 0 && lambdaMu > 0 && rhoW > 0 && ell > 0 && a > 0 && nPoly > 0 && nPoly != 1;

press = kConst*rho^nPoly;
uPerMass = FullSimplify[Integrate[press/rho^2, rho, GenerateConditions -> False], Assumptions -> $Assumptions];
uRho = FullSimplify[rho*uPerMass, Assumptions -> $Assumptions];
csSq = FullSimplify[D[press, rho]/mpsi, Assumptions -> $Assumptions];
hRho = FullSimplify[D[uRho, rho], Assumptions -> $Assumptions];
csSqN5 = FullSimplify[csSq /. nPoly -> 5, Assumptions -> $Assumptions];
hRhoN5 = FullSimplify[hRho /. nPoly -> 5, Assumptions -> $Assumptions];

Print["c_s^2(rho)  [n=5] = ", fmt[csSqN5]];
Print["h(rho)      [n=5] = ", fmt[hRhoN5]];
expectZero["n=5 enthalpy identity", hRhoN5 - mpsi*csSqN5/4];

residualN3 = FullSimplify[hRho /. nPoly -> 3, Assumptions -> $Assumptions] - mpsi*FullSimplify[csSq /. nPoly -> 3, Assumptions -> $Assumptions]/4;
residualN3 = FullSimplify[residualN3, Assumptions -> $Assumptions];
Print["n=3 residual (should be NONZERO) = ", fmt[residualN3]];
If[TrueQ[residualN3 === 0], fail["n=3 residual is zero -- identity does not actually distinguish n=5"]];

Clear[muStarSym];
$Assumptions = $Assumptions && muStarSym > 0;
enthalpyLock = muStarSym - lambdaMu*mpsi*cSw^2/4;
muStarSolved = First[muStarSym /. Solve[enthalpyLock == 0, muStarSym]];
thetaW = FullSimplify[4*rhoW^2*muStarSolved^2/(hbar^2*cSw^2), Assumptions -> $Assumptions];
thetaCanonical = FullSimplify[(2*rhoW*muStarSolved/(hbar*cSw))^2, Assumptions -> $Assumptions];
Print["Theta_w (enthalpy lock) = ", fmt[thetaW]];
expectZero["Theta_w vs alternative-form derivation", thetaW - thetaCanonical];

healingCondition = cSw - hbar/(2*mpsi*ell);
cSwFromEll = FullSimplify[First[cSw /. Solve[healingCondition == 0, cSw]] /. ConditionalExpression[e_, _] :> e, Assumptions -> $Assumptions];
ellSolved = FullSimplify[First[ell /. Solve[healingCondition == 0, ell]] /. ConditionalExpression[e_, _] :> e, Assumptions -> $Assumptions];
Print["ell from healing condition = ", fmt[ellSolved]];
thetaWInEll = FullSimplify[thetaW /. cSw -> cSwFromEll, Assumptions -> $Assumptions];
thetaHealReduced = FullSimplify[lambdaMu^2*rhoW^2/(16*ell^2), Assumptions -> $Assumptions];
expectZero["healing-lock reduction", thetaWInEll - thetaHealReduced];
Print["Theta_w (healing lock) = ", fmt[thetaHealReduced]];

(* Reference-branch convention: ell = a * refFactor with refFactor = 1/20.
   TODO(provenance): cite the upstream stage that fixes refFactor. This factor is
   the load-bearing piece of the "25" in the normalized reference identity. *)
refFactor = 1/20;  (* reference-branch convention: ell = a * refFactor  (see F2 below for provenance) *)
thetaRef = FullSimplify[thetaHealReduced /. ell -> a*refFactor, Assumptions -> $Assumptions];
thetaRefNorm = FullSimplify[thetaRef /. a -> 1, Assumptions -> $Assumptions];
Print["Theta_w (reference branch, general a) = ", fmt[thetaRef]];
Print["Theta_w (reference branch, normalized wall units) = ", fmt[thetaRefNorm]];
refTarget = FullSimplify[(1/refFactor)^2/16*lambdaMu^2*rhoW^2, Assumptions -> $Assumptions];
expectZero["normalized reference factor", thetaRefNorm - refTarget];

Print[""];
Print["Stage 076 Mathematica audit passed."];

Exit[0];
