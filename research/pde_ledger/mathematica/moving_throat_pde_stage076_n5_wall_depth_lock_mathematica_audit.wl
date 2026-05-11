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

Clear[kConst, rho, mpsi, hbar, cSw, lambdaMu, rhoW, ell, a];
$Assumptions =
  Element[{kConst, rho, mpsi, hbar, cSw, lambdaMu, rhoW, ell, a}, Reals] &&
  kConst > 0 && rho > 0 && mpsi > 0 && hbar > 0 && cSw > 0 && lambdaMu > 0 && rhoW > 0 && ell > 0 && a > 0;

press = kConst*rho^5;
uRho = kConst*rho^5/4;
csSq = FullSimplify[D[press, rho]/mpsi, Assumptions -> $Assumptions];
hRho = FullSimplify[D[uRho, rho], Assumptions -> $Assumptions];

Print["c_s^2(rho) = ", fmt[csSq]];
Print["h(rho) = ", fmt[hRho]];
expectZero["n=5 enthalpy identity", hRho - mpsi*csSq/4];

muStar = FullSimplify[lambdaMu*mpsi*cSw^2/4, Assumptions -> $Assumptions];
thetaW = FullSimplify[4*rhoW^2*muStar^2/(hbar^2*cSw^2), Assumptions -> $Assumptions];
thetaExpected = FullSimplify[lambdaMu^2*mpsi^2*rhoW^2*cSw^2/(4*hbar^2), Assumptions -> $Assumptions];

Print["Theta_w (enthalpy lock) = ", fmt[thetaW]];
expectZero["Theta_w - expected", thetaW - thetaExpected];

thetaHeal = FullSimplify[
  thetaW /. cSw -> hbar/(2*mpsi*ell),
  Assumptions -> $Assumptions
];
thetaHealExpected = FullSimplify[lambdaMu^2*rhoW^2/(16*ell^2), Assumptions -> $Assumptions];
expectZero["healing-lock reduction", thetaHeal - thetaHealExpected];
Print["Theta_w (healing lock) = ", fmt[thetaHealExpected]];

thetaRef = FullSimplify[thetaHealExpected /. ell -> a/20, Assumptions -> $Assumptions];
thetaRefNorm = FullSimplify[thetaRef /. a -> 1, Assumptions -> $Assumptions];
Print["Theta_w (reference branch, general a) = ", fmt[thetaRef]];
Print["Theta_w (reference branch, normalized wall units) = ", fmt[thetaRefNorm]];
expectZero["normalized reference factor", thetaRefNorm - 25*lambdaMu^2*rhoW^2];

Print[""];
Print["Stage 076 Mathematica audit passed."];

Exit[0];
