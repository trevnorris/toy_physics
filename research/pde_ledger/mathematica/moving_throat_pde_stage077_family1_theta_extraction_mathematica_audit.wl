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

expectApprox[name_String, value_, target_, tol_] := Module[{diff},
  diff = Abs[N[value, 40] - N[target, 40]];
  Print[name, " diff = ", fmt[diff]];
  If[TrueQ[diff <= tol], pass[name], fail[name, diff]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 060 — SHELL-WEIGHTED THETA EXTRACTION"];

Clear[xi, alphaR];
$Assumptions = Element[{xi, alphaR}, Reals] && alphaR > 0;

sXi = (1 + Tanh[xi])/2;
chi = FullSimplify[D[sXi, xi], Assumptions -> $Assumptions];
ifMom = FullSimplify[Integrate[chi^2, {xi, -Infinity, Infinity}], Assumptions -> $Assumptions];
xiStar = FullSimplify[ArcTanh[2/Sqrt[alphaR] - 1], Assumptions -> $Assumptions];

Print["chi_phi(xi) = ", fmt[chi]];
Print["I_f = ", fmt[ifMom]];
expectZero["I_f - 1/3", ifMom - 1/3];
Print["xi_* = ", fmt[xiStar]];
Print["xi_*(alpha_r=10) = ", fmt[N[xiStar /. alphaR -> 10, 30]]];
sAtStar = (1 + Tanh[xiStar])/2;
rhoQuarticAtStar = FullSimplify[1 - alphaR*sAtStar^2, Assumptions -> $Assumptions];
expectZero["1 - alphaR*S[xi_*]^2", rhoQuarticAtStar];

banner["NUMERICAL FAMILY-1 EXTRACTION"];

alphaNum = 10;
xiCut = N[xiStar /. alphaR -> alphaNum, 50];

sNum[x_?NumericQ] := (1 + Tanh[x])/2;
chiNum[x_?NumericQ] := 1/(2*Cosh[x]^2);
rhoNum[x_?NumericQ] := (1 - alphaNum*sNum[x]^2)^(1/4);
rhoSqNum[x_?NumericQ] := Sqrt[1 - alphaNum*sNum[x]^2];

denNum = Quiet[
  NIntegrate[chiNum[x]^2, {x, -Infinity, Infinity}, WorkingPrecision -> 60, AccuracyGoal -> 28, PrecisionGoal -> 28],
  NIntegrate::precw
];
num1 = Quiet[
  NIntegrate[chiNum[x]^2*rhoNum[x], {x, -Infinity, xiCut}, WorkingPrecision -> 60, AccuracyGoal -> 28, PrecisionGoal -> 28],
  NIntegrate::precw
];
num2 = Quiet[
  NIntegrate[chiNum[x]^2*rhoSqNum[x], {x, -Infinity, xiCut}, WorkingPrecision -> 60, AccuracyGoal -> 28, PrecisionGoal -> 28],
  NIntegrate::precw
];

r1 = N[num1/denNum, 50];
r2 = N[num2/denNum, 50];
thetaChi = N[25*r2, 50];
thetaJ = N[25*r1^2, 50];

Print["numeric xi_* = ", fmt[xiCut]];
Print["<rho>_chi = ", fmt[r1]];
Print["<rho^2>_chi = ", fmt[r2]];
Print["denominator = ", fmt[denNum]];
Print["Theta_w^(chi) / lambda_mu^2 = ", fmt[thetaChi]];
Print["Theta_w^(J) / lambda_mu^2 = ", fmt[thetaJ]];

expectApprox["<rho>_chi numeric check", r1, ToExpression["0.19261900555649309777068139356018510792903510747507`50"], 10^-28];
expectApprox["<rho^2>_chi numeric check", r2, ToExpression["0.16274529400326462037087418498629868328210821103971`50"], 10^-28];
expectApprox["Theta_w^(chi) numeric check", thetaChi, ToExpression["4.0686323500816155092718546246574670820527052759928`50"], 10^-26];
expectApprox["Theta_w^(J) numeric check", thetaJ, ToExpression["0.92755203253930797183993260663904217023332624032789`50"], 10^-27];
expectTrue["Theta_w^(chi) >= Theta_w^(J) > 0", thetaChi >= thetaJ && thetaJ > 0];

Print[""];
Print["Stage 077 Mathematica audit passed."];

Exit[0];
