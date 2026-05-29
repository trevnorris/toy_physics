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

banner["STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q"];

Clear[z, sNorm, beta, sigma0, sigma2, sigma4, sigma5];
$Assumptions =
  Element[{z, sNorm, beta, sigma0, sigma2, sigma4, sigma5}, Reals] &&
  sNorm != 0 && beta != 0 && 3*sNorm - sigma0 != 0;

lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9;

yScale = Expand[Normal[Series[(-3*sNorm)/(sNorm*lambdaOut), {z, 0, 5}]]];
(* anchor to literal canonical fingerprint (appendix eq:app-part04-Yout-dtn) *)
yCanLiteral = 1 + z^2/9 + (4*z^4)/81 + (I/27)*z^5;
expectZero["pure scale invariance", yScale - yCanLiteral];

yArg = Expand[Normal[Series[(-3*sNorm)/(sNorm*(lambdaOut /. z -> beta*z)), {z, 0, 5}]]];
m2Arg = FullSimplify[Coefficient[yArg, z, 2], Assumptions -> $Assumptions];
m4Arg = FullSimplify[Coefficient[yArg, z, 4], Assumptions -> $Assumptions];
chiArg = FullSimplify[(Coefficient[yArg, z, 5]/I)/(1/27), Assumptions -> $Assumptions];

Print["Y_scale+arg(z) = ", fmt[yArg]];
Print["m2_arg = ", fmt[m2Arg]];
Print["m4_arg = ", fmt[m4Arg]];
Print["chi_arg = ", fmt[chiArg]];
expectZero["m2_arg - beta^2/9", m2Arg - beta^2/9];
expectZero["m4_arg - 4 beta^4/81", m4Arg - 4*beta^4/81];
expectZero["chi_arg - beta^5", chiArg - beta^5];

betaRoots = Sort[beta /. Solve[{beta^2 == 1, beta^4 == 1}, beta, Reals]];
Print["solutions preserving canonical even fingerprint = ", fmt[betaRoots]];
If[betaRoots =!= {-1, 1}, fail["expected beta roots { -1, 1 }", betaRoots]];
expectZero["chi_arg(beta=1) - 1", (chiArg /. beta -> 1) - 1];

lambdaAdd = Expand[sNorm*lambdaOut + sigma0 + sigma2*z^2 + sigma4*z^4 + I*sigma5*z^5];
l0 = FullSimplify[lambdaAdd /. z -> 0, Assumptions -> $Assumptions];
l2 = FullSimplify[Coefficient[lambdaAdd, z, 2], Assumptions -> $Assumptions];
l4 = FullSimplify[Coefficient[lambdaAdd, z, 4], Assumptions -> $Assumptions];
l5 = FullSimplify[Coefficient[lambdaAdd, z, 5]/I, Assumptions -> $Assumptions];

m2 = FullSimplify[-l2/l0, Assumptions -> $Assumptions];
m4 = FullSimplify[l2^2/l0^2 - l4/l0, Assumptions -> $Assumptions];
sol = Solve[{m2 == 1/9, m4 == 4/81}, {sigma2, sigma4}, Reals];
If[Length[sol] =!= 1, fail["unique additive even-match solution", sol]];
sol = First[sol];

Print["Sigma2(beta=1) = ", fmt[sigma2 /. sol]];
Print["Sigma4(beta=1) = ", fmt[sigma4 /. sol]];
expectZero["Sigma2(beta=1) + sigma0/9", (sigma2 /. sol) + sigma0/9];
expectZero["Sigma4(beta=1) + sigma0/27", (sigma4 /. sol) + sigma0/27];

chiAdd = FullSimplify[((-l5/l0)/(1/27)) /. sol, Assumptions -> $Assumptions];
Print["chi_add = ", fmt[Factor[chiAdd]]];
expectZero["chi_add - 3(sNorm + 9 sigma5)/(3 sNorm - sigma0)", chiAdd - 3*(sNorm + 9*sigma5)/(3*sNorm - sigma0)];

sigma5Pres = FullSimplify[sigma5 /. First[Solve[chiAdd == 1, sigma5, Reals]], Assumptions -> $Assumptions];
Print["Sigma5 preservation locus = ", fmt[sigma5Pres]];
expectZero["Sigma5 locus (Class C) = -sigma0/27", sigma5Pres - (-sigma0/27)];

(* Class D: general scale + argument + additive (β-parameterized preservation submanifold). *)
(* Notes box: Σ_5 = S(1 - β^5)/9 - Σ_0/27 (general locus); Class C is the β=1 reduction. *)
lambdaGen = Expand[
  sNorm*(lambdaOut /. z -> beta*z)
  + sigma0 + sigma2*z^2 + sigma4*z^4 + I*sigma5*z^5
];
l0g = FullSimplify[lambdaGen /. z -> 0, Assumptions -> $Assumptions];
l2g = FullSimplify[Coefficient[lambdaGen, z, 2], Assumptions -> $Assumptions];
l4g = FullSimplify[Coefficient[lambdaGen, z, 4], Assumptions -> $Assumptions];
l5g = FullSimplify[Coefficient[lambdaGen, z, 5]/I, Assumptions -> $Assumptions];
m2g = FullSimplify[-l2g/l0g, Assumptions -> $Assumptions];
m4g = FullSimplify[l2g^2/l0g^2 - l4g/l0g, Assumptions -> $Assumptions];
solG = Solve[{m2g == 1/9, m4g == 4/81}, {sigma2, sigma4}, Reals];
If[Length[solG] =!= 1, fail["unique β-parameterized even-match solution", solG]];
solG = First[solG];
Print["Sigma2(beta) = ", fmt[Simplify[sigma2 /. solG]]];
Print["Sigma4(beta) = ", fmt[Simplify[sigma4 /. solG]]];
chiGen = FullSimplify[((-l5g/l0g)/(1/27)) /. solG, Assumptions -> $Assumptions];
Print["chi_gen(beta) = ", fmt[Factor[chiGen]]];
(* Independent route: build chiGen from raw L-coefficients without Series. *)
L0raw = -3*sNorm + sigma0;
L2raw = sNorm*beta^2/3 + sigma2;
L4raw = sNorm*beta^4/9 + sigma4;
L5raw = sNorm*beta^5/9 + sigma5;
solAlt = First[Solve[{-L2raw/L0raw == 1/9, L2raw^2/L0raw^2 - L4raw/L0raw == 4/81},
                     {sigma2, sigma4}, Reals]];
chiGenAlt = FullSimplify[(27*(-L5raw/L0raw)) /. solAlt, Assumptions -> $Assumptions];
expectZero["chiGen independent-route agreement", chiGenAlt - chiGen];
sigma5PresGen = FullSimplify[sigma5 /. First[Solve[chiGen == 1, sigma5, Reals]], Assumptions -> $Assumptions];
Print["Sigma5 general preservation locus = ", fmt[sigma5PresGen]];
expectZero[
  "general preservation submanifold = S(1 - beta^5)/9 - sigma0/27",
  sigma5PresGen - (sNorm*(1 - beta^5)/9 - sigma0/27)
];
Print["general preservation locus check = ", fmt[FullSimplify[(chiGen /. sigma5 -> sigma5PresGen) - 1, Assumptions -> $Assumptions]]];
expectZero[
  "general locus reduces to beta=1 (Class C)",
  (sigma5PresGen - (-sigma0/27)) /. beta -> 1
];

Print[""];
Print["Stage 108 Mathematica audit passed."];

Exit[0];
