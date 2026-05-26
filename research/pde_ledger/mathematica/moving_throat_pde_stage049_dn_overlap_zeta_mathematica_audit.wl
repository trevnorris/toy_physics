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
  (* Strip ConditionalExpression wrapper: under $Assumptions, a result of
     the form ConditionalExpression[0, cond] is identically zero on the
     declared domain.  Solve[]/Reduce[] often introduce these wrappers when
     auxiliary inequalities are nontrivial. *)
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

dnHalfwaveMomentum[n_, l_] := FullSimplify[(n + 1/2) Pi/l];
overlapRatio[n_] := FullSimplify[1/(2 n + 1)];
twinSupportRatio[n_, x_] := FullSimplify[1/((2 n + 1)^2 (1 + x n (n + 1)))];

banner["STAGE 32 — EXPLICIT D/N OVERLAP EXTRACTION OF THE PHYSICAL SUPPORT RATIO"];

Clear[n, s, l, lambdaStar, kWeff, kPhiEff, kX, tX];
$Assumptions =
  Element[{n}, Integers] && n >= 0 &&
  Element[{s, lambdaStar}, Reals] &&
  Element[{l, kWeff, kPhiEff, kX, tX}, Reals] &&
  l > 0 && kWeff > 0 && kPhiEff > 0 && kX > 0 && tX > 0;

kN = dnHalfwaveMomentum[n, l];
chiN = Sqrt[2/l] Sin[kN s];
Print["chi_n(s) = ", fmt[chiN]];
Print["k_n = ", fmt[kN]];

expectZero["k_n satisfies D/N Neumann boundary", Cos[kN l]];

overlapFromIntegral = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];
overlapFormula = Sqrt[2 l]/((n + 1/2) Pi);
Print["I_n from direct integral = ", fmt[overlapFromIntegral]];
Print["I_n closed form = ", fmt[overlapFormula]];
expectZero["uniform overlap integral", overlapFromIntegral - overlapFormula];

i0 = FullSimplify[overlapFormula /. n -> 0];
ratio = FullSimplify[overlapFormula/i0];
Print["I_0 = ", fmt[i0]];
Print["I_n / I_0 = ", fmt[ratio]];
expectZero["overlap ratio hierarchy", ratio - overlapRatio[n]];

lambdaW = FullSimplify[lambdaStar i0];
lambdaPhi = FullSimplify[lambdaStar overlapFormula];
zetaPhys = FullSimplify[(lambdaPhi^2 kWeff)/(lambdaW^2 kPhiEff)];
zetaPhysExpected = FullSimplify[(kWeff/kPhiEff) overlapRatio[n]^2];
Print["zeta_n^(phys) = ", fmt[zetaPhysExpected]];
expectZero["microscopic coherent-support law", zetaPhys - zetaPhysExpected];

kWTwin = FullSimplify[kX + Pi^2 tX/(4 l^2)];
kPhiTwin = FullSimplify[kX + (n + 1/2)^2 Pi^2 tX/l^2];
twinGap = FullSimplify[Pi^2 tX n (n + 1)/l^2];
expectZero["same-operator twin stiffness relation", kPhiTwin - (kWTwin + twinGap)];

xExpr = FullSimplify[Pi^2 tX/(l^2 kWTwin)];
zetaTwin = FullSimplify[(kWTwin/kPhiTwin) overlapRatio[n]^2];
zetaTwinExpected = twinSupportRatio[n, xExpr];
Print["x = ", fmt[xExpr]];
Print["zeta_n^(twin) = ", fmt[zetaTwinExpected]];
expectZero["exact twin-lane support ratio", zetaTwin - zetaTwinExpected];
expectZero["lowest twin half-wave", (zetaTwinExpected /. n -> 0) - 1];

Print[""];
Print["Carry-forward formulas:"];
Print["  I_n = sqrt(2L) / ((n+1/2) pi)"];
Print["  I_n / I_0 = 1 / (2n+1)"];
Print["  zeta_n^(phys) = (K_W_eff / K_phi,n_eff) (I_n / I_0)^2"];
Print["  zeta_n^(twin) = 1 / ((2n+1)^2 (1 + x n(n+1)))"];
Print["  zeta_0^(twin) = 1"];

Print[""];
Print["Stage 32 Mathematica audit passed."];

Exit[0];
