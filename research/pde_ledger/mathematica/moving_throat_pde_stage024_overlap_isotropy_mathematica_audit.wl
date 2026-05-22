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
  If[ListQ[expr],
    res = Map[FullSimplify[Together[Expand[#]], Assumptions -> $Assumptions] &, expr, {ArrayDepth[expr]}];
    Print[name, " = ", fmt[res]];
    If[TrueQ[res === ConstantArray[0, Dimensions[expr]]], pass[name], fail[name, res]],
    res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
    Print[name, " = ", fmt[res]];
    If[TrueQ[res === 0], pass[name], fail[name, res]]
  ]
];

n[1] = Sin[theta] Cos[phi];
n[2] = Sin[theta] Sin[phi];
n[3] = Cos[theta];

i4[i_, j_, k_, l_] := Integrate[
  n[i] n[j] n[k] n[l] Sin[theta],
  {theta, 0, Pi}, {phi, 0, 2 Pi}
];

i6[i_, j_, k_, l_, m_, nn_] := Integrate[
  n[i] n[j] n[k] n[l] n[m] n[nn] Sin[theta],
  {theta, 0, Pi}, {phi, 0, 2 Pi}
];

quadOverlap[aMat_, bMat_] := FullSimplify[
  Sum[aMat[[i, j]]*bMat[[k, l]]*i4[i, j, k, l], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}, {l, 1, 3}],
  Assumptions -> $Assumptions
];

tripleOverlap[aMat_, qMat_, bMat_] := FullSimplify[
  Sum[
    aMat[[i, j]]*qMat[[k, l]]*bMat[[m, nn]]*i6[i, j, k, l, m, nn],
    {i, 1, 3}, {j, 1, 3}, {k, 1, 3}, {l, 1, 3}, {m, 1, 3}, {nn, 1, 3}
  ],
  Assumptions -> $Assumptions
];

banner["STAGE 007 — OVERLAP ISOTROPY"];

$Assumptions = True;
s5 = Sqrt[5];
nrm = Sqrt[15/(8*Pi)];

e20 = {{-1/Sqrt[6], 0, 0}, {0, -1/Sqrt[6], 0}, {0, 0, 2/Sqrt[6]}};
e21c = {{0, 0, 1/Sqrt[2]}, {0, 0, 0}, {1/Sqrt[2], 0, 0}};
e21s = {{0, 0, 0}, {0, 0, 1/Sqrt[2]}, {0, 1/Sqrt[2], 0}};
e22c = {{1/Sqrt[2], 0, 0}, {0, -1/Sqrt[2], 0}, {0, 0, 0}};
e22s = {{0, 1/Sqrt[2], 0}, {1/Sqrt[2], 0, 0}, {0, 0, 0}};
basis = nrm*# & /@ {e20, e21c, e21s, e22c, e22s};

banner["SECTION I — NORMALIZED HARMONICS AND SOURCE MAP"];
gram = Table[quadOverlap[basis[[i]], basis[[j]]], {i, 1, 5}, {j, 1, 5}];
expectZero["Gram - I5", gram - IdentityMatrix[5]];

Clear[s20, s21c, s21s, s22c, s22s];
$Assumptions = Element[{s20, s21c, s21s, s22c, s22s}, Reals];
svec = {s20, s21c, s21s, s22c, s22s};
expectZero["projected coefficients - source coefficients", gram.svec - svec];

banner["SECTION II — ISOTROPIC GROUPED-BUNDLE COLLAPSE"];
Clear[x0, x20, x21, x22];
$Assumptions = Element[{x0, x20, x21, x22}, Reals];
xbar = (x20 + 2*x21 + 2*x22)/5;
ax = (2*x20 - x21 - x22)/10;
bx = (x21 - x22)/2;
expectZero["xbar - x0", (xbar /. {x20 -> x0, x21 -> x0, x22 -> x0}) - x0];
expectZero["a_x on equal lanes", ax /. {x20 -> x0, x21 -> x0, x22 -> x0}];
expectZero["b_x on equal lanes", bx /. {x20 -> x0, x21 -> x0, x22 -> x0}];
expectZero["a_x witness - 1/5", (ax /. {x20 -> x0 + 1, x21 -> x0, x22 -> x0}) - 1/5];
expectZero["b_x witness", bx /. {x20 -> x0 + 1, x21 -> x0, x22 -> x0}];
expectZero["a_x second witness + 1/10", (ax /. {x20 -> x0, x21 -> x0 + 1, x22 -> x0}) + 1/10];
expectZero["b_x second witness - 1/2", (bx /. {x20 -> x0, x21 -> x0 + 1, x22 -> x0}) - 1/2];

banner["SECTION III — ISOTROPIC RADIAL/AXIAL OVERLAP MOMENTS"];
Clear[omega, k, m, lambdaB1, lambdaB2, iEta1, iEta2, varpi1, varpi2];
$Assumptions =
  Element[{omega, k, m, lambdaB1, lambdaB2, iEta1, iEta2, varpi1, varpi2}, Reals] &&
  varpi1 > 0 && varpi2 > 0;

c1 = lambdaB1*iEta1;
c2 = lambdaB2*iEta2;
expectZero["C_alpha - lambda_B,alpha I_etaalpha", {c1, c2} - {lambdaB1*iEta1, lambdaB2*iEta2}];

bResp = FullSimplify[c1^2/(varpi1^2 - omega^2) + c2^2/(varpi2^2 - omega^2), Assumptions -> $Assumptions];
bSeries = Expand[Normal[Series[bResp, {omega, 0, 4}]]];
b0 = FullSimplify[Coefficient[bSeries, omega, 0], Assumptions -> $Assumptions];
b2 = FullSimplify[Coefficient[bSeries, omega, 2], Assumptions -> $Assumptions];
b4 = FullSimplify[Coefficient[bSeries, omega, 4], Assumptions -> $Assumptions];
expectZero["B0 sum formula", b0 - (c1^2/varpi1^2 + c2^2/varpi2^2)];
expectZero["B2 sum formula", b2 - (c1^2/varpi1^4 + c2^2/varpi2^4)];
expectZero["B4 sum formula", b4 - (c1^2/varpi1^6 + c2^2/varpi2^6)];

Clear[lambdaU, lambdaW, lambdaR, iEtaU, iEtaW, iUW, omegaU, omegaW, rMix];
$Assumptions =
  Element[{omega, k, m, lambdaU, lambdaW, lambdaR, iEtaU, iEtaW, iUW, omegaU, omegaW, rMix}, Reals];

gU = lambdaU*iEtaU;
gW = lambdaW*iEtaW;
rPair = lambdaR*iUW;

zResp = Expand[Normal[Series[
  ((gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2) - (gU^2 + gW^2)*omega^2) /
  ((omegaU^2*omegaW^2 - rPair^2) - (omegaU^2 + omegaW^2)*omega^2 + omega^4),
  {omega, 0, 4}
]]];
z0 = FullSimplify[Coefficient[zResp, omega, 0], Assumptions -> $Assumptions];
z2 = FullSimplify[Coefficient[zResp, omega, 2], Assumptions -> $Assumptions];
z4 = FullSimplify[Coefficient[zResp, omega, 4], Assumptions -> $Assumptions];
nResp = Expand[Normal[Series[
  ((omegaU^2*gW + rPair*gU) - gW*omega^2)^2 /
  ((omegaU^2*omegaW^2 - rPair^2) - (omegaU^2 + omegaW^2)*omega^2 + omega^4)^2,
  {omega, 0, 4}
]]];
n0 = FullSimplify[Coefficient[nResp, omega, 0], Assumptions -> $Assumptions];
n2 = FullSimplify[Coefficient[nResp, omega, 2], Assumptions -> $Assumptions];
n4 = FullSimplify[Coefficient[nResp, omega, 4], Assumptions -> $Assumptions];
dResp = Expand[Normal[Series[k - m*omega^2 - bResp - zResp, {omega, 0, 4}]]];
d0 = FullSimplify[Coefficient[dResp, omega, 0], Assumptions -> $Assumptions];
d2 = FullSimplify[Coefficient[dResp, omega, 2], Assumptions -> $Assumptions];
d4 = FullSimplify[Coefficient[dResp, omega, 4], Assumptions -> $Assumptions];

expectZero["Z0 formula", z0 - (gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2)/(omegaU^2*omegaW^2 - rPair^2)];
expectZero["Z2 formula", z2 - ((gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2)*(omegaU^2 + omegaW^2) - (gU^2 + gW^2)*(omegaU^2*omegaW^2 - rPair^2))/(omegaU^2*omegaW^2 - rPair^2)^2];
expectZero["Z4 formula", z4 - ((gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2)*((omegaU^2 + omegaW^2)^2 - (omegaU^2*omegaW^2 - rPair^2)) - (omegaU^2 + omegaW^2)*(gU^2 + gW^2)*(omegaU^2*omegaW^2 - rPair^2))/(omegaU^2*omegaW^2 - rPair^2)^3];
expectZero["N0 formula", n0 - (omegaU^2*gW + rPair*gU)^2/(omegaU^2*omegaW^2 - rPair^2)^2];
expectZero["N2 formula", n2 - 2*(omegaU^2*gW + rPair*gU)*((omegaU^2*gW + rPair*gU)*(omegaU^2 + omegaW^2) - (omegaU^2*omegaW^2 - rPair^2)*gW)/(omegaU^2*omegaW^2 - rPair^2)^3];
expectZero["N4 formula", n4 - ((omegaU^2*omegaW^2 - rPair^2)^2*gW^2 - 2*(omegaU^2*omegaW^2 - rPair^2)*(omegaU^2*gW + rPair*gU)^2 - 4*(omegaU^2*omegaW^2 - rPair^2)*(omegaU^2*gW + rPair*gU)*(omegaU^2 + omegaW^2)*gW + 3*(omegaU^2*gW + rPair*gU)^2*(omegaU^2 + omegaW^2)^2)/(omegaU^2*omegaW^2 - rPair^2)^4];
expectZero["D0 formula", d0 - (k - b0 - z0)];
expectZero["D2 formula", d2 + (m + b2 + z2)];
expectZero["D4 formula", d4 + (b4 + z4)];

banner["SECTION IV — AXISYMMETRIC SPLITTING MATRIX"];
$Assumptions = True;
qMat = basis[[1]];
m20 = Table[tripleOverlap[basis[[i]], qMat, basis[[j]]], {i, 1, 5}, {j, 1, 5}];
kappaStar = Sqrt[5]/(7*Sqrt[Pi]);
mtarget = DiagonalMatrix[{kappaStar, kappaStar/2, kappaStar/2, -kappaStar, -kappaStar}];
expectZero["M - M_target", m20 - mtarget];

Clear[eps, x1];
$Assumptions = Element[{x0, eps, x1}, Reals];
xLane[lam_] := x0 + eps*lam*x1;
x20ax = xLane[1];
x21ax = xLane[1/2];
x22ax = xLane[-1];
xbarAx = FullSimplify[(x20ax + 2*x21ax + 2*x22ax)/5, Assumptions -> $Assumptions];
axAx = FullSimplify[(2*x20ax - x21ax - x22ax)/10, Assumptions -> $Assumptions];
bxAx = FullSimplify[(x21ax - x22ax)/2, Assumptions -> $Assumptions];
expectZero["xbar axisymmetric - x0", xbarAx - x0];
expectZero["a_x - eps*x1/4", axAx - eps*x1/4];
expectZero["b_x - 3 eps*x1/4", bxAx - 3*eps*x1/4];
expectZero["b_x - 3 a_x", bxAx - 3*axAx];

banner["SECTION V — FIRST-ORDER TRANSPORT LAW"];
Clear[d0, d1, n0, n1];
$Assumptions = Element[{eps, d0, d1, n0, n1}, Reals] && d0 != 0;
laneRatio[lam_] := Expand[Normal[Series[(n0 + eps*lam*n1)/(d0 + eps*lam*d1), {eps, 0, 1}]]];
p20 = laneRatio[1];
p21 = laneRatio[1/2];
p22 = laneRatio[-1];
p0 = n0/d0;
p1 = FullSimplify[(n1*d0 - n0*d1)/d0^2, Assumptions -> $Assumptions];
expectZero["P20 - (P0 + eps P1)", p20 - (p0 + eps*p1)];
expectZero["P21 - (P0 + eps P1/2)", p21 - (p0 + eps*p1/2)];
expectZero["P22 - (P0 - eps P1)", p22 - (p0 - eps*p1)];
pbar = FullSimplify[(p20 + 2*p21 + 2*p22)/5, Assumptions -> $Assumptions];
aP = FullSimplify[(2*p20 - p21 - p22)/10, Assumptions -> $Assumptions];
bP = FullSimplify[(p21 - p22)/2, Assumptions -> $Assumptions];
expectZero["Pbar - P0", pbar - p0];
expectZero["a_P - eps*P1/4", aP - eps*p1/4];
expectZero["b_P - 3 eps*P1/4", bP - 3*eps*p1/4];
expectZero["b_P - 3 a_P", bP - 3*aP];

Print[""];
Print["FINAL STAGE-007 LEDGER:"];
Print["  Verified the normalized real STF harmonic orthonormality, the exact angular"];
Print["  source-map identity, the O(3) isotropy collapse, representative isotropic"];
Print["  overlap formulas for C_alpha/B_n/Z_n/N_n/D(omega), the axisymmetric Y20"];
Print["  triple-overlap matrix, the grouped splitting law b = 3 a, and the"];
Print["  corresponding first-order transport law for P_A = N_A / D_A."];

Exit[0];
