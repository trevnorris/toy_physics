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

(* memoized canonical sphere moments *)
i4[i_, j_, k_, l_] /; {i, j, k, l} =!= Sort[{i, j, k, l}] := i4 @@ Sort[{i, j, k, l}];
i4[i_, j_, k_, l_] := i4[i, j, k, l] = Integrate[
  n[i] n[j] n[k] n[l] Sin[theta],
  {theta, 0, Pi}, {phi, 0, 2 Pi}
];

i6[i_, j_, k_, l_, m_, nn_] /; {i, j, k, l, m, nn} =!= Sort[{i, j, k, l, m, nn}] := i6 @@ Sort[{i, j, k, l, m, nn}];
i6[i_, j_, k_, l_, m_, nn_] := i6[i, j, k, l, m, nn] = Integrate[
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
expectZero["a_x witness - 1/5", (ax /. {x20 -> x0 + 1, x21 -> x0, x22 -> x0}) - 1/5];
expectZero["b_x witness", bx /. {x20 -> x0 + 1, x21 -> x0, x22 -> x0}];
expectZero["a_x second witness + 1/10", (ax /. {x20 -> x0, x21 -> x0 + 1, x22 -> x0}) + 1/10];
expectZero["b_x second witness - 1/2", (bx /. {x20 -> x0, x21 -> x0 + 1, x22 -> x0}) - 1/2];
Clear[p, q, rr];
$Assumptions = Element[{p, q, rr}, Reals];
xbarMix = xbar /. {x20 -> p, x21 -> q, x22 -> rr};
axMix = ax /. {x20 -> p, x21 -> q, x22 -> rr};
bxMix = bx /. {x20 -> p, x21 -> q, x22 -> rr};
x20Re = FullSimplify[xbarMix + 4*axMix, Assumptions -> $Assumptions];
x21Re = FullSimplify[xbarMix - axMix + bxMix, Assumptions -> $Assumptions];
x22Re = FullSimplify[xbarMix - axMix - bxMix, Assumptions -> $Assumptions];
expectZero["x20 reassembled - p", x20Re - p];
expectZero["x21 reassembled - q", x21Re - q];
expectZero["x22 reassembled - rr", x22Re - rr];

banner["SECTION III — ISOTROPIC RADIAL/AXIAL OVERLAP MOMENTS"];
Clear[omega, k, m, lambdaB1, lambdaB2, iEta1, iEta2, varpi1, varpi2];
$Assumptions =
  Element[{omega, k, m, lambdaB1, lambdaB2, iEta1, iEta2, varpi1, varpi2}, Reals] &&
  varpi1 > 0 && varpi2 > 0;

c1 = lambdaB1*iEta1;
c2 = lambdaB2*iEta2;

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

(* Per-pair conservative 2x2 matrix for the (U,W) modes. *)
mPair = {{omegaU^2 - omega^2, -rPair}, {-rPair, omegaW^2 - omega^2}};
coupling = {gU, gW};

(* Z response: contract coupling^T . mPair^{-1} . coupling. *)
zFromMatrix = FullSimplify[
  First[First[{coupling}.Inverse[mPair].Transpose[{coupling}]]],
  Assumptions -> $Assumptions
];

(* N response: square of the W-component projection from the matrix inverse. *)
nFromMatrix = FullSimplify[
  ((Inverse[mPair].coupling)[[2]])^2,
  Assumptions -> $Assumptions
];

(* Reference closed forms from the paper card (Eqs B-moments, ZN-moments). *)
qRef = gU^2*omegaW^2 + 2*gU*gW*rPair + gW^2*omegaU^2;
hRef = gU^2 + gW^2;
pRef = omegaU^2*gW + rPair*gU;
deltaRef = omegaU^2*omegaW^2 - rPair^2;
sRef = omegaU^2 + omegaW^2;
zRefRational = (qRef - hRef*omega^2)/(deltaRef - sRef*omega^2 + omega^4);
nRefRational = (pRef - gW*omega^2)^2/(deltaRef - sRef*omega^2 + omega^4)^2;

(* Anchor: matrix-inverse derivation matches the paper's closed-form rational. *)
expectZero["Z_full from matrix inverse matches paper rational", zFromMatrix - zRefRational];
expectZero["N_full from matrix inverse matches paper rational", nFromMatrix - nRefRational];

zSeries = Expand[Normal[Series[zFromMatrix, {omega, 0, 4}]]];
z0 = FullSimplify[Coefficient[zSeries, omega, 0], Assumptions -> $Assumptions];
z2 = FullSimplify[Coefficient[zSeries, omega, 2], Assumptions -> $Assumptions];
z4 = FullSimplify[Coefficient[zSeries, omega, 4], Assumptions -> $Assumptions];

nSeries = Expand[Normal[Series[nFromMatrix, {omega, 0, 4}]]];
n0 = FullSimplify[Coefficient[nSeries, omega, 0], Assumptions -> $Assumptions];
n2 = FullSimplify[Coefficient[nSeries, omega, 2], Assumptions -> $Assumptions];
n4 = FullSimplify[Coefficient[nSeries, omega, 4], Assumptions -> $Assumptions];

dSeries = Expand[Normal[Series[k - m*omega^2 - bResp - zFromMatrix, {omega, 0, 4}]]];
d0 = FullSimplify[Coefficient[dSeries, omega, 0], Assumptions -> $Assumptions];
d2 = FullSimplify[Coefficient[dSeries, omega, 2], Assumptions -> $Assumptions];
d4 = FullSimplify[Coefficient[dSeries, omega, 4], Assumptions -> $Assumptions];

expectZero["Z0 formula", z0 - qRef/deltaRef];
expectZero["Z2 formula", z2 - (qRef*sRef - hRef*deltaRef)/deltaRef^2];
expectZero["Z4 formula", z4 - (qRef*(sRef^2 - deltaRef) - sRef*hRef*deltaRef)/deltaRef^3];
expectZero["N0 formula", n0 - pRef^2/deltaRef^2];
expectZero["N2 formula", n2 - 2*pRef*(pRef*sRef - deltaRef*gW)/deltaRef^3];
expectZero["N4 formula", n4 - (deltaRef^2*gW^2 - 2*deltaRef*pRef^2 - 4*deltaRef*pRef*sRef*gW + 3*pRef^2*sRef^2)/deltaRef^4];
expectZero["D0 formula", d0 - (k - b0 - z0)];
expectZero["D2 formula", d2 + (m + b2 + z2)];
expectZero["D4 formula", d4 + (b4 + z4)];

banner["SECTION III.5 — LANE COLLAPSE UNDER O(3) INVARIANCE"];
Clear[gU20, gU21, gU22, gW20, gW21, gW22, rr20, rr21, rr22, cc20, cc21, cc22, varpi, omU, omW, gUiso, gWiso, rrIso, cIso, delta];
$Assumptions =
  Element[{omega, k, m, gU20, gU21, gU22, gW20, gW21, gW22, rr20, rr21, rr22, cc20, cc21, cc22, varpi, omU, omW, gUiso, gWiso, rrIso, cIso, delta}, Reals] &&
  varpi > 0 && omU > 0 && omW > 0;

perLaneD[gUA_, gWA_, rrA_, cA_] := Module[{bA, qA, hA, sA, deltaA, zA},
  bA = cA^2/(varpi^2 - omega^2);
  qA = gUA^2*omW^2 + 2*gUA*gWA*rrA + gWA^2*omU^2;
  hA = gUA^2 + gWA^2;
  sA = omU^2 + omW^2;
  deltaA = omU^2*omW^2 - rrA^2;
  zA = (qA - hA*omega^2)/(deltaA - sA*omega^2 + omega^4);
  FullSimplify[k - m*omega^2 - bA - zA, Assumptions -> $Assumptions]
];
perLaneN[gUA_, gWA_, rrA_] := Module[{pA, sA, deltaA},
  pA = omU^2*gWA + rrA*gUA;
  sA = omU^2 + omW^2;
  deltaA = omU^2*omW^2 - rrA^2;
  FullSimplify[(pA - gWA*omega^2)^2/(deltaA - sA*omega^2 + omega^4)^2, Assumptions -> $Assumptions]
];
d20Lane = perLaneD[gU20, gW20, rr20, cc20];
d21Lane = perLaneD[gU21, gW21, rr21, cc21];
d22Lane = perLaneD[gU22, gW22, rr22, cc22];
n20Lane = perLaneN[gU20, gW20, rr20];
n21Lane = perLaneN[gU21, gW21, rr21];
n22Lane = perLaneN[gU22, gW22, rr22];
isoSubs = {gU20 -> gUiso, gU21 -> gUiso, gU22 -> gUiso, gW20 -> gWiso, gW21 -> gWiso, gW22 -> gWiso, rr20 -> rrIso, rr21 -> rrIso, rr22 -> rrIso, cc20 -> cIso, cc21 -> cIso, cc22 -> cIso};
expectZero["D_20 - D_21 (isotropic)", FullSimplify[(d20Lane - d21Lane) /. isoSubs, Assumptions -> $Assumptions]];
expectZero["D_21 - D_22 (isotropic)", FullSimplify[(d21Lane - d22Lane) /. isoSubs, Assumptions -> $Assumptions]];
expectZero["D_20 - D_22 (isotropic)", FullSimplify[(d20Lane - d22Lane) /. isoSubs, Assumptions -> $Assumptions]];
expectZero["N_20 - N_21 (isotropic)", FullSimplify[(n20Lane - n21Lane) /. isoSubs, Assumptions -> $Assumptions]];
expectZero["N_21 - N_22 (isotropic)", FullSimplify[(n21Lane - n22Lane) /. isoSubs, Assumptions -> $Assumptions]];
expectZero["N_20 - N_22 (isotropic)", FullSimplify[(n20Lane - n22Lane) /. isoSubs, Assumptions -> $Assumptions]];
breakSubs = {gU20 -> gUiso + delta, gU21 -> gUiso, gU22 -> gUiso, gW20 -> gWiso, gW21 -> gWiso, gW22 -> gWiso, rr20 -> rrIso, rr21 -> rrIso, rr22 -> rrIso, cc20 -> cIso, cc21 -> cIso, cc22 -> cIso};
diffLin = Normal[Series[(d20Lane - d21Lane) /. breakSubs, {delta, 0, 1}]];
coeffDelta = FullSimplify[Coefficient[diffLin, delta, 1], Assumptions -> $Assumptions];
Print["linear coefficient of delta in (D_20 - D_21) = ", fmt[coeffDelta]];
If[TrueQ[coeffDelta === 0], fail["Lane-breaking witness produced no defect"], pass["Lane-breaking witness: collapse check is non-tautological"]];

banner["SECTION IV — AXISYMMETRIC SPLITTING MATRIX"];
(* Reset symbol context before sphere-integral table — F3/F4 additions leaked symbols. *)
ClearAll[gU, gW, rPair, omegaU, omegaW, mPair, coupling, zFromMatrix, nFromMatrix,
  qRef, hRef, pRef, deltaRef, sRef, zRefRational, nRefRational,
  zSeries, z0, z2, z4, nSeries, n0, n2, n4, dSeries, d0, d2, d4,
  gUiso, gWiso, rrIso, omU, omW];
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
Clear[d0sym, d1sym, n0sym, n1sym, eps, lambdaSym];
$Assumptions = Element[{eps, d0sym, d1sym, n0sym, n1sym, lambdaSym}, Reals] && d0sym != 0;

(* Direct quotient-rule derivation: P_A(eps) = (n0 + eps*lambda*n1)/(d0 + eps*lambda*d1).
   Compute P_A and its first eps-derivative symbolically and assemble the
   first-order expansion from the derivative, independent of SymPy's series-helper. *)
pAFull[lam_] := (n0sym + eps*lam*n1sym)/(d0sym + eps*lam*d1sym);
pAAt0[lam_] := pAFull[lam] /. eps -> 0;
pAFirst[lam_] := FullSimplify[D[pAFull[lam], eps] /. eps -> 0, Assumptions -> $Assumptions];

p0Closed = n0sym/d0sym;
p1Closed = FullSimplify[(n1sym*d0sym - n0sym*d1sym)/d0sym^2, Assumptions -> $Assumptions];

(* Anchor each lane's first-order expansion derived from the quotient rule
   against the paper-card formula. *)
expectZero["P20 first-order from quotient rule",
  pAFirst[1] - 1*p1Closed
];
expectZero["P21 first-order from quotient rule",
  pAFirst[1/2] - (1/2)*p1Closed
];
expectZero["P22 first-order from quotient rule",
  pAFirst[-1] - (-1)*p1Closed
];
expectZero["P0 from quotient rule", pAAt0[1] - p0Closed];

(* Defect map applied to the first-order expansion (1, 2, 2) weighting from Stage 022. *)
p20Lin = pAAt0[1] + eps*pAFirst[1];
p21Lin = pAAt0[1/2] + eps*pAFirst[1/2];
p22Lin = pAAt0[-1] + eps*pAFirst[-1];
pbarChk = FullSimplify[(p20Lin + 2*p21Lin + 2*p22Lin)/5, Assumptions -> $Assumptions];
aPChk = FullSimplify[(2*p20Lin - p21Lin - p22Lin)/10, Assumptions -> $Assumptions];
bPChk = FullSimplify[(p21Lin - p22Lin)/2, Assumptions -> $Assumptions];

expectZero["Pbar - P0", pbarChk - p0Closed];
expectZero["a_P - eps*P1/4", aPChk - eps*p1Closed/4];
expectZero["b_P - 3 eps*P1/4", bPChk - 3*eps*p1Closed/4];
expectZero["b_P - 3 a_P", bPChk - 3*aPChk];

Print[""];
Print["FINAL STAGE-007 LEDGER:"];
Print["  Verified the normalized real STF harmonic orthonormality, the exact angular"];
Print["  source-map identity, the O(3) isotropy collapse, representative isotropic"];
Print["  overlap formulas for C_alpha/B_n/Z_n/N_n/D(omega), the axisymmetric Y20"];
Print["  triple-overlap matrix, the grouped splitting law b = 3 a, and the"];
Print["  corresponding first-order transport law for P_A = N_A / D_A."];

Exit[0];
