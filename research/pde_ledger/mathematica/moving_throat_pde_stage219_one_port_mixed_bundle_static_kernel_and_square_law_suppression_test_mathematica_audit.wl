ClearAll["Global`*"];
$HistoryLength = 0;

Print["Stage 219 Mathematica audit: one-port mixed-bundle static kernel"];

zeroQ[e_] := TrueQ[Simplify[Together[e]] == 0];

Kred = {
  {Kstar, -GU, -GW},
  {-GU, OmegaU^2, -R},
  {-GW, -R, OmegaW^2}
};

Delta = OmegaU^2 OmegaW^2 - R^2;
Q = GU^2 OmegaW^2 + 2 GU GW R + GW^2 OmegaU^2;
P = OmegaU^2 GW + R GU;
PU = GU OmegaW^2 + R GW;
D0 = Kstar - Q/Delta;

m1Residual = Det[Kred] - Delta*D0;
Print["M1 residual = ", Simplify[Together[m1Residual]]];
If[!zeroQ[m1Residual], Print["FAIL M1"]; Exit[1], Print["PASS M1"]];

Kint = {
  {OmegaU^2, -R},
  {-R, OmegaW^2}
};
cpl = {-GU, -GW};
schur = Kstar - cpl . Inverse[Kint] . cpl;
m2Residual = schur - D0;
Print["M2 residual = ", Simplify[Together[m2Residual]]];
If[!zeroQ[m2Residual], Print["FAIL M2"]; Exit[1], Print["PASS M2"]];

Kinv = Inverse[Kred];
m3Residuals = {
  Kinv[[1, 1]] - 1/D0,
  Kinv[[1, 2]] - PU/(Delta*D0),
  Kinv[[1, 3]] - P/(Delta*D0),
  Kinv[[2, 2]] - (Kstar*OmegaW^2 - GW^2)/(Delta*D0),
  Kinv[[2, 3]] - (Kstar*R + GU*GW)/(Delta*D0),
  Kinv[[3, 3]] - (Kstar*OmegaU^2 - GU^2)/(Delta*D0)
};
Print["M3 residuals = ", Simplify[Together /@ m3Residuals]];
If[!TrueQ[And @@ (zeroQ /@ m3Residuals)], Print["FAIL M3"]; Exit[1], Print["PASS M3"]];

J = S*{sq, sU, sW};
dV = -1/2*(J . Kinv . J);
Ns = (
  Delta*sq^2
  + 2 PU sq sU
  + 2 P sq sW
  + (Kstar OmegaW^2 - GW^2) sU^2
  + 2 (Kstar R + GU GW) sU sW
  + (Kstar OmegaU^2 - GU^2) sW^2
);
chiS = Ns/(Delta*D0);
m4Residual = dV - (-1/2 chiS S^2);
Print["M4 residual = ", Simplify[Together[m4Residual]]];
If[!zeroQ[m4Residual], Print["FAIL M4"]; Exit[1], Print["PASS M4"]];

Lambda = P/Delta;
N0 = Lambda^2;
P0 = N0/D0;
chiqW = Kinv[[1, 3]];
m5Residuals = {
  chiqW - Lambda/D0,
  chiqW^2 - P0/D0
};
Print["M5 residuals = ", Simplify[Together /@ m5Residuals]];
If[!TrueQ[And @@ (zeroQ /@ m5Residuals)], Print["FAIL M5"]; Exit[1], Print["PASS M5"]];

SQ = x^(-3);
SY = Exp[-2 kappa x]/x;
Jp = {betaQ*SQ, betaU*SY, betaW*SY};
dVp = -1/2*(Jp . Kinv . Jp);

C6 = (1/D0) betaQ^2;
C4 = (PU/(Delta*D0)) betaQ betaU + (P/(Delta*D0)) betaQ betaW;
C2 = (
  ((Kstar OmegaW^2 - GW^2)/(Delta*D0)) betaU^2
  + 2 ((Kstar R + GU GW)/(Delta*D0)) betaU betaW
  + ((Kstar OmegaU^2 - GU^2)/(Delta*D0)) betaW^2
);
m6Residual = dVp - (-1/2 (C6/x^6 + 2 C4 Exp[-2 kappa x]/x^4 + C2 Exp[-4 kappa x]/x^2));

scaledFamily = Simplify[Together[Expand[-2*dVp*x^6*Exp[4 kappa x]]]];
familyZ = Collect[
  Expand[scaledFamily /. {Exp[4 kappa x] -> z^2, Exp[2 kappa x] -> z}],
  {z, x},
  Simplify
];
zDegrees = Exponent[Coefficient[familyZ, z, #], x] & /@ {0, 1, 2};
x5Coefficients = Coefficient[familyZ, x, 5];
zCoefficientLists = CoefficientList[Coefficient[familyZ, z, #], x] & /@ {0, 1, 2};
m6StructuralOK = TrueQ[
  PolynomialQ[familyZ, {x, z}]
    && Exponent[familyZ, z] == 2
    && Exponent[familyZ, x] <= 4
    && zDegrees === {4, 2, 0}
    && zeroQ[x5Coefficients]
    && And @@ ((!zeroQ[Coefficient[familyZ, z, #]]) & /@ {0, 1, 2})
];

Print["M6 residual = ", Simplify[Together[m6Residual]]];
Print["M6 scaled family = ", familyZ];
Print["M6 x-degrees by {1, Exp[2 kappa x], Exp[4 kappa x]} = ", zDegrees];
Print["M6 coefficient-list lengths = ", Length /@ zCoefficientLists];
If[
  !TrueQ[zeroQ[m6Residual] && m6StructuralOK],
  Print["FAIL M6"]; Exit[1],
  Print["PASS M6"]
];

subs = {Kstar -> 11, OmegaU -> 3, OmegaW -> 4, R -> 2, GU -> 1, GW -> 2};
Knum = Kred /. subs;
deltaNum = Delta /. subs;
d0Num = Simplify[D0 /. subs];
detNum = Det[Knum];
minorNums = Det /@ Table[Knum[[1 ;; n, 1 ;; n]], {n, 1, 3}];
eigNums = Eigenvalues[Knum];
eigRealOK = And @@ Thread[Chop[Im[N[eigNums]]] == 0];
eigPositiveOK = And @@ Thread[Re[N[eigNums]] > 0];
m7Residuals = {
  deltaNum - 140,
  d0Num - 74/7,
  detNum - 1480,
  minorNums[[1]] - 11,
  minorNums[[2]] - 98,
  minorNums[[3]] - 1480
};
m7PositiveOK = TrueQ[
  And @@ (zeroQ /@ m7Residuals)
    && And @@ Thread[minorNums > 0]
    && PositiveDefiniteMatrixQ[Knum]
    && eigRealOK
    && eigPositiveOK
];

Print["M7 residuals = ", Simplify[Together /@ m7Residuals]];
Print["M7 leading principal minors = ", minorNums];
Print["M7 eigenvalues = ", N[eigNums]];
Print["M7 PositiveDefiniteMatrixQ = ", PositiveDefiniteMatrixQ[Knum]];
If[!m7PositiveOK, Print["FAIL M7"]; Exit[1], Print["PASS M7"]];

Print["Stage 219 Mathematica audit passed."];
Exit[0];
