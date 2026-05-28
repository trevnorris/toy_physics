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

banner["STAGE 174 — STATIC SELF-SIMILARITY DECOMPOSITION"];

Clear[k, b0, z0, k1, b01, z01, n0, n01];
$Assumptions = Element[{k, b0, z0, k1, b01, z01, n0, n01}, Reals] &&
  k != 0 && b0 != 0 && z0 != 0 && n0 != 0;

d0 = k - b0 - z0;
d01 = k1 - b01 - z01;

deltaK = k1/k;
deltaB = b01/b0;
deltaZ = z01/z0;
deltaN = n01/n0;
deltaD = d01/d0;

omegaK = k/d0;
omegaB = b0/d0;
omegaZ = z0/d0;

expectZero["weight identity", omegaK - omegaB - omegaZ - 1];
expectZero[
  "delta_D weighted decomposition",
  deltaD - (omegaK*deltaK - omegaB*deltaB - omegaZ*deltaZ)
];

xiLoad = FullSimplify[deltaN - deltaD, Assumptions -> $Assumptions];
xiWallRef = FullSimplify[
  (deltaN - deltaK) + omegaB*(deltaB - deltaK) + omegaZ*(deltaZ - deltaK),
  Assumptions -> $Assumptions
];
expectZero["Xi_load wall-referenced form", xiLoad - xiWallRef];
Print["Xi_load = ", fmt[xiLoad]];

banner["BdG bundle weighted logarithmic slope"];
Clear[c1, c2, w1, w2, dc1, dc2, dw1, dw2];
$Assumptions = Element[{c1, c2, w1, w2, dc1, dc2, dw1, dw2}, Reals] &&
  c1 != 0 && c2 != 0 && w1 != 0 && w2 != 0;

b1 = c1^2/w1^2;
b2Term = c2^2/w2^2;
b0Two = FullSimplify[b1 + b2Term, Assumptions -> $Assumptions];
b0Eps = (c1 + eps*dc1)^2/(w1 + eps*dw1)^2 + (c2 + eps*dc2)^2/(w2 + eps*dw2)^2;
b01Two = FullSimplify[D[b0Eps, eps] /. eps -> 0, Assumptions -> $Assumptions];

rhoB1 = FullSimplify[b1/b0Two, Assumptions -> $Assumptions];
rhoB2 = FullSimplify[b2Term/b0Two, Assumptions -> $Assumptions];
deltaBWeighted = FullSimplify[
  rhoB1*(2*dc1/c1 - 2*dw1/w1) + rhoB2*(2*dc2/c2 - 2*dw2/w2),
  Assumptions -> $Assumptions
];
expectZero["BdG weighted-average formula", b01Two/b0Two - deltaBWeighted];

banner["Conservative Maxwell/mixed weighted logarithmic slope"];
Clear[q1, q2, delta1, delta2, q1p, q2p, delta1p, delta2p];
$Assumptions = Element[{q1, q2, delta1, delta2, q1p, q2p, delta1p, delta2p}, Reals] &&
  q1 != 0 && q2 != 0 && delta1 != 0 && delta2 != 0;

z1 = q1/delta1;
z2Term = q2/delta2;
z0Two = FullSimplify[z1 + z2Term, Assumptions -> $Assumptions];
z0Eps = (q1 + eps*q1p)/(delta1 + eps*delta1p) + (q2 + eps*q2p)/(delta2 + eps*delta2p);
z01Two = FullSimplify[D[z0Eps, eps] /. eps -> 0, Assumptions -> $Assumptions];

rhoZ1 = FullSimplify[z1/z0Two, Assumptions -> $Assumptions];
rhoZ2 = FullSimplify[z2Term/z0Two, Assumptions -> $Assumptions];
deltaZWeighted = FullSimplify[
  rhoZ1*(q1p/q1 - delta1p/delta1) + rhoZ2*(q2p/q2 - delta2p/delta2),
  Assumptions -> $Assumptions
];
expectZero["Z weighted-average formula", z01Two/z0Two - deltaZWeighted];

banner["Outgoing-transfer weighted logarithmic slope"];
Clear[p1, p2, p1p, p2p];
$Assumptions = Element[{p1, p2, p1p, p2p, delta1, delta2, delta1p, delta2p}, Reals] &&
  p1 != 0 && p2 != 0 && delta1 != 0 && delta2 != 0;

n1s = p1^2/delta1^2;
n2s = p2^2/delta2^2;
n0Two = FullSimplify[n1s + n2s, Assumptions -> $Assumptions];
n0Eps = (p1 + eps*p1p)^2/(delta1 + eps*delta1p)^2 + (p2 + eps*p2p)^2/(delta2 + eps*delta2p)^2;
n01Two = FullSimplify[D[n0Eps, eps] /. eps -> 0, Assumptions -> $Assumptions];

rhoN1 = FullSimplify[n1s/n0Two, Assumptions -> $Assumptions];
rhoN2 = FullSimplify[n2s/n0Two, Assumptions -> $Assumptions];
deltaNWeighted = FullSimplify[
  rhoN1*(2*p1p/p1 - 2*delta1p/delta1) + rhoN2*(2*p2p/p2 - 2*delta2p/delta2),
  Assumptions -> $Assumptions
];
expectZero["N weighted-average formula", n01Two/n0Two - deltaNWeighted];

banner["Static self-similarity theorem"];
Clear[deltaStar, sigmaB, sigmaZ, sigmaN];
$Assumptions = Element[{k, b0, z0, n0, deltaStar, sigmaB, sigmaZ, sigmaN}, Reals] &&
  k != 0 && b0 != 0 && z0 != 0 && n0 != 0;

xiSelfSimilar = FullSimplify[
  xiLoad /. {k1 -> deltaStar*k, b01 -> deltaStar*b0, z01 -> deltaStar*z0, n01 -> deltaStar*n0},
  Assumptions -> $Assumptions
];
expectZero["Xi_load under static self-similarity", xiSelfSimilar];

xiSigma = FullSimplify[
  xiLoad /. {b01 -> (deltaK + sigmaB)*b0, z01 -> (deltaK + sigmaZ)*z0, n01 -> (deltaK + sigmaN)*n0},
  Assumptions -> $Assumptions
];
expectZero["Xi_load mismatch-field form", xiSigma - (sigmaN + omegaB*sigmaB + omegaZ*sigmaZ)];
Print["Xi_load mismatch-field prototype = ", fmt[FullSimplify[sigmaN + omegaB*sigmaB + omegaZ*sigmaZ, Assumptions -> $Assumptions]]];

Print[""];
Print["Carry-forward formulas:"];
Print["  delta_D = omegaK deltaK - omegaB deltaB - omegaZ deltaZ"];
Print["  omegaK - omegaB - omegaZ = 1"];
Print["  Xi_load = (deltaN-deltaK) + omegaB(deltaB-deltaK) + omegaZ(deltaZ-deltaK)"];
Print["  delta_B = weighted average of 2 dln(c_alpha/varpi_alpha)"];
Print["  delta_Z = weighted average of dln(Q_r/Delta_r)"];
Print["  delta_N = weighted average of 2 dln(P_r/Delta_r)"];
Print["  static self-similarity => Xi_load = 0"];

Print[""];
Print["Stage 174 Mathematica audit passed."];

Exit[0];
