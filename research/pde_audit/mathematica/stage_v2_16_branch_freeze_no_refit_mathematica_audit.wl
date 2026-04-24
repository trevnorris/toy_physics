Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-16 branch-freeze/no-refit Mathematica audit"];

Clear[K, M, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4, theta, G, cs, a, c, mhat, Sport];
$Assumptions = Element[{K, M, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4, theta}, Reals] &&
  G > 0 && cs > 0 && a > 0 && c > 0 && mhat > 0 && Sport > 0;

D0 = K - B0 - Z0;
D2 = -(M + B2 + Z2);
D4 = -(B4 + Z4);
A = M + B2 + Z2;
C = B4 + Z4;
Tgr = 54 G cs^5/(5 a^5 c^5);

residuals = {
  D0 C - 3 A^2,
  mhat^2 Sport N0 - Tgr D0,
  D0 N2 - 2 D2 N0,
  D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0,
  theta (c/cs)^3 - 1
};
vars = {K, N0, N2, N4, theta};
jac = D[residuals, {vars}];
jdet = Cancel[Together[Det[jac]]];

subbanner["Target residual slots"];
checkTrue["five residuals", Length[residuals] == 5];
checkTrue["five branch-output knobs", Length[vars] == 5];
checkTrue["Jacobian is generically nonzero", !TrueQ[jdet == 0], jdet];
Print["target_residual_jacobian_det = ", fmt[jdet]];

subbanner["No-refit discipline"];
checkZero["constant-prefactor polynomial P2 residual", residuals[[3]] /. N2 -> 2 D2 N0/D0];
checkZero[
  "constant-prefactor polynomial P4 residual",
  residuals[[4]] /. {N2 -> 2 D2 N0/D0, N4 -> N0 (D2^2 + 2 D0 D4)/D0^2}
];

Print["PASS: target algebra is tuneable unless branch data are frozen before residual evaluation."];
Print["All Stage V2-16 Mathematica checks passed."];
