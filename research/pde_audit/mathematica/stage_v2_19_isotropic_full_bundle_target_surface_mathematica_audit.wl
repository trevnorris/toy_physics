Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-19 isotropic full-bundle target-surface Mathematica audit"];

Clear[K, M, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4, theta, G, cs, a, c, mhat, Sport];
$Assumptions = Element[{K, M, B0, B2, B4, Z0, Z2, Z4, N0, N2, N4, theta}, Reals] &&
  G > 0 && cs > 0 && a > 0 && c > 0 && mhat > 0 && Sport > 0 &&
  B4 + Z4 != 0 && K - B0 - Z0 != 0;

D0 = K - B0 - Z0;
D2 = -(M + B2 + Z2);
D4 = -(B4 + Z4);
A = M + B2 + Z2;
Ciso = B4 + Z4;
Tgr = 54 G cs^5/(5 a^5 c^5);

resPole = D0 Ciso - 3 A^2;
resNorm = mhat^2 Sport N0 - Tgr D0;
resP2 = D0 N2 - 2 D2 N0;
resP4 = D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0;
resTail = theta (c/cs)^3 - 1;

D0surf = 3 A^2/Ciso;
N0surf = Tgr D0surf/(mhat^2 Sport);
N2surf = 2 D2 N0surf/D0surf;
N4surf = N0surf (D2^2 + 2 D0surf D4)/D0surf^2;
surface = {
  K -> B0 + Z0 + D0surf,
  N0 -> N0surf,
  N2 -> N2surf,
  N4 -> N4surf,
  theta -> (cs/c)^3
};

subbanner["Target residuals on surface"];
checkZero["one-pole residual", resPole /. surface];
checkZero["normalization residual", resNorm /. surface];
checkZero["P2 residual", resP2 /. surface];
checkZero["P4 residual", resP4 /. surface];
checkZero["tail transport residual", resTail /. surface];

subbanner["Local codimension"];
vars = {K, N0, N2, N4, theta};
jdet = Cancel[Together[Det[D[{resPole, resNorm, resP2, resP4, resTail}, {vars}]]];
checkTrue["target Jacobian nonzero generically", !TrueQ[jdet == 0], jdet];
Print["target_jacobian_det = ", fmt[jdet]];

subbanner["Stability implication"];
checkEqual["stable nondegenerate D0 surface", D0 /. surface, 3 A^2/Ciso];
Print["For D0>0 and A != 0, the one-pole surface requires C = B4+Z4 > 0."];
Print["All Stage V2-19 Mathematica checks passed."];
