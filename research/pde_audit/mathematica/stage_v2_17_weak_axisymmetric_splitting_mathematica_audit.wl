Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-17 weak-axisymmetric splitting Mathematica audit"];

Clear[x0, x1, eps, D0, D2, D4, D01, D21, D41, lam];
$Assumptions = D0 != 0 && Element[{x0, x1, eps, D0, D2, D4, D01, D21, D41, lam}, Reals];

subbanner["Grouped P2 axisymmetric line"];
x20 = x0 + eps x1;
x21 = x0 + eps x1/2;
x22 = x0 - eps x1;
checkEqual["trace unchanged", weightedBar[x20, x21, x22], x0];
checkZero["axisymmetric b equals 3a", groupB[x20, x21, x22] - 3 groupA[x20, x21, x22]];

subbanner["Hidden-even response condition"];
D0lane = D0 + eps lam D01;
D2lane = D2 + eps lam D21;
D4lane = D4 + eps lam D41;
u2lane = -D2lane/D0lane;
u4lane = (D2lane^2 - D0lane D4lane)/D0lane^2;
canonical = {D2 -> -D0/9, D4 -> -D0/27};
u2slope = Cancel[Together[D[u2lane, eps] /. eps -> 0 /. canonical]];
u4slope = Cancel[Together[D[u4lane, eps] /. eps -> 0 /. canonical]];
hiddenResidual = Cancel[Together[u4slope - (8/9) u2slope]];
hiddenExpected = lam (D01 + 18 D21 - 27 D41)/(27 D0);
checkZero["hidden-even residual formula", hiddenResidual - hiddenExpected];
checkZero["hidden-even condition", hiddenResidual /. D41 -> D01/27 + 2 D21/3];

Print["u2_slope = ", fmt[u2slope]];
Print["u4_slope = ", fmt[u4slope]];
Print["All Stage V2-17 Mathematica checks passed."];
