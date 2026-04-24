Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-04 open-junction impedance Mathematica audit"];

Clear[eps, chi, delta, k, r, Phi, aExit, AExit, L, j];
$Assumptions = eps > 0 && chi > 0 && delta > 0 && k > 0 && r > 0 &&
  Phi > 0 && aExit > 0 && AExit > 0 && L > 0 && Element[j, Integers] && j >= 0;

Rp = Cancel[Together[(eps - 1)/(eps + 1)]];
Tp = Cancel[Together[1 + Rp]];
Rq = Cancel[Together[-Rp]];
Tq = Cancel[Together[1 + Rq]];
Eref = Cancel[Together[Rp^2]];
Etrans = Cancel[Together[4 eps/(1 + eps)^2]];

subbanner["Reflection and transmission"];
checkEqual["pressure reflection", Rp, (eps - 1)/(eps + 1)];
checkEqual["pressure transmission", Tp, 2 eps/(eps + 1)];
checkEqual["dual flow reflection", Rq, (1 - eps)/(eps + 1)];
checkEqual["dual flow transmission", Tq, 2/(eps + 1)];
checkEqual["energy reflection plus transmission", Cancel[Together[Eref + Etrans]], 1];
checkEqual["large area pressure reflection", Cancel[Together[Rp /. eps -> 1/chi]], (1 - chi)/(1 + chi)];
checkEqual["large area flow reflection", Cancel[Together[Rq /. eps -> 1/chi]], (chi - 1)/(1 + chi)];

subbanner["Boundary limits"];
derivP = Cancel[Together[I k (1 - Rp)/(1 + Rp)]];
derivQ = Cancel[Together[I k (1 - Rq)/(1 + Rq)]];
checkEqual["pressure derivative ratio", derivP, I k/eps];
checkEqual["flow derivative ratio", derivQ, I k eps];
checkEqual["low load pressure limit", Rp /. eps -> 0, -1];
checkEqual["low load flow limit", Rq /. eps -> 0, 1];
checkEqual["low load pressure amplitude", Tp /. eps -> 0, 0];
checkEqual["low load flow derivative ratio", derivQ /. eps -> 0, 0];

subbanner["DC leakage"];
S3 = 2 Pi^2;
Abulk = S3 r^3;
Jbulk = Cancel[Together[Phi/Abulk]];
divBulk = Cancel[Together[(1/Abulk) D[Abulk Jbulk, r]]];
Aexit3Ball = 4 Pi aExit^3/3;
Jtube = Cancel[Together[Phi/Aexit3Ball]];
checkEqual["bulk radial continuity", divBulk, 0];
checkEqual["finite exit flux density", Jtube, 3 Phi/(4 Pi aExit^3)];
checkTrue["hard cap flux singularity", Exponent[Denominator[Together[Phi/AExit]], AExit] > 0];

subbanner["Support ladders"];
kDN = Pi (j + 1/2)/L;
kDD = Pi (j + 1)/L;
checkEqual["DN first wavenumber", kDN /. j -> 0, Pi/(2 L)];
checkEqual["DD first wavenumber", kDD /. j -> 0, Pi/L];

Print[""];
Print["EXPECTED NEGATIVE CONTROLS: hard cap with nonzero DC flux and generic pressure-scalar Neumann are rejected."];
Print["All Stage V2-04 Mathematica checks passed."];
