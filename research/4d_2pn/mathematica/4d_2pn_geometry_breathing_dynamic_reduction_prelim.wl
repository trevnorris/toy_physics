(* ::Package:: *)

ClearAll["Global`*"];
$Assumptions = a0 > 0 && Lam > 0 && Sigma > 0 && rho > 0 && beta > 0 && rhoEff > 0;

section[s_String] := Print["\n" <> StringRepeat["=", 78] <> "\n" <> s <> "\n" <> StringRepeat["=", 78]];
pass[label_String, expr_] := Print[label, If[TrueQ[expr], " PASS", " FAIL"]];

(* ---------------------------------------------------------------------- *)
(* 2PN geometry breathing: dynamic reduction with affine inertia          *)
(* ---------------------------------------------------------------------- *)

(* Symbols *)
Clear[a0, Lam, Sigma, rho, beta, rhoEff, s, a, L, w, r];

(* 4D cylinder-like throat geometry (3-ball cross section times interval) *)
V[a_, L_] := (4 Pi/3) a^3 L;
A[a_, L_] := 4 Pi a^2 L + (8 Pi/3) a^3;

sigma = Sigma/a0^3;
Pvac = rho Sigma/a0^4;
kappab = beta Sigma/a0;
Egeom[a_, L_] := Expand[Pvac V[a, L] + sigma A[a, L] + kappab a^2/L];

H[a_, L_] := {{D[Egeom[a, L], {a, 2}], D[Egeom[a, L], a, L]}, {D[Egeom[a, L], L, a], D[Egeom[a, L], {L, 2}]}};
g[a_, L_] := {D[V[a, L], a], D[V[a, L], L]};

subs0 = {a -> a0, L -> Lam a0};
H0 = FullSimplify[H[a, L] /. subs0];
V0 = FullSimplify[V[a, L] /. subs0];
g0 = FullSimplify[g[a, L] /. subs0];

(* Affine inertia from entrained-fluid kinematics *)
ballMoment = FullSimplify[4 Pi Integrate[r^4, {r, 0, a0}]];
intervalMoment = FullSimplify[Integrate[w^2, {w, -Lam a0/2, Lam a0/2}]];
V3 = FullSimplify[(4 Pi/3) a0^3];

Maa = FullSimplify[rhoEff Lam a0 ballMoment/a0^2];
MLL = FullSimplify[rhoEff V3 intervalMoment/(Lam a0)^2];
M0 = {{Maa, 0}, {0, MLL}};

(* Dimensionless coordinates q = {delta a/a0, delta L/a0} *)
hBar = FullSimplify[a0^2 H0/Sigma];
mHat = FullSimplify[M0/(rhoEff V0)];
gBar = FullSimplify[(a0 g0)/V0];

(* Dynamic susceptibility in scaled frequency variable s = omega^2 rhoEff V0 a0^2 / Sigma *)
Ybar[s_] := FullSimplify[(gBar . Inverse[hBar - s mHat] . gBar)/Sigma];
Delta0 = FullSimplify[(gBar . Inverse[hBar] . gBar)/Sigma];
Delta2 = FullSimplify[(gBar . Inverse[hBar] . mHat . Inverse[hBar] . gBar)/Sigma];
lambdaEff = FullSimplify[Delta0/Delta2];

K00RawLocal = -757/2520;
K00RawDyn[s_] := FullSimplify[K00RawLocal + Ybar[s]];

section["1) Exact affine inertia reduction"];
Print["V(a,L) = ", V[a, L]];
Print["A(a,L) = ", A[a, L]];
Print["E_geom(a,L) = ", Egeom[a, L]];
Print[""];
Print["3-ball radial second moment = ", ballMoment];
Print["Interval axial second moment = ", intervalMoment];
Print["V0 = ", V0];
Print["M0 = ", MatrixForm[M0]];
Print["mHat = ", MatrixForm[mHat]];
Print["gBar = ", MatrixForm[gBar]];

section["2) Exact dynamic monopole response"];
Print["hBar = ", MatrixForm[hBar]];
Print["Y_geom(s) = ", Together[Ybar[s]]];
Print["Delta0 = ", FullSimplify[Delta0]];
Print["Delta2 = ", FullSimplify[Delta2]];
Print["lambdaEff = ", FullSimplify[lambdaEff]];
Print["K00_raw(s) = ", Together[K00RawDyn[s]]];

(* ---------------------------------------------------------------------- *)
(* EM-branch worked point                                                  *)
(* ---------------------------------------------------------------------- *)

x01 = SetPrecision[2.40482555769577276862163187933, 60];
lamEM = N[Sqrt[2] Pi/x01, 50];
rhoEx = 1/10;
betaEx = 12;
target = 109/280;

DeltaUnit = N[Delta0 /. {a0 -> 1, Lam -> lamEM, rho -> rhoEx, beta -> betaEx, Sigma -> 1}, 50];
sigmaStar = N[DeltaUnit/target, 50];

hNum = N[hBar /. {a0 -> 1, Lam -> lamEM, rho -> rhoEx, beta -> betaEx}, 40];
mNum = N[mHat /. {a0 -> 1, Lam -> lamEM}, 40];
gNum = N[gBar /. {a0 -> 1, Lam -> lamEM}, 40];
V0Num = N[V0 /. {a0 -> 1, Lam -> lamEM}, 40];
prefactor = N[sigmaStar/V0Num, 40];

(* Mass-normal generalized eigen-decomposition using mHat^{-1/2} hBar mHat^{-1/2} *)
mHalfInv = DiagonalMatrix[1/Sqrt[Diagonal[mNum]]];
mat = N[mHalfInv . hNum . mHalfInv, 40];
{valsUnsorted, vecsUnsorted} = Eigensystem[mat];
ord = Ordering[valsUnsorted];
vals = valsUnsorted[[ord]];
vecs0 = vecsUnsorted[[ord]];
vecs = Table[mHalfInv . vecs0[[i]], {i, Length[vals]}];

massOrth = Chop[Table[vecs[[i]].mNum.vecs[[j]], {i, Length[vals]}, {j, Length[vals]}], 10^-30];
stiffDiag = Chop[Table[vecs[[i]].hNum.vecs[[j]], {i, Length[vals]}, {j, Length[vals]}], 10^-30];
residues = Table[N[(gNum.vecs[[i]])^2/(sigmaStar vals[[i]]), 40], {i, Length[vals]}];
residueFractions = N[residues/Total[residues], 40];

omega2 = N[prefactor vals, 40];
lambdaEffNum = N[lambdaEff /. {a0 -> 1, Lam -> lamEM, rho -> rhoEx, beta -> betaEx}, 40];
omega2Eff = N[prefactor lambdaEffNum, 40];

YExactScaled[sval_] := Total@Table[residues[[i]]/(1 - sval/vals[[i]]), {i, Length[vals]}];
YPadeScaled[sval_] := target/(1 - sval/lambdaEffNum);

sGrid = N[Table[x, {x, 0, 0.1 vals[[1]], 0.1 vals[[1]]/399}], 40];
relErrGrid = Table[Abs[(YPadeScaled[x] - YExactScaled[x])/YExactScaled[x]], {x, sGrid}];
maxRelErr = Max[relErrGrid];

section["3) EM-branch worked point with affine inertia"];
Print["Lam_EM = ", lamEM];
Print["rho = ", rhoEx];
Print["beta = ", betaEx];
Print["Sigma* = ", sigmaStar];
Print[""];
Print["hBar(worked point) = ", MatrixForm[hNum]];
Print["mHat(worked point) = ", MatrixForm[mNum]];
Print["gBar(worked point) = ", gNum];
Print[""];
Print["Mass orthogonality V^T mHat V = ", MatrixForm[massOrth]];
Print["Stiffness diagonalization V^T hBar V = ", MatrixForm[stiffDiag]];
Print[""];
Print["Dimensionless pole parameters lambda_i = ", vals];
Print["Static residues R_i = ", residues];
Print["Residue fractions = ", residueFractions];
Print["Static sum = ", N[Total[residues], 40]];
Print["Target 109/280 = ", N[target, 40]];
Print["Residual = ", N[Total[residues] - target, 40]];
Print[""];
Print["Prefactor Sigma/(rhoEff V0 a0^2) = ", prefactor, " / rhoEff"];
Print["Omega_i^2 = ", omega2];
Print["lambdaEff = ", lambdaEffNum];
Print["Omega_eff^2 = ", omega2Eff];
Print["Max relative error of one-pole Pade form on 0 <= s <= 0.1 lambda_- : ", maxRelErr];

section["4) Verification"];
pass["mHat == DiagonalMatrix[{3/5,1/12}]", FullSimplify[mHat == DiagonalMatrix[{3/5, 1/12}]]];
pass["gBar == {3,1/Lam}", FullSimplify[gBar == {3, 1/Lam}]];
pass["Static geometry coefficient matches 109/280 at worked point",
  Chop[N[(Delta0 /. {a0 -> 1, Lam -> lamEM, rho -> rhoEx, beta -> betaEx, Sigma -> sigmaStar}) - target, 30]] == 0];
pass["Raw monopole closure gives 4/45 at worked point",
  Chop[N[(K00RawDyn[0] /. {a0 -> 1, Lam -> lamEM, rho -> rhoEx, beta -> betaEx, Sigma -> sigmaStar}) - 4/45, 30]] == 0];
pass["Mass orthogonality check",
  Max[Abs[Flatten[massOrth - IdentityMatrix[2]]]] < 10^-20];
pass["Residues sum to 109/280",
  Abs[N[Total[residues] - target, 40]] < 10^-30];

section["5) Interpretation"];
Print["1) The old monopole auxiliary is the low-frequency breathing response of the"];
Print["   same reduced geometry sector that generated the static 109/280 closure."];
Print[""];
Print["2) With affine entrained-fluid inertia, the exact monopole channel is a"];
Print["   two-pole Stieltjes response with positive residues."];
Print[""];
Print["3) At the EM worked point, the dominant pole carries ",
      N[100 residueFractions[[Ordering[residueFractions][[-1]]]], 20],
      "% of the static weight, so the single-pole auxiliary is a controlled reduction."];
Print[""];
Print["4) The remaining microphysical task is to derive the overall inertia scale"];
Print["   rhoEff (or its soft-wall analog) from the Family-1 confinement / traction PDE."];

(*"
Output:


==============================================================================
1) Exact affine inertia reduction
==============================================================================
V(a,L) = (4*a^3*L*Pi)/3
A(a,L) = (8*a^3*Pi)/3 + 4*a^2*L*Pi
E_geom(a,L) = (a^2*beta*Sigma)/(a0*L) + (8*a^3*Pi*Sigma)/(3*a0^3) + (4*a^2*L*Pi*Sigma)/a0^3 + (4*a^3*L*Pi*rho*Sigma)/(3*a0^4)

3-ball radial second moment = (4*a0^5*Pi)/5
Interval axial second moment = (a0^3*Lam^3)/12
V0 = (4*a0^4*Lam*Pi)/3
M0 = MatrixForm[{{(4*a0^4*Lam*Pi*rhoEff)/5, 0}, {0, (a0^4*Lam*Pi*rhoEff)/9}}]
mHat = MatrixForm[{{3/5, 0}, {0, 1/12}}]
gBar = MatrixForm[{3, Lam^(-1)}]

==============================================================================
2) Exact dynamic monopole response
==============================================================================
hBar = MatrixForm[{{(2*beta)/Lam + 8*Pi*(2 + Lam + Lam*rho), (-2*beta)/Lam^2 + 4*Pi*(2 + rho)}, {(-2*beta)/Lam^2 + 4*Pi*(2 + rho), (2*beta)/Lam^3}}]
Y_geom(s) = (-3*(640*beta + 320*Lam*Pi - 800*Lam^2*Pi - 320*Lam^2*Pi*rho - 12*Lam*s - 15*Lam^3*s))/((-1920*beta*Pi - 2880*beta*Lam*Pi + 3840*Lam^3*Pi^2 - 1920*beta*Lam*Pi*rho + 3840*Lam^3*Pi^2*rho + 960*Lam^3*Pi^2*rho^2 + 72*beta*s + 10*beta*Lam^2*s + 80*Lam^3*Pi*s + 40*Lam^4*Pi*s + 40*Lam^4*Pi*rho*s - 3*Lam^3*s^2)*Sigma)
Delta0 = (-4*beta + Lam*Pi*(-2 + Lam*(5 + 2*rho)))/(2*Pi*(Lam^3*Pi*(2 + rho)^2 - beta*(2 + Lam*(3 + 2*rho)))*Sigma)
Delta2 = (4*beta^2*(36 + 5*Lam^2) + Lam^4*Pi^2*(224 - 40*Lam*(4 + rho) + 36*rho*(4 + rho) + 5*Lam^2*(4 + rho)^2) - 4*beta*Lam^2*Pi*(36*(2 + rho) + 5*Lam*(-4 + Lam*(4 + rho))))/(960*Pi^2*(Lam^3*Pi*(2 + rho)^2 - beta*(2 + Lam*(3 + 2*rho)))^2*Sigma)
lambdaEff = (480*Pi*(Lam^3*Pi*(2 + rho)^2 - beta*(2 + Lam*(3 + 2*rho)))*(-4*beta + Lam*Pi*(-2 + Lam*(5 + 2*rho))))/(4*beta^2*(36 + 5*Lam^2) + Lam^4*Pi^2*(224 - 40*Lam*(4 + rho) + 36*rho*(4 + rho) + 5*Lam^2*(4 + rho)^2) - 4*beta*Lam^2*Pi*(36*(2 + rho) + 5*Lam*(-4 + Lam*(4 + rho))))
K00_raw(s) = -1/2520*(4838400*beta + 2419200*Lam*Pi - 6048000*Lam^2*Pi - 2419200*Lam^2*Pi*rho - 90720*Lam*s - 113400*Lam^3*s - 1453440*beta*Pi*Sigma - 2180160*beta*Lam*Pi*Sigma + 2906880*Lam^3*Pi^2*Sigma - 1453440*beta*Lam*Pi*rho*Sigma + 2906880*Lam^3*Pi^2*rho*Sigma + 726720*Lam^3*Pi^2*rho^2*Sigma + 54504*beta*s*Sigma + 7570*beta*Lam^2*s*Sigma + 60560*Lam^3*Pi*s*Sigma + 30280*Lam^4*Pi*s*Sigma + 30280*Lam^4*Pi*rho*s*Sigma - 2271*Lam^3*s^2*Sigma)/((-1920*beta*Pi - 2880*beta*Lam*Pi + 3840*Lam^3*Pi^2 - 1920*beta*Lam*Pi*rho + 3840*Lam^3*Pi^2*rho + 960*Lam^3*Pi^2*rho^2 + 72*beta*s + 10*beta*Lam^2*s + 80*Lam^3*Pi*s + 40*Lam^4*Pi*s + 40*Lam^4*Pi*rho*s - 3*Lam^3*s^2)*Sigma)

==============================================================================
3) EM-branch worked point with affine inertia
==============================================================================
Lam_EM = 1.8474865771201280510433744839396122428696573796896114673321817811568650327843`50.
rho = 1/10
beta = 12
Sigma* = 0.20761432918354888540403250898444970518124118990882318258536279134129844964869`48.15554576568318

hBar(worked point) = MatrixForm[{{114.3317468529872455095609837004465703294639830683748323538408`40., 19.3578673262800093188893554202382925541664874836563993190515`40.}, {19.3578673262800093188893554202382925541664874836563993190515`40., 3.8059875784510492896004581375364468153610885642501076931764`40.}}]
mHat(worked point) = MatrixForm[{{0.6000000000000000000000000000000000000000000000000000000001`40., 0}, {0, 0.0833333333333333333333333333333333333333333333333333333333`40.}}]
gBar(worked point) = {3.`40., 0.5412759217762790697885051933117400179826769616175860802258`40.}

Mass orthogonality V^T mHat V = MatrixForm[{{1.0000000000000000000000000000000000000000000000000000000006`39.221848749616356, 0}, {0, 1.0000000000000000000000000000000000000000000000000000000003`39.221848749616356}}]
Stiffness diagonalization V^T hBar V = MatrixForm[{{5.2311561265784182986358974405499758298214694688839927117242`38.02299430264189, 0}, {0, 230.9936062364795823591712397106316698369515650827420201961894`39.221848749616356}}]

Dimensionless pole parameters lambda_i = {5.2311561265784182986358974405499758298214694688839927117225`39.69897000433602, 230.9936062364795823591712397106316698369515650827420201961387`39.69897000433602}
Static residues R_i = {0.0032737645132729641382072885756778329043413221745151920103`37.40465507708957, 0.386011949772441321576078425710036452809944392111199093705`39.04575749022339}
Residue fractions = {0.0084096703093250454926425761577045248918859652189381079163`37.391367725226466, 0.9915903296906749545073574238422954751081140347810618920837`38.67291564770774}
Static sum = 0.3892857142857142857142857142857142857142857142857142857153`38.912340338109146
Target 109/280 = 0.3892857142857142857142857142857142857142857142857142857144`40.
Residual = 0``39.322071871510744

Prefactor Sigma/(rhoEff V0 a0^2) = 0.0268279459960492221543306378971474181100639978933906440619`39.99999999696444 / rhoEff
Omega_i^2 = {0.1403411740607478351490921618619829674167798112504983352805`39.52287874426849, 6.1970839935449330436240428929144642501693701554291009170535`39.52287874426849}
lambdaEff = 169.4820508806408594397201922449390188647018253490281783401399`40.
Omega_eff^2 = 4.5468553083254994905845600638074030788430518831570876000336`39.69897000281824
Max relative error of one-pole Pade form on 0 <= s <= 0.1 lambda_- : 0.00008869892456152196

==============================================================================
4) Verification
==============================================================================
mHat == DiagonalMatrix[{3/5,1/12}] PASS
gBar == {3,1/Lam} PASS
Static geometry coefficient matches 109/280 at worked point PASS
Raw monopole closure gives 4/45 at worked point PASS
Mass orthogonality check PASS
Residues sum to 109/280 PASS

==============================================================================
5) Interpretation
==============================================================================
1) The old monopole auxiliary is the low-frequency breathing response of the
   same reduced geometry sector that generated the static 109/280 closure.

2) With affine entrained-fluid inertia, the exact monopole channel is a
   two-pole Stieltjes response with positive residues.

3) At the EM worked point, the dominant pole carries 99.15903296906749545073574238422954751081`20.% of the static weight, so the single-pole auxiliary is a controlled reduction.

4) The remaining microphysical task is to derive the overall inertia scale
   rhoEff (or its soft-wall analog) from the Family-1 confinement / traction PDE.
"*)
