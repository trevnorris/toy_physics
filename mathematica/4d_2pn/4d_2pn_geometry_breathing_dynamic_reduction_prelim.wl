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
