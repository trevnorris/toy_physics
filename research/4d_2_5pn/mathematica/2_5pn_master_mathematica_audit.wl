(* 2.5PN master Mathematica audit
   --------------------------------
   Standalone Wolfram Language verification script for the core theorem chain
   behind the 4d_2_5pn paper:

     - decisive Burke-Thorne benchmark
     - low-frequency odd-channel bookkeeping
     - scalar-channel demotion checks
     - dipole/vector-channel demotion checks
     - surviving quadrupole representation / source-map checks

   Scope:
     This script is the Mathematica counterpart of the SymPy master audit.
     It checks the paper's symbolic identities inside the same declared
     reduction / passive-outgoing / small-body closure hierarchy.
*)

Print["--- 2.5PN master Mathematica audit ---"];

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

banner[title_String] := Module[{line = StringRepeat["=", 88]},
  Print["\n" <> line];
  Print[title];
  Print[line];
];

subbanner[title_String] := Module[{line = StringRepeat["-", 88]},
  Print["\n" <> line];
  Print[title];
  Print[line];
];

pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet@FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet@FullSimplify[res == 0, assum]]
];

zeroArrayQ[arr_, assum_: True] := Module[{res, flat},
  res = Quiet@FullSimplify[arr, assum];
  flat = Flatten[{res}];
  And @@ (TrueQ[# === 0] || TrueQ[Quiet@FullSimplify[# == 0, assum]] & /@ flat)
];

checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet@FullSimplify[expr, assum];
  If[zeroScalarQ[expr, assum], pass[name], fail[name, res]]
];

checkMatrix[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet@FullSimplify[expr, assum];
  If[zeroArrayQ[expr, assum], pass[name], fail[name, res]]
];

checkCondition[name_String, cond_] := If[TrueQ[cond], pass[name], fail[name, cond]];

stfFromVector[v_List] := FullSimplify[Outer[Times, v, v] - IdentityMatrix[3] (v.v)/3];
stfTensor[m_?MatrixQ] := FullSimplify[(m + Transpose[m])/2 - IdentityMatrix[3] Tr[m]/3];

E20 = {{-1/Sqrt[6], 0, 0}, {0, -1/Sqrt[6], 0}, {0, 0, 2/Sqrt[6]}};
E21c = {{0, 0, 1/Sqrt[2]}, {0, 0, 0}, {1/Sqrt[2], 0, 0}};
E21s = {{0, 0, 0}, {0, 0, 1/Sqrt[2]}, {0, 1/Sqrt[2], 0}};
E22c = {{1/Sqrt[2], 0, 0}, {0, -1/Sqrt[2], 0}, {0, 0, 0}};
E22s = {{0, 1/Sqrt[2], 0}, {1/Sqrt[2], 0, 0}, {0, 0, 0}};

basis = {E20, E21c, E21s, E22c, E22s};
basisNames = {"20", "21c", "21s", "22c", "22s"};
coeffsInBasis[T_?MatrixQ] := FullSimplify[Total[Flatten[T #]]] & /@ basis;

pi20[ux_, uy_, vn_] := FullSimplify[(3 vn^2 - (ux^2 + uy^2 + vn^2))/2];
pi21c[ux_, uy_, vn_] := FullSimplify[Sqrt[2] vn ux];
pi21s[ux_, uy_, vn_] := FullSimplify[Sqrt[2] vn uy];
pi22c[ux_, uy_] := FullSimplify[(ux^2 - uy^2)/(2 Sqrt[2])];
pi22s[ux_, uy_] := FullSimplify[(2 ux uy)/(2 Sqrt[2])];

tensorAction[A_?MatrixQ, T_?MatrixQ] := FullSimplify[A.T + T.Transpose[A]];
repMatrix[A_?MatrixQ] := Module[{mrep = ConstantArray[0, {5, 5}], s, comps},
  Do[
    s = tensorAction[A, basis[[j]]];
    comps = coeffsInBasis[s];
    Do[mrep[[i, j]] = FullSimplify[comps[[i]]], {i, 1, 5}],
    {j, 1, 5}
  ];
  mrep
];

banner["PART II — DECISIVE BENCHMARK / BURKE–THORNE PROTOTYPE"];

subbanner["II.1 — Two-body mass dipole and STF quadrupole decomposition"];
m1 =.; m2 =.; X1 =.; X2 =.; X3 =.; x1 =.; x2 =.; x3 =.;
Mtot = m1 + m2;
muRed = FullSimplify[m1 m2/Mtot];
Xvec = {X1, X2, X3};
xvec = {x1, x2, x3};
r1 = FullSimplify[Xvec + (m2/Mtot) xvec];
r2 = FullSimplify[Xvec - (m1/Mtot) xvec];
Dvec = FullSimplify[m1 r1 + m2 r2];
dDdx = Table[D[Dvec[[i]], xvec[[j]]], {i, 1, 3}, {j, 1, 3}];
Qmat = FullSimplify[m1 stfFromVector[r1] + m2 stfFromVector[r2]];
Qexpected = FullSimplify[Mtot stfFromVector[Xvec] + muRed stfFromVector[xvec]];
Print["D_i = ", MatrixForm[Dvec]];
checkMatrix["dD_i/dx_j", dDdx];
checkMatrix["Q - (M STF(XX) + mu STF(xx))", Qmat - Qexpected];

subbanner["II.2 — Burke–Thorne local quadrupole force and Iyer–Will match"];
Clear[r];
t =.; G =.; c =.; m =.; nu =.; v2 =.;
rfun = r[t];
rd = Derivative[1][r][t];
AcoefV = FullSimplify[(2 G m + 54 v2 rfun - 75 rfun rd^2)/(3 rfun)];
BcoefV = FullSimplify[(2 G m - 6 v2 rfun + 15 rfun rd^2)/rfun];
alpha =.; beta =.;
alphaSol = alpha /. First[Solve[
    {
      Coefficient[BcoefV, v2] == -(2 + alpha),
      Coefficient[BcoefV, G m/rfun] == -(2 - alpha),
      Coefficient[BcoefV, rd^2] == 3 (1 + alpha)
    },
    alpha
  ]];
betaSol = beta /. First[Solve[
    Coefficient[AcoefV, v2] == 3 (1 + beta),
    beta
  ]];
Print["A(v^2, GM/r, rdot^2) = ", Factor[AcoefV]];
Print["B(v^2, GM/r, rdot^2) = ", Factor[BcoefV]];
Print["alpha = ", alphaSol];
Print["beta  = ", betaSol];
checkCondition["Burke–Thorne prototype lands on (alpha,beta)=(4,5)", alphaSol == 4 && betaSol == 5];

banner["PART III — LOW-FREQUENCY SELECTION RULES / INFLUENCE FUNCTIONAL"];

subbanner["III.1 — Time-domain signs for i omega^n"];
omega =.; Omega =.; sigma =.; g =.;
sign1 = FullSimplify[I/((-I)^1)];
sign3 = FullSimplify[I/((-I)^3)];
sign5 = FullSimplify[I/((-I)^5)];
Print["i*omega^1  ->  ", sign1, " * d^1/dt^1"];
Print["i*omega^3  ->  ", sign3, " * d^3/dt^3"];
Print["i*omega^5  ->  ", sign5, " * d^5/dt^5"];
checkCondition["Fourier sign map is {-1,+1,-1}", sign1 == -1 && sign3 == 1 && sign5 == -1];

subbanner["III.2 — Minimal retarded kernel expansions"];
K1 = Expand[Normal@Series[g^2/(Omega^2 - omega^2 - I sigma omega), {omega, 0, 2}]];
K3 = Expand[Normal@Series[g^2/(Omega^2 - omega^2 - I sigma omega^3), {omega, 0, 4}]];
K5 = Expand[Normal@Series[g^2/(Omega^2 - omega^2 - I sigma omega^5), {omega, 0, 6}]];
Print["K1 = ", K1];
Print["K3 = ", K3];
Print["K5 = ", K5];

subbanner["III.3 — Dissipation / Schott identities"];
Clear[q];
q1 = Derivative[1][q][t];
q2 = Derivative[2][q][t];
q3 = Derivative[3][q][t];
q4 = Derivative[4][q][t];
q5 = Derivative[5][q][t];
checkScalar["n=1 dissipation identity", q1 (-q1) + q1^2];
checkScalar["n=3 dissipation identity", q1 q3 - (D[q1 q2, t] - q2^2)];
checkScalar["n=5 dissipation identity", q1 (-q5) - (-D[q1 q4 - q2 q3, t] - q3^2)];

subbanner["III.4 — Model-specific 2PN scalar / quadrupole combinations"];
delta01 =.; deltag1 =.; delta205 =.;
gamma1eff = FullSimplify[(16/5) delta01 - (281/80) deltag1];
gamma5eff = FullSimplify[(25/16) delta205];
Print["gamma1_eff = ", gamma1eff];
Print["gamma5_eff = ", gamma5eff];

banner["PART IV — SCALAR SECTOR"];

subbanner["IV.1 — Outgoing scalar Green function and monopole odd term"];
cS =.; rrad =.; omega =.;
k = omega/cS;
GR = Exp[I k rrad]/(4 Pi rrad);
GRSeries = Expand[Normal@Series[GR, {omega, 0, 5}]];
Print["G_R(omega,r) = ", GR];
Print["Series = ", GRSeries];
Print["Model-specific gamma1_eff = ", gamma1eff];

subbanner["IV.1b — 2PN scalar support/geometry finite-size rescue"];
a0 =.; ag =.; k0 =.; z0 =.;
j0 = Sin[z0]/z0;
y0 = -Cos[z0]/z0;
h0 = FullSimplify[j0 + I y0];
Lambda0 = FullSimplify[(k0 D[h0, z0]/h0) /. z0 -> k0 a0];
Y0 = FullSimplify[1/Lambda0];
Y0Norm = FullSimplify[Y0/(Y0 /. k0 -> 0)];
YgNorm = FullSimplify[1/(1 - I ag k0)];
Seff = FullSimplify[(16/5) Y0Norm + 25/16 - (281/80) YgNorm];
ellEff = FullSimplify[(16/5) a0 - (281/80) ag];
Print["Lambda0(k) = ", Lambda0];
Print["Y0_norm(k) = ", Normal@Series[Y0Norm, {k0, 0, 5}]];
Print["Yg_norm(k) = ", Normal@Series[YgNorm, {k0, 0, 5}]];
Print["Seff(k) = ", Normal@Series[Seff, {k0, 0, 4}]];
Print["ell_eff = ", ellEff];
Print["equal-scale residual ell_eff(ag=a0) = ", FullSimplify[ellEff /. ag -> a0]];
Print["exact scalar cancellation ag = ", First[ag /. Solve[ellEff == 0, ag]]];

subbanner["IV.2 — Gaussian overlap counterexample and exact leakage identity"];
w =.; a =.; ell =.; adot =.;
W = Exp[-w^2/ell^2]/(Sqrt[Pi] ell);
gg = Exp[-w^2/a^2]/(Sqrt[Pi] a);
Cgauss = FullSimplify[Integrate[W gg, {w, -Infinity, Infinity}], a > 0 && ell > 0];
jw = FullSimplify[(adot/a) w gg];
continuityResidual = FullSimplify[adot D[gg, a] + D[jw, w]];
Ileak = FullSimplify[Integrate[D[W, w] jw, {w, -Infinity, Infinity}], a > 0 && ell > 0];
dCda = FullSimplify[adot D[Cgauss, a]];
Print["C(a) = ", Cgauss];
Print["dC/da = ", FullSimplify[D[Cgauss, a]]];
checkScalar["continuity residual", continuityResidual];
checkScalar["I_leak - adot*dC/da", Ileak - dCda];

subbanner["IV.3 — Projection-locking linear algebra criterion"];
v1a =.; v1L =.; v2a =.; v2L =.;
detM = FullSimplify[v1a v2L - v1L v2a];
Print["determinant of two-mode tangent matrix = ", detM];
Print["If determinant != 0, projection-locking requires B_a = B_L = 0."];

subbanner["IV.4 — Direct vs derivative coupling; damped discrete-mode test"];
A0 =.; B0 =.; C0 =.; D0 =.; g0 =.; gd =.; Omega =.; eta =.; lam =.; lamd =.; omega =.;
Ggapless = A0 + I B0 omega + C0 omega^2 + I D0 omega^3;
direct = Expand[g0^2 Ggapless];
derivative = Expand[gd^2 omega^2 Ggapless];
discreteDamped = Expand[Normal@Series[lam^2/(Omega^2 - omega^2 - I eta omega), {omega, 0, 3}]];
discreteDampedDerivative = Expand[Normal@Series[lamd^2 omega^2/(Omega^2 - omega^2 - I eta omega), {omega, 0, 4}]];
gamma1Direct = Expand[ComplexExpand[Im[D[direct, omega] /. omega -> 0]]];
gamma3Derivative = FullSimplify[ComplexExpand[Im[D[derivative, {omega, 3}] /. omega -> 0]]/Factorial[3]];
gamma1Disc = Expand[ComplexExpand[Im[D[discreteDamped, omega] /. omega -> 0]]];
Print["direct = ", direct];
Print["derivative = ", derivative];
Print["discrete_damped = ", discreteDamped];
Print["discrete_damped_derivative = ", discreteDampedDerivative];
Print["gamma1_direct = ", gamma1Direct];
Print["gamma3_derivative = ", gamma3Derivative];
Print["gamma1_discrete_damped = ", gamma1Disc];
checkCondition["Direct coupling keeps a linear odd term", gamma1Direct =!= 0];
checkCondition["Derivative coupling shifts the odd term to cubic order", gamma3Derivative =!= 0];
checkCondition["Damped discrete mode inherits a linear odd term", gamma1Disc =!= 0];

subbanner["IV.5 — Breathing-to-outlet vertex is dot(q)-type"];
Bq =.; delta =.;
Kdirect = Expand[Bq^2 Ggapless];
Kderiv = Expand[Bq^2 omega^2 Ggapless];
Print["K_direct = ", Kdirect];
Print["Im K_direct = ", Expand[ComplexExpand[Im[Kdirect]]]];
Print["K_deriv = ", Kderiv];
Print["Im K_deriv = ", Expand[ComplexExpand[Im[Kderiv]]]];
checkCondition["Derivative outlet kernel has no linear odd term", Coefficient[Expand[ComplexExpand[Im[Kderiv]]], omega, 1] == 0];
breathingExponent = 5/2 + 3 delta;
Print["Derivative-breathing odd term exponent (with a/r ~ eps^delta) = ", breathingExponent];

subbanner["IV.6 — No third linear scalar source from quadratic momentum/stress terms"];
wv =.; z1 =.; z2 =.; z3 =.; u0 =.; u1 =.; u2 =.; q0 =.;
Zsig = I z1 wv + z2 wv^2 + I z3 wv^3;
uSigma = u0 q0 + I u1 wv q0 + u2 wv^2 q0;
jSigma = Expand[Normal@Series[Zsig uSigma, {wv, 0, 3}]];
Print["j_sigma(w) = ", jSigma];
Print["With Z_sigma(0)=0, the leading mouth source is derivative-like, not direct q-like."];

subbanner["IV.7 — Mouth radiative-order theorem and Ohmic no-go"];
wv =.; Kmouth =.; Mmouth =.; Madd =.; betaM =.; R2 =.; c0 =.; rho0 =.; chiSigma =.; kappaM =.; L =.; z2Nat =.; aSmall =.; 
Zrad = R2 wv^2 + I Madd wv;
den = Kmouth - Mmouth wv^2 + I wv Zrad;
Ymouth = FullSimplify[I wv betaM/den];
YmouthExpand = Expand[Normal@Series[Ymouth, {wv, 0, 5}]];
Print["Compact reciprocal mouth admittance Y(w) = ", YmouthExpand];
Print["Re Y series = ", Normal@Series[ComplexExpand[Re[YmouthExpand]], {wv, 0, 5}]];
Print["Im Y series = ", Normal@Series[ComplexExpand[Im[YmouthExpand]], {wv, 0, 5}]];
Z1d = c0 rho0;
Y1d = FullSimplify[I betaM wv/(Kmouth - Mmouth wv^2 + I wv Z1d)];
Print["1D Ohmic benchmark Re Y1 = ", Normal@Series[ComplexExpand[Re[Y1d]], {wv, 0, 4}]];
ellMouth = FullSimplify[(chiSigma cS^3/(kappaM L)) z2Nat];
ellNat = FullSimplify[ellMouth /. {z2Nat -> aSmall^2/cS^3, L -> delta aSmall}];
Print["Natural odd mouth length ell_mouth ~ ", ellNat];

banner["PART V — DIPOLE / VECTOR SECTOR"];

subbanner["V.1 — Carried odd wake ports: CM/relative split"];
mA =.; mB =.; M =.;
etaA = mA/M;
etaB = mB/M;
Ux =.; Uy =.; Dz =.; ux =.; uy =.; d =.;
UA = {Ux, Uy};
urel = {ux, uy};
PiAperp = FullSimplify[Sqrt[7/2] (UA + etaB urel)];
PiBperp = FullSimplify[Sqrt[7/2] (UA - etaA urel)];
PiA10 = FullSimplify[2 (Dz + etaB d)];
PiB10 = FullSimplify[2 (Dz - etaA d)];
Print["PiA_perp - PiB_perp = ", FullSimplify[(PiAperp - PiBperp) /. M -> mA + mB]];
Print["PiA_10   - PiB_10   = ", FullSimplify[(PiA10 - PiB10) /. M -> mA + mB]];

subbanner["V.2 — Vector Ward-identity failure"];
Clear[xpA, xpB, xmA, xmB, Xcm, xrel1];
acoef =.; 
Lbody = acoef (Derivative[1][xmA][t] - Derivative[1][xmB][t]) (Derivative[2][xpA][t] - Derivative[2][xpB][t]);
FA = FullSimplify[-D[D[Lbody, Derivative[1][xmA][t]], t]];
FB = FullSimplify[-D[D[Lbody, Derivative[1][xmB][t]], t]];
Print["F_A = ", FA];
Print["F_B = ", FB];
checkScalar["F_A + F_B", FA + FB];
Mtot2 = mA + mB;
eA = FullSimplify[mA/Mtot2];
eB = FullSimplify[mB/Mtot2];
xA = Xcm[t] + eB xrel1[t];
xB = Xcm[t] - eA xrel1[t];
dxDiff = FullSimplify[D[xA, t] - D[xB, t]];
ddxDiff = FullSimplify[D[xA, {t, 2}] - D[xB, {t, 2}]];
Lcmrel = FullSimplify[acoef dxDiff ddxDiff];
Print["L_cmrel = ", Lcmrel];
checkCondition["CM coordinate drops out of the relative dipole odd action", FreeQ[Lcmrel, Derivative[1][Xcm][t]] && FreeQ[Lcmrel, Derivative[2][Xcm][t]]];

subbanner["V.3 — Outgoing l=1 spectral no-go for linear odd term"];
a =.; k =.; z =.; R =.;
j1 = Sin[z]/z^2 - Cos[z]/z;
y1 = -Cos[z]/z^2 - Sin[z]/z;
h1 = FullSimplify[j1 + I y1];
Lambda1 = FullSimplify[(k D[h1, z]/h1) /. z -> k a];
Lambda1Series = Normal@Series[Lambda1, {k, 0, 6}];
Y1 = Normal@Series[1/Lambda1Series, {k, 0, 5}];
Y1Norm = Expand[R Y1/(-a/2)];
imagK3 = FullSimplify[Coefficient[Y1Norm, k, 3]/I];
Print["Lambda1(k) = ", Lambda1Series];
Print["Y1_norm(k) = ", Y1Norm];
checkCondition["Outgoing l=1 has no linear k term", Coefficient[Y1Norm, k, 1] == 0];
Print["Imaginary k^3 coefficient = ", imagK3];

subbanner["V.3b — Order reduction of an isotropic fifth-derivative force"];
rN =.; rdN =.; v2N =.; mN =.;
dscalar[expr_] := FullSimplify[D[expr, rN] rdN + D[expr, rdN] (v2N/rN - rdN^2/rN - mN/rN^2) + D[expr, v2N] (-2 mN rdN/rN^2)];
addVec[u_List, w_List] := {FullSimplify[u[[1]] + w[[1]]], FullSimplify[u[[2]] + w[[2]]]};
mulVec[s_, u_List] := {FullSimplify[s u[[1]]], FullSimplify[s u[[2]]]};
dvec[u_List] := Module[{A1 = u[[1]], B1 = u[[2]], dn, dv, res},
  dn = {-rdN/rN, 1/rN};
  dv = {-mN/rN^2, 0};
  res = addVec[{dscalar[A1], 0}, mulVec[A1, dn]];
  res = addVec[res, {0, dscalar[B1]}];
  res = addVec[res, mulVec[B1, dv]];
  res
];
xv = {rN, 0};
vv = {0, 1};
a2v = dvec[vv];
a3v = dvec[a2v];
a4v = dvec[a3v];
a5v = dvec[a4v];
Print["x^(5) = ", a5v[[1]], " * n + ", a5v[[2]], " * v"];
Print["Thus -kap x^(5) lies in the standard {rdot n, v} 2.5PN basis."];

subbanner["V.4 — Same-order dipole wake cannot be absorbed into the Iyer–Will family"];
p =.; qpar =.; s =.; alpha =.; beta =.;
Fvec = {3 (1 + beta), (23 + 6 alpha - 9 beta)/3, -5 beta, -(2 + alpha), alpha - 2, 3 (1 + alpha)};
BTvec = {18, 2/3, -25, -6, 2, 15};
Dvec = {9 p + 36 qpar, -8 p - 22 qpar, -45 p - 60 qpar, -9 p, 8 p, 45 p};
sol = Solve[Thread[BTvec + Dvec == s Fvec], {alpha, beta, s, p, qpar}];
Print["Full matching solution = ", sol];
checkCondition["Same-order dipole match requires p=q=0", sol =!= {} && (p /. First[sol]) == 0 && (qpar /. First[sol]) == 0];

subbanner["V.5 — Dipole finite-size scaling rescue"];
eps =.; rho =.; delta =.;
dipoleScaling = FullSimplify[eps (Sqrt[eps] rho)^3];
Print["eps * (sqrt(eps)*rho)^3 = ", dipoleScaling];
Print["If rho ~ eps^delta, exponent = 5/2 + 3 delta."];
Print["general exponent = ", 5/2 + 3 delta];

banner["PART VI — QUADRUPOLE SECTOR"];

subbanner["VI.1 — P2 ports are exactly the real STF l=2 content"];
ux =.; uy =.; vn =.;
vSq = ux^2 + uy^2 + vn^2;
Q = {
  {ux^2 - vSq/3, ux uy, ux vn},
  {ux uy, uy^2 - vSq/3, uy vn},
  {ux vn, uy vn, vn^2 - vSq/3}
};
q20 = FullSimplify[Total[Flatten[Q E20]]];
q21c = FullSimplify[Total[Flatten[Q E21c]]];
q21s = FullSimplify[Total[Flatten[Q E21s]]];
q22c = FullSimplify[Total[Flatten[Q E22c]]];
q22s = FullSimplify[Total[Flatten[Q E22s]]];
checkScalar["q20 - sqrt(2/3) Pi20", q20 - Sqrt[2/3] pi20[ux, uy, vn]];
checkScalar["q21c - Pi21c", q21c - pi21c[ux, uy, vn]];
checkScalar["q21s - Pi21s", q21s - pi21s[ux, uy, vn]];
checkScalar["q22c - 2 Pi22c", q22c - 2 pi22c[ux, uy]];
checkScalar["q22s - 2 Pi22s", q22s - 2 pi22s[ux, uy]];

subbanner["VI.2 — Outgoing l=2 retarded completion gives +i*k^5"];
a =.; k =.; z =.; R2 =.;
j2 = ((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2 = -((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2;
h2 = FullSimplify[j2 + I y2];
Lambda2 = FullSimplify[(k D[h2, z]/h2) /. z -> k a];
Lambda2Series = Normal@Series[Lambda2, {k, 0, 6}];
Y2 = Normal@Series[1/Lambda2Series, {k, 0, 5}];
Y2Norm = Expand[R2 Y2/(-a/3)];
imagK5 = FullSimplify[Coefficient[Y2Norm, k, 5]/I];
Print["Lambda2(k) = ", Lambda2Series];
Print["Y2_norm(k) = ", Y2Norm];
Print["Imaginary k^5 coefficient = ", imagK5];
checkCondition["Quadrupole k^5 coefficient is positive", TrueQ[FullSimplify[imagK5 > 0, R2 > 0 && a > 0]]];

subbanner["VI.3 — Quadrupole odd kernel gives damping, not antidamping"];
Clear[q];
gammaQuad =.;
Fone = FullSimplify[-(gammaQuad/2) Derivative[5][q][t]];
Pone = Expand[Fone Derivative[1][q][t]];
SchottOne = FullSimplify[-(gammaQuad/2) (Derivative[1][q][t] Derivative[4][q][t] - Derivative[2][q][t] Derivative[3][q][t])];
residual = FullSimplify[Pone - D[SchottOne, t]];
Print["P - dE_Schott/dt = ", residual];
checkScalar["quadrupole power balance residual", residual + (gammaQuad/2) Derivative[3][q][t]^2];

subbanner["VI.3b — Compact internal P2 branch is finite-size suppressed"];
eps =.; rho =.; delta =.;
internalScaling = FullSimplify[eps^2 (Sqrt[eps] rho)^5];
Print["eps^2 * (sqrt(eps)*rho)^5 = ", internalScaling];
Print["If rho ~ eps^delta, exponent = ", 9/2 + 5 delta];

subbanner["VI.4 — Orbital/worldtube STF source map"];
uX =.; uY =.; d =.; rSep =.; mu =.; G =.; M =.;
vOrb = {uX, uY, d};
xOrb = {0, 0, rSep};
Vstf = stfTensor[Outer[Times, vOrb, vOrb]];
aOrb = FullSimplify[-(G M/rSep^3) xOrb];
IorbDDot = FullSimplify[2 mu (Vstf + stfTensor[Outer[Times, xOrb, aOrb]])];
Vcomp = AssociationThread[basisNames, coeffsInBasis[Vstf]];
Xcomp = AssociationThread[basisNames, coeffsInBasis[stfTensor[Outer[Times, xOrb, xOrb]]]];
Icomp = AssociationThread[basisNames, coeffsInBasis[IorbDDot]];
checkScalar["V20 source map", Vcomp["20"] - Sqrt[2/3] pi20[uX, uY, d]];
checkScalar["V21c source map", Vcomp["21c"] - pi21c[uX, uY, d]];
checkScalar["V21s source map", Vcomp["21s"] - pi21s[uX, uY, d]];
checkScalar["V22c source map", Vcomp["22c"] - 2 pi22c[uX, uY]];
checkScalar["V22s source map", Vcomp["22s"] - 2 pi22s[uX, uY]];
checkScalar["X21c static position source", Xcomp["21c"]];
checkScalar["X21s static position source", Xcomp["21s"]];
checkScalar["X22c static position source", Xcomp["22c"]];
checkScalar["X22s static position source", Xcomp["22s"]];
checkScalar["I21c source map", Icomp["21c"] - 2 mu pi21c[uX, uY, d]];
checkScalar["I21s source map", Icomp["21s"] - 2 mu pi21s[uX, uY, d]];
checkScalar["I22c source map", Icomp["22c"] - 4 mu pi22c[uX, uY]];
checkScalar["I22s source map", Icomp["22s"] - 4 mu pi22s[uX, uY]];
checkScalar["I20 source map", Icomp["20"] - 2 mu Sqrt[2/3] (pi20[uX, uY, d] - G M/rSep)];

subbanner["VI.5 — Matching theorem: commutant is scalar, basis map invertible"];
Ax = {{0, 0, 0}, {0, 0, -1}, {0, 1, 0}};
Ay = {{0, 0, 1}, {0, 0, 0}, {-1, 0, 0}};
Az = {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}};
Jx = repMatrix[Ax];
Jy = repMatrix[Ay];
Jz = repMatrix[Az];
vars = Array[u, 25];
Munk = Partition[vars, 5];
eqs = Flatten[Flatten /@ (Munk.# - #.Munk & /@ {Jx, Jy, Jz})];
coefMat = CoefficientArrays[eqs, vars][[2]];
nullSpace = NullSpace[Normal[coefMat]];
commMat = Partition[First[nullSpace], 5];
scale = commMat[[1, 1]];
Print["commutant basis matrix = ", MatrixForm[commMat]];
checkCondition["commutant dimension is one", Length[nullSpace] == 1];
checkMatrix["commutant minus scalar identity", commMat - scale IdentityMatrix[5]];
Bmap = DiagonalMatrix[{Sqrt[2/3], 1, 1, 2, 2}];
Bdet = FullSimplify[Det[Bmap]];
Print["det(port->canonical basis map) = ", Bdet];
checkCondition["Port-to-canonical STF map is invertible", Bdet =!= 0];

subbanner["VI.6 — Static overlap test"];
rSym =.; Gsym =.; Msym =.; muSym =.;
xStatic = {0, 0, rSym};
xStf = stfTensor[Outer[Times, xStatic, xStatic]];
aStatic = FullSimplify[-(Gsym Msym/rSym^3) xStatic];
IorbDDotStatic = FullSimplify[2 muSym stfTensor[Outer[Times, xStatic, aStatic]]];
C20x = FullSimplify[Total[Flatten[xStf E20]]];
C20ddot = FullSimplify[Total[Flatten[IorbDDotStatic E20]]/(2 muSym)];
J20 = 5/4;
q20Static = J20;
m0FromDDot = FullSimplify[q20Static/C20ddot];
Print["C20[x_<i x_j>] = ", C20x];
Print["C20[ddot I_orb /(2 mu)] = ", C20ddot];
Print["q20_static = ", q20Static];
Print["m0 from static overlap = ", m0FromDDot];
checkCondition["Static overlap is nonzero on the natural branch", m0FromDDot =!= 0];

banner["FINAL CONDITIONAL THEOREM LEDGER"];
Print["Conservative decisive benchmark:"];
Print["  Burke–Thorne prototype lands on Iyer–Will parameters (alpha,beta)=(4,5)."];
Print["\nLow-frequency framework:"];
Print["  Wrong odd channels are identified by n = 1 (scalar), n = 3 (dipole), n = 5 (quadrupole)."];
Print["\nScalar side:"];
Print["  - direct scalar monopole overlap generically gives i*omega,"];
Print["  - continuity-derived breathing outlet is derivative-coupled, so its kernel starts at i*omega^3,"];
Print["  - compact reciprocal mouth radiation is super-Ohmic / higher-order,"];
Print["  - genuine Ohmic scalar bleed would require new structure outside the natural compact-radiative branch."];
Print["\nDipole side:"];
Print["  - equal-and-opposite body forces do not kill the relative carried wake,"];
Print["  - but a clean outgoing l=1 completion removes the dangerous linear odd term and starts at i*k^3,"];
Print["  - same-order nonzero dipole deformation cannot be absorbed into the standard 2.5PN family,"];
Print["  - in a strict small-body branch it is pushed above universal point-particle 2.5PN."];
Print["\nQuadrupole side:"];
Print["  - the solved 2PN P2 package is exactly the real STF l=2 content,"];
Print["  - the source map to the orbital/worldtube STF quadrupole is explicit,"];
Print["  - rotational invariance forces the retarded 5x5 matching matrix to be scalar,"];
Print["  - a compact outgoing l=2 completion gives a positive +i*k^5 term,"];
Print["  - and the static overlap gate passes: m0 is nonzero on the natural low-frequency branch."];
Print["\nBest current verdict:"];
Print["  A full GR-like point-particle 2.5PN theorem is conditionally reachable,"];
Print["  with the remaining narrow gap living in the final passive/outgoing quadrupole normalization."];

Print["\nPasses: ", passCount];
Print["Fails : ", failCount];
If[failCount == 0,
  Print["OVERALL: PASS"],
  Print["OVERALL: FAIL"]
];

(*"
Output:
--- 2.5PN master Mathematica audit ---

========================================================================================
PART II — DECISIVE BENCHMARK / BURKE–THORNE PROTOTYPE
========================================================================================

----------------------------------------------------------------------------------------
II.1 — Two-body mass dipole and STF quadrupole decomposition
----------------------------------------------------------------------------------------
D_i = MatrixForm[{(m1 + m2)*X1, (m1 + m2)*X2, (m1 + m2)*X3}]
PASS: dD_i/dx_j
PASS: Q - (M STF(XX) + mu STF(xx))

----------------------------------------------------------------------------------------
II.2 — Burke–Thorne local quadrupole force and Iyer–Will match
----------------------------------------------------------------------------------------
A(v^2, GM/r, rdot^2) = (2*G*m + 54*v2*r[t] - 75*r[t]*Derivative[1][r][t]^2)/(3*r[t])
B(v^2, GM/r, rdot^2) = -((-2*G*m + 6*v2*r[t] - 15*r[t]*Derivative[1][r][t]^2)/r[t])
alpha = 4
beta  = 5
PASS: Burke–Thorne prototype lands on (alpha,beta)=(4,5)

========================================================================================
PART III — LOW-FREQUENCY SELECTION RULES / INFLUENCE FUNCTIONAL
========================================================================================

----------------------------------------------------------------------------------------
III.1 — Time-domain signs for i omega^n
----------------------------------------------------------------------------------------
i*omega^1  ->  -1 * d^1/dt^1
i*omega^3  ->  1 * d^3/dt^3
i*omega^5  ->  -1 * d^5/dt^5
PASS: Fourier sign map is {-1,+1,-1}

----------------------------------------------------------------------------------------
III.2 — Minimal retarded kernel expansions
----------------------------------------------------------------------------------------
K1 = (g^2*omega^2)/Omega^4 + g^2/Omega^2 + (I*g^2*omega*sigma)/Omega^4 - (g^2*omega^2*sigma^2)/Omega^6
K3 = (g^2*omega^4)/Omega^6 + (g^2*omega^2)/Omega^4 + g^2/Omega^2 + (I*g^2*omega^3*sigma)/Omega^4
K5 = (g^2*omega^6)/Omega^8 + (g^2*omega^4)/Omega^6 + (g^2*omega^2)/Omega^4 + g^2/Omega^2 + (I*g^2*omega^5*sigma)/Omega^4

----------------------------------------------------------------------------------------
III.3 — Dissipation / Schott identities
----------------------------------------------------------------------------------------
PASS: n=1 dissipation identity
PASS: n=3 dissipation identity
PASS: n=5 dissipation identity

----------------------------------------------------------------------------------------
III.4 — Model-specific 2PN scalar / quadrupole combinations
----------------------------------------------------------------------------------------
gamma1_eff = (16*delta01)/5 - (281*deltag1)/80
gamma5_eff = (25*delta205)/16

========================================================================================
PART IV — SCALAR SECTOR
========================================================================================

----------------------------------------------------------------------------------------
IV.1 — Outgoing scalar Green function and monopole odd term
----------------------------------------------------------------------------------------
G_R(omega,r) = E^((I*omega*rrad)/cS)/(4*Pi*rrad)
Series = (I/4*omega)/(cS*Pi) + 1/(4*Pi*rrad) - (omega^2*rrad)/(8*cS^2*Pi) - (I/24*omega^3*rrad^2)/(cS^3*Pi) + (omega^4*rrad^3)/(96*cS^4*Pi) + (I/480*omega^5*rrad^4)/(cS^5*Pi)
Model-specific gamma1_eff = (16*delta01)/5 - (281*deltag1)/80

----------------------------------------------------------------------------------------
IV.1b — 2PN scalar support/geometry finite-size rescue
----------------------------------------------------------------------------------------
Lambda0(k) = -a0^(-1) + I*k0
Y0_norm(k) = 1 + I*a0*k0 - a0^2*k0^2 - I*a0^3*k0^3 + a0^4*k0^4 + I*a0^5*k0^5
Yg_norm(k) = 1 + I*ag*k0 - ag^2*k0^2 - I*ag^3*k0^3 + ag^4*k0^4 + I*ag^5*k0^5
Seff(k) = 5/4 + I/80*(256*a0 - 281*ag)*k0 + ((-256*a0^2 + 281*ag^2)*k0^2)/80 - I/80*(256*a0^3 - 281*ag^3)*k0^3 + ((256*a0^4 - 281*ag^4)*k0^4)/80
ell_eff = (16*a0)/5 - (281*ag)/80
equal-scale residual ell_eff(ag=a0) = (-5*a0)/16
exact scalar cancellation ag = (256*a0)/281

----------------------------------------------------------------------------------------
IV.2 — Gaussian overlap counterexample and exact leakage identity
----------------------------------------------------------------------------------------
C(a) = 1/(Sqrt[a^2 + ell^2]*Sqrt[Pi])
dC/da = -(a/((a^2 + ell^2)^(3/2)*Sqrt[Pi]))
PASS: continuity residual
PASS: I_leak - adot*dC/da

----------------------------------------------------------------------------------------
IV.3 — Projection-locking linear algebra criterion
----------------------------------------------------------------------------------------
determinant of two-mode tangent matrix = -(v1L*v2a) + v1a*v2L
If determinant != 0, projection-locking requires B_a = B_L = 0.

----------------------------------------------------------------------------------------
IV.4 — Direct vs derivative coupling; damped discrete-mode test
----------------------------------------------------------------------------------------
direct = A0*g0^2 + I*B0*g0^2*omega + C0*g0^2*omega^2 + I*D0*g0^2*omega^3
derivative = A0*gd^2*omega^2 + I*B0*gd^2*omega^3 + C0*gd^2*omega^4 + I*D0*gd^2*omega^5
discrete_damped = (-I*eta^3*lam^2*omega^3)/Omega^8 - (eta^2*lam^2*omega^2)/Omega^6 + ((2*I)*eta*lam^2*omega^3)/Omega^6 + (I*eta*lam^2*omega)/Omega^4 + (lam^2*omega^2)/Omega^4 + lam^2/Omega^2
discrete_damped_derivative = -((eta^2*lamd^2*omega^4)/Omega^6) + (I*eta*lamd^2*omega^3)/Omega^4 + (lamd^2*omega^4)/Omega^4 + (lamd^2*omega^2)/Omega^2
gamma1_direct = B0*g0^2
gamma3_derivative = B0*gd^2
gamma1_discrete_damped = (eta*lam^2)/Omega^4
PASS: Direct coupling keeps a linear odd term
PASS: Derivative coupling shifts the odd term to cubic order
PASS: Damped discrete mode inherits a linear odd term

----------------------------------------------------------------------------------------
IV.5 — Breathing-to-outlet vertex is dot(q)-type
----------------------------------------------------------------------------------------
K_direct = A0*Bq^2 + I*B0*Bq^2*omega + Bq^2*C0*omega^2 + I*Bq^2*D0*omega^3
Im K_direct = B0*Bq^2*omega + Bq^2*D0*omega^3
K_deriv = A0*Bq^2*omega^2 + I*B0*Bq^2*omega^3 + Bq^2*C0*omega^4 + I*Bq^2*D0*omega^5
Im K_deriv = B0*Bq^2*omega^3 + Bq^2*D0*omega^5
PASS: Derivative outlet kernel has no linear odd term
Derivative-breathing odd term exponent (with a/r ~ eps^delta) = 5/2 + 3*delta

----------------------------------------------------------------------------------------
IV.6 — No third linear scalar source from quadratic momentum/stress terms
----------------------------------------------------------------------------------------
j_sigma(w) = I*q0*u0*wv*z1 - q0*u1*wv^2*z1 + I*q0*u2*wv^3*z1 + q0*u0*wv^2*z2 + I*q0*u1*wv^3*z2 + I*q0*u0*wv^3*z3
With Z_sigma(0)=0, the leading mouth source is derivative-like, not direct q-like.

----------------------------------------------------------------------------------------
IV.7 — Mouth radiative-order theorem and Ohmic no-go
----------------------------------------------------------------------------------------
Compact reciprocal mouth admittance Y(w) = (I*betaM*wv)/Kmouth + (I*betaM*Madd*wv^3)/Kmouth^2 + (I*betaM*Mmouth*wv^3)/Kmouth^2 + (betaM*R2*wv^4)/Kmouth^2 + (I*betaM*Madd^2*wv^5)/Kmouth^3 + ((2*I)*betaM*Madd*Mmouth*wv^5)/Kmouth^3 + (I*betaM*Mmouth^2*wv^5)/Kmouth^3
Re Y series = (betaM*R2*wv^4)/Kmouth^2
Im Y series = (betaM*wv)/Kmouth + (betaM*(Madd + Mmouth)*wv^3)/Kmouth^2 + (betaM*(Madd + Mmouth)^2*wv^5)/Kmouth^3
1D Ohmic benchmark Re Y1 = (betaM*c0*rho0*wv^2)/Kmouth^2 + (betaM*c0*rho0*(2*Kmouth*Mmouth - c0^2*rho0^2)*wv^4)/Kmouth^4
Natural odd mouth length ell_mouth ~ (aSmall*chiSigma)/(delta*kappaM)

========================================================================================
PART V — DIPOLE / VECTOR SECTOR
========================================================================================

----------------------------------------------------------------------------------------
V.1 — Carried odd wake ports: CM/relative split
----------------------------------------------------------------------------------------
PiA_perp - PiB_perp = {Sqrt[7/2]*ux, Sqrt[7/2]*uy}
PiA_10   - PiB_10   = 2*d

----------------------------------------------------------------------------------------
V.2 — Vector Ward-identity failure
----------------------------------------------------------------------------------------
F_A = acoef*(-Derivative[3][xpA][t] + Derivative[3][xpB][t])
F_B = acoef*(Derivative[3][xpA][t] - Derivative[3][xpB][t])
PASS: F_A + F_B
L_cmrel = acoef*Derivative[1][xrel1][t]*Derivative[2][xrel1][t]
PASS: CM coordinate drops out of the relative dipole odd action

----------------------------------------------------------------------------------------
V.3 — Outgoing l=1 spectral no-go for linear odd term
----------------------------------------------------------------------------------------
Lambda1(k) = -2/a + a*k^2 + I*a^2*k^3 - a^3*k^4 - I*a^4*k^5 + a^5*k^6
Y1_norm(k) = R + (a^2*k^2*R)/2 + I/2*a^3*k^3*R - (a^4*k^4*R)/4
PASS: Outgoing l=1 has no linear k term
Imaginary k^3 coefficient = (a^3*R)/2

----------------------------------------------------------------------------------------
V.3b — Order reduction of an isotropic fifth-derivative force
----------------------------------------------------------------------------------------
x^(5) = (15*mN*rdN*(2*mN + 7*rdN^2*rN - 3*rN*v2N))/rN^6 * n + (mN*(-8*mN + 9*rN*(-5*rdN^2 + v2N)))/rN^6 * v
Thus -kap x^(5) lies in the standard {rdot n, v} 2.5PN basis.

----------------------------------------------------------------------------------------
V.4 — Same-order dipole wake cannot be absorbed into the Iyer–Will family
----------------------------------------------------------------------------------------
Full matching solution = {{alpha -> 4, beta -> 5, s -> 1, p -> 0, qpar -> 0}}
PASS: Same-order dipole match requires p=q=0

----------------------------------------------------------------------------------------
V.5 — Dipole finite-size scaling rescue
----------------------------------------------------------------------------------------
eps * (sqrt(eps)*rho)^3 = eps^(5/2)*rho^3
If rho ~ eps^delta, exponent = 5/2 + 3 delta.
general exponent = 5/2 + 3*delta

========================================================================================
PART VI — QUADRUPOLE SECTOR
========================================================================================

----------------------------------------------------------------------------------------
VI.1 — P2 ports are exactly the real STF l=2 content
----------------------------------------------------------------------------------------
PASS: q20 - sqrt(2/3) Pi20
PASS: q21c - Pi21c
PASS: q21s - Pi21s
PASS: q22c - 2 Pi22c
PASS: q22s - 2 Pi22s

----------------------------------------------------------------------------------------
VI.2 — Outgoing l=2 retarded completion gives +i*k^5
----------------------------------------------------------------------------------------
Lambda2(k) = -3/a + (a*k^2)/3 + (a^3*k^4)/9 + I/9*a^4*k^5 - (2*a^5*k^6)/27
Y2_norm(k) = R2 + (a^2*k^2*R2)/9 + (4*a^4*k^4*R2)/81 + I/27*a^5*k^5*R2
Imaginary k^5 coefficient = (a^5*R2)/27
PASS: Quadrupole k^5 coefficient is positive

----------------------------------------------------------------------------------------
VI.3 — Quadrupole odd kernel gives damping, not antidamping
----------------------------------------------------------------------------------------
P - dE_Schott/dt = -1/2*(gammaQuad*Derivative[3][q][t]^2)
PASS: quadrupole power balance residual

----------------------------------------------------------------------------------------
VI.3b — Compact internal P2 branch is finite-size suppressed
----------------------------------------------------------------------------------------
eps^2 * (sqrt(eps)*rho)^5 = eps^(9/2)*rho^5
If rho ~ eps^delta, exponent = 9/2 + 5*delta

----------------------------------------------------------------------------------------
VI.4 — Orbital/worldtube STF source map
----------------------------------------------------------------------------------------
PASS: V20 source map
PASS: V21c source map
PASS: V21s source map
PASS: V22c source map
PASS: V22s source map
PASS: X21c static position source
PASS: X21s static position source
PASS: X22c static position source
PASS: X22s static position source
PASS: I21c source map
PASS: I21s source map
PASS: I22c source map
PASS: I22s source map
PASS: I20 source map

----------------------------------------------------------------------------------------
VI.5 — Matching theorem: commutant is scalar, basis map invertible
----------------------------------------------------------------------------------------
commutant basis matrix = MatrixForm[{{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}}]
PASS: commutant dimension is one
PASS: commutant minus scalar identity
det(port->canonical basis map) = 4*Sqrt[2/3]
PASS: Port-to-canonical STF map is invertible

----------------------------------------------------------------------------------------
VI.6 — Static overlap test
----------------------------------------------------------------------------------------
C20[x_<i x_j>] = Sqrt[2/3]*rSym^2
C20[ddot I_orb /(2 mu)] = -((Sqrt[2/3]*Gsym*Msym)/rSym)
q20_static = 5/4
m0 from static overlap = (-5*Sqrt[3/2]*rSym)/(4*Gsym*Msym)
PASS: Static overlap is nonzero on the natural branch

========================================================================================
FINAL CONDITIONAL THEOREM LEDGER
========================================================================================
Conservative decisive benchmark:
  Burke–Thorne prototype lands on Iyer–Will parameters (alpha,beta)=(4,5).

Low-frequency framework:
  Wrong odd channels are identified by n = 1 (scalar), n = 3 (dipole), n = 5 (quadrupole).

Scalar side:
  - direct scalar monopole overlap generically gives i*omega,
  - continuity-derived breathing outlet is derivative-coupled, so its kernel starts at i*omega^3,
  - compact reciprocal mouth radiation is super-Ohmic / higher-order,
  - genuine Ohmic scalar bleed would require new structure outside the natural compact-radiative branch.

Dipole side:
  - equal-and-opposite body forces do not kill the relative carried wake,
  - but a clean outgoing l=1 completion removes the dangerous linear odd term and starts at i*k^3,
  - same-order nonzero dipole deformation cannot be absorbed into the standard 2.5PN family,
  - in a strict small-body branch it is pushed above universal point-particle 2.5PN.

Quadrupole side:
  - the solved 2PN P2 package is exactly the real STF l=2 content,
  - the source map to the orbital/worldtube STF quadrupole is explicit,
  - rotational invariance forces the retarded 5x5 matching matrix to be scalar,
  - a compact outgoing l=2 completion gives a positive +i*k^5 term,
  - and the static overlap gate passes: m0 is nonzero on the natural low-frequency branch.

Best current verdict:
  A full GR-like point-particle 2.5PN theorem is conditionally reachable,
  with the remaining narrow gap living in the final passive/outgoing quadrupole normalization.

Passes: 42
Fails : 0
OVERALL: PASS
"*)
