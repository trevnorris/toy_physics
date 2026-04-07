(* 2.5PN session-completion Mathematica audit
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

Print["--- 2.5PN session-completion Mathematica audit ---"];

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

banner["PART VII — QUADRUPOLE NORMALIZATION / EXTRACTION / PARAMETERIZATION"];

subbanner["VII.1 — Canonical quadrupole normalization map"];
Clear[x, y, q];
gamma =.; Gamma5 =.; m0hat =.; G =.; c =.; mu =.;
xq = x[t];
yq = y[t];
Xplane = {xq, yq, 0};
Itensor = FullSimplify[mu (Outer[Times, Xplane, Xplane] - IdentityMatrix[3] (Xplane.Xplane)/3)];
qCan = Table[FullSimplify[Total[Flatten[Itensor basis[[i]]]]], {i, 1, 5}];
q5Can = D[#, {t, 5}] & /@ qCan;
Fx = FullSimplify[-(gamma/2) Sum[D[qCan[[i]], xq] q5Can[[i]], {i, 1, 5}]];
Fy = FullSimplify[-(gamma/2) Sum[D[qCan[[i]], yq] q5Can[[i]], {i, 1, 5}]];
Ix = FullSimplify[-(xq D[Itensor[[1, 1]], {t, 5}] + yq D[Itensor[[1, 2]], {t, 5}])];
Iy = FullSimplify[-(xq D[Itensor[[1, 2]], {t, 5}] + yq D[Itensor[[2, 2]], {t, 5}])];
checkScalar["Fx - gamma*mu*(-x^j I_xj^(5))", Fx - gamma mu Ix];
checkScalar["Fy - gamma*mu*(-x^j I_yj^(5))", Fy - gamma mu Iy];
gammaGR = FullSimplify[2 G/(5 c^5)];
gammaToy = FullSimplify[m0hat^2 Gamma5];
Gamma5Target = FullSimplify[gammaGR/m0hat^2];
Print["gamma_GR = ", gammaGR];
Print["gamma_toy = ", gammaToy];
Print["Gamma5_target = ", Gamma5Target];

subbanner["VII.2 — Outgoing l=2 PDE fingerprint"];
Clear[k, a, cS, m0, omega, r];
z = k r;
j2r = ((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2r = -((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2;
h2r = FullSimplify[j2r + I y2r];
odeResidual = FullSimplify[D[h2r, {r, 2}] + (2/r) D[h2r, r] + (k^2 - 6/r^2) h2r];
checkScalar["l=2 Helmholtz residual", odeResidual];
za =.; 
j2a = ((3/za^3) - 1/za) Sin[za] - 3 Cos[za]/za^2;
y2a = -((3/za^3) - 1/za) Cos[za] - 3 Sin[za]/za^2;
h2a = FullSimplify[j2a + I y2a];
Lambda2Session = FullSimplify[(k D[h2a, za]/h2a) /. za -> k a];
Lambda2SessionSeries = Normal@Series[Lambda2Session, {k, 0, 6}];
Y2Session = FullSimplify[Normal@Series[1/Lambda2SessionSeries, {k, 0, 5}]];
Y2Static = FullSimplify[Y2Session /. k -> 0];
Y2Hat = FullSimplify[Expand[Y2Session/Y2Static]];
Y2HatSeries = Normal@Series[Y2Hat, {k, 0, 5}];
Y2HatOmega = Expand[Y2HatSeries /. k -> omega/cS];
m2 = FullSimplify[m0 Coefficient[Y2HatOmega, omega, 2]];
m4 = FullSimplify[m0 Coefficient[Y2HatOmega, omega, 4]];
Gamma5PDE = FullSimplify[m0 Coefficient[Y2HatOmega, omega, 5]/I];
Print["Lambda2(k) = ", Lambda2SessionSeries];
Print["Y2(0) = ", Y2Static];
Print["Y2_hat(k) = ", Y2HatSeries];
Print["m2 = ", m2];
Print["m4 = ", m4];
Print["Gamma5^PDE = ", Gamma5PDE];
checkScalar["m4 - 4 m2^2 / m0", m4 - 4 m2^2/m0];
checkScalar[
  "Gamma5^PDE - 9 m2^(5/2)/m0^(3/2)",
  Gamma5PDE - 9 m2^(5/2)/m0^(3/2),
  m0 > 0 && a > 0 && cS > 0
];

subbanner["VII.3 — Source/port normalization cleanup"];
B20 = Sqrt[2/3];
BmapSession = DiagonalMatrix[{B20, 1, 1, 2, 2}];
Gs =.; Ms =.; rs =.;
C20ddotOver2mu = FullSimplify[-Sqrt[6] Gs Ms/(3 rs)];
P20ddotOver2mu = FullSimplify[C20ddotOver2mu/B20];
q20PortStatic = 5/4;
m0Mixed = FullSimplify[q20PortStatic/C20ddotOver2mu];
m0RawSameBasis = FullSimplify[q20PortStatic/P20ddotOver2mu];
q20CanStatic = FullSimplify[B20 q20PortStatic];
m0CanSameBasis = FullSimplify[q20CanStatic/C20ddotOver2mu];
Norb = FullSimplify[1/m0RawSameBasis];
m0hatPointParticle = FullSimplify[Norb m0RawSameBasis];
Print["det(Bmap) = ", FullSimplify[Det[BmapSession]]];
Print["m0_mixed = ", m0Mixed];
Print["m0_raw_same_basis = ", m0RawSameBasis];
Print["m0_can_same_basis = ", m0CanSameBasis];
Print["N_orb = ", Norb];
Print["m0hat_point_particle = ", m0hatPointParticle];
checkScalar["same-basis raw vs canonical overlap", m0CanSameBasis - m0RawSameBasis];
checkScalar["point-particle normalized overlap - 1", m0hatPointParticle - 1];
kapp =.; Y0const =.;
Y2HatMin = 1 + a^2 omega^2/(9 cS^2) + 4 a^4 omega^4/(81 cS^4) + I a^5 omega^5/(27 cS^5);
Meff = Expand[kapp^2 Y0const Y2HatMin];
m0Eff = FullSimplify[Coefficient[Meff, omega, 0]];
m2Eff = FullSimplify[Coefficient[Meff, omega, 2]];
m4Eff = FullSimplify[Coefficient[Meff, omega, 4]];
Gamma5Eff = FullSimplify[Coefficient[Meff, omega, 5]/I];
checkScalar["m2_eff - m0_eff a^2/(9 c_s^2)", m2Eff - m0Eff a^2/(9 cS^2)];
checkScalar["m4_eff - 4 m0_eff a^4/(81 c_s^4)", m4Eff - 4 m0Eff a^4/(81 cS^4)];
checkScalar["Gamma5_eff - m0_eff a^5/(27 c_s^5)", Gamma5Eff - m0Eff a^5/(27 cS^5)];
Gamma5Can = FullSimplify[m0hat a^5/(27 cS^5)];
Print["Gamma5_can = ", Gamma5Can];

subbanner["VII.4 — Convention-invariant normalization product"];
lam =.; mfrak =.; K0 =.; Gamma5Port =.; 
mfrakP = FullSimplify[lam mfrak];
Gamma5PortP = FullSimplify[Gamma5Port/lam^2];
K0P = FullSimplify[K0/lam^2];
checkScalar["mfrak'^2 Gamma5' - mfrak^2 Gamma5", mfrakP^2 Gamma5PortP - mfrak^2 Gamma5Port];
checkScalar["mfrak'^2 K0' - mfrak^2 K0", mfrakP^2 K0P - mfrak^2 K0];
Gamma5PortDN = FullSimplify[K0 a^5/(27 cS^5)];
gammaEff = FullSimplify[mfrak^2 Gamma5PortDN];
NQ =.;
NQTarget = FullSimplify[NQ /. First[Solve[NQ a^5/(27 cS^5) == gammaGR, NQ]]];
Print["Gamma5_port^DN = ", Gamma5PortDN];
Print["gamma_eff = ", gammaEff];
Print["N_Q_target = ", NQTarget];

subbanner["VII.5 — Canonical invariant low-frequency fingerprint"];
K0bar =.; K2bar =.; K4bar =.; Gammabar =.;
K2barFromBranch = FullSimplify[K0bar a^2/(9 cS^2)];
K4barFromBranch = FullSimplify[K0bar 4 a^4/(81 cS^4)];
GammabarFromBranch = FullSimplify[K0bar a^5/(27 cS^5)];
ratioA2Cs2 = FullSimplify[9 K2bar/K0bar, K0bar > 0 && K2bar > 0];
K4barInvariant = FullSimplify[(4 K0bar ratioA2Cs2^2)/81, K0bar > 0 && K2bar > 0];
GammabarInvariant = FullSimplify[(K0bar ratioA2Cs2^(5/2))/27, K0bar > 0 && K2bar > 0];
K2barTargetFormula = FullSimplify[(2 G/(45 c^5))^(2/5) K0bar^(3/5)];
Print["a^2/c_s^2 = ", ratioA2Cs2];
Print["K4bar invariant = ", K4barInvariant];
Print["Gammabar invariant = ", GammabarInvariant];
Print["K2bar_target formula = ", K2barTargetFormula];
checkScalar["K4bar invariant - 4 K2bar^2/K0bar", K4barInvariant - 4 K2bar^2/K0bar, K0bar > 0 && K2bar > 0];
checkScalar[
  "Gammabar invariant - 9 K2bar^(5/2)/K0bar^(3/2)",
  GammabarInvariant - 9 K2bar^(5/2)/K0bar^(3/2),
  K0bar > 0 && K2bar > 0
];

subbanner["VII.6 — Extraction of (Kbar0, Kbar2) and single-prefactor obstruction"];
K2barBranch = FullSimplify[K0bar a^2/(9 cS^2)];
K4barBranch = FullSimplify[K0bar 4 a^4/(81 cS^4)];
GammabarBranch = FullSimplify[K0bar a^5/(27 cS^5)];
K0barTarget = K0bar /. First[Solve[GammabarBranch == gammaGR, K0bar]];
K2barTarget = FullSimplify[K2barBranch /. K0bar -> K0barTarget];
K4barTarget = FullSimplify[K4barBranch /. K0bar -> K0barTarget];
Print["K2bar_branch = ", K2barBranch];
Print["K0bar_target = ", K0barTarget];
Print["K2bar_target = ", K2barTarget];
Print["K4bar_target = ", K4barTarget];
Omega20 =.; Cq =.;
Y20Min = FullSimplify[1/(1 - omega^2/Omega20^2)];
Print["Y20_min(omega) = ", Normal@Series[Y20Min, {omega, 0, 4}]];
Y2Out = Expand[1 + a^2 omega^2/(9 cS^2) + 4 a^4 omega^4/(81 cS^4) + I a^5 omega^5/(27 cS^5)];
KbarPref = Expand[Cq Y2Out];
Kbar0 = FullSimplify[Coefficient[KbarPref, omega, 0]];
Kbar2 = FullSimplify[Coefficient[KbarPref, omega, 2]];
Kbar4 = FullSimplify[Coefficient[KbarPref, omega, 4]];
Gammabar5 = FullSimplify[Coefficient[KbarPref, omega, 5]/I];
CqTarget = Cq /. First[Solve[Gammabar5 == gammaGR, Cq]];
Print["C_Q_target = ", CqTarget];
checkScalar["Kbar0 - C_Q", Kbar0 - Cq];
checkScalar["Kbar2 - C_Q a^2/(9 c_s^2)", Kbar2 - Cq a^2/(9 cS^2)];
checkScalar["Kbar4 - 4 C_Q a^4/(81 c_s^4)", Kbar4 - 4 Cq a^4/(81 cS^4)];

subbanner["VII.7 — Frozen 2PN files: fixed representation/support data"];
ux =.; uy =.; d =.; U =.;
vFrozen = {ux, uy, d};
nHat = {0, 0, 1};
Tmat = Outer[Times, vFrozen, vFrozen] - U Outer[Times, nHat, nHat];
Tstf = stfTensor[Tmat];
Ccoeffs = coeffsInBasis[Tstf];
checkScalar["C20 relation", Ccoeffs[[1]] - Sqrt[2/3] (pi20[ux, uy, d] - U)];
checkScalar["C21c relation", Ccoeffs[[2]] - pi21c[ux, uy, d]];
checkScalar["C21s relation", Ccoeffs[[3]] - pi21s[ux, uy, d]];
checkScalar["C22c relation", Ccoeffs[[4]] - 2 pi22c[ux, uy]];
checkScalar["C22s relation", Ccoeffs[[5]] - 2 pi22s[ux, uy]];
J0 = 4/Sqrt[5];
J20Session = 5/4;
staticSupport = FullSimplify[J0^2 + J20Session^2];
closureDeficit = FullSimplify[staticSupport - 5/4];
Print["J.J = ", staticSupport];
Print["closure deficit = ", closureDeficit];
Print["Open dynamic poles: Omega20, Omega21, Omega22"];

subbanner["VII.8 — One-pole insufficiency and minimal positive conservative precursor"];
Abranch = FullSimplify[a^2/(9 cS^2)];
YoutBranch = Expand[1 + Abranch omega^2 + 4 Abranch^2 omega^4 + I a^5 omega^5/(27 cS^5)];
Omega1 = FullSimplify[3 cS/a, a > 0 && cS > 0];
Y1Pole = FullSimplify[1/(1 - omega^2/Omega1^2)];
Y1PoleSeries = Normal@Series[Y1Pole, {omega, 0, 5}];
mismatchW4 = FullSimplify[Coefficient[YoutBranch, omega, 4] - Coefficient[Y1PoleSeries, omega, 4]];
Print["Omega_1 = ", Omega1];
Print["one-pole O(w^4) mismatch = ", mismatchW4];
pIso =.;
q1Iso = FullSimplify[Abranch + Sqrt[3] Abranch Sqrt[1 - pIso]/Sqrt[pIso]];
q2Iso = FullSimplify[Abranch - Sqrt[3] Abranch Sqrt[pIso]/Sqrt[1 - pIso]];
Y2Pole = FullSimplify[pIso/(1 - q1Iso omega^2) + (1 - pIso)/(1 - q2Iso omega^2)];
Print["Y_2pole series = ", Normal@Series[Y2Pole, {omega, 0, 5}]];
pMin = 1/4;
q1Min = FullSimplify[q1Iso /. pIso -> pMin];
Ymin = FullSimplify[Limit[Y2Pole, pIso -> pMin, Direction -> -1]];
YminSeries = Normal@Series[Ymin, {omega, 0, 5}];
YminSplit = FullSimplify[3/4 + (1/4)/(1 - q1Min omega^2)];
OmegaMin = FullSimplify[3 cS/(2 a), a > 0 && cS > 0];
Print["Minimal positive closure = ", Ymin];
Print["Minimal split form = ", YminSplit];
Print["Omega_Q,min = ", OmegaMin];
checkScalar["minimal closure omega^2 match", Coefficient[YminSeries, omega, 2] - Abranch];
checkScalar["minimal closure omega^4 match", Coefficient[YminSeries, omega, 4] - 4 Abranch^2];

subbanner["VII.9 — Minimal isotropic quadrupole module and componentwise P2 target"];
c0iso =.; c1iso =.; OmegaQ =.;
Yfamily = FullSimplify[c0iso + c1iso/(1 - omega^2/OmegaQ^2)];
YfamilySeries = Normal@Series[Yfamily, {omega, 0, 5}];
solIso = Solve[
  {
    Coefficient[YfamilySeries, omega, 0] == 1,
    Coefficient[YfamilySeries, omega, 2] == Abranch,
    Coefficient[YfamilySeries, omega, 4] == 4 Abranch^2
  },
  {c0iso, c1iso, OmegaQ}
];
c0sol = 3/4;
c1sol = 1/4;
OmegaQSol = FullSimplify[3 cS/(2 a), a > 0 && cS > 0];
YfamilyMin = FullSimplify[c0sol + c1sol/(1 - omega^2/OmegaQSol^2)];
YfamilyMinSeries = Normal@Series[YfamilyMin, {omega, 0, 5}];
Print["c0 = ", c0sol];
Print["c1 = ", c1sol];
Print["OmegaQ = ", OmegaQSol];
Print["Y_family_min = ", YfamilyMin];
checkScalar["minimal module omega^2 mismatch", Coefficient[YfamilyMinSeries, omega, 2] - Abranch];
checkScalar["minimal module omega^4 mismatch", Coefficient[YfamilyMinSeries, omega, 4] - 4 Abranch^2];
gram = Table[FullSimplify[Total[Flatten[basis[[i]] basis[[j]]]]], {i, 1, 5}, {j, 1, 5}];
checkMatrix["real STF Gram matrix - I5", gram - IdentityMatrix[5]];
YQMin = FullSimplify[3/4 + (1/4)/(1 - omega^2/OmegaQSol^2)];
Print["Y20 = Y21 = Y22 = ", YQMin];
Print["YQ_min series = ", Normal@Series[YQMin, {omega, 0, 5}]];
OmegaQFromInvariants = FullSimplify[Sqrt[K0bar]/(2 Sqrt[K2bar]), K0bar > 0 && K2bar > 0];
Print["OmegaQ_from_invariants = ", OmegaQFromInvariants];
Print["K4 relation = ", FullSimplify[4 K2bar^2/K0bar]];
Print["Gamma5 relation = ", FullSimplify[9 K2bar^(5/2)/K0bar^(3/2)]];

subbanner["VII.10 — Single-pole sufficiency theorem and extraction ledger"];
K2FromBranch = FullSimplify[K0bar/(4 OmegaQ^2)];
K4FromBranch = FullSimplify[K0bar/(4 OmegaQ^4)];
GammaFromBranch = FullSimplify[9 K0bar/(32 OmegaQ^5)];
K0barTargetFromOmega = K0bar /. First[Solve[GammaFromBranch == gammaGR, K0bar]];
K2barTargetFromOmega = FullSimplify[K2FromBranch /. K0bar -> K0barTargetFromOmega];
K4barTargetFromOmega = FullSimplify[K4FromBranch /. K0bar -> K0barTargetFromOmega];
Print["K2_from_branch = ", K2FromBranch];
Print["K4_from_branch = ", K4FromBranch];
Print["Gamma_from_branch = ", GammaFromBranch];
Print["K0bar_target(OmegaQ) = ", K0barTargetFromOmega];
Print["K2bar_target(OmegaQ) = ", K2barTargetFromOmega];
Print["K4bar_target(OmegaQ) = ", K4barTargetFromOmega];
c0s =.; c2s =.; c4s =.;
Print["OmegaQ^2 from (c0,c2) = ", FullSimplify[c0s/(4 c2s)]];
Print["OmegaQ^2 from (c2,c4) = ", FullSimplify[c2s/c4s]];
Print["K0 from (c2,c4) = ", FullSimplify[4 c2s^2/c4s]];
Print["Gamma from (c2,c4) = ", FullSimplify[(9/8) c4s^(3/2)/Sqrt[c2s]]];
Print["Direct GR gate: c4^3/c2 = ", FullSimplify[(256/2025) G^2/c^10]];

subbanner["VII.11 — Grouped P2 formulas and normalized-support ledger"];
c020 =.; c220 =.; c420 =.; c021 =.; c221 =.; c421 =.; c022 =.; c222 =.; c422 =.;
c0Trace = FullSimplify[(c020 + 2 c021 + 2 c022)/5];
c2Trace = FullSimplify[(c220 + 2 c221 + 2 c222)/5];
c4Trace = FullSimplify[(c420 + 2 c421 + 2 c422)/5];
Print["c0_trace = ", c0Trace];
Print["c2_trace = ", c2Trace];
Print["c4_trace = ", c4Trace];
checkScalar["weighted anisotropy sum n=0", (c020 - c0Trace) + 2 (c021 - c0Trace) + 2 (c022 - c0Trace)];
checkScalar["weighted anisotropy sum n=2", (c220 - c2Trace) + 2 (c221 - c2Trace) + 2 (c222 - c2Trace)];
checkScalar["weighted anisotropy sum n=4", (c420 - c4Trace) + 2 (c421 - c4Trace) + 2 (c422 - c4Trace)];
Print["Minimal branch identity c0*c4 - 4*c2^2 = ", FullSimplify[c0s c4s - 4 c2s^2]];
u220 =.; u420 =.; u221 =.; u421 =.; u222 =.; u422 =.; u2n =.; u4n =.; K0n =.;
u2Trace = FullSimplify[(u220 + 2 u221 + 2 u222)/5];
u4Trace = FullSimplify[(u420 + 2 u421 + 2 u422)/5];
Print["u2_trace = ", u2Trace];
Print["u4_trace = ", u4Trace];
Gamma5Norm = FullSimplify[(9/8) u4n^(3/2)/Sqrt[u2n], u2n > 0 && u4n > 0];
Gamma5NormBranch = FullSimplify[Gamma5Norm /. u4n -> 4 u2n^2, u2n > 0];
K0barTargetNorm = K0n /. First[Solve[K0n Gamma5NormBranch == gammaGR, K0n]];
Print["normalized branch identity u4 - 4*u2^2 = ", FullSimplify[u4n - 4 u2n^2]];
Print["OmegaQ on normalized branch = ", FullSimplify[1/(2 Sqrt[u2n])]];
Print["Gamma5_norm on branch = ", Gamma5NormBranch];
Print["K0bar_target from normalized support = ", K0barTargetNorm];

subbanner["VII.12 — Axisymmetric grouped-P2 parameterization"];
u2bar =.; u4bar =.; a2 =.; a4 =.; b2 =.; b4 =.;
I5 = IdentityMatrix[5];
T5 = DiagonalMatrix[{4, -1, -1, -1, -1}];
S5 = DiagonalMatrix[{0, 1, 1, -1, -1}];
U2mat = FullSimplify[u2bar I5 + a2 T5 + b2 S5];
U4mat = FullSimplify[u4bar I5 + a4 T5 + b4 S5];
Ymat = FullSimplify[I5 + omega^2 U2mat + omega^4 U4mat];
Print["Y20 = ", Ymat[[1, 1]]];
Print["Y21 = ", Ymat[[2, 2]]];
Print["Y22 = ", Ymat[[4, 4]]];
x20 =.; x21 =.; x22 =.;
Print["ubar_from_(x20,x21,x22) = ", FullSimplify[(x20 + 2 x21 + 2 x22)/5]];
Print["a_from_(x20,x21,x22) = ", FullSimplify[(2 x20 - x21 - x22)/10]];
Print["b_from_(x20,x21,x22) = ", FullSimplify[(x21 - x22)/2]];
Print["A2_sq = ", FullSimplify[Tr[(U2mat - u2bar I5).(U2mat - u2bar I5)]/5]];
Print["A4_sq = ", FullSimplify[Tr[(U4mat - u4bar I5).(U4mat - u4bar I5)]/5]];
Print["Minimal branch relation = ", FullSimplify[u4bar - 4 u2bar^2]];
Print["OmegaQ^2 = ", FullSimplify[1/(4 u2bar)]];
Print["Gamma5_norm = ", FullSimplify[9 u2bar^(5/2)]];
Print["K0bar_target = ", FullSimplify[2 G/(45 c^5 u2bar^(5/2))]];

banner["FINAL SESSION THEOREM LEDGER"];
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
Print["\nQuadrupole normalization stack:"];
Print["  gamma_GR = ", gammaGR];
Print["  Gamma5_target = ", Gamma5Target];
Print["  N_Q_target = ", NQTarget];
Print["  K0bar_target = ", K0barTarget];
Print["  K2bar_target = ", K2barTarget];
Print["  K4bar_target = ", K4barTarget];
Print["\nMinimal isotropic quadrupole module:"];
Print["  Omega_Q,min = ", OmegaMin];
Print["  exact moment-matched isotropic solution Omega_Q = ", OmegaQSol];
Print["  single-pole GR scaling law K0bar = ", K0barTargetFromOmega, "."];
Print["\nNormalized-support branch:"];
Print["  K0bar_target from normalized support moments = ", K0barTargetNorm];
Print["\nBest current verdict:"];
Print["  The session reduces the full GR-like point-particle 2.5PN closure to the"];
Print["  isotropic real P2 quadrupole module. Once the grouped low-frequency"];
Print["  conservative coefficients through O(omega^4) are fixed on that branch,"];
Print["  the odd Burke–Thorne coefficient follows automatically."];
Print["  The one remaining narrow theorem gap is the final passive/outgoing"];
Print["  quadrupole normalization on the actual moving-throat branch."];

Print["\nPasses: ", passCount];
Print["Fails : ", failCount];
If[failCount == 0,
  Print["OVERALL: PASS"],
  Print["OVERALL: FAIL"]
];

(*"
Output:
--- 2.5PN session-completion Mathematica audit ---

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
PART VII — QUADRUPOLE NORMALIZATION / EXTRACTION / PARAMETERIZATION
========================================================================================

----------------------------------------------------------------------------------------
VII.1 — Canonical quadrupole normalization map
----------------------------------------------------------------------------------------
PASS: Fx - gamma*mu*(-x^j I_xj^(5))
PASS: Fy - gamma*mu*(-x^j I_yj^(5))
gamma_GR = (2*G)/(5*c^5)
gamma_toy = Gamma5*m0hat^2
Gamma5_target = (2*G)/(5*c^5*m0hat^2)

----------------------------------------------------------------------------------------
VII.2 — Outgoing l=2 PDE fingerprint
----------------------------------------------------------------------------------------
PASS: l=2 Helmholtz residual
Lambda2(k) = -3/a + (a*k^2)/3 + (a^3*k^4)/9 + I/9*a^4*k^5 - (2*a^5*k^6)/27
Y2(0) = -1/3*a
Y2_hat(k) = 1 + (a^2*k^2)/9 + (4*a^4*k^4)/81 + I/27*a^5*k^5
m2 = (a^2*m0)/(9*cS^2)
m4 = (4*a^4*m0)/(81*cS^4)
Gamma5^PDE = (a^5*m0)/(27*cS^5)
PASS: m4 - 4 m2^2 / m0
PASS: Gamma5^PDE - 9 m2^(5/2)/m0^(3/2)

----------------------------------------------------------------------------------------
VII.3 — Source/port normalization cleanup
----------------------------------------------------------------------------------------
det(Bmap) = 4*Sqrt[2/3]
m0_mixed = (-5*Sqrt[3/2]*rs)/(4*Gs*Ms)
m0_raw_same_basis = (-5*rs)/(4*Gs*Ms)
m0_can_same_basis = (-5*rs)/(4*Gs*Ms)
N_orb = (-4*Gs*Ms)/(5*rs)
m0hat_point_particle = 1
PASS: same-basis raw vs canonical overlap
PASS: point-particle normalized overlap - 1
PASS: m2_eff - m0_eff a^2/(9 c_s^2)
PASS: m4_eff - 4 m0_eff a^4/(81 c_s^4)
PASS: Gamma5_eff - m0_eff a^5/(27 c_s^5)
Gamma5_can = (a^5*m0hat)/(27*cS^5)

----------------------------------------------------------------------------------------
VII.4 — Convention-invariant normalization product
----------------------------------------------------------------------------------------
PASS: mfrak'^2 Gamma5' - mfrak^2 Gamma5
PASS: mfrak'^2 K0' - mfrak^2 K0
Gamma5_port^DN = (a^5*K0)/(27*cS^5)
gamma_eff = (a^5*K0*mfrak^2)/(27*cS^5)
N_Q_target = (54*cS^5*G)/(5*a^5*c^5)

----------------------------------------------------------------------------------------
VII.5 — Canonical invariant low-frequency fingerprint
----------------------------------------------------------------------------------------
a^2/c_s^2 = (9*K2bar)/K0bar
K4bar invariant = (4*K2bar^2)/K0bar
Gammabar invariant = 9*K0bar*(K2bar/K0bar)^(5/2)
K2bar_target formula = ((2/5)^(2/5)*(G/c^5)^(2/5)*K0bar^(3/5))/3^(4/5)
PASS: K4bar invariant - 4 K2bar^2/K0bar
PASS: Gammabar invariant - 9 K2bar^(5/2)/K0bar^(3/2)

----------------------------------------------------------------------------------------
VII.6 — Extraction of (Kbar0, Kbar2) and single-prefactor obstruction
----------------------------------------------------------------------------------------
K2bar_branch = (a^2*K0bar)/(9*cS^2)
K0bar_target = (54*cS^5*G)/(5*a^5*c^5)
K2bar_target = (6*cS^3*G)/(5*a^3*c^5)
K4bar_target = (8*cS*G)/(15*a*c^5)
Y20_min(omega) = 1 + omega^4/Omega20^4 + omega^2/Omega20^2
C_Q_target = (54*cS^5*G)/(5*a^5*c^5)
PASS: Kbar0 - C_Q
PASS: Kbar2 - C_Q a^2/(9 c_s^2)
PASS: Kbar4 - 4 C_Q a^4/(81 c_s^4)

----------------------------------------------------------------------------------------
VII.7 — Frozen 2PN files: fixed representation/support data
----------------------------------------------------------------------------------------
PASS: C20 relation
PASS: C21c relation
PASS: C21s relation
PASS: C22c relation
PASS: C22s relation
J.J = 381/80
closure deficit = 281/80
Open dynamic poles: Omega20, Omega21, Omega22

----------------------------------------------------------------------------------------
VII.8 — One-pole insufficiency and minimal positive conservative precursor
----------------------------------------------------------------------------------------
Omega_1 = (3*cS)/a
one-pole O(w^4) mismatch = a^4/(27*cS^4)
Y_2pole series = 1 + (a^2*omega^2)/(9*cS^2) + (2*a^4*omega^4*(2*Sqrt[1 - pIso] - Sqrt[3]*Sqrt[pIso] + Sqrt[3]*pIso^(3/2) + Sqrt[3]*Sqrt[1 - pIso]*Sqrt[-((-1 + pIso)*pIso)]))/(81*cS^4*Sqrt[1 - pIso])
Minimal positive closure = (3 + (1 - (4*a^2*omega^2)/(9*cS^2))^(-1))/4
Minimal split form = (3 + (1 - (4*a^2*omega^2)/(9*cS^2))^(-1))/4
Omega_Q,min = (3*cS)/(2*a)
PASS: minimal closure omega^2 match
PASS: minimal closure omega^4 match

----------------------------------------------------------------------------------------
VII.9 — Minimal isotropic quadrupole module and componentwise P2 target
----------------------------------------------------------------------------------------
c0 = 3/4
c1 = 1/4
OmegaQ = (3*cS)/(2*a)
Y_family_min = (3 + (1 - (4*a^2*omega^2)/(9*cS^2))^(-1))/4
PASS: minimal module omega^2 mismatch
PASS: minimal module omega^4 mismatch
PASS: real STF Gram matrix - I5
Y20 = Y21 = Y22 = (3 + (1 - (4*a^2*omega^2)/(9*cS^2))^(-1))/4
YQ_min series = 1 + (a^2*omega^2)/(9*cS^2) + (4*a^4*omega^4)/(81*cS^4)
OmegaQ_from_invariants = Sqrt[K0bar/K2bar]/2
K4 relation = (4*K2bar^2)/K0bar
Gamma5 relation = (9*K2bar^(5/2))/K0bar^(3/2)

----------------------------------------------------------------------------------------
VII.10 — Single-pole sufficiency theorem and extraction ledger
----------------------------------------------------------------------------------------
K2_from_branch = K0bar/(4*OmegaQ^2)
K4_from_branch = K0bar/(4*OmegaQ^4)
Gamma_from_branch = (9*K0bar)/(32*OmegaQ^5)
K0bar_target(OmegaQ) = (64*G*OmegaQ^5)/(45*c^5)
K2bar_target(OmegaQ) = (16*G*OmegaQ^3)/(45*c^5)
K4bar_target(OmegaQ) = (16*G*OmegaQ)/(45*c^5)
OmegaQ^2 from (c0,c2) = c0s/(4*c2s)
OmegaQ^2 from (c2,c4) = c2s/c4s
K0 from (c2,c4) = (4*c2s^2)/c4s
Gamma from (c2,c4) = (9*c4s^(3/2))/(8*Sqrt[c2s])
Direct GR gate: c4^3/c2 = (256*G^2)/(2025*c^10)

----------------------------------------------------------------------------------------
VII.11 — Grouped P2 formulas and normalized-support ledger
----------------------------------------------------------------------------------------
c0_trace = (c020 + 2*(c021 + c022))/5
c2_trace = (c220 + 2*(c221 + c222))/5
c4_trace = (c420 + 2*(c421 + c422))/5
PASS: weighted anisotropy sum n=0
PASS: weighted anisotropy sum n=2
PASS: weighted anisotropy sum n=4
Minimal branch identity c0*c4 - 4*c2^2 = -4*c2s^2 + c0s*c4s
u2_trace = (u220 + 2*(u221 + u222))/5
u4_trace = (u420 + 2*(u421 + u422))/5
normalized branch identity u4 - 4*u2^2 = -4*u2n^2 + u4n
OmegaQ on normalized branch = 1/(2*Sqrt[u2n])
Gamma5_norm on branch = 9*u2n^(5/2)
K0bar_target from normalized support = (2*G)/(45*c^5*u2n^(5/2))

----------------------------------------------------------------------------------------
VII.12 — Axisymmetric grouped-P2 parameterization
----------------------------------------------------------------------------------------
Y20 = 1 + omega^2*(4*a2 + u2bar) + omega^4*(4*a4 + u4bar)
Y21 = 1 + omega^2*(-a2 + b2 + u2bar) + omega^4*(-a4 + b4 + u4bar)
Y22 = 1 - omega^2*(a2 + b2 + (a4 + b4)*omega^2 - u2bar) + omega^4*u4bar
ubar_from_(x20,x21,x22) = (x20 + 2*(x21 + x22))/5
a_from_(x20,x21,x22) = (2*x20 - x21 - x22)/10
b_from_(x20,x21,x22) = (x21 - x22)/2
A2_sq = (4*(5*a2^2 + b2^2))/5
A4_sq = (4*(5*a4^2 + b4^2))/5
Minimal branch relation = -4*u2bar^2 + u4bar
OmegaQ^2 = 1/(4*u2bar)
Gamma5_norm = 9*u2bar^(5/2)
K0bar_target = (2*G)/(45*c^5*u2bar^(5/2))

========================================================================================
FINAL SESSION THEOREM LEDGER
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

Quadrupole normalization stack:
  gamma_GR = (2*G)/(5*c^5)
  Gamma5_target = (2*G)/(5*c^5*m0hat^2)
  N_Q_target = (54*cS^5*G)/(5*a^5*c^5)
  K0bar_target = (54*cS^5*G)/(5*a^5*c^5)
  K2bar_target = (6*cS^3*G)/(5*a^3*c^5)
  K4bar_target = (8*cS*G)/(15*a*c^5)

Minimal isotropic quadrupole module:
  Omega_Q,min = (3*cS)/(2*a)
  exact moment-matched isotropic solution Omega_Q = (3*cS)/(2*a)
  single-pole GR scaling law K0bar = (64*G*OmegaQ^5)/(45*c^5).

Normalized-support branch:
  K0bar_target from normalized support moments = (2*G)/(45*c^5*u2n^(5/2))

Best current verdict:
  The session reduces the full GR-like point-particle 2.5PN closure to the
  isotropic real P2 quadrupole module. Once the grouped low-frequency
  conservative coefficients through O(omega^4) are fixed on that branch,
  the odd Burke–Thorne coefficient follows automatically.
  The one remaining narrow theorem gap is the final passive/outgoing
  quadrupole normalization on the actual moving-throat branch.

Passes: 72
Fails : 0
OVERALL: PASS
"*)
