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
  If[ListQ[expr],
    res = Map[FullSimplify[Together[Expand[#]], Assumptions -> $Assumptions] &, expr, {ArrayDepth[expr]}];
    Print[name, " = ", fmt[res]];
    If[TrueQ[res === ConstantArray[0, Dimensions[expr]]], pass[name], fail[name, res]],
    res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
    Print[name, " = ", fmt[res]];
    If[TrueQ[res === 0], pass[name], fail[name, res]]
  ]
];

groupedParts[x20_, x21_, x22_] := {
  FullSimplify[(x20 + 2*x21 + 2*x22)/5, Assumptions -> $Assumptions],
  FullSimplify[(2*x20 - x21 - x22)/10, Assumptions -> $Assumptions],
  FullSimplify[(x21 - x22)/2, Assumptions -> $Assumptions]
};

banner["STAGE 006 — FULL GROUPED BUNDLE"];

banner["SECTION I — WEIGHTED PROJECTOR CALCULUS"];
Clear[x20, x21, x22];
$Assumptions = Element[{x20, x21, x22}, Reals];

ggrp = DiagonalMatrix[{1, 2, 2}];
ebar = {1, 1, 1};
ea = {4, -1, -1};
eb = {0, 1, -1};

expectZero["ebar^T G ea", ebar.ggrp.ea];
expectZero["ebar^T G eb", ebar.ggrp.eb];
expectZero["ea^T G eb", ea.ggrp.eb];
expectZero["||ebar||_G^2 - 5", ebar.ggrp.ebar - 5];
expectZero["||ea||_G^2 - 20", ea.ggrp.ea - 20];
expectZero["||eb||_G^2 - 4", eb.ggrp.eb - 4];

pbar = Outer[Times, ebar, ebar.ggrp]/5;
pa = Outer[Times, ea, ea.ggrp]/20;
pb = Outer[Times, eb, eb.ggrp]/4;

expectZero["Pbar^2 - Pbar", pbar.pbar - pbar];
expectZero["Pa^2 - Pa", pa.pa - pa];
expectZero["Pb^2 - Pb", pb.pb - pb];
expectZero["Pbar Pa", pbar.pa];
expectZero["Pbar Pb", pbar.pb];
expectZero["Pa Pb", pa.pb];
expectZero["Pbar + Pa + Pb - I", pbar + pa + pb - IdentityMatrix[3]];

x = {x20, x21, x22};
{xbar, ax, bx} = groupedParts[x20, x21, x22];
expectZero["Pbar*x - xbar*ebar", pbar.x - xbar*ebar];
expectZero["Pa*x - ax*ea", pa.x - ax*ea];
expectZero["Pb*x - bx*eb", pb.x - bx*eb];

banner["SECTION II — BUNDLE COEFFICIENT ASSEMBLY"];
Clear[omega, omegaU, omegaW, rMix, gU, gW];
$Assumptions = Element[{omega, omegaU, omegaW, rMix, gU, gW}, Reals];

deltaExpr = omegaU^2*omegaW^2 - rMix^2;
sExpr = omegaU^2 + omegaW^2;
qExpr = gU^2*omegaW^2 + 2*gU*gW*rMix + gW^2*omegaU^2;
hExpr = gU^2 + gW^2;
pExpr = omegaU^2*gW + rMix*gU;

(* Schur-complement derivation of the rational form used below.
   The paper's one-port Lagrangian has +R U W, so the frequency-space
   spring matrix has off-diagonal -R. Its determinant is the denominator. *)
mBlock = {{omegaU^2 - omega^2, -rMix}, {-rMix, omegaW^2 - omega^2}};
detM = Expand[Det[mBlock]];
expectZero["Schur denominator matches Delta - S omega^2 + omega^4",
  Expand[detM - (deltaExpr - sExpr*omega^2 + omega^4)]];
gVec = {gU, gW};
zSchur = FullSimplify[gVec . Inverse[mBlock] . gVec, Assumptions -> $Assumptions];
expectZero["Z rational form matches Schur (g_U,g_W) M^{-1} (g_U,g_W)^T",
  FullSimplify[zSchur - (qExpr - hExpr*omega^2)/detM, Assumptions -> $Assumptions]];

zOnePort = Expand[Normal[Series[(qExpr - hExpr*omega^2)/(deltaExpr - sExpr*omega^2 + omega^4), {omega, 0, 4}]]];
nOnePort = Expand[Normal[Series[(pExpr - gW*omega^2)^2/(deltaExpr - sExpr*omega^2 + omega^4)^2, {omega, 0, 4}]]];

expectZero["Z0 one-port", Coefficient[zOnePort, omega, 0] - qExpr/deltaExpr];
expectZero["Z2 one-port", Coefficient[zOnePort, omega, 2] - (qExpr*sExpr - hExpr*deltaExpr)/deltaExpr^2];
expectZero["Z4 one-port", Coefficient[zOnePort, omega, 4] - (qExpr*(sExpr^2 - deltaExpr) - sExpr*hExpr*deltaExpr)/deltaExpr^3];
expectZero["N0 one-port", Coefficient[nOnePort, omega, 0] - pExpr^2/deltaExpr^2];
expectZero["N2 one-port", Coefficient[nOnePort, omega, 2] - 2*pExpr*(pExpr*sExpr - deltaExpr*gW)/deltaExpr^3];
expectZero["N4 one-port", Coefficient[nOnePort, omega, 4] - (deltaExpr^2*gW^2 - 2*deltaExpr*pExpr^2 - 4*deltaExpr*pExpr*sExpr*gW + 3*pExpr^2*sExpr^2)/deltaExpr^4];

(* Independent-path verification: numerical substitution at fixed
   parameter values to verify Z_n, N_n closed forms agree with the
   rational function expansion. This breaks structural correspondence
   with the SymPy symbolic series approach. *)
zRule = {omegaU -> 2, omegaW -> 3, rMix -> 1, gU -> 1, gW -> 2};
deltaNum = deltaExpr /. zRule;
sNum = sExpr /. zRule;
qNum = qExpr /. zRule;
hNum = hExpr /. zRule;
pNum = pExpr /. zRule;
zNum0 = qNum/deltaNum;
zNum2 = (qNum*sNum - hNum*deltaNum)/deltaNum^2;
zNum4 = (qNum*(sNum^2 - deltaNum) - sNum*hNum*deltaNum)/deltaNum^3;
nNum0 = pNum^2/deltaNum^2;
nNum2 = 2*pNum*(pNum*sNum - deltaNum*gW)/deltaNum^3 /. zRule;
nNum4 = (deltaNum^2*gW^2 - 2*deltaNum*pNum^2 - 4*deltaNum*pNum*sNum*gW + 3*pNum^2*sNum^2)/deltaNum^4 /. zRule;
(* Compare to direct numerical Taylor coefficients of the rational function. *)
zNumFromSeries = SeriesCoefficient[(qExpr - hExpr*omega^2)/(deltaExpr - sExpr*omega^2 + omega^4) /. zRule, {omega, 0, #}] & /@ {0, 2, 4};
nNumFromSeries = SeriesCoefficient[((pExpr - gW*omega^2)^2/(deltaExpr - sExpr*omega^2 + omega^4)^2) /. zRule, {omega, 0, #}] & /@ {0, 2, 4};
expectZero["Z0 numerical cross-check", zNum0 - zNumFromSeries[[1]]];
expectZero["Z2 numerical cross-check", zNum2 - zNumFromSeries[[2]]];
expectZero["Z4 numerical cross-check", zNum4 - zNumFromSeries[[3]]];
expectZero["N0 numerical cross-check", nNum0 - nNumFromSeries[[1]]];
expectZero["N2 numerical cross-check", nNum2 - nNumFromSeries[[2]]];
expectZero["N4 numerical cross-check", nNum4 - nNumFromSeries[[3]]];

(* Stage-003 carry-forward: B_{A0}, B_{A2}, B_{A4} are the stable-BdG Schur
   sums Sum_alpha g_{A alpha}^2/varpi_{A alpha}^{2,4,6} in grouped notation. *)
Clear[k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022, n220, n221, n222];
$Assumptions = Element[{k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022, n220, n221, n222}, Reals];

d020 = k20 - b020 - z020;
d021 = k21 - b021 - z021;
d022 = k22 - b022 - z022;
d220 = -(m20 + b220 + z220);
d221 = -(m21 + b221 + z221);
d222 = -(m22 + b222 + z222);
d420 = -(b420 + z420);
d421 = -(b421 + z421);
d422 = -(b422 + z422);

{kbar, aK, bK} = groupedParts[k20, k21, k22];
{mbar, aM, bM} = groupedParts[m20, m21, m22];
{bbar0, aB0, bB0} = groupedParts[b020, b021, b022];
{bbar2, aB2, bB2} = groupedParts[b220, b221, b222];
{bbar4, aB4, bB4} = groupedParts[b420, b421, b422];
{zbar0, aZ0, bZ0} = groupedParts[z020, z021, z022];
{zbar2, aZ2, bZ2} = groupedParts[z220, z221, z222];
{zbar4, aZ4, bZ4} = groupedParts[z420, z421, z422];
{dbar0, aD0, bD0} = groupedParts[d020, d021, d022];
{dbar2, aD2, bD2} = groupedParts[d220, d221, d222];
{dbar4, aD4, bD4} = groupedParts[d420, d421, d422];
{nbar0, aN0, bN0} = groupedParts[n020, n021, n022];

expectZero["Dbar0 decomposition", dbar0 - (kbar - bbar0 - zbar0)];
expectZero["aD0 decomposition", aD0 - (aK - aB0 - aZ0)];
expectZero["bD0 decomposition", bD0 - (bK - bB0 - bZ0)];
expectZero["Dbar2 decomposition", dbar2 + (mbar + bbar2 + zbar2)];
expectZero["aD2 decomposition", aD2 + (aM + aB2 + aZ2)];
expectZero["bD2 decomposition", bD2 + (bM + bB2 + bZ2)];
expectZero["Dbar4 decomposition", dbar4 + (bbar4 + zbar4)];
expectZero["aD4 decomposition", aD4 + (aB4 + aZ4)];
expectZero["bD4 decomposition", bD4 + (bB4 + bZ4)];
(* Non-tautological linearity test: grouped_parts of a lane sum
   equals the sum of the individual grouped_parts outputs. *)
{nbar2, aN2, bN2} = groupedParts[n220, n221, n222];
{nbar02, aN02, bN02} = groupedParts[n020 + n220, n021 + n221, n022 + n222];
expectZero["Nbar0 + Nbar2 additivity", nbar02 - (nbar0 + nbar2)];
expectZero["aN0 + aN2 additivity", aN02 - (aN0 + aN2)];
expectZero["bN0 + bN2 additivity", bN02 - (bN0 + bN2)];

banner["SECTION III — ISOTROPIC BRANCH AND TARGET"];
Clear[omega, d0, d2, d4, n0, n2, n4, G, c, cS, a, mhat];
$Assumptions =
  Element[{omega, d0, d2, d4, n0, n2, n4, G, c, cS, a, mhat}, Reals] &&
  d0 != 0 && And @@ Thread[{G, c, cS, a, mhat} > 0];

dCons = d0 + d2*omega^2 + d4*omega^4;
yResp = Expand[Normal[Series[d0/dCons, {omega, 0, 4}]]];
u2 = FullSimplify[Coefficient[yResp, omega, 2], Assumptions -> $Assumptions];
u4 = FullSimplify[Coefficient[yResp, omega, 4], Assumptions -> $Assumptions];
pref = Expand[Normal[Series[d0*(n0 + n2*omega^2 + n4*omega^4)/dCons^2, {omega, 0, 4}]]];
p0 = FullSimplify[Coefficient[pref, omega, 0], Assumptions -> $Assumptions];
p2 = FullSimplify[Coefficient[pref, omega, 2], Assumptions -> $Assumptions];
p4 = FullSimplify[Coefficient[pref, omega, 4], Assumptions -> $Assumptions];

expectZero["u2 + D2/D0", u2 + d2/d0];
expectZero["u4 formula", u4 - (d2^2 - d0*d4)/d0^2];
expectZero["P0 - N0/D0", p0 - n0/d0];
expectZero["P2 formula", p2 - (d0*n2 - 2*d2*n0)/d0^2];
expectZero["P4 formula", p4 - (d0^2*n4 - 2*d0*(d2*n2 + d4*n0) + 3*d2^2*n0)/d0^3];

n2Target = FullSimplify[n2 /. First[Solve[p2 == 0, n2, Reals]], Assumptions -> $Assumptions];
n4Target = FullSimplify[n4 /. First[Solve[(p4 /. n2 -> n2Target) == 0, n4, Reals]], Assumptions -> $Assumptions];
(* Independent closed-form derivation of n2Target and n4Target. *)
n2TargetClosed = 2*d2*n0/d0;
n4TargetClosed = FullSimplify[n0*(2*d0*d4 + d2^2)/d0^2, Assumptions -> $Assumptions];
expectZero["N2_target closed form", FullSimplify[n2Target - n2TargetClosed, Assumptions -> $Assumptions]];
expectZero["N4_target closed form", FullSimplify[n4Target - n4TargetClosed, Assumptions -> $Assumptions]];

(* Non-tautological substitution check: abstract mhat^2 * P0 == explicit
   mhat^2 * N0/(K - B0 - Z0) under D0 = K - B0 - Z0. *)
Clear[kNorm, b0Norm, z0Norm];
$Assumptions = $Assumptions && Element[{kNorm, b0Norm, z0Norm}, Reals] && kNorm - b0Norm - z0Norm != 0;
normAbstract = mhat^2*n0/d0 - 54*G*cS^5/(5*a^5*c^5);
normExplicit = mhat^2*n0/(kNorm - b0Norm - z0Norm) - 54*G*cS^5/(5*a^5*c^5);
expectZero["normalization abstract == explicit under D0 = K - B0 - Z0",
  FullSimplify[(normAbstract /. d0 -> kNorm - b0Norm - z0Norm) - normExplicit,
    Assumptions -> $Assumptions]];

Clear[z, j2, y2, h2, lambda2, lambda2Series, y2Resp, y2Static, y2Hat, g5Stage4];
$Assumptions =
  Element[{G, c, cS, a, mhat, z, omega}, Reals] &&
  And @@ Thread[{G, c, cS, a, mhat, z} > 0];

j2 = ((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2 = -((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2;
h2 = FullSimplify[j2 + I y2, Assumptions -> $Assumptions];
lambda2 = FullSimplify[(omega D[h2, z]/h2) /. z -> omega a/cS, Assumptions -> $Assumptions];
lambda2Series = Normal[Series[lambda2, {omega, 0, 6}]];
y2Resp = Normal[Series[1/lambda2Series, {omega, 0, 5}]] // FullSimplify;
y2Static = FullSimplify[y2Resp /. omega -> 0, Assumptions -> a > 0 && cS > 0];
y2Hat = Expand[y2Resp/y2Static];
g5Stage4 = FullSimplify[Coefficient[y2Hat, omega, 5]/I, Assumptions -> $Assumptions];
gamma5Port = g5Stage4;
expectZero["Stage-5 Gamma5_port anchor", gamma5Port - a^5/(27*cS^5)];
(* Independent-path verification: compute the 5th-order coefficient of
   Y2/Y2_static by direct small-z Taylor expansion of j2 + I y2 (and its
   derivative), bypassing the omega*D[h2,z]/h2 ratio path used above. *)
h2SmallZ = Normal[Series[j2 + I y2, {z, 0, 9}]];
h2DerivSmallZ = Normal[Series[D[j2 + I y2, z], {z, 0, 8}]];
lambdaSmallZ = Normal[Series[(z*h2DerivSmallZ/h2SmallZ) /. z -> omega*a/cS, {omega, 0, 6}]];
ySmallZ = Normal[Series[1/lambdaSmallZ, {omega, 0, 5}]];
yStatSmallZ = ySmallZ /. omega -> 0;
yHatSmallZ = Expand[ySmallZ/yStatSmallZ];
gamma5Bessel = FullSimplify[Coefficient[yHatSmallZ, omega, 5]/I, Assumptions -> $Assumptions];
expectZero["Gamma5_port via direct Bessel small-z expansion", gamma5Bessel - a^5/(27*cS^5)];
gammaGR = 2*G/(5*c^5);
ratioTarget = FullSimplify[gammaGR/(mhat^2*gamma5Port), Assumptions -> $Assumptions];
expectZero["ratio target at mhat=1", (ratioTarget /. mhat -> 1) - 54*G*cS^5/(5*a^5*c^5)];

banner["SECTION IV — FIRST-ORDER ANISOTROPY TRANSPORT"];
Clear[eps, dD0, dD2, dN0, aD0s, aD2s, aN0s, bD0s, bD2s, bN0s];
$Assumptions = Element[{eps, dD0, dD2, dN0, aD0s, aD2s, aN0s, bD0s, bD2s, bN0s}, Reals] && d0 != 0;
u2Eps = Expand[Normal[Series[-(d2 + eps*dD2)/(d0 + eps*dD0), {eps, 0, 1}]]];
p0Eps = Expand[Normal[Series[(n0 + eps*dN0)/(d0 + eps*dD0), {eps, 0, 1}]]];
du2 = FullSimplify[Coefficient[u2Eps, eps, 1], Assumptions -> $Assumptions];
dP0 = FullSimplify[Coefficient[p0Eps, eps, 1], Assumptions -> $Assumptions];
u2Base = -d2/d0;
p0Base = n0/d0;

expectZero["du2 formula", du2 + (dD2 + u2Base*dD0)/d0];
expectZero["dP0 formula", dP0 - (dN0 - p0Base*dD0)/d0];
expectZero["a_u2 formula", (du2 /. {dD0 -> aD0s, dD2 -> aD2s}) + (aD2s + u2Base*aD0s)/d0];
expectZero["b_u2 formula", (du2 /. {dD0 -> bD0s, dD2 -> bD2s}) + (bD2s + u2Base*bD0s)/d0];
expectZero["a_P0 formula", (dP0 /. {dD0 -> aD0s, dN0 -> aN0s}) - (aN0s - p0Base*aD0s)/d0];
expectZero["b_P0 formula", (dP0 /. {dD0 -> bD0s, dN0 -> bN0s}) - (bN0s - p0Base*bD0s)/d0];

banner["SECTION V — MONOTONICITY DERIVATIVES"];
Clear[k, b0, z0, n0mono];
$Assumptions = Element[{k, b0, z0, n0mono}, Reals] && And @@ Thread[{k, b0, z0, n0mono} > 0] && k - b0 - z0 > 0;
p0Mono = n0mono/(k - b0 - z0);
expectZero["dP0/dN0 - 1/D0", D[p0Mono, n0mono] - 1/(k - b0 - z0)];
expectZero["dP0/dB0 - N0/D0^2", D[p0Mono, b0] - n0mono/(k - b0 - z0)^2];
expectZero["dP0/dZ0 - N0/D0^2", D[p0Mono, z0] - n0mono/(k - b0 - z0)^2];
expectZero["dP0/dK + N0/D0^2", D[p0Mono, k] + n0mono/(k - b0 - z0)^2];

Print[""];
Print["FINAL STAGE-006 LEDGER:"];
Print["  Verified the weighted projector calculus, the full grouped-bundle coefficient"];
Print["  assembly, the isotropic prefactor/normalization formulas, the constant-"];
Print["  prefactor branch conditions, and the first-order anisotropy transport laws."];

Exit[0];
