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

zOnePort = Expand[Normal[Series[(qExpr - hExpr*omega^2)/(deltaExpr - sExpr*omega^2 + omega^4), {omega, 0, 4}]]];
nOnePort = Expand[Normal[Series[(pExpr - gW*omega^2)^2/(deltaExpr - sExpr*omega^2 + omega^4)^2, {omega, 0, 4}]]];

expectZero["Z0 one-port", Coefficient[zOnePort, omega, 0] - qExpr/deltaExpr];
expectZero["Z2 one-port", Coefficient[zOnePort, omega, 2] - (qExpr*sExpr - hExpr*deltaExpr)/deltaExpr^2];
expectZero["Z4 one-port", Coefficient[zOnePort, omega, 4] - (qExpr*(sExpr^2 - deltaExpr) - sExpr*hExpr*deltaExpr)/deltaExpr^3];
expectZero["N0 one-port", Coefficient[nOnePort, omega, 0] - pExpr^2/deltaExpr^2];
expectZero["N2 one-port", Coefficient[nOnePort, omega, 2] - 2*pExpr*(pExpr*sExpr - deltaExpr*gW)/deltaExpr^3];
expectZero["N4 one-port", Coefficient[nOnePort, omega, 4] - (deltaExpr^2*gW^2 - 2*deltaExpr*pExpr^2 - 4*deltaExpr*pExpr*sExpr*gW + 3*pExpr^2*sExpr^2)/deltaExpr^4];

Clear[k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022];
$Assumptions = Element[{k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022}, Reals];

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
expectZero["Nbar0 formula", nbar0 - (n020 + 2*n021 + 2*n022)/5];
expectZero["aN0 formula", aN0 - (2*n020 - n021 - n022)/10];
expectZero["bN0 formula", bN0 - (n021 - n022)/2];

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
expectZero["P2 under N2 target", p2 /. n2 -> n2Target];
expectZero["P4 under N2,N4 targets", p4 /. {n2 -> n2Target, n4 -> n4Target}];

gamma5Port = a^5/(27*cS^5);
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
