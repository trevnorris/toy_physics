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
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 005 — GROUPED P2 NORMALIZATION BRIDGE"];

Clear[omega, d0, d2, d4, n0, n2, n4, aa, bb, g5];
$Assumptions =
  Element[{omega, d0, d2, d4, n0, n2, n4, aa, bb, g5}, Reals] &&
  d0 != 0;

dCons = d0 + d2*omega^2 + d4*omega^4;
yResp = Expand[Normal[Series[d0/dCons, {omega, 0, 4}]]];
u2 = FullSimplify[Coefficient[yResp, omega, 2], Assumptions -> $Assumptions];
u4 = FullSimplify[Coefficient[yResp, omega, 4], Assumptions -> $Assumptions];

banner["SECTION I — CONSERVATIVE OPERATOR TO RESPONSE"];
Print["Y_resp(omega) = ", fmt[yResp]];
expectZero["u2 formula", u2 + d2/d0];
expectZero["u4 formula", u4 - (d2^2 - d0*d4)/d0^2];

nFac = n0 + n2*omega^2 + n4*omega^4;
pref = Expand[Normal[Series[d0*nFac/dCons^2, {omega, 0, 4}]]];
p0 = FullSimplify[Coefficient[pref, omega, 0], Assumptions -> $Assumptions];
p2 = FullSimplify[Coefficient[pref, omega, 2], Assumptions -> $Assumptions];
p4 = FullSimplify[Coefficient[pref, omega, 4], Assumptions -> $Assumptions];

y2Out = 1 + aa*omega^2 + bb*omega^4 + I*g5*omega^5;
yBranch = Expand[pref*y2Out];
k0 = FullSimplify[Coefficient[yBranch, omega, 0], Assumptions -> $Assumptions];
k2 = FullSimplify[Coefficient[yBranch, omega, 2], Assumptions -> $Assumptions];
k4 = FullSimplify[Coefficient[yBranch, omega, 4], Assumptions -> $Assumptions];
gamma5 = FullSimplify[Coefficient[yBranch, omega, 5]/I, Assumptions -> $Assumptions];

banner["SECTION II — OUTGOING PREFACTOR AND BRANCH COEFFICIENTS"];
expectZero["P0 formula", p0 - n0/d0];
expectZero["P2 formula", p2 - (d0*n2 - 2*d2*n0)/d0^2];
expectZero[
  "P4 formula",
  p4 - (d0^2*n4 - 2*d0*(d2*n2 + d4*n0) + 3*d2^2*n0)/d0^3
];
expectZero["K0 formula", k0 - p0];
expectZero["K2 formula", k2 - (p2 + aa*p0)];
expectZero["K4 formula", k4 - (p4 + aa*p2 + bb*p0)];
expectZero["Gamma5 formula", gamma5 - g5*p0];

Clear[x20, x21, x22, xQ];
$Assumptions = Element[{x20, x21, x22, xQ}, Reals];
xbar = FullSimplify[(x20 + 2*x21 + 2*x22)/5, Assumptions -> $Assumptions];
ax = FullSimplify[(2*x20 - x21 - x22)/10, Assumptions -> $Assumptions];
bx = FullSimplify[(x21 - x22)/2, Assumptions -> $Assumptions];

banner["SECTION III — GROUPED REAL P2 INVERSE MAP"];
expectZero["x20 recovered", (xbar + 4*ax) - x20];
expectZero["x21 recovered", (xbar - ax + bx) - x21];
expectZero["x22 recovered", (xbar - ax - bx) - x22];
expectZero["xbar isotropic", (xbar /. {x20 -> xQ, x21 -> xQ, x22 -> xQ}) - xQ];
expectZero["ax isotropic", ax /. {x20 -> xQ, x21 -> xQ, x22 -> xQ}];
expectZero["bx isotropic", bx /. {x20 -> xQ, x21 -> xQ, x22 -> xQ}];

Clear[delta0, s2, p0proto, gW, omegaA, omegaW, rMix, gA];
$Assumptions =
  Element[{omega, delta0, s2, p0proto, gW, omegaA, omegaW, rMix, gA}, Reals] &&
  delta0 != 0;

nProto = (p0proto - gW*omega^2)^2/(delta0 - s2*omega^2 + omega^4)^2;
nSeries = Expand[Normal[Series[nProto, {omega, 0, 4}]]];
n0Proto = FullSimplify[Coefficient[nSeries, omega, 0], Assumptions -> $Assumptions];
n2Proto = FullSimplify[Coefficient[nSeries, omega, 2], Assumptions -> $Assumptions];
n4Proto = FullSimplify[Coefficient[nSeries, omega, 4], Assumptions -> $Assumptions];

banner["SECTION IV — STAGE-4 ONE-PORT PROTOTYPE"];
expectZero["N0 prototype", n0Proto - p0proto^2/delta0^2];
expectZero["N2 prototype", n2Proto - 2*p0proto*(p0proto*s2 - delta0*gW)/delta0^3];
expectZero[
  "N4 prototype",
  n4Proto - (delta0^2*gW^2 - 2*delta0*p0proto^2 - 4*delta0*p0proto*s2*gW + 3*p0proto^2*s2^2)/delta0^4
];

deltaBack = omegaA^2*omegaW^2 - rMix^2;
s2Back = omegaA^2 + omegaW^2;
p0Back = omegaA^2*gW + rMix*gA;
nStage4 = Expand[
  Normal[Series[((p0Back - gW*omega^2)^2)/(deltaBack - s2Back*omega^2 + omega^4)^2, {omega, 0, 4}]]
];
expectZero["N0 round-trip", Coefficient[nStage4, omega, 0] - (n0Proto /. {delta0 -> deltaBack, s2 -> s2Back, p0proto -> p0Back})];
expectZero["N2 round-trip", Coefficient[nStage4, omega, 2] - (n2Proto /. {delta0 -> deltaBack, s2 -> s2Back, p0proto -> p0Back})];
expectZero["N4 round-trip", Coefficient[nStage4, omega, 4] - (n4Proto /. {delta0 -> deltaBack, s2 -> s2Back, p0proto -> p0Back})];

Clear[G, c, cS, a, mhat, p0Static];
$Assumptions =
  Element[{G, c, cS, a, mhat, p0Static}, Reals] &&
  And @@ Thread[{G, c, cS, a, mhat, p0Static} > 0];

Clear[z, j2, y2, h2, lambda2, lambda2Series, y2, y2Static, y2Hat, y2HatOmega, aStage4, bStage4, g5Stage4];
$Assumptions =
  Element[{G, c, cS, a, mhat, p0Static, z, omega}, Reals] &&
  And @@ Thread[{G, c, cS, a, mhat, p0Static, z} > 0];

j2 = ((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2 = -((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2;
h2 = FullSimplify[j2 + I y2, Assumptions -> $Assumptions];
lambda2 = FullSimplify[(omega D[h2, z]/h2) /. z -> omega a/cS, Assumptions -> $Assumptions];
lambda2Series = Normal[Series[lambda2, {omega, 0, 6}]];
y2Resp = Normal[Series[1/lambda2Series, {omega, 0, 5}]] // FullSimplify;
y2Static = FullSimplify[y2Resp /. omega -> 0, Assumptions -> a > 0 && cS > 0];
y2Hat = Expand[y2Resp/y2Static];
aStage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 2], Assumptions -> $Assumptions];
bStage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 4], Assumptions -> $Assumptions];
g5Stage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 5]/I, Assumptions -> $Assumptions];

banner["SECTION V — STAGE-4 OUTGOING l=2 FINGERPRINT ANCHOR"];
expectZero["Stage-4 A coefficient", aStage4 - a^2/(9*cS^2)];
expectZero["Stage-4 B coefficient", bStage4 - 4*a^4/(81*cS^4)];
expectZero["Stage-4 G5 coefficient", g5Stage4 - a^5/(27*cS^5)];

gamma5Port = g5Stage4;
gammaGR = FullSimplify[2*G/(5*c^5), Assumptions -> $Assumptions];
ratioTarget = FullSimplify[p0Static /. First[Solve[mhat^2*p0Static*gamma5Port == gammaGR, p0Static, Reals]], Assumptions -> $Assumptions];
k2Target = FullSimplify[ratioTarget*aStage4, Assumptions -> $Assumptions];
k4Target = FullSimplify[ratioTarget*bStage4, Assumptions -> $Assumptions];

banner["SECTION VI — QUADRUPOLE NORMALIZATION PRODUCT"];
Print["Gamma5_port = ", fmt[gamma5Port]];
Print["Required mhat^2 * P0 = ", fmt[ratioTarget]];
expectZero["mhat=1 K0 target", (ratioTarget /. mhat -> 1) - 54*G*cS^5/(5*a^5*c^5)];
expectZero["mhat=1 K2 target", (k2Target /. mhat -> 1) - 6*G*cS^3/(5*a^3*c^5)];
expectZero["mhat=1 K4 target", (k4Target /. mhat -> 1) - 8*G*cS/(15*a*c^5)];

Print[""];
Print["FINAL STAGE-005 LEDGER:"];
Print["  Verified the exact grouped P2 bridge from conservative operator moments to"];
Print["  normalized response moments, outgoing prefactor coefficients, the one-port"];
Print["  prototype N0/N2/N4 formulas, and the universal quadrupole normalization product."];

Exit[0];
