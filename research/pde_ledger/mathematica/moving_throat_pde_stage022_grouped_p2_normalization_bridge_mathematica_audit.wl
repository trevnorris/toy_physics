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

banner["STAGE 022 — GROUPED P2 NORMALIZATION BRIDGE"];

Clear[omega, d0, d2, d4, n0, n2, n4, aa, bb, g5];
$Assumptions =
  Element[{omega, d0, d2, d4, n0, n2, n4, aa, bb, g5}, Reals] &&
  d0 != 0;

dCons = d0 + d2*omega^2 + d4*omega^4;
(* Inverse-relation route: assume Y_resp = 1 + u2Sym omega^2 + u4Sym omega^4 + O(omega^6)
   and impose Y_resp * dCons - d0 = 0 modulo omega^6. *)
Clear[u2Sym, u4Sym];
yRespCand = 1 + u2Sym*omega^2 + u4Sym*omega^4;
prod = Expand[yRespCand*dCons - d0];
coeffEqs = Table[Coefficient[prod, omega, k] == 0, {k, 0, 4}];
sol = First[Solve[coeffEqs, {u2Sym, u4Sym}]];
u2 = FullSimplify[u2Sym /. sol, Assumptions -> $Assumptions];
u4 = FullSimplify[u4Sym /. sol, Assumptions -> $Assumptions];

banner["SECTION I — CONSERVATIVE OPERATOR TO RESPONSE"];
Print["u2 (inverse-relation route) = ", fmt[u2]];
Print["u4 (inverse-relation route) = ", fmt[u4]];
expectZero["u2 formula", u2 + d2/d0];
expectZero["u4 formula", u4 - (d2^2 - d0*d4)/d0^2];

nFac = n0 + n2*omega^2 + n4*omega^4;
(* Inverse-relation route: assume pref = p0Sym + p2Sym omega^2 + p4Sym omega^4 + O(omega^6)
   and impose pref * dCons^2 - d0 * nFac = 0 modulo omega^6. *)
Clear[p0Sym, p2Sym, p4Sym];
prefCand = p0Sym + p2Sym*omega^2 + p4Sym*omega^4;
prodPref = Expand[prefCand*dCons^2 - d0*nFac];
coeffEqsPref = Table[Coefficient[prodPref, omega, k] == 0, {k, 0, 4}];
solPref = First[Solve[coeffEqsPref, {p0Sym, p2Sym, p4Sym}]];
p0 = FullSimplify[p0Sym /. solPref, Assumptions -> $Assumptions];
p2 = FullSimplify[p2Sym /. solPref, Assumptions -> $Assumptions];
p4 = FullSimplify[p4Sym /. solPref, Assumptions -> $Assumptions];

(* Branch coefficients: direct polynomial expansion (Expand only, no Series). *)
y2Out = 1 + aa*omega^2 + bb*omega^4 + I*g5*omega^5;
prefTrunc = p0 + p2*omega^2 + p4*omega^4;
yBranch = Expand[prefTrunc*y2Out];
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
(* Intentional parallel check: the 3x3 inverse-map identity admits no
   engine-independent route. Both engines verify the same algebra as a
   sanity cross-check. *)
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

dProto = delta0 - s2*omega^2 + omega^4;
pProto = p0proto - gW*omega^2;
(* Inverse-relation route: assume nSeries = n0Cand + n2Cand omega^2 + n4Cand omega^4
   and impose nSeries * dProto^2 - pProto^2 = 0 modulo omega^6. *)
Clear[n0Cand, n2Cand, n4Cand];
nCand = n0Cand + n2Cand*omega^2 + n4Cand*omega^4;
prodN = Expand[nCand*dProto^2 - pProto^2];
coeffEqsN = Table[Coefficient[prodN, omega, k] == 0, {k, 0, 4}];
solN = First[Solve[coeffEqsN, {n0Cand, n2Cand, n4Cand}]];
n0Proto = FullSimplify[n0Cand /. solN, Assumptions -> $Assumptions];
n2Proto = FullSimplify[n2Cand /. solN, Assumptions -> $Assumptions];
n4Proto = FullSimplify[n4Cand /. solN, Assumptions -> $Assumptions];

banner["SECTION IV — STAGE-021 ONE-PORT PROTOTYPE"];
expectZero["N0 prototype", n0Proto - p0proto^2/delta0^2];
expectZero["N2 prototype", n2Proto - 2*p0proto*(p0proto*s2 - delta0*gW)/delta0^3];
expectZero[
  "N4 prototype",
  n4Proto - (delta0^2*gW^2 - 2*delta0*p0proto^2 - 4*delta0*p0proto*s2*gW + 3*p0proto^2*s2^2)/delta0^4
];

deltaBack = omegaA^2*omegaW^2 - rMix^2;
s2Back = omegaA^2 + omegaW^2;
p0Back = omegaA^2*gW + rMix*gA;

Clear[G, c, cS, a, mhat, p0Static];
$Assumptions =
  Element[{G, c, cS, a, mhat, p0Static}, Reals] &&
  And @@ Thread[{G, c, cS, a, mhat, p0Static} > 0];

Clear[z, j2, y2, h2, lambda2, lambda2Series, y2, y2Static, y2Hat, y2HatOmega, aStage4, bStage4, g5Stage4];
$Assumptions =
  Element[{G, c, cS, a, mhat, p0Static, z, omega}, Reals] &&
  And @@ Thread[{G, c, cS, a, mhat, p0Static, z} > 0];

(* Use Mathematica's built-in spherical Hankel function instead of the
   explicit polynomial-rational form, so the derivation of A, B, G5 is
   independent of the SymPy script's choice of j2, y2 expressions. *)
h2 = SphericalHankelH1[2, z];
lambda2 = FullSimplify[(omega D[h2, z]/h2) /. z -> omega a/cS, Assumptions -> $Assumptions];
lambda2Series = Normal[Series[lambda2, {omega, 0, 6}]];
y2Resp = Normal[Series[1/lambda2Series, {omega, 0, 5}]] // FullSimplify;
y2Static = FullSimplify[y2Resp /. omega -> 0, Assumptions -> a > 0 && cS > 0];
y2Hat = Expand[y2Resp/y2Static];
aStage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 2], Assumptions -> $Assumptions];
bStage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 4], Assumptions -> $Assumptions];
g5Stage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 5]/I, Assumptions -> $Assumptions];

banner["SECTION V — STAGE-021 OUTGOING l=2 FINGERPRINT ANCHOR"];
expectZero["Stage-021 A coefficient", aStage4 - a^2/(9*cS^2)];
expectZero["Stage-021 B coefficient", bStage4 - 4*a^4/(81*cS^4)];
expectZero["Stage-021 G5 coefficient", g5Stage4 - a^5/(27*cS^5)];

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
Print["FINAL STAGE-022 LEDGER:"];
Print["  Verified the exact grouped P2 bridge from conservative operator moments to"];
Print["  normalized response moments, outgoing prefactor coefficients, the one-port"];
Print["  prototype N0/N2/N4 formulas, and the universal quadrupole normalization product."];

Exit[0];
