ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
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

Needs["VariationalMethods`"];

banner["STAGE 004 — MAXWELL + MIXED-SECTOR REDUCTION"];

subbanner["I. Mixed 4+1 Maxwell fields are gauge invariant"];

Clear[t, x, w, chiFun, a0Fun, aaFun, awFun];
$Assumptions = Element[{t, x, w}, Reals];

chi = chiFun[t, x, w];
a0 = a0Fun[t, x, w];
aa = aaFun[t, x, w];
aw = awFun[t, x, w];

eW = -D[aw, t] - D[a0, w];
cA = D[aw, x] - D[aa, w];
a0p = a0 - D[chi, t];
aap = aa + D[chi, x];
awp = aw + D[chi, w];
eWp = -D[awp, t] - D[a0p, w];
cAp = D[awp, x] - D[aap, w];

expectZero["E_w gauge variation", eWp - eW];
expectZero["C_a gauge variation", cAp - cA];

subbanner["II. Conservative wall + brane-MAXWELL + mixed-mode reduction"];

Clear[omega, qFun, aFun, wFun, m, k, oA, oW, r, gA, gW, d0, s2, n0, g2];
$Assumptions =
  Element[{t, omega, m, k, oA, oW, r, gA, gW, d0, s2, n0, g2}, Reals] &&
  m > 0 && k > 0 && oA > 0 && oW > 0 && d0 > 0 && s2 > 0;

q = qFun[t];
a = aFun[t];
ww = wFun[t];

lRed = (
  1/2 m D[q, t]^2 - 1/2 k q^2
  + 1/2 D[a, t]^2 - 1/2 oA^2 a^2
  + 1/2 D[ww, t]^2 - 1/2 oW^2 ww^2
  + r a ww + gA q a + gW q ww
);

elList = EulerEquations[lRed, {qFun[t], aFun[t], wFun[t]}, t];
expectZero["Q equation", (elList[[1]] /. Equal[lhs_, rhs_] :> lhs - rhs) + (m qFun''[t] + k qFun[t] - gA aFun[t] - gW wFun[t])];
expectZero["A equation", (elList[[2]] /. Equal[lhs_, rhs_] :> lhs - rhs) + (aFun''[t] + oA^2 aFun[t] - r wFun[t] - gA qFun[t])];
expectZero["W equation", (elList[[3]] /. Equal[lhs_, rhs_] :> lhs - rhs) + (wFun''[t] + oW^2 wFun[t] - r aFun[t] - gW qFun[t])];

aKer = oA^2 - omega^2;
wKer = oW^2 - omega^2;
delta = FullSimplify[aKer wKer - r^2, Assumptions -> $Assumptions];
sigmaCons = FullSimplify[(gA^2 wKer + 2 gA gW r + gW^2 aKer)/delta, Assumptions -> $Assumptions];
dCons = FullSimplify[k - m omega^2 - sigmaCons, Assumptions -> $Assumptions];
matEAW = {{aKer, -r}, {-r, wKer}};
solAW = LinearSolve[matEAW, {gA, gW}];
aSol = FullSimplify[solAW[[1]], Assumptions -> $Assumptions];
wSol = FullSimplify[solAW[[2]], Assumptions -> $Assumptions];
sigmaConsDerived = FullSimplify[gA aSol + gW wSol, Assumptions -> $Assumptions];
expectZero["sigmaCons from LinearSolve matches closed form", sigmaConsDerived - sigmaCons];

Print["Sigma_EM+mix^cons(omega) = ", fmt[sigmaCons]];
Print["D_cons(omega) = ", fmt[dCons]];
expectZero["A exact solution residual", aKer aSol - r wSol - gA];
expectZero["W exact solution residual", wKer wSol - r aSol - gW];

toy = (n0 - g2 omega^2)/(d0 - s2 omega^2 + omega^4);
toySeries = Expand[Normal[Series[toy, {omega, 0, 4}]]];
z0Toy = FullSimplify[Coefficient[toySeries, omega, 0], Assumptions -> $Assumptions];
z2Toy = FullSimplify[Coefficient[toySeries, omega, 2], Assumptions -> $Assumptions];
z4Toy = FullSimplify[Coefficient[toySeries, omega, 4], Assumptions -> $Assumptions];

expectZero["z0 formula", z0Toy - n0/d0];
expectZero["z2 formula", z2Toy - (n0 s2 - g2 d0)/d0^2];
expectZero["z4 formula", z4Toy - (n0 (s2^2 - d0) - s2 g2 d0)/d0^3];

subsDict = {d0 -> oA^2 oW^2 - r^2, s2 -> oA^2 + oW^2, n0 -> gA^2 oW^2 + 2 gA gW r + gW^2 oA^2, g2 -> gA^2 + gW^2};
sigmaSeries = Expand[Normal[Series[sigmaCons, {omega, 0, 4}]]];
z0 = FullSimplify[Coefficient[sigmaSeries, omega, 0], Assumptions -> $Assumptions];
z2 = FullSimplify[Coefficient[sigmaSeries, omega, 2], Assumptions -> $Assumptions];
z4 = FullSimplify[Coefficient[sigmaSeries, omega, 4], Assumptions -> $Assumptions];

Print["z0^(EM+mix) = ", fmt[z0]];
Print["z2^(EM+mix) = ", fmt[z2]];
Print["z4^(EM+mix) = ", fmt[z4]];
expectZero["Sigma z0", z0 - (n0/d0 /. subsDict)];
expectZero["Sigma z2", z2 - ((n0 s2 - g2 d0)/d0^2 /. subsDict)];
expectZero["Sigma z4", z4 - ((n0 (s2^2 - d0) - s2 g2 d0)/d0^3 /. subsDict)];

subbanner["III. Outgoing dressing of the mixed block"];

Clear[piOut, gammaPort];
$Assumptions =
  Element[{omega, oA, oW, r, gA, gW, piOut, gammaPort}, Reals] &&
  oA > 0 && oW > 0 && gammaPort > 0;

aKer = oA^2 - omega^2;
wKer = oW^2 - omega^2;
delta = FullSimplify[aKer wKer - r^2, Assumptions -> $Assumptions];
sigmaCons = FullSimplify[(gA^2 wKer + 2 gA gW r + gW^2 aKer)/delta, Assumptions -> $Assumptions];
sigmaFull = FullSimplify[(gA^2 (wKer - piOut) + 2 gA gW r + gW^2 aKer)/(aKer (wKer - piOut) - r^2), Assumptions -> $Assumptions];
nOmega = FullSimplify[D[sigmaFull, piOut] /. piOut -> 0, Assumptions -> $Assumptions];
n0 = FullSimplify[nOmega /. omega -> 0, Assumptions -> oA > 0 && oW > 0];
dCorr = FullSimplify[-I gammaPort omega^5 n0, Assumptions -> $Assumptions];

Print["Sigma_full(omega) = ", fmt[sigmaFull]];
Print["N(omega) = ", fmt[nOmega]];
Print["N(0) = ", fmt[n0]];
Print["delta D_wall^(odd) = ", fmt[dCorr]];
expectZero["N(omega) compact formula", nOmega - (aKer gW + r gA)^2/delta^2];
expectZero["N(0) positive-square form", n0 - (oA^2 gW + r gA)^2/(oA^2 oW^2 - r^2)^2];

subbanner["IV. Compact outgoing l=2 fingerprint"];

Clear[kWave, radius, cS, za];
$Assumptions =
  Element[{kWave, radius, cS, omega, za}, Reals] && kWave > 0 && radius > 0 && cS > 0;

h2a = SphericalHankelH1[2, za];
lambda2 = FullSimplify[(kWave D[h2a, za]/h2a) /. za -> kWave radius, Assumptions -> $Assumptions];
lambda2Series = Normal[Series[lambda2, {kWave, 0, 6}]];
y2 = Normal[Series[1/lambda2Series, {kWave, 0, 5}]] // FullSimplify;
y2Static = FullSimplify[y2 /. kWave -> 0, Assumptions -> radius > 0 && cS > 0];
y2Hat = Expand[y2/y2Static];
y2HatOmega = Expand[y2Hat /. kWave -> omega/cS];
gamma5Port = FullSimplify[Coefficient[y2HatOmega, omega, 5]/I, Assumptions -> radius > 0 && cS > 0];

Print["Lambda2(k) = ", fmt[Expand[lambda2Series]]];
Print["Y2_hat(omega) = ", fmt[y2HatOmega]];
Print["Gamma5_port = ", fmt[gamma5Port]];
expectZero[
  "Y2_hat minimal branch",
  y2HatOmega - (1 + radius^2 omega^2/(9 cS^2) + 4 radius^4 omega^4/(81 cS^4) + I radius^5 omega^5/(27 cS^5))
];
expectZero["Gamma5_port - a^5/(27 c_s^5)", gamma5Port - radius^5/(27 cS^5)];

subbanner["V. Derivative-coupled scalar mixed outlet starts at i omega^3"];

Clear[eta, gamma1];
$Assumptions =
  Element[{omega, oA, oW, r, eta, gamma1}, Reals] && oA > 0 && oW > 0 && eta > 0 && gamma1 > 0;

gA = 0;
gW = eta omega;
aKer = oA^2 - omega^2;
wKer = oW^2 - omega^2;
delta = FullSimplify[aKer wKer - r^2, Assumptions -> $Assumptions];
nOmega = FullSimplify[(aKer gW + r gA)^2/delta^2, Assumptions -> $Assumptions];
nSeries = Expand[Normal[Series[nOmega, {omega, 0, 2}]]];
pi0 = I gamma1 omega;
deltaD0 = Expand[pi0 nSeries];

Print["N_scalar(omega) = ", fmt[nSeries]];
Print["Pi0_out * N_scalar = ", fmt[deltaD0]];
expectZero["N_scalar leading term", nSeries - eta^2 oA^4 omega^2/(oA^2 oW^2 - r^2)^2];
expectZero["scalar odd order", deltaD0 - I gamma1 eta^2 oA^4 omega^3/(oA^2 oW^2 - r^2)^2];

Print[""];
Print["Stage 004 Mathematica audit passed."];

Exit[0];
