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
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

cleanScalar[expr_, assumptions_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Together[res], Assumptions -> assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> assumptions]
];

expectZeroUnder[name_String, expr_, assumptions_] := Module[{res},
  res = cleanScalar[expr, assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectZero[name_String, expr_] := expectZeroUnder[name, expr, $Assumptions];

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

totalDegreeRS[poly_] := Module[{rules, powers},
  rules = CoefficientRules[Expand[poly], {r, s}];
  powers = rules[[All, 1]];
  Max[Total /@ powers]
];

banner["STAGE 211 -- FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION"];

Clear[
  r, s, ki, kj, kk, H0, A, B, Ccoef, Dcoef, Ecoef, Fcoef,
  u, k, ud, ux
];

den = 1 + r^2 + s^2;
linearK = ki + kj r + kk s;
Delta = A + B r + Ccoef s + Dcoef r^2 + Ecoef r s + Fcoef s^2;
sqrtDelta = Sqrt[Delta];
Phi = linearK/Sqrt[den] + sqrtDelta/Sqrt[den];
tau = 2 H0/Phi;

$Assumptions = (
  Element[
    {r, s, ki, kj, kk, H0, A, B, Ccoef, Dcoef, Ecoef, Fcoef, u, k, ud, ux},
    Reals
  ]
  && r >= 0 && s >= 0
  && ki > 0 && kj > 0 && kk > 0 && H0 > 0 && k > 0
  && Delta > 0
);

subbanner["M1. Exact stationary numerator law"];

stationaryDenominator = 2 den^(3/2) sqrtDelta;
directNumerator[var_] := FullSimplify[
  Numerator[Together[D[Phi, var]]],
  Assumptions -> $Assumptions
];

Mr = den kj - r linearK;
Ms = den kk - s linearK;
Lr = den D[Delta, r] - 2 r Delta;
Ls = den D[Delta, s] - 2 s Delta;
Nr = 2 Mr sqrtDelta + Lr;
Ns = 2 Ms sqrtDelta + Ls;

Print["M1 direct numerator from D[Phi, r] = ", fmt[directNumerator[r]]];
Print["M1 direct numerator from D[Phi, s] = ", fmt[directNumerator[s]]];
expectZero["M1 D[Phi,r] numerator minus paper N_r", directNumerator[r] - Nr];
expectZero["M1 D[Phi,s] numerator minus paper N_s", directNumerator[s] - Ns];
expectZero["M1 derivative law r", D[Phi, r] - Nr/stationaryDenominator];
expectZero["M1 derivative law s", D[Phi, s] - Ns/stationaryDenominator];

subbanner["M2. Quartic cross-consistency"];

Ccross = Expand[Ms Lr - Mr Ls];
crossDegree = totalDegreeRS[Ccross];

expectZero["M2 square-root-free cross identity", Ms Nr - Mr Ns - Ccross];
Print["M2 total degree C_cross = ", crossDegree];
expectTrue["M2 C_cross has total degree 4", crossDegree == 4];

subbanner["M3. Sextic square eliminants"];

Sr = Expand[Lr^2 - 4 Mr^2 Delta];
Ss = Expand[Ls^2 - 4 Ms^2 Delta];
srDegree = totalDegreeRS[Sr];
ssDegree = totalDegreeRS[Ss];

expectZero["M3 square eliminant identity r", Nr (Nr - 4 Mr sqrtDelta) - Sr];
expectZero["M3 square eliminant identity s", Ns (Ns - 4 Ms sqrtDelta) - Ss];
Print["M3 total degree S_r = ", srDegree];
Print["M3 total degree S_s = ", ssDegree];
expectTrue["M3 S_r has total degree 6", srDegree == 6];
expectTrue["M3 S_s has total degree 6", ssDegree == 6];

subbanner["M4. Bezout bound"];

bezoutBound = crossDegree srDegree;
Print["M4 computed Bezout product = ", bezoutBound];
expectTrue["M4 Bezout product equals 24", bezoutBound == 24];

subbanner["M5. Diagonal-isotropic reduction"];

isoSubstitution = {
  A -> ki^2 - 2 H0 u,
  B -> 2 ki kj,
  Ccoef -> 2 ki kk,
  Dcoef -> kj^2 - 2 H0 u,
  Ecoef -> 2 kj kk,
  Fcoef -> kk^2 - 2 H0 u
};
DeltaIso = FullSimplify[Delta /. isoSubstitution, Assumptions -> $Assumptions];
DeltaIsoExpected = linearK^2 - 2 H0 u den;
krs = linearK/Sqrt[den];
tauIso = tau /. isoSubstitution;
tauIsoExpected = 2 H0/(krs + Sqrt[krs^2 - 2 H0 u]);
isoAssumptions = $Assumptions && DeltaIsoExpected > 0 && krs^2 - 2 H0 u > 0;

expectZero["M5 Delta_iso reduction", DeltaIso - DeltaIsoExpected];
expectZeroUnder["M5 tau_iso reduction", tauIso - tauIsoExpected, isoAssumptions];

subbanner["M6. Full-symmetry equal-mix stationarity"];

symSubstitution = {
  ki -> k,
  kj -> k,
  kk -> k,
  A -> k^2 - 2 H0 ud,
  B -> 2 k^2 - 4 H0 ux,
  Ccoef -> 2 k^2 - 4 H0 ux,
  Dcoef -> k^2 - 2 H0 ud,
  Ecoef -> 2 k^2 - 4 H0 ux,
  Fcoef -> k^2 - 2 H0 ud
};
NrSymEqualMix = Nr /. symSubstitution /. {r -> 1, s -> 1};
NsSymEqualMix = Ns /. symSubstitution /. {r -> 1, s -> 1};

expectZero["M6 symmetric N_r(1,1)", NrSymEqualMix];
expectZero["M6 symmetric N_s(1,1)", NsSymEqualMix];

Print[""];
Print["All Stage 211 identities verified."];
Exit[0];
