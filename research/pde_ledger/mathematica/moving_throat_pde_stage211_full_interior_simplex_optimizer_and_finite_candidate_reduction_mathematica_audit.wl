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
primitivePartRS[poly_] := Module[{parts},
  parts = FactorTermsList[Expand[poly], {r, s}];
  Last[parts]
];

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

stationaryNumerator[var_] := Numerator[Together[D[Phi, var]]];
stationaryDenominator[var_] := Denominator[Together[D[Phi, var]]];

numR = stationaryNumerator[r];
numS = stationaryNumerator[s];

Print["M1 direct numerator from D[Phi, r] = ", fmt[numR]];
Print["M1 direct numerator from D[Phi, s] = ", fmt[numS]];
expectTrue["M1 D[Phi,r] numerator is nonzero", ! TrueQ[Expand[numR] === 0]];
expectTrue["M1 D[Phi,s] numerator is nonzero", ! TrueQ[Expand[numS] === 0]];
expectZero["M1 reconstruct D[Phi,r] from derived numerator", D[Phi, r] - numR/stationaryDenominator[r]];
expectZero["M1 reconstruct D[Phi,s] from derived numerator", D[Phi, s] - numS/stationaryDenominator[s]];

subbanner["M2. Quartic cross-consistency"];

numRq = Expand[numR /. Sqrt[Delta] -> q];
numSq = Expand[numS /. Sqrt[Delta] -> q];
crossDerived = primitivePartRS[Resultant[numRq, numSq, q]];
crossDegree = totalDegreeRS[crossDerived];

expectTrue["M2 numerator r is polynomial in q", PolynomialQ[numRq, q]];
expectTrue["M2 numerator s is polynomial in q", PolynomialQ[numSq, q]];
expectTrue["M2 resultant removed q", FreeQ[crossDerived, q]];
Print["M2 derived cross polynomial = ", fmt[crossDerived]];
Print["M2 total degree derived cross = ", crossDegree];
expectTrue["M2 derived cross has total degree 4", crossDegree == 4];

(* SymPy comparison target only. *)
mRTarget = Expand[den kj - r linearK];
mSTarget = Expand[den kk - s linearK];
lRTarget = Expand[den D[Delta, r] - 2 r Delta];
lSTarget = Expand[den D[Delta, s] - 2 s Delta];
crossTarget = Expand[mSTarget lRTarget - mRTarget lSTarget];
crossRatio = cleanScalar[crossDerived/crossTarget, $Assumptions];

Print["M2 derived cross / SymPy target = ", fmt[crossRatio]];
expectTrue[
  "M2 derived cross ratio is nonzero constant",
  FreeQ[crossRatio, r] && FreeQ[crossRatio, s] && TrueQ[crossRatio != 0]
];
expectZero["M2 derived cross minus scaled SymPy target", crossDerived - crossRatio crossTarget];

subbanner["M3. Sextic square eliminants"];

srDerived = primitivePartRS[Resultant[numRq, q^2 - Delta, q]];
ssDerived = primitivePartRS[Resultant[numSq, q^2 - Delta, q]];
srDegree = totalDegreeRS[srDerived];
ssDegree = totalDegreeRS[ssDerived];

expectTrue["M3 r square eliminant removed q", FreeQ[srDerived, q]];
expectTrue["M3 s square eliminant removed q", FreeQ[ssDerived, q]];
Print["M3 derived S_r polynomial = ", fmt[srDerived]];
Print["M3 derived S_s polynomial = ", fmt[ssDerived]];
Print["M3 total degree derived S_r = ", srDegree];
Print["M3 total degree derived S_s = ", ssDegree];
expectTrue["M3 derived S_r has total degree 6", srDegree == 6];
expectTrue["M3 derived S_s has total degree 6", ssDegree == 6];

squareTargetR = Expand[lRTarget^2 - 4 mRTarget^2 Delta];
squareTargetS = Expand[lSTarget^2 - 4 mSTarget^2 Delta];
srRatio = cleanScalar[srDerived/squareTargetR, $Assumptions];
ssRatio = cleanScalar[ssDerived/squareTargetS, $Assumptions];

Print["M3 derived S_r / SymPy target = ", fmt[srRatio]];
Print["M3 derived S_s / SymPy target = ", fmt[ssRatio]];
expectTrue[
  "M3 derived S_r ratio is nonzero constant",
  FreeQ[srRatio, r] && FreeQ[srRatio, s] && TrueQ[srRatio != 0]
];
expectTrue[
  "M3 derived S_s ratio is nonzero constant",
  FreeQ[ssRatio, r] && FreeQ[ssRatio, s] && TrueQ[ssRatio != 0]
];
expectZero["M3 derived S_r minus scaled SymPy target", srDerived - srRatio squareTargetR];
expectZero["M3 derived S_s minus scaled SymPy target", ssDerived - ssRatio squareTargetS];

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
tauIsoExpected = 2 H0/(krs + Sqrt[krs^2 - 2 H0 u]);
isoAssumptions = $Assumptions && DeltaIsoExpected > 0 && krs^2 - 2 H0 u > 0;
tauIso = FullSimplify[
  PowerExpand[tau /. isoSubstitution],
  Assumptions -> isoAssumptions
];

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
numRSymEqualMix = FullSimplify[numR /. symSubstitution /. {r -> 1, s -> 1}];
numSSymEqualMix = FullSimplify[numS /. symSubstitution /. {r -> 1, s -> 1}];

expectZero["M6 symmetric derived numerator r(1,1)", numRSymEqualMix];
expectZero["M6 symmetric derived numerator s(1,1)", numSSymEqualMix];

Print[""];
Print["All Stage 211 identities verified."];
Exit[0];
