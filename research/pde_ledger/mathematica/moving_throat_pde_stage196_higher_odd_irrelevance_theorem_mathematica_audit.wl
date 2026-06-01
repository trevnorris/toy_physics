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
stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[stripCE[expr]]], Assumptions -> $Assumptions];
  res = stripCE[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectZeroList[name_String, exprs_List] := Module[{res},
  res = FullSimplify[Together[Expand[stripCE /@ exprs]], Assumptions -> $Assumptions];
  res = stripCE[res];
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[And @@ Thread[res === ConstantArray[0, Length[res]]]],
    pass[name],
    fail[name, res]
  ];
];

banner["STAGE 196 - HIGHER-ODD IRRELEVANCE THEOREM"];

Clear[w, pole, chi, tau7, z, l0, l2, l4, l5, l7];
Clear[x, scale, stretch, sig0, sig2, sig4, sig5];
Clear[radius, soundSpeed, lightSpeed, bigG];

$Assumptions =
  Element[
    {w, pole, chi, tau7, z, l0, l2, l4, l5, l7,
     x, scale, stretch, sig0, sig2, sig4, sig5,
     radius, soundSpeed, lightSpeed, bigG},
    Reals
  ] &&
  pole != 0 && l0 != 0 &&
  scale != 0 && stretch != 0 && 3*scale - sig0 != 0 &&
  scale*stretch^5 + 9*sig5 != 0 &&
  radius > 0 && soundSpeed > 0 && lightSpeed > 0 && bigG > 0;

(* I. Response-side tail: construct the tailful module and extract coefficients. *)
subbanner["I. Response-side higher-odd difference identity"];

sigmaCan = FullSimplify[9/(8*pole^5), Assumptions -> $Assumptions];
frontCore = w^2/pole^2 + I*chi*sigmaCan*w^5;
retResponse[tail_] := 3/4 + 1/(4*(1 - frontCore - tail));

retNoTail = retResponse[0];
retWithTail = retResponse[I*tau7*w^7];
retDifference = FullSimplify[retWithTail - retNoTail, Assumptions -> $Assumptions];
retDifferenceTarget = I*tau7*w^7/(4*(1 - frontCore)*(1 - frontCore - I*tau7*w^7));
retThrough5 = Expand[Normal[Series[retWithTail, {w, 0, 5}]]];
retDifferenceThrough7 = Expand[Normal[Series[retDifference, {w, 0, 7}]]];
retDifferenceLowCoeffs = Table[Coefficient[retDifferenceThrough7, w, k], {k, 0, 6}];
retSympyLowTarget = (
  1
  + w^2/(4*pole^2)
  + w^4/(4*pole^4)
  + I*chi*9*w^5/(32*pole^5)
);

Print["sigma_Q^can = ", fmt[sigmaCan]];
Print["retarded response through w^5 with explicit w^7 tail = ", fmt[retThrough5]];
Print["retarded tail difference through w^7 = ", fmt[retDifferenceThrough7]];

expectZero["exact response-side difference identity - SymPy target", retDifference - retDifferenceTarget];
expectZeroList["response-side tail absent through w^6", retDifferenceLowCoeffs];
expectZero[
  "Yhat_Q^(ret,>=7) through O(w^5) - SymPy one-pole form",
  retThrough5 - retSympyLowTarget
];
expectZero[
  "first response-side higher-odd correction - I tau7 w^7/4",
  retDifferenceThrough7 - I*tau7*w^7/4
];

(* II. Generic DtN-side tail: use denominator algebra and coefficient extraction. *)
subbanner["II. DtN-side higher-odd difference identity"];

dtnDenominator[tail_] := l0 + l2*z^2 + l4*z^4 + I*l5*z^5 + tail;
dtnResponse[tail_] := l0/dtnDenominator[tail];

dtnNoTail = dtnResponse[0];
dtnWithTail = dtnResponse[I*l7*z^7];
dtnDifference = FullSimplify[dtnWithTail - dtnNoTail, Assumptions -> $Assumptions];
dtnDifferenceTarget = -l0*(I*l7*z^7)/(dtnDenominator[0]*dtnDenominator[I*l7*z^7]);
dtnThrough5 = Expand[Normal[Series[dtnWithTail, {z, 0, 5}]]];
dtnThrough7 = Expand[Normal[Series[dtnWithTail, {z, 0, 7}]]];
dtnDifferenceLowCoeffs = Table[Coefficient[Expand[Normal[Series[dtnDifference, {z, 0, 7}]]], z, k], {k, 0, 6}];
dtnSympyLowTarget = (
  1
  - l2*z^2/l0
  + (l2^2/l0^2 - l4/l0)*z^4
  - I*l5*z^5/l0
);

Print["generic DtN response through z^7 = ", fmt[dtnThrough7]];

expectZero["exact DtN-side difference identity - SymPy target", dtnDifference - dtnDifferenceTarget];
expectZeroList["DtN tail absent through z^6", dtnDifferenceLowCoeffs];
expectZero["DtN normalized response through O(z^5) - SymPy target", dtnThrough5 - dtnSympyLowTarget];
expectZero[
  "L7 coefficient in z^7 DtN response - SymPy target",
  Coefficient[Coefficient[dtnThrough7, z, 7], l7] + I/l0
];

(* III. Stage 194 matching: derive the slots from the native outgoing l=2 operator. *)
subbanner["III. Canonical-even matching with an explicit z^7 DtN tail"];

lambdaOut = FunctionExpand[x*D[SphericalHankelH1[2, x], x]/SphericalHankelH1[2, x]];
lambdaWindow5 = Expand[Normal[Series[lambdaOut, {x, 0, 5}]]];
deformedDtn = Expand[
  scale*(lambdaWindow5 /. x -> stretch*z)
  + sig0
  + sig2*z^2
  + sig4*z^4
  + I*sig5*z^5
  + I*l7*z^7
];
normalDtn = Expand[Normal[Series[Coefficient[deformedDtn, z, 0]/deformedDtn, {z, 0, 7}]]];

evenRules = Solve[
  {
    Coefficient[normalDtn, z, 2] == 1/9,
    Coefficient[normalDtn, z, 4] == 4/81
  },
  {sig2, sig4},
  Reals
];
If[Length[evenRules] =!= 1, fail["unique canonical-even matching rule", evenRules]];
evenRule = First[evenRules];

sig2Target = -(3*scale*stretch^2 - 3*scale + sig0)/9;
sig4Target = -(3*scale*stretch^4 - 3*scale + sig0)/27;
matchedDtn = FullSimplify[normalDtn /. evenRule, Assumptions -> $Assumptions];
matchedThrough5 = Expand[Normal[Series[matchedDtn, {z, 0, 5}]]];
z5Slot = FullSimplify[Coefficient[matchedDtn, z, 5], Assumptions -> $Assumptions];
chiFromSeries = FullSimplify[z5Slot/(I/27), Assumptions -> $Assumptions];
chiSympyTarget = 3*(scale*stretch^5 + 9*sig5)/(3*scale - sig0);
matchedTargetThrough5 = (
  1
  + z^2/9
  + 4*z^4/81
  + I*chiSympyTarget*z^5/27
);

Print["outgoing l=2 DtN window through z^5 = ", fmt[lambdaWindow5]];
Print["Sigma_2 matching law = ", fmt[FullSimplify[sig2 /. evenRule, Assumptions -> $Assumptions]]];
Print["Sigma_4 matching law = ", fmt[FullSimplify[sig4 /. evenRule, Assumptions -> $Assumptions]]];
Print["matched normalized DtN response through z^5 = ", fmt[matchedThrough5]];
Print["chi_Q extracted from z^5 slot = ", fmt[chiFromSeries]];

expectZero["Sigma_2 matching law - SymPy target", (sig2 /. evenRule) - sig2Target];
expectZero["Sigma_4 matching law - SymPy target", (sig4 /. evenRule) - sig4Target];
expectZero["z^5 coefficient - deformation slot target", z5Slot - I*(scale*stretch^5/9 + sig5)/(3*scale - sig0)];
expectZero["canonical-even plus chi_Q matching through O(z^5)", matchedThrough5 - matchedTargetThrough5];
expectZero["chi_Q extractor - SymPy deformation formula", chiFromSeries - chiSympyTarget];
expectZero["d chi_Q^(series extractor) / dL7", D[chiFromSeries, l7]];
expectZeroList[
  "matched coefficients z^0..z^5 independent of L7",
  Table[D[Coefficient[matchedDtn, z, k], l7], {k, 0, 5}]
];

(* IV. Stage 195 source-map reduction inherits the same independence. *)
subbanner["IV. Source-map reduction is unchanged"];

p0Target = 54*bigG*soundSpeed^5/(5*radius^5*lightSpeed^5);
nNatural = FullSimplify[1/chiFromSeries, Assumptions -> $Assumptions];
deltaNatural = FullSimplify[p0Target*(nNatural - 1), Assumptions -> $Assumptions];
nSympyTarget = (3*scale - sig0)/(3*(scale*stretch^5 + 9*sig5));

Print["N_Q on the natural point-particle source-map branch = ", fmt[nNatural]];
Print["Delta_norm on the natural point-particle source-map branch = ", fmt[deltaNatural]];

expectZero["N_Q natural branch - SymPy target", nNatural - nSympyTarget];
expectZero["d N_Q / dL7", D[nNatural, l7]];
expectZero["d Delta_norm / dL7", D[deltaNatural, l7]];
expectZero[
  "Delta_norm - P0_target*(1/chi_Q - 1)",
  deltaNatural - p0Target*(1/chiFromSeries - 1)
];

banner["STAGE 196 MATHEMATICA LEDGER"];
Print["1. An explicit I tau7 w^7 retarded tail changes the grouped one-pole response first at w^7."];
Print["2. An explicit I L7 z^7 DtN tail changes the normalized DtN response first at z^7."];
Print["3. Native outgoing l=2 coefficient extraction gives the same canonical-even laws."];
Print["4. The extracted chi_Q and the natural source-map residual are independent of L7."];
Print["5. Therefore no higher odd datum beginning at order seven re-enters the 2.5PN finish line."];

Exit[0];
