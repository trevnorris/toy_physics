(* PathA-34 Gate 5 cross-ell unification, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
boolText[x_] := If[TrueQ[x], "true", "false"];
numText[x_] := StringReplace[
  ToString[N[x, 17], InputForm],
  {RegularExpression["`[0-9.]*"] -> "", "*^" -> "e"}
];
quoteText[x_] := Module[{s = If[StringQ[x], x, ToString[x, InputForm]]},
  "'" <> StringReplace[s, {"'" -> "''", "\n" -> " "}] <> "'"
];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_34_cross_l_unification.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
yamlOut = FileNameJoin[{scratchDir, "pathA_34_mathematica_results.yaml"}];
If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];

$Assumptions = Element[
    {z, omega, a, cs, M0, D1, N0, D0, N2, D2, N4, D4, K0c, Keta, TOmega,
     Z0ret, Z1ret, OmegaU, OmegaW, Rmix, gU, gW, etaNull, gain0, gain1},
    Reals
  ] && a > 0 && cs > 0 && D0 != 0 && K0c > 0 && Keta > 0 && TOmega > 0 &&
  Z0ret > 0 && Z1ret > 0 && OmegaU > 0 && OmegaW > 0 && gU > 0 && gW > 0 &&
  etaNull > 0;

serZ[expr_, order_] := Collect[Normal[Series[expr, {z, 0, order}]], z, FullSimplify];
serW[expr_, order_] := Collect[Normal[Series[expr, {omega, 0, order}]], omega, FullSimplify];

j0 = Sin[z]/z;
y0 = -Cos[z]/z;
j1 = Sin[z]/z^2 - Cos[z]/z;
y1 = -Cos[z]/z^2 - Sin[z]/z;
j2 = (3/z^3 - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2 = (1/z - 3/z^3) Cos[z] - 3 Sin[z]/z^2;

branchData[ell_, h_, kind_] := Module[
  {lam, yout, yser, lamser, pow, rad},
  lam = FullSimplify[z D[h, z]/h, $Assumptions];
  yout = FullSimplify[-(ell + 1)/lam, $Assumptions];
  lamser = serZ[lam, 8];
  yser = serZ[yout, 8];
  pow = 2 ell + 1;
  rad = FullSimplify[Coefficient[yser, z, pow]/I, $Assumptions];
  <|
    "ell" -> ell,
    "kind" -> kind,
    "normalizationFactor" -> -(ell + 1),
    "lambdaSeries" -> lamser,
    "YSeries" -> yser,
    "static" -> FullSimplify[Coefficient[yser, z, 0], $Assumptions],
    "radiativePower" -> pow,
    "radiativeCoeff" -> rad,
    "rawOrder" -> If[TrueQ[rad == 0], "null", pow],
    "u2z" -> FullSimplify[Coefficient[yser, z, 2], $Assumptions],
    "u4z" -> FullSimplify[Coefficient[yser, z, 4], $Assumptions],
    "v5z" -> FullSimplify[Coefficient[yser, z, 5]/I, $Assumptions]
  |>
];

out0 = branchData[0, j0 + I y0, "outgoing_hankel1"];
out1 = branchData[1, j1 + I y1, "outgoing_hankel1"];
out2 = branchData[2, j2 + I y2, "outgoing_hankel1"];
in0 = branchData[0, j0 - I y0, "incoming_hankel2"];
in1 = branchData[1, j1 - I y1, "incoming_hankel2"];
in2 = branchData[2, j2 - I y2, "incoming_hankel2"];

fingerprintOk = And[
  TrueQ[out0["normalizationFactor"] == -1],
  TrueQ[out1["normalizationFactor"] == -2],
  TrueQ[out2["normalizationFactor"] == -3],
  TrueQ[FullSimplify[out0["static"] == 1, $Assumptions]],
  TrueQ[FullSimplify[out1["static"] == 1, $Assumptions]],
  TrueQ[FullSimplify[out2["static"] == 1, $Assumptions]],
  TrueQ[FullSimplify[out0["radiativeCoeff"] == 1, $Assumptions]],
  TrueQ[FullSimplify[out1["radiativeCoeff"] == 1/2, $Assumptions]],
  TrueQ[FullSimplify[out2["radiativeCoeff"] == 1/27, $Assumptions]],
  TrueQ[FullSimplify[in0["radiativeCoeff"] == -out0["radiativeCoeff"], $Assumptions]],
  TrueQ[FullSimplify[in1["radiativeCoeff"] == -out1["radiativeCoeff"], $Assumptions]],
  TrueQ[FullSimplify[in2["radiativeCoeff"] == -out2["radiativeCoeff"], $Assumptions]]
];

Pport = OmegaU^2 gW + Rmix gU;
DeltaPort = OmegaU^2 OmegaW^2 - Rmix^2;
N0port = FullSimplify[Pport^2/DeltaPort^2, $Assumptions];
P0port = FullSimplify[N0port/D0, $Assumptions];

Nomega = N0port + N2 omega^2 + N4 omega^4;
Dcons = D0 + D2 omega^2 + D4 omega^4;
prefObj = D0 Nomega/Dcons^2;
prefSeries = serW[prefObj, 6];
p0 = FullSimplify[Coefficient[prefSeries, omega, 0], $Assumptions];
p2 = FullSimplify[Coefficient[prefSeries, omega, 2], $Assumptions];
p4 = FullSimplify[Coefficient[prefSeries, omega, 4], $Assumptions];
expectedP0 = N0port/D0;
expectedP2 = (D0 N2 - 2 D2 N0port)/D0^2;
expectedP4 = (D0^2 N4 - 2 D0 (D2 N2 + D4 N0port) + 3 D2^2 N0port)/D0^3;
resP0 = FullSimplify[p0 - expectedP0, $Assumptions];
resP2 = FullSimplify[p2 - expectedP2, $Assumptions];
resP4 = FullSimplify[p4 - expectedP4, $Assumptions];
chiQ = FullSimplify[out2["v5z"]/(1/27), $Assumptions];
gate4Ok = And[
  TrueQ[FullSimplify[out2["u2z"] == 1/9, $Assumptions]],
  TrueQ[FullSimplify[out2["u4z"] == 4/81, $Assumptions]],
  TrueQ[FullSimplify[out2["v5z"] == 1/27, $Assumptions]],
  TrueQ[FullSimplify[chiQ == 1, $Assumptions]],
  TrueQ[resP0 == 0], TrueQ[resP2 == 0], TrueQ[resP4 == 0]
];

K1dc = Keta + 2 TOmega;
T0dc = FullSimplify[K0c/(K0c + Z0ret), $Assumptions];
T1dc = FullSimplify[K1dc/(K1dc + Z1ret), $Assumptions];
eps0 = FullSimplify[Z0ret/K0c, $Assumptions];
eps1 = FullSimplify[Z1ret/K1dc, $Assumptions];
A0lead = FullSimplify[I out0["radiativeCoeff"] (a omega/cs) M0 (1 - T0dc), $Assumptions];
A1lead = FullSimplify[I out1["radiativeCoeff"] (a omega/cs)^3 D1 (1 - T1dc), $Assumptions];
expectedA0 = FullSimplify[I a omega M0 eps0/(cs (1 + eps0)), $Assumptions];
expectedA1 = FullSimplify[I a^3 omega^3 D1 eps1/(2 cs^3 (1 + eps1)), $Assumptions];
resA0 = FullSimplify[A0lead - expectedA0, $Assumptions];
resA1 = FullSimplify[A1lead - expectedA1, $Assumptions];
residualOk = TrueQ[resA0 == 0] && TrueQ[resA1 == 0];

rankDofs = {OmegaU, OmegaW, Rmix, gU, gW, D0, K0c, Keta, TOmega, Z0ret, Z1ret};
rankRow[expr_, dofs_] := FullSimplify[D[expr, #], $Assumptions] & /@ dofs;
rankAuditFor[extraRows_, t0expr_, t1expr_, dofs_] := Module[
  {rows, r0, nullity, augRank, movingNullity},
  rows = rankRow[#, dofs] & /@ Join[{P0port, K0c, Keta + 2 TOmega}, extraRows];
  r0 = MatrixRank[rows];
  nullity = Length[dofs] - r0;
  augRank = MatrixRank[Join[rows, {rankRow[t0expr, dofs], rankRow[t1expr, dofs]}]];
  movingNullity = augRank - r0;
  <|
    "rows" -> rows,
    "rank" -> r0,
    "nullity" -> nullity,
    "returnMovingNullity" -> movingNullity,
    "movesReturn" -> TrueQ[movingNullity > 0]
  |>
];
rankBaseline = rankAuditFor[{}, T0dc, T1dc, rankDofs];
nativeNullMovesReturn = rankBaseline["movesReturn"];
selectorT0 = FullSimplify[K0c/(K0c + K0c), $Assumptions];
selectorT1 = FullSimplify[K1dc/(K1dc + K1dc), $Assumptions];
rankSelector = rankAuditFor[{Z0ret - K0c, Z1ret - (Keta + 2 TOmega)}, selectorT0, selectorT1, rankDofs];

positiveBoundedQ[t_, eps_] := TrueQ[FullSimplify[eps > 0 && t > 0 && t < 1, $Assumptions]];
overcancelQ[epsA_, epsB_] := TrueQ[FullSimplify[epsA == 0 && epsB == 0, $Assumptions]];
noConsistentQ[epsA_, epsB_, tA_, tB_] := TrueQ[
  FullSimplify[epsA == -2 && epsB == -2 && ! (tA > 0 && tA < 1 && tB > 0 && tB < 1), $Assumptions]
];
gateVerdictFor[z0_, z1_, extraRows_, dimOk_, tautological_] := Module[
  {k1, t0, t1, e0, e1, rankA, positive, over, noConsistent, epsMismatch},
  k1 = Keta + 2 TOmega;
  t0 = FullSimplify[K0c/(K0c + z0), $Assumptions];
  t1 = FullSimplify[k1/(k1 + z1), $Assumptions];
  e0 = FullSimplify[z0/K0c, $Assumptions];
  e1 = FullSimplify[z1/k1, $Assumptions];
  rankA = rankAuditFor[extraRows, t0, t1, rankDofs];
  positive = positiveBoundedQ[t0, e0] && positiveBoundedQ[t1, e1];
  over = overcancelQ[e0, e1];
  noConsistent = noConsistentQ[e0, e1, t0, t1];
  epsMismatch = ! positive;
  Which[
    tautological, "FAIL_TAUTOLOGICAL",
    ! dimOk, "FAIL_DIMENSIONAL",
    noConsistent, "FAIL_NO_CONSISTENT_RETURN",
    over, "FAIL_OVERCANCEL",
    epsMismatch, "FAIL_EPSILON_MISMATCH",
    TrueQ[rankA["movesReturn"]], "FAIL_UNDERDETERMINED_NOT_PREDICTIVE",
    True, "CROSS_L_RESIDUAL_PREDICTION"
  ]
];

zeroDim = {0, 0, 0};
dimScale[x_, q_] := q x;
dimOf[expr_, dims_] := Module[{args, ds, base, pow},
  Which[
    TrueQ[expr == 0] || NumericQ[expr], zeroDim,
    AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr],
    AtomQ[expr], fail["missing dimension for " <> ToString[Unevaluated[expr], InputForm]],
    Head[expr] === Times, Total[dimOf[#, dims] & /@ (List @@ expr)],
    Head[expr] === Power,
      base = expr[[1]];
      pow = expr[[2]];
      If[! NumericQ[pow], fail["non-numeric dimension exponent"]];
      dimScale[dimOf[base, dims], pow],
    Head[expr] === Plus,
      args = List @@ expr;
      ds = dimOf[#, dims] & /@ Select[args, ! TrueQ[# == 0] &];
      If[Length[ds] == 0, zeroDim,
        If[Length[DeleteDuplicates[ds]] != 1, fail["dimension mismatch in sum"]];
        First[ds]
      ],
    True, fail["unsupported dimension expression " <> ToString[expr, InputForm]]
  ]
];

baseDims = <|
  a -> {1, 0, 0},
  cs -> {1, 0, -1},
  omega -> {0, 0, -1},
  M0 -> {0, 1, -1},
  D1 -> {1, 1, -1},
  R0 -> {0, 1, -1},
  R1 -> {1, 1, -1},
  N0 -> {-1, 1, 0},
  D0 -> {-1, 1, -2},
  K0c -> {0, 1, -2},
  Keta -> {0, 1, -2},
  TOmega -> {0, 1, -2},
  Z0ret -> {0, 1, -2},
  Z1ret -> {0, 1, -2},
  OmegaU -> {0, 0, -1},
  OmegaW -> {0, 0, -1},
  Rmix -> {0, 0, -2},
  gU -> {-1/2, 1/2, -2},
  gW -> {-1/2, 1/2, -2},
  etaNull -> zeroDim,
  gain0 -> zeroDim,
  gain1 -> zeroDim,
  qfree -> zeroDim
|>;

expectedDims = <|
  "M0" -> {0, 1, -1},
  "D1" -> {1, 1, -1},
  "R0" -> {0, 1, -1},
  "R1" -> {1, 1, -1},
  "A0" -> {0, 1, -1},
  "A1" -> {1, 1, -1},
  "T0" -> zeroDim,
  "T1" -> zeroDim,
  "epsilon0_eff" -> zeroDim,
  "epsilon1_eff" -> zeroDim,
  "N0" -> {-1, 1, 0},
  "D0" -> {-1, 1, -2},
  "P0_physical" -> zeroDim
|>;

dimEntries = {
  {"M0", M0}, {"D1", D1}, {"R0", R0}, {"R1", R1},
  {"A0", A0lead}, {"A1", A1lead}, {"T0", T0dc}, {"T1", T1dc},
  {"epsilon0_eff", eps0}, {"epsilon1_eff", eps1}, {"N0", N0port}, {"D0", D0},
  {"P0_physical", (cs/a)^2 P0port}
};

dimCheck[dims_] := Module[{computed},
  computed = Table[
    With[{name = dimEntries[[i, 1]], expr = dimEntries[[i, 2]]},
      If[TrueQ[expr == 0], expectedDims[name], dimOf[expr, dims]] == expectedDims[name]
    ],
    {i, Length[dimEntries]}
  ];
  And @@ computed
];
dimensionalOk = TrueQ[dimCheck[baseDims]];
corruptSourcedDims = Join[KeyDrop[baseDims, M0], <|M0 -> {1, 1, -1}|>];
corruptSourcedOk = TrueQ[dimCheck[corruptSourcedDims]];
corruptFreeDims = Join[KeyDrop[baseDims, qfree], <|qfree -> {7, 0, 0}|>];
corruptFreeOk = TrueQ[dimCheck[corruptFreeDims]];

baselineVerdict = gateVerdictFor[Z0ret, Z1ret, {}, dimensionalOk, False];
selectorVerdict = gateVerdictFor[
  K0c,
  Keta + 2 TOmega,
  {Z0ret - K0c, Z1ret - (Keta + 2 TOmega)},
  dimensionalOk,
  False
];
probe3cWith = gateVerdictFor[-Z0ret, -Z1ret, {}, dimensionalOk, False];
probe3dWith = gateVerdictFor[0, 0, {}, dimensionalOk, False];
probe3hWith = gateVerdictFor[-2 K0c, -2 (Keta + 2 TOmega), {}, dimensionalOk, False];
probe3fWith = gateVerdictFor[Z0ret, Z1ret, {}, corruptSourcedOk, False];

sampleRules = {
  a -> 3.0, cs -> 2.0, omega -> 0.07, M0 -> 5.0, D1 -> 11.0,
  K0c -> 13.0, Keta -> 17.0, TOmega -> 19.0, Z0ret -> 23.0, Z1ret -> 29.0,
  N0 -> 31.0, D0 -> 37.0, N2 -> 41.0, D2 -> 43.0, N4 -> 47.0, D4 -> 53.0,
  OmegaU -> 2.0, OmegaW -> 3.0, Rmix -> 1.0, gU -> 5.0, gW -> 7.0,
  etaNull -> 0.25, gain0 -> 1.3, gain1 -> 0.7
};
evalSample[expr_] := N[expr /. sampleRules, 17];

headlineOk = fingerprintOk && gate4Ok && residualOk && dimensionalOk &&
  TrueQ[rankBaseline["nullity"] == 8] && TrueQ[nativeNullMovesReturn] &&
  TrueQ[rankBaseline["returnMovingNullity"] == 2] &&
  TrueQ[baselineVerdict == "FAIL_UNDERDETERMINED_NOT_PREDICTIVE"] &&
  TrueQ[selectorVerdict == "CROSS_L_RESIDUAL_PREDICTION"] &&
  TrueQ[probe3cWith == "FAIL_EPSILON_MISMATCH"] &&
  TrueQ[probe3dWith == "FAIL_OVERCANCEL"] &&
  TrueQ[probe3hWith == "FAIL_NO_CONSISTENT_RETURN"] &&
  TrueQ[probe3fWith == "FAIL_DIMENSIONAL"] &&
  ! corruptSourcedOk && corruptFreeOk;

lines = {
  "schema: pathA_34_mathematica_scratch/v1",
  "engine: mathematica",
  "fingerprints:",
  "  ell0:",
  "    Y_series_z: " <> quoteText[out0["YSeries"]],
  "    normalization_factor: " <> ToString[out0["normalizationFactor"]],
  "    radiative_coeff_z: " <> quoteText[out0["radiativeCoeff"]],
  "    raw_order: " <> ToString[out0["rawOrder"]],
  "    incoming_radiative_coeff_z: " <> quoteText[in0["radiativeCoeff"]],
  "  ell1:",
  "    Y_series_z: " <> quoteText[out1["YSeries"]],
  "    normalization_factor: " <> ToString[out1["normalizationFactor"]],
  "    radiative_coeff_z: " <> quoteText[out1["radiativeCoeff"]],
  "    raw_order: " <> ToString[out1["rawOrder"]],
  "    incoming_radiative_coeff_z: " <> quoteText[in1["radiativeCoeff"]],
  "  ell2:",
  "    Y_series_z: " <> quoteText[out2["YSeries"]],
  "    normalization_factor: " <> ToString[out2["normalizationFactor"]],
  "    radiative_coeff_z: " <> quoteText[out2["radiativeCoeff"]],
  "    raw_order: " <> ToString[out2["rawOrder"]],
  "    incoming_radiative_coeff_z: " <> quoteText[in2["radiativeCoeff"]],
  "residuals:",
  "  T0_dc: " <> quoteText[T0dc],
  "  T1_dc: " <> quoteText[T1dc],
  "  epsilon0_eff: " <> quoteText[eps0],
  "  epsilon1_eff: " <> quoteText[eps1],
  "  A0_leading: " <> quoteText[A0lead],
  "  A1_leading: " <> quoteText[A1lead],
  "  A0_expected: " <> quoteText[expectedA0],
  "  A1_expected: " <> quoteText[expectedA1],
  "  A0_residual_to_expected: " <> quoteText[resA0],
  "  A1_residual_to_expected: " <> quoteText[resA1],
  "  A0_sample:",
  "    re: " <> numText[Re[evalSample[A0lead]]],
  "    im: " <> numText[Im[evalSample[A0lead]]],
  "  A1_sample:",
  "    re: " <> numText[Re[evalSample[A1lead]]],
  "    im: " <> numText[Im[evalSample[A1lead]]],
  "  T0_sample: " <> numText[evalSample[T0dc]],
  "  T1_sample: " <> numText[evalSample[T1dc]],
  "  epsilon0_sample: " <> numText[evalSample[eps0]],
  "  epsilon1_sample: " <> numText[evalSample[eps1]],
  "gate4_non_regression:",
  "  ok: " <> boolText[gate4Ok],
  "  chi_Q: " <> quoteText[chiQ],
  "  P0_residual: " <> quoteText[resP0],
  "  P2_residual: " <> quoteText[resP2],
  "  P4_residual: " <> quoteText[resP4],
  "  P0_sample: " <> numText[evalSample[p0]],
  "  P2_sample: " <> numText[evalSample[p2]],
  "  P4_sample: " <> numText[evalSample[p4]],
  "rank:",
  "  native_nullspace_dimension: " <> ToString[rankBaseline["nullity"]],
  "  native_null_moves_return: " <> boolText[nativeNullMovesReturn],
  "  return_moving_nullity: " <> ToString[rankBaseline["returnMovingNullity"]],
  "  determined_prediction: " <> boolText[! nativeNullMovesReturn],
  "  baseline_verdict: " <> quoteText[baselineVerdict],
  "  selector_control_verdict: " <> quoteText[selectorVerdict],
  "dimension:",
  "  dimensional_ok: " <> boolText[dimensionalOk],
  "  corrupt_sourced_verdict: " <> quoteText[If[corruptSourcedOk, "NO_FAIL", "FAIL_DIMENSIONAL"]],
  "  corrupt_free_carrier_verdict: " <> quoteText[If[corruptFreeOk, "NO_FAIL", "FAIL_DIMENSIONAL"]],
  "headline_booleans:",
  "  fingerprints_ok: " <> boolText[fingerprintOk],
  "  residual_form_ok: " <> boolText[residualOk],
  "  gate4_ok: " <> boolText[gate4Ok],
  "  dimension_ok: " <> boolText[dimensionalOk],
  "verdict_probes:",
  "  rank_baseline: " <> quoteText[baselineVerdict],
  "  rank_selector_equation: " <> quoteText[selectorVerdict],
  "  3c_wrong_sign:",
  "    with_mutation: " <> quoteText[probe3cWith],
  "    without_mutation: " <> quoteText[baselineVerdict],
  "  3d_perfect_return:",
  "    with_mutation: " <> quoteText[probe3dWith],
  "    without_mutation: " <> quoteText[baselineVerdict],
  "  3h_no_consistent_return:",
  "    with_mutation: " <> quoteText[probe3hWith],
  "    without_mutation: " <> quoteText[baselineVerdict],
  "  3f_corrupt_dimension:",
  "    with_mutation: " <> quoteText[probe3fWith],
  "    without_mutation: " <> quoteText[baselineVerdict]
};

Export[yamlOut, StringRiffle[lines, "\n"] <> "\n", "Text"];
If[! TrueQ[headlineOk], fail["headline checks failed"]];
Print["Wrote Mathematica scratch: ", yamlOut];
Exit[0];
