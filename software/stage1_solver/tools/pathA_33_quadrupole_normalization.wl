(* PathA-33 Gate 4 quadrupole normalization, Mathematica engine. *)

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
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_33_quadrupole_normalization.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
yamlOut = FileNameJoin[{scratchDir, "pathA_33_mathematica_results.yaml"}];
If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];

$Assumptions = Element[{z, omega, a, cs, c, G, D0, D2, D4, N0, N2, N4}, Reals] &&
  a > 0 && cs > 0 && c > 0 && D0 != 0;

serZ[expr_, order_] := Collect[Normal[Series[expr, {z, 0, order}]], z, FullSimplify];
serW[expr_, order_] := Collect[Normal[Series[expr, {omega, 0, order}]], omega, FullSimplify];

j2 = (3/z^3 - 1/z) Sin[z] - 3 Cos[z]/z^2;
y2 = (1/z - 3/z^3) Cos[z] - 3 Sin[z]/z^2;

branchData[name_, h_, source_] := Module[
  {lam, yhat, hser, lamser, yser, u2z, u4z, v5z, static},
  lam = FullSimplify[z D[h, z]/h, $Assumptions];
  yhat = FullSimplify[-3/lam, $Assumptions];
  hser = serZ[h, 7];
  lamser = serZ[lam, 6];
  yser = serZ[yhat, 6];
  u2z = FullSimplify[Coefficient[yser, z, 2], $Assumptions];
  u4z = FullSimplify[Coefficient[yser, z, 4], $Assumptions];
  v5z = FullSimplify[Coefficient[yser, z, 5]/I, $Assumptions];
  static = FullSimplify[Coefficient[yser, z, 0], $Assumptions];
  <|
    "name" -> name,
    "source" -> source,
    "hSeries" -> hser,
    "lambda" -> lam,
    "lambdaSeries" -> lamser,
    "yhat" -> yhat,
    "yhatSeries" -> yser,
    "static" -> static,
    "u2z" -> u2z,
    "u4z" -> u4z,
    "v5z" -> v5z,
    "u2" -> FullSimplify[u2z a^2/cs^2, $Assumptions],
    "u4" -> FullSimplify[u4z a^4/cs^4, $Assumptions],
    "v5" -> FullSimplify[v5z a^5/cs^5, $Assumptions]
  |>
];

out = branchData["outgoing_hankel1", j2 + I y2, "dtn_hankel1"];
incoming = branchData["incoming_hankel2", j2 - I y2, "dtn_hankel2"];
standing = branchData["standing_j2", j2, "regular_standing_j2"];

u2Match = TrueQ[FullSimplify[out["u2z"] == 1/9, $Assumptions]];
u4Match = TrueQ[FullSimplify[out["u4z"] == 4/81, $Assumptions]];
v5Match = TrueQ[FullSimplify[out["v5z"] == 1/27, $Assumptions]];
chiQ = FullSimplify[out["v5"]/(a^5/(27 cs^5)), $Assumptions];
chiQIncoming = FullSimplify[incoming["v5"]/(a^5/(27 cs^5)), $Assumptions];
chiMatch = TrueQ[FullSimplify[chiQ == 1, $Assumptions]];

Nomega = N0 + N2 omega^2 + N4 omega^4;
Dcons = D0 + D2 omega^2 + D4 omega^4;
prefObj = D0 Nomega/Dcons^2;
plainObj = Nomega/Dcons;
prefSeries = serW[prefObj, 6];
plainSeries = serW[plainObj, 6];
p0 = FullSimplify[Coefficient[prefSeries, omega, 0], $Assumptions];
p2 = FullSimplify[Coefficient[prefSeries, omega, 2], $Assumptions];
p4 = FullSimplify[Coefficient[prefSeries, omega, 4], $Assumptions];
expectedP0 = N0/D0;
expectedP2 = (D0 N2 - 2 D2 N0)/D0^2;
expectedP4 = (D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0)/D0^3;
resP0 = FullSimplify[p0 - expectedP0, $Assumptions];
resP2 = FullSimplify[p2 - expectedP2, $Assumptions];
resP4 = FullSimplify[p4 - expectedP4, $Assumptions];
p0Match = TrueQ[resP0 == 0];
p2Match = TrueQ[resP2 == 0];
p4Match = TrueQ[resP4 == 0];
plainP2 = FullSimplify[Coefficient[plainSeries, omega, 2], $Assumptions];
plainEqualsCorrectP2 = TrueQ[FullSimplify[plainP2 == expectedP2, $Assumptions]];

sampleRules = {
  a -> 3.0, cs -> 2.0, c -> 5.0, G -> 7.0,
  D0 -> 19.0, D2 -> 23.0, D4 -> 29.0,
  N0 -> 11.0, N2 -> 13.0, N4 -> 17.0
};
evalSample[expr_] := N[expr /. sampleRules, 17];

zeroDim = {0, 0, 0};
dimAdd[x_, y_] := x + y;
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
expText[e_] := Module[{r = Rationalize[e]},
  If[Denominator[r] == 1, ToString[Numerator[r]], ToString[Numerator[r]] <> "/" <> ToString[Denominator[r]]]
];
dimText[d_] := Module[{pairs, parts},
  pairs = {{"L", d[[1]]}, {"T", d[[3]]}, {"M", d[[2]]}};
  parts = (If[TrueQ[#[[2]] == 1], #[[1]], #[[1]] <> "^" <> expText[#[[2]]]] &) /@
    Select[pairs, ! TrueQ[#[[2]] == 0] &];
  If[Length[parts] == 0, "1", StringRiffle[parts, " "]]
];

rawDims = <|
  a -> {1, 0, 0},
  cs -> {1, 0, -1},
  c -> {1, 0, -1},
  G -> {3, -1, -2},
  omega -> {0, 0, -1},
  D0 -> {-1, 1, -2},
  N0 -> {-1, 1, 0},
  chiQsym -> zeroDim,
  mtilde0 -> zeroDim
|>;
P0Raw = N0/D0;
frequencyNormalization = (cs/a)^2;
P0Physical = frequencyNormalization P0Raw;
YhatPhysical = 1 + out["u2"] omega^2 + out["u4"] omega^4 + I out["v5"] omega^5;
Gamma5 = chiQsym P0Physical a^5/(27 cs^5);
targetRHS = 54 G cs^5/(5 a^5 c^5);
p0RawDim = dimOf[P0Raw, rawDims];
frequencyNormDim = dimOf[frequencyNormalization, rawDims];
p0Dim = dimOf[P0Physical, rawDims];
dimensionalOk = TrueQ[p0Dim == zeroDim];
dropNormDim = dimOf[P0Raw, rawDims];
dropNormOk = TrueQ[dropNormDim == zeroDim];
dropNormVerdict = If[dropNormOk, "NO_FAIL", "FAIL_DIMENSIONAL"];
corruptN0Dims = Join[KeyDrop[rawDims, N0], <|N0 -> zeroDim|>];
corruptN0RawDim = dimOf[P0Raw, corruptN0Dims];
corruptN0P0Dim = dimOf[P0Physical, corruptN0Dims];
corruptN0Ok = TrueQ[corruptN0P0Dim == zeroDim];
corruptN0Verdict = If[corruptN0Ok, "NO_FAIL", "FAIL_DIMENSIONAL"];
rhsDim = dimOf[targetRHS, rawDims];
muDim = (rhsDim - p0Dim)/2;
dims = Join[rawDims, <|muHat0 -> muDim|>];
lhs = (muHat0 mtilde0)^2 P0Physical;
lhsRawMutation = (muHat0 mtilde0)^2 P0Raw;
lhsDim = dimOf[lhs, dims];
lhsRawDim = dimOf[lhsRawMutation, dims];
requiredP0Dim = rhsDim - 2 muDim;
gamma5Dim = dimOf[Gamma5, dims];
yhatDim = dimOf[YhatPhysical, dims];
homogeneityPass = TrueQ[lhsDim == rhsDim && p0Dim == requiredP0Dim && yhatDim == zeroDim];

headlineOk = And[
  u2Match, u4Match, v5Match, chiMatch,
  p0Match, p2Match, p4Match,
  dimensionalOk, ! dropNormOk, ! corruptN0Ok
];

lines = {
  "schema: pathA_33_mathematica_scratch/v1",
  "engine: mathematica",
  "fingerprint:",
  "  lambda_out_series_z: " <> quoteText[out["lambdaSeries"]],
  "  Yhat_out_series_z: " <> quoteText[out["yhatSeries"]],
  "  coefficients_z:",
  "    u2: " <> quoteText[out["u2z"]],
  "    u4: " <> quoteText[out["u4z"]],
  "    v5: " <> quoteText[out["v5z"]],
  "  coefficients_physical:",
  "    u2: " <> quoteText[out["u2"]],
  "    u4: " <> quoteText[out["u4"]],
  "    v5: " <> quoteText[out["v5"]],
  "  coefficients_z_numeric:",
  "    u2: " <> numText[out["u2z"]],
  "    u4: " <> numText[out["u4z"]],
  "    v5: " <> numText[out["v5z"]],
  "  coefficient_matches:",
  "    u2_z: " <> boolText[u2Match],
  "    u4_z: " <> boolText[u4Match],
  "    v5_z: " <> boolText[v5Match],
  "  chi_Q: " <> quoteText[chiQ],
  "  chi_Q_numeric: " <> numText[chiQ],
  "  incoming_chi_Q: " <> quoteText[chiQIncoming],
  "  incoming_chi_Q_numeric: " <> numText[chiQIncoming],
  "  standing_lambda_static: " <> quoteText[Coefficient[standing["lambdaSeries"], z, 0]],
  "  standing_Yhat_static: " <> quoteText[standing["static"]],
  "  standing_v5_z: " <> quoteText[standing["v5z"]],
  "prefactor:",
  "  correct_series: " <> quoteText[prefSeries],
  "  coefficients:",
  "    P0: " <> quoteText[p0],
  "    P2: " <> quoteText[p2],
  "    P4: " <> quoteText[p4],
  "  residuals_to_formula:",
  "    P0: " <> quoteText[resP0],
  "    P2: " <> quoteText[resP2],
  "    P4: " <> quoteText[resP4],
  "  matches:",
  "    P0: " <> boolText[p0Match],
  "    P2: " <> boolText[p2Match],
  "    P4: " <> boolText[p4Match],
  "  plain_P2: " <> quoteText[plainP2],
  "  plain_equals_correct_P2: " <> boolText[plainEqualsCorrectP2],
  "  sample_values:",
  "    P0: " <> numText[evalSample[p0]],
  "    P2: " <> numText[evalSample[p2]],
  "    P4: " <> numText[evalSample[p4]],
  "    plain_P2: " <> numText[evalSample[plainP2]],
  "    plain_minus_correct_P2: " <> numText[evalSample[plainP2 - expectedP2]],
  "dimension:",
  "  P0_raw_dimension: " <> quoteText[dimText[p0RawDim]],
  "  frequency_normalization_dimension: " <> quoteText[dimText[frequencyNormDim]],
  "  P0_dimension: " <> quoteText[dimText[p0Dim]],
  "  P0_physical_dimension: " <> quoteText[dimText[p0Dim]],
  "  dimensional_ok: " <> boolText[dimensionalOk],
  "  mu_hat0_dimension: " <> quoteText[dimText[muDim]],
  "  Gamma5_dimension: " <> quoteText[dimText[gamma5Dim]],
  "  lhs_dimension: " <> quoteText[dimText[lhsDim]],
  "  rhs_dimension: " <> quoteText[dimText[rhsDim]],
  "  mu_hat0_homogeneity_pass: " <> boolText[homogeneityPass],
  "  mu_hat0_homogeneity_label: " <> quoteText["non-able-to-fail (mu_hat0 free carrier)"],
  "  drop_norm_verdict: " <> quoteText[dropNormVerdict],
  "  corrupt_N0_verdict: " <> quoteText[corruptN0Verdict],
  "headline_booleans:",
  "  fingerprint_ok: " <> boolText[u2Match && u4Match && v5Match && chiMatch],
  "  prefactor_ok: " <> boolText[p0Match && p2Match && p4Match],
  "  dimension_ok: " <> boolText[dimensionalOk]
};

Export[yamlOut, StringRiffle[lines, "\n"] <> "\n", "Text"];
If[! TrueQ[headlineOk], fail["headline checks failed"]];
Print["Wrote Mathematica scratch: ", yamlOut];
Exit[0];
