(* pathA_30 frozen-wall D/N unit test, Mathematica transfer-matrix engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
boolText[x_] := If[TrueQ[x], "true", "false"];
quoteText[x_] := "'" <> StringReplace[If[StringQ[x], x, ToString[x, InputForm]], "'" -> "''"] <> "'";

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_30_dn_unit_test.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyExprFile = FileNameJoin[{scratchDir, "pathA_30_sympy_exprs.wl"}];
yamlOut = FileNameJoin[{scratchDir, "pathA_30_mathematica_results.yaml"}];

If[! FileExistsQ[sympyExprFile], fail["missing SymPy expression export: " <> sympyExprFile]];
Get[sympyExprFile];

$Assumptions = L0 > 0 && cS > 0 && omega > 0 && Element[alpha, Reals] &&
  Element[{psi0, p0}, Reals] && psi0 != 0;

k = omega/cS;
transferMatrix = {
  {Cos[k L0], Sin[k L0]/k},
  {-k Sin[k L0], Cos[k L0]}
};
mouthState = {psi0, p0};
capState = FullSimplify[transferMatrix . mouthState];

neumannCapEquation = FullSimplify[capState[[2]] == 0];
neumannResidual = FullSimplify[capState[[2]]];
p0FromNeumann = FullSimplify[-Coefficient[neumannResidual, psi0] psi0/Coefficient[neumannResidual, p0]];
p0OverPsi0N = FullSimplify[p0FromNeumann/psi0];
dtnTransfer = FullSimplify[-p0OverPsi0N];
poleDenominatorTransfer = FullSimplify[Cos[k L0]];

robinCapEquation = FullSimplify[capState[[2]] == alpha capState[[1]]];
robinResidual = FullSimplify[capState[[2]] - alpha capState[[1]]];
p0FromRobin = FullSimplify[-Coefficient[robinResidual, psi0] psi0/Coefficient[robinResidual, p0]];
p0OverPsi0Robin = FullSimplify[p0FromRobin/psi0];
robinDtnTransfer = FullSimplify[-p0OverPsi0Robin];
robinDenominatorTransfer = FullSimplify[k Cos[k L0] - alpha Sin[k L0]];
robinAlpha0 = FullSimplify[robinDtnTransfer /. alpha -> 0];
robinAlphaInf = FullSimplify[Limit[robinDtnTransfer, alpha -> Infinity]];
ddTransfer = FullSimplify[k Cot[k L0]];
staticSeriesTransfer = FullSimplify[Normal[Series[dtnTransfer, {omega, 0, 5}]]];

dtnAgree = TrueQ[FullSimplify[dtnTransfer - sympyDtn == 0]];
poleDenAgree = TrueQ[FullSimplify[poleDenominatorTransfer - sympyPoleDenominator == 0]];
robinDtnAgree = TrueQ[FullSimplify[robinDtnTransfer - sympyRobinDtn == 0]];
robinDenAgree = TrueQ[FullSimplify[robinDenominatorTransfer - sympyRobinDenominatorCore == 0]];
staticAgree = TrueQ[FullSimplify[staticSeriesTransfer - sympyStaticSeriesPoly == 0]];
ddAgree = TrueQ[FullSimplify[robinAlphaInf - sympyDDBranchDtn == 0]] &&
  TrueQ[FullSimplify[robinAlphaInf - ddTransfer == 0]];
alpha0Agree = TrueQ[FullSimplify[robinAlpha0 - dtnTransfer == 0]];

If[! And[dtnAgree, poleDenAgree, robinDtnAgree, robinDenAgree, staticAgree, ddAgree, alpha0Agree],
  fail["Mathematica transfer route disagrees with SymPy dsolve artifacts"]
];

zeroDim = {0, 0, 0};
dimAdd[x_, y_] := x + y;
dimSub[x_, y_] := x - y;
dimScale[x_, q_] := q x;
dimOf[expr_, dims_] := Module[{args, ds, base, pow, argDims},
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
    MemberQ[{Sin, Cos, Tan, Cot, Sinh, Cosh, Tanh, Coth, Sech, Csch}, Head[expr]],
      argDims = dimOf[#, dims] & /@ (List @@ expr);
      If[AnyTrue[argDims, # =!= zeroDim &], fail["dimensionful argument in dimensionless function"]];
      zeroDim,
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
dimTupleText[d_] := "(" <> StringRiffle[expText /@ d, ","] <> ")";
verdictFromDim[ok_] := If[TrueQ[ok], "DIMENSIONAL_OK", "DN_UNITTEST_FAIL_DIMENSIONAL"];

spatialDimensions = 4;
lengthDim = {1, 0, 0};
energyDim = {2, 1, -2};
fourVolumeDim = dimScale[lengthDim, spatialDimensions];
Pdim = dimSub[energyDim, fourVolumeDim];
rhoDim = dimScale[lengthDim, -spatialDimensions];
Kdim = dimSub[Pdim, 5 rhoDim];

dimRules = <|
  L0 -> lengthDim,
  omega -> {0, 0, -1},
  cS -> {1, 0, -1},
  K -> Kdim,
  rhoStar -> rhoDim,
  m -> {0, 1, 0}
|>;
csSquaredFromEOS = 5 K rhoStar^4/m;
tanArgument = omega L0/cS;
dtnPrefactor = -omega/cS;
csDim = dimOf[csSquaredFromEOS, dimRules];
tanArgumentDim = dimOf[tanArgument, dimRules];
z00Dim = dimOf[dtnTransfer, dimRules];
dimensionalOk = TrueQ[csDim == {2, 0, -2} && tanArgumentDim == zeroDim && z00Dim == {-1, 0, 0}];
corruptKRules = Join[KeyDrop[dimRules, K], <|K -> dimAdd[dimRules[K], {1, 0, 0}]|>];
corruptCsDim = dimOf[csSquaredFromEOS, corruptKRules];
corruptTanArgumentDim = dimOf[tanArgument, corruptKRules];
corruptZ00Dim = dimOf[dtnTransfer, corruptKRules];
corruptKOk = TrueQ[corruptCsDim == {2, 0, -2} && corruptTanArgumentDim == zeroDim && corruptZ00Dim == {-1, 0, 0}];
corruptKProbeVerdict = If[corruptKOk, "NO_FAIL", "DN_UNITTEST_FAIL_DIMENSIONAL"];

If[! TrueQ[dimensionalOk && corruptKProbeVerdict === "DN_UNITTEST_FAIL_DIMENSIONAL"],
  fail["Mathematica dimensional gate did not pass baseline and fail corrupted K probe"]
];

yaml = StringRiffle[{
  "schema: pathA_30_dn_unit_test_mathematica/v1",
  "route: transfer_matrix_resolvent",
  "sympy_expression_digest: " <> quoteText[sympyExpressionDigest],
  "transfer_matrix: " <> quoteText[transferMatrix],
  "neumann_cap_equation: " <> quoteText[neumannCapEquation],
  "robin_cap_equation: " <> quoteText[robinCapEquation],
  "dtn_transfer: " <> quoteText[dtnTransfer],
  "pole_denominator_transfer: " <> quoteText[poleDenominatorTransfer],
  "robin_dtn_transfer: " <> quoteText[robinDtnTransfer],
  "robin_denominator_transfer: " <> quoteText[robinDenominatorTransfer],
  "static_series_transfer: " <> quoteText[staticSeriesTransfer],
  "dimensional:",
  "  dimension_order:",
  "  - L",
  "  - M",
  "  - T",
  "  headline_quantities_walked:",
  "    cs_squared_from_EOS: " <> quoteText[csSquaredFromEOS],
  "    tan_argument: " <> quoteText[tanArgument],
  "    Z00: " <> quoteText[dtnTransfer],
  "  symbol_dimensions:",
  "    L0: " <> quoteText["L"],
  "    omega: " <> quoteText["T^-1"],
  "    c_s: " <> quoteText["L T^-1"],
  "    m: " <> quoteText["M"],
  "    energy: " <> quoteText[dimText[energyDim]],
  "    four_volume_L4: " <> quoteText[dimText[fourVolumeDim]],
  "    P: " <> quoteText[dimText[Pdim]],
  "    rho_star: " <> quoteText[dimText[rhoDim]],
  "    K: " <> quoteText[dimText[Kdim]],
  "  sourcing_note:",
  "    K_source: " <> quoteText["P=K*rho^5"],
  "    spatial_dimensions: " <> ToString[spatialDimensions],
  "    no_cs_dependency: true",
  "    energy_dimension_LMT:",
  "    - " <> quoteText[expText[energyDim[[1]]]],
  "    - " <> quoteText[expText[energyDim[[2]]]],
  "    - " <> quoteText[expText[energyDim[[3]]]],
  "    four_volume_dimension_LMT:",
  "    - " <> quoteText[expText[fourVolumeDim[[1]]]],
  "    - " <> quoteText[expText[fourVolumeDim[[2]]]],
  "    - " <> quoteText[expText[fourVolumeDim[[3]]]],
  "    P_dimension_LMT:",
  "    - " <> quoteText[expText[Pdim[[1]]]],
  "    - " <> quoteText[expText[Pdim[[2]]]],
  "    - " <> quoteText[expText[Pdim[[3]]]],
  "    rho_dimension_LMT:",
  "    - " <> quoteText[expText[rhoDim[[1]]]],
  "    - " <> quoteText[expText[rhoDim[[2]]]],
  "    - " <> quoteText[expText[rhoDim[[3]]]],
  "    K_derivation: " <> quoteText["[K]=[P]-5[rho]"],
  "    K_dimension_LMT:",
  "    - " <> quoteText[expText[Kdim[[1]]]],
  "    - " <> quoteText[expText[Kdim[[2]]]],
  "    - " <> quoteText[expText[Kdim[[3]]]],
  "    derived_chain: " <> quoteText[
    "[P]=" <> dimTupleText[Pdim] <> ", [rho]=" <> dimTupleText[rhoDim] <>
      ", [K]=[P]-5[rho]=" <> dimTupleText[Kdim]
  ],
  "  computed_dimensions:",
  "    cs_squared_from_EOS: " <> quoteText[dimText[csDim]],
  "    tan_argument: " <> quoteText[dimText[tanArgumentDim]],
  "    Z00: " <> quoteText[dimText[z00Dim]],
  "  dimensional_ok: " <> boolText[dimensionalOk],
  "  dimensional_status: " <> quoteText[verdictFromDim[dimensionalOk]],
  "  DN_UNITTEST_FAIL_DIMENSIONAL_probe:",
  "    mutation: " <> quoteText["corrupt sourced [K] by one extra power of L"],
  "    sourced_K_dimension: " <> quoteText[dimText[dimRules[K]]],
  "    corrupted_K_dimension: " <> quoteText[dimText[corruptKRules[K]]],
  "    with_mutation_dimensional_ok: " <> boolText[corruptKOk],
  "    without_mutation_dimensional_ok: " <> boolText[dimensionalOk],
  "    probe_verdict: " <> quoteText[corruptKProbeVerdict],
  "    mutated_dimensions:",
  "      cs_squared_from_EOS: " <> quoteText[dimText[corruptCsDim]],
  "      tan_argument: " <> quoteText[dimText[corruptTanArgumentDim]],
  "      Z00: " <> quoteText[dimText[corruptZ00Dim]],
  "engine_agreement:",
  "  dtn: " <> boolText[dtnAgree],
  "  pole_denominator: " <> boolText[poleDenAgree],
  "  robin_dtn: " <> boolText[robinDtnAgree],
  "  robin_denominator: " <> boolText[robinDenAgree],
  "  static_series: " <> boolText[staticAgree],
  "  dd_limit: " <> boolText[ddAgree],
  "  alpha0_limit: " <> boolText[alpha0Agree],
  "  dimensional_ok: " <> boolText[dimensionalOk],
  "  dimension_probe_verdict: " <> boolText[corruptKProbeVerdict === "DN_UNITTEST_FAIL_DIMENSIONAL"]
}, "\n"] <> "\n";

Export[yamlOut, yaml, "Text"];
Print["verdict: mathematica_transfer_engine_pass"];
Print["yaml: ", yamlOut];
Exit[0];
