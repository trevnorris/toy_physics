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
  "engine_agreement:",
  "  dtn: " <> boolText[dtnAgree],
  "  pole_denominator: " <> boolText[poleDenAgree],
  "  robin_dtn: " <> boolText[robinDtnAgree],
  "  robin_denominator: " <> boolText[robinDenAgree],
  "  static_series: " <> boolText[staticAgree],
  "  dd_limit: " <> boolText[ddAgree],
  "  alpha0_limit: " <> boolText[alpha0Agree]
}, "\n"] <> "\n";

Export[yamlOut, yaml, "Text"];
Print["verdict: mathematica_transfer_engine_pass"];
Print["yaml: ", yamlOut];
Exit[0];
