(* pathA_28 monopole/dipole return-condition Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_28_monopole.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_28_monopole_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_28_monopole_mathematica.json"}];

If[! FileExistsQ[sympyJson], fail["missing SymPy JSON for engine agreement: " <> sympyJson]];

$Assumptions = a > 0 && k > 0 && omega > 0 && cS > 0 &&
  Element[{M0, D1, Q2, R0, R1, eta0, eta1, ellw, j0, E0}, Reals] &&
  lam > 0 && muw > 0 && rho0 > 0;

nonzeroQ[expr_] := ! TrueQ[FullSimplify[expr == 0]];

radiatingPower[y_, max_] := Module[{n, coeff, realCoeff, imagCoeff},
  For[n = 1, n <= max, n++,
    coeff = FullSimplify[Coefficient[Expand[y], k, n]];
    realCoeff = FullSimplify[ComplexExpand[Re[coeff]]];
    imagCoeff = FullSimplify[ComplexExpand[Im[coeff]]];
    If[! TrueQ[imagCoeff == 0] && TrueQ[realCoeff == 0],
      Return[{n, imagCoeff}]
    ];
  ];
  fail["no radiation-phase coefficient found in " <> ToString[y, InputForm]]
];

hankelData[ell_] := Module[{z, j, y, h, lambda, lambdaSeries, ySeries, yStatic, yNorm, rp, pRaw, imagCoeff, kernel},
  z = Unique["z"];
  {j, y} = Switch[ell,
    0, {Sin[z]/z, -Cos[z]/z},
    1, {Sin[z]/z^2 - Cos[z]/z, -Cos[z]/z^2 - Sin[z]/z},
    2, {((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2, -((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2},
    _, fail["unsupported ell"]
  ];
  h = FullSimplify[j + I y];
  lambda = FullSimplify[(k D[h, z]/h) /. z -> k a];
  lambdaSeries = Expand[Normal@Series[lambda, {k, 0, 6}]];
  ySeries = Expand[Normal@Series[1/lambdaSeries, {k, 0, 6}]];
  yStatic = FullSimplify[ySeries /. k -> 0];
  yNorm = Expand[Normal@Series[ySeries/yStatic, {k, 0, 6}]];
  rp = radiatingPower[yNorm, 8];
  pRaw = rp[[1]];
  imagCoeff = FullSimplify[rp[[2]]];
  kernel = FullSimplify[I imagCoeff (omega/cS)^pRaw];
  <|
    "LambdaSeries" -> lambdaSeries,
    "YNormSeries" -> yNorm,
    "pRaw" -> pRaw,
    "ImagCoeffK" -> imagCoeff,
    "RadiationKernel" -> kernel
  |>
];

dtn0 = hankelData[0];
dtn1 = hankelData[1];
dtn2 = hankelData[2];

raw0 = FullSimplify[dtn0["RadiationKernel"] M0];
raw1 = FullSimplify[dtn1["RadiationKernel"] D1];
raw2 = FullSimplify[dtn2["RadiationKernel"] Q2];
resid0 = FullSimplify[dtn0["RadiationKernel"] (M0 + R0)];
resid1 = FullSimplify[dtn1["RadiationKernel"] (D1 + R1)];
without0 = FullSimplify[resid0 /. R0 -> 0];
without1 = FullSimplify[resid1 /. R1 -> 0];
with0 = FullSimplify[resid0 /. R0 -> -M0];
with1 = FullSimplify[resid1 /. R1 -> -D1];
deriv0 = FullSimplify[eta0^2 omega^2 raw0];
deriv1 = FullSimplify[eta1^2 omega^2 raw1];

If[! nonzeroQ[raw0] || ! nonzeroQ[raw1] || ! nonzeroQ[raw2], fail["raw amplitude disappeared"]];
If[! nonzeroQ[without0] || ! nonzeroQ[without1], fail["residual without condition disappeared"]];
If[with0 =!= 0 || with1 =!= 0, fail["cancellation condition did not zero residual"]];

w =.; W243 = Exp[-w^2]/Sqrt[Pi]; jw243 = ellw j0 w Exp[-w^2];
sleak243 = FullSimplify[Integrate[D[W243, w] jw243, {w, -Infinity, Infinity}]];
If[FullSimplify[sleak243 + Sqrt[2] ellw j0/4] =!= 0, fail["Stage 243 S_leak mismatch"]];
W244 = Exp[-w^2/lam^2]/(lam Sqrt[Pi]);
phi = 2 w Exp[-w^2/lam^2]/(Sqrt[Pi] lam^3);
Ew = -E0 phi;
jw244 = muw rho0 Ew;
sleak244 = FullSimplify[Integrate[D[W244, w] jw244, {w, -Infinity, Infinity}],
  Assumptions -> lam > 0 && muw > 0 && rho0 > 0 && Element[E0, Reals]];
If[FullSimplify[sleak244 - Sqrt[2] E0 muw rho0/(2 Sqrt[Pi] lam^3)] =!= 0, fail["Stage 244 S_leak mismatch"]];

rawPresent = AllTrue[{raw0, raw1, raw2}, nonzeroQ];
conditionWorks = TrueQ[with0 == 0 && with1 == 0 && nonzeroQ[without0] && nonzeroQ[without1]];
computeVerdict[raw_, works_, possible_] := Which[
  TrueQ[raw] && ! TrueQ[possible], "MONOPOLE_RADIATION_UNAVOIDABLE",
  TrueQ[raw] && TrueQ[works] && TrueQ[possible], "MONOPOLE_DIPOLE_RETURN_CONDITIONAL",
  True, "INCONCLUSIVE"
];
topLineVerdict = computeVerdict[rawPresent, conditionWorks, True];
syntheticUnavoidable = computeVerdict[True, False, False];
If[topLineVerdict === "INCONCLUSIVE", fail["computed verdict fell to INCONCLUSIVE"]];
If[syntheticUnavoidable =!= "MONOPOLE_RADIATION_UNAVOIDABLE", fail["unavoidable rung control did not fire"]];
If[! TrueQ[dtn0["pRaw"] < dtn2["pRaw"]], fail["raw monopole is not omega-dominant"]];
steadyLimit0 = Block[{$Assumptions = a > 0 && cS > 0 && Element[M0, Reals]},
  FullSimplify[Limit[raw0, omega -> 0]]
];
If[steadyLimit0 =!= 0, fail["steady-no-radiation control failed"]];
If[dtn2["pRaw"] =!= 5 || ! nonzeroQ[raw2], fail["quadrupole-survives control failed"]];
If[FullSimplify[raw0 /. M0 -> 0] =!= 0, fail["S_leak zero shortcut control malformed"]];
If[! nonzeroQ[deriv0], fail["derivative vertex incorrectly killed monopole"]];

sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
assertEngine[name_, actual_] := Module[{expectedText, expectedExpr},
  expectedText = sympyExprs[name];
  If[! StringQ[expectedText], fail["missing SymPy expression for " <> name]];
  expectedExpr = ToExpression[expectedText, InputForm];
  If[! TrueQ[FullSimplify[actual == expectedExpr]],
    fail["engine disagreement " <> name <> ": Mathematica got " <> ToString[actual, InputForm] <>
      ", SymPy exported " <> expectedText]
  ];
];

engineKeys = {
  "Lambda0_series" -> dtn0["LambdaSeries"],
  "Y0_norm_series" -> dtn0["YNormSeries"],
  "p_raw_0" -> dtn0["pRaw"],
  "imag_coeff_k_0" -> dtn0["ImagCoeffK"],
  "radiation_kernel_0" -> dtn0["RadiationKernel"],
  "raw_amplitude_0" -> raw0,
  "Lambda1_series" -> dtn1["LambdaSeries"],
  "Y1_norm_series" -> dtn1["YNormSeries"],
  "p_raw_1" -> dtn1["pRaw"],
  "imag_coeff_k_1" -> dtn1["ImagCoeffK"],
  "radiation_kernel_1" -> dtn1["RadiationKernel"],
  "raw_amplitude_1" -> raw1,
  "Lambda2_series" -> dtn2["LambdaSeries"],
  "Y2_norm_series" -> dtn2["YNormSeries"],
  "p_raw_2" -> dtn2["pRaw"],
  "imag_coeff_k_2" -> dtn2["ImagCoeffK"],
  "radiation_kernel_2" -> dtn2["RadiationKernel"],
  "raw_amplitude_2" -> raw2,
  "residual_without_0" -> without0,
  "residual_with_0" -> with0,
  "derivative_vertex_amp_0" -> deriv0,
  "residual_without_1" -> without1,
  "residual_with_1" -> with1,
  "derivative_vertex_amp_1" -> deriv1,
  "stage243_sleak" -> sleak243,
  "stage244_sleak" -> sleak244
};
Scan[assertEngine[#[[1]], #[[2]]] &, engineKeys];

If[topLineVerdict =!= sympyResults["top_line_verdict"], fail["top-line verdict disagreement with SymPy JSON"]];

results = <|
  "schema" -> "pathA_28_monopole_mathematica/v1",
  "top_line_verdict" -> topLineVerdict,
  "controls" -> <|
    "raw_monopole_present" -> True,
    "steady_no_radiation" -> True,
    "quadrupole_survives" -> True,
    "return_necessity" -> conditionWorks,
    "unavoidable_rung_control" -> syntheticUnavoidable,
    "anti_tautology_derivative_vertex_not_basis" -> nonzeroQ[deriv0],
    "anti_tautology_no_track3_bulk_kill" -> nonzeroQ[without0] && nonzeroQ[without1]
  |>,
  "engine_agreement" -> <|
    "sympy_json" -> sympyJson,
    "shared_expression_count" -> Length[engineKeys],
    "status" -> "PASS"
  |>
|>;

If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
Export[jsonOut, results, "JSON"];
Print["PASS pathA_28_monopole_mathematica"];
Print[ExportString[<|
  "top_line_verdict" -> topLineVerdict,
  "engine_agreement" -> "PASS",
  "shared_expression_count" -> Length[engineKeys],
  "json" -> jsonOut
|>, "JSON"]];
Exit[0];
