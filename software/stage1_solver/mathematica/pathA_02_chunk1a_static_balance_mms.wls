#!/usr/bin/env wolframscript

(* Path-A chunk 1a static-balance MMS verifier.
   Mathematica is used only as an independent continuum-source engine. *)

$HistoryLength = 0;
ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/pathA_02_chunk1a_static_balance_mms.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
runDir = FileNameJoin[{scriptDir, "runs", "pathA_02_chunk1a_static_balance_mms"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

diagnosticsPath = FileNameJoin[{runDir, "pathA_02_chunk1a_static_balance_mms_diagnostics.json"}];

wMin = 0.0;
wMax = 1.4;
length = wMax - wMin;

payload = Module[
  {
    z, rho, x, r0, tw, u, twR, uR, r0w, fluxDivergence,
    gradientSquare, potentialGradient, source, sampleW, sourceSamples,
    fluxSamples, gradientSamples, potentialSamples, checks
  },
  x = (z - wMin)/length;
  r0 = 1.08
    + 0.12*256*x^4*(1 - x)^4
    + 0.025*Sin[Pi*x]^4;
  tw = 1.32
    + 0.18*Sin[2*Pi*x]
    + 0.11*rho
    + 0.04*rho^2
    + 0.03*rho*Cos[Pi*x];
  u = 0.42
    + 0.20*rho^2
    + 0.03*rho^3
    + 0.05*rho*Sin[Pi*x];

  twR = D[tw, rho];
  uR = D[u, rho];
  r0w = D[r0, z];
  fluxDivergence = -D[(tw /. rho -> r0)*r0w, z];
  gradientSquare = (1/2)*(twR /. rho -> r0)*r0w^2;
  potentialGradient = uR /. rho -> r0;
  source = fluxDivergence + gradientSquare + potentialGradient;

  sampleW = N[wMin + length*{0.07, 0.19, 0.37, 0.53, 0.71, 0.89}, 17];
  sourceSamples = Chop[Re[N[source /. z -> #, 17]]] & /@ sampleW;
  fluxSamples = Chop[Re[N[fluxDivergence /. z -> #, 17]]] & /@ sampleW;
  gradientSamples = Chop[Re[N[gradientSquare /. z -> #, 17]]] & /@ sampleW;
  potentialSamples = Chop[Re[N[potentialGradient /. z -> #, 17]]] & /@ sampleW;

  checks = {
    <|"check" -> "manufactured_radius_nonconstant", "pass" -> TrueQ[Max[Abs[Chop[Re[N[r0w /. z -> #]]] & /@ sampleW]] > 0]|>,
    <|"check" -> "tw_depends_on_radius", "pass" -> TrueQ[Max[Abs[Chop[Re[N[twR /. {z -> #, rho -> 1.1}]]] & /@ sampleW]] > 0]|>,
    <|"check" -> "u_gradient_nonzero", "pass" -> TrueQ[Max[Abs[potentialSamples]] > 0]|>,
    <|"check" -> "gradient_square_nonzero", "pass" -> TrueQ[Max[Abs[gradientSamples]] > 0]|>
  };

  <|
    "engine" -> "Mathematica",
    "script" -> "pathA_02_chunk1a_static_balance_mms.wls",
    "domain" -> <|"w_min" -> wMin, "w_max" -> wMax|>,
    "manufactured_radius_expression_input_form" -> ToString[r0, InputForm],
    "tw_expression_input_form" -> ToString[tw, InputForm],
    "u_expression_input_form" -> ToString[u, InputForm],
    "source_expression_input_form" -> ToString[source, InputForm],
    "sample_w" -> sampleW,
    "source_samples" -> sourceSamples,
    "term_samples" -> <|
      "flux_divergence" -> fluxSamples,
      "gradient_square" -> gradientSamples,
      "potential_gradient" -> potentialSamples
    |>,
    "genuine_checks" -> checks,
    "passed" -> And @@ (TrueQ[#["pass"]] & /@ checks)
  |>
];

exportResult = Export[diagnosticsPath, payload, "JSON"];
Print["wrote diagnostics: ", diagnosticsPath];
If[exportResult === $Failed, Exit[1]];
If[! TrueQ[payload["passed"]], Exit[1]];
Exit[0];
