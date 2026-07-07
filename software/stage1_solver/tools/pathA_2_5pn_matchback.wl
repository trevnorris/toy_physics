(* A3 2.5PN match-back verification, Wolfram Language engine.

   Self-contained by design: this script restates the calibrated literals
   locally, does not read or write repo artifacts, and does not consume the
   SymPy peer in any way.
*)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

$Assumptions = G > 0 && cs > 0 && a > 0 && c > 0;

residualOrder = {
  "INV1_moment_invariant",
  "INV2_pathA43_form_to_BT",
  "INV3_corpus_form_to_BT",
  "INV4_redundant_cross_form",
  "INV5_Kbar0_coeff_54_5",
  "INV5_Kbar2_coeff_6_5",
  "INV5_Kbar4_coeff_8_15",
  "INV5_BT_coeff_2_5",
  "INV5_pathA43_denominator_27",
  "INV5_corpus_factor_9",
  "INV5_exp_Kbar2_5_2",
  "INV5_exp_Kbar0_3_2"
};

baseConfig = <|
  "k0Coeff" -> 54/5,
  "k2Coeff" -> 6/5,
  "k4Coeff" -> 8/15,
  "btCoeff" -> 2/5,
  "pathDenominator" -> 27,
  "corpusFactor" -> 9,
  "expK2" -> 5/2,
  "expK0" -> 3/2,
  "k0Scale" -> 1,
  "k2Scale" -> 1,
  "k4Scale" -> 1,
  "btScale" -> 1
|>;

mergeConfig[updates_Association] := Association[Normal[baseConfig], Normal[updates]];

mutations = {
  <|"name" -> "Kbar4_coeff_8_15_to_8_14", "config" -> mergeConfig[<|"k4Coeff" -> 8/14|>]|>,
  <|"name" -> "Kbar4_sign_flip", "config" -> mergeConfig[<|"k4Coeff" -> -8/15|>]|>,
  <|"name" -> "Kbar2_coeff_6_5_to_7_5", "config" -> mergeConfig[<|"k2Coeff" -> 7/5|>]|>,
  <|"name" -> "Kbar0_coeff_54_5_to_55_5", "config" -> mergeConfig[<|"k0Coeff" -> 55/5|>]|>,
  <|"name" -> "pathA43_denominator_27_to_26", "config" -> mergeConfig[<|"pathDenominator" -> 26|>]|>,
  <|"name" -> "corpus_factor_9_to_8", "config" -> mergeConfig[<|"corpusFactor" -> 8|>]|>,
  <|"name" -> "exp_Kbar2_5_2_to_3_2", "config" -> mergeConfig[<|"expK2" -> 3/2|>]|>,
  <|"name" -> "exp_Kbar0_3_2_to_1", "config" -> mergeConfig[<|"expK0" -> 1|>]|>,
  <|"name" -> "BT_coeff_2_5_to_3_5", "config" -> mergeConfig[<|"btCoeff" -> 3/5|>]|>,
  <|
    "name" -> "coherent_scale_all_moments_and_BT_x2",
    "config" -> mergeConfig[<|"k0Scale" -> 2, "k2Scale" -> 2, "k4Scale" -> 2, "btScale" -> 2|>]
  |>,
  <|
    "name" -> "coupled_moments_x2_BT_fixed",
    "config" -> mergeConfig[<|"k0Scale" -> 2, "k2Scale" -> 2, "k4Scale" -> 2|>]
  |>
};

expectedCaughtBy = <|
  "Kbar4_coeff_8_15_to_8_14" -> {
    "INV1_moment_invariant",
    "INV5_Kbar4_coeff_8_15"
  },
  "Kbar4_sign_flip" -> {
    "INV1_moment_invariant",
    "INV5_Kbar4_coeff_8_15"
  },
  "Kbar2_coeff_6_5_to_7_5" -> {
    "INV1_moment_invariant",
    "INV3_corpus_form_to_BT",
    "INV4_redundant_cross_form",
    "INV5_Kbar2_coeff_6_5"
  },
  "Kbar0_coeff_54_5_to_55_5" -> {
    "INV1_moment_invariant",
    "INV2_pathA43_form_to_BT",
    "INV3_corpus_form_to_BT",
    "INV4_redundant_cross_form",
    "INV5_Kbar0_coeff_54_5"
  },
  "pathA43_denominator_27_to_26" -> {
    "INV2_pathA43_form_to_BT",
    "INV4_redundant_cross_form",
    "INV5_pathA43_denominator_27"
  },
  "corpus_factor_9_to_8" -> {
    "INV3_corpus_form_to_BT",
    "INV4_redundant_cross_form",
    "INV5_corpus_factor_9"
  },
  "exp_Kbar2_5_2_to_3_2" -> {
    "INV3_corpus_form_to_BT",
    "INV4_redundant_cross_form",
    "INV5_exp_Kbar2_5_2"
  },
  "exp_Kbar0_3_2_to_1" -> {
    "INV3_corpus_form_to_BT",
    "INV4_redundant_cross_form",
    "INV5_exp_Kbar0_3_2"
  },
  "BT_coeff_2_5_to_3_5" -> {
    "INV2_pathA43_form_to_BT",
    "INV3_corpus_form_to_BT",
    "INV5_BT_coeff_2_5"
  },
  "coherent_scale_all_moments_and_BT_x2" -> {
    "INV5_Kbar0_coeff_54_5",
    "INV5_Kbar2_coeff_6_5",
    "INV5_Kbar4_coeff_8_15",
    "INV5_BT_coeff_2_5"
  },
  "coupled_moments_x2_BT_fixed" -> {
    "INV2_pathA43_form_to_BT",
    "INV3_corpus_form_to_BT",
    "INV5_Kbar0_coeff_54_5",
    "INV5_Kbar2_coeff_6_5",
    "INV5_Kbar4_coeff_8_15"
  }
|>;

stripConditionalZero[expr_] := expr /. ConditionalExpression[0, _] :> 0;
compact[expr_] := stripConditionalZero[FullSimplify[Cancel[Together[expr]], $Assumptions]];

assertNoFloat[label_, expr_] := Module[{floats},
  floats = Cases[HoldComplete[expr], _Real, Infinity];
  If[floats =!= {}, fail[label <> " contains machine real atoms: " <> ToString[floats, InputForm]]];
];

buildResiduals[cfg_Association] := Module[
  {k0, k2, k4, bt, pathForm, corpusForm, raw, reduced},
  k0 = cfg["k0Scale"] cfg["k0Coeff"] G cs^5/(a^5 c^5);
  k2 = cfg["k2Scale"] cfg["k2Coeff"] G cs^3/(a^3 c^5);
  k4 = cfg["k4Scale"] cfg["k4Coeff"] G cs/(a c^5);
  bt = cfg["btScale"] cfg["btCoeff"] G/c^5;
  pathForm = k0 a^5/(cfg["pathDenominator"] cs^5);
  corpusForm = cfg["corpusFactor"] k2^cfg["expK2"]/k0^cfg["expK0"];

  raw = <|
    "INV1_moment_invariant" -> k4 k0 - 4 k2^2,
    "INV2_pathA43_form_to_BT" -> pathForm - bt,
    "INV3_corpus_form_to_BT" -> corpusForm - bt,
    "INV4_redundant_cross_form" -> pathForm - corpusForm,
    "INV5_Kbar0_coeff_54_5" -> k0 a^5 c^5/(G cs^5) - 54/5,
    "INV5_Kbar2_coeff_6_5" -> k2 a^3 c^5/(G cs^3) - 6/5,
    "INV5_Kbar4_coeff_8_15" -> k4 a c^5/(G cs) - 8/15,
    "INV5_BT_coeff_2_5" -> bt c^5/G - 2/5,
    "INV5_pathA43_denominator_27" -> cfg["pathDenominator"] - 27,
    "INV5_corpus_factor_9" -> cfg["corpusFactor"] - 9,
    "INV5_exp_Kbar2_5_2" -> cfg["expK2"] - 5/2,
    "INV5_exp_Kbar0_3_2" -> cfg["expK0"] - 3/2
  |>;

  KeyValueMap[assertNoFloat[#1 <> " raw", #2] &, raw];
  reduced = Association[KeyValueMap[(#1 -> compact[#2]) &, raw]];
  KeyValueMap[assertNoFloat[#1 <> " reduced", #2] &, reduced];
  reduced
];

exprText[expr_] := ToString[expr, InputForm];
vectorText[res_Association] := "[" <> StringRiffle[(# <> "=" <> exprText[res[#]]) & /@ residualOrder, ", "] <> "]";
firedLabels[res_Association] := Select[residualOrder, ! TrueQ[res[#] === 0] &];

printProvenance[] := (
  Print["PROVENANCE:"];
  Print["  Kbar0,Kbar2,Kbar4 are CALIBRATED closure inputs; moments are not derived here."];
  Print["  External calibrated literals: 2/5, 54/5, and G; G=GENUINE_BLOCKED."];
  Print["  The pathA_43 denominator 27 is upstream-EARNED density-Hankel fingerprint, imported here."];
  Print["  Corpus form 9*Kbar2^(5/2)/Kbar0^(3/2) is imported from 4d_2_5pn.tex, not re-derived."];
  Print["  Runtime isolation: no peer-engine import/subprocess/Get, no report/source/note/_scratch reads, no writes."]
);

printProvenance[];

baseline = buildResiduals[baseConfig];
Print["BASELINE RESIDUAL VECTOR (Wolfram):"];
Print["  " <> vectorText[baseline]];
If[! AllTrue[Values[baseline], TrueQ[# === 0] &], fail["baseline residual vector is not all zero"]];

Print["MUTATED RESIDUAL VECTORS (Wolfram):"];
caughtMatrix = <||>;
Do[
  name = mutation["name"];
  residuals = buildResiduals[mutation["config"]];
  caught = firedLabels[residuals];
  caughtMatrix[name] = caught;
  Print["  " <> name <> ":"];
  Print["    residuals: " <> vectorText[residuals]];
  Print["    caught_by: [" <> StringRiffle[caught, ", "] <> "]"],
  {mutation, mutations}
];

Print["EXPECTED CAUGHT-BY MATRIX (Wolfram):"];
Do[
  name = mutation["name"];
  expected = expectedCaughtBy[name];
  actual = caughtMatrix[name];
  Print["  " <> name <> ": expected=[" <> StringRiffle[expected, ", "] <>
    "] actual=[" <> StringRiffle[actual, ", "] <> "]"];
  If[! TrueQ[actual === expected], fail["caught-by mismatch for " <> name]];
  If[Length[actual] == 0, fail["mutation was not caught: " <> name]],
  {mutation, mutations}
];

Print["NO-FLOAT GUARD: PASS"];
Print["BASELINE ALL-ZERO: PASS"];
Print["MUTATION PROBE: PASS"];
Print["WOLFRAM_MATCHBACK: PASS"];
Exit[0];
