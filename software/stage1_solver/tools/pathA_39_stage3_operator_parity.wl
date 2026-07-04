(* pathA_39 Stage 3 operator-parity gate, Mathematica engine. *)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_39_stage3_operator_parity.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
reportsDir = FileNameJoin[{stage1Root, "reports"}];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
sympyJson = FileNameJoin[{scratchDir, "pathA_39_stage3_operator_parity_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_39_stage3_operator_parity_mathematica.json"}];
path38Yaml = FileNameJoin[{reportsDir, "pathA_38_results.yaml"}];
stage01Report = FileNameJoin[{reportsDir, "pathA_39_scalar_admixture_screen.md"}];
stage2Report = FileNameJoin[{reportsDir, "pathA_39_magnetic_force.md"}];

If[! FileExistsQ[sympyJson], fail["missing SymPy JSON: " <> sympyJson]];
If[! FileExistsQ[path38Yaml], fail["missing pathA_38 YAML: " <> path38Yaml]];
If[! FileExistsQ[stage01Report], fail["missing Stage 0+1 report: " <> stage01Report]];
If[! FileExistsQ[stage2Report], fail["missing Stage 2 report: " <> stage2Report]];

$Assumptions = ell > 0;

path38Text = ReadString[path38Yaml];
stage01Text = ReadString[stage01Report];
stage2Text = ReadString[stage2Report];

If[! StringContainsQ[path38Text, "fields:\n  - eta\n  - n\n  - tau"], fail["pathA_38 fields import mismatch"]];
If[! StringContainsQ[path38Text, "parity_matrix:"], fail["pathA_38 Sigma import missing"]];
If[! StringContainsQ[path38Text, "name: translation_goldstone"], fail["missing f0 mode"]];
If[! StringContainsQ[path38Text, "name: wall_shape_partner"], fail["missing wall partner mode"]];
If[! StringContainsQ[path38Text, "name: relative_lock_partner"], fail["missing relative-lock mode"]];
If[! StringContainsQ[path38Text, "name: tau_source_tilt"], fail["missing tau mode"]];
If[! StringContainsQ[path38Text, "status: ENGINE_AGREE"], fail["pathA_38 engine agreement missing"]];
If[! StringContainsQ[stage01Text, "FAIL_OBSERVABLE_SCALAR_ADMIXTURE"] || ! StringContainsQ[stage01Text, "ENGINE_AGREE"],
  fail["Stage 0+1 import mismatch"]
];
If[! StringContainsQ[stage2Text, "MAGNETIC_FORCE_DERIVED"] || ! StringContainsQ[stage2Text, "ENGINE_AGREE"],
  fail["Stage 2 import mismatch"]
];

verdictCodes = <|
  "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING" -> 1,
  "PASS_CONDITIONAL_ON_NO_OPERATOR_PARITY_MIXING" -> 2,
  "PASS_TO_O(V)_CONTROLLED_O(V2)_DEFERRED" -> 3,
  "WEAK_BREAKING_EPSILON_V" -> 4,
  "ACCIDENTAL_ZERO_SIM_DEFERRED" -> 5,
  "BASIS_INCOMPLETE" -> 6,
  "CLEAN_PARITY_PRESERVING" -> 7,
  "REJECT_ILLEGAL_INPUT" -> 8,
  "REJECT_UNDECLARED_SPURION" -> 9,
  "PARITY_TWO_WAYS_AGREE" -> 10,
  "CLEAN_STATIC_SPLIT" -> 11,
  "A1_IDS_RESTORED" -> 12
|>;
verdictCode[label_] := verdictCodes[label];

fields = {"eta", "n", "tau"};
fieldSigma = <|"eta" -> -1, "n" -> -1, "tau" -> 1|>;
xFactors = {"iVdotk", "omegaVdotk"};
xAdjoint = <|"iVdotk" -> -1, "omegaVdotk" -> 1|>;
coeffTypes = {"C", "B", "K"};
wModule = <|
  "C" -> <|"N_w" -> 0, "adjoint_sign" -> 1, "module" -> "C(w)"|>,
  "B" -> <|"N_w" -> 1, "adjoint_sign" -> -1, "module" -> "A_B=1/2[B D_w + D_w B]"|>,
  "K" -> <|"N_w" -> 2, "adjoint_sign" -> 1, "module" -> "-D_w K D_w"|>
|>;
pfSigns = <|"even" -> 1, "odd" -> -1|>;
expectedRawCountFull = 216;

matrixDefs = {
  <|"name" -> "etaeta", "kind" -> "symmetric", "left" -> "eta", "right" -> "eta", "matrix" -> {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}}|>,
  <|"name" -> "nn", "kind" -> "symmetric", "left" -> "n", "right" -> "n", "matrix" -> {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}|>,
  <|"name" -> "tautau", "kind" -> "symmetric", "left" -> "tau", "right" -> "tau", "matrix" -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}|>,
  <|"name" -> "etan", "kind" -> "symmetric", "left" -> "eta", "right" -> "n", "matrix" -> {{0, 1, 0}, {1, 0, 0}, {0, 0, 0}}|>,
  <|"name" -> "etatau", "kind" -> "symmetric", "left" -> "eta", "right" -> "tau", "matrix" -> {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}}|>,
  <|"name" -> "ntau", "kind" -> "symmetric", "left" -> "n", "right" -> "tau", "matrix" -> {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}|>,
  <|"name" -> "[etan]", "kind" -> "antisymmetric", "left" -> "eta", "right" -> "n", "matrix" -> {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}|>,
  <|"name" -> "[etatau]", "kind" -> "antisymmetric", "left" -> "eta", "right" -> "tau", "matrix" -> {{0, 0, 1}, {0, 0, 0}, {-1, 0, 0}}|>,
  <|"name" -> "[ntau]", "kind" -> "antisymmetric", "left" -> "n", "right" -> "tau", "matrix" -> {{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}|>
};

sigmaAB[m_] := fieldSigma[m["left"]] fieldSigma[m["right"]];
matrixAdjointSign[m_] := If[m["kind"] === "symmetric", 1, -1];
matrixAdjointClass[m_] := If[m["kind"] === "symmetric", "self-adjoint", "anti-self-adjoint"];

makeSlot[x_, a_, m_, coeff_, pf_] := <|
  "X" -> x,
  "a" -> a,
  "Mdef" -> m,
  "M" -> m["name"],
  "coeff_type" -> coeff,
  "p_f" -> pf,
  "id" -> StringJoin[x, "|", ToString[a], "|", m["name"], "|", coeff, "|", pf]
|>;

generateSlots[aValues_] := Module[{slots = {}},
  Do[
    AppendTo[slots, makeSlot[x, a, m, coeff, pf]],
    {x, xFactors}, {a, aValues}, {m, matrixDefs}, {coeff, coeffTypes}, {pf, {"even", "odd"}}
  ];
  slots
];

concreteCoefficientProfile[pf_] := Switch[pf,
  "even", Sech[z]^2,
  "odd", Tanh[z] Sech[z]^2,
  _, fail["unknown coefficient parity: " <> ToString[pf]]
];

representativeDerivativeTokens[coeffType_] := Switch[coeffType,
  "C", {},
  "B", {dWB},
  "K", {dWKLeft, dWKRight},
  _, fail["unknown coeff_type: " <> ToString[coeffType]]
];

productExpr[list_] := If[Length[list] == 0, 1, Times @@ list];

matrixSameQ[left_, right_] := And @@ Flatten[Map[TrueQ[FullSimplify[# == 0, ell > 0]] &, left - right, {2}]];

parityConjugationSign[slot_] := Module[
  {sigmaMat, mat, derivativeTokens, sSym, scalar, operator, substitutions, conjugatedScalar, conjugated},
  sigmaMat = DiagonalMatrix[{-1, -1, 1}];
  mat = slot["Mdef"]["matrix"];
  derivativeTokens = representativeDerivativeTokens[slot["coeff_type"]];
  sSym = Unique["s"];
  scalar = sSym^slot["a"] concreteCoefficientProfile[slot["p_f"]] productExpr[derivativeTokens];
  operator = scalar mat;
  substitutions = Join[{z -> -z, sSym -> -sSym}, (# -> -#) & /@ derivativeTokens];
  conjugatedScalar = FullSimplify[scalar /. substitutions, ell > 0];
  conjugated = conjugatedScalar (sigmaMat . mat . sigmaMat);
  Which[
    matrixSameQ[conjugated, operator], 1,
    matrixSameQ[conjugated, -operator], -1,
    True, fail["literal P O P^-1 is not a parity eigen-operator: " <> slot["id"]]
  ]
];

classifySlot[slot_] := Module[
  {xSign, mSign, wSign, product, nW, pSign, parityRule, parityConj, c0Sign},
  xSign = xAdjoint[slot["X"]];
  mSign = matrixAdjointSign[slot["Mdef"]];
  wSign = wModule[slot["coeff_type"]]["adjoint_sign"];
  product = xSign mSign wSign;
  nW = wModule[slot["coeff_type"]]["N_w"];
  pSign = pfSigns[slot["p_f"]];
  parityRule = (-1)^slot["a"] sigmaAB[slot["Mdef"]] pSign (-1)^nW;
  parityConj = parityConjugationSign[slot];
  c0Sign = (-1)^slot["a"];
  <|
    "id" -> slot["id"],
    "X" -> slot["X"],
    "a" -> slot["a"],
    "M" -> slot["M"],
    "coeff_type" -> slot["coeff_type"],
    "p_f" -> slot["p_f"],
    "sigma_Asigma_B" -> sigmaAB[slot["Mdef"]],
    "N_w" -> nW,
    "M_adjoint_class" -> matrixAdjointClass[slot["Mdef"]],
    "X_adjoint_sign" -> xSign,
    "M_adjoint_sign" -> mSign,
    "w_module_adjoint_sign" -> wSign,
    "adjoint_product" -> product,
    "hermitian_retained" -> TrueQ[product == 1],
    "pi_P_rule" -> parityRule,
    "pi_P_conjugation" -> parityConj,
    "parity_methods_agree" -> TrueQ[parityRule == parityConj],
    "P_w_invariant" -> TrueQ[parityRule == 1],
    "C0_sign" -> c0Sign,
    "ibp_representative" -> slot["id"],
    "ibp_note" -> "singleton canonical operator class"
  |>
];

buildTables[aValues_] := Module[
  {slots, rows, manifestRows, adjointRows, hermitianRows, parityRows, physicalRows},
  slots = generateSlots[aValues];
  rows = classifySlot /@ slots;
  manifestRows = (<|
      "id" -> #["id"],
      "sigma_Asigma_B" -> #["sigma_Asigma_B"],
      "N_w" -> #["N_w"],
      "M_adjoint_class" -> #["M_adjoint_class"],
      "p_f" -> #["p_f"]
    |>) & /@ rows;
  adjointRows = (<|
      "id" -> #["id"],
      "X_adjoint_sign" -> #["X_adjoint_sign"],
      "M_adjoint_sign" -> #["M_adjoint_sign"],
      "w_module_adjoint_sign" -> #["w_module_adjoint_sign"],
      "adjoint_product" -> #["adjoint_product"],
      "hermitian_retained" -> #["hermitian_retained"]
    |>) & /@ rows;
  hermitianRows = Select[rows, #["hermitian_retained"] &];
  parityRows = (<|
      "id" -> #["id"],
      "pi_P_rule" -> #["pi_P_rule"],
      "pi_P_conjugation" -> #["pi_P_conjugation"],
      "parity_methods_agree" -> #["parity_methods_agree"],
      "P_w_invariant" -> #["P_w_invariant"],
      "a_class" -> If[#["a"] == 1 && #["P_w_invariant"], "combined-parity-mixing", "parity-preserving"]
    |>) & /@ hermitianRows;
  physicalRows = Select[hermitianRows, #["P_w_invariant"] &];
  <|
    "slots" -> slots,
    "classified_rows" -> rows,
    "manifest_rows" -> manifestRows,
    "adjoint_rows" -> adjointRows,
    "hermitian_rows" -> hermitianRows,
    "parity_rows" -> parityRows,
    "physical_rows" -> physicalRows
  |>
];

dw[expr_] := D[expr, z]/ell;
moduleOnVector[coeffType_, coeff_, vec_] := Switch[coeffType,
  "C", coeff vec,
  "B", coeff (dw /@ vec) + (1/2) dw[coeff] vec,
  "K", (-dw[coeff dw[#]] &) /@ vec,
  _, fail["unknown coeff_type: " <> coeffType]
];

integrateW[expr_] := Module[
  {sSym, tSym, raw, rules, total = 0, powers, coeff, m, n, left, right},
  sSym = Unique["S"];
  tSym = Unique["T"];
  raw = Expand[(expr ell) /. {Sech[z] -> sSym, Tanh[z] -> tSym}];
  If[TrueQ[raw == 0], Return[0]];
  rules = CoefficientRules[raw, {sSym, tSym}];
  Do[
    powers = rule[[1]];
    coeff = rule[[2]];
    m = powers[[1]];
    n = powers[[2]];
    If[OddQ[n], Continue[]];
    If[m <= 0, fail["non-decaying witness monomial"]];
    left = (n + 1)/2;
    right = m/2;
    total = total + coeff Gamma[left] Gamma[right]/Gamma[left + right],
    {rule, rules}
  ];
  FullSimplify[total, ell > 0]
];

slotById[slots_] := Association[(#["id"] -> #) & /@ slots];

buildWitnesses[physicalRows_, slotsById_] := Module[
  {f0, wall, tau, profiles, profileText, pairNames, pairs, a1Rows, witnessRows = {}, exprs = <||>,
    slot, coeff, mat, finiteValues, staticValues, genericFunctionals, nonzeroPairs, finiteNonzeroPairs,
    burden, witnessStatus, bra, ket, acted, overlap, staticValue, pairName, decisive},
  f0 = {Sech[z]^2/ell, Sech[z]^2/ell, 0};
  wall = {Sech[z] Tanh[z], Sech[z] Tanh[z], 0};
  tau = {0, 0, 1};
  profiles = <|"even" -> Sech[z]^2, "odd" -> Tanh[z] Sech[z]^2|>;
  profileText = <|"even" -> "sech(w/ell)^2", "odd" -> "tanh(w/ell)*sech(w/ell)^2"|>;
  pairNames = {"f0_to_wall_partner", "f0_to_tau", "wall_partner_to_f0"};
  pairs = <|"f0_to_wall_partner" -> {f0, wall}, "f0_to_tau" -> {f0, tau}, "wall_partner_to_f0" -> {wall, f0}|>;
  a1Rows = Select[physicalRows, #["a"] == 1 &];
  Do[
    slot = slotsById[row["id"]];
    coeff = profiles[slot["p_f"]];
    mat = slot["Mdef"]["matrix"];
    finiteValues = <||>;
    staticValues = <||>;
    genericFunctionals = <||>;
    Do[
      {bra, ket} = pairs[pairName];
      acted = mat . moduleOnVector[slot["coeff_type"], coeff, ket];
      overlap = integrateW[bra . acted];
      staticValue = If[slot["X"] === "omegaVdotk", 0, overlap];
      finiteValues[pairName] = overlap;
      staticValues[pairName] = staticValue;
      exprs[StringJoin["witness_", IntegerString[index - 1, 10, 3], "_", pairName, "_finite"]] = overlap;
      exprs[StringJoin["witness_", IntegerString[index - 1, 10, 3], "_", pairName, "_static"]] = staticValue;
      genericFunctionals[pairName] = StringJoin[
        slot["X"], "*s^", ToString[slot["a"]], "*Int[bra=", First[StringSplit[pairName, "_to_"]],
        ", M=", slot["M"], ", module=", slot["coeff_type"], ", coeff_parity=", slot["p_f"],
        ", ket=", Last[StringSplit[pairName, "_to_"]], "]"
      ],
      {pairName, pairNames}
    ];
    nonzeroPairs = Select[pairNames, ! TrueQ[FullSimplify[staticValues[#] == 0, ell > 0]] &];
    finiteNonzeroPairs = Select[pairNames, ! TrueQ[FullSimplify[finiteValues[#] == 0, ell > 0]] &];
    burden = slot["X"] === "iVdotk";
    witnessStatus = Which[
      burden && Length[nonzeroPairs] > 0, "NONZERO_STATIC_WITNESS",
      burden, "STRUCTURAL_ZERO_FOR_IMPORTED_TARGETS",
      True, "OMEGA_STATIC_MODE_NO_WITNESS_BURDEN"
    ];
    AppendTo[witnessRows, <|
      "slot_id" -> row["id"],
      "X" -> slot["X"],
      "a" -> slot["a"],
      "M" -> slot["M"],
      "coeff_type" -> slot["coeff_type"],
      "p_f" -> slot["p_f"],
      "witness_profile" -> profileText[slot["p_f"]],
      "functionals" -> genericFunctionals,
      "finite_frequency_w_overlap_values" -> Association[(# -> ToString[finiteValues[#], InputForm]) & /@ pairNames],
      "static_omega0_values" -> Association[(# -> ToString[staticValues[#], InputForm]) & /@ pairNames],
      "static_nonzero_pairs" -> nonzeroPairs,
      "finite_frequency_nonzero_pairs" -> finiteNonzeroPairs,
      "witness_status" -> witnessStatus,
      "continuum_threshold_flag" -> "not_excited_by_native_discrete_witness"
    |>],
    {index, Length[a1Rows]}, {row, {a1Rows[[index]]}}
  ];
  decisive = Select[witnessRows, #["X"] === "iVdotk" && #["witness_status"] === "NONZERO_STATIC_WITNESS" &];
  <|"a1_witness_rows" -> witnessRows, "decisive_static_witness_rows" -> decisive, "exprs" -> exprs|>
];

neutralOverlapControlExpr[slotsById_] := Module[
  {slot, f0, wall, coeff, acted, value},
  slot = slotsById["iVdotk|0|etaeta|B|odd"];
  f0 = {Sech[z]^2/ell, Sech[z]^2/ell, 0};
  wall = {Sech[z] Tanh[z], Sech[z] Tanh[z], 0};
  coeff = Tanh[z] Sech[z]^2;
  acted = slot["Mdef"]["matrix"] . moduleOnVector[slot["coeff_type"], coeff, wall];
  value = integrateW[f0 . acted];
  <|"slot_id" -> "iVdotk|0|etaeta|B|odd", "profile" -> "tanh(w/ell)*sech(w/ell)^2", "overlap" -> value|>
];

illegalAttemptSpec[illegalInput_] := Module[
  {base},
  base = <|
    "X" -> "iVdotk",
    "a" -> 1,
    "M" -> "etaeta",
    "coeff_type" -> "C",
    "p_f" -> "even",
    "odd_spurions" -> {"s"}
  |>;
  Switch[illegalInput,
    "extra_odd_spurion",
      Join[base, <|
        "attempted_structure" -> "iVdotk*s*q_odd_extra*C_even(w)*etaeta",
        "odd_spurions" -> {"s", "q_odd_extra"}
      |>],
    "V_i n_i",
      Join[base, <|
        "attempted_structure" -> "V_i n_i*s*C_even(w)*etaeta",
        "X" -> "V_i n_i",
        "requires_declared_brane_vector" -> "d_i"
      |>],
    "epsilon_ijk V_i k_j",
      Join[base, <|
        "attempted_structure" -> "epsilon_ijk V_i k_j*s*C_even(w)*etaeta",
        "X" -> "epsilon_ijk V_i k_j",
        "requires_declared_axial_spurion" -> "J_k"
      |>],
    _, fail["unknown illegal control input: " <> ToString[illegalInput]]
  ]
];

grammarRefusal[verdict_, reason_, observation_] := <|
  "accepted" -> False,
  "verdict" -> verdict,
  "reason" -> reason,
  "observation" -> observation
|>;

validateOperatorAttempt[attempt_] := Module[
  {declaredOddSpurions, declaredBraneVectors, declaredAxialSpurions, observations, oddSpurions,
    undeclaredOdd, requiredBraneVector, requiredAxial, matrixDef, slot, generatedIds},
  declaredOddSpurions = {"s"};
  declaredBraneVectors = {};
  declaredAxialSpurions = {};
  observations = <|
    "declared_external_scalars" -> xFactors,
    "declared_odd_spurions" -> declaredOddSpurions,
    "declared_brane_vectors" -> declaredBraneVectors,
    "declared_axial_spurions" -> declaredAxialSpurions,
    "attempted_structure" -> attempt["attempted_structure"]
  |>;

  oddSpurions = Lookup[attempt, "odd_spurions", {}];
  undeclaredOdd = Sort[Complement[oddSpurions, declaredOddSpurions]];
  If[Length[undeclaredOdd] > 0,
    Return[grammarRefusal[
      "REJECT_ILLEGAL_INPUT",
      "undeclared P_w-odd spurion distinct from s",
      Join[observations, <|"undeclared_odd_spurions" -> undeclaredOdd|>]
    ]]
  ];

  requiredBraneVector = Lookup[attempt, "requires_declared_brane_vector", None];
  If[requiredBraneVector =!= None && ! MemberQ[declaredBraneVectors, requiredBraneVector],
    Return[grammarRefusal[
      "REJECT_UNDECLARED_SPURION",
      attempt["X"] <> " requires an undeclared brane vector",
      Join[observations, <|"missing_declared_brane_vector" -> requiredBraneVector|>]
    ]]
  ];

  requiredAxial = Lookup[attempt, "requires_declared_axial_spurion", None];
  If[requiredAxial =!= None && ! MemberQ[declaredAxialSpurions, requiredAxial],
    Return[grammarRefusal[
      "REJECT_UNDECLARED_SPURION",
      attempt["X"] <> " requires an undeclared axial spurion",
      Join[observations, <|"missing_declared_axial_spurion" -> requiredAxial|>]
    ]]
  ];

  If[! MemberQ[xFactors, attempt["X"]],
    Return[grammarRefusal[
      "REJECT_ILLEGAL_INPUT",
      attempt["X"] <> " is not a declared external scalar",
      Join[observations, <|"undeclared_external_scalar" -> attempt["X"]|>]
    ]]
  ];
  If[! MemberQ[{0, 1}, attempt["a"]],
    Return[grammarRefusal[
      "REJECT_ILLEGAL_INPUT",
      "a=" <> ToString[attempt["a"]] <> " is outside the leading O(s^1) grammar",
      Join[observations, <|"invalid_a" -> attempt["a"]|>]
    ]]
  ];
  If[! MemberQ[coeffTypes, attempt["coeff_type"]],
    Return[grammarRefusal[
      "REJECT_ILLEGAL_INPUT",
      attempt["coeff_type"] <> " is not a declared coefficient module",
      Join[observations, <|"invalid_coeff_type" -> attempt["coeff_type"]|>]
    ]]
  ];
  If[! KeyExistsQ[pfSigns, attempt["p_f"]],
    Return[grammarRefusal[
      "REJECT_ILLEGAL_INPUT",
      attempt["p_f"] <> " is not a declared coefficient parity",
      Join[observations, <|"invalid_p_f" -> attempt["p_f"]|>]
    ]]
  ];

  matrixDef = SelectFirst[matrixDefs, #["name"] === attempt["M"] &, Missing["not-found"]];
  If[MissingQ[matrixDef],
    Return[grammarRefusal[
      "REJECT_ILLEGAL_INPUT",
      attempt["M"] <> " is not in the declared field-matrix basis",
      Join[observations, <|"invalid_matrix" -> attempt["M"]|>]
    ]]
  ];

  slot = makeSlot[attempt["X"], attempt["a"], matrixDef, attempt["coeff_type"], attempt["p_f"]];
  generatedIds = (#["id"] &) /@ generateSlots[{0, 1}];
  If[! MemberQ[generatedIds, slot["id"]],
    Return[grammarRefusal[
      "REJECT_ILLEGAL_INPUT",
      slot["id"] <> " was not generated by the canonical grammar",
      Join[observations, <|"candidate_slot_id" -> slot["id"]|>]
    ]]
  ];
  <|"accepted" -> True, "slot" -> slot|>
];

illegalControlResult[name_, illegalInput_] := Module[
  {attempt, validation},
  attempt = illegalAttemptSpec[illegalInput];
  validation = validateOperatorAttempt[attempt];
  If[! validation["accepted"],
    Return[<|
      "name" -> name,
      "status" -> "FIRED",
      "verdict" -> validation["verdict"],
      "verdict_code" -> verdictCode[validation["verdict"]],
      "reason" -> validation["reason"],
      "attempted_structure" -> attempt["attempted_structure"],
      "validation_observation" -> "REFUSED_BY_DECLARED_GRAMMAR",
      "validator_observation" -> validation["observation"]
    |>]
  ];
  <|
    "name" -> name,
    "status" -> "NOT_FIRED",
    "verdict" -> "ILLEGAL_INPUT_ACCEPTED",
    "reason" -> "declared grammar admitted an illegal control structure",
    "attempted_structure" -> attempt["attempted_structure"],
    "accepted_slot_id" -> validation["slot"]["id"]
  |>
];

decideCase[name_, opts___Rule] := Module[
  {allowedA, c0Symmetry, deleteId, illegalInput, vZero, witnessById, tables, generatedIds, expectedIds,
    missing, extra, hermitian, parityOk, c0Ok, active, a0, a1, nonzeroA1, verdict, result},
  allowedA = "allowed_a" /. {opts} /. "allowed_a" -> {0, 1};
  c0Symmetry = "c0_symmetry" /. {opts} /. "c0_symmetry" -> False;
  deleteId = "delete_id" /. {opts} /. "delete_id" -> None;
  illegalInput = "illegal_input" /. {opts} /. "illegal_input" -> None;
  vZero = "v_zero" /. {opts} /. "v_zero" -> False;
  witnessById = "witness_by_id" /. {opts} /. "witness_by_id" -> <||>;

  If[illegalInput =!= None, Return[illegalControlResult[name, illegalInput]]];

  tables = buildTables[allowedA];
  generatedIds = (#["id"] &) /@ tables["classified_rows"];
  expectedIds = generatedIds;
  If[deleteId =!= None, generatedIds = Select[generatedIds, # =!= deleteId &]];
  missing = Sort[Complement[expectedIds, generatedIds]];
  extra = Sort[Complement[generatedIds, expectedIds]];
  If[Length[missing] > 0 || Length[extra] > 0,
    Return[<|"name" -> name, "status" -> "FIRED", "verdict" -> "BASIS_INCOMPLETE",
      "expected_raw_count" -> Length[expectedIds], "generated_raw_count" -> Length[generatedIds],
      "missing_ids" -> missing, "extra_ids" -> extra|>]
  ];
  hermitian = Select[tables["hermitian_rows"], MemberQ[generatedIds, #["id"]] &];
  parityOk = Select[hermitian, #["P_w_invariant"] &];
  c0Ok = Select[parityOk, (! c0Symmetry) || #["C0_sign"] == 1 &];
  active = If[vZero, {}, c0Ok];
  a0 = Select[active, #["a"] == 0 &];
  a1 = Select[active, #["a"] == 1 &];
  nonzeroA1 = Select[a1, KeyExistsQ[witnessById, #["id"]] && witnessById[#["id"]]["witness_status"] === "NONZERO_STATIC_WITNESS" &];
  verdict = Which[
    vZero, "CLEAN_STATIC_SPLIT",
    c0Symmetry && Length[a1] == 0, "PASS_CONDITIONAL_ON_NO_OPERATOR_PARITY_MIXING",
    Length[a1] == 0, "CLEAN_PARITY_PRESERVING",
    Length[nonzeroA1] > 0, "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING",
    True, "ACCIDENTAL_ZERO_SIM_DEFERRED"
  ];
  result = <|
    "name" -> name,
    "status" -> "FIRED",
    "verdict" -> verdict,
    "verdict_code" -> verdictCode[verdict],
    "raw_count" -> Length[generatedIds],
    "hermitian_count" -> Length[hermitian],
    "physical_count" -> Length[c0Ok],
    "active_operator_count" -> Length[active],
    "a0_count" -> Length[a0],
    "a1_count" -> Length[a1],
    "a1_nonzero_static_witness_count" -> Length[nonzeroA1],
    "a1_ids" -> (#["id"] &) /@ a1,
    "a1_nonzero_static_witness_ids" -> (#["id"] &) /@ nonzeroA1
  |>;
  If[c0Symmetry,
    result = Join[result, <|"c0_label_basis" -> If[Length[a1] == 0,
      "observed_post_C0_filter_a1_count_eq_0",
      "observed_post_C0_filter_a1_count_nonzero"
    ]|>]
  ];
  result
];

controlLabel[item_] := ToString[Lookup[item, "verdict", Lookup[item, "status", "UNKNOWN"]]];

tables = buildTables[{0, 1}];
slotsById = slotById[tables["slots"]];
witnesses = buildWitnesses[tables["physical_rows"], slotsById];
witnessById = Association[(#["slot_id"] -> #) & /@ witnesses["a1_witness_rows"]];
neutralControl = neutralOverlapControlExpr[slotsById];
parityMismatches = Select[tables["hermitian_rows"], ! #["parity_methods_agree"] &];
If[Length[parityMismatches] > 0, fail["parity two-method mismatch"]];

main = decideCase["charged_moving_throat", "witness_by_id" -> witnessById];
neutral = decideCase["neutral_moving_throat", "allowed_a" -> {0}, "witness_by_id" -> witnessById];
c0 = decideCase["injected_C0", "c0_symmetry" -> True, "witness_by_id" -> witnessById];
deleted = decideCase["deleted_operator", "delete_id" -> tables["classified_rows"][[1]]["id"], "witness_by_id" -> witnessById];
extraOdd = decideCase["injected_extra_odd_spurion", "illegal_input" -> "extra_odd_spurion"];
smuggledVector = decideCase["smuggled_V_i_n_i", "illegal_input" -> "V_i n_i"];
smuggledAxial = decideCase["smuggled_epsilon_ijk_V_i_k_j", "illegal_input" -> "epsilon_ijk V_i k_j"];
vZero = decideCase["V_equals_zero", "v_zero" -> True, "witness_by_id" -> witnessById];
c0RemovedRestores = <|
  "name" -> "C0_removed_restores_a1",
  "status" -> If[main["a1_count"] > 0 && c0["a1_count"] == 0, "FIRED", "NOT_FIRED"],
  "verdict" -> If[main["a1_count"] > 0 && c0["a1_count"] == 0, "A1_IDS_RESTORED", "RESTORE_FAILED"],
  "without_C0_a1_count" -> main["a1_count"],
  "with_C0_a1_count" -> c0["a1_count"]
|>;
parityTwoWays = <|
  "name" -> "parity_two_ways",
  "status" -> If[Length[parityMismatches] == 0, "FIRED", "NOT_FIRED"],
  "verdict" -> If[Length[parityMismatches] == 0, "PARITY_TWO_WAYS_AGREE", "PARITY_TWO_WAYS_DISAGREE"],
  "checked_retained_count" -> Length[tables["hermitian_rows"]]
|>;
controls = <|
  "neutral_moving_throat" -> Join[neutral, <|"required_overlap_exercised" -> <|
      "slot_id" -> neutralControl["slot_id"], "profile" -> neutralControl["profile"],
      "overlap" -> ToString[neutralControl["overlap"], InputForm],
      "status" -> If[TrueQ[FullSimplify[neutralControl["overlap"] == 0, ell > 0]], "FIRED", "NOT_FIRED"]|>|>],
  "charged_moving_throat" -> main,
  "injected_C0" -> c0,
  "C0_removed_restores_a1" -> c0RemovedRestores,
  "injected_extra_odd_spurion" -> extraOdd,
  "deleted_operator" -> deleted,
  "smuggled_V_i_n_i" -> smuggledVector,
  "smuggled_epsilon_ijk_V_i_k_j" -> smuggledAxial,
  "parity_two_ways" -> parityTwoWays,
  "V_equals_zero" -> vZero
|>;
If[! And @@ ((#["status"] === "FIRED") & /@ Values[controls]), fail["a control did not fire"]];

physicalRows = tables["physical_rows"];
a0Rows = Select[physicalRows, #["a"] == 0 &];
a1Rows = Select[physicalRows, #["a"] == 1 &];
countSummary = <|
  "raw_manifest_count" -> Length[tables["classified_rows"]],
  "expected_raw_manifest_count" -> expectedRawCountFull,
  "hermitian_retained_pre_ibp" -> Length[tables["hermitian_rows"]],
  "hermitian_retained_post_ibp" -> Length[DeleteDuplicates[(#["ibp_representative"] &) /@ tables["hermitian_rows"]]],
  "parity_retained_pre_ibp" -> Length[physicalRows],
  "parity_retained_post_ibp" -> Length[DeleteDuplicates[(#["ibp_representative"] &) /@ physicalRows]],
  "a0_physical_count" -> Length[a0Rows],
  "a1_physical_count" -> Length[a1Rows],
  "a1_static_nonzero_witness_count" -> Length[witnesses["decisive_static_witness_rows"]]
|>;
If[countSummary["raw_manifest_count"] =!= expectedRawCountFull, fail["raw manifest count mismatch"]];

expectedLandingObservation = <|
  "expected" -> "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING",
  "computed" -> main["verdict"],
  "status" -> If[
    main["verdict"] === "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING",
    "OBSERVED_EXPECTED_FAIL",
    "COMPUTED_NONFAIL_EMITTED"
  ]
|>;

actuals = Join[
  <|
    "raw_manifest_count" -> countSummary["raw_manifest_count"],
    "hermitian_retained_post_ibp" -> countSummary["hermitian_retained_post_ibp"],
    "parity_retained_post_ibp" -> countSummary["parity_retained_post_ibp"],
    "a0_physical_count" -> countSummary["a0_physical_count"],
    "a1_physical_count" -> countSummary["a1_physical_count"],
    "a1_static_nonzero_witness_count" -> countSummary["a1_static_nonzero_witness_count"],
    "parity_mismatch_count" -> Length[parityMismatches],
    "neutral_control_overlap" -> neutralControl["overlap"],
    "main_verdict_code" -> verdictCode[main["verdict"]],
    "neutral_verdict_code" -> verdictCode[neutral["verdict"]],
    "c0_verdict_code" -> verdictCode[c0["verdict"]],
    "deleted_verdict_code" -> verdictCode[deleted["verdict"]],
    "extra_odd_verdict_code" -> verdictCode[extraOdd["verdict"]],
    "smuggled_vector_verdict_code" -> verdictCode[smuggledVector["verdict"]],
    "smuggled_axial_verdict_code" -> verdictCode[smuggledAxial["verdict"]],
    "v_zero_verdict_code" -> verdictCode[vZero["verdict"]]
  |>,
  witnesses["exprs"]
];

sympyResults = Import[sympyJson, "RawJSON"];
sympyExprs = sympyResults["engine_agreement"]["mathematica_exprs"];
sympyDigest = sympyResults["engine_agreement"]["sympy_expression_digest"];

assertEngine[name_] := Module[{expectedText, expectedExpr, actual},
  expectedText = sympyExprs[name];
  If[! StringQ[expectedText], fail["missing SymPy expression for " <> name]];
  If[! KeyExistsQ[actuals, name], fail["missing Mathematica actual for " <> name]];
  expectedExpr = ToExpression[expectedText, InputForm];
  actual = actuals[name];
  If[! TrueQ[FullSimplify[actual == expectedExpr, ell > 0]],
    fail["engine disagreement " <> name <> ": Mathematica got " <>
      ToString[actual, InputForm] <> ", SymPy exported " <> expectedText]
  ];
];

engineKeys = Keys[sympyExprs];
Scan[assertEngine, engineKeys];

topLineVerdict = main["verdict"];
If[topLineVerdict =!= sympyResults["top_line_verdict"], fail["top-line verdict disagreement"]];

controlVerdicts = Association @ KeyValueMap[#1 -> controlLabel[#2] &, controls];
agreementPayload = <|
  "top_line_verdict" -> topLineVerdict,
  "main_verdict_code" -> verdictCode[topLineVerdict],
  "counts" -> countSummary,
  "manifest_ids" -> (#["id"] &) /@ tables["classified_rows"],
  "hermitian_ids" -> (#["id"] &) /@ tables["hermitian_rows"],
  "physical_ids" -> (#["id"] &) /@ physicalRows,
  "a0_ids" -> (#["id"] &) /@ a0Rows,
  "a1_ids" -> (#["id"] &) /@ a1Rows,
  "adjoint_table" -> tables["adjoint_rows"],
  "parity_table" -> tables["parity_rows"],
  "witness_status_table" -> ((<|
      "slot_id" -> #["slot_id"],
      "witness_status" -> #["witness_status"],
      "static_nonzero_pairs" -> #["static_nonzero_pairs"],
      "finite_frequency_nonzero_pairs" -> #["finite_frequency_nonzero_pairs"]
    |>) & /@ witnesses["a1_witness_rows"]),
  "control_verdicts" -> controlVerdicts,
  "checked_expression_count" -> Length[engineKeys],
  "expr_digest" -> sympyDigest
|>;

(* RawJSON structural equality is enforced by the SymPy --compare runner after
   this engine exports.  Mathematica has already checked every expression in
   sympyExprs and the top-line verdict independently above. *)

out = <|
  "schema" -> "pathA_39_stage3_operator_parity_mathematica/v1",
  "status" -> "OK",
  "headline" -> topLineVerdict,
  "checked_quantities" -> engineKeys,
  "sympy_expression_digest" -> sympyDigest,
  "expected_landing_observation" -> expectedLandingObservation,
  "agreement_payload" -> agreementPayload,
  "computed" -> <|
    "counts" -> countSummary,
    "control_verdicts" -> controlVerdicts,
    "decisive_witness_count" -> Length[witnesses["decisive_static_witness_rows"]]
  |>
|>;

Export[jsonOut, out, "RawJSON"];
Print["OK pathA_39_stage3_operator_parity_mathematica"];
Print[ExportString[<|"json" -> jsonOut, "headline" -> topLineVerdict|>, "RawJSON"]];
