(* Independent Wolfram Language reader and static gate for reduction-registry-v1. *)

ClearAll[registryFail, require, inputFormString, qidSetString];

registryFail[msg_] := (
  Print["WL_REGISTRY_FAIL: ", msg];
  Exit[1];
);

require[condition_, msg_] := If[! TrueQ[condition], registryFail[msg]];

inputFormString[value_] := ToString[value, InputForm];

qidSetString[qids_List] := "{" <> StringRiffle[qids, ", "] <> "}";

scriptPath = ExpandFileName[$InputFileName];
defaultRegistryDirectory = ExpandFileName[
  FileNameJoin[{DirectoryName[scriptPath], "..", "reduction"}]
];
(* `math -script f.wl` leaves $ScriptCommandLine EMPTY in this environment, so
   Rest[] would fire Rest::norest and trip the usage gate before any registry is
   read.  Drop the script name only when one is actually present. *)
scriptArguments = If[Length[$ScriptCommandLine] > 0, Rest[$ScriptCommandLine], {}];

registryDirectory = Which[
  scriptArguments === {},
    defaultRegistryDirectory,
  Length[scriptArguments] === 2 && scriptArguments[[1]] === "--registry-dir",
    ExpandFileName[scriptArguments[[2]]],
  True,
    registryFail[
      "usage: math -script registry_dimensional_gate.wl " <>
      "[--registry-dir DIRECTORY]"
    ]
];

ClearAll[loadYAMLDocument];
loadYAMLDocument[path_String] := Module[{document},
  require[FileExistsQ[path], "registry document not found: " <> path];
  document = Quiet[Check[Import[path, {"YAML", "Data"}], $Failed]];
  require[document =!= $Failed, "could not parse YAML document: " <> path];
  require[AssociationQ[document], "YAML document is not a mapping: " <> path];
  document
];

quantitiesPath = FileNameJoin[{registryDirectory, "quantities.yaml"}];
relationsPath = FileNameJoin[{registryDirectory, "relations.yaml"}];
quantitiesDocument = loadYAMLDocument[quantitiesPath];
relationsDocument = loadYAMLDocument[relationsPath];

require[
  Lookup[quantitiesDocument, "schema_version", Missing["KeyAbsent"]] ===
    "reduction-registry-v1",
  "unsupported quantities schema_version"
];
require[
  Lookup[relationsDocument, "schema_version", Missing["KeyAbsent"]] ===
    "reduction-registry-v1",
  "unsupported relations schema_version"
];
require[
  Lookup[quantitiesDocument, "document_kind", Missing["KeyAbsent"]] ===
    "quantities",
  "quantities.yaml has the wrong document_kind"
];
require[
  Lookup[relationsDocument, "document_kind", Missing["KeyAbsent"]] ===
    "relations",
  "relations.yaml has the wrong document_kind"
];
require[
  Lookup[relationsDocument, "expression_language", Missing["KeyAbsent"]] ===
    "prefix-v1",
  "unsupported relation expression language"
];

dimensionConvention = Lookup[
  quantitiesDocument, "dimension_convention", Missing["KeyAbsent"]
];
require[AssociationQ[dimensionConvention], "missing dimension_convention mapping"];
require[
  Lookup[dimensionConvention, "name", Missing["KeyAbsent"]] ===
    "LTM-exponent-vector-v1",
  "unsupported dimension convention"
];
require[
  Lookup[dimensionConvention, "ordered_bases", Missing["KeyAbsent"]] ===
    {"L", "T", "M"},
  "dimension ordered_bases must be [L, T, M]"
];

activeRegime = Lookup[quantitiesDocument, "active_regime", Missing["KeyAbsent"]];
require[StringQ[activeRegime], "missing active_regime"];

quantities = Lookup[quantitiesDocument, "quantities", Missing["KeyAbsent"]];
relations = Lookup[relationsDocument, "relations", Missing["KeyAbsent"]];
require[ListQ[quantities] && AllTrue[quantities, AssociationQ],
  "quantities must be a list of mappings"];
require[ListQ[relations] && AllTrue[relations, AssociationQ],
  "relations must be a list of mappings"];

qids = Lookup[quantities, "qid", Missing["KeyAbsent"]];
require[AllTrue[qids, StringQ], "every quantity must have a string qid"];
require[DuplicateFreeQ[qids], "quantity qids are not unique"];
quantityByQID = Association @ Map[Function[q, q["qid"] -> q], quantities];

ClearAll[validateQuantity];
validateQuantity[quantity_Association] := Module[
  {qid, dimension, exponents, state, axis, regime},
  qid = Lookup[quantity, "qid", Missing["KeyAbsent"]];
  dimension = Lookup[quantity, "dimension", Missing["KeyAbsent"]];
  state = Lookup[quantity, "state", Missing["KeyAbsent"]];
  axis = Lookup[quantity, "counting_axis", Missing["KeyAbsent"]];
  regime = Lookup[quantity, "regime", Missing["KeyAbsent"]];
  require[AssociationQ[dimension], "missing dimension mapping for " <> qid];
  require[
    Lookup[dimension, "convention", Missing["KeyAbsent"]] ===
      "LTM-exponent-vector-v1",
    "dimension convention mismatch for " <> qid
  ];
  exponents = Lookup[dimension, "exponents", Missing["KeyAbsent"]];
  require[
    ListQ[exponents] && Length[exponents] === 3 && AllTrue[exponents, IntegerQ],
    "dimension exponents must be three exact integers for " <> qid
  ];
  require[MemberQ[{"live", "retired"}, state], "invalid state for " <> qid];
  require[
    MemberQ[{"continuous-model", "convention-orbit", "discrete-structural"}, axis],
    "invalid counting_axis for " <> qid
  ];
  require[ListQ[regime] && AllTrue[regime, StringQ], "invalid regime for " <> qid];
];
Scan[validateQuantity, quantities];

dimensionsByQID = Association @ Map[
  Function[q, q["qid"] -> q["dimension"]["exponents"]],
  quantities
];
qidSymbols = AssociationThread[qids, Table[Unique["registryQ$"], {Length[qids]}]];

ClearAll[parsePrefix];
parsePrefix[node_] := Module[{operator, arguments, id, p, q},
  If[IntegerQ[node], Return[node]];
  require[ListQ[node] && Length[node] >= 1,
    "invalid prefix-v1 node: " <> inputFormString[node]];
  operator = First[node];
  arguments = Rest[node];
  Switch[operator,
    "Q",
      require[Length[arguments] === 1 && StringQ[arguments[[1]]],
        "Q requires one string QID"];
      id = arguments[[1]];
      require[KeyExistsQ[qidSymbols, id], "unknown QID in expression: " <> id];
      qidSymbols[id],
    "Rat",
      require[Length[arguments] === 2 && AllTrue[arguments, IntegerQ],
        "Rat requires two integer arguments"];
      {p, q} = arguments;
      require[q =!= 0, "Rat denominator is zero"];
      Rational[p, q],
    "Add",
      require[Length[arguments] >= 2, "Add requires at least two arguments"];
      Plus @@ (parsePrefix /@ arguments),
    "Mul",
      require[Length[arguments] >= 2, "Mul requires at least two arguments"];
      Times @@ (parsePrefix /@ arguments),
    "Sub",
      require[Length[arguments] === 2, "Sub requires two arguments"];
      parsePrefix[arguments[[1]]] - parsePrefix[arguments[[2]]],
    "Div",
      require[Length[arguments] === 2, "Div requires two arguments"];
      parsePrefix[arguments[[1]]]/parsePrefix[arguments[[2]]],
    "Pow",
      require[Length[arguments] === 2, "Pow requires two arguments"];
      parsePrefix[arguments[[1]]]^parsePrefix[arguments[[2]]],
    "Neg",
      require[Length[arguments] === 1, "Neg requires one argument"];
      -parsePrefix[arguments[[1]]],
    "Sqrt",
      require[Length[arguments] === 1, "Sqrt requires one argument"];
      Sqrt[parsePrefix[arguments[[1]]]],
    _,
      registryFail["unsupported prefix-v1 operator: " <> inputFormString[operator]]
  ]
];

relationIDs = Lookup[relations, "relation_id", Missing["KeyAbsent"]];
require[AllTrue[relationIDs, StringQ], "every relation must have a string relation_id"];
require[DuplicateFreeQ[relationIDs], "relation ids are not unique"];

ClearAll[qidsIn];
qidsIn[tree_] := DeleteDuplicates[Cases[tree, {"Q", id_String} :> id, Infinity]];

ClearAll[auditRelationDataflow];
auditRelationDataflow[relation_Association] := Module[
  {relationID, residual, output, declaredInputs, rhs, actualInputs, issues = {}, parsed},
  relationID = Lookup[relation, "relation_id", Missing["KeyAbsent"]];
  residual = Lookup[relation, "residual", Missing["KeyAbsent"]];
  output = Lookup[relation, "designated_output", Missing["KeyAbsent"]];
  declaredInputs = Lookup[relation, "input_qids", Missing["KeyAbsent"]];

  require[residual =!= Null && ListQ[residual],
    relationID <> " has no executable prefix-v1 residual"];
  require[StringQ[output] && KeyExistsQ[quantityByQID, output],
    relationID <> " has an unknown designated_output"];
  require[ListQ[declaredInputs] && AllTrue[declaredInputs, StringQ],
    relationID <> " input_qids must be a string list"];
  require[AllTrue[declaredInputs, KeyExistsQ[quantityByQID, #] &],
    relationID <> " input_qids contains an unknown QID"];
  require[StringQ[Lookup[relation, "provenance_status", Missing["KeyAbsent"]]],
    relationID <> " has no provenance_status"];
  require[ListQ[Lookup[relation, "regime", Missing["KeyAbsent"]]] &&
      AllTrue[relation["regime"], StringQ],
    relationID <> " regime must be a string list"];
  require[ListQ[Lookup[relation, "applied_functions", Missing["KeyAbsent"]]],
    relationID <> " applied_functions must be a list"];

  parsed = parsePrefix[residual];
  If[!(Length[residual] === 3 && residual[[1]] === "Sub" &&
      MatchQ[residual[[2]], {"Q", _String}]),
    AppendTo[issues, "residual is not explicit [Sub,[Q,output],rhs] transport"]
  ];

  If[issues === {},
    rhs = residual[[3]];
    actualInputs = qidsIn[rhs];
    If[residual[[2, 2]] =!= output,
      AppendTo[issues, "designated_output is not the left QID"]
    ];
    If[Count[residual, {"Q", output}, Infinity] =!= 1,
      AppendTo[issues, "designated_output is not isolated to the left leaf"]
    ];
    If[! DuplicateFreeQ[declaredInputs] || Sort[declaredInputs] =!= Sort[actualInputs],
      AppendTo[issues,
        "declared inputs " <> qidSetString[declaredInputs] <>
        " differ from RHS QIDs " <> qidSetString[actualInputs]]
    ];
  ];

  <|
    "relation_id" -> relationID,
    "parsed_residual" -> parsed,
    "issues" -> issues
  |>
];

dataflowAudits = auditRelationDataflow /@ relations;
parsedResiduals = AssociationThread[
  relationIDs, Lookup[dataflowAudits, "parsed_residual"]
];
dataflowMatchCount = Count[Lookup[dataflowAudits, "issues"], {}];
dataflowMismatchCount = Length[relations] - dataflowMatchCount;

designatedOutputs = Lookup[relations, "designated_output", Missing["KeyAbsent"]];
require[DuplicateFreeQ[designatedOutputs],
  "more than one relation has the same designated_output"];

dimensionAuditTag = Unique["dimensionAuditTag$"];
ClearAll[dimensionIssue, dimensionOf];
dimensionIssue[status_String, detail_String] := Throw[
  <|"status" -> status, "detail" -> detail|>,
  dimensionAuditTag
];

dimensionOf[node_] := Module[{operator, arguments, childDimensions, p, q},
  If[IntegerQ[node], Return[{0, 0, 0}]];
  If[! ListQ[node] || Length[node] < 1,
    dimensionIssue["UNDETERMINED", "invalid prefix node"]
  ];
  operator = First[node];
  arguments = Rest[node];
  Switch[operator,
    "Q",
      If[Length[arguments] =!= 1 || ! StringQ[arguments[[1]]] ||
          ! KeyExistsQ[dimensionsByQID, arguments[[1]]],
        dimensionIssue["UNDETERMINED", "invalid or unknown Q leaf"]
      ];
      dimensionsByQID[arguments[[1]]],
    "Rat",
      If[Length[arguments] =!= 2 || ! AllTrue[arguments, IntegerQ],
        dimensionIssue["UNDETERMINED", "invalid Rat node"]
      ];
      {p, q} = arguments;
      If[q === 0, dimensionIssue["UNDETERMINED", "zero Rat denominator"]];
      {0, 0, 0},
    "Add" | "Sub",
      childDimensions = dimensionOf /@ arguments;
      If[! SameQ @@ childDimensions,
        dimensionIssue[
          "INHOMOGENEOUS",
          operator <> " term dimensions are " <> inputFormString[childDimensions]
        ]
      ];
      First[childDimensions],
    "Mul",
      Total[dimensionOf /@ arguments],
    "Div",
      dimensionOf[arguments[[1]]] - dimensionOf[arguments[[2]]],
    "Pow",
      If[! IntegerQ[arguments[[2]]],
        dimensionIssue["UNDETERMINED", "Pow exponent is not a bare integer"]
      ];
      arguments[[2]] dimensionOf[arguments[[1]]],
    "Neg",
      dimensionOf[arguments[[1]]],
    "Sqrt",
      dimensionOf[arguments[[1]]]/2,
    _,
      dimensionIssue[
        "UNDETERMINED", "unsupported operator " <> inputFormString[operator]
      ]
  ]
];

ClearAll[auditRelationDimension];
auditRelationDimension[relation_Association] := Module[{result},
  result = Catch[
    <|
      "status" -> "HOMOGENEOUS",
      "dimension" -> dimensionOf[relation["residual"]]
    |>,
    dimensionAuditTag
  ];
  Join[<|"relation_id" -> relation["relation_id"]|>, result]
];

dimensionAudits = auditRelationDimension /@ relations;
homogeneousCount = Count[Lookup[dimensionAudits, "status"], "HOMOGENEOUS"];
inhomogeneousCount = Count[Lookup[dimensionAudits, "status"], "INHOMOGENEOUS"];
undeterminedCount = Count[Lookup[dimensionAudits, "status"], "UNDETERMINED"];
homogeneousDimensionStrings = Cases[
  dimensionAudits,
  audit_ /; audit["status"] === "HOMOGENEOUS" :>
    audit["relation_id"] <> "->" <> inputFormString[audit["dimension"]]
];

Print[
  "DIMENSIONAL_HOMOGENEITY: HOMOGENEOUS=", homogeneousCount,
  "; INHOMOGENEOUS=", inhomogeneousCount,
  "; UNDETERMINED=", undeterminedCount,
  "; TERM_DIMENSIONS={", StringRiffle[homogeneousDimensionStrings, ", "], "}"
];
Scan[
  Function[audit,
    If[audit["status"] =!= "HOMOGENEOUS",
      Print[
        "DIMENSIONAL_DETAIL: ", audit["relation_id"], " ", audit["status"],
        " - ", audit["detail"]
      ]
    ]
  ],
  dimensionAudits
];

Print[
  "RELATION_DECLARATIONS: MATCH=", dataflowMatchCount,
  "; MISMATCH=", dataflowMismatchCount,
  "; RELATIONS=", qidSetString[relationIDs]
];
Scan[
  Function[audit,
    Scan[
      Function[issue,
        Print["RELATION_DETAIL: ", audit["relation_id"], " - ", issue]
      ],
      audit["issues"]
    ]
  ],
  dataflowAudits
];

If[inhomogeneousCount > 0 || undeterminedCount > 0 || dataflowMismatchCount > 0,
  registryFail[
    "registry rejected by dimensional homogeneity or relation declaration gate"
  ]
];

ClearAll[activeLiveQuantityQ];
activeLiveQuantityQ[qid_String] := Module[{quantity = quantityByQID[qid]},
  quantity["state"] === "live" && MemberQ[quantity["regime"], activeRegime]
];

terminalKinds = {"parameter", "boundary-datum", "control", "discrete-choice"};
preliminaryAdmittedRelations = Select[
  relations,
  Function[relation,
    relation["provenance_status"] === "DERIVED-EXECUTED" &&
    MemberQ[relation["regime"], activeRegime] &&
    relation["applied_functions"] === {} &&
    activeLiveQuantityQ[relation["designated_output"]] &&
    AllTrue[relation["input_qids"], activeLiveQuantityQ]
  ]
];
preliminaryOutputToRelation = Association @ Map[
  Function[relation, relation["designated_output"] -> relation],
  preliminaryAdmittedRelations
];

ClearAll[inputClosureQ];
inputClosureQ[qid_String, stack_List] := Module[{quantity, producer},
  quantity = quantityByQID[qid];
  If[MemberQ[terminalKinds, quantity["kind"]], Return[True]];
  If[quantity["kind"] =!= "intermediate", Return[False]];
  If[MemberQ[stack, qid], Return[False]];
  If[! KeyExistsQ[preliminaryOutputToRelation, qid], Return[False]];
  producer = preliminaryOutputToRelation[qid];
  AllTrue[
    producer["input_qids"],
    inputClosureQ[#, Append[stack, qid]] &
  ]
];

admittedRelations = Select[
  preliminaryAdmittedRelations,
  Function[relation,
    AllTrue[
      relation["input_qids"],
      Function[inputQID,
        inputClosureQ[inputQID, {relation["designated_output"]}]
      ]
    ]
  ]
];
admittedOutputs = Lookup[admittedRelations, "designated_output", {}];
activeContinuousQIDs = Map[
  #1["qid"] &,
  Select[
    quantities,
    #1["state"] === "live" &&
      #1["counting_axis"] === "continuous-model" &&
      MemberQ[#1["regime"], activeRegime] &
  ]
];
residue = Select[activeContinuousQIDs, ! MemberQ[admittedOutputs, #] &];

Print["RESIDUE: SIZE=", Length[residue], "; QIDS=", qidSetString[residue]];
Print["WL_REGISTRY_PASS"];
Exit[0];
