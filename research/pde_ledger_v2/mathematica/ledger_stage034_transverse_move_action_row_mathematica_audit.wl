(* Ledger stage034 Mathematica audit: moving-throat transverse action row.

   Standalone, print-only, assert-zero, exact, and file-I/O-free.
   The stage003/pathA_36 transverse row is instantiated as tagged provenance;
   stage034 earns only the finite-profile moving coupling.

   Independent route: this engine constructs a Fourier-space field on a
   NullSpace-derived transverse basis, extracts quadratic forms with native
   CoefficientArrays, verifies positivity with PositiveDefiniteMatrixQ, and
   solves the two polarization dispersion equations natively.  It does not
   port the SymPy engine's explicit three-component Hessian/eigenvalue route.
   Tooth-local ablation uses LEDGER_STAGE034_MUTATION.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE034_MUTATION";
activeMutation = Environment[mutationEnvironment];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

verdictToken = "TRANSVERSE_MOVE_ACTION_ROW";
noInconsistency = "none";
manifestDigestExpected = "4343bd60cd974f653a0a8ac2eeced6c7aca15b1831c81be20bd358f449c454af";

toothOrder = {
  "ACTION_KINETIC",
  "ACTION_COUPLING",
  "ACTION_STABILITY",
  "G0_DAMAGE",
  "LEDGER_READY_ROW",
  "FIELD_IDENTITY_UNITS",
  "GUARD_IMPORT_VS_EARN",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS",
  "UNITS_RESTORED",
  "VERDICT_REDERIVATION",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "ACTION_KINETIC" ->
    "rho_br/2 -> rho_br/3 in the Fourier action before CoefficientArrays",
  "ACTION_COUPLING" ->
    "q_T=lambda_T*tau_d -> lambda_T in the differentiated source",
  "ACTION_STABILITY" ->
    "rho_br -> -rho_br before PositiveDefiniteMatrixQ",
  "G0_DAMAGE" ->
    "absorb the moving row into the parsed active F_flux ledger row",
  "LEDGER_READY_ROW" ->
    "make the local moving coupling quadratic in u_T",
  "FIELD_IDENTITY_UNITS" ->
    "[b_T]=1 -> L^-1 in the formal unit-rescaling object",
  "GUARD_IMPORT_VS_EARN" ->
    "misclassify imported u_T_kinetic as earned by stage034",
  "TARGET_BLINDNESS" ->
    "inject electric A_E into the live action dependency graph",
  "DUAL_ENGINE_TERMS" ->
    "drop q_T_relation from the computed term inventory",
  "UNITS_RESTORED" ->
    "[q_T]=M T^-1 -> M T^-2 in the whole-density rescaling firewall",
  "VERDICT_REDERIVATION" ->
    "make the verdict pipeline's Fourier kinetic form non-PD",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "drop one deferred source tooth from the canonical partition"
|>;

raise[name_] := Throw[name, "ledgerStage034Failure"];

assertExact[name_, expression_] := Module[{reals},
  reals = Cases[Unevaluated[expression], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FIRST_FAILURE=", name];
    Print["FAIL  ", name, ": machine-real atom(s): ", InputForm[reals]];
    raise[name]
  ]
];

cleanZero[expression_] :=
  FullSimplify[expression] /. ConditionalExpression[0, _] -> 0;

expectZero[name_, residual_, evidence_: None] := Module[{clean},
  assertExact[name, residual];
  clean = cleanZero[residual];
  assertExact[name, clean];
  If[TrueQ[clean === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FIRST_FAILURE=", name];
    If[activeMutation === name, Print["FIRED_AT_OWN_ASSERT=", name]];
    Print["FAIL  ", name, ": residual = ", InputForm[clean]];
    If[evidence =!= None, Print["      evidence = ", InputForm[evidence]]];
    raise[name]
  ]
];

expectBool[name_, condition_, evidence_: None] :=
  expectZero[name, If[TrueQ[condition], 0, 1], evidence];

section[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

uComponents = {u1, u2, u3};
uDots = {ud1, ud2, ud3};
bComponents = {b1, b2, b3};
velocities = {v1, v2, v3};

actionBuild[
    kineticFactor_: 1/2, rhoSign_: 1, muSign_: 1,
    nonlinearCoupling_: False, qDefinition_: Automatic,
    electricLeak_: False] := Module[
  {kinetic, gradient, coupling, density, qDef, linearForm},
  qDef = If[qDefinition === Automatic, lambdaT tauD, qDefinition];
  kinetic = Expand[kineticFactor rhoSign rhoBr Total[uDots^2]];
  gradient = Expand[-muSign muR Total[bComponents^2]/2];
  linearForm = velocities.uComponents;
  coupling = qT s etaA linearForm;
  If[TrueQ[nonlinearCoupling], coupling = Expand[coupling u1]];
  density = Expand[kinetic + gradient + coupling];
  If[TrueQ[electricLeak], density = Expand[density + aElectric u1]];
  <|
    "Density" -> density,
    "Kinetic" -> kinetic,
    "Gradient" -> gradient,
    "Coupling" -> coupling,
    "QDefinition" -> qDef
  |>
];

quadraticMatrix[polynomial_, variables_List] := Module[{arrays},
  arrays = CoefficientArrays[Expand[polynomial], variables];
  2 Normal[arrays[[3]]]
];

fourierData[
    kineticFactor_: 1/2, rhoSign_: 1, muSign_: 1] := Module[
  {direction, basis, amplitude, amplitudeDot, curlAmplitude,
   fourierKinetic, fourierGradientEnergy, kineticMatrix,
   stiffnessMatrix, waveOperator, dispersionRoots, transverseChecks},
  direction = {0, 0, 1};
  basis = NullSpace[{direction}];
  amplitude = a1 basis[[1]] + a2 basis[[2]];
  amplitudeDot = d1 basis[[1]] + d2 basis[[2]];
  curlAmplitude = Cross[k direction, amplitude];
  fourierKinetic =
    Expand[kineticFactor rhoSign rhoBr Total[amplitudeDot^2]];
  fourierGradientEnergy =
    Expand[muSign muR Total[curlAmplitude^2]/2];
  kineticMatrix = FullSimplify[
    quadraticMatrix[fourierKinetic, {d1, d2}]];
  stiffnessMatrix = FullSimplify[
    quadraticMatrix[fourierGradientEnergy, {a1, a2}]];
  waveOperator = FullSimplify[
    stiffnessMatrix - omegaSquared kineticMatrix];
  dispersionRoots = Table[
    FullSimplify[
      omegaSquared /. First@Solve[
        waveOperator[[index, index]] == 0, omegaSquared, Reals]],
    {index, 2}];
  transverseChecks = FullSimplify[direction.#] & /@ basis;
  <|
    "Basis" -> basis,
    "TransverseChecks" -> transverseChecks,
    "KineticMatrix" -> kineticMatrix,
    "StiffnessMatrix" -> stiffnessMatrix,
    "WaveOperator" -> waveOperator,
    "DispersionRoots" -> dispersionRoots
  |>
];

stabilityObject[rhoSign_: 1, muSign_: 1] := Module[
  {data, witnessRules, kineticWitness, stiffnessWitness,
   rootsWitness, kineticPD, stiffnessPD, dispersionPositive},
  data = fourierData[1/2, rhoSign, muSign];
  witnessRules = {rhoBr -> 2, muR -> 3, k -> 5};
  kineticWitness = data["KineticMatrix"] /. witnessRules;
  stiffnessWitness = data["StiffnessMatrix"] /. witnessRules;
  rootsWitness = FullSimplify[data["DispersionRoots"] /. witnessRules];
  kineticPD = TrueQ[PositiveDefiniteMatrixQ[kineticWitness]];
  stiffnessPD = TrueQ[PositiveDefiniteMatrixQ[stiffnessWitness]];
  dispersionPositive = AllTrue[rootsWitness, TrueQ[# > 0] &];
  Join[data, <|
    "KineticWitness" -> kineticWitness,
    "StiffnessWitness" -> stiffnessWitness,
    "RootsWitness" -> rootsWitness,
    "KineticPD" -> kineticPD,
    "StiffnessPD" -> stiffnessPD,
    "DispersionPositive" -> dispersionPositive,
    "Stable" -> (kineticPD && stiffnessPD && dispersionPositive)
  |>]
];

g0RowTranscript = StringRiffle[{
  "bulk_scalar|rho,theta retained scalar action",
  "electric_scalar|localized H/h scalar row",
  "longitudinal|u_L retained scalar-longitudinal row",
  "drain|active drain throughput Gamma_0",
  "return_F_flux|remote return momentum ledger ACTIVE_DEFERRED",
  "wall_r_B|holonomic r_B wall and passive reaction",
  "geon|held-out geon rest constant",
  "zero_bulk_scalar_couplings|r_BH,r_B^2H^2,Hrho,Hdelta_rho,Hdt_theta,Hgrad_theta=0",
  "zero_source_modulation|delta_J_m and neighbor response=0",
  "zero_scalar_transverse_mixing|r_Bu_T,Hu_T,u_Lu_T,two-gradient mixing=0",
  "zero_cross_kinetic|dt_u_L_dt_h and Berry rows=0",
  "zero_scalar_masses|u_L^2,h^2,higher gradients and nonlinearities=0",
  "zero_bulk_modulus|independent B(divu)^2=0",
  "zero_brane_phase|theta_B and brane-phase drain=0",
  "zero_wall_dynamics|bending,anchoring,surface storage,dissipation=0",
  "zero_drain_response|dynamic drain and return responses=0",
  "zero_direct_drain_sources|direct h,u_L drain sources=0",
  "zero_geon_derivatives|field-dependent geon derivatives=0",
  "zero_viscosity|viscosity,drag,no-slip,permeability,phase-jump=0",
  "zero_other_branches|E4,E5,E1 and mixture terms=0",
  "prohibited_ancestry|Maxwell/gauge fields,point sources,native current law,Coulomb prior absent"
}, "\n"];

parseG0Rows[transcript_String] := Module[{pairs, keys},
  pairs = StringSplit[#, "|", 2] & /@
    Select[StringSplit[StringTrim[transcript], "\n"], StringLength[#] > 0 &];
  keys = pairs[[All, 1]];
  If[Length[keys] =!= Length[DeleteDuplicates[keys]],
    raise["G0_DAMAGE"]];
  Association[(#[[1]] -> #[[2]]) & /@ pairs]
];

g0DamageObject[absorbFlux_: False] := Module[
  {before, amendment, after, preexisting, changed, overlap,
   fluxUntouched, newRows, internal},
  before = parseG0Rows[g0RowTranscript];
  amendment = <|
    "delta_transverse_action" -> "imported pathA_36 u_T row",
    "delta_moving_coupling" -> "q_T sum_i s_i eta_a V_i.u_T"
  |>;
  If[TrueQ[absorbFlux],
    AssociateTo[amendment,
      "return_F_flux" -> "SILENTLY_ABSORBED_INTO_CONSERVATIVE_ROW"]];
  after = Join[before, amendment];
  preexisting = Keys[before];
  changed = Sort@Select[preexisting, after[#] =!= before[#] &];
  overlap = Sort@Intersection[preexisting, Keys[amendment]];
  fluxUntouched = after["return_F_flux"] === before["return_F_flux"];
  newRows = Sort@Complement[Keys[after], preexisting];
  internal = If[
    changed === {} && overlap === {} && TrueQ[fluxUntouched],
    noInconsistency, "g0-damage"];
  <|
    "ParsedRowCount" -> Length[before],
    "ChangedPreexisting" -> changed,
    "AmendmentOverlap" -> overlap,
    "FluxUntouched" -> fluxUntouched,
    "NewRows" -> newRows,
    "InternalInconsistency" -> internal
  |>
];

ledgerReadyObject[nonlinear_: False] := Module[
  {build, source, rules, degree, local, sourceLinear,
   momentum, curlResponse, variational, wellFormed},
  build = actionBuild[1/2, 1, 1, nonlinear];
  source = D[build["Coupling"], #] & /@ uComponents;
  rules = CoefficientRules[Expand[build["Coupling"]], uComponents];
  degree = Max[Total /@ (First /@ rules)];
  local = FreeQ[
    build["Density"], Alternatives[Integrate, Sum, Inactive, Inverse]];
  sourceLinear = FreeQ[source, Alternatives @@ uComponents];
  momentum = D[build["Density"], #] & /@ uDots;
  curlResponse = D[build["Density"], #] & /@ bComponents;
  variational =
    momentum === rhoBr uDots &&
    curlResponse === -muR bComponents &&
    AnyTrue[source, !TrueQ[# === 0] &];
  wellFormed =
    TrueQ[local] && degree === 1 && TrueQ[sourceLinear] && TrueQ[variational];
  <|
    "Local" -> local,
    "CouplingDegree" -> degree,
    "SourceLinear" -> sourceLinear,
    "Variational" -> variational,
    "Responses" -> <|
      "Momentum" -> momentum,
      "CurlResponse" -> curlResponse,
      "FieldSource" -> source|>,
    "WellFormed" -> wellFormed
  |>
];

unitScaleRules[badB_: False, badQ_: False] := Join[
  {
    rhoBr -> rhoBr massScale lengthScale^-3,
    muR -> muR massScale lengthScale^-1 timeScale^-2,
    qT -> qT massScale timeScale^If[TrueQ[badQ], -2, -1],
    lambdaT -> lambdaT massScale timeScale^-1,
    tauD -> tauD,
    s -> s,
    etaA -> etaA lengthScale^-3
  },
  Thread[uComponents -> uComponents lengthScale],
  Thread[uDots -> uDots lengthScale timeScale^-1],
  Thread[bComponents ->
    bComponents If[TrueQ[badB], lengthScale^-1, 1]],
  Thread[velocities -> velocities lengthScale timeScale^-1]
];

homogeneousUnderScalingQ[expression_, rules_, targetScale_] :=
  TrueQ[FullSimplify[
    Expand[(expression /. rules) - targetScale expression] == 0,
    lengthScale > 0 && timeScale > 0 && massScale > 0]];

dimensionUnderScaling[expression_, rules_] := Module[{scaled},
  scaled = Expand[expression /. rules];
  FullSimplify[
      # D[scaled, #]/scaled,
      lengthScale > 0 && timeScale > 0 && massScale > 0] & /@
    {lengthScale, timeScale, massScale}
];

dimensionObject[badB_: False, badQ_: False] := Module[
  {rules, build, densityScale, actionScale, identity, termChecks,
   actionCheck, densityDimension, actionDimension},
  rules = unitScaleRules[badB, badQ];
  build = actionBuild[];
  densityScale = massScale lengthScale^-1 timeScale^-2;
  actionScale = massScale lengthScale^2 timeScale^-1;
  identity = <|
    "u_T" -> {1, 0, 0},
    "u_dot_T" -> {1, -1, 0},
    "curl_u_T" -> If[TrueQ[badB], {-1, 0, 0}, {0, 0, 0}],
    "rho_br" -> {-3, 0, 1},
    "mu_R" -> {-1, -2, 1},
    "q_T" -> If[TrueQ[badQ], {0, -2, 1}, {0, -1, 1}],
    "eta_a" -> {-3, 0, 0},
    "b_T" -> If[TrueQ[badB], {-1, 0, 0}, {0, 0, 0}]
  |>;
  termChecks = Association[
    (# -> homogeneousUnderScalingQ[build[#], rules, densityScale]) & /@
      {"Kinetic", "Gradient", "Coupling", "Density"}];
  actionCheck =
    homogeneousUnderScalingQ[
      build["Density"] measureSymbol,
      Append[rules,
        measureSymbol -> measureSymbol timeScale lengthScale^3],
      actionScale];
  densityDimension = dimensionUnderScaling[build["Density"], rules];
  actionDimension = dimensionUnderScaling[
    build["Density"] measureSymbol,
    Append[rules,
      measureSymbol -> measureSymbol timeScale lengthScale^3]];
  <|
    "FieldIdentity" -> identity,
    "TermChecks" -> termChecks,
    "ActionCheck" -> actionCheck,
    "DensityDimension" -> densityDimension,
    "ActionDimension" -> actionDimension
  |>
];

computedTermInventory[includeQRelation_: True] := Module[
  {build, terms, groups = {}, containsAny},
  build = actionBuild[];
  terms = If[Head[Expand[build["Density"]]] === Plus,
    List @@ Expand[build["Density"]], {Expand[build["Density"]]}];
  containsAny[term_, variables_] :=
    AnyTrue[variables, !FreeQ[term, #] &];
  Do[
    Which[
      containsAny[term, uDots], AppendTo[groups, "u_T_kinetic"],
      containsAny[term, bComponents], AppendTo[groups, "u_T_gradient"],
      containsAny[term, uComponents], AppendTo[groups, "moving_coupling"],
      True, AppendTo[groups, "unclassified"]
    ],
    {term, terms}];
  groups = Join[groups, {"c_gamma_relation", "transverse_constraint"}];
  If[TrueQ[includeQRelation], AppendTo[groups, "q_T_relation"]];
  Sort@DeleteDuplicates[groups]
];

expectedTermInventory = Sort@{
  "c_gamma_relation",
  "moving_coupling",
  "q_T_relation",
  "transverse_constraint",
  "u_T_gradient",
  "u_T_kinetic"
};

accountingObject[overclaim_: False] := Module[
  {tags, earned, imported},
  tags = <|
    "u_T_kinetic" -> "IMPORTED_STAGE003_PATHA36",
    "u_T_gradient" -> "IMPORTED_STAGE003_PATHA36",
    "c_gamma_relation" -> "IMPORTED_STAGE003_PATHA36",
    "moving_coupling" -> "EARNED_STAGE034"
  |>;
  If[TrueQ[overclaim], tags["u_T_kinetic"] = "EARNED_STAGE034"];
  earned = Sort@Keys@Select[tags, # === "EARNED_STAGE034" &];
  imported = Sort@Keys@Select[tags, # === "IMPORTED_STAGE003_PATHA36" &];
  <|"EarnedNewTerms" -> earned, "ImportedTerms" -> imported, "Tags" -> tags|>
];

dependencyObject[electricLeak_: False] := Module[
  {build, liveSymbols, forbidden, allowed, forbiddenLive, unknown},
  build = actionBuild[1/2, 1, 1, False, Automatic, electricLeak];
  liveSymbols = DeleteDuplicates@Cases[
    build["Density"],
    symbol_Symbol /; Context[symbol] === "Global`",
    Infinity];
  forbidden = {aElectric, qElectric, gElectric};
  allowed = Join[
    {rhoBr, muR, qT, s, etaA},
    uComponents, uDots, bComponents, velocities];
  forbiddenLive = Sort@Intersection[liveSymbols, forbidden];
  unknown = Sort@Complement[liveSymbols, allowed];
  <|
    "ForbiddenLive" -> forbiddenLive,
    "UnknownLive" -> unknown,
    "SourceSideOnly" -> (forbiddenLive === {} && unknown === {})
  |>
];

deriveVerdict[stability_Association, g0_Association, row_Association] := Module[
  {sectors = {}, internal, token},
  If[!TrueQ[stability["Stable"]], AppendTo[sectors, "stability-fail"]];
  If[g0["InternalInconsistency"] =!= noInconsistency,
    AppendTo[sectors, "g0-damage"]];
  If[!TrueQ[row["WellFormed"]], AppendTo[sectors, "row-malformed"]];
  internal = If[sectors === {}, noInconsistency, StringRiffle[sectors, "+"]];
  token = Which[
    MemberQ[sectors, "stability-fail"], "ROW_UNSTABLE",
    MemberQ[sectors, "g0-damage"], "G0_DAMAGED",
    sectors =!= {}, "AMENDMENT_INCONSISTENT",
    True, verdictToken
  ];
  {token, internal}
];

(* Exact source-build order: all 35 source teeth, no wildcard families. *)
sourceToothIDs = {
  "SOURCE_TRANSLATION_CONTINUITY",
  "SOURCE_NOT_IMPORTED",
  "SOURCE_BASIS",
  "PARITY_RW",
  "PARITY_PW",
  "PARITY_ROTATION",
  "PARITY_TIME_REVERSAL",
  "FIELD_IDENTITY_UNITS",
  "ACTION_KINETIC",
  "ACTION_COUPLING",
  "ACTION_STABILITY",
  "G0_DAMAGE",
  "ROUTE_INDEPENDENCE",
  "BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER",
  "BOOST_COMMON_VELOCITY",
  "DIRECT_SOURCE",
  "DIRECT_PROJECTOR",
  "DIRECT_EXCHANGE_SIGN",
  "DIRECT_FALLOFF",
  "DIRECT_VELOCITY_ORDER",
  "COMPARE_COMPUTED",
  "DELTA_RATIO",
  "CONE_RATIO",
  "QMAG_R1",
  "UNITS_RESTORED",
  "ACTIVE_FLUX_CAVEAT",
  "HOOK_LORENTZ",
  "LEDGER_READY_ROW",
  "TRUTH_TOTALITY",
  "TRUTH_PRECEDENCE",
  "LANDING_OWNERSHIP",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS"
};

sourceManifest = {
  {"SOURCE_TRANSLATION_CONTINUITY", "SCOPED_OUT", "STAGE035_V2"},
  {"SOURCE_NOT_IMPORTED", "SCOPED_OUT", "STAGE035_V2"},
  {"SOURCE_BASIS", "SCOPED_OUT", "STAGE035_V2"},
  {"PARITY_RW", "SCOPED_OUT", "STAGE035_V2"},
  {"PARITY_PW", "SCOPED_OUT", "STAGE035_V2"},
  {"PARITY_ROTATION", "SCOPED_OUT", "STAGE035_V2"},
  {"PARITY_TIME_REVERSAL", "SCOPED_OUT", "STAGE035_V2"},
  {"FIELD_IDENTITY_UNITS", "REPLACED_BY_STRONGER", "STAGE034_DIMENSION_OBJECT"},
  {"ACTION_KINETIC", "REPLACED_BY_STRONGER", "STAGE034_HESSIAN_DISPERSION"},
  {"ACTION_COUPLING", "REPLACED_BY_STRONGER", "STAGE034_DIFFERENTIATED_SOURCE"},
  {"ACTION_STABILITY", "REPLACED_BY_STRONGER", "STAGE034_PD_HESSIANS"},
  {"G0_DAMAGE", "REPLACED_BY_STRONGER", "STAGE034_PARSED_G0_DIFF"},
  {"ROUTE_INDEPENDENCE", "SCOPED_OUT", "STAGE037_V4"},
  {"BOOST_PROJECTOR", "SCOPED_OUT", "STAGE036_V3"},
  {"BOOST_GENERAL_VELOCITIES", "SCOPED_OUT", "STAGE036_V3"},
  {"BOOST_NEXT_ORDER", "SCOPED_OUT", "STAGE036_V3"},
  {"BOOST_COMMON_VELOCITY", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_SOURCE", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_PROJECTOR", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_EXCHANGE_SIGN", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_FALLOFF", "SCOPED_OUT", "STAGE037_V4"},
  {"DIRECT_VELOCITY_ORDER", "SCOPED_OUT", "STAGE037_V4"},
  {"COMPARE_COMPUTED", "SCOPED_OUT", "STAGE037_V4"},
  {"DELTA_RATIO", "SCOPED_OUT", "STAGE037_V4"},
  {"CONE_RATIO", "SCOPED_OUT", "STAGE037_V4"},
  {"QMAG_R1", "SCOPED_OUT", "STAGE037_V4"},
  {"UNITS_RESTORED", "REPLACED_BY_STRONGER", "STAGE034_WHOLE_DENSITY_FIREWALL"},
  {"ACTIVE_FLUX_CAVEAT", "SCOPED_OUT", "STAGE038_V5"},
  {"HOOK_LORENTZ", "SCOPED_OUT", "STAGE038_V5"},
  {"LEDGER_READY_ROW", "REPLACED_BY_STRONGER", "STAGE034_LOCAL_VARIATIONAL_ROW"},
  {"TRUTH_TOTALITY", "SCOPED_OUT", "STAGE038_V5"},
  {"TRUTH_PRECEDENCE", "SCOPED_OUT", "STAGE038_V5"},
  {"LANDING_OWNERSHIP", "SCOPED_OUT", "STAGE038_V5"},
  {"TARGET_BLINDNESS", "PRESERVED", "STAGE034_SOURCE_SIDE_DEPENDENCIES"},
  {"DUAL_ENGINE_TERMS", "REPLACED_BY_STRONGER", "STAGE034_CANONICAL_TERM_INVENTORY"}
};

expectedManifestCounts = KeySort@<|
  "PRESERVED" -> 1,
  "REPLACED_BY_STRONGER" -> 8,
  "SCOPED_OUT" -> 26
|>;

lexicographicCodeLess[left_List, right_List] := Module[{limit, index},
  limit = Min[Length[left], Length[right]];
  index = SelectFirst[
    Range[limit], left[[#]] =!= right[[#]] &, Missing["NotFound"]];
  If[MissingQ[index], Length[left] < Length[right],
    left[[index]] < right[[index]]]
];

canonicalManifestText[manifest_] := Module[{rows},
  rows = StringRiffle[#, "|"] & /@ manifest;
  StringRiffle[
    Sort[rows,
      lexicographicCodeLess[
        ToCharacterCode[#1], ToCharacterCode[#2]] &],
    "\n"]
];

manifestSHA256[manifest_] :=
  IntegerString[Hash[canonicalManifestText[manifest], "SHA256"], 16, 64];

ok = Catch[
  If[activeMutation =!= "" && !MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print["ledger_stage034_transverse_move_action_row Mathematica audit"];
  Print["ROUTE=NullSpace transverse basis + Fourier CoefficientArrays + PositiveDefiniteMatrixQ dispersion"];
  Print["FILE_IO=none; CROSS_ENGINE_COMPARE=none"];
  Print["PROGRESS MATHEMATICA STAGE034 START"];
  If[activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print["MUTATED_PRIMITIVE=", ablationDescriptions[activeMutation]]
  ];

  section[
    "Imported pathA_36 transverse row: native Fourier forms and two polarizations"];
  kineticFactor = If[activeMutation === "ACTION_KINETIC", 1/3, 1/2];
  kineticBuild = actionBuild[kineticFactor];
  kineticData = fourierData[kineticFactor];
  kineticCoefficient = Coefficient[kineticBuild["Kinetic"], ud1^2];
  gradientCoefficient = Coefficient[kineticBuild["Gradient"], b1^2];
  kineticOK =
    TrueQ[FullSimplify[kineticCoefficient == rhoBr/2]] &&
    TrueQ[FullSimplify[gradientCoefficient == -muR/2]] &&
    kineticData["KineticMatrix"] === rhoBr IdentityMatrix[2] &&
    kineticData["StiffnessMatrix"] === muR k^2 IdentityMatrix[2] &&
    And @@ (TrueQ[FullSimplify[# == (muR/rhoBr) k^2]] & /@
      kineticData["DispersionRoots"]) &&
    Length[kineticData["Basis"]] === 2 &&
    kineticData["TransverseChecks"] === {0, 0};
  expectBool["ACTION_KINETIC", kineticOK, <|
    "KineticCoefficient" -> kineticCoefficient,
    "GradientCoefficient" -> gradientCoefficient,
    "Dispersion" -> kineticData["DispersionRoots"],
    "cGammaSquared" -> muR/rhoBr|>];
  Print[
    "      imported provenance=stage003/pathA_36; omega^2=(mu_R/rho_br) k^2; polarizations=2"];
  Print["PROGRESS MATHEMATICA ACTION_KINETIC COMPLETE"];

  section["Magnetism-new finite-profile moving coupling"];
  qDefinition = If[
    activeMutation === "ACTION_COUPLING", lambdaT, lambdaT tauD];
  couplingBuild = actionBuild[1/2, 1, 1, False, qDefinition];
  sourceVector =
    FullSimplify[
      D[couplingBuild["Coupling"] /. qT -> couplingBuild["QDefinition"], #]
        & /@ uComponents];
  expectedSource = lambdaT tauD s etaA velocities;
  normalizedProfile =
    Exp[-radialCoordinate^2/mouthRadius^2]/
      (Pi^(3/2) mouthRadius^3);
  profileIntegral = FullSimplify[
    Integrate[
      4 Pi radialCoordinate^2 normalizedProfile,
      {radialCoordinate, 0, Infinity},
      Assumptions -> mouthRadius > 0],
    mouthRadius > 0];
  couplingOK =
    sourceVector === expectedSource &&
    TrueQ[FullSimplify[
      couplingBuild["QDefinition"] == lambdaT tauD]] &&
    profileIntegral === 1 &&
    And @@ (!FreeQ[couplingBuild["Coupling"], #] & /@
      Join[{qT, s, etaA}, velocities]);
  expectBool["ACTION_COUPLING", couplingOK, <|
    "dLduT" -> sourceVector,
    "QDefinition" -> couplingBuild["QDefinition"],
    "ProfileIntegral" -> profileIntegral|>];
  Print[
    "      earned_new=moving_coupling; q_T=lambda_T*tau_d; integral eta_a d^3x=1"];
  Print["PROGRESS MATHEMATICA ACTION_COUPLING COMPLETE"];

  section["Fourier dispersion and native positive-definite forms"];
  productionStability = stabilityObject[
    If[activeMutation === "ACTION_STABILITY", -1, 1], 1];
  ghostControl = stabilityObject[-1, 1];
  tachyonControl = stabilityObject[1, -1];
  stabilityOK =
    TrueQ[productionStability["Stable"]] &&
    !TrueQ[ghostControl["Stable"]] &&
    !TrueQ[tachyonControl["Stable"]] &&
    productionStability["KineticMatrix"] === rhoBr IdentityMatrix[2] &&
    productionStability["StiffnessMatrix"] ===
      muR k^2 IdentityMatrix[2];
  expectBool["ACTION_STABILITY", stabilityOK, <|
    "KineticWitness" -> productionStability["KineticWitness"],
    "StiffnessWitness" -> productionStability["StiffnessWitness"],
    "RootsWitness" -> productionStability["RootsWitness"],
    "GhostControl" -> ghostControl["Stable"],
    "TachyonControl" -> tachyonControl["Stable"]|>];
  Print[
    "      PositiveDefiniteMatrixQ kinetic=True; gradient=True; no_ghost; no_tachyon"];
  Print["PROGRESS MATHEMATICA ACTION_STABILITY COMPLETE"];

  section["Parsed G0 row diff and active-flux preservation"];
  damage = g0DamageObject[activeMutation === "G0_DAMAGE"];
  expectBool["G0_DAMAGE",
    damage["ParsedRowCount"] === 21 &&
    damage["ChangedPreexisting"] === {} &&
    damage["AmendmentOverlap"] === {} &&
    TrueQ[damage["FluxUntouched"]] &&
    damage["NewRows"] ===
      {"delta_moving_coupling", "delta_transverse_action"} &&
    damage["InternalInconsistency"] === noInconsistency,
    damage];
  Print[
    "      scalar/drain/return_F_flux/wall_r_B/geon/declared-zero rows unchanged"];
  Print["      active F_flux=untouched,deferred(V-5/Part-VII)"];

  section["Local variational G0+delta row"];
  ledger = ledgerReadyObject[activeMutation === "LEDGER_READY_ROW"];
  expectBool["LEDGER_READY_ROW",
    TrueQ[ledger["WellFormed"]] &&
    ledger["CouplingDegree"] === 1 &&
    TrueQ[ledger["SourceLinear"]] &&
    TrueQ[ledger["Local"]] &&
    TrueQ[ledger["Variational"]],
    ledger];
  Print[
    "      one local density; linear current.field source; imported+new pieces variational"];

  section["Field identity and restored units via formal rescaling"];
  fieldDimensions = dimensionObject[
    activeMutation === "FIELD_IDENTITY_UNITS", False];
  expectedFieldIdentity = <|
    "u_T" -> {1, 0, 0},
    "u_dot_T" -> {1, -1, 0},
    "curl_u_T" -> {0, 0, 0},
    "rho_br" -> {-3, 0, 1},
    "mu_R" -> {-1, -2, 1},
    "q_T" -> {0, -1, 1},
    "eta_a" -> {-3, 0, 0},
    "b_T" -> {0, 0, 0}
  |>;
  expectBool["FIELD_IDENTITY_UNITS",
    fieldDimensions["FieldIdentity"] === expectedFieldIdentity &&
    AllTrue[Values[fieldDimensions["TermChecks"]], TrueQ] &&
    TrueQ[fieldDimensions["ActionCheck"]] &&
    fieldDimensions["DensityDimension"] === {-1, -2, 1} &&
    fieldDimensions["ActionDimension"] === {2, -1, 1},
    fieldDimensions];
  Print[
    "      [u_T]=L; [u_dot_T]=L/T; [curl u_T]=[b_T]=1; [S]=E*T"];

  section["Provenance accounting and build-global guards"];
  accounting = accountingObject[
    activeMutation === "GUARD_IMPORT_VS_EARN"];
  expectBool["GUARD_IMPORT_VS_EARN",
    accounting["EarnedNewTerms"] === {"moving_coupling"} &&
    accounting["ImportedTerms"] ===
      {"c_gamma_relation", "u_T_gradient", "u_T_kinetic"},
    accounting];
  Print[
    "      earned_new_terms={moving_coupling}; imported=stage003/pathA_36"];

  dependencies = dependencyObject[
    activeMutation === "TARGET_BLINDNESS"];
  expectBool["TARGET_BLINDNESS",
    TrueQ[dependencies["SourceSideOnly"]] &&
    dependencies["ForbiddenLive"] === {} &&
    dependencies["UnknownLive"] === {},
    dependencies];
  Print[
    "      no A_E; no electric q/g knob; no downstream sign/landing token"];

  inventory = computedTermInventory[
    activeMutation =!= "DUAL_ENGINE_TERMS"];
  expectBool["DUAL_ENGINE_TERMS",
    inventory === expectedTermInventory &&
    !MemberQ[inventory, "unclassified"],
    <|"Computed" -> inventory, "Required" -> expectedTermInventory|>];
  Print["      TERM_INVENTORY=", StringRiffle[inventory, ","]];

  restored = dimensionObject[
    False, activeMutation === "UNITS_RESTORED"];
  expectBool["UNITS_RESTORED",
    AllTrue[Values[restored["TermChecks"]], TrueQ] &&
    TrueQ[restored["ActionCheck"]] &&
    restored["FieldIdentity"]["q_T"] === {0, -1, 1} &&
    restored["DensityDimension"] === {-1, -2, 1} &&
    restored["ActionDimension"] === {2, -1, 1},
    restored];
  Print[
    "      whole action density=[E L^-3]; dt d^3x density=[E T]"];

  section["Computed verdict re-derivation"];
  verdictStability = stabilityObject[
    If[activeMutation === "VERDICT_REDERIVATION", -1, 1], 1];
  {liveVerdict, liveInternal} = deriveVerdict[
    verdictStability, g0DamageObject[], ledgerReadyObject[]];
  {namedAlternative, alternativeInternal} = deriveVerdict[
    stabilityObject[-1, 1], g0DamageObject[], ledgerReadyObject[]];
  expectBool["VERDICT_REDERIVATION",
    liveVerdict === verdictToken &&
    liveInternal === noInconsistency &&
    namedAlternative === "ROW_UNSTABLE" &&
    alternativeInternal === "stability-fail",
    <|
      "Rederived" -> liveVerdict,
      "InternalInconsistency" -> liveInternal,
      "NamedNegativeControl" -> namedAlternative,
      "NegativeInternal" -> alternativeInternal|>];

  section["Canonical source-to-stage predicate manifest"];
  manifestTest = If[
    activeMutation === "SOURCE_TO_STAGE_MANIFEST",
    Most[sourceManifest], sourceManifest];
  identifiers = manifestTest[[All, 1]];
  partitionCounts = KeySort@Counts[manifestTest[[All, 2]]];
  deferred = Cases[
    manifestTest, {identifier_, "SCOPED_OUT", _} :> identifier];
  expectedDeferred = Cases[
    sourceManifest, {identifier_, "SCOPED_OUT", _} :> identifier];
  manifestDigest = manifestSHA256[manifestTest];
  manifestOK =
    identifiers === sourceToothIDs &&
    Length[identifiers] === Length[DeleteDuplicates[identifiers]] === 35 &&
    partitionCounts === expectedManifestCounts &&
    deferred === expectedDeferred &&
    Length[deferred] === 26 &&
    AllTrue[
      manifestTest[[All, 3]],
      StringStartsQ[
        #,
        Alternatives[
          "STAGE034_", "STAGE035_", "STAGE036_", "STAGE037_", "STAGE038_"]]
        &] &&
    manifestDigest === manifestDigestExpected;
  expectBool["SOURCE_TO_STAGE_MANIFEST", manifestOK, <|
    "Entries" -> Length[manifestTest],
    "Partition" -> partitionCounts,
    "Deferred" -> deferred,
    "Digest" -> manifestDigest|>];
  Print[
    "      entries=", Length[manifestTest],
    "; partition=", partitionCounts,
    "; deferred=", Length[deferred],
    "; digest=", manifestDigest];

  Print[""];
  Print[
    "ACTION_ROW=S_T+move:int dt d^3x [rho_br/2 |u_dot_T|^2 - mu_R/2 |curl u_T|^2 + q_T sum_i s_i eta_a V_i.u_T]"];
  Print["TRANSVERSE_CONSTRAINT=div u_T=0; POLARIZATIONS=2"];
  Print["C_GAMMA_SQUARED=mu_R/rho_br"];
  Print["Q_T=lambda_T*tau_d"];
  Print["PROVENANCE=stage003/pathA_36 imported-not-earned"];
  Print["internal_inconsistency=", liveInternal];
  Print["VERDICT_TOKEN: ", liveVerdict];

  If[activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage034Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[TrueQ[ok],
  Print[
    "OVERALL PASS: Mathematica independently reached TRANSVERSE_MOVE_ACTION_ROW"];
  Exit[0],
  Print["OVERALL FAIL: Mathematica stage034 audit did not close"];
  Exit[1]
]
