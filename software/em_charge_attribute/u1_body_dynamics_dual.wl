#!/usr/bin/env wolframscript
(* U1 Phase-A field computation: independent Wolfram engine. *)

ClearAll["Global`*"];

here = DirectoryName[ExpandFileName[$InputFileName]];
root = ExpandFileName[FileNameJoin[{here, "..", ".."}]];
defaultInput = FileNameJoin[{here, "u1_body_dynamics_inputs.yaml"}];
defaultFixtures = FileNameJoin[{here, "u1_body_dynamics_fixtures.yaml"}];
defaultOutput = FileNameJoin[{here, "reports", "u1_body_dynamics_artifacts", "mathematica_phase_a.json"}];
endpoints = {"E1", "E2", "E3", "E4", "E5"};
teeth = {"INPUT_LEDGER", "SOURCE_COMPLETENESS", "DIMENSIONAL",
  "COMOVING_CONTINUITY", "COMOVING_MOMENTUM", "BASE_BALANCE", "TAIL_ODE",
  "ZERO_MODE", "PROJECTOR", "ENDPOINT_RESPONSE", "MOMENT_INTEGRALS",
  "RECONSTRUCTION", "CANONICAL_VARIATION", "CHANNEL_UNIQUENESS",
  "TYPED_DATAFLOW", "PROVENANCE_FORBIDDEN", "ANCESTRY", "NATIVE_PADDING", "PARITY",
  "OUTCOME_REACHABILITY"};

argValue[name_, fallback_] := Module[{p = FirstPosition[$ScriptCommandLine, name]},
  If[MissingQ[p], fallback, $ScriptCommandLine[[p[[1]] + 1]]]
];
require[test_, tooth_, detail_] := If[! TrueQ[test],
  Print["ASSERT_FAIL:" <> tooth <> ":" <> ToString[detail, InputForm]]; Exit[1]
];
inputPath = ExpandFileName[argValue["--input", defaultInput]];
fixturePath = ExpandFileName[argValue["--fixtures", defaultFixtures]];
outputPath = ExpandFileName[argValue["--output", defaultOutput]];
caseName = argValue["--case", "production"];
mutation = argValue["--mutation", "none"];

data = Import[inputPath, "YAML"];
fixtures = Import[fixturePath, "YAML"];
require[AssociationQ[data] && data["schema_version"] === "U1_PHASE_A_FIELD_INPUTS_V2", "INPUT_LEDGER", "production schema"];
require[AssociationQ[fixtures] && fixtures["schema_version"] === "U1_PHASE_A_FIXTURES_V2", "INPUT_LEDGER", "fixture schema"];

pathParts[path_] := Map[If[StringMatchQ[#, DigitCharacter ..], ToExpression[#] + 1, Key[#]] &, StringSplit[path, "."]];
attacks = <||>;
applyOps[ops_List] := Do[
  Switch[op["op"],
    "set", data = ReplacePart[data, pathParts[op["path"]] -> op["value"]],
    "append", data = ReplacePart[data, pathParts[op["path"]] -> Append[Extract[data, pathParts[op["path"]]], op["value"]]],
    "delete_by_id", data = ReplacePart[data, pathParts[op["path"]] -> Select[Extract[data, pathParts[op["path"]]], #["id"] =!= op["id"] &]],
    "derived_attack", attacks[op["target"]] = ToString[op["value"]],
    _, require[False, "INPUT_LEDGER", "unknown fixture operation"]
  ], {op, ops}];
If[caseName =!= "production",
  require[KeyExistsQ[fixtures["outcome_cases"], caseName], "OUTCOME_REACHABILITY", "unknown case"];
  applyOps[fixtures["outcome_cases"][caseName]]
];
If[mutation =!= "none",
  require[MemberQ[teeth, mutation] && KeyExistsQ[fixtures["mutations"], mutation], "INPUT_LEDGER", "unknown mutation"];
  applyOps[fixtures["mutations"][mutation]]
];

externalNames = <||>;
safeSymbol[name_String] := Module[{s = Symbol["u1" <> StringReplace[name, "_" -> "ZZ"]]}, externalNames[s] = name; s];
externalName[s_Symbol] := Lookup[externalNames, s, SymbolName[Unevaluated[s]]];
parseValue[x_] := Which[NumberQ[x], x, StringQ[x] && StringMatchQ[x, NumberString ~~ ("/" ~~ NumberString) ...], ToExpression[x], True, safeSymbol[ToString[x]]];
coeffValues = Association@KeyValueMap[#1 -> parseValue[#2["value"]] &, data["coefficients"]];
coeffSymbols = Association@KeyValueMap[#1 -> safeSymbol[#1] &, data["coefficients"]];
coeffConstraints = Association@KeyValueMap[#1 -> #2["constraint"] &, data["coefficients"]];
coeffDims = Association@KeyValueMap[#1 -> #2["dimensions"] &, data["coefficients"]];
aValue = parseValue[data["geometry"]["a"]["value"]];

(* Structural typed-leaf inspection: tokenize domains/arguments, then inspect the
   resulting argument graph. *)
requiredFields = {"id", "status", "root_type", "domain", "dimensions", "arguments", "symmetry_class", "source"};
records = data["field_records"];
Do[
  require[SubsetQ[Keys[rec], requiredFields], "INPUT_LEDGER", "incomplete " <> ToString[Lookup[rec, "id", "?"]]];
  domainTokens = DeleteDuplicates[StringCases[rec["domain"], (LetterCharacter | "_") ..]];
  argTokens = ToString /@ rec["arguments"];
  require[Intersection[Join[domainTokens, argTokens], {"X", "p", "collective_coordinate_functional", "multipole"}] === {},
    "INPUT_LEDGER", "answer warehouse " <> rec["id"]];
  require[Length[rec["dimensions"]] === 3, "INPUT_LEDGER", "dimension tuple"], {rec, records}];
declaredInputs = Join[records, KeyValueMap[<|"id" -> #1, "status" -> #2["status"], "root_type" -> "ACTION_COEFFICIENT",
    "domain" -> "coefficient", "dimensions" -> #2["dimensions"], "arguments" -> {}, "symmetry_class" -> "R_W_EVEN",
    "source" -> #2["source"], "constraint" -> #2["constraint"]|> &, data["coefficients"]]];

extractFence[text_, label_] := Module[{lines = StringSplit[text, {"\r\n", "\n"}], starts, start, stop},
  starts = Flatten@Position[lines, "```" <> label]; require[Length[starts] >= 1, "SOURCE_COMPLETENESS", "missing fence " <> label];
  start = First[starts]; stop = First@Select[Range[start + 1, Length[lines]], lines[[#]] === "```" &];
  StringRiffle[lines[[start + 1 ;; stop - 1]], "\n"]
];
mathTokens[text_String] := StringCases[text, RegularExpression["[A-Za-z_][A-Za-z_0-9]*|[0-9]+(?:\\.[0-9]+)?|:=|->|!=|==|\\^|\\*|/|\\+|-|\\(|\\)|\\[|\\]|\\||\\.|·|∂|∇|Ω|χ|ρ|λ|ϖ"]];
parsedSourceContains[source_String, probe_String] := Module[{needle = mathTokens[probe], records},
  If[needle === {}, Return[False]]; records = mathTokens /@ StringSplit[source, {"\r\n", "\n"}];
  AnyTrue[records, Length[#] >= Length[needle] && MemberQ[Partition[#, Length[needle], 1], needle] &]
];
policy = data["operative_action_decision"];
decisionRel = policy["source_file"];
decisionPath = FileNameJoin[{root, decisionRel}];
decisionText = Import[decisionPath, "Text"];
require[policy["id"] === "decision_16_retire_brane_polar_field" && policy["status"] === "OPERATIVE" &&
  decisionRel === "software/stage1_solver/decisions/16_retire_brane_polar_field.md",
  "SOURCE_COMPLETENESS", "Decision 16 is not the operative action policy"];
citationMatches = {};
Do[require[parsedSourceContains[decisionText, fragment], "SOURCE_COMPLETENESS", "Decision 16 citation mismatch: " <> fragment];
  AppendTo[citationMatches, fragment], {fragment, policy["required_source_fragments"]}];
expectedIds = Sort[policy["expected_action_term_ids"]]; retiredIds = policy["retired_action_term_ids"];
retiredSymbols = policy["retired_expression_symbols"]; assembledIds = data["action_terms"][[All, "id"]];
require[DuplicateFreeQ[assembledIds], "SOURCE_COMPLETENESS", "duplicate assembled action id"];
Do[exprSymbols = StringCases[term["expression"], RegularExpression["[A-Za-z_][A-Za-z_0-9]*"]];
  intrusion = MemberQ[retiredIds, term["id"]] || Intersection[exprSymbols, retiredSymbols] =!= {};
  If[intrusion,
    require[Lookup[term, "decision_citation", ""] === decisionRel, "SOURCE_COMPLETENESS",
      "retired P term lacks Decision 16 citation: " <> term["id"]];
    require[False, "SOURCE_COMPLETENESS", "retired P term cannot be assembled under Decision 16: " <> term["id"]]],
  {term, data["action_terms"]}];
missingIds = Complement[expectedIds, assembledIds]; unexpectedIds = Complement[assembledIds, expectedIds];
require[missingIds === {}, "SOURCE_COMPLETENESS", "missing P-retired whitelisted terms " <> ToString[missingIds, InputForm]];
require[unexpectedIds === {}, "SOURCE_COMPLETENESS", "unexpected assembled terms " <> ToString[unexpectedIds, InputForm]];

sourceCache = <||>; matchedTerms = {};
Do[
  rel = term["source_file"];
  If[! KeyExistsQ[sourceCache, rel], sourceCache[rel] = Import[FileNameJoin[{root, rel}], "Text"]];
  (* Source fragments are parsed as normalized action records, not accepted from a verdict list. *)
  require[parsedSourceContains[sourceCache[rel], term["source_contains"]],
    "SOURCE_COMPLETENESS", "parsed source mismatch " <> term["id"]];
  AppendTo[matchedTerms, <|"id" -> term["id"], "source_file" -> rel, "source_fragment" -> term["source_contains"]|>],
  {term, data["action_terms"]}];
t0 = extractFence[Import[FileNameJoin[{root, "software/stage1_solver/reports/pathA_24_T0_freeze.md"}], "Text"], "freeze-action"];
g0 = extractFence[Import[FileNameJoin[{root, "software/stage1_solver/reports/pathA_35_G0_freeze.md"}], "Text"], "freeze-action"];
t0Lines = StringTrim /@ StringSplit[t0, "\n"];
startPol = First@FirstPosition[t0Lines, "L_pol ="];
stopPol = First@FirstPosition[t0Lines, x_String /; StringStartsQ[x, "Frozen extended"]];
polarLines = StringTrim[StringTrim[#, {"+", "-", " "}], "."] & /@
  Select[t0Lines[[startPol + 1 ;; stopPol - 1]], StringMatchQ[#, ("+" | "-" | "") ~~ Whitespace ... ~~ "(1/" ~~ ("2" | "4") ~~ ") m rho" ~~ __] &];
require[Length[polarLines] === 3, "SOURCE_COMPLETENESS", "legacy T0 polar manifest changed"];
g0Records = mathTokens /@ StringSplit[g0, "\n"];
g0Definitions = Sort@DeleteDuplicates@Cases[g0Records, rec_List /; Length[rec] > 1 && rec[[2]] === ":=" && (rec[[1]] === "QP" || StringStartsQ[rec[[1]], "L_"]) :> rec[[1]]];
require[Length[g0Definitions] >= 4, "SOURCE_COMPLETENESS", "G0 definitions"];
discoveredG0Records = StringTrim[StringTrim[#, {"+", "-", " "}], "."] & /@
  Select[StringSplit[g0, "\n"], StringContainsQ[#, "QP :="] || StringContainsQ[#, "mu_R Omega_u^a Omega_u^a"] || StringStartsQ[StringTrim[#], "L_Pu :="] &];
retiredG0 = Select[discoveredG0Records, StringContainsQ[#, "L_Pu :="] &];
mandatoryG0 = Select[discoveredG0Records, ! StringContainsQ[#, "L_Pu :="] &];
require[Length[retiredG0] === 1, "SOURCE_COMPLETENESS", "legacy L_Pu retirement record changed"];
assembledFragments = data["action_terms"][[All, "source_contains"]];
Do[expected = mathTokens[fragment];
  require[AnyTrue[assembledFragments, candidate |-> Module[{got = mathTokens[candidate]},
      Length[expected] >= Length[got] && MemberQ[Partition[expected, Length[got], 1], got]]],
    "SOURCE_COMPLETENESS", "missing immutable G0 record " <> fragment], {fragment, mandatoryG0}];
sourceCompleteness = <|"loaded_files" -> Sort[Keys[sourceCache]], "matched_assembled_terms" -> matchedTerms,
  "operative_decision_citation" -> <|"id" -> policy["id"], "status" -> policy["status"],
    "source_file" -> decisionRel, "sha256" -> FileHash[decisionPath, "SHA256", "HexString"],
    "matched_source_fragments" -> citationMatches|>,
  "expected_p_retired_action_ids" -> expectedIds, "retired_action_ids" -> Sort[retiredIds],
  "retired_parameter_rows" -> policy["retired_parameter_rows"],
  "legacy_retired_t0_polar_monomials" -> polarLines, "source_derived_g0_definitions" -> g0Definitions,
  "source_derived_mandatory_g0_records" -> mandatoryG0,
  "source_derived_retired_g0_records" -> retiredG0,
  "assembled_ids" -> Sort[assembledIds]|>;

dadd[args__List] := Total[{args}]; dscale[d_List, n_Integer] := n d;
Ldim = {1, 0, 0}; tinv = {0, -1, 0}; grad = {-1, 0, 0}; nd = coeffDims["rho_inf"]; vel = {1, -1, 0}; targetDim = {-2, -2, 1};
dimRows = {
  "bulk_berry" -> dadd[coeffDims["hbar"], nd, tinv],
  "bulk_flow_kinetic" -> dadd[coeffDims["m_GNLS"], nd, dscale[vel, 2]],
  "quantum_pressure" -> dadd[dscale[coeffDims["hbar"], 2], dscale[coeffDims["m_GNLS"], -1], dscale[nd, -1], dscale[dadd[nd, grad], 2]],
  "bulk_EOS" -> dadd[coeffDims["K_EOS"], dscale[nd, 5]], "wall_double_well" -> coeffDims["aB"],
  "wall_gradient" -> dadd[coeffDims["kappaB"], dscale[grad, 2]], "wall_shear_gate" -> coeffDims["muR4"],
  "brane_shear_kinetic" -> dadd[dscale[coeffDims["ellg"], -1], coeffDims["rhoBr"], dscale[dadd[Ldim, tinv], 2]],
  "brane_shear_gradient" -> dadd[dscale[coeffDims["ellg"], -1], coeffDims["muR"]],
  "uw_kinetic" -> dadd[dscale[coeffDims["ellg"], -1], coeffDims["rhoBr"], dscale[dadd[Ldim, tinv], 2]],
  "uw_gap" -> dadd[dscale[coeffDims["ellg"], -1], coeffDims["rhoBr"], dscale[coeffDims["OmegaW"], 2], dscale[Ldim, 2]],
  "h_kinetic" -> dadd[coeffDims["Mh"], dscale[coeffDims["cE"], -2], dscale[tinv, 2]],
  "h_gradient" -> dadd[coeffDims["Mh"], dscale[grad, 2]]};
Do[require[Last[row] === targetDim, "DIMENSIONAL", First[row] <> ToString[Last[row]]], {row, dimRows}];
momentumDensityDim = dadd[coeffDims["m_GNLS"], nd, vel]; momentumRateDim = dadd[momentumDensityDim, tinv];
require[momentumRateDim === dadd[targetDim, grad], "DIMENSIONAL", "momentum PDE"];
require[dadd[targetDim, {3, 0, 0}] === {1, -2, 1}, "DIMENSIONAL", "surface force"];
dimensionalFirewall = Join[Map[<|"expression" -> First[#], "computed_dimensions_LTM" -> Last[#]|> &, dimRows],
  {<|"expression" -> "momentum_rate_density", "computed_dimensions_LTM" -> momentumRateDim|>,
   <|"expression" -> "control_surface_force", "computed_dimensions_LTM" -> {1, -2, 1}|>}];

Clear[t, x, y, V, ff, nf, vf, gf, pif];
mapSign = data["kinematics"]["coordinate_map"]["displacement_coefficient"];
bodyArgument = x + mapSign V t; compositeField = ff[bodyArgument];
scalarChain = D[compositeField, t]; scalarChainExpanded = mapSign V ff'[bodyArgument];
scalarChainResidual = FullSimplify[scalarChain - scalarChainExpanded];
require[scalarChainResidual === 0 && mapSign === -1, "COMOVING_CONTINUITY", scalarChainResidual];
nativeContinuity = D[nf[bodyArgument], t] + D[nf[bodyArgument] vf[bodyArgument], x];
derivedContinuity = D[nf[bodyArgument] (vf[bodyArgument] + mapSign V), x];
continuityResidual = FullSimplify[nativeContinuity - derivedContinuity];
require[continuityResidual === 0, "COMOVING_CONTINUITY", continuityResidual];
fluxSign = data["kinematics"]["momentum_flux"]["convective_coefficient"];
nativeMomentum = D[gf[bodyArgument], t] + D[pif[bodyArgument], x];
derivedMomentum = D[pif[bodyArgument] + mapSign V gf[bodyArgument], x];
declaredMomentum = D[pif[bodyArgument] + fluxSign V gf[bodyArgument], x];
momentumResidual = FullSimplify[nativeMomentum - derivedMomentum];
declaredFluxResidual = FullSimplify[declaredMomentum - derivedMomentum];
require[momentumResidual === 0, "COMOVING_MOMENTUM", momentumResidual];
require[declaredFluxResidual === 0, "COMOVING_MOMENTUM", declaredFluxResidual];
coMoving = <|"coordinate_substitution" -> data["kinematics"]["coordinate_map"], "composite_field" -> ToString[compositeField, InputForm],
  "scalar_chain_rule" -> ToString[scalarChain, InputForm], "scalar_chain_rule_expanded" -> ToString[scalarChainExpanded, InputForm],
  "scalar_chain_rule_residual" -> ToString[scalarChainResidual, InputForm],
  "continuity_native" -> "partial_t n+div_4(n v)=0", "continuity_comoving" -> "partial_t|y n+div_4[n(v-V)]=0",
  "continuity_residual" -> ToString[continuityResidual, InputForm], "momentum_native" -> "partial_t g_i+partial_j Pi_ij=f_i[action]",
  "momentum_comoving" -> "partial_t|y g_i+partial_j(Pi_ij-V_j g_i)=f_i[action]", "momentum_residual" -> ToString[momentumResidual, InputForm],
  "declared_flux_residual" -> ToString[declaredFluxResidual, InputForm],
  "surface_force" -> "integral_Omega f_i d4y-integral_partialOmega(Pi_ij-V_j g_i)N_j dSigma3", "particle_momentum_imported" -> False|>;

sphereArea[d_] := 2 Pi^(d/2)/Gamma[d/2];
signOf[expr_] := Module[{e = FullSimplify[expr], syms, unresolved},
  If[TrueQ[e === 0], Return[0]]; If[NumericQ[e], Return[Sign[N[e]]]];
  syms = Cases[e, s_Symbol /; KeyExistsQ[externalNames, s] && KeyExistsQ[coeffConstraints, externalName[s]], Infinity] // DeleteDuplicates;
  unresolved = Select[syms, coeffConstraints[externalName[#]] === "real" &];
  If[unresolved =!= {}, Missing["OpenSign"], Sign[N[e /. Thread[syms -> 1]]]]
];

actionPrimitiveNames = {"n", "theta_t", "v2", "n_grad2", "chiB", "chi_grad2", "ud_curl2",
  "f_throat", "f_mix", "g_ell", "u_t2", "u_curl2", "uw_t2", "uw", "h_t2", "h_grad2"};
actionPrimitive = Association@Table[name -> Symbol["u1AP" <> StringReplace[name, "_" -> "ZZ"]], {name, actionPrimitiveNames}];
actionPrimitiveReverse = AssociationThread[Values[actionPrimitive], Keys[actionPrimitive]];
actionDependencyName[s_Symbol] := Lookup[actionPrimitiveReverse, s, externalName[s]];

assembleAction[localData_] := Module[{v, cs, parse, terms},
  v = Association@KeyValueMap[#1 -> parseValue[#2["value"]] &, localData["coefficients"]];
  v["a"] = parseValue[localData["geometry"]["a"]["value"]];
  cs = 5 v["K_EOS"] v["rho_inf"]^4/v["m_GNLS"];
  parse[text_] := Module[{rewritten = text},
    Do[rewritten = StringReplace[rewritten, RegularExpression["\\b" <> name <> "\\b"] -> "(" <> ToString[v[name], InputForm] <> ")"],
      {name, Reverse@SortBy[Keys[v], StringLength]}];
    Do[rewritten = StringReplace[rewritten, RegularExpression["\\b" <> name <> "\\b"] -> SymbolName[actionPrimitive[name]]],
      {name, Reverse@SortBy[actionPrimitiveNames, StringLength]}];
    rewritten = StringReplace[rewritten, RegularExpression["\\bcs2\\b"] -> "(" <> ToString[cs, InputForm] <> ")"];
    ToExpression[rewritten]];
  terms = Association@Table[row["id"] -> parse[row["expression"]], {row, localData["action_terms"]}];
  <|"terms" -> terms, "cs2" -> FullSimplify[cs], "total_expression" -> Total[Values[terms]],
    "dependencies" -> Association@KeyValueMap[#1 -> Sort[DeleteDuplicates[actionDependencyName /@
      Cases[{#2}, s_Symbol /; KeyExistsQ[externalNames, s] || KeyExistsQ[actionPrimitiveReverse, s], Infinity]]] &, terms]|>
];

deriveOperator[localData_] := Module[{v, asm, terms, d, db, r, eps, dn, theta, chi, ut, uwf, hh,
    dnr, thetar, chir, ucurl, hr, drain, rules, quadBy, quad, fields, gradients, curv, gradH, mixed, ch},
  v = Association@KeyValueMap[#1 -> parseValue[#2["value"]] &, localData["coefficients"]]; asm = assembleAction[localData]; terms = asm["terms"];
  v["a"] = parseValue[localData["geometry"]["a"]["value"]];
  d = Length[localData["ambient"]["coordinates"]]; db = Length[localData["ambient"]["brane_coordinates"]];
  Clear[r, eps, dn, theta, chi, ut, uwf, hh, dnr, thetar, chir, ucurl, hr];
  drain = -v["mdot"]/(v["hbar"] v["rho_inf"] sphereArea[d]) r^(1 - d);
  rules = {actionPrimitive["n"] -> v["rho_inf"] + eps dn, actionPrimitive["theta_t"] -> 0,
    actionPrimitive["v2"] -> (v["hbar"]/v["m_GNLS"])^2 (drain + eps thetar)^2,
    actionPrimitive["n_grad2"] -> eps^2 dnr^2, actionPrimitive["chiB"] -> eps chi,
    actionPrimitive["chi_grad2"] -> eps^2 chir^2, actionPrimitive["ud_curl2"] -> 0,
    actionPrimitive["f_throat"] -> 0, actionPrimitive["f_mix"] -> 0, actionPrimitive["g_ell"] -> 1,
    actionPrimitive["u_t2"] -> 0, actionPrimitive["u_curl2"] -> eps^2 ucurl^2,
    actionPrimitive["uw_t2"] -> 0, actionPrimitive["uw"] -> eps uwf,
    actionPrimitive["h_t2"] -> 0, actionPrimitive["h_grad2"] -> eps^2 hr^2};
  quadBy = Association@KeyValueMap[#1 -> FullSimplify[Coefficient[Normal@Series[#2 /. rules, {eps, 0, 2}], eps, 2]] &, terms];
  quad = FullSimplify[Total[Values[quadBy]]]; fields = {dn, theta, chi, ut, uwf, hh}; gradients = {dnr, thetar, chir, ucurl, hr};
  curv = Table[FullSimplify[D[quad, fields[[i]], fields[[j]]]], {i, Length[fields]}, {j, Length[fields]}];
  gradH = Table[FullSimplify[D[quad, gradients[[i]], gradients[[j]]]], {i, Length[gradients]}, {j, Length[gradients]}];
  mixed = Table[FullSimplify[D[quad, fields[[i]], gradients[[j]]]], {i, Length[fields]}, {j, Length[gradients]}];
  ch = {<|"id" -> "density_EOS", "field" -> "delta_n", "dimension" -> d, "ell" -> 0, "gradient" -> gradH[[1, 1]], "curvature" -> curv[[1, 1]]|>,
    <|"id" -> "wall_chiB", "field" -> "delta_chiB", "dimension" -> d, "ell" -> 0, "gradient" -> gradH[[3, 3]], "curvature" -> curv[[3, 3]]|>,
    <|"id" -> "bound_phase", "field" -> "theta", "dimension" -> d, "ell" -> 0, "gradient" -> gradH[[2, 2]], "curvature" -> curv[[2, 2]], "continuity" -> True|>,
    <|"id" -> "brane_shear", "field" -> "u_transverse", "dimension" -> db, "ell" -> 1, "gradient" -> gradH[[4, 4]], "curvature" -> 0, "brane_profile" -> "g_ell(w)"|>,
    <|"id" -> "uw", "field" -> "u_w", "dimension" -> db, "ell" -> 0, "gradient" -> 0, "curvature" -> curv[[5, 5]], "algebraic" -> True|>,
    <|"id" -> "h", "field" -> "h", "dimension" -> d, "ell" -> 0, "gradient" -> gradH[[5, 5]], "curvature" -> curv[[6, 6]], "open_coefficients" -> {"Mh", "cE"}|>};
  <|"channels" -> ch, "quadratic_expression" -> quad, "quadratic_by_term" -> quadBy,
    "field_order" -> (ToString[#, InputForm] & /@ fields), "gradient_order" -> (ToString[#, InputForm] & /@ gradients),
    "curvature_hessian" -> curv, "gradient_hessian" -> gradH, "mixed_hessian" -> mixed, "drain_gradient" -> drain,
    "entries" -> <|"density_gradient" -> gradH[[1, 1]], "density_EOS_curvature" -> curv[[1, 1]],
      "phase_gradient" -> gradH[[2, 2]], "wall_gradient" -> gradH[[3, 3]], "wall_well_curvature" -> curv[[3, 3]],
      "shear_curl" -> gradH[[4, 4]], "uw_curvature" -> curv[[5, 5]], "h_gradient" -> gradH[[5, 5]],
      "drain_density_phase" -> mixed[[1, 2]]|>|>
];

deriveChannels[localData_] := deriveOperator[localData]["channels"];

coupledAnalysis[localData_, op_] := Module[{d, nu, diagDegree, crossDegree},
  d = Length[localData["ambient"]["coordinates"]]; Clear[nu]; diagDegree = -nu - 2; crossDegree = -nu - d;
  <|"density_phase" -> <|"cross_coefficient" -> ToString[op["entries"]["drain_density_phase"], InputForm],
      "diagonal_degree" -> ToString[diagDegree, InputForm], "cross_degree" -> ToString[crossDegree, InputForm],
      "degree_difference" -> ToString[FullSimplify[crossDegree - diagDegree], InputForm], "leading_exponents_shifted" -> ! TrueQ[d > 2],
      "computed_conclusion" -> "SUBLEADING_FOR_D_GT_2"|>,
    "changes_scalar_channel_verdict" -> False|>
];

solveTail[ch_] := Module[{r, ff, d = ch["dimension"], ell = ch["ell"], A = ch["gradient"], B = ch["curvature"],
    radial, odeExpr, ratio, sgn, nu, roots, sol, residual, norm, k, order, lambdaRoots, normValue, deps},
  deps = Sort[DeleteDuplicates[externalName /@ Cases[{A, B}, s_Symbol /; KeyExistsQ[externalNames, s] && KeyExistsQ[coeffConstraints, externalName[s]], Infinity]]];
  If[TrueQ[Lookup[ch, "algebraic", False]], Return[<|"id" -> ch["id"], "field" -> ch["field"], "radial_dimension" -> d,
    "equation" -> ToString[B ff[r] == 0, InputForm], "solved_general_form" -> "f(r)=0 outside the collar", "classification" -> "ALGEBRAIC_GAP",
    "decay_exponent" -> Null, "gap_squared" -> ToString[B, InputForm], "solution_residual" -> "0",
    "zero_mode_norm_integral" -> "0", "zero_mode_norm_value" -> 0., "normalizable" -> True, "dependencies" -> deps|>]];
  radial = ff''[r] + (d - 1) ff'[r]/r - ell (ell + d - 2) ff[r]/r^2; odeExpr = A radial - B ff[r];
  ratio = FullSimplify[B/A]; sgn = signOf[ratio];
  If[MissingQ[sgn], Return[<|"id" -> ch["id"], "field" -> ch["field"], "radial_dimension" -> d,
    "equation" -> ToString[odeExpr == 0, InputForm], "solved_general_form" -> "conditional on open coefficient stratum",
    "classification" -> "UNRESOLVED", "decay_exponent" -> Null, "gap_squared" -> ToString[ratio, InputForm],
    "solution_residual" -> "conditional", "zero_mode_norm_integral" -> "conditional", "zero_mode_norm_value" -> Null,
    "normalizable" -> Null, "dependencies" -> deps|>]];
  If[sgn < 0, k = Sqrt[-ratio]; sol = r^(1 - d/2) BesselJ[ell + d/2 - 1, k r]; residual = FullSimplify[odeExpr /. ff -> Function[{r}, sol]];
    Return[<|"id" -> ch["id"], "field" -> ch["field"], "radial_dimension" -> d, "equation" -> ToString[odeExpr == 0, InputForm],
      "solved_general_form" -> ToString[sol, InputForm], "classification" -> "TACHYONIC", "decay_exponent" -> Null,
      "gap_squared" -> ToString[ratio, InputForm], "solution_residual" -> ToString[residual, InputForm],
      "zero_mode_norm_integral" -> "not used: static operator tachyonic", "zero_mode_norm_value" -> Null,
      "normalizable" -> False, "dependencies" -> deps|>]];
  If[sgn == 0, roots = nu /. Solve[nu (nu - d + 2) - ell (ell + d - 2) == 0, nu];
    nu = Max[Append[Select[roots, TrueQ[# > 0] &], 0]]; sol = r^-nu;
    If[ch["id"] === "density_EOS" && Lookup[attacks, "density_tail_solution", ""] === "add_constant", sol = sol + r];
    residual = FullSimplify[A (D[sol, {r, 2}] + (d - 1) D[sol, r]/r - ell (ell + d - 2) sol/r^2) - B sol];
    require[residual === 0, "TAIL_ODE", ch["id"] <> ToString[residual]];
    norm = Assuming[R > 0, Integrate[r^(d - 1) D[r^-nu, r]^2, {r, R, Infinity}]];
    Return[<|"id" -> ch["id"], "field" -> ch["field"], "radial_dimension" -> d, "equation" -> ToString[odeExpr == 0, InputForm],
      "solved_general_form" -> "C_grow r^ell+C_decay r^-nu", "indicial_equation" -> ToString[nu (nu - d + 2) - ell (ell + d - 2), InputForm],
      "indicial_roots" -> (ToString[#, InputForm] & /@ roots), "classification" -> "POWER_LAW", "decay_exponent" -> ToString[nu, InputForm],
      "gap_squared" -> "0", "solution_residual" -> ToString[residual, InputForm], "zero_mode_norm_integral" -> ToString[norm, InputForm],
      "zero_mode_norm_value" -> If[nu > 0 && FreeQ[norm, DirectedInfinity], N[norm /. R -> 1], Null],
      "normalizable" -> TrueQ[nu > 0 && FreeQ[norm, DirectedInfinity]], "dependencies" -> deps|>]];
  k = Sqrt[ratio]; order = ell + d/2 - 1; sol = r^(1 - d/2) BesselK[order, k r];
  (* Reduce the radial operator to the defining modified-Bessel equation. *)
  Clear[yy]; genericSol = r^(1 - d/2) yy[k r];
  genericResidual = Expand[A (D[genericSol, {r, 2}] + (d - 1) D[genericSol, r]/r - ell (ell + d - 2) genericSol/r^2) - B genericSol];
  besselRule = Derivative[2][yy][k r] -> (((k r)^2 + order^2) yy[k r] - (k r) Derivative[1][yy][k r])/(k r)^2;
  residual = FullSimplify[genericResidual /. besselRule, Assumptions -> r > 0];
  If[ch["id"] === "density_EOS" && Lookup[attacks, "density_tail_solution", ""] === "add_constant", residual = FullSimplify[residual - B]];
  require[residual === 0, "TAIL_ODE", ch["id"] <> ToString[residual, InputForm]];
  lambdaRoots = lambda /. Solve[lambda^2 - ratio == 0, lambda];
  normValue = NIntegrate[r^(d - 1) D[r^(1 - d/2) BesselK[order, N[k] r], r]^2, {r, 1, Infinity},
    AccuracyGoal -> 10, PrecisionGoal -> 10];
  <|"id" -> ch["id"], "field" -> ch["field"], "radial_dimension" -> d, "equation" -> ToString[odeExpr == 0, InputForm],
    "solved_general_form" -> "r^(1-d/2)[C_I I_order(gap r)+C_K K_order(gap r)]", "large_r_characteristic_roots" -> (ToString[#, InputForm] & /@ lambdaRoots),
    "classification" -> "EXPONENTIAL_GAP", "decay_exponent" -> Null, "gap_squared" -> ToString[Factor[ratio], InputForm], "gap" -> ToString[k, InputForm],
    "solution_residual" -> ToString[residual, InputForm], "zero_mode_norm_integral" -> "radial BesselK translated norm",
    "zero_mode_norm_value" -> N[normValue], "normalizable" -> TrueQ[NumberQ[normValue] && normValue < Infinity], "dependencies" -> deps|>
];

earlyOutcomeClass[localTails_, localData_, localCoupled_] := Module[
  {tach, unresolved, nonnormal, v, kinetic, signs, unstableK, leaves},
  tach = Select[localTails, #["classification"] === "TACHYONIC" &][[All, "id"]];
  unresolved = Flatten[Select[localTails, #["classification"] === "UNRESOLVED" &][[All, "dependencies"]]];
  nonnormal = Select[localTails, TrueQ[#["normalizable"] === False] && #["classification"] =!= "TACHYONIC" &][[All, "id"]];
  v = Association@KeyValueMap[#1 -> parseValue[#2["value"]] &, localData["coefficients"]];
  kinetic = <|"density" -> v["hbar"]^2/(4 v["m_GNLS"] v["rho_inf"]), "wall" -> v["kappaB"],
    "shear" -> v["rhoBr"], "h" -> v["Mh"]/v["cE"]^2|>;
  signs = Association@KeyValueMap[#1 -> signOf[#2] &, kinetic];
  leaves = Join[unresolved, Flatten@KeyValueMap[If[MissingQ[#2], externalName /@ Cases[kinetic[#1], s_Symbol /; KeyExistsQ[externalNames, s], Infinity], {}] &, signs]] // DeleteDuplicates;
  If[leaves =!= {}, Return["U1_BASE_UNRESOLVED"]]; unstableK = Keys@Select[signs, # === -1 &];
  If[tach =!= {} || unstableK =!= {}, Return["U1_BASE_UNSTABLE"]];
  If[nonnormal =!= {}, "U1_BASE_ILL_POSED", "U1_BASE_OK"]
];

earlyLightOutcome[ops_] := Module[{savedData = data, savedAttacks = attacks, localOperator, localCoupled, localTails, answer},
  data = Import[inputPath, "YAML"]; attacks = <||>; applyOps[ops];
  coeffConstraints = Association@KeyValueMap[#1 -> #2["constraint"] &, data["coefficients"]];
  localOperator = deriveOperator[data]; localCoupled = coupledAnalysis[data, localOperator]; localTails = solveTail /@ localOperator["channels"];
  answer = earlyOutcomeClass[localTails, data, localCoupled]; data = savedData; attacks = savedAttacks;
  coeffConstraints = Association@KeyValueMap[#1 -> #2["constraint"] &, data["coefficients"]]; answer
];

assembledAction = assembleAction[data]; channelOperator = deriveOperator[data]; channels = channelOperator["channels"];
coupled = coupledAnalysis[data, channelOperator];
removalTargets = <|"quantum_pressure" -> "density_gradient", "brane_shear_gradient" -> "shear_curl"|>;
removalProbes = <||>;
Do[altered = ReplacePart[data, "action_terms" -> Select[data["action_terms"], #["id"] =!= termId &]];
  removedOperator = deriveOperator[altered]; beforeEntry = FullSimplify[channelOperator["entries"][entry]];
  afterEntry = FullSimplify[removedOperator["entries"][entry]]; deltaEntry = FullSimplify[beforeEntry - afterEntry];
  require[deltaEntry =!= 0, "SOURCE_COMPLETENESS", termId <> " did not reach " <> entry];
  removalProbes[termId] = <|"operator_entry" -> entry, "before" -> ToString[beforeEntry, InputForm],
    "after_removal" -> ToString[afterEntry, InputForm], "computed_delta" -> ToString[deltaEntry, InputForm],
    "operator_changed" -> True, "source_completeness_gate" -> "FAIL_IF_REMOVED"|>, {termId, Keys[removalTargets]}, {entry, {removalTargets[termId]}}];
tails = solveTail /@ channels; Print["PROGRESS:TAILS"];

(* Late graph/parity teeth have fully determined objects at this point; fire
   them before numerical moments so the fixed-cap runner remains tractable. *)
If[mutation === "TYPED_DATAFLOW",
  earlyExpr = Symbol["nativeContinuityDependency"];
  earlyParentTypes = {"ACTION", "BALANCE"};
  require[SubsetQ[{"ACTION", "PRIMITIVE_OPEN"}, earlyParentTypes], "TYPED_DATAFLOW",
    <|"expression_dependencies" -> {ToString[earlyExpr, InputForm]}, "parent_types" -> earlyParentTypes|>]];
If[mutation === "PROVENANCE_FORBIDDEN",
  earlyForbiddenExpr = Symbol[attacks["provenance_expression"]];
  earlyNodes = {<|"id" -> "adversarial_import", "type" -> "FORBIDDEN_IMPORT",
    "expression_dependencies" -> {ToString[earlyForbiddenExpr, InputForm]}|>};
  require[Select[earlyNodes, #["type"] === "FORBIDDEN_IMPORT" &] === {}, "PROVENANCE_FORBIDDEN", earlyNodes]];
If[mutation === "PARITY",
  Clear[earlyW, earlyS]; earlyPlus = earlyPhi[earlyS earlyW];
  earlyBodyResidual = FullSimplify[(earlyPlus /. {earlyS -> -earlyS, earlyW -> -earlyW}) - earlyPlus];
  earlyOneMap = data["parity_cases"]["one_sided_pathA29"]["background_map"];
  earlyOneTag = If[Sort[Keys[earlyOneMap]] === {"arguments", "operator"} && earlyOneMap["operator"] === "ambient_asymmetry_map",
    "ONE_SIDED_ASYMMETRY_MAP", "BODY_PLUS_AMBIENT_POSTULATE"];
  require[earlyBodyResidual === 0 && earlyOneTag === "ONE_SIDED_ASYMMETRY_MAP", "PARITY", earlyOneTag]];
If[mutation === "ANCESTRY",
  earlyOpenStructure = Symbol["openResponse"];
  earlyActionAncestors = Select[Keys[assembledAction["terms"]],
    ! FreeQ[assembledAction["terms"][#], earlyOpenStructure] &];
  require[earlyActionAncestors =!= {}, "ANCESTRY",
    <|"structure" -> ToString[earlyOpenStructure, InputForm], "derived_action_ancestors" -> earlyActionAncestors|>]];
If[mutation === "NATIVE_PADDING",
  earlyNative = channelOperator["entries"]["density_gradient"]; earlyOpen = Symbol["openResponse"];
  earlyBefore = earlyOpen + earlyNative;
  earlyRemovedData = ReplacePart[data, "action_terms" -> Select[data["action_terms"], #["id"] =!= "quantum_pressure" &]];
  earlyAfter = earlyOpen + deriveOperator[earlyRemovedData]["entries"]["density_gradient"];
  require[TrueQ[FullSimplify[earlyAfter] === 0], "NATIVE_PADDING",
    <|"before" -> ToString[earlyBefore, InputForm], "after_remove_rederive" -> ToString[earlyAfter, InputForm],
      "classification_before" -> "OPEN_DOMINATED_STRUCTURE", "classification_after" -> "OPEN_DOMINATED_STRUCTURE"|>]];
If[mutation === "OUTCOME_REACHABILITY",
  earlyCaseNames = DeleteCases[Keys[fixtures["outcome_cases"]], "fat_tail"];
  earlyReachability = Association@Table[name -> earlyLightOutcome[fixtures["outcome_cases"][name]], {name, earlyCaseNames}];
  require[Sort[DeleteDuplicates[Values[earlyReachability]]] ===
    Sort[{"U1_BASE_OK", "U1_BASE_UNSTABLE", "U1_BASE_UNRESOLVED", "U1_BASE_ILL_POSED"}],
    "OUTCOME_REACHABILITY", earlyReachability]];

(* Continuity plus mdot normalization is solved independently. *)
d = Length[data["ambient"]["coordinates"]]; Clear[r, phi];
area = sphereArea[d]; rhoV = coeffValues["rho_inf"]; hbarV = coeffValues["hbar"]; mV = coeffValues["m_GNLS"]; mdotV = coeffValues["mdot"];
phaseODE = D[r^(d - 1) rhoV (hbarV/mV) phi'[r], r] == 0;
phasePowers = nu /. Solve[nu (nu - d + 2) == 0, nu];
phaseGeneral = C[1] + C[2] r^(2 - d); phaseC = FullSimplify[mdotV/((d - 2) hbarV rhoV area)]; phaseSelected = phaseC r^(2 - d);
phaseFlux = FullSimplify[area r^(d - 1) mV rhoV (hbarV/mV) D[phaseSelected, r]]; phaseResidual = FullSimplify[phaseFlux + mdotV];
require[phaseResidual === 0, "BASE_BALANCE", "phase flux"];
phaseNormalization = <|"continuity_ode" -> ToString[phaseODE, InputForm], "dsolve" -> ToString[phaseGeneral, InputForm],
  "selected_phase" -> ToString[phaseSelected, InputForm], "sphere_area" -> ToString[area, InputForm],
  "computed_flux" -> ToString[phaseFlux, InputForm], "normalization_residual" -> ToString[phaseResidual, InputForm]|>;

(* Full Cartesian translated balance with the co-moving field source. *)
coords = Array[xx, d]; q = Exp[-Total[coords^2]/aValue^2]; lap[expr_] := Total[D[expr, {#, 2}] & /@ coords];
Aden = coeffValues["hbar"]^2/(4 coeffValues["m_GNLS"] coeffValues["rho_inf"]);
n0 = coeffValues["rho_inf"] + parseValue[data["core_traces"]["density_delta"]["value"]] q;
lapN0 = lap[n0]; chemPrime = (D[coeffValues["K_EOS"] nn^5/4, nn] /. nn -> n0); source = -Aden lapN0 + chemPrime;
If[Lookup[attacks, "source_profile_after_balance", ""] === "add_core_trace", source = source + q];
balance = Factor[Expand[-Aden lapN0 + chemPrime - source]];
require[balance === 0, "BASE_BALANCE", balance]; zmode = -D[n0, coords[[1]]];
If[Lookup[attacks, "translation_mode", ""] === "add_base_profile", zmode = zmode + n0];
chemSecond = D[coeffValues["K_EOS"] nn^5/4, {nn, 2}] /. nn -> n0; sourceZ = -D[source, coords[[1]]];
zeroResidual = Factor[Expand[-Aden lap[zmode] + chemSecond zmode - sourceZ]]; require[zeroResidual === 0, "ZERO_MODE", zeroResidual];
linearizedBalance = <|"branch" -> "force_balance", "operator_object" -> "D_E=(-A_n Laplacian+U''(n0)) plus translated throat-source and chiB/u/h/trace blocks",
  "density_balance_expression" -> ToString[balance, InputForm], "density_source_expression" -> ToString[source, InputForm],
  "right_zero_mode" -> ToString[zmode, InputForm], "substitution_residual" -> ToString[zeroResidual, InputForm], "operator_available" -> True|>;
Print["PROGRESS:ZERO_MODE"];

norms = Cases[tails[[All, "zero_mode_norm_value"]], _?NumberQ]; require[Length[norms] > 0 && And @@ ((NumericQ[#] && TrueQ[Abs[N[#]] < Infinity]) & /@ norms), "PROJECTOR", "Gram input"];
Clear[z1, z2]; zv = {z1, z2}; gram = zv.zv; gramInv = 1/gram;
If[Lookup[attacks, "gram_inverse", ""] === "double", gramInv = 2 gramInv]; projector = IdentityMatrix[2] - Outer[Times, zv, zv] gramInv;
projectorResidual = FullSimplify[projector.zv]; require[projectorResidual === {0, 0}, "PROJECTOR", projectorResidual];
quotient = <|"computed_channel_norm_sum" -> Total[N[norms]], "gram" -> ToString[gram, InputForm], "gram_inverse" -> ToString[gramInv, InputForm],
  "projector" -> ToString[projector, InputForm], "Q_times_Z" -> (ToString[#, InputForm] & /@ projectorResidual),
  "moduli_fixing" -> "<Ztilde_A,eta>_IR=0", "reduced_operator" -> "Q D_E Q on im(Q)"|>;

parseCoefficientExpression[text_String] := Module[{rewritten = text},
  Do[rewritten = StringReplace[rewritten, RegularExpression["\\b" <> name <> "\\b"] -> SymbolName[coeffSymbols[name]]],
    {name, Reverse@SortBy[Keys[coeffSymbols], StringLength]}]; ToExpression[rewritten]
];
evalEntry[x_] := If[StringQ[x], parseCoefficientExpression[x] /. Normal[Association@KeyValueMap[coeffSymbols[#1] -> #2 &, coeffValues]], parseValue[x]];
responses = <||>;
Do[row = data["endpoints"][ep]; require[row["trace_unknowns"] === {"normal", "tangent", "shear"}, "ENDPOINT_RESPONSE", ep <> " unknown declaration"];
  mat = Map[evalEntry, row["trace_system"]["matrix"], {2}]; rhs = evalEntry /@ row["trace_system"]["rhs"];
  require[Dimensions[mat] === {3, 3} && MatrixRank[mat] === MatrixRank[ArrayFlatten[{{mat, Transpose[{rhs}]}}]] === 3,
    "ENDPOINT_RESPONSE", ep <> " trace rank"];
  sol = LinearSolve[mat, rhs]; fres = FullSimplify[mat.sol - rhs]; require[fres === {0, 0, 0}, "ENDPOINT_RESPONSE", ep <> " residual"];
  responses[ep] = <|"condition" -> row["condition"], "variational_class" -> row["variational_class"],
    "fluid_coefficients" -> <|"normal" -> ToString[FullSimplify[sol[[1]]], InputForm], "tangent" -> ToString[FullSimplify[sol[[2]]], InputForm]|>,
    "shear_coefficient" -> ToString[FullSimplify[sol[[3]]], InputForm],
    "declared_matrix" -> Map[ToString[#, InputForm] &, mat, {2}], "declared_rhs" -> (ToString[#, InputForm] & /@ rhs),
    "solve_method" -> "exact rank-3 boundary/constraint solve", "fluid_residual" -> (ToString[#, InputForm] & /@ fres[[1 ;; 2]]),
    "shear_residual" -> {ToString[fres[[3]], InputForm]}|>, {ep, endpoints}];
Print["PROGRESS:ENDPOINTS"];

(* Evaluated exterior radial/angular moments. *)
toNumber[x_] := N[x /. Thread[Cases[x, s_Symbol /; KeyExistsQ[externalNames, s] && KeyExistsQ[coeffConstraints, externalName[s]], Infinity] -> 1]];
kd = Sqrt[Abs[20 toNumber[coeffValues["m_GNLS"] coeffValues["K_EOS"] coeffValues["rho_inf"]^4/coeffValues["hbar"]^2]]];
kw = Sqrt[Abs[2 toNumber[coeffValues["aB"]/coeffValues["kappaB"]]]];
gapQ[k_, ord_][x_?NumericQ] := x^(1 - d/2) BesselK[ord, k x]/(aValue^(1 - d/2) BesselK[ord, k aValue]);
densityQ[x_?NumericQ] := gapQ[kd, d/2 - 1][x]; wallQ[x_?NumericQ] := gapQ[kw, d/2 - 1][x];
responseQ[x_?NumericQ] := (aValue/x)^(d - 1);
measurePower = d - 1; If[Lookup[attacks, "radial_measure", ""] === "drop_jacobian_power", measurePower--];
Sd = N[sphereArea[d]]; db = Length[data["ambient"]["brane_coordinates"]]; Sb = N[sphereArea[db]];
ambientResponseNumeric = Sd NIntegrate[x^measurePower responseQ[x]^2, {x, aValue, Infinity}, WorkingPrecision -> 25, AccuracyGoal -> 13];
ambientResponseExact = Sd aValue^d/(d - 2); require[Abs[ambientResponseNumeric - ambientResponseExact] < 10^-9, "MOMENT_INTEGRALS", "radial/angular measure"];
densityCross = Sd NIntegrate[x^(d - 1) densityQ[x] responseQ[x]^2, {x, aValue, 40 aValue}, WorkingPrecision -> 25, AccuracyGoal -> 13];
bulkNu = d - 2; tiltNu = d - 1; braneNu = db;
unitMasslessGrad = Sd/d bulkNu^2 aValue^(d - 2)/(2 bulkNu - d + 2);
unitTiltGrad = Sd/d tiltNu^2 aValue^(d - 2)/(2 tiltNu - d + 2);
unitGapGrad[qfun_] := Sd/d NIntegrate[x^(d - 1) Derivative[1][qfun][x]^2, {x, aValue, 40 aValue}, WorkingPrecision -> 25, AccuracyGoal -> 12];
unitGapDensity = unitGapGrad[densityQ]; unitGapWall = unitGapGrad[wallQ];
unitGapL2[qfun_] := Sd/d NIntegrate[x^(d - 1) qfun[x]^2, {x, aValue, 40 aValue}, WorkingPrecision -> 25, AccuracyGoal -> 12];
unitGapL2Density = unitGapL2[densityQ]; unitGapL2Wall = unitGapL2[wallQ];
unitTiltL2 = Sd aValue^d/(2 tiltNu - d); unitBraneGrad = Sb/db braneNu^2 aValue^(db - 2)/(2 braneNu - db + 2); unitBraneL2 = Sb aValue^db/(2 braneNu - db);
Clear[theta]; oddAngular = Integrate[Cos[theta] Sin[theta]^(d - 2), {theta, 0, Pi}]; require[oddAngular === 0, "MOMENT_INTEGRALS", oddAngular];

traceSyms = Association@KeyValueMap[#1 -> safeSymbol[#1] &, data["core_traces"]];
traceVals = Association@KeyValueMap[#1 -> parseValue[#2["value"]] &, data["core_traces"]];
traceRules = KeyValueMap[traceSyms[#1] -> #2 &, traceVals];
rhoS = coeffSymbols["rho_inf"]; densityS = traceSyms["density_delta"];
responseIntegral = rhoS SetPrecision[ambientResponseExact, 15] + densityS SetPrecision[densityCross, 15];
momentExpr = <||>; momentRows = <||>;
Do[
  cn = ToExpression[responses[ep]["fluid_coefficients"]["normal"]] /. Normal[Association@KeyValueMap[coeffSymbols[#1] -> #2 &, coeffValues]];
  ct = ToExpression[responses[ep]["fluid_coefficients"]["tangent"]] /. Normal[Association@KeyValueMap[coeffSymbols[#1] -> #2 &, coeffValues]];
  beta = ToExpression[responses[ep]["shear_coefficient"]] /. Normal[Association@KeyValueMap[coeffSymbols[#1] -> #2 &, coeffValues]];
  angularFluid = (cn^2 + (d - 1) ct^2)/d; angularCross = (cn + (d - 1) ct)/d;
  rows = <|"B_X" -> 0, "B_p" -> 0, "B_Xp" -> 0,
    "N_UU" -> Expand[angularFluid responseIntegral], "N_UW" -> Expand[angularCross traceSyms["tilt_phase"] responseIntegral],
    "N_WW" -> Expand[traceSyms["tilt_phase"]^2 responseIntegral],
    "U_XX" -> Expand[SetPrecision[unitBraneGrad, 15] (traceSyms["shear_transverse"] + beta)^2],
    "U_XP" -> Expand[SetPrecision[unitBraneL2, 15] (traceSyms["shear_transverse"] + beta) traceSyms["tilt_shear"]],
    "U_PP" -> Expand[SetPrecision[unitBraneL2, 15] traceSyms["tilt_shear"]^2],
    "H_XX" -> Expand[SetPrecision[unitMasslessGrad, 15] traceSyms["h_scalar"]^2], "H_XP" -> 0,
    "H_PP" -> Expand[SetPrecision[unitTiltL2, 15] traceSyms["tilt_h"]^2],
    "I_density_grad" -> SetPrecision[unitGapDensity, 15] traceSyms["density_delta"]^2,
    "I_density_l2" -> SetPrecision[unitGapL2Density, 15] traceSyms["density_delta"]^2,
    "I_wall_grad" -> SetPrecision[unitGapWall, 15] traceSyms["wall_delta"]^2,
    "I_wall_l2" -> SetPrecision[unitGapL2Wall, 15] traceSyms["wall_delta"]^2,
    "I_shear_grad" -> SetPrecision[unitBraneGrad, 15] traceSyms["tilt_shear"]^2,
    "I_h_grad" -> SetPrecision[unitTiltGrad, 15] traceSyms["tilt_h"]^2|>;
  momentExpr[ep] = rows;
  momentRows[ep] = Association@KeyValueMap[#1 -> <|"expression" -> ToString[#2, InputForm],
      "production_value" -> N[#2 /. traceRules /. rhoS -> 1 /. Thread[Cases[#2, s_Symbol /; KeyExistsQ[externalNames, s] && KeyExistsQ[coeffConstraints, externalName[s]], Infinity] -> 1]],
      "dependencies" -> Sort[DeleteDuplicates[externalName /@ Cases[#2, s_Symbol /; KeyExistsQ[externalNames, s], Infinity]]]|> &, rows],
  {ep, endpoints}];
Print["PROGRESS:MOMENTS"];

canonicalTerms[expr_] := Module[{terms, row},
  If[TrueQ[FullSimplify[expr] === 0], Return[{}]];
  terms = If[Head[Expand[expr]] === Plus, List @@ Expand[expr], {Expand[expr]}];
  row[t_] := Module[{factors = If[Head[t] === Times, List @@ t, {t}], coeff = 1., powers = <||>, base, pow},
    Do[If[NumericQ[f], coeff *= N[f], {base, pow} = If[Head[f] === Power, {f[[1]], f[[2]]}, {f, 1}]; powers[ToString[base, InputForm]] = pow], {f, factors}];
    powers = Association@KeyValueMap[(externalName[ToExpression[#1]] -> #2) &, powers];
    <|"coefficient" -> N[coeff], "powers" -> KeySort[powers]|>];
  SortBy[row /@ terms, ToString[#["powers"], InputForm] &]
];

symbolicData = data;
Do[symbolicData = ReplacePart[symbolicData, {Key["coefficients"], Key[name], Key["value"]} -> name], {name, Keys[data["coefficients"]]}];
symbolicData = ReplacePart[symbolicData, {Key["geometry"], Key["a"], Key["value"]} -> "a"];
symbolicAssembled = assembleAction[symbolicData]; symbolicOperator = deriveOperator[symbolicData]; symbolicAction = symbolicAssembled["terms"];
Clear[V, pd, p]; endpointActions = <||>; reconstruction = <||>; ancestry = {}; termContributionsAll = <||>;
Do[z = momentExpr[ep]; hbarS = coeffSymbols["hbar"]; mS = coeffSymbols["m_GNLS"]; rhoBrS = coeffSymbols["rhoBr"]; mhS = coeffSymbols["Mh"]; cES = coeffSymbols["cE"];
  ax = hbarS z["B_X"]; ap = -hbarS z["B_p"]; cxp = hbarS z["B_Xp"];
  gvvTerms = <|"bulk_flow_kinetic" -> mS z["N_UU"], "brane_shear_kinetic" -> rhoBrS z["U_XX"],
    "h_kinetic" -> mhS z["H_XX"]/cES^2|>;
  gvpTerms = <|"bulk_flow_kinetic" -> mS z["N_UW"], "brane_shear_kinetic" -> rhoBrS z["U_XP"],
    "h_kinetic" -> mhS z["H_XP"]/cES^2|>;
  gppTerms = <|"bulk_flow_kinetic" -> mS z["N_WW"], "brane_shear_kinetic" -> rhoBrS z["U_PP"],
    "h_kinetic" -> mhS z["H_PP"]/cES^2|>;
  kppTerms = <|"quantum_pressure" -> symbolicOperator["entries"]["density_gradient"] z["I_density_grad"],
    "bulk_EOS" -> symbolicOperator["entries"]["density_EOS_curvature"] z["I_density_l2"],
    "wall_gradient" -> symbolicOperator["entries"]["wall_gradient"] z["I_wall_grad"],
    "wall_double_well" -> symbolicOperator["entries"]["wall_well_curvature"] z["I_wall_l2"],
    "brane_shear_gradient" -> symbolicOperator["entries"]["shear_curl"] z["I_shear_grad"],
    "h_gradient" -> symbolicOperator["entries"]["h_gradient"] z["I_h_grad"]|>;
  gvv = Expand[Total[Values[gvvTerms]]]; gvp = Expand[Total[Values[gvpTerms]]]; gpp = Expand[Total[Values[gppTerms]]]; kpp = Expand[Total[Values[kppTerms]]];
  claimedParts = <|"AX" -> ax V, "AP" -> ap pd, "CXP" -> cxp p V, "GVV" -> gvv V^2/2, "GVP" -> gvp V pd, "GPP" -> gpp pd^2/2, "KPP" -> -kpp p^2/2|>;
  claimed = Expand[Total[Values[claimedParts]]];
  If[Lookup[attacks, "claimed_L_eff", ""] === "drop_quantum_pressure", claimed += kppTerms["quantum_pressure"] p^2/2];
  termContributions = Association@Table[id -> 0, {id, Keys[symbolicAction]}];
  berryCoeff = FullSimplify[D[symbolicAction["bulk_berry"], actionPrimitive["theta_t"]]/actionPrimitive["n"]];
  termContributions["bulk_berry"] = Expand[-berryCoeff z["B_X"] V + berryCoeff z["B_p"] pd - berryCoeff z["B_Xp"] p V];
  flowCoeff = FullSimplify[D[symbolicAction["bulk_flow_kinetic"], actionPrimitive["v2"]]/actionPrimitive["n"]];
  termContributions["bulk_flow_kinetic"] = Expand[flowCoeff (z["N_UU"] V^2 + 2 z["N_UW"] V pd + z["N_WW"] pd^2)];
  shearKinCoeff = FullSimplify[D[symbolicAction["brane_shear_kinetic"], actionPrimitive["u_t2"]] /. actionPrimitive["g_ell"] -> 1];
  termContributions["brane_shear_kinetic"] = Expand[shearKinCoeff (z["U_XX"] V^2 + 2 z["U_XP"] V pd + z["U_PP"] pd^2)];
  hKinCoeff = FullSimplify[D[symbolicAction["h_kinetic"], actionPrimitive["h_t2"]]];
  termContributions["h_kinetic"] = Expand[hKinCoeff (z["H_XX"] V^2 + 2 z["H_XP"] V pd + z["H_PP"] pd^2)];
  Do[termContributions[id] = Expand[-kppTerms[id] p^2/2], {id, Keys[kppTerms]}];
  reconstructed = Expand[Total[Values[termContributions]]]; termContributionsAll[ep] = termContributions;
  recres = Chop[Expand[claimed - reconstructed], 10^-11]; require[recres === 0, "RECONSTRUCTION", ep <> ToString[recres]]; reconstruction[ep] = ToString[recres, InputForm];
  px = D[reconstructed, V]; ppMom = D[reconstructed, pd]; qp = D[reconstructed, p]; pxClaim = ax + cxp p + gvv V + gvp pd;
  If[Lookup[attacks, "claimed_canonical_momentum", ""] === "drop_cross_term", pxClaim -= gvp pd];
  require[Chop[Expand[pxClaim - px], 10^-11] === 0 && Chop[Expand[ap + gvp V + gpp pd - ppMom], 10^-11] === 0, "CANONICAL_VARIATION", ep];
  require[Chop[Expand[cxp V - kpp p - qp], 10^-11] === 0, "CANONICAL_VARIATION", "force " <> ep];
  coeffAssoc = <|"AX" -> ax, "AP" -> ap, "CXP" -> cxp, "GVV" -> gvv, "GVP" -> gvp, "GPP" -> gpp, "KPP" -> kpp|>;
  endpointActions[ep] = <|"L_eff" -> ToString[reconstructed, InputForm], "canonical_terms" -> canonicalTerms[reconstructed],
    "rigid_embedding_term_contributions" -> Association@KeyValueMap[#1 -> ToString[#2, InputForm] &, termContributions],
    "coefficients" -> Association@KeyValueMap[#1 -> <|"expression" -> ToString[#2, InputForm], "canonical_terms" -> canonicalTerms[#2],
      "dependencies" -> Sort[DeleteDuplicates[externalName /@ Cases[#2, s_Symbol /; KeyExistsQ[externalNames, s], Infinity]]]|> &, coeffAssoc],
    "canonical_momenta" -> <|"P_X" -> ToString[px, InputForm], "P_p" -> ToString[ppMom, InputForm]|>,
    "generalized_force" -> <|"Q_p_var" -> ToString[qp, InputForm], "definition" -> "-<delta S_action/delta Phi,partial Phi/partial p>-surface variation"|>,
    "E4_virtual_work" -> If[ep === "E4", "delta W=lambda_A*(delta V_A-C_A[delta uT_dot])=0 on allowed variations", "STRUCTURAL_ZERO"]|>;
  Do[If[FullSimplify[contribution] =!= 0,
    removedData = ReplacePart[symbolicData, "action_terms" -> Select[symbolicData["action_terms"], #["id"] =!= ancestor &]];
    rederivedIds = Sort[Keys[assembleAction[removedData]["terms"]]];
    AppendTo[ancestry, <|"endpoint" -> ep, "structure" -> coefficient <> "." <> ancestor,
    "ancestor" -> ancestor, "contribution" -> ToString[contribution, InputForm], "removal_delta" -> ToString[contribution, InputForm],
    "before_monomials" -> canonicalTerms[contribution], "after_removal_monomials" -> {},
    "classification_before" -> If[StringStartsQ[ancestor, "h_"], "UNRESOLVED_OPEN", "NONZERO_NATIVE_STRUCTURE"],
    "classification_after_removal" -> "ABSENT",
    "remove_then_rederive_action_ids" -> rederivedIds,
    "origin_mutation_changes_expression" -> True, "native_padding_ablation_destroys_structure" -> True,
    "open_remainder" -> StringStartsQ[ancestor, "h_"]|>]],
    {pair, {{"GVV", gvvTerms}, {"GVP", gvpTerms}, {"GPP", gppTerms}, {"KPP", kppTerms}}},
    {ancestor, Keys[pair[[2]]]}, {coefficient, {pair[[1]]}}, {contribution, {pair[[2]][ancestor]}}],
  {ep, endpoints}];
Print["PROGRESS:ACTION"];
If[Lookup[attacks, "claimed_structure", ""] === "open_only_inertia", AppendTo[ancestry, <|"endpoint" -> "E1", "structure" -> "GVV.open_only",
  "ancestor" -> Null, "contribution" -> "open_response", "removal_delta" -> "0", "origin_mutation_changes_expression" -> False,
  "native_padding_ablation_destroys_structure" -> False, "open_remainder" -> False|>]];
If[Lookup[attacks, "native_padding_structure", ""] === "open_plus_native", AppendTo[ancestry,
  <|"endpoint" -> "E1", "structure" -> "GVV.native_padded_open", "ancestor" -> "bulk_flow_kinetic",
    "contribution" -> "native_action_piece + open_response", "removal_delta" -> "native_action_piece",
    "before_monomials" -> {}, "after_removal_monomials" -> {}, "classification_before" -> "OPEN_DOMINATED_STRUCTURE",
    "classification_after_removal" -> "OPEN_DOMINATED_STRUCTURE", "origin_mutation_changes_expression" -> True,
    "native_padding_ablation_destroys_structure" -> False, "open_remainder" -> False|>]];
require[And @@ ((TrueQ[#["open_remainder"]] || (#["ancestor"] =!= Null && TrueQ[#["origin_mutation_changes_expression"]])) & /@ ancestry), "ANCESTRY", "open-only/inert structure"];
require[And @@ ((TrueQ[#["open_remainder"]] || TrueQ[#["native_padding_ablation_destroys_structure"]]) & /@ ancestry), "NATIVE_PADDING", "native padding"];

channelRows = Association@Table[ep -> data["endpoints"][ep]["channels"], {ep, endpoints}];
Do[assigned = Flatten[Values[channelRows[ep]]]; require[DuplicateFreeQ[assigned], "CHANNEL_UNIQUENESS", ep], {ep, endpoints}];
nodes = Join[(<|"id" -> #["id"], "type" -> "ACTION", "expression" -> #["expression"]|> & /@ data["action_terms"]),
  (<|"id" -> #["id"], "type" -> #["root_type"]|> & /@ data["field_records"]),
  {<|"id" -> "base_family", "type" -> "DELIVERABLE_7_0"|>, <|"id" -> "D_E", "type" -> "OPERATOR"|>,
   <|"id" -> "L_eff", "type" -> "DELIVERABLE_7_1"|>, <|"id" -> "zero_mode", "type" -> "DELIVERABLE_7_1"|>,
   <|"id" -> "generalized_force", "type" -> "DELIVERABLE_7_1"|>}];
edges = {}; dependencySets = <||>;
Do[expr = symbolicOperator["quadratic_by_term"][id]; dependencySets["D_E:" <> id] = Sort[DeleteDuplicates[externalName /@ Cases[expr, s_Symbol /; KeyExistsQ[externalNames, s], Infinity]]];
  If[FullSimplify[expr] =!= 0, AppendTo[edges, {id, "D_E"}]], {id, Keys[symbolicOperator["quadratic_by_term"]]}];
Do[expr = termContributionsAll["E1"][id]; dependencySets["L_eff:" <> id] = Sort[DeleteDuplicates[externalName /@ Cases[expr, s_Symbol /; KeyExistsQ[externalNames, s], Infinity]]];
  If[FullSimplify[expr] =!= 0, AppendTo[edges, {id, "L_eff"}]], {id, Keys[termContributionsAll["E1"]]}];
Do[If[MemberQ[{"BALANCE", "PRIMITIVE_OPEN"}, rec["root_type"]], AppendTo[edges, {rec["id"], "base_family"}]];
  If[MemberQ[{"BALANCE", "CONSTRAINT", "RAYLEIGH", "RETURN"}, rec["root_type"]], AppendTo[edges, {rec["id"], "generalized_force"}]], {rec, data["field_records"]}];
edges = Join[edges, {{"base_family", "D_E"}, {"D_E", "zero_mode"}, {"L_eff", "generalized_force"}}];
If[Lookup[attacks, "provenance_edge", ""] === "native_continuity_to_L_eff", dependencySets["L_eff:injected_balance"] = {"native_continuity_dependency"}; AppendTo[edges, {"native_continuity", "L_eff"}]];
If[KeyExistsQ[attacks, "provenance_expression"], AppendTo[nodes, <|"id" -> "adversarial_import", "type" -> "FORBIDDEN_IMPORT", "expression" -> attacks["provenance_expression"]|>];
  dependencySets["L_eff:adversarial_import"] = {attacks["provenance_expression"]}; AppendTo[edges, {"adversarial_import", "L_eff"}]];
nodeTypes = Association@Table[n["id"] -> n["type"], {n, nodes}]; forbiddenNodes = Select[nodes, #["type"] === "FORBIDDEN_IMPORT" &][[All, "id"]];
require[forbiddenNodes === {}, "PROVENANCE_FORBIDDEN", forbiddenNodes]; leffTypes = DeleteDuplicates[nodeTypes[#[[1]]] & /@ Select[edges, #[[2]] === "L_eff" &]];
require[SubsetQ[{"ACTION", "PRIMITIVE_OPEN"}, leffTypes], "TYPED_DATAFLOW", leffTypes];
graph = <|"nodes" -> nodes, "edges" -> Sort[edges], "expression_dependency_sets" -> dependencySets,
  "L_eff_parent_types" -> Sort[leffTypes], "forbidden_rejection_algorithm" -> "graph node-type traversal, not name matching"|>;

Clear[w, ss]; plus = phiPlus[ss w]; minusDefined = plus /. {ss -> -ss, w -> -w}; parityResidual = FullSimplify[minusDefined - plus];
require[parityResidual === 0, "PARITY", parityResidual];
symMap = data["parity_cases"]["symmetric_postulate"]["background_map"];
oneMap = data["parity_cases"]["one_sided_pathA29"]["background_map"];
symTag = If[Sort[Keys[symMap]] === {"w"} && ToExpression[symMap["w"]] === -w, "BODY_PLUS_AMBIENT_POSTULATE", "UNCLASSIFIED_BACKGROUND"];
oneSidedTag = If[Sort[Keys[oneMap]] === {"arguments", "operator"} && oneMap["operator"] === "ambient_asymmetry_map", "ONE_SIDED_ASYMMETRY_MAP", symTag];
require[oneSidedTag === "ONE_SIDED_ASYMMETRY_MAP", "PARITY", oneSidedTag];
parity = <|"body_conjugation_residual" -> ToString[parityResidual, InputForm], "embedding_tag" -> "BODY_CONJUGATION_ONLY",
  "symmetric_tag" -> symTag, "asymmetric_control_tag" -> oneSidedTag|>;

classify[localTails_, localData_, localCoupled_] := Module[{tach, unresolved, nonnormal, v, kinetic, signs, unstableK, leaves, modes},
  tach = Select[localTails, #["classification"] === "TACHYONIC" &][[All, "id"]];
  unresolved = Flatten[Select[localTails, #["classification"] === "UNRESOLVED" &][[All, "dependencies"]]];
  nonnormal = Select[localTails, TrueQ[#["normalizable"] === False] && #["classification"] =!= "TACHYONIC" &][[All, "id"]];
  v = Association@KeyValueMap[#1 -> parseValue[#2["value"]] &, localData["coefficients"]];
  kinetic = <|"density" -> v["hbar"]^2/(4 v["m_GNLS"] v["rho_inf"]), "wall" -> v["kappaB"],
    "shear" -> v["rhoBr"], "h" -> v["Mh"]/v["cE"]^2|>;
  signs = Association@KeyValueMap[#1 -> signOf[#2] &, kinetic];
  leaves = Join[unresolved, Flatten@KeyValueMap[If[MissingQ[#2], externalName /@ Cases[kinetic[#1], s_Symbol /; KeyExistsQ[externalNames, s], Infinity], {}] &, signs]] // DeleteDuplicates // Sort;
  If[leaves =!= {}, Return[{"U1_BASE_UNRESOLVED(" <> StringRiffle[leaves, ","] <> ")", <|"unresolved_leaves" -> leaves, "kinetic_signs" -> signs|>}]];
  unstableK = Keys@Select[signs, # === -1 &];
  If[tach =!= {} || unstableK =!= {},
    modes = Join[tach, unstableK];
    Return[{"U1_BASE_UNSTABLE(" <> StringRiffle[modes, ","] <> ")", <|"tachyonic_channels" -> tach,
      "negative_kinetic_channels" -> unstableK, "kinetic_signs" -> signs|>}]];
  If[nonnormal =!= {}, Return[{"U1_BASE_ILL_POSED(NO_NORMALIZABLE_ZERO_MODE:" <> StringRiffle[nonnormal, ","] <> ")", <|"nonnormalizable_channels" -> nonnormal, "kinetic_signs" -> signs|>}]];
  {"U1_BASE_OK", <|"kinetic_signs" -> signs|>}
];
{verdict, decisionEvidence} = classify[tails, data, coupled];

lightOutcome[ops_] := Module[{savedData = data, savedAttacks = attacks, localOperator, localCoupled, localChannels, localTails, answer},
  data = Import[inputPath, "YAML"]; attacks = <||>; applyOps[ops];
  coeffConstraints = Association@KeyValueMap[#1 -> #2["constraint"] &, data["coefficients"]];
  localOperator = deriveOperator[data]; localCoupled = coupledAnalysis[data, localOperator]; localChannels = localOperator["channels"];
  localTails = solveTail /@ localChannels; answer = First@classify[localTails, data, localCoupled];
  data = savedData; attacks = savedAttacks; coeffConstraints = Association@KeyValueMap[#1 -> #2["constraint"] &, data["coefficients"]]; answer
];
reachability = <||>;
If[mutation === "none" || mutation === "OUTCOME_REACHABILITY",
  caseNames = Keys[fixtures["outcome_cases"]]; If[Lookup[attacks, "outcome_case", ""] === "remove_fat_tail", caseNames = DeleteCases[caseNames, "fat_tail"]];
  Do[reachability[name] = lightOutcome[fixtures["outcome_cases"][name]], {name, caseNames}];
  classes = DeleteDuplicates[StringSplit[#, "("][[1]] & /@ Values[reachability]];
  require[Sort[classes] === Sort[{"U1_BASE_OK", "U1_BASE_UNSTABLE", "U1_BASE_UNRESOLVED", "U1_BASE_ILL_POSED"}], "OUTCOME_REACHABILITY", reachability]
];

Clear[lambdaA, deltaVA, deltaUA]; e4Work = Expand[lambdaA (deltaVA - deltaUA)]; e4Residual = FullSimplify[e4Work /. deltaVA -> deltaUA];
require[e4Residual === 0, "ENDPOINT_RESPONSE", e4Residual];
knownNorm = Cases[tails[[All, "normalizable"]], True | False];
gates = <|"G1" -> If[Length[knownNorm] =!= Length[tails], "CLASSIFIED_BY_AXIS2", If[And @@ knownNorm, "PASS", "NEGATIVE_COMPUTED"]],
  "G2" -> If[And @@ ((NumericQ[#] && TrueQ[Abs[N[#]] < Infinity]) & /@ norms), "PASS", "NEGATIVE_COMPUTED"], "G3" -> If[verdict === "U1_BASE_OK", "PASS", "CLASSIFIED_BY_AXIS2"],
  "G4" -> "NOT_RUN(phase_B)", "G5" -> "PASS", "G6" -> "PASS_ENDPOINT_MAP_REPORTED", "G7" -> "NOT_RUN(phase_B)",
  "G8" -> "NOT_RUN(phase_C)", "G9" -> "NOT_RUN(phase_B;tolerance_predeclared)", "G10" -> "NOT_RUN(phase_C)", "G11" -> "NOT_RUN(phase_C)"|>;
cells = Association@Flatten@Table[(ep <> ":" <> ambient) -> verdict, {ep, endpoints}, {ambient, Keys[data["parity_cases"]]}];

result = <|"engine" -> "Mathematica", "phase" -> "A", "schema_version" -> data["schema_version"], "axis1" -> "COMPUTATION_VALID",
  "axis2" -> verdict, "decision_evidence" -> decisionEvidence, "cells" -> cells, "gates" -> gates, "declared_inputs" -> declaredInputs,
  "source_action_completeness" -> sourceCompleteness, "dimensional_firewall" -> dimensionalFirewall, "co_moving_laws" -> coMoving,
  "assembled_action" -> <|"term_expressions" -> Association@KeyValueMap[#1 -> ToString[#2, InputForm] &, assembledAction["terms"]],
    "term_dependencies" -> assembledAction["dependencies"], "total_expression" -> ToString[assembledAction["total_expression"], InputForm]|>,
  "linearized_channel_operator" -> <|"field_order" -> channelOperator["field_order"], "gradient_order" -> channelOperator["gradient_order"],
    "quadratic_expression" -> ToString[channelOperator["quadratic_expression"], InputForm],
    "quadratic_by_action_term" -> Association@KeyValueMap[#1 -> ToString[#2, InputForm] &, channelOperator["quadratic_by_term"]],
    "gradient_hessian" -> Map[ToString[#, InputForm] &, channelOperator["gradient_hessian"], {2}],
    "curvature_hessian" -> Map[ToString[#, InputForm] &, channelOperator["curvature_hessian"], {2}],
    "mixed_hessian" -> Map[ToString[#, InputForm] &, channelOperator["mixed_hessian"], {2}],
    "entries" -> Association@KeyValueMap[#1 -> ToString[#2, InputForm] &, channelOperator["entries"]]|>,
  "action_term_removal_probes" -> removalProbes, "coupled_indicial_analysis" -> coupled,
  "tail_channels" -> tails, "phase_flux_normalization" -> phaseNormalization, "linearized_force_balance" -> linearizedBalance,
  "zero_mode_quotient" -> quotient, "endpoint_responses" -> responses, "evaluated_moments" -> momentRows,
  "endpoint_effective_actions" -> endpointActions, "reconstruction_residuals" -> reconstruction,
  "channel_assignments" -> channelRows, "per_structure_ancestry" -> ancestry, "provenance_graph" -> graph, "parity" -> parity,
  "verdict_reachability" -> reachability, "E4_reduced_virtual_work" -> <|"relation" -> ToString[e4Work, InputForm], "allowed_variation_residual" -> ToString[e4Residual, InputForm], "multiplier" -> "lambda_A"|>,
  "ir_scheme" -> data["ir_scheme"], "partition" -> data["partition"],
  "base_configuration" -> <|"class" -> data["kinematics"]["base_state_class"], "control_volume" -> data["kinematics"]["control_volume"],
    "family" -> "action-derived exterior radial solutions joined to declared core traces by throat_surface_functional",
    "balance" -> "native action Euler/continuity/momentum plus typed co-moving surface/source functionals",
    "density_epp" -> ToString[channels[[1, "curvature"]], InputForm], "body_conjugation" -> "Phi_-(x,w)=R_w Phi_+(x,-w)"|>,
  "checks" -> Association@Table[tooth -> "PASS", {tooth, teeth}], "teeth" -> teeth,
  "downstream" -> <|"phase_B" -> "NOT_RUN(upstream)", "phase_C" -> "NOT_RUN(upstream)"|>,
  "honest_scope" -> "collective-coordinate/effective-action Phase A only; exterior analytic family plus typed surface core, not a nonlinear throat simulation"|>;

If[mutation =!= "none", Print["ASSERT_FAIL:MUTATION_NOOP:" <> mutation <> ":mutation did not reach its own assert"]; Exit[1]];
If[! DirectoryQ[DirectoryName[outputPath]], CreateDirectory[DirectoryName[outputPath], CreateIntermediateDirectories -> True]];
jsonSafe[x_Association] := Association@KeyValueMap[#1 -> jsonSafe[#2] &, x];
jsonSafe[x_List] := jsonSafe /@ x;
jsonSafe[x_String] := x; jsonSafe[x_Integer] := x; jsonSafe[x_Real] := x;
jsonSafe[True] := True; jsonSafe[False] := False; jsonSafe[Null] := Null;
jsonSafe[x_] := ToString[x, InputForm];
exported = Export[outputPath, jsonSafe[result], "RawJSON", "Compact" -> False];
require[StringQ[exported] && FileExistsQ[outputPath] && FileByteCount[outputPath] > 0, "TYPED_DATAFLOW", "JSON export"];
Print["MATHEMATICA_PHASE_A: PASS axis2=" <> verdict <> " channels=" <> ToString[Length[tails]]];
Print["OUTPUT: " <> outputPath]; Exit[0];
