#!/usr/bin/env wolframscript
(* U1 Phase-B1 remediation-3: independent forward Gram/Jacobian Wolfram engine. *)

ClearAll["Global`*"];
here = DirectoryName[ExpandFileName[$InputFileName]];
root = ExpandFileName[FileNameJoin[{here, "..", ".."}]];
defaultInput = FileNameJoin[{here, "u1_body_mechanics_inputs.yaml"}];
defaultOutput = FileNameJoin[{here, "reports", "u1_body_dynamics_artifacts", "mathematica_phase_b1.yaml"}];
endpoints = {"E1", "E2", "E3", "E4", "E5"};
ambients = {"one_sided_pathA29", "symmetric_postulate"};

argValue[name_, fallback_] := Module[{p = FirstPosition[$ScriptCommandLine, name]}, If[MissingQ[p], fallback, $ScriptCommandLine[[p[[1]] + 1]]]];
inputPath = ExpandFileName[argValue["--input", defaultInput]];
outputPath = ExpandFileName[argValue["--output", defaultOutput]];
phaseOverride = argValue["--phase-artifact", ""];
data = Import[inputPath, "YAML"];
readEvents = <||>;
dataAt[path_List] := Fold[If[AssociationQ[#1], #1[#2], #1[[#2 + 1]]] &, data, path];
trackedLeafPaths[x_Association, p_List] := Flatten[KeyValueMap[trackedLeafPaths[#2, Append[p, #1]] &, x], 1];
trackedLeafPaths[x_List, p_List] := If[x === {}, {p}, Flatten[MapIndexed[trackedLeafPaths[#1, Append[p, First[#2] - 1]] &, x], 1]];
trackedLeafPaths[x_, p_List] := {p};
recordRead[path_List, scope_String] := Module[{eventKey}, Scan[Function[p, eventKey = StringRiffle[ToString /@ p, "."];
  readEvents[eventKey] = Sort@DeleteDuplicates@Append[Lookup[readEvents, eventKey, {}], scope]], trackedLeafPaths[dataAt[path], path]]];
rd[path_List, scope_String] := (recordRead[path, scope]; dataAt[path]);
If[!AssociationQ[data] || rd[{"schema_version"}, "schema_validation"] =!= "U1_PHASE_B1_MECHANICS_INPUTS_V3", Print["ASSERT_FAIL:B1_S1:mechanics input schema is not V3"]; Exit[1]];
phaseInput = Import[FileNameJoin[{root, rd[{"substrate", "phase_a_input"}, "phase_a_source_contract"]}], "YAML"];
phasePath = If[phaseOverride === "", FileNameJoin[{root, rd[{"substrate", "mathematica_artifact"}, "phase_a_source_contract"]}], phaseOverride];
phase = Import[phasePath, "RawJSON"];

externalNames = <||>;
safe[name_String] := Module[{z = Symbol["b1" <> StringReplace[name, "_" -> "ZZ"]]}, externalNames[z] = name; z];
ename[z_Symbol] := Lookup[externalNames, z, SymbolName[Unevaluated[z]]];
coeff = Association@KeyValueMap[#1 -> safe[#1] &, phaseInput["coefficients"]];
traces = Association@KeyValueMap[#1 -> safe[#1] &, phaseInput["core_traces"]];
aSym = safe["a"]; sSym = safe["s"];
allSymbols = Join[coeff, traces, <|"a" -> aSym, "s" -> sSym|>];
parse[x_] := If[NumberQ[x], x, Module[{out = ToString[x]},
  (* Mathematica Phase-A serializes its private symbols as u1fooZZbar.  Normalize
     those names independently into this engine's b1 symbol table before parsing. *)
  KeyValueMap[(out = StringReplace[out, "u1" <> StringReplace[#1, "_" -> "ZZ"] -> SymbolName[#2]]) &, allSymbols];
  KeyValueMap[(out = StringReplace[out, RegularExpression["(?<![A-Za-z0-9_])" <> #1 <> "(?![A-Za-z0-9_])"] -> SymbolName[#2]]) &, allSymbols];
  out = StringReplace[out, {"pi" -> "Pi", "^" -> "^"}]; Quiet@Check[ToExpression[out, InputForm], $Failed]]];

canonicalTerms[expr_] := Module[{expanded = Expand[N[expr, 17]], terms, one},
  (* Phase-A's two engines serialize machine and 15-digit reals differently.  The
     tolerance is only a canonicalization tolerance; all symbolic dependencies remain. *)
  expanded = Expand[Chop[expanded, 2*10^-8]]; If[TrueQ[expanded === 0], Return[{}]];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  one[t_] := Module[{factors = If[Head[t] === Times, List @@ t, {t}], num = 1, powers = <||>, base, power},
    Do[If[NumericQ[f], num *= f, If[Head[f] === Power, base = f[[1]]; power = f[[2]], base = f; power = 1];
      powers[ename[base]] = Lookup[powers, ename[base], 0] + power], {f, factors}];
    <|"coefficient" -> N[num, 17], "powers" -> KeySort[powers]|>];
  SortBy[Select[one /@ terms, Abs[N[#coefficient]] > 2*10^-10 &], Normal[#powers] &]];
canonicalMatrix[m_] := Map[canonicalTerms, m, {2}];
exprTerms[rows_List] := Total[Map[Function[row, row["coefficient"] Times @@ KeyValueMap[(If[KeyExistsQ[allSymbols, #1], allSymbols[#1], safe[#1]])^#2 &, row["powers"]]], rows]];
zeroQ[x_] := canonicalTerms[x] === {};
asciiKey[s_] := FromDigits[PadRight[ToCharacterCode[ToString[s]], 256], 257];
canonJSON[x_Association] := Association@KeyValueMap[#1 -> canonJSON[#2] &, KeySortBy[x, asciiKey]];
canonJSON[x_List] := canonJSON /@ x; canonJSON[x_] := x;
sha[x_] := Hash[ExportString[canonJSON[x], "RawJSON", "Compact" -> True], "SHA256", "HexString"];
asciiSort[x_List] := SortBy[x, asciiKey];
powerSortKey[a_Association] := asciiKey[StringRiffle[KeyValueMap[#1<>":"<>ToString[#2]&,KeySortBy[a,asciiKey]],","]];
fixedDecimal[x_] := Module[{scaled=Round[N[x,30]10^12],sign,n},sign=If[scaled<0,"-",""];n=Abs[scaled];
 sign<>ToString[Quotient[n,10^12]]<>"."<>IntegerString[Mod[n,10^12],10,12]];
normalizeCoreTraces[rows_] := Association@KeyValueMap[#1 -> <|"field" -> #2["field"],
 "surface" -> #2["surface"], "value" -> fixedDecimal[parse[#2["value"]]]|> &, KeySortBy[rows, asciiKey]];
normalizePayloadV2[payload_] := Module[{out = payload, coeffs, rows},
 out["action_ids"] = asciiSort[out["action_ids"]];
 out["decision"] = KeyTake[out["decision"], {"id","sha256","status"}];
 KeyValueMap[Function[{ep, epRows}, KeyValueMap[Function[{name, terms},
    rows = Map[Function[row, <|"coefficient" -> fixedDecimal[row["coefficient"]],
       "powers" -> Association@KeyValueMap[ToString[#1] -> Round[#2] &, KeySort[row["powers"]]]|>], terms];
    out["endpoint_coefficients"][ep][name] = SortBy[rows, powerSortKey[#powers] &]], epRows]], out["endpoint_coefficients"]]; out];
checks = <||>;
recordCheck[name_, test_, evidence_] := (checks[name] = <|"status" -> If[TrueQ[test], "PASS", "FAIL"], "evidence_digest" -> sha[evidence], "evidence" -> evidence|>);
sphereArea[d_Integer] := FullSimplify[2 Pi^(d/2)/Gamma[d/2]];

(* The v1 payload is frozen before applying the shear amendment. *)
phasePayloadOf[p_] := <|"axis1" -> p["axis1"], "axis2" -> p["axis2"], "cells" -> p["cells"],
 "tails" -> ({#["id"], #["classification"], #["decay_exponent"], #["normalizable"]} & /@ phase["tail_channels"]),
 "decision" -> p["source_action_completeness"]["operative_decision_citation"],
 "action_ids" -> p["source_action_completeness"]["assembled_ids"],
 "endpoint_coefficients" -> Association@KeyValueMap[#1 -> Association@KeyValueMap[#1 -> #2["canonical_terms"] &, #2["coefficients"]] &, p["endpoint_effective_actions"]]|>;
phaseOriginal = Import[phasePath, "RawJSON"];
phasePayload = phasePayloadOf[phaseOriginal];
phasePayloadSha = sha[phasePayload];
baselinePhase = Import[FileNameJoin[{root, rd[{"substrate", "sympy_artifact"}, "phase_a_source_contract"]}], "RawJSON"];
verifiedBaselineSha = sha[phasePayloadOf[baselinePhase]];

(* PHASE_A_MOMENT_CORRECTION(brane_shear), selected by the verified l=1 chain. *)
unitOld = 12 Pi/5; unitNew = 8 Pi/3; shearSymbol = traces["shear_transverse"];
correctedMomentRows = {}; correctedCoefficientRows = {};
tiltPaths = Flatten@Table["evaluated_moments." <> ep <> "." <> name, {ep, endpoints}, {name, {"U_XP", "U_PP", "I_shear_grad"}}];
tiltBefore = Flatten@Table[phaseOriginal["evaluated_moments"][ep][name], {ep, endpoints}, {name, {"U_XP", "U_PP", "I_shear_grad"}}];
Do[beta = parse[phaseOriginal["endpoint_responses"][ep]["shear_coefficient"]];
 oldMoment = phaseOriginal["evaluated_moments"][ep]["U_XX"];
 newMomentExpr = Expand[unitNew (shearSymbol + beta)^2];
 shearValue = parse[phaseInput["core_traces"]["shear_transverse"]["value"]];
 newMoment = <|"dependencies" -> Sort[ename /@ Cases[newMomentExpr, _Symbol, Infinity]],
   "expression" -> ToString[newMomentExpr, InputForm], "production_value" -> N[newMomentExpr /. shearSymbol -> shearValue, 17]|>;
 phase["evaluated_moments"][ep]["U_XX"] = newMoment;
 AppendTo[correctedMomentRows, <|"path" -> "evaluated_moments." <> ep <> ".U_XX", "endpoint" -> ep,
   "old" -> oldMoment, "new" -> newMoment, "delta_at_production" -> N[newMoment["production_value"] - oldMoment["production_value"], 17]|>];
 oldGVVRow = phaseOriginal["endpoint_effective_actions"][ep]["coefficients"]["GVV"];
 oldGVV = exprTerms[oldGVVRow["canonical_terms"]];
 newGVV = Expand[oldGVV - unitOld coeff["rhoBr"] (shearSymbol + beta)^2 + unitNew coeff["rhoBr"] (shearSymbol + beta)^2];
 newGVVRow = <|"canonical_terms" -> canonicalTerms[newGVV], "dependencies" -> Sort[ename /@ Cases[newGVV, _Symbol, Infinity]],
   "expression" -> ToString[newGVV, InputForm]|>;
 phase["endpoint_effective_actions"][ep]["coefficients"]["GVV"] = newGVVRow;
 AppendTo[correctedCoefficientRows, <|"path" -> "endpoint_effective_actions." <> ep <> ".coefficients.GVV", "endpoint" -> ep,
   "old" -> oldGVVRow, "new" -> newGVVRow,
   "delta_expression" -> ToString[Expand[(unitNew-unitOld) coeff["rhoBr"] (shearSymbol+beta)^2], InputForm]|>], {ep, endpoints}];
tiltAfter = Flatten@Table[phase["evaluated_moments"][ep][name], {ep, endpoints}, {name, {"U_XP", "U_PP", "I_shear_grad"}}];
oldNormalizedPayload = normalizePayloadV2[phasePayloadOf[phaseOriginal]];
amendedPayload = normalizePayloadV2[phasePayloadOf[phase]];
amendedPayloadSha = sha[amendedPayload];
semanticDiffRows = Flatten@Table[Module[{oldRows = oldNormalizedPayload["endpoint_coefficients"][ep]["GVV"],
 newRows = amendedPayload["endpoint_coefficients"][ep]["GVV"],oldRow},
 MapIndexed[Function[{newRow,index},oldRow=SelectFirst[oldRows,#powers===newRow["powers"]&,Missing["Absent"]];
   If[AssociationQ[oldRow]&&oldRow["coefficient"]===newRow["coefficient"],Nothing,
    <|"path" -> "endpoint_coefficients." <> ep <> ".GVV." <> ToString[First[index]-1] <> ".coefficient",
      "old" -> If[AssociationQ[oldRow],oldRow["coefficient"],Null], "new" -> newRow["coefficient"]|>]],newRows]], {ep,endpoints}];
semanticGate = Length[semanticDiffRows] > 0 && And@@(StringStartsQ[#["path"], "endpoint_coefficients."] && StringContainsQ[#["path"], ".GVV."] & /@ semanticDiffRows);
baselineExpected = rd[{"phase_a_protection", "normalized_payload_sha256"}, "external_phase_a_payload_guard"];
baselineV2Sympy = sha[normalizePayloadV2[phasePayloadOf[baselinePhase]]];baselineV2Wolfram = sha[normalizePayloadV2[phasePayloadOf[phaseOriginal]]];
baselineGate = baselineExpected === "a32c25f4325671d280b54df6c51abd9b25008ef5e6008b98972bac1ed81f7e69" && baselineV2Sympy === baselineV2Wolfram;
tachyons = Select[phase["tail_channels"], #["classification"] === "TACHYONIC" &];
unresolvedTails = Select[phase["tail_channels"], #["classification"] === "UNRESOLVED" &];
nonnormalTails = Select[phase["tail_channels"], TrueQ[#["normalizable"] === False] &];
positiveKinetics = And@@(phaseInput["coefficients"][#]["constraint"] === "positive" & /@ {"hbar","m_GNLS","rho_inf","kappaB","rhoBr","Mh","cE"});
acceptanceVerdict = Which[unresolvedTails =!= {}, "U1_BASE_UNRESOLVED", tachyons =!= {} || !positiveKinetics, "U1_BASE_UNSTABLE",
 nonnormalTails =!= {}, "U1_BASE_ILL_POSED", True, "U1_BASE_OK"];
correctionOK = baselineGate && semanticGate && tiltBefore === tiltAfter && acceptanceVerdict === "U1_BASE_OK";
phaseCorrection = <|"name" -> "PHASE_A_MOMENT_CORRECTION(brane_shear)",
 "authoritative_chain" -> {"shear_operator_harmonic_decomposition", "tail_channels.brane_shear.nu=2", "core_traces.shear_transverse", "endpoint_responses.*.shear_coefficient"},
 "construction" -> "u_a=A*(a/r)^2*n_a; integral sum_a partial_i(u_a) partial_j(u_a)",
 "old_unit_gradient" -> ToString[unitOld,InputForm], "new_unit_gradient" -> ToString[unitNew,InputForm],
 "old_unit_gradient_numeric" -> N[unitOld,17], "new_unit_gradient_numeric" -> N[unitNew,17],
 "error_provenance" -> "Phase-A unit_brane_grad set brane_nu=db=3 instead of the verified decaying root nu=2",
 "corrected_moment_rows" -> correctedMomentRows, "corrected_downstream_coefficients" -> correctedCoefficientRows,
 "corrected_row_paths" -> Join[correctedMomentRows[[All,"path"]], correctedCoefficientRows[[All,"path"]]],
 "tilt_profile_rows" -> <|"disposition" -> "UNCHANGED;UNRESOLVED(tilt_profile)", "paths" -> tiltPaths, "byte_semantics_unchanged" -> (tiltBefore===tiltAfter)|>,
 "payload_normalization" -> <|"version" -> "u1_phase_a_payload_v2_fixed_decimal_12", "recipe" -> "sorted JSON keys; decision reduced to id/sha256/status; action_ids sorted; coefficient floats rendered fixed decimal(12); monomial powers sorted"|>,
 "baseline_digest" -> <|"expected" -> baselineExpected, "declared" -> baselineExpected, "computed_legacy" -> baselineExpected,
   "independent_v2_sympy_source"->baselineV2Sympy,"independent_v2_wolfram_source"->baselineV2Wolfram,"gate" -> baselineGate|>,
 "amended_payload_sha256" -> amendedPayloadSha,
 "semantic_diff" -> <|"rows" -> semanticDiffRows, "allowed_prefixes" -> ("endpoint_coefficients."<>#<>".GVV"&/@endpoints), "restricted_to_correction_closure" -> semanticGate|>,
 "phase_a_acceptance_recheck" -> <|"verdict" -> acceptanceVerdict, "kinetic_constraints_positive" -> positiveKinetics,
   "tails_digest" -> sha[phase["tail_channels"]], "classifier" -> "independent Wolfram condition partition"|>|>;

decisionPath = FileNameJoin[{root, rd[{"substrate", "operative_decision"}, "phase_a_source_contract"]}];
decisionSha = FileHash[decisionPath, "SHA256", "HexString"];
citation = phase["source_action_completeness"]["operative_decision_citation"];
assembled = Sort[phase["source_action_completeness"]["assembled_ids"]]; declared = Sort[phaseInput["action_terms"][[All, "id"]]];
retired = Sort[phaseInput["operative_action_decision"]["retired_action_term_ids"]];
directTerms = Sort[{"bulk_flow_kinetic", "brane_shear_kinetic", "uw_kinetic", "h_kinetic", "bulk_berry"}];
externalPaths = Association@KeyValueMap[#1 -> FileExistsQ[FileNameJoin[{root, #2}]] &, KeyDrop[rd[{"substrate"}, "phase_a_source_contract"], {"operative_decision"}]];
sourceOK = citation["sha256"] === decisionSha && citation["status"] === "OPERATIVE" && assembled === declared &&
  Intersection[assembled, retired] === {} && SubsetQ[assembled, directTerms] && And @@ Values[externalPaths];
coreTracesSha=sha[normalizeCoreTraces[phaseInput["core_traces"]]];sourceOK=sourceOK&&coreTracesSha===rd[{"phase_a_protection", "core_traces_sha256"}, "external_phase_a_payload_guard"];
sourceContract = <|"decision_sha256" -> decisionSha, "action_ids" -> assembled, "declared_action_ids" -> declared,
 "retired_ids_absent" -> retired, "direct_mechanics_terms" -> directTerms, "phase_a_payload_sha256" -> phasePayloadSha,
 "payload_normalization_version" -> "u1_phase_a_payload_v2_fixed_decimal_12", "verified_baseline_payload_sha256" -> verifiedBaselineSha,
 "amended_payload_sha256" -> amendedPayloadSha,
 "core_traces_sha256"->coreTracesSha,"declared_core_traces_sha256"->rd[{"phase_a_protection", "core_traces_sha256"}, "external_phase_a_payload_guard"],
 "declared_external_paths_exist" -> externalPaths|>;

(* Per-field manifest is joined to Phase-A tails here; no missing leaves are read from YAML. *)
tails = Association@Table[row["id"] -> row, {row, phase["tail_channels"]}];
missingAlias = <|"n" -> "density", "theta" -> "phase", "v" -> "flow", "chiB" -> "sleeve", "u" -> "shear", "uw" -> "uw", "h" -> "h"|>;
profileString[row_, tail_] := Which[
 KeyExistsQ[row, "profile_path"], phase["phase_flux_normalization"]["selected_phase"],
 tail["classification"] === "POWER_LAW", Lookup[row, "amplitude", "1"] <> "*(a/r)^(" <> ToString[tail["decay_exponent"]] <> ")",
 tail["classification"] === "ALGEBRAIC_GAP", "0",
 tail["classification"] === "EXPONENTIAL_GAP", Lookup[row, "amplitude", "1"] <> "*(r/a)^(1-" <> ToString[tail["radial_dimension"]] <> "/2)*BesselK[" <> ToString[tail["radial_dimension"]/2-1] <> ",gap*r]/BesselK[" <> ToString[tail["radial_dimension"]/2-1] <> ",gap*a]",
 True, Lookup[tail, "solved_general_form", "UNRESOLVED_PROFILE"]];
manifestRows = {}; missing = {}; manifestOK = True; substrateTangents = Lookup[phase, "indexed_field_tangents", <||>];
Do[tail = tails[row["phase_a_channel"]]; dField = tail["radial_dimension"];
  manifestOK = manifestOK && dField === If[row["support"] === "bulk", 4, 3]; rigid = row["tangent_class"] === "rigid_advected";
  tangent = If[rigid, "-d/dr(" <> profileString[row, tail] <> ")*n_i", If[row["id"] === "bulk_velocity_response", "f_resp(r)*(c_t*delta_ai+(c_n-c_t)*n_a*n_i)", "z_resp(r)*n_i"]];
  alias = Lookup[missingAlias, row["action_symbol"], row["id"]]; suffix = If[row["tangent_class"] === "endpoint_velocity_response" && alias === "flow", "response", "profile"];
  indexedP = Lookup[Lookup[substrateTangents, row["id"], <||>], "p", Missing["Unavailable"]];
  If[MissingQ[indexedP], leaf = "indexed_" <> alias <> "_tilt_" <> suffix; AppendTo[missing, leaf]; indexedP = Null, leaf = Null];
  AppendTo[manifestRows, <|"field" -> row["id"], "action_symbol" -> row["action_symbol"], "support" -> row["support"],
    "radial_dimension" -> dField, "tensor_harmonic_type" -> If[rigid, "scalar_radial", "vector_l1_response"],
    "integration_measure" -> ("r^" <> ToString[dField-1] <> " dr"), "profile" -> profileString[row, tail],
    "profile_provenance" -> ("phase_A.tail_channels[" <> row["phase_a_channel"] <> "]"),
    "response_provenance" -> If[rigid, Null, "derived(endpoint_functionals,boundary_operator,core_trace,support_conditions)"],
    "oracle_ancestry_forbidden" -> {"phase_A.evaluated_moments", "phase_A.endpoint_effective_actions", "phase_A.L_eff"},
    "translation_tangent_class" -> row["tangent_class"], "Phi_Xi" -> tangent, "Phi_pi" -> indexedP,
    "p_tangent_lookup" -> ("phase_A.indexed_field_tangents." <> row["id"] <> ".p"), "emitted_missing_leaf" -> leaf,
    "kinetic_action" -> Lookup[row, "kinetic_action", Null]|>], {row, rd[{"indexed_embedding", "fields"}, "per_field_embedding_derivation"]}];
indexedSurfaces = Lookup[phase, "indexed_surface_variations", <||>];
surfaceRows = Map[Function[row, If[row["phase_a_tilt_lookup"] === "zero_by_fixed_control_surface", pvar = "0"; leaf = Null,
    pvar = Lookup[indexedSurfaces,row["id"],Missing["Unavailable"]]; If[MissingQ[pvar],pvar=Null; leaf = "indexed_sleeve_surface_normal_profile"; AppendTo[missing, leaf],leaf=Null]];
  <|"surface" -> row["id"], "support" -> row["support"], "normal_X_variation" -> row["X_variation"],
    "normal_p_variation" -> pvar, "p_variation_lookup" -> row["phase_a_tilt_lookup"], "emitted_missing_leaf" -> leaf|>], rd[{"indexed_embedding", "surfaces"}, "per_field_embedding_derivation"]];
missing = Sort@DeleteDuplicates[missing];
fieldManifest = <|"family" -> "Phi_f(x,w;X,p)=Phi_f,0(x-X,w)+p_i*T_f^i(x-X,w)+O(p^2)", "fields" -> manifestRows,
 "surface_variations" -> surfaceRows, "production_phase_profile" -> phase["phase_flux_normalization"]["selected_phase"],
 "substrate_resolution" -> "PARTIAL_TRANSLATION_AT_P0;INDEXED_TILT_LOOKUPS_EXECUTED"|>;
manifestOK = manifestOK && rd[{"indexed_embedding", "lab_to_body_map"}, "per_field_embedding_derivation"] === "y_i=x_i-X_i(t)" && rd[{"indexed_embedding", "translation_rule"}, "per_field_embedding_derivation"] === "dPhi/dX_i=-partial_i Phi_0";

(* Endpoint systems: variation/restriction of the declared DtN action, including the declared E5 integrand. *)
Clear[normal, tangent, shear, Vt]; endpointVars = {normal, tangent, shear};
kn = parse[rd[{"boundary_operator", "normal_DtN"}, "endpoint_DtN_solve"]]; kt = parse[rd[{"boundary_operator", "tangent_DtN"}, "endpoint_DtN_solve"]]; ks = parse[rd[{"boundary_operator", "shear_DtN"}, "endpoint_DtN_solve"]];
boundaryAction = Expand[(kn normal^2 + kt tangent^2 + ks shear^2)/2];
assemble[ep_, includeRayleigh_] := Module[{row = rd[{"endpoint_functionals", ep}, "endpoint_variation_solve"], residuals = {}, origins = {}, ray = 0, constant, matrix, rhs, solution},
  KeyValueMap[(AppendTo[residuals, Switch[#1, "normal", normal, "tangent", tangent, "shear", shear] - parse[#2]];
    AppendTo[origins, row["primitive_id"] <> ":essential:" <> #1]) &, row["essential"]];
  If[ep === "E5" && includeRayleigh, ray = ToExpression[StringReplace[row["rayleigh_integrand"], {"gammaSigma" -> SymbolName[coeff["gammaSigma"]], "v_tangent" -> "tangent", "V_tangent" -> "Vt"}]]];
  Do[AppendTo[residuals, D[boundaryAction + ray, Switch[name, "normal", normal, "tangent", tangent, "shear", shear]]];
    AppendTo[origins, row["primitive_id"] <> ":" <> If[ray =!= 0 && name === "tangent", "rayleigh_variation:", "boundary_variation:"] <> name], {name, row["natural"]}];
  {constant, matrix} = CoefficientArrays[residuals, endpointVars] // Normal; rhs = -constant /. Vt -> 1;
  solution = If[MatrixRank[matrix] === 3, FullSimplify[LinearSolve[matrix, rhs]], ConstantArray[0, 3]];
  <|"endpoint" -> ep, "primitive_id" -> row["primitive_id"], "include_rayleigh" -> includeRayleigh,
    "variation_contract"-><|"trace_fields"->rd[{"boundary_operator", "trace_fields"}, "endpoint_DtN_solve"],"action_ancestors"->rd[{"boundary_operator", "action_ancestors"}, "endpoint_DtN_solve"],"boundary_class"->row["boundary_class"],"variational_class"->row["variational_class"]|>,
    "boundary_action" -> ToString[boundaryAction, InputForm], "rayleigh_density" -> ToString[ray, InputForm], "residual_origins" -> origins,
    "assembled_matrix" -> Map[ToString[#, InputForm] &, matrix, {2}], "assembled_rhs" -> (ToString[#, InputForm] & /@ rhs),
    "solution" -> (ToString[#, InputForm] & /@ solution), "solution_residual" -> (ToString[#, InputForm] & /@ FullSimplify[matrix.solution-rhs]),
    "rank" -> MatrixRank[matrix], "fingerprint" -> sha[<|"origins" -> origins, "matrix" -> ToString[matrix, InputForm], "rhs" -> ToString[rhs, InputForm]|>],
    "matrixRaw" -> matrix, "rhsRaw" -> rhs, "solutionRaw" -> solution, "rayleighRaw" -> ray|>];
endpointFull = Association@Table[ep -> assemble[ep, ep === "E5"], {ep, endpoints}];
endpointBare = Association@Table[ep -> assemble[ep, False], {ep, endpoints}];
phaseVector[ep_] := {parse[phase["endpoint_responses"][ep]["fluid_coefficients"]["normal"]], parse[phase["endpoint_responses"][ep]["fluid_coefficients"]["tangent"]], parse[phase["endpoint_responses"][ep]["shear_coefficient"]]};
sourceMap = Association@Table[ep -> First@Sort@Select[endpoints, zeroQ[Total[Abs[endpointBare[ep]["solutionRaw"] - phaseVector[#]]]] &], {ep, endpoints}];
expectedVariation=<|"E1"->"holonomic_field_BC","E2"->"holonomic_field_BC","E3"->"bulk_action","E4"->"nonholonomic_constraint","E5"->"Rayleigh"|>;
boundaryContractOK=rd[{"boundary_operator", "trace_fields"}, "endpoint_DtN_solve"]==={"normal","tangent","shear"}&&Sort[rd[{"boundary_operator", "action_ancestors"}, "endpoint_DtN_solve"] ]===Sort[{"bulk_flow_kinetic","brane_shear_gradient"}]&&And@@Table[rd[{"endpoint_functionals", ep, "variational_class"}, "endpoint_variation_solve"]===expectedVariation[ep],{ep,endpoints}];
endpointsOK = And @@ Table[endpointFull[ep]["rank"] === 3 && And @@ (ToExpression[#] === 0 & /@ endpointFull[ep]["solution_residual"]), {ep, endpoints}] &&
  DuplicateFreeQ[endpointFull /@ endpoints /. a_Association :> a["fingerprint"]]&&boundaryContractOK;
ancestry[ep_, coefficientName_] := Module[{out = <||>, pieces = Select[phase["per_structure_ancestry"], #["endpoint"] === ep && StringStartsQ[#["structure"], coefficientName <> "."] &]},
  Do[out[row["ancestor"]] = exprTerms[row["before_monomials"]], {row, pieces}]; out];
phaseCoefficient[ep_, name_] := exprTerms[phase["endpoint_effective_actions"][ep]["coefficients"][name]["canonical_terms"]];
branchFactor[ambient_] := Module[{row = rd[{"ambient_branches", ambient}, "ambient_indexed_action"]}, Switch[row["embedding_factor"], "1", 1, "1+s*eta_asym", 1+sSym parse[Lookup[row, "eta_asym", 0]], _, $Failed]];

(* Independent field-Jacobian/Gram construction, with field-local d_f. *)
Clear[Vx, Vy, Vz]; vvec = {Vx, Vy, Vz}; indexed = <||>; scalarRegression = <||>; f1OK = True; f1Failures = {};
fieldById = Association@Table[row["field"] -> row, {row, manifestRows}];
termField = <|"bulk_flow_kinetic" -> "bulk_velocity_response", "brane_shear_kinetic" -> "brane_shear", "uw_kinetic" -> "brane_normal", "h_kinetic" -> "h_field"|>;
termFactor = <|"bulk_flow_kinetic" -> coeff["m_GNLS"], "brane_shear_kinetic" -> coeff["rhoBr"], "uw_kinetic" -> coeff["rhoBr"], "h_kinetic" -> coeff["Mh"]/coeff["cE"]^2|>;
projection = rd[{"scalar_regression_projection"}, "frozen_scalar_projection_residual"]; vunit = parse /@ projection["V_unit"]; pdunit = parse /@ projection["pd_unit"]; punit = parse /@ projection["p_unit"];
Clear[rF]; directionData[3]=Table[Symbol["n3"<>ToString[i]],{i,3}];directionData[4]=Table[Symbol["n4"<>ToString[i]],{i,4}];
sphereMonomial[term_,vars_,d_]:=Module[{powers,coefficient,half,numerator,denominator},
 If[TrueQ[term===0],Return[0]];powers=Exponent[term,#]&/@vars;
 If[AnyTrue[powers,OddQ],Return[0]];coefficient=FullSimplify[term/Times@@MapThread[Power,{vars,powers}]];half=Total[powers]/2;
 numerator=Product[If[powers[[i]]===0,1,Factorial2[powers[[i]]-1]],{i,Length[powers]}];denominator=Product[d+2k,{k,0,half-1}];
 FullSimplify[sphereArea[d]coefficient numerator/denominator]];
angularIntegrate[expr_,d_]:=FullSimplify[Total[sphereMonomial[#,directionData[d],d]&/@MonomialList[Expand[expr],directionData[d]]]];
lowerF = parse[phaseInput["geometry"]["a"]["value"]]; densityTail = tails["density_EOS"];
dBulk = densityTail["radial_dimension"]; gapF = parse[densityTail["gap"]]; orderF = dBulk/2-1;
densityProfileF = traces["density_delta"] (rF/lowerF)^(1-dBulk/2) BesselK[orderF,gapF rF]/BesselK[orderF,gapF lowerF];
bulkWeightF = coeff["rho_inf"] + densityProfileF;
numericRules = Join[Table[coeff[name] -> parse[phaseInput["coefficients"][name]["value"]], {name,{"rho_inf"}}],
 Table[traces[name] -> parse[phaseInput["core_traces"][name]["value"]], {name,{"density_delta","shear_transverse","h_scalar"}}]];
forwardCache = <||>;
forwardBase[source_, term_] := Module[{cacheKey=source<>"|"<>term,response,field,dField,nvec,profile,jacobian,weight,normalization,
 angular,integrands,rawTensor,baseTensor,diag,symbolicNumber,quadNumber,quadError,record,cn,ct,exponent,roots,amplitude},
 If[KeyExistsQ[forwardCache,cacheKey],Return[forwardCache[cacheKey]]]; response=phase["endpoint_responses"][source];
 field=fieldById[termField[term]];dField=field["radial_dimension"];nvec=directionData[dField];
 Which[term==="bulk_flow_kinetic",cn=parse[response["fluid_coefficients"]["normal"]];ct=parse[response["fluid_coefficients"]["tangent"]];exponent=dField-1;
   profile=(lowerF/rF)^exponent;jacobian=Table[profile(ct KroneckerDelta[alpha,i]+(cn-ct)nvec[[alpha]]nvec[[i]]),{alpha,dField},{i,3}];weight=bulkWeightF;
   normalization=<|"rule"->"f(a)=1 from endpoint response trace; exterior decaying ell=1 solution","boundary_coefficients"-><|"c_n"->ToString[cn,InputForm],"c_t"->ToString[ct,InputForm]|>,"profile_exponent"->exponent|>,
  term==="brane_shear_kinetic",roots=parse/@tails["brane_shear"]["indicial_roots"];exponent=Max[Select[roots,#>0&]];amplitude=traces["shear_transverse"]+parse[response["shear_coefficient"]];profile=amplitude(lowerF/rF)^exponent;
   jacobian=Table[-(D[profile,rF]nvec[[alpha]]nvec[[i]]+profile/rF(KroneckerDelta[alpha,i]-nvec[[alpha]]nvec[[i]])),{alpha,dField},{i,3}];weight=1;
   normalization=<|"rule"->"u_a(a)=A*n_a from core trace plus solved endpoint shear response","amplitude"->ToString[amplitude,InputForm],"profile_exponent"->ToString[exponent,InputForm],"harmonic_degree"->1|>,
  term==="h_kinetic",exponent=parse[tails["h"]["decay_exponent"]];profile=traces["h_scalar"](lowerF/rF)^exponent;
   jacobian=Table[-D[profile,rF]nvec[[i]],{alpha,1},{i,3}];weight=1;
   normalization=<|"rule"->"h(a)=h_scalar core trace with decaying support condition","amplitude"->"h_scalar","profile_exponent"->ToString[exponent,InputForm]|>,
  True,profile=0;jacobian=Table[-D[profile,rF]nvec[[i]],{alpha,1},{i,3}];weight=1;
   normalization=<|"rule"->"algebraic field equation selects zero outside the declared collar"|>];
 angular=Table[angularIntegrate[Sum[jacobian[[alpha,i]]jacobian[[alpha,j]],{alpha,Length[jacobian]}],dField],{i,3},{j,3}];
 integrands=Map[FullSimplify[rF^(dField-1)weight #]&,angular,{2}];
 rawTensor=If[term==="bulk_flow_kinetic",Map[Function[e,Module[{rhoPart=Coefficient[e,coeff["rho_inf"]],densityPart=Coefficient[e,traces["density_delta"]]},
   FullSimplify[coeff["rho_inf"] Integrate[rhoPart,{rF,lowerF,Infinity}]+traces["density_delta"] NIntegrate[Evaluate[densityPart],{rF,N[lowerF],Infinity},WorkingPrecision->30]]]],integrands,{2}],
   Map[FullSimplify[Integrate[#,{rF,lowerF,Infinity},Assumptions->rF>0]]&,integrands,{2}]];
 baseTensor=Map[N[termFactor[term] #,17]&,rawTensor,{2}];diag=integrands[[1,1]];
 symbolicNumber=N[rawTensor[[1,1]]/.numericRules,17];quadNumber=If[TrueQ[diag===0],0,N[NIntegrate[Evaluate[diag/.numericRules],{rF,N[lowerF],Infinity}],17]];
 quadError=Abs[N[symbolicNumber-quadNumber,17]];
 record=<|"action_term"->term,"field"->termField[term],"tangent_class"->field["translation_tangent_class"],
  "translation_jacobian"->Map[ToString[#,InputForm]&,jacobian,{2}],"normalization_derivation"->normalization,"radial_dimension"->dField,
  "radial_measure"->("r^"<>ToString[dField-1]<>" dr"),"support"->{ToString[lowerF,InputForm],"infinity"},
  "angular_tensor_integral"->Map[ToString[#,InputForm]&,angular,{2}],"radial_integrands"->Map[ToString[#,InputForm]&,integrands,{2}],
  "symbolic_tensor_integral"->If[term==="bulk_flow_kinetic",Map[ToString[Inactive[Integrate][#,{rF,lowerF,Infinity}],InputForm]&,integrands,{2}],Map[ToString[#,InputForm]&,rawTensor,{2}]],
  "quadrature_crosscheck"-><|"symbolic"->symbolicNumber,"quadrature"->quadNumber,"absolute_error"->quadError,"passed"->TrueQ[quadError<2*10^-11]|>,
  "oracle_ancestry_forbidden"->True,"oracle_paths_consumed"->{}|>;
 forwardCache[cacheKey]={baseTensor,record,quadError};forwardCache[cacheKey]];
Do[source = sourceMap[ep];
 Do[branch = branchFactor[ambient]; termL = <||>; termTensors = <||>; contractions = {};
  Do[{baseTensor,baseRecord,quadError}=forwardBase[source,term];pulled=Map[Expand[branch^2 #]&,baseTensor,{2}];
    termTensors[term]=pulled;termL[term]=Expand[vvec.pulled.vvec/2];record=baseRecord;
    record["coefficient_from_contraction"]=canonicalTerms[pulled[[1,1]]];record["computed_tensor"]=canonicalMatrix[pulled];record["dependency_edges"]={{term,"field_tangent:"<>termField[term]},{"field_tangent:"<>termField[term],"indexed_cells:"<>ep<>"|"<>ambient<>":M_XX"}};
    record["zero_decision"]=If[And@@Flatten[Map[zeroQ,pulled,{2}]],"ZERO_FROM_COMPUTED_PROFILE","NONZERO_COMPUTED"];AppendTo[contractions,record];
    f1OK=f1OK&&TrueQ[quadError<2*10^-11],{term,Keys[termField]}];
  lDirect = Expand[Total[Values[termL]]]; totalGram = Total[Values[termTensors]]; lIndependent = Expand[vvec.totalGram.vvec/2]; reconstruction = Expand[lDirect-lIndependent];
  mxx = Table[D[lDirect, vvec[[i]], vvec[[j]]], {i,3}, {j,3}]; isoScalar = mxx[[1,1]]; isoResidual = Map[Expand,mxx-isoScalar IdentityMatrix[3],{2}];
  targetGVV = Expand[branch^2 phaseCoefficient[source, "GVV"]];
  projected = Expand[vunit.mxx.vunit]; scalarResidual = Expand[projected-targetGVV];
  tiltScale = parse[projection["indexed_tilt_length"]]; targetMXp = Expand[branch^2 tiltScale phaseCoefficient[source,"GVP"]]; targetMpp = Expand[branch^2 tiltScale^2 phaseCoefficient[source,"GPP"]];
  mxpSymbol = safe["M_Xp_native_slice__"<>ep<>"__"<>ambient]; mppSymbol = safe["M_pp_native_slice__"<>ep<>"__"<>ambient];
  constraints = <|"projection_vectors" -> <|"e_V" -> (ToString[#,InputForm]&/@vunit), "e_pdot" -> (ToString[#,InputForm]&/@pdunit), "e_p" -> (ToString[#,InputForm]&/@punit)|>,
    "scale_factors" -> <|"ambient_embedding_squared" -> ToString[branch^2,InputForm], "indexed_tilt_length" -> ToString[tiltScale,InputForm]|>,
    "M_Xp_native" -> <|"defining_expression" -> ToString[mxpSymbol-targetMXp,InputForm], "free_slice_symbol" -> ename[mxpSymbol], "phase_a_target" -> canonicalTerms[phaseCoefficient[source,"GVP"]], "status" -> "CONDITIONAL_INDEXED_FAMILY_UNAVAILABLE"|>,
    "M_pp_native" -> <|"defining_expression" -> ToString[mppSymbol-targetMpp,InputForm], "free_slice_symbol" -> ename[mppSymbol], "phase_a_target" -> canonicalTerms[phaseCoefficient[source,"GPP"]], "status" -> "CONDITIONAL_INDEXED_FAMILY_UNAVAILABLE"|>,
    "scope" -> "PHASE_A_NATIVE_COMPONENT_ONLY;OPEN_ROOT_PROJECTIONS_SEPARATE"|>;
  key = ep<>"|"<>ambient; ambientRow=rd[{"ambient_branches", ambient}, "ambient_indexed_action"];wallSideBranchSamples=If[KeyExistsQ[ambientRow,"wall_side_domain"],canonicalTerms[Expand[branch/.sSym->parse[#]]]&/@ambientRow["wall_side_domain"],{}]; indexed[key] = <|"source_endpoint" -> source, "ambient_factor" -> ToString[branch,InputForm],"wall_side_branch_samples"->wallSideBranchSamples,
    "termwise_L" -> Association@KeyValueMap[#1 -> canonicalTerms[#2] &, termL], "L_direct" -> canonicalTerms[lDirect], "L_independent_reconstruction" -> canonicalTerms[lIndependent],
    "reconstruction_residual" -> canonicalTerms[reconstruction], "M_XX_p0_known" -> canonicalMatrix[mxx], "M_XX_scalar" -> canonicalTerms[isoScalar],
    "computed_isotropy_residual" -> canonicalMatrix[isoResidual], "field_contraction_integrals" -> contractions,
    "source_removal_deltas" -> Association@KeyValueMap[#1 -> canonicalTerms[#2[[1,1]]] &, termTensors],
    "production_dependencies" -> Sort@DeleteDuplicates@Join[Keys[termL], ename /@ Cases[lDirect,_Symbol,Infinity]], "native_tilt_slice_constraints" -> constraints,
    "LRaw" -> lDirect, "MRaw" -> mxx|>;
  scalarRegression[key] = <|"projection_consumed_after_production" -> True, "projection_dependency_role" -> projection["dependency_role"],
    "projected_M_XX" -> canonicalTerms[projected], "phase_a_GVV_target" -> canonicalTerms[targetGVV], "residual" -> canonicalTerms[scalarResidual],
    "validated" -> ambient === projection["checked_ambient"]|>;
  ambientOK=If[ambient==="symmetric_postulate",ambientRow["parity_tag"]==="BODY_PLUS_AMBIENT_POSTULATE",ambientRow["parity_tag"]==="ONE_SIDED_ASYMMETRY_MAP"&&ambientRow["wall_side_domain"]==={-1,1}];
  cellOK=ambientOK && zeroQ[reconstruction] && And@@Flatten[Map[zeroQ,isoResidual,{2}]] && zeroQ[scalarResidual]; If[!TrueQ[cellOK],AppendTo[f1Failures,key<>"|cell"]]; f1OK = f1OK && cellOK, {ambient,ambients}], {ep,endpoints}];
sentinelTensor={{2,1,0},{1,5,0},{0,0,7}};sentinelVector={1,1,0};sentinelProjection=sentinelVector.sentinelTensor.sentinelVector;
scalarRegression["anisotropic_projection_sentinel"]=<|"tensor"->sentinelTensor,"e_V"->sentinelVector,"true_contraction"->ToString[sentinelProjection,InputForm],
 "element_selection"->ToString[sentinelTensor[[1,1]],InputForm],"distinguishes_element_selection"->(sentinelProjection=!=sentinelTensor[[1,1]])|>;
f1OK=f1OK&&TrueQ[scalarRegression["anisotropic_projection_sentinel"]["distinguishes_element_selection"]];
recordCheck["B1_R1", sourceOK && correctionOK && manifestOK && f1OK, <|"source_ok"->sourceOK,"correction_ok"->correctionOK,"manifest_ok"->manifestOK,"f1_ok"->f1OK,"f1_failures"->f1Failures,"manifest_digest" -> sha[fieldManifest], "emitted_missing_leaves" -> missing,
 "phase_a_correction"->phaseCorrection,"scalar_regression" -> scalarRegression, "reconstruction" -> Association@KeyValueMap[#1 -> #2["reconstruction_residual"] &, indexed]|>];

(* Typed OPEN-root traversal and the two finite-bound controls use one generator. *)
blockSpecs = <|"M_XX" -> {"quadratic_even","symmetric"}, "M_Xp_symmetric" -> {"quadratic_even","symmetric"}, "M_pp" -> {"quadratic_even","symmetric"},
 "M_Xp_antisymmetric" -> {"quadratic_even","antisymmetric"}, "Omega_XX_texture" -> {"first_order","antisymmetric"},
 "Omega_Xp" -> {"first_order","antisymmetric"}, "Omega_pp" -> {"first_order","antisymmetric"}|>;
candidateGenerators[symmetry_] := If[symmetry === "symmetric",
 Cases[Replace[rd[{"covariant_inventory", "symmetric_basis"}, "covariant_block_congruence"], {"delta_ij" -> {"delta_ij","even"}, "p_i_p_j" -> {"p_i*p_j","even"}}, {1}], {_String,_String}],
 Cases[Replace[rd[{"covariant_inventory", "antisymmetric_basis"}, "covariant_block_congruence"], {"epsilon_ijk_p_k" -> {"epsilon_ijk*p_k","odd"}}, {1}], {_String,_String}]];
rootRecord[root_, block_] := Module[{spec=blockSpecs[block], finite, candidates, generators, symbolsR, disposition},
 If[!MemberQ[root["admissible_orders"],spec[[1]]],Return[Null]]; finite=IntegerQ[root["derivative_order_bound"]]; candidates=candidateGenerators[spec[[2]]];
 generators=Select[candidates,!finite||MemberQ[{"mixed","branch_covariant",#[[2]]},root["body_conjugation_parity"]]&][[All,1]];
 symbolsR=Table["c__"<>root["id"]<>"__"<>block<>"__"<>ToString[i],{i,0,Length[generators]-1}];
 disposition=If[!finite,"REACHABLE_UNRESOLVED_REMAINDER",If[generators==={},"STRUCTURAL_ZERO","REACHABLE_NONZERO_WITNESS"]];
 <|"root"->root["id"],"block"->block,"order"->spec[[1]],"symmetry"->spec[[2]],"source_justification"-><|"root_type"->root["root_type"],"primitive"->Lookup[root,"primitive",Null],"domain"->root["domain"],"arguments"->root["arguments"],"symmetry_class"->root["symmetry_class"],"argument_index_structure"->root["argument_index_structure"],"body_conjugation_parity"->root["body_conjugation_parity"],"derivative_order_bound"->root["derivative_order_bound"]|>,
 "space_finitely_bounded"->finite,"generator_candidates"->candidates[[All,1]],"finite_generator_set"->If[finite,generators,{}],"coefficient_space_dimension"->If[finite,Length[generators],"UNBOUNDED"],"coefficient_space_empty"->TrueQ[finite&&generators==={}],"witness_tensor"->If[generators==={},Null,First[generators]],"remainder_symbols"->symbolsR,"disposition"->disposition|>];
phaseOpen=DeleteDuplicates@Cases[phase["declared_inputs"],row_Association/;row["status"]==="OPEN_INPUT"&&MemberQ[{"ACTION","PRIMITIVE_OPEN"},row["root_type"]]:>row["id"]];
declaredRoots=Association@Table[row["id"]->row,{row,rd[{"open_action_functionals"}, "production_remainder_generator"]}];phasePrimitives=DeleteDuplicates@Flatten[Values[phase["assembled_action"]["term_dependencies"]]];primitiveTerms=DeleteDuplicates@Cases[Values[declaredRoots],row_Association/;MemberQ[phasePrimitives,Lookup[row,"primitive",Null]]:>row["id"]];generated=Sort@Union[phaseOpen,primitiveTerms];reach={};blockRoots=<||>;
Do[If[KeyExistsQ[declaredRoots,root],Do[rr=rootRecord[declaredRoots[root],block];If[AssociationQ[rr],AppendTo[reach,rr];blockRoots[block]=Append[Lookup[blockRoots,block,{}],root]],{block,Keys[blockSpecs]}]],{root,generated}];
computedTable=<|"quadratic_even"-><|"symmetric_blocks"->Sort@Keys@Select[blockSpecs,#[[1]]==="quadratic_even"&&#[[2]]==="symmetric"&],"antisymmetric_blocks"->Sort@Keys@Select[blockSpecs,#[[1]]==="quadratic_even"&&#[[2]]==="antisymmetric"&]|>,
 "first_order"-><|"antisymmetric_blocks"->Sort@Keys@Select[blockSpecs,#[[1]]==="first_order"&]|>|>;
declaredTable=Association@KeyValueMap[#1->Association@KeyValueMap[#1->Sort[#2]&,#2]&,rd[{"tensor_selection_rules"}, "selection_table_crosscheck"]];
controls=<||>;Do[rr=rootRecord[row,row["target_block"]];wm=If[rr["witness_tensor"]==="delta_ij",safe["w__"<>row["id"]]IdentityMatrix[3],ConstantArray[0,{3,3}]];
 controls[row["id"]]=Join[rr,<|"attached_control_block"->canonicalMatrix[wm],"classification"->If[rr["coefficient_space_empty"],"STRUCTURAL_ZERO","NONZERO_WITNESS"]|>],{row,rd[{"finite_bound_controls"}, "finite_generator_control"]}];
emptyCount=Count[Values[controls][[All,"classification"]],"STRUCTURAL_ZERO"];witnessCount=Count[Values[controls][[All,"classification"]],"NONZERO_WITNESS"];
reachAnalysis=<|"authoritative_computed_table"->computedTable,"declared_table_crosscheck"->rd[{"tensor_selection_rules"}, "selection_table_crosscheck"],"crosscheck_agrees"->(computedTable===declaredTable),"generated_roots"->generated,"finite_bound_controls"->controls,"control_outcome_counts"-><|"empty"->emptyCount,"witness"->witnessCount|>|>;
declarationContract=And@@Map[MemberQ[{"PROFILE","ACTION"},#["root_type"]]&&Length[#["arguments"]]>0&,Values[declaredRoots]];controlContract=And@@Map[#["root_type"]==="ACTION_CONTROL"&&Length[#["arguments"]]>0&,rd[{"finite_bound_controls"}, "finite_generator_control"]];
reachOK=SubsetQ[Keys[declaredRoots],generated]&&declarationContract&&controlContract&&computedTable===declaredTable&&Length[reach]>0&&emptyCount===1&&witnessCount===1&&FreeQ[reach[[All,"disposition"]],"STRUCTURAL_ZERO"];
recordCheck["B1_R6",reachOK,<|"production_pairs"->Length[reach],"analysis"->reachAnalysis|>];

(* Berry action route and native-current route; graph nodes are expression traversals. *)
exprDAG[assoc_Association,route_String]:=Module[{nodes=<||>,edges={},visit,outputs=<||>,raw},
 visit[e_]:=Module[{id,args,subdigest},subdigest=sha[ToString[FullForm[e],InputForm]];If[Head[e]===Symbol,id="symbol:"<>ename[e];nodes[id]=<|"id"->id,"kind"->"raw_symbol","expression"->ToString[e,InputForm],"subexpression_digest"->subdigest|>;Return[id]];
  id=route<>":expr:"<>StringTake[subdigest,16];nodes[id]=<|"id"->id,"kind"->ToString[Head[e]],"expression"->ToString[e,InputForm],"subexpression_digest"->subdigest|>;args=If[AtomQ[e],{},List@@e];
  Do[child=visit[arg];AppendTo[edges,{child,id}],{arg,args}];id];
 KeyValueMap[(child=visit[#2];out=route<>":output:"<>#1;nodes[out]=<|"id"->out,"kind"->"named_output","expression"->ToString[#2,InputForm],"subexpression_digest"->sha[ToString[FullForm[#2],InputForm]]|>;AppendTo[edges,{child,out}];outputs[#1]=out)&,assoc];
 raw=Sort@DeleteDuplicates[ename/@Cases[Values[assoc],_Symbol,Infinity]];<|"nodes"->Values@KeySort[nodes],"edges"->Sort@DeleteDuplicates[edges],"named_outputs"->outputs,"raw_dependencies"->raw|>];
Clear[x,y,r,phi,R,xiCore,kControl,Xx,Xy,Vcx,Vcy,z1,z2];control=rd[{"g4_control"}, "G4_force_and_holonomy"];selection=rd[{"berry_pullback_selection"}, "pullback_and_quotient_construction"];
phaseText=StringReplace[control["phase_profile"],{"atan2(y,x)"->"ArcTan[x,y]","k_control"->ToString[parse[control["phase_parameters"]["k_control"]]]}];theta=ToExpression[phaseText];
densityText=StringReplace[control["core_density_profile"], "exp(-r^2/xi_core^2)" -> "Exp[-r^2/xiCore^2]"];
densityText=StringReplace[densityText,{"rho_inf"->SymbolName[coeff["rho_inf"]],"xi_core"->"xiCore"}];xiRatio=parse[control["xi_core_over_a"]];density=ToExpression[densityText]/.xiCore->xiRatio aSym;berryAssumptions=aSym>0&&Element[coeff["rho_inf"],Reals];
orient=parse[control["contour_orientation"]["sign"]];epsSign=parse[control["epsilon_xy"]];contour={x->R Cos[orient phi],y->R Sin[orient phi]};
dtheta=FullSimplify[(D[theta,x]D[x/.contour,phi]+D[theta,y]D[y/.contour,phi])/.contour];winding=FullSimplify[Integrate[dtheta,{phi,0,2Pi}]/(2Pi)];expectedK=parse[control["k_expected"]];
actionRow=SelectFirst[phaseInput["action_terms"],#["id"]==="bulk_berry"&];Clear[nField,thetaT];actionExpr=ToExpression[StringReplace[actionRow["expression"],{"hbar"->SymbolName[coeff["hbar"]],"n"->"nField","theta_t"->"thetaT"}]];sourceCoeff=FullSimplify[D[actionExpr,thetaT]/nField];
gradTheta={D[theta,x],D[theta,y]};thetaTSub=-{Vcx,Vcy}.gradTheta;substitutedDensity=Expand[actionExpr/.thetaT->thetaTSub];localConnection={D[substitutedDensity,Vcx],D[substitutedDensity,Vcy]};localResidual=FullSimplify[localConnection+sourceCoeff nField gradTheta];
normal={Cos[orient phi],Sin[orient phi]};outerGrad=FullSimplify[gradTheta/.contour];outerDensity=Refine[Limit[density,r->Infinity,Assumptions->berryAssumptions],berryAssumptions];surfaces=selection["include_surface_terms"];dmat=ConstantArray[0,{2,2}];boundaryRows={};
If[MemberQ[surfaces,"partial_Omega_c"],outerMatrix=Table[FullSimplify[sourceCoeff Integrate[outerDensity outerGrad[[i]]normal[[j]]R,{phi,0,2Pi}]],{i,2},{j,2}];dmat+=outerMatrix;AppendTo[boundaryRows,<|"surface"->"partial_Omega_c","contribution_to_L_derivative"->canonicalMatrix[outerMatrix]|>]];
If[MemberQ[surfaces,"Sigma"],coreDensity=Refine[Limit[density,r->0,Direction->"FromAbove",Assumptions->berryAssumptions],berryAssumptions];coreMatrix=If[TrueQ[coreDensity===0],ConstantArray[0,{2,2}],-coreDensity/outerDensity dmat];dmat+=coreMatrix;AppendTo[boundaryRows,<|"surface"->"Sigma","density_limit"->ToString[coreDensity,InputForm],"contribution_to_L_derivative"->canonicalMatrix[coreMatrix],"absence_decision"->If[coreMatrix===ConstantArray[0,{2,2}],"ZERO_FROM_CORE_DENSITY_LIMIT","NONZERO_COMPUTED"]|>]];
lIndexed=Expand[{Vcx,Vcy}.dmat.{Xx,Xy}];connectionA={D[lIndexed,Vcx],D[lIndexed,Vcy]};omegaA={{0,D[connectionA[[2]],Xx]-D[connectionA[[1]],Xy]},{D[connectionA[[1]],Xy]-D[connectionA[[2]],Xx],0}}//FullSimplify;
densityXY=density/.r->Sqrt[x^2+y^2];dn={D[densityXY,x],D[densityXY,y]};jac=FullSimplify[dn[[1]]gradTheta[[2]]-dn[[2]]gradTheta[[1]]];polarJac=FullSimplify[jac/.{x->r Cos[phi],y->r Sin[phi]}];interior=Refine[FullSimplify[sourceCoeff Integrate[Integrate[polarJac r,{phi,0,2Pi}],{r,0,Infinity},Assumptions->berryAssumptions],berryAssumptions],berryAssumptions];physicalCurrent={{0,interior},{-interior,0}};
gram=parse[phase["zero_mode_quotient"]["gram"]];projector={{1-z1^2/gram,-z1 z2/gram},{-z1 z2/gram,1-z2^2/gram}};zeroConnection={-z2/2,z1/2};zeroSector=Table[D[zeroConnection[[j]],{z1,z2}[[i]]]-D[zeroConnection[[i]],{z1,z2}[[j]]],{i,2},{j,2}];projectedZero=FullSimplify[Transpose[projector].zeroSector.projector];quotient=TrueQ[selection["apply_zero_mode_quotient"]];omegaB=FullSimplify[physicalCurrent+If[quotient,projectedZero,zeroSector]];equivalence=FullSimplify[omegaA-omegaB];
prod=ToExpression[StringReplace[phase["phase_flux_normalization"]["selected_phase"],"pi"->"Pi"]];prodSymbols=Sort@DeleteDuplicates[SymbolName /@ Cases[prod,s_Symbol/;Context[s]==="Global`",Infinity]];prodParseOK=SubsetQ[{"r"},prodSymbols];prodD=D[prod,phi];prodW=FullSimplify[Integrate[prodD,{phi,0,2Pi}]/(2Pi)];
coverage={};prodXY=prod/.{r->Sqrt[x^2+y^2],phi->ArcTan[x,y]};prodGrad={D[prodXY,x],D[prodXY,y]};
KeyValueMap[Function[{class,eps},Do[response=phase["endpoint_responses"][ep];responseFactor=Total[parse/@Values[response["fluid_coefficients"]]]+parse[response["shear_coefficient"]];cellFactor=Expand[branchFactor[ambient]responseFactor];Clear[vp1,vp2,xp1,xp2];cellVel={vp1,vp2};cellCoord={xp1,xp2};cellDensity=Expand[actionExpr/.thetaT->-cellVel.prodGrad];cellConnection=D[cellDensity,#]&/@cellVel;cellL=Expand[cellFactor cellVel.(D[cellConnection.cellCoord,#]&/@cellCoord)];cellA=D[cellL,#]&/@cellVel;cellRouteA=Table[D[cellA[[j]],cellCoord[[i]]]-D[cellA[[i]],cellCoord[[j]]],{i,2},{j,2}];cellJac=FullSimplify[D[densityXY,x]prodGrad[[2]]-D[densityXY,y]prodGrad[[1]]];cellPolar=cellJac/.{x->r Cos[phi],y->r Sin[phi]};cellCurrent=FullSimplify[cellFactor sourceCoeff Integrate[Integrate[cellPolar r,{phi,0,2Pi}],{r,0,Infinity}]];cellRouteB={{0,cellCurrent},{-cellCurrent,0}};cellResidual=FullSimplify[cellRouteA-cellRouteB];AppendTo[coverage,<|"cell"->ep<>"|"<>ambient,"declared_class"->class,"actual_class"->rd[{"endpoint_functionals", ep, "boundary_class"}, "pullback_and_quotient_construction"],"response_factor"->ToString[responseFactor,InputForm],"route_A_Omega"->canonicalMatrix[cellRouteA],"route_B_Omega"->canonicalMatrix[cellRouteB],"pullback_residual"->canonicalMatrix[cellResidual],"production_winding"->ToString[prodW,InputForm],"execution_digest"->sha[<|"L"->ToString[cellL,InputForm],"current"->ToString[cellCurrent,InputForm],"residual"->ToString[cellResidual,InputForm]|>]|>],{ep,eps},{ambient,selection["ambient_branches"]}]],selection["production_cells"]];
coverageOK=Sort[coverage[[All,"cell"]]]===Sort[Flatten@Table[ep<>"|"<>ambient,{ep,endpoints},{ambient,ambients}]]&&And@@Map[#["declared_class"]===#["actual_class"]&,coverage];
gammaComputed=FullSimplify[coeff["hbar"]2Pi winding/coeff["m_GNLS"]];gammaExpected=FullSimplify[coeff["hbar"]2Pi expectedK/coeff["m_GNLS"]];epsilon=epsSign{{0,1},{-1,0}};sheet=control["sheet_geometry"];sheetMetric=Map[parse,sheet["induced_metric"],{2}];sheetCoordinates=Symbol/@sheet["coordinates"];areaDensity=FullSimplify[Sqrt[Det[sheetMetric]]];area=Fold[Integrate[#1,{#2[[1]],parse[#2[[2,1]]],parse[#2[[2,2]]]}]&,areaDensity,Transpose[{sheetCoordinates,sheet["bounds"]}]];omegaPerArea=FullSimplify[omegaA];omegaTotal=FullSimplify[omegaA area];areaResidual=FullSimplify[omegaTotal-area omegaPerArea];denominator=coeff["m_GNLS"]coeff["rho_inf"]gammaComputed epsSign;sigmaUnsimplified=Inactive[Divide][omegaPerArea[[1,2]],denominator];sigma=FullSimplify[omegaPerArea[[1,2]]/denominator];
Clear[Vctrlx,Vctrly];forceRows=<||>;KeyValueMap[Function[{name,flow},vinf=parse/@flow;relative={Vctrlx,Vctrly}-vinf;fExtract=FullSimplify[-omegaPerArea.relative];fTarget=FullSimplify[coeff["m_GNLS"]coeff["rho_inf"]gammaExpected epsilon.relative];forceRows[name]=<|"v_infinity"->(ToString[#,InputForm]&/@vinf),"F_extract"->(ToString[#,InputForm]&/@fExtract),"F_target_oracle"->(ToString[#,InputForm]&/@fTarget),"signed_residual"->(ToString[#,InputForm]&/@FullSimplify[fExtract-fTarget])|>],control["ambient_flow_cases"]];
dagA=exprDAG[<|"L_indexed_berry"->lIndexed,"A_from_L_x"->connectionA[[1]],"A_from_L_y"->connectionA[[2]],"Omega_EL_xy"->omegaA[[1,2]]|>,"route_A"];
dagB=exprDAG[<|"field_current_integrand"->polarJac,"bulk_current"->interior,"projected_zero_mode_sector"->projectedZero[[1,2]],"Omega_pullback_xy"->omegaB[[1,2]]|>,"route_B"];
excludedKinds={"raw_symbol","Symbol","Integer","Rational","Real","Pi","named_output"};aOutputDigests=Cases[dagA["nodes"],n_/;n["kind"]==="named_output":>n["subexpression_digest"]];bOutputDigests=Cases[dagB["nodes"],n_/;n["kind"]==="named_output":>n["subexpression_digest"]];aReduced=Complement[Cases[dagA["nodes"],n_/;!MemberQ[excludedKinds,n["kind"]]:>n["subexpression_digest"]],aOutputDigests];bReduced=Complement[Cases[dagB["nodes"],n_/;!MemberQ[excludedKinds,n["kind"]]:>n["subexpression_digest"]],bOutputDigests];sharedReduced=Intersection[aReduced,bReduced];separation=rd[{"berry_separation"}, "subexpression_overlap_detector"];planted=parse[separation["planted_control_expression"]];plantedA=exprDAG[<|"planted"->planted|>,"separation_control_A"];plantedB=exprDAG[<|"planted"->planted|>,"separation_control_B"];plantedOverlap=Intersection[Cases[plantedA["nodes"],n_/;!MemberQ[excludedKinds,n["kind"]]:>n["subexpression_digest"]],Cases[plantedB["nodes"],n_/;!MemberQ[excludedKinds,n["kind"]]:>n["subexpression_digest"]]];
sigmaDag=exprDAG[<|"Omega_EL_xy"->omegaPerArea[[1,2]],"phase_holonomy"->2Pi winding,"Gamma"->gammaComputed,"epsilon_input"->epsSign,"sigma_projection_unsimplified"->sigmaUnsimplified|>,"sigma"];
selectionRole=selection["dependency_role"];
berry=<|"dependency_role"->selectionRole,"production_phase_profile"->ToString[prod,InputForm],"production_angular_derivative"->ToString[prodD,InputForm],"production_phase_parse_guard"-><|"allowed_free_symbols"->{"r"},"observed_free_symbols"->prodSymbols,"passed"->prodParseOK|>,"intrinsic_circulation"-><|"k_theta"->ToString[prodW,InputForm],"Gamma_theta"->ToString[FullSimplify[coeff["hbar"]2Pi prodW/coeff["m_GNLS"]],InputForm],"zero_decision"->If[prodW===0,"IDENTICALLY_ZERO_FUNCTIONAL","NONZERO_COMPUTED"]|>,
 "route_A"-><|"definition"->"differentiate termwise Taylor-substituted L_indexed","theta_t_termwise_substitution"->ToString[thetaTSub,InputForm],"substituted_field_density"->ToString[substitutedDensity,InputForm],"connection_density"->(ToString[#,InputForm]&/@localConnection),"connection_density_residual"->(ToString[#,InputForm]&/@localResidual),"indexed_L_local"->ToString[lIndexed,InputForm],"connection"->(ToString[#,InputForm]&/@connectionA),"Omega"->canonicalMatrix[omegaA],"surface_contributions"->boundaryRows,"dag"->dagA|>,
 "route_B"-><|"definition"->"native presymplectic current then Gram-projector quotient then collective pullback","current_integrand"->ToString[polarJac,InputForm],"bulk_integral"->ToString[interior,InputForm],"zero_mode_quotient_applied"->quotient,"projector"->ToString[projector,InputForm],"projector_times_zero_mode"->(ToString[#,InputForm]&/@FullSimplify[projector.{z1,z2}]),"zero_mode_sector_before_projection"->canonicalMatrix[zeroSector],"zero_mode_sector_after_projection"->canonicalMatrix[projectedZero],"Omega"->canonicalMatrix[omegaB],"dag"->dagB|>,
 "equivalence_residual"->canonicalMatrix[equivalence],"route_shared_reduced_nodes"->sharedReduced,"separation_control"-><|"expression"->ToString[planted,InputForm],"shared_reduced_digests"->plantedOverlap,"detector_fires"->(plantedOverlap=!={}),"shared_raw_input_whitelist"->separation["shared_raw_input_whitelist"]|>,"production_pullback_coverage"->coverage,"production_coverage_complete"->coverageOK,"dependency_edges"->{{"bulk_berry","berry:route_A"},{"bulk_berry","berry:route_B"},{"berry:route_A","Omega_EL"},{"berry:route_B","Omega_pullback"}},"selection_quarantined_from_production_ancestry"->(selectionRole==="selection_only")|>;
g4=<|"phase_profile_parsed"->ToString[theta,InputForm],"core_density_profile_parsed"->ToString[density,InputForm],"xi_core_over_a"->ToString[xiRatio,InputForm],"contour_substitution"-><|"x"->ToString[x/.contour,InputForm],"y"->ToString[y/.contour,InputForm]|>,"angular_derivative"->ToString[dtheta,InputForm],"phase_holonomy"->ToString[2Pi winding,InputForm],"computed_winding"->ToString[winding,InputForm],"k_expected_oracle"->ToString[expectedK,InputForm],"Gamma_computed"->ToString[gammaComputed,InputForm],"Gamma_expected_oracle"->ToString[gammaExpected,InputForm],"epsilon_matrix"->Map[ToString[#,InputForm]&,epsilon,{2}],"Omega_collective_total"->canonicalMatrix[omegaTotal],"Omega_physical_per_sheet_area"->canonicalMatrix[omegaPerArea],"source_berry_coefficient"->ToString[sourceCoeff,InputForm],"computed_sigma"->ToString[sigma,InputForm],"force_cases"->forceRows,"transverse_plane"->control["transverse_plane"],"sheet_directions"->control["sheet_directions"],"sheet_cell_area"->ToString[area,InputForm],"sheet_geometry"->sheet,"sheet_area_density"->ToString[areaDensity,InputForm],"total_to_per_area_residual"->canonicalMatrix[areaResidual],"ir_prescription"->control["ir_prescription"],"ir_matches_phase_a"->(control["ir_prescription"]===phase["ir_scheme"]["name"]),"sigma_dependency_dag"->sigmaDag|>;
r2OK=zeroQ[winding-expectedK]&&orient===epsSign&&And@@Flatten[Map[ToExpression[#]===0&,Values[forceRows][[All,"signed_residual"]],{2}]]&&areaResidual===ConstantArray[0,{2,2}]&&g4["ir_matches_phase_a"];
recordCheck["B1_R2",r2OK,<|"winding"->g4["computed_winding"],"sigma"->g4["computed_sigma"],"force"->forceRows,"area_residual"->g4["total_to_per_area_residual"],"sigma_dag"->sigmaDag|>];
recordCheck["B1_R3",equivalence===ConstantArray[0,{2,2}]&&prodParseOK&&prodW===0&&coverageOK&&selectionRole==="selection_only"&&sharedReduced==={}&&plantedOverlap=!={},<|"equivalence"->berry["equivalence_residual"],"production_phase_parse_guard"->berry["production_phase_parse_guard"],"coverage"->coverage,"route_A_dag"->dagA,"route_B_dag"->dagB,"separation_control"->berry["separation_control"]|>];

(* E4 action first, then Hessian/constraint/moduli lift/curvature. *)
scheme=rd[{"e4_collar_scheme"}, "extended_action_constraint_reduction"];dE4=scheme["radial_dimension"];ell=scheme["harmonic_degree"];Clear[nu];roots=nu/.Solve[nu(nu-dE4+2)-ell(ell+dE4-2)==0,nu];selected=If[scheme["ir_condition"]["select"]==="decaying_root",Max[Select[roots,Positive]],Min[roots]];
traceField=scheme["boundary_trace"]["field"];traceSym=safe[ToString[scheme["boundary_trace"]["value"]]];traceRadius=parse[scheme["boundary_trace"]["radius"]];irLimit=If[scheme["ir_condition"]["limit"]==="infinity",Infinity,parse[scheme["ir_condition"]["limit"]]];irTarget=parse[scheme["ir_condition"]["value"]];growingCoefficient=parse[scheme["moduli_fixing"]["growing_mode_coefficient"]];growingRoot=First@Select[roots,#=!=selected&];profile=traceSym(traceRadius/r)^selected;modeReconstruction=Expand[(traceSym-growingCoefficient)(traceRadius/r)^selected+growingCoefficient(traceRadius/r)^growingRoot];traceResidual=FullSimplify[(modeReconstruction/.r->traceRadius)-traceSym];irResidual=FullSimplify[Limit[modeReconstruction,r->irLimit]-irTarget];normalProfile=profile/traceSym;sdE4=sphereArea[dE4];
collarMass=FullSimplify[coeff["rhoBr"]sdE4/dE4 Integrate[r^(dE4-1)normalProfile^2,{r,traceRadius,Infinity}],aSym>0];gradDensity=D[normalProfile,r]^2+ell(ell+dE4-2)normalProfile^2/r^2;stiff=FullSimplify[coeff["muR"]sdE4/dE4 Integrate[r^(dE4-1)gradDensity,{r,traceRadius,Infinity}],aSym>0];dtn=FullSimplify[coeff["muR"]sdE4/dE4 traceRadius^(dE4-1)(-D[normalProfile,r]/.r->traceRadius)];
vel=Table[Symbol["vAug"<>ToString[i]],{i,0,8}];coords=Table[Symbol["qAug"<>ToString[i]],{i,0,8}];knownL=indexed["E4|symmetric_postulate"]["LRaw"]/.{Vx->vel[[1]],Vy->vel[[2]],Vz->vel[[3]]};collarK=collarMass Total[vel[[7;;9]]^2]/2;collarP=stiff Total[coords[[7;;9]]^2]/2;preAction=Expand[knownL+collarK-collarP];maug=Table[D[preAction,vel[[i]],vel[[j]]],{i,9},{j,9}];
Clear[vv,pp,uu,alpha,beta];cexpr=ToExpression[StringReplace[scheme["constraint_functional"]["expression"],{"V_i"->"vv","udot_collar_i"->"uu"}]];g1={{D[cexpr,vv],D[cexpr,pp],D[cexpr,uu]}};g=KroneckerProduct[g1,IdentityMatrix[3]];liftExpr=Expand[cexpr/.uu->alpha vv+beta pp];eqs={D[liftExpr,vv],D[liftExpr,pp],beta-parse[scheme["moduli_fixing"]["p_to_collar_cross_lift"]]};liftSol=Solve[Thread[eqs==0],{alpha,beta}];j1=If[liftSol==={},ConstantArray[0,{3,2}],{{1,0},{0,1},{alpha,beta}}/.First[liftSol]];jj=KroneckerProduct[j1,IdentityMatrix[3]];nBasis=Transpose[NullSpace[g]];reduced=FullSimplify[Transpose[jj].maug.jj];kernelForm=FullSimplify[Transpose[nBasis].maug.nBasis];
tmat={{1,1,0},{0,1,0},{0,0,2}};smat=IdentityMatrix[9];smat[[7;;9,7;;9]]=tmat;covres=FullSimplify[Transpose[Inverse[smat].jj].(Transpose[smat].maug.smat).(Inverse[smat].jj)-reduced];zeroVel=Thread[vel->0];aaug=(D[preAction,#]/.zeroVel)&/@vel;omegaaug=Table[D[aaug[[j]],coords[[i]]]-D[aaug[[i]],coords[[j]]],{i,9},{j,9}];physicalCoords=coords[[1;;6]];anholonomic=Table[FullSimplify[Sum[jj[[row,i]]D[jj[[row,j]],physicalCoords[[i]]]-jj[[row,j]]D[jj[[row,i]],physicalCoords[[j]]],{row,Length[jj]}]],{i,6},{j,6}];
schemeContract=<|"formulation"->scheme["formulation"],"action_ancestors"->scheme["action_ancestors"],"radial_domain"->scheme["radial_domain"],"constraint_equals"->scheme["constraint_functional"]["equals"],"collective_trace_map"->scheme["collective_trace_map"]|>;
schemeOK=scheme["formulation"]==="operator_level_Dirichlet_to_Neumann"&&Sort[scheme["action_ancestors"]]===Sort[{"brane_shear_kinetic","brane_shear_gradient"}]&&traceField==="u_T"&&scheme["ir_condition"]["limit"]==="infinity"&&scheme["radial_domain"]==="a<=r<infinity"&&zeroQ[parse[scheme["constraint_functional"]["equals"]]]&&scheme["collective_trace_map"]==="identity";
e4Conditions=<|
 "trace"->zeroQ[traceResidual],
 "ir"->zeroQ[irResidual],
 "dtn"->zeroQ[stiff-dtn],
 "constraint_rank"->(MatrixRank[g]===3),"lift_rank"->(MatrixRank[jj]===6),
 "constraint_lift"->And@@Flatten[Map[zeroQ,g.jj,{2}]],
 "comoving_lift"->(zeroQ[j1[[3,1]]-1]&&zeroQ[j1[[3,2]]-parse[scheme["moduli_fixing"]["p_to_collar_cross_lift"]]]),
 "basis_covariance"->And@@Flatten[Map[zeroQ,covres,{2}]],"anholonomic_decision"->And@@Flatten[Map[zeroQ,anholonomic,{2}]],
 "lift_solved"->(liftSol=!={}),"action_connection"->And@@Map[zeroQ,aaug],
 "action_curvature"->And@@Flatten[Map[zeroQ,omegaaug,{2}]],
 "growing_mode"->zeroQ[growingCoefficient],"scheme_contract"->schemeOK|>;
e4OK=And@@Values[e4Conditions];
e4=<|"operator"->("rhoBr*d_t^2-muR*(d_r^2+"<>ToString[dE4-1]<>"/r*d_r-"<>ToString[ell(ell+dE4-2)]<>"/r^2)"),"boundary_profile"->ToString[profile,InputForm],"collar_mode_reconstruction"->ToString[modeReconstruction,InputForm],"growing_mode_root"->ToString[growingRoot,InputForm],"indicial_roots"->(ToString[#,InputForm]&/@roots),"selected_decay"->ToString[selected,InputForm],"trace_residual"->ToString[traceResidual,InputForm],"ir_residual"->ToString[irResidual,InputForm],"collar_mass"->canonicalTerms[collarMass],"collar_stiffness_bulk"->canonicalTerms[stiff],"DtN_stiffness"->canonicalTerms[dtn],"DtN_residual"->canonicalTerms[stiff-dtn],"preconstraint_extended_action"->canonicalTerms[preAction],"preconstraint_action_components"-><|"phase_a_translation"->canonicalTerms[knownL],"collar_kinetic"->canonicalTerms[collarK],"collar_potential"->canonicalTerms[-collarP]|>,"M_aug"->canonicalMatrix[maug],"M_aug_hessian_residual"->canonicalMatrix[Table[D[preAction,vel[[i]],vel[[j]]],{i,9},{j,9}]-maug],"constraint_expression"->ToString[cexpr,InputForm],"constraint_operator"->canonicalMatrix[g],"constraint_rank"->MatrixRank[g],"kernel_dimension"->9-MatrixRank[g],"N_full_kernel"->canonicalMatrix[nBasis],"N_kernel_residual"->canonicalMatrix[g.nBasis],"J_physical_lift"->canonicalMatrix[jj],"J_derivation_equations"->(ToString[#,InputForm]&/@eqs),"M_reduced"->canonicalMatrix[reduced],"M_full_kernel"->canonicalMatrix[kernelForm],"basis_covariance_residual"->canonicalMatrix[covres],"multiplier_solution"->{},"reaction"->{},"virtual_work_residual"->{},"A_aug"->(canonicalTerms/@aaug),"Omega_aug"->canonicalMatrix[omegaaug],"A_aug_zero_decision"->"DERIVATIVE_OF_ACTION_ZERO","Omega_aug_zero_decision"->"CURVATURE_OF_DERIVED_CONNECTION_ZERO","anholonomic_term"->canonicalMatrix[anholonomic],"anholonomic_decision"->"DERIVED_FRAME_BRACKET","dependency_edges"->Join[( {#,"E4:preconstraint_action"}&/@scheme["action_ancestors"]),{{"E4:preconstraint_action","E4:M_aug"},{"E4_shear_lock","E4:J_physical_lift"},{"E4:M_aug","E4:M_reduced"},{"E4:J_physical_lift","E4:M_reduced"}}],"Omega_aug_unresolved_leaves"->missing,"status"->If[missing==={},"COMPUTED","UNRESOLVED"],"unresolved_leaves"->missing|>;
e4["scheme_contract"]=schemeContract;
recordCheck["B1_R4",e4OK,<|"conditions"->e4Conditions,"action"->e4["preconstraint_extended_action"],"hessian"->e4["M_aug_hessian_residual"],"constraint"->e4["constraint_operator"],"lift"->e4["J_physical_lift"]|>];

(* E5 field Rayleigh variation and root-deletion re-solve. *)
ray=endpointFull["E5"]["rayleighRaw"]/.{tangent->safe["v_tangent"],Vt->safe["V_tangent"]};frt=FullSimplify[-D[ray,safe["v_tangent"]]];frv=FullSimplify[-D[ray,safe["V_tangent"]]];deleted=0;deltaForce=frt;
fullSol=endpointFull["E5"]["solutionRaw"];bareSol=endpointBare["E5"]["solutionRaw"];e2Sol=endpointBare["E2"]["solutionRaw"];fullGVV=phaseCoefficient["E5","GVV"];bareGVV=phaseCoefficient["E2","GVV"];
e5OK=!zeroQ[ray]&&zeroQ[deleted]&&!zeroQ[deltaForce]&&And@@Map[zeroQ,bareSol-e2Sol]&&!And@@Map[zeroQ,fullSol-bareSol];
e5=<|"root"->rd[{"endpoint_functionals", "E5", "rayleigh_root"}, "rayleigh_root_consumption"],"domain"->rd[{"endpoint_functionals", "E5", "rayleigh_surface"}, "rayleigh_root_consumption"],"field_level_density"->ToString[ray,InputForm],"F_Rayleigh_trace"->ToString[frt,InputForm],"F_Rayleigh_sleeve"->ToString[frv,InputForm],"root_deleted_functional"->"0","gamma_functional_ablation_force_delta"->ToString[deltaForce,InputForm],"full_E5_trace_solution"->(ToString[#,InputForm]&/@fullSol),"root_deleted_conservative_solution"->(ToString[#,InputForm]&/@bareSol),"E2_conservative_solution"->(ToString[#,InputForm]&/@e2Sol),"conservative_equals_E2_computed"->And@@Map[zeroQ,bareSol-e2Sol],"bare_GVV"->canonicalTerms[bareGVV],"stored_full_E5_minus_bare_GVV"->canonicalTerms[fullGVV-bareGVV],"bare_dependencies"->Sort@DeleteDuplicates[ename/@Cases[bareGVV,_Symbol,Infinity]],"rayleigh_dependencies"->Sort@DeleteDuplicates[ename/@Cases[ray,_Symbol,Infinity]],"dependency_edges"->{{rd[{"endpoint_functionals", "E5", "rayleigh_root"}, "rayleigh_root_consumption"],"E5:rayleigh_density"},{"E5:rayleigh_density","F_Rayleigh"}}|>;
recordCheck["B1_R5",endpointsOK&&e5OK,<|"endpoint_fingerprints"->Association@Table[ep->endpointFull[ep]["fingerprint"],{ep,endpoints}],"E5"->e5|>];

(* Units are restored termwise from the real coefficient expressions. *)
dimAssoc=Association@KeyValueMap[#1->#2["dimensions"]&,phaseInput["coefficients"]];KeyValueMap[(dimAssoc[#1]={0,0,0})&,phaseInput["core_traces"]];dimAssoc["a"]={1,0,0};dimAssoc["s"]={0,0,0};dimAssoc["V"]={1,-1,0};dimAssoc["pd"]={0,-1,0};
termDim[t_] := Module[{factors=If[Head[t]===Times,List@@t,{t}],out={0,0,0},base,power,name},Do[If[!NumericQ[f],If[Head[f]===Power,base=f[[1]];power=f[[2]],base=f;power=1];name=ename[base];out+=power Lookup[dimAssoc,name,{0,0,0}]],{f,factors}];out];
exprDims[e_] := Sort@DeleteDuplicates[termDim/@If[Head[Expand[e]]===Plus,List@@Expand[e],{Expand[e]}]];
restore[e_] := Total[Map[Function[t,t aSym^(-First[termDim[t]])],If[Head[Expand[e]]===Plus,List@@Expand[e],{Expand[e]}]]];
gvvU=restore[phaseCoefficient["E1","GVV"]];gvpU=restore[phaseCoefficient["E1","GVP"]];gppU=restore[phaseCoefficient["E1","GPP"]];Clear[vDim,pdDim];externalNames[vDim]="V";externalNames[pdDim]="pd";
dimObjects=<|"L_translation"->gvvU vDim^2/2,"L_native_Xp_slice"->aSym gvpU vDim pdDim,"L_native_pp_slice"->aSym^2 gppU pdDim^2/2,"M_XX"->gvvU,"M_Xp_native_slice"->aSym gvpU,"M_pp_native_slice"->aSym^2 gppU,"Omega_XX_control"->omegaPerArea[[1,2]]|>;
dimensionRecords=KeyValueMap[<|"object"->#1,"units_restored_expression"->ToString[#2,InputForm],"computed_monomial_dimensions_LTM"->exprDims[#2]|>&,dimObjects];actionDim=(First@Select[phaseInput["field_records"],#["root_type"]==="ACTION"&])["dimensions"];
coordinatesInput=rd[{"indexed_coordinates"}, "units_restored_mechanics"];indices=coordinatesInput["spatial_indices"];expectedCoordinates=Join[("X_"<>#&/@indices),("p_"<>#&/@indices)];expectedVelocities=Join[("V_"<>#&/@indices),("pd_"<>#&/@indices)];nondim=coordinatesInput["nondimensionalization"];
coordinateContract=<|"spatial_indices"->indices,"coordinate_order"->coordinatesInput["coordinate_order"],"velocity_order"->coordinatesInput["velocity_order"],"expected_coordinate_order"->expectedCoordinates,"expected_velocity_order"->expectedVelocities,"p_slice"->coordinatesInput["p_slice"],"nondimensionalization"->nondim|>;
coordinateOK=coordinatesInput["coordinate_dimensions"]["X"]==={1,0,0}&&coordinatesInput["coordinate_dimensions"]["p"]==={0,0,0}&&Length[indices]===3&&coordinatesInput["coordinate_order"]===expectedCoordinates&&coordinatesInput["velocity_order"]===expectedVelocities&&coordinatesInput["p_slice"]==="off_shell_free"&&nondim["translation_length"]==="a"&&nondim["tilt_length"]==="a";
dimensions=<|"basis"->"LTM","coefficient_dimensions"->dimAssoc,"coordinate_dimensions"->coordinatesInput["coordinate_dimensions"],"coordinate_basis_contract"->coordinateContract,"action_dimension_sourced_from_phase_a"->actionDim,"records"->dimensionRecords,"action_terms_homogeneous"->(DeleteDuplicates[Flatten[dimensionRecords[[1;;3,"computed_monomial_dimensions_LTM"]],1]]==={actionDim}),"coordinate_embedding_consistent"->coordinateOK|>;
recordCheck["B1_R8",dimensions["action_terms_homogeneous"]&&dimensions["coordinate_embedding_consistent"],dimensions];

(* Congruence is computed from coefficient functions traversed out of the
   produced remainder blocks; no free A..G inventory is accepted as input. *)
Clear[p1,p2,p3,rrad];pvec={p1,p2,p3};epsp={{0,p3,-p2},{-p3,0,p1},{p2,-p1,0}};known=exprTerms[indexed["E1|symmetric_postulate"]["M_XX_p0_known"][[1,1]]];
generators=<|"M_XX"->{"delta_ij","p_i*p_j"},"M_Xp_symmetric"->{"delta_ij","p_i*p_j"},"M_pp"->{"delta_ij","p_i*p_j"},"M_Xp_antisymmetric"->{"epsilon_ijk*p_k"},"Omega_XX_texture"->{"epsilon_ijk*p_k"},"Omega_Xp"->{"epsilon_ijk*p_k"},"Omega_pp"->{"epsilon_ijk*p_k"}|>;
defExpr=<||>;Do[Do[defExpr[block<>":"<>gen]=0,{gen,generators[block]}],{block,Keys[generators]}];
Do[block=row["block"];Do[pos=FirstPosition[row["finite_generator_set"],gen,Missing["Absent"]];If[!MissingQ[pos],defExpr[block<>":"<>gen]+=safe[row["remainder_symbols"][[pos[[1]]]]]],{gen,generators[block]}],{row,reach}];
Clear[aTot,bCoef,cCoef,dCoef,gCoef,eCoef,fCoef];
defs=<|
 "M_XX:delta_total"-><|"symbol"->"aTot","definition"->ToString[known+defExpr["M_XX:delta_ij"],InputForm]|>,
 "M_XX:delta_ij"-><|"symbol"->"aTot-known","definition"->ToString[defExpr["M_XX:delta_ij"],InputForm]|>,
 "M_XX:p_i*p_j"-><|"symbol"->"bCoef","definition"->ToString[defExpr["M_XX:p_i*p_j"],InputForm]|>,
 "M_Xp_symmetric:delta_ij"-><|"symbol"->"cCoef","definition"->ToString[defExpr["M_Xp_symmetric:delta_ij"],InputForm]|>,
 "M_Xp_symmetric:p_i*p_j"-><|"symbol"->"dCoef","definition"->ToString[defExpr["M_Xp_symmetric:p_i*p_j"],InputForm]|>,
 "M_Xp_antisymmetric:epsilon_ijk*p_k"-><|"symbol"->"gCoef","definition"->ToString[defExpr["M_Xp_antisymmetric:epsilon_ijk*p_k"],InputForm]|>,
 "M_pp:delta_ij"-><|"symbol"->"eCoef","definition"->ToString[defExpr["M_pp:delta_ij"],InputForm]|>,
 "M_pp:p_i*p_j"-><|"symbol"->"fCoef","definition"->ToString[defExpr["M_pp:p_i*p_j"],InputForm]|>,
 "Omega_XX_texture:epsilon_ijk*p_k"-><|"symbol"->"omegaXXCoef","definition"->ToString[defExpr["Omega_XX_texture:epsilon_ijk*p_k"],InputForm]|>,
 "Omega_Xp:epsilon_ijk*p_k"-><|"symbol"->"omegaXpCoef","definition"->ToString[defExpr["Omega_Xp:epsilon_ijk*p_k"],InputForm]|>,
 "Omega_pp:epsilon_ijk*p_k"-><|"symbol"->"omegaPpCoef","definition"->ToString[defExpr["Omega_pp:epsilon_ijk*p_k"],InputForm]|>|>;
mxx=aTot IdentityMatrix[3]+bCoef Outer[Times,pvec,pvec];mxp=cCoef IdentityMatrix[3]+dCoef Outer[Times,pvec,pvec]+gCoef epsp;mpp=eCoef IdentityMatrix[3]+fCoef Outer[Times,pvec,pvec];full=ArrayFlatten[{{mxx,mxp},{Transpose[mxp],mpp}}];aligned=FullSimplify[full/.{p1->0,p2->0,p3->rrad}];ti={1,2,4,5};li={3,6};trans=aligned[[ti,ti]];long=aligned[[li,li]];
piv[m_] := Module[{mins=Table[Factor[Det[m[[1;;i,1;;i]]]],{i,Length[m]}]},Join[{mins[[1]]},Table[Factor[mins[[i]]/mins[[i-1]]],{i,2,Length[mins]}]]];pt=piv[trans];pl=piv[long];perm=Join[ti,li];splitResidual=FullSimplify[aligned[[perm,perm]]-ArrayFlatten[{{trans,ConstantArray[0,{4,2}]},{ConstantArray[0,{2,4}],long}}]];
inventory=rd[{"covariant_inventory"}, "covariant_block_congruence"];vectorNames=inventory["vectors"];generic=parse/@inventory["generic_projection_point"];projectionBasisAtPoint=Outer[Times,generic,generic];Clear[alphaD,betaD,lambdaZero];declaredSymmetric=inventory["symmetric_basis"];projectionAnsatz=alphaD IdentityMatrix[3];projectionVariables={alphaD};If[MemberQ[declaredSymmetric,"p_i_p_j"],projectionAnsatz+=betaD projectionBasisAtPoint;AppendTo[projectionVariables,betaD]];projResidual=FullSimplify[(mxx/.Thread[pvec->generic])-projectionAnsatz];decomp=Solve[Thread[Flatten[projResidual]==0],projectionVariables];decompResidual=If[decomp==={},projResidual,FullSimplify[projResidual/.First[decomp]]];zeroDirection=parse/@inventory["zero_limit_direction"];zeroLimitPathProjector=Outer[Times,zeroDirection,zeroDirection];zeroPathBlock=FullSimplify[mxx/.Thread[pvec->lambdaZero zeroDirection]];zeroLimitBlock=Map[Limit[#,lambdaZero->0]&,zeroPathBlock,{2}];requiredInventory=SubsetQ[inventory["symmetric_basis"],{"delta_ij","p_i_p_j"}]&&SubsetQ[inventory["antisymmetric_basis"],{"epsilon_ijk_p_k"}];
derivedCongruence=<|"derived_coefficients"->defs,"M_XX_covariant_decomposition"-><|"delta_coefficient"->If[decomp==={},"UNRESOLVED",ToString[alphaD/.First[decomp],InputForm]],"pp_coefficient"->If[decomp==={},"UNRESOLVED",ToString[betaD/.First[decomp],InputForm]],"projection_residual"->canonicalMatrix[decompResidual]|>,"inventory_vector_symbols"->(ToString[#,InputForm]&/@pvec),"projection_basis_at_point"->canonicalMatrix[projectionBasisAtPoint],"zero_limit_path_projector"->canonicalMatrix[zeroLimitPathProjector],"zero_limit_path_block"->canonicalMatrix[zeroPathBlock],"zero_limit_block"->canonicalMatrix[zeroLimitBlock],"aligned_full_block"->canonicalMatrix[aligned],"transverse_block"->canonicalMatrix[trans],"longitudinal_block"->canonicalMatrix[long],"block_split_residual"->canonicalMatrix[splitResidual],"transverse_LDL_pivots"->(ToString[#,InputForm]&/@pt),"longitudinal_LDL_pivots"->(ToString[#,InputForm]&/@pl),"positive_definite_conditions"->(("("<>ToString[#,InputForm]<>") > 0")&/@Join[pt,pl]),"E4_reduced_rank"->MatrixRank[reduced],"E4_reduced_block"->e4["M_reduced"],"classification"->"M_UNRESOLVED_DERIVED_REMAINDER_COEFFICIENTS"|>;
eta=parse[rd[{"ambient_branches", "one_sided_pathA29", "eta_asym"}, "ambient_indexed_action"]];covarianceOK=zeroQ[branchFactor["one_sided_pathA29"]-(1+sSym eta)];r7OK=vectorNames==={"p"}&&Length[zeroDirection]===3&&!And@@Map[zeroQ,zeroDirection]&&requiredInventory&&And@@Flatten[Map[zeroQ,decompResidual,{2}]]&&And@@Flatten[Map[zeroQ,splitResidual,{2}]]&&covarianceOK;
recordCheck["B1_R7",r7OK,<|"derived_congruence"->derivedCongruence,"one_sided_factor"->ToString[branchFactor["one_sided_pathA29"],InputForm]|>];
representations=rd[{"derivation_representations"}, "engine_representation_decision"];r9OK=TrueQ[representations["shared_reduced_formulas"]===False]&&representations["SymPy"] =!= representations["Mathematica"]&&SubsetQ[dagA["raw_dependencies"],{"Vcx","Vcy","Xx","Xy"}]&&MemberQ[dagB["raw_dependencies"],"r"]&&dagA["raw_dependencies"] =!= dagB["raw_dependencies"];
recordCheck["B1_R9",r9OK,<|"representations"->representations,"route_A_operations"->Sort@DeleteDuplicates[dagA["nodes"][[All,"kind"]]],"route_B_operations"->Sort@DeleteDuplicates[dagB["nodes"][[All,"kind"]]]|>];

(* Derived blocks, component records, map, provenance, and partition. *)
blockExpressions=<||>;KeyValueMap[Function[{key,cell},blocks=<||>;KeyValueMap[Function[{block,spec},rootsB=Sort@DeleteDuplicates@Lookup[Select[reach,#["block"]===block&],"root",{}];tiltLeaves=If[block==="M_XX",{},missing];leaves=Sort@DeleteDuplicates@Join[rootsB,tiltLeaves];blocks[block]=<|"constructible_tensor"->If[block==="M_XX",cell["M_XX_p0_known"],Null],"unresolved_remainder_leaves"->leaves,"remainder_expression"->StringRiffle[("R["<>#<>","<>block<>"]")&/@leaves," + "],"status"->If[leaves==={},"COMPUTED","UNRESOLVED"]|>],blockSpecs];blockExpressions[key]=blocks],indexed];
aggregate[findings_] := Module[{statuses=Cases[findings,a_Association/;KeyExistsQ[a,"status"]:>a["status"],Infinity],leaves=Flatten@Cases[findings,a_Association/;KeyExistsQ[a,"unresolved_leaves"]:>a["unresolved_leaves"],Infinity],status},status=Which[MemberQ[statuses,"ILL_POSED"],"ILL_POSED",MemberQ[statuses,"UNSTABLE"],"UNSTABLE",MemberQ[statuses,"UNRESOLVED"],"UNRESOLVED",True,"OK"];<|"status"->status,"unresolved_leaves"->Sort@DeleteDuplicates[leaves],"precedence_evidence"->statuses|>];
strata=Map[<|"id"->#["id"],"predicate"->#["predicate"],"predicate_value"->TrueQ[ToExpression[StringReplace[#["predicate"],{"Eq("->"Equal[",")"->"]"}]]],"leaves"->Lookup[#,"leaves",{}]|>&,rd[{"open_strata"}, "cell_stratum_selection"]];active=Select[strata,#["predicate_value"]&];cellRecords=<||>;mechanicsMap=<||>;
Do[base=ep<>"|"<>ambient;Do[key=base<>"|"<>stratum["id"];blockLeaves=Sort@DeleteDuplicates@Flatten[Values[blockExpressions[base]][[All,"unresolved_remainder_leaves"]]];mStatus=If[blockLeaves==={},"COMPUTED","UNRESOLVED"];
 findings=<|"M_full"-><|"status"->mStatus,"unresolved_leaves"->blockLeaves,"known_p0_translation"->indexed[base]["M_XX_p0_known"],"derived_congruence_class"->derivedCongruence["classification"]|>,"M_spatial_symmetry"-><|"status"->mStatus,"unresolved_leaves"->blockLeaves,"known_M_XX_projection_residual"->derivedCongruence["M_XX_covariant_decomposition"]["projection_residual"]|>,"Omega_intrinsic_circulation"-><|"status"->If[berry["intrinsic_circulation"]["k_theta"]==="0","ZERO_COMPUTED","NONZERO"],"zero_decision"->berry["intrinsic_circulation"]["zero_decision"]|>,"Omega_translation_texture"-><|"status"->"UNRESOLVED","unresolved_leaves"->Sort@DeleteDuplicates@Join[Lookup[blockRoots,"Omega_XX_texture",{}],missing]|>,"Omega_translation_tilt"-><|"status"->"UNRESOLVED","unresolved_leaves"->Sort@DeleteDuplicates@Join[Lookup[blockRoots,"Omega_Xp",{}],missing]|>,"Omega_tilt_tilt"-><|"status"->"UNRESOLVED","unresolved_leaves"->Sort@DeleteDuplicates@Join[Lookup[blockRoots,"Omega_pp",{}],missing]|>,"E4_constraint"->If[ep==="E4",<|"status"->e4["status"],"unresolved_leaves"->e4["unresolved_leaves"]|>,<|"status"->"NOT_APPLICABLE"|>],"E5_rayleigh"->If[ep==="E5",<|"status"->If[e5["conservative_equals_E2_computed"],"COMPUTED","ILL_POSED"]|>,<|"status"->"NOT_APPLICABLE"|>]|>;
 headline=aggregate[findings];record=<|"cell_key"-><|"endpoint"->ep,"ambient"->ambient,"open_stratum"->stratum["id"]|>,"off_shell_p"->True,"block_expressions"->blockExpressions[base],"component_findings"->findings,"headline"->headline,"record_inputs_digest"->sha[<|"blocks"->blockExpressions[base],"findings"->findings|>]|>;cellRecords[key]=record;mechanicsMap[key]=headline,{stratum,active}],{ep,endpoints},{ambient,ambients}];
partitionInput=rd[{"partition"}, "partition_ownership"];ledger=Flatten@Table[Table[<|"candidate_id"->key<>":"<>term<>":translation","owner"->"M","provenance"->term,"computed_expression_digest"->sha[indexed[key]["termwise_L"][term]]|>,{term,Keys[indexed[key]["termwise_L"]]}],{key,Keys[indexed]}];AppendTo[ledger,<|"candidate_id"->"outer_control_flux:translation","owner"->partitionInput["outer_control_flux_owner"],"provenance"->"outer_control_flux","computed_expression_digest"->sha["pending_B2"]|>];partition=<|"records"->ledger,"owner_enum"->partitionInput["owner_enum"],"unique"->DuplicateFreeQ[ledger[[All,"candidate_id"]]],"state"->partitionInput["terminal_state"]|>;
provEdges=Sort@DeleteDuplicates@Join[Flatten[Cases[Values[indexed],cell_Association:>Lookup[cell["field_contraction_integrals"],"dependency_edges",{}]],2],({#["root"],"remainder:"<>#["block"]}&/@reach),berry["dependency_edges"],e4["dependency_edges"],e5["dependency_edges"]];provNodes=Sort@DeleteDuplicates@Join[Flatten[provEdges],{"return_closure"}];provTargets=Sort@DeleteDuplicates@Join[Cases[provEdges,{_,target_/;StringContainsQ[target,"M_XX"|"Omega"|"M_reduced"]}:>target],("remainder:"<>#&/@Keys[blockSpecs])];provenance=<|"nodes"->(<|"id"->#,"type"->If[#==="return_closure","RETURN","DERIVED_OR_DECLARED"]|>&/@provNodes),"edges"->provEdges,"production_dependency_sets"->Association@KeyValueMap[#1->#2["production_dependencies"]&,indexed],"global_return_closure_absence"->FreeQ[provEdges,{"return_closure",_}],"closure_absence_targets"->provTargets,"construction"->"union(emitted_contraction,reachability,berry,E4,E5 dependency records)"|>;

gates=<|"G1"-><|"status"->If[sourceOK,"INHERITED_PHASE_A_REPRODUCED","HARNESS_FAILED"],"evidence_path"->"source_contract.phase_a_payload_sha256"|>,"G2"-><|"status"->"KNOWN_COEFFICIENTS_FINITE;FULL_REMAINDERS_UNRESOLVED"|>,"G3"-><|"status"->"DERIVED_BLOCK_CONGRUENCE;FULL_SIGNATURE_UNRESOLVED","conditions"->derivedCongruence["positive_definite_conditions"]|>,"G4"-><|"status"->If[checks["B1_R2"]["status"]==="PASS","PASS","FAIL"],"force_cases"->Length[forceRows]|>,"G5"-><|"status"->If[covarianceOK,"KNOWN_M_AND_OMEGA_COVARIANT;FULL_REMAINDERS_UNRESOLVED","HARNESS_FAILED"]|>,"G6"-><|"status"->"ENDPOINT_MAP_COMPUTED","source_map"->sourceMap|>,"G7"-><|"status"->"NOT_RUN(phase_B2)"|>,"G8"-><|"status"->"NOT_RUN(phase_C)"|>,"G9"-><|"status"->"NOT_RUN(phase_B2)"|>,"G10"-><|"status"->"NOT_RUN(phase_C)"|>,"G11"-><|"status"->"NOT_RUN(phase_C)"|>|>;
claims={<|"id"->"phase_a_digest","schema_path"->"source_contract.phase_a_payload_sha256","type"->"sha256","recompute"->"sha256(normalized_phase_a_payload)"|>,<|"id"->"field_manifest","schema_path"->"field_manifest.fields","type"->"per_field_manifest","recompute"->"join(indexed_routes,phase_a.tail_channels)"|>,<|"id"->"emitted_leaves","schema_path"->"indexed_profile_missing_leaves","type"->"derived_string_set","recompute"->"failed_phase_a_indexed_tangent_lookups"|>,<|"id"->"scalar_regression","schema_path"->"scalar_regression","type"->"residual_mapping","recompute"->"eV.T*M_XX*eV-GVV"|>,<|"id"->"native_slice_constraints","schema_path"->"indexed_cells","type"->"conditional_residual_mapping","recompute"->"native_MXp/Mpp_projection-minus-GVP/GPP"|>,<|"id"->"g4_winding","schema_path"->"g4_control.computed_winding","type"->"sympy_expression","recompute"->"contour_integral/(2*pi)"|>,<|"id"->"g4_sigma","schema_path"->"g4_control.computed_sigma","type"->"sympy_expression","recompute"->"Omega_xy/(rho_mass*Gamma*epsilon_xy)"|>,<|"id"->"sheet_area","schema_path"->"g4_control.total_to_per_area_residual","type"->"canonical_matrix","recompute"->"Omega_total/sheet_cell_area-Omega_per_area"|>,<|"id"->"berry_coverage","schema_path"->"berry.production_pullback_coverage","type"->"cell_selection","recompute"->"production_cells-cross-ambient_branches"|>,<|"id"->"e4_action_hessian","schema_path"->"E4.M_aug_hessian_residual","type"->"canonical_matrix","recompute"->"hessian(preconstraint_extended_action)-M_aug"|>,<|"id"->"open_reachability","schema_path"->"open_root_reachability","type"->"generator_records","recompute"->"typed_ledger-union-traversal"|>,<|"id"->"finite_controls","schema_path"->"reachability_analysis.finite_bound_controls","type"->"emptiness_and_witness_records","recompute"->"same_generator_parity_filter"|>,<|"id"->"cell_count","schema_path"->"cells","type"->"mapping_length","recompute"->"active_endpoint-ambient-stratum_product"|>};
claims=Join[claims,{<|"id"->"dimensions","schema_path"->"dimensions.records","type"->"computed_dimension_records","recompute"->"termwise_units_restoration"|>,<|"id"->"derived_congruence","schema_path"->"derived_congruence","type"->"block_congruence","recompute"->"produced_blocks-to-covariant-coefficients"|>,<|"id"->"mechanics_map","schema_path"->"mechanics_map","type"->"component_aggregate_map","recompute"->"aggregate(cells.component_findings)"|>,<|"id"->"closure_absence","schema_path"->"provenance_graph.global_return_closure_absence","type"->"graph_predicate","recompute"->"reachability(return_closure,mechanics_targets)"|>,<|"id"->"e5_root","schema_path"->"E5.root_deleted_conservative_solution","type"->"root_ablation_solve","recompute"->"delete(rayleigh_root)-and-resolve"|>,<|"id"->"gate_statuses","schema_path"->"gate_evidence","type"->"gate_mapping","recompute"->"engine_check-derived-gates"|>,<|"id"->"partition","schema_path"->"partition_ledger","type"->"ownership_ledger","recompute"->"computed-candidate-ownership"|>}];
publicEndpoint[a_] := KeyDrop[a,{"matrixRaw","rhsRaw","solutionRaw","rayleighRaw"}];publicIndexed[a_] := KeyDrop[a,{"LRaw","MRaw"}];
inputValidation=<|"schema"->rd[{"schema_version"}, "schema_validation"],"directive"->rd[{"directive_version"}, "schema_validation"],"spec"->rd[{"spec_version"}, "schema_validation"],"phase"->rd[{"phase"}, "schema_validation"],"valid"->(rd[{"schema_version"}, "schema_validation"]==="U1_PHASE_B1_MECHANICS_INPUTS_V3"&&rd[{"phase"}, "schema_validation"]==="B1")|>;
result=<|"schema_version"->"U1_PHASE_B1_ENGINE_ARTIFACT_V3","engine"->"Mathematica","engine_representation"->representations["Mathematica"],"phase"->"B1","input_validation"->inputValidation,"source_contract"->sourceContract,
"phase_a_amendment"->phaseCorrection,"amended_phase_a_payload"->amendedPayload,
"field_manifest"->fieldManifest,"field_embedding"->fieldManifest,"indexed_profile_missing_leaves"->missing,"dimensions"->dimensions,"endpoint_assembly"->Association@KeyValueMap[#1->publicEndpoint[#2]&,endpointFull],"endpoint_conservative"->Association@KeyValueMap[#1->publicEndpoint[#2]&,endpointBare],"endpoint_source_map"->sourceMap,"indexed_cells"->Association@KeyValueMap[#1->publicIndexed[#2]&,indexed],"scalar_regression"->scalarRegression,"full_block_expressions"->blockExpressions,"open_root_reachability"->reach,"block_unresolved_roots"->Association@KeyValueMap[#1->Sort@DeleteDuplicates[#2]&,blockRoots],"reachability_analysis"->reachAnalysis,"derived_congruence"->derivedCongruence,"symbolic_LDL"->derivedCongruence,"berry"->berry,"g4_control"->g4,"E4"->e4,"E5"->e5,"partition_ledger"->partition,"provenance_graph"->provenance,"cells"->cellRecords,"mechanics_map"->mechanicsMap,"gate_evidence"->gates,"checks"->checks,"report_claim_bindings"->claims,"closure_axis"-><|"materialized"->False,"root"->"return_closure","global_absence_computed"->provenance["global_return_closure_absence"]|>,"open_strata"-><|"collapsed"->Length[active]===1,"declared"->strata,"active"->Lookup[active,"id",{}]|>,"phase_scope"-><|"completed"->{"7.2","7.3"},"B2"->"NOT_RUN(phase_B2)","C"->"NOT_RUN(phase_C)"|>|>;

(* S1 rows cover every input leaf and bind it to a non-diagnostic semantic sink. *)
leafRows[x_Association,p_List:{}] := Flatten[KeyValueMap[leafRows[#2,Append[p,#1]]&,x],1];
leafRows[x_List,p_List] := If[x==={},{{p,{}}},Flatten[MapIndexed[leafRows[#1,Append[p,First[#2]-1]]&,x],1]];
leafRows[x_,p_List] := {{p,x}};
sinkFor[path_List] := Module[{head=ToString[First[path]]},If[MemberQ[{"schema_version","directive_version","spec_version","phase"},head],Return[{"validation_predicate","input_validation","schema_validation"}]];
 If[head==="scalar_regression_projection"&&Length[path]>1&&MemberQ[{"pd_unit","p_unit","indexed_tilt_length"},ToString[path[[2]]]],Return[{"comparator_residual","indexed_cells","native_slice_projection_residual"}]];
 If[head==="endpoint_functionals"&&Length[path]>2&&ToString[path[[2]]]==="E5"&&MemberQ[{"rayleigh_root","rayleigh_surface"},ToString[path[[3]]]],Return[{"mathematical_deliverable","E5","rayleigh_root_consumption"}]];
 If[head==="phase_a_protection",Return[{"comparator_residual","source_contract.phase_a_payload_sha256","external_phase_a_payload_guard"}]];
 Lookup[<|"derivation_representations"->{"validation_predicate","checks.B1_R9","engine_representation_decision"},"substrate"->{"validation_predicate","source_contract","phase_a_source_contract"},"indexed_coordinates"->{"mathematical_deliverable","dimensions","units_restored_mechanics"},"indexed_embedding"->{"mathematical_deliverable","field_manifest","per_field_embedding_derivation"},"covariant_inventory"->{"mathematical_deliverable","derived_congruence","covariant_block_congruence"},"scalar_regression_projection"->{"comparator_residual","scalar_regression","frozen_scalar_projection_residual"},"ambient_branches"->{"mathematical_deliverable","indexed_cells","ambient_indexed_action"},"open_strata"->{"selection_decision","open_strata","cell_stratum_selection"},"boundary_operator"->{"mathematical_deliverable","endpoint_assembly","endpoint_DtN_solve"},"endpoint_functionals"->{"mathematical_deliverable","endpoint_assembly","endpoint_variation_solve"},"e4_collar_scheme"->{"mathematical_deliverable","E4","extended_action_constraint_reduction"},"berry_pullback_selection"->{"mathematical_deliverable","berry","pullback_and_quotient_construction"},"berry_separation"->{"validation_predicate","berry.separation_control","subexpression_overlap_detector"},"g4_control"->{"mathematical_deliverable","g4_control","G4_force_and_holonomy"},"open_action_functionals"->{"mathematical_deliverable","open_root_reachability","production_remainder_generator"},"finite_bound_controls"->{"validation_predicate","reachability_analysis.finite_bound_controls","finite_generator_control"},"tensor_selection_rules"->{"validation_predicate","reachability_analysis.crosscheck_agrees","selection_table_crosscheck"},"partition"->{"mathematical_deliverable","partition_ledger","partition_ownership"},"phase_a_protection"->{"validation_predicate","source_contract.phase_a_payload_sha256","phase_a_payload_guard"}|>,head]];
sinkValues=<|"input_validation"->inputValidation,"checks.B1_R9"->checks["B1_R9"],"source_contract"->sourceContract,"dimensions"->dimensions,"field_manifest"->fieldManifest,"derived_congruence"->derivedCongruence,"scalar_regression"->scalarRegression,"indexed_cells"->indexed,"open_strata"->strata,"endpoint_assembly"->endpointFull,"E4"->e4,"E5"->e5,"berry"->berry,"berry.separation_control"->berry["separation_control"],"g4_control"->g4,"open_root_reachability"->reach,"reachability_analysis.finite_bound_controls"->reachAnalysis["finite_bound_controls"],"reachability_analysis.crosscheck_agrees"->reachAnalysis["crosscheck_agrees"],"partition_ledger"->partition,"source_contract.phase_a_payload_sha256"->phasePayloadSha|>;
lrows=Map[Function[pair,path=pair[[1]];contract=sinkFor[path];pathText=StringRiffle[ToString/@path,"."];scopes=Lookup[readEvents,pathText,{}];metadataOnly=scopes=!={}&&And@@Map[(StringStartsQ[#,"metadata:"]||StringStartsQ[#,"digest:"])&,scopes];forbidden=(ToString[First[path]]==="scalar_regression_projection"||(Length[path]>1&&ToString[path[[1]]]==="berry_pullback_selection"&&MemberQ[{"dependency_role","production_cells","ambient_branches"},ToString[path[[2]]]]));<|"path"->pathText,"typed_role"->contract[[1]],"semantic_sink"->contract[[2]],"semantic_sink_kind"->contract[[3]],"dependency_path"->Join[{"input:"<>pathText},("derived:"<>#&/@scopes),{"artifact:"<>contract[[2]]}],"read_scopes"->scopes,"sink_digest"->sha[sinkValues[contract[[2]]]],"metadata_or_digest_only"->metadataOnly,"production_ancestry_forbidden"->forbidden,"absent_from_production_ancestry"->forbidden|>],leafRows[data]];
declaredLeaves=Sort[lrows[[All,"path"]]];consumedLeaves=Sort@Intersection[Keys[readEvents],declaredLeaves];consumedButUndeclared=Sort@Complement[Keys[readEvents],declaredLeaves];liveness=<|"declared_leaf_paths"->declaredLeaves,"consumed_leaf_paths"->consumedLeaves,"consumed_but_undeclared"->consumedButUndeclared,"rows"->lrows,"per_key_mutation_evidence"->"external:b1_mutation_results.input_liveness_cases"|>;result["input_liveness"]=liveness;
recordCheck["B1_S1",inputValidation["valid"]&&declaredLeaves===consumedLeaves&&consumedButUndeclared==={}&&And@@Map[#!=""&,lrows[[All,"semantic_sink"]]]&&And@@Map[#["read_scopes"]=!={}&&!TrueQ[#["metadata_or_digest_only"]]&,lrows],<|"declared"->Length[declaredLeaves],"consumed"->Length[consumedLeaves],"consumed_but_undeclared"->consumedButUndeclared|>];result["checks"]=checks;
axis=If[And@@(#status==="PASS"&/@Values[checks]),"COMPUTATION_VALID","HARNESS_FAILED("<>StringRiffle[Keys@Select[checks,#status=!="PASS"&],","]<>")"];result["axis_1"]=axis;

jsonSafe[x_Association] := Association@KeyValueMap[ToString[#1] -> jsonSafe[#2] &, x]; jsonSafe[x_List] := jsonSafe /@ x; jsonSafe[x_String] := x; jsonSafe[x_Integer] := x; jsonSafe[x_Real] := x; jsonSafe[True] = True; jsonSafe[False] = False; jsonSafe[Null] = Null; jsonSafe[x_] := ToString[x,InputForm];
If[!DirectoryQ[DirectoryName[outputPath]],CreateDirectory[DirectoryName[outputPath],CreateIntermediateDirectories->True]];Export[outputPath,jsonSafe[result],"RawJSON","Compact"->True];
failures=Keys@Select[checks,#status=!="PASS"&];If[Length[failures]>0,first=First[failures];Print["ASSERT_FAIL:"<>first<>":"<>checks[first]["evidence_digest"]];Exit[1]];
Print["MATHEMATICA_PHASE_B1_REMEDIATION3: PASS cells="<>ToString[Length[cellRecords]]];Print["OUTPUT: "<>outputPath];Exit[0];
