SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
(* S12 decision list, independent Wolfram construction. No external imports.
   Supplied input: the two conservation laws, signed conversion identity,
   and U1's live interior field between two live faces. All subsequent physical
   expressions are differentiation, elimination, restriction, or extraction.
   Distribution convention: n and u retain their original differentiable
   bulk-field class. J has independent differentiable restrictions on the
   three open regions, independent one-sided traces, and an arbitrary
   distribution supported at each face. Derivatives of those distributions
   stay formal; no face-flux form, continuity, or boundary condition is imposed.
   The order can jump. No point value at a face is selected.
   This script computes and prints; it has no expected-value acceptance guard. *)
$HistoryLength = 0;
transcriptPath = "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S12_order_source.out";
stream = OpenWrite[transcriptPath, PageWidth -> Infinity];
emit[tag_String, object_] := Module[{line},
  line = tag <> ": " <> ToString[object, InputForm, PageWidth -> Infinity];
  Print[line]; WriteString[stream, line, "\n"];
];
object[tag_String, value_] := (
  emit[tag, value]; emit[tag <> "/comparison", "there is nothing to compare it with"];
);
working[tag_String] := emit["WORKING_ON", tag];

(* Dimension vectors use {length,time}; material count is dimensionless.
   Only primitive meanings are assigned. The order-flux units are solved
   from the supplied order balance, rather than prescribed. *)
unit[fun_] := Switch[fun,
  n, {-4, 0}, chi | f | g, {0, 0},
  um | u1 | u2 | u3, {1, -1},
  wm | wp, {1, 0}, gb | gr | gd, {0, -1},
  j1 | j2 | j3 | jm, fluxUnit,
  _, If[MemberQ[fluxFunctionHeads, fun], fluxUnit,
    Throw[{"unassigned_function_dimension", HoldForm[fun]}]]
];
dimension[e_?NumberQ] := {0, 0};
dimension[t] := {0, 1};
dimension[x1 | x2 | x3 | w] := {1, 0};
dimension[e_Plus] := Module[{parts, dims, ref},
  parts = List @@ e; dims = dimension /@ parts; ref = First[dims];
  AppendTo[dimensionLedger, {parts, dims, (# - ref & /@ dims)}]; ref
];
dimension[e_Times] := Total[dimension /@ (List @@ e)];
dimension[Power[e_, p_?NumberQ]] := p dimension[e];
dimension[HeavisideTheta[e_]] := (dimension[e]; {0, 0});
dimension[DiracDelta[e_]] := -dimension[e];
dimension[Derivative[orders__][fun_][args___]] :=
  unit[fun] - Total[MapThread[Times, {{orders}, dimension /@ {args}}]];
dimension[fun_[args___]] := unit[fun];
dimension[e_] := Throw[{"unassigned_expression_dimension", HoldForm[e]}];

dimReport[tag_String, expression_, reference_:None] := Module[{pieces, dims, ref},
  working[tag <> "/dimensions"]; dimensionLedger = {};
  pieces = If[Head[Expand[expression]] === Plus, List @@ Expand[expression], {Expand[expression]}];
  dims = dimension /@ pieces;
  emit[tag <> "/dimension_operands", Transpose[{pieces, dims}]];
  emit[tag <> "/nested_dimension_operands", Take[#, 2] & /@ dimensionLedger];
  emit[tag <> "/nested_dimension_residual", (#[[3]] & /@ dimensionLedger)];
  If[reference === None,
    emit[tag <> "/dimension_comparison", "there is nothing to compare it with"],
    ref = dimension[reference];
    emit[tag <> "/dimension_reference_operand", {reference, ref}];
    emit[tag <> "/dimension_residual", (# - ref & /@ dims)]
  ]; dims
];

main[] := Module[{space, coordinates, density, velocity, flux, fraction,
    faces, levels, assumptions, regions, material, order, conversion,
    timeRule, fluxUnknowns, dimTerms, dimVectors, dimReference, dimResidual,
    fluxRules, packages, label, orderField, raw, source, regular,
    restricted, delta, coefficient, trace, faceTerm, fullFaces, branchAssumptions,
    regionIndicators, regionFluxHeads, faceFluxHeads, regionFlux, faceFlux,
    markers, markerRules, markedFlux, fluxReplacement, toField, markedSource,
    supportedSources},
  working["supplied_inputs"];
  emit["engine", {$Version, $SystemID}];
  emit["governing_input", "S12_decision_list.md; supplied equations and U1 encoding are unfalsifiable inputs"];
  emit["calculus_scope", {"distributional order; n and u retain differentiable bulk fields with traces",
    "J has independent restrictions and one-sided traces in all three open regions",
    "J also retains arbitrary face-supported distributions; their derivatives are distributional and remain formal",
    "no continuity, exterior vanishing, or boundary condition on J",
    "pointwise face values unspecified", "interior range and strict face order are supplied domain conditions"}];
  coordinates = {t, x1, x2, x3, w}; space = Rest[coordinates];
  density = n @@ coordinates;
  velocity = (#[Sequence @@ coordinates] & /@ {u1, u2, u3, um});
  flux = (#[Sequence @@ coordinates] & /@ {j1, j2, j3, jm});
  fraction = chi @@ coordinates;
  faces = {wm[t, x1, x2, x3], wp[t, x1, x2, x3]};
  (* BEGIN SUPPLIED BALANCES *)
  material = D[density, t] + Total[MapThread[D[#1, #2] &, {density velocity, space}]] == 0;
  order = D[fraction density, t] + Total[MapThread[D[#1, #2] &, {fraction density velocity + flux, space}]] == density (gb @@ coordinates);
  conversion = (gb @@ coordinates) == (gr @@ coordinates) - (gd @@ coordinates);
  (* END SUPPLIED BALANCES *)
  (* BEGIN U1 ENCODING. The difference of steps describes the same two
     strictly ordered faces without selecting a point value at either face. *)
  levels = (w - # & /@ faces);
  assumptions = Element[coordinates, Reals] && (Less @@ faces);
  regions = {w < First[faces], First[faces] < w < Last[faces], w > Last[faces]};
  packages = {
    {"general_order", (f @@ coordinates) (Subtract @@ (HeavisideTheta /@ levels))},
    {"conditional_fast_ordering", g[w - First[faces], Last[faces] - w] (Subtract @@ (HeavisideTheta /@ levels))}
  };
  (* END U1 ENCODING *)
  (* Partition the unrestricted flux using the already computed ordered-face
     indicators. These are restrictions of J, not a constitutive profile.
     The additional face fields denote arbitrary supported distributions,
     not values or traces of the regional fields. Markers isolate their
     differentiated contributions algebraically, without setting a field to 0. *)
  regionIndicators = -Differences[Join[{1}, HeavisideTheta /@ levels, {0}]];
  regionFluxHeads = Table[Symbol["j" <> ToString[c] <> "Region" <> ToString[r]],
    {r, Length[regions]}, {c, Length[space]}];
  faceFluxHeads = Table[Symbol["j" <> ToString[c] <> "Face" <> ToString[k]],
    {k, Length[faces]}, {c, Length[space]}];
  fluxFunctionHeads = Join[Flatten[regionFluxHeads], Flatten[faceFluxHeads]];
  regionFlux = Map[(# @@ coordinates) &, regionFluxHeads, {2}];
  faceFlux = Map[(# @@ coordinates) &, faceFluxHeads, {2}];
  markers = Table[Unique["faceFluxMarker$"], {Length[faces]}];
  markerRules = Thread[markers -> ConstantArray[1, Length[faces]]];
  markedFlux = regionIndicators . regionFlux + markers . faceFlux;
  toField[expr_] := Function[{tt, xx1, xx2, xx3, ww},
    Evaluate[expr /. Thread[coordinates -> {tt, xx1, xx2, xx3, ww}]]];
  fluxReplacement = Thread[{j1, j2, j3, jm} -> (toField /@ markedFlux)];
  object["supplied/material_balance", material];
  object["supplied/order_balance", order];
  object["supplied/conversion_identity", conversion];
  object["supplied/faces", faces];
  object["supplied/domain", {assumptions, 0 <= (f @@ coordinates) <= 1}];
  object["live_fields", {density, velocity, flux, fraction, faces}];
  object["region_domains", regions];
  object["flux/region_indicators", regionIndicators];
  object["flux/region_fields", regionFlux];
  object["flux/face_supported_fields", faceFlux];
  object["flux/face_support_domains", Table[{faceFlux[[k]], w == faces[[k]]}, {k, Length[faces]}]];
  object["flux/distribution_field", markedFlux /. markerRules];
  working["material_elimination"];
  timeRule = First[Solve[material, D[density, t]]];
  object["material/time_derivative_rule", timeRule];

  working["order_flux_dimensions"];
  emit["dimension_basis", {"length", "time", "number count dimensionless"}];
  fluxUnknowns = {jLength, jTime}; fluxUnit = fluxUnknowns;
  dimensionLedger = {};
  dimTerms = List @@ Expand[First[order]];
  dimVectors = dimension /@ dimTerms;
  dimReference = dimension[Last[order]];
  emit["order_flux_dimension_inference/operands", {Transpose[{dimTerms, dimVectors}], {Last[order], dimReference}}];
  dimResidual = (# - dimReference & /@ dimVectors);
  emit["order_flux_dimension_inference/residual", dimResidual];
  fluxRules = Solve[Thread[DeleteDuplicates[Flatten[dimResidual]] == 0], fluxUnknowns];
  object["order_flux_dimension_inference/rules", fluxRules];
  fluxUnit = fluxUnknowns /. Flatten[fluxRules];
  object["order_flux_dimension_inference/object", fluxUnit];
  dimReport["supplied/material_balance_lhs", First[material], D[density, t]];
  object["supplied/material_balance_rhs", Last[material]];
  emit["supplied/material_balance_rhs/dimension", "additive identity has no intrinsic unit"];
  dimReport["supplied/order_balance_lhs", First[order], Last[order]];
  dimReport["supplied/order_balance_rhs", Last[order]];
  dimReport["supplied/conversion_identity_rhs", Last[conversion], First[conversion]];

  Do[
    label = package[[1]]; orderField = package[[2]];
    working[label <> "/order_field"];
    object[label <> "/order_field", orderField];
    If[label === "conditional_fast_ordering",
      emit[label <> "/premise", "conditional specialization: ordering keeps pace with the slab; G is general in both offsets"];
      object[label <> "/interior_range", 0 <= g[w - First[faces], Last[faces] - w] <= 1]
    ];
    working[label <> "/source"];
    (* Substitution into the supplied balance precedes all eliminations. *)
    raw = Expand[First[order] /. Join[{chi -> toField[orderField]}, fluxReplacement]];
    object[label <> "/balance_source_operand", raw];
    object[label <> "/material_elimination_operand", timeRule];
    markedSource = Expand[raw /. timeRule];
    object[label <> "/source_with_flux_markers", markedSource];
    supportedSources = (Coefficient[markedSource, #] & /@ markers);
    object[label <> "/face_supported_flux_sources", supportedSources];
    source = Expand[markedSource /. markerRules];
    object[label <> "/source_distribution", source];
    regular = Expand[source - Total[supportedSources]] /. DiracDelta[_] -> 0;
    object[label <> "/regular_distribution_part", regular];
    dimReport[label <> "/order_field", orderField];
    dimReport[label <> "/source_distribution", source, Last[order]];
    dimReport[label <> "/regular_distribution_part", regular, Last[order]];
    Do[
      working[label <> "/region/" <> ToString[k]];
      branchAssumptions = assumptions && regions[[k]];
      object[label <> "/region/" <> ToString[k] <> "/domain", branchAssumptions];
      restricted = Refine[regular, branchAssumptions];
      object[label <> "/region/" <> ToString[k] <> "/source", restricted];
      dimReport[label <> "/region/" <> ToString[k] <> "/source", restricted, Last[order]],
      {k, Length[regions]}];
    fullFaces = {};
    Do[
      working[label <> "/face/" <> ToString[k]];
      delta = DiracDelta[levels[[k]]];
      coefficient = Coefficient[source, delta];
      trace = Expand[coefficient /. w -> faces[[k]]];
      faceTerm = trace delta + supportedSources[[k]];
      object[label <> "/face/" <> ToString[k] <> "/locus", w == faces[[k]]];
      object[label <> "/face/" <> ToString[k] <> "/delta", delta];
      object[label <> "/face/" <> ToString[k] <> "/coefficient_operand", coefficient];
      object[label <> "/face/" <> ToString[k] <> "/trace_coefficient", trace];
      object[label <> "/face/" <> ToString[k] <> "/supported_flux_operand", faceFlux[[k]]];
      object[label <> "/face/" <> ToString[k] <> "/supported_flux_source", supportedSources[[k]]];
      object[label <> "/face/" <> ToString[k] <> "/source_contribution", faceTerm];
      dimReport[label <> "/face/" <> ToString[k] <> "/trace_coefficient", trace];
      dimReport[label <> "/face/" <> ToString[k] <> "/supported_flux_source", supportedSources[[k]], Last[order]];
      dimReport[label <> "/face/" <> ToString[k] <> "/source_contribution", faceTerm, Last[order]];
      AppendTo[fullFaces, faceTerm],
      {k, Length[faces]}];
    object[label <> "/source_with_face_traces", regular + Total[fullFaces]],
    {package, packages}];
  emit["completion", {"general_order", "conditional_fast_ordering", "dimensional_objects"}];
];
result = Catch[Check[main[]; "completed", {"wolfram_messages", $MessageList}]];
emit["execution", result];
Close[stream];
If[result =!= "completed", Exit[1]];
Exit[0];
