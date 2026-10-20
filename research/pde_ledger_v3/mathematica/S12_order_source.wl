SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
(* S12 decision list, independent Wolfram construction. No external imports.
   Supplied input: the two conservation laws, signed conversion identity,
   and U1's live interior field between two live faces. All subsequent physical
   expressions are differentiation, elimination, restriction, or extraction.
   Revision 4: n, u and J have independent differentiable restrictions and
   one-sided traces on the three open regions. The order is live inside,
   zero outside. No input field has a face-concentrated part. Conserved
   products are formed within a region before applying the region indicators.
   Both balances are differentiated distributionally, including the material
   balance at each face. No point value or jump is selected.
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
  _, Which[
    MemberQ[densityFunctionHeads, fun], {-4, 0},
    MemberQ[velocityFunctionHeads, fun], {1, -1},
    MemberQ[fluxFunctionHeads, fun], fluxUnit,
    True, Throw[{"unassigned_function_dimension", HoldForm[fun]}]]
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
    materialDensity, materialCurrent, orderDensity, orderCurrent,
    fluxUnknowns, dimTerms, dimVectors, dimReference, dimResidual, fluxRules,
    packages, label, orderField, interiorOrder, regionalOrders, raw, source,
    regular, restricted, delta, coefficient, trace, faceTerm, fullFaces,
    regionIndicators, regionDensityHeads, regionVelocityHeads, regionFluxHeads,
    regionDensity, regionVelocity, regionFlux, fieldRules, toField,
    materialDensities, materialCurrents, localMaterial, timeRules,
    materialDensityField, materialCurrentField, materialRaw, materialRequired,
    materialFaceEquations, materialFaceTerms, orderDensities, orderCurrents,
    orderDensityField, orderCurrentField, localOrder, fieldTraces},
  working["supplied_inputs"];
  emit["engine", {$Version, $SystemID}];
  emit["governing_input", "S12_decision_list.md revision 4, b08673a2; supplied equations and face rules are unfalsifiable inputs"];
  emit["calculus_scope", {
    "n, u and J are independent live functions in all three open regions with unrestricted one-sided traces",
    "order is a live interior field, zero outside by U1",
    "input fields have no face-concentrated parts; those possibilities remain OPEN outside this object",
    "the computed source may have face-concentrated parts",
    "pointwise face values and jumps are unspecified"}];
  coordinates = {t, x1, x2, x3, w}; space = Rest[coordinates];
  density = n @@ coordinates;
  velocity = (#[Sequence @@ coordinates] & /@ {u1, u2, u3, um});
  flux = (#[Sequence @@ coordinates] & /@ {j1, j2, j3, jm});
  fraction = chi @@ coordinates;
  faces = {wm[t, x1, x2, x3], wp[t, x1, x2, x3]};
  (* BEGIN SUPPLIED BALANCES. Retain the undifferentiated conserved operands
     so products can be restricted before distributional differentiation. *)
  materialDensity = density; materialCurrent = density velocity;
  orderDensity = fraction density; orderCurrent = fraction density velocity + flux;
  material = D[materialDensity, t] + Total[MapThread[D[#1, #2] &, {materialCurrent, space}]] == 0;
  order = D[orderDensity, t] + Total[MapThread[D[#1, #2] &, {orderCurrent, space}]] == density (gb @@ coordinates);
  conversion = (gb @@ coordinates) == (gr @@ coordinates) - (gd @@ coordinates);
  (* END SUPPLIED BALANCES *)
  (* BEGIN U1 ENCODING *)
  levels = (w - # & /@ faces);
  assumptions = Element[coordinates, Reals] && (Less @@ faces);
  regions = {w < First[faces], First[faces] < w < Last[faces], w > Last[faces]};
  packages = {
    {"general_order", f @@ coordinates},
    {"conditional_fast_ordering", g[w - First[faces], Last[faces] - w]}
  };
  (* END U1 ENCODING *)
  regionIndicators = -Differences[Join[{1}, HeavisideTheta /@ levels, {0}]];
  regionDensityHeads = Table[Symbol["nRegion" <> ToString[r]], {r, Length[regions]}];
  regionVelocityHeads = Table[Symbol["u" <> ToString[c] <> "Region" <> ToString[r]],
    {r, Length[regions]}, {c, Length[space]}];
  regionFluxHeads = Table[Symbol["j" <> ToString[c] <> "Region" <> ToString[r]],
    {r, Length[regions]}, {c, Length[space]}];
  densityFunctionHeads = regionDensityHeads;
  velocityFunctionHeads = Flatten[regionVelocityHeads];
  fluxFunctionHeads = Flatten[regionFluxHeads];
  regionDensity = (# @@ coordinates & /@ regionDensityHeads);
  regionVelocity = Map[(# @@ coordinates) &, regionVelocityHeads, {2}];
  regionFlux = Map[(# @@ coordinates) &, regionFluxHeads, {2}];
  toField[expr_] := Function[{tt, xx1, xx2, xx3, ww},
    Evaluate[expr /. Thread[coordinates -> {tt, xx1, xx2, xx3, ww}]]];
  fieldRules = Table[
    Thread[Join[{n}, {u1, u2, u3, um}, {j1, j2, j3, jm}] ->
      (toField /@ Join[{regionDensity[[r]]}, regionVelocity[[r]], regionFlux[[r]]])],
    {r, Length[regions]}];
  object["supplied/material_balance", material];
  object["supplied/order_balance", order];
  object["supplied/conversion_identity", conversion];
  object["supplied/faces", faces];
  object["supplied/domain", {assumptions, 0 <= (f @@ coordinates) <= 1}];
  object["live_fields", {density, velocity, flux, fraction, faces}];
  object["region_domains", regions];
  object["fields/region_indicators", regionIndicators];
  object["density/region_fields", regionDensity];
  object["velocity/region_fields", regionVelocity];
  object["flux/region_fields", regionFlux];
  object["density/distribution_field", regionIndicators . regionDensity];
  object["velocity/distribution_field", regionIndicators . regionVelocity];
  object["flux/distribution_field", regionIndicators . regionFlux];

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

  working["material/region_balances"];
  materialDensities = Table[materialDensity /. fieldRules[[r]], {r, Length[regions]}];
  materialCurrents = Table[materialCurrent /. fieldRules[[r]], {r, Length[regions]}];
  localMaterial = Table[material /. fieldRules[[r]], {r, Length[regions]}];
  timeRules = Flatten[Table[First[Solve[localMaterial[[r]], D[regionDensity[[r]], t]]],
    {r, Length[regions]}]];
  object["material/region_density_operands", materialDensities];
  object["material/region_current_operands", materialCurrents];
  object["material/region_balance_operands", localMaterial];
  object["material/time_derivative_rules", timeRules];
  materialDensityField = regionIndicators . materialDensities;
  materialCurrentField = regionIndicators . materialCurrents;
  object["material/conserved_density_operand", materialDensityField];
  object["material/conserved_current_operand", materialCurrentField];
  materialRaw = Expand[D[materialDensityField, t] +
    Total[MapThread[D[#1, #2] &, {materialCurrentField, space}]]];
  object["material/balance_distribution_operand", materialRaw];
  materialRequired = Expand[materialRaw /. timeRules];
  object["material/balance_distribution", materialRequired];
  dimReport["material/balance_distribution_operand", materialRaw, D[density, t]];
  dimReport["material/balance_distribution", materialRequired, D[density, t]];
  Do[
    object["material/region/" <> ToString[r] <> "/balance", localMaterial[[r]]];
    dimReport["material/region/" <> ToString[r] <> "/balance_lhs", First[localMaterial[[r]]], D[density, t]],
    {r, Length[regions]}];
  materialFaceEquations = {}; materialFaceTerms = {};
  Do[
    working["material/face/" <> ToString[k]];
    delta = DiracDelta[levels[[k]]];
    coefficient = Coefficient[materialRaw, delta];
    trace = Expand[coefficient /. w -> faces[[k]]];
    faceTerm = trace delta;
    object["material/face/" <> ToString[k] <> "/coefficient_operand", coefficient];
    object["material/face/" <> ToString[k] <> "/trace_coefficient", trace];
    emit["material/face/" <> ToString[k] <> "/balance_operands", {trace, Last[material]}];
    emit["material/face/" <> ToString[k] <> "/balance_residual", trace - Last[material]];
    AppendTo[materialFaceEquations, trace == Last[material]];
    object["material/face/" <> ToString[k] <> "/requirement", Last[materialFaceEquations]];
    object["material/face/" <> ToString[k] <> "/balance_contribution", faceTerm];
    dimReport["material/face/" <> ToString[k] <> "/trace_coefficient", trace];
    dimReport["material/face/" <> ToString[k] <> "/balance_contribution", faceTerm, D[density, t]];
    AppendTo[materialFaceTerms, faceTerm],
    {k, Length[faces]}];
  object["material/face_requirements", materialFaceEquations];
  object["material/balance_with_face_traces", Total[materialFaceTerms]];

  Do[
    label = package[[1]]; interiorOrder = package[[2]];
    (* Regional order values are exactly U1: only the interior is live. *)
    regionalOrders = {0, interiorOrder, 0};
    orderField = regionIndicators . regionalOrders;
    working[label <> "/order_field"];
    object[label <> "/order_field", orderField];
    object[label <> "/region_order_operands", regionalOrders];
    If[label === "conditional_fast_ordering",
      emit[label <> "/premise", "conditional specialization: ordering keeps pace with the slab; G is general in both offsets"];
      object[label <> "/interior_range", 0 <= interiorOrder <= 1]
    ];
    working[label <> "/conserved_operands"];
    orderDensities = Table[orderDensity /. Join[fieldRules[[r]], {chi -> toField[regionalOrders[[r]]]}],
      {r, Length[regions]}];
    orderCurrents = Table[orderCurrent /. Join[fieldRules[[r]], {chi -> toField[regionalOrders[[r]]]}],
      {r, Length[regions]}];
    object[label <> "/region_density_operands", orderDensities];
    object[label <> "/region_current_operands", orderCurrents];
    orderDensityField = regionIndicators . orderDensities;
    orderCurrentField = regionIndicators . orderCurrents;
    object[label <> "/conserved_density_operand", orderDensityField];
    object[label <> "/conserved_current_operand", orderCurrentField];
    working[label <> "/source"];
    raw = Expand[D[orderDensityField, t] +
      Total[MapThread[D[#1, #2] &, {orderCurrentField, space}]]];
    object[label <> "/balance_source_operand", raw];
    object[label <> "/material_elimination_operand", timeRules];
    source = Expand[raw /. timeRules];
    object[label <> "/source_distribution", source];
    regular = source /. DiracDelta[_] -> 0;
    object[label <> "/regular_distribution_part", regular];
    dimReport[label <> "/order_field", orderField];
    dimReport[label <> "/balance_source_operand", raw, Last[order]];
    dimReport[label <> "/source_distribution", source, Last[order]];
    dimReport[label <> "/regular_distribution_part", regular, Last[order]];
    localOrder = Table[Expand[D[orderDensities[[r]], t] +
      Total[MapThread[D[#1, #2] &, {orderCurrents[[r]], space}]]], {r, Length[regions]}];
    Do[
      working[label <> "/region/" <> ToString[r]];
      object[label <> "/region/" <> ToString[r] <> "/domain", assumptions && regions[[r]]];
      object[label <> "/region/" <> ToString[r] <> "/field_operands",
        {regionDensity[[r]], regionVelocity[[r]], regionFlux[[r]], regionalOrders[[r]]}];
      restricted = Expand[localOrder[[r]] /. timeRules];
      object[label <> "/region/" <> ToString[r] <> "/source", restricted];
      dimReport[label <> "/region/" <> ToString[r] <> "/source", restricted, Last[order]],
      {r, Length[regions]}];
    fullFaces = {};
    Do[
      working[label <> "/face/" <> ToString[k]];
      delta = DiracDelta[levels[[k]]];
      coefficient = Coefficient[source, delta];
      trace = Expand[coefficient /. w -> faces[[k]]];
      faceTerm = trace delta;
      fieldTraces = Table[
        {regionDensity[[r]], regionVelocity[[r]], regionFlux[[r]], regionalOrders[[r]]} /. w -> faces[[k]],
        {r, k, k + 1}];
      object[label <> "/face/" <> ToString[k] <> "/locus", w == faces[[k]]];
      object[label <> "/face/" <> ToString[k] <> "/field_trace_operands", fieldTraces];
      object[label <> "/face/" <> ToString[k] <> "/delta", delta];
      object[label <> "/face/" <> ToString[k] <> "/coefficient_operand", coefficient];
      object[label <> "/face/" <> ToString[k] <> "/trace_coefficient", trace];
      object[label <> "/face/" <> ToString[k] <> "/source_contribution", faceTerm];
      object[label <> "/face/" <> ToString[k] <> "/material_requirement", materialFaceEquations[[k]]];
      dimReport[label <> "/face/" <> ToString[k] <> "/trace_coefficient", trace];
      dimReport[label <> "/face/" <> ToString[k] <> "/source_contribution", faceTerm, Last[order]];
      AppendTo[fullFaces, faceTerm],
      {k, Length[faces]}];
    object[label <> "/source_with_face_traces", regular + Total[fullFaces]],
    {package, packages}];
  emit["completion", {"material_face_requirements", "general_order", "conditional_fast_ordering", "dimensional_objects"}];
];
result = Catch[Check[main[]; "completed", {"wolfram_messages", $MessageList}]];
emit["execution", result];
Close[stream];
If[result =!= "completed", Exit[1]];
Exit[0];
