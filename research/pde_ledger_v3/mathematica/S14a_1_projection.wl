SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
$HistoryLength = 0;
(* S14a part 1. Supplied inputs are confined to suppliedData. No imports.
   Construction: pull back the spacetime current to a fixed slab coordinate;
   integrate its transverse derivative; differentiate HeavisideTheta regional
   currents separately, counting each boundary delta with unit weight.
   Formal integrals retain arbitrary native fields. No symbolic-value guards.
   Equalities printed with Inactive[Equal] denote the supplied balance imposed
   on computed operands; they are not asserted identities. *)
ClearAll[emit, working, alone, operands];
emit[tag_, value_] := Print[tag, ": ", ToString[value, InputForm, PageWidth -> Infinity]];
working[tag_] := emit["WORKING_ON", tag];
alone[tag_, value_] := (emit[tag, value]; emit[tag <> "/comparison", "There is nothing to compare this object with."]);
operands[tag_, a_, b_] := (emit[tag <> "/operand_a", a]; emit[tag <> "/operand_b", b]; emit[tag <> "/residual", Cancel[Together[Expand[normalizeArguments[a - b]]]]]);

(* Dimensions are exponent vectors in {material number, length, time}.
   They are assignments, not measured results. *)
fieldDimensions = <||>;
coordinateDimensions = <|t -> {0, 0, 1}, x1 -> {0, 1, 0}, x2 -> {0, 1, 0}, x3 -> {0, 1, 0}, w -> {0, 1, 0}, r -> {0, 1, 0}, ri -> {0, 1, 0}, ro -> {0, 1, 0}, z -> {0, 0, 0}|>;
register[heads_, dim_] := Scan[(AssociateTo[fieldDimensions, # -> dim]) &, Flatten[{heads}]];
makeHeads[stem_, count_] := Table[Symbol[stem <> ToString[k]], {k, count}];
ClearAll[suppliedData];
suppliedData[spherical_] := Module[{args, faceArgs, spatial, ns, us, js, cs, gs, grs, gds, lower, upper, density, flux, rhs, balance, orderIdentity, factor, count},
 (* SUPPLIED EQUATIONS, U1 ENCODING, FIELD ARGUMENTS AND FACE FACTOR.
    Each region has independent one-sided fields; the only absent order
    fields are those set absent outside the slab by U1 itself. *)
 If[spherical,
  args = {r, w, t}; faceArgs = {r, t}; spatial = {r, w}; count = 2;
  ns = makeHeads["sn", 3]; us = Table[makeHeads["su" <> ToString[k] <> "c", count], {k, 3}];
  js = Table[makeHeads["sj" <> ToString[k] <> "c", count], {k, 3}];
  cs = schi; gs = makeHeads["sgamma", 3]; grs = makeHeads["sreturn", 3]; gds = makeHeads["sdrain", 3];
  lower = swLower @@ faceArgs; upper = swUpper @@ faceArgs,
  args = {x1, x2, x3, w, t}; faceArgs = {x1, x2, x3, t}; spatial = {x1, x2, x3, w}; count = 4;
  ns = makeHeads["n", 3]; us = Table[makeHeads["u" <> ToString[k] <> "c", count], {k, 3}];
  js = Table[makeHeads["j" <> ToString[k] <> "c", count], {k, 3}];
  cs = chi; gs = makeHeads["gamma", 3]; grs = makeHeads["return", 3]; gds = makeHeads["drain", 3];
  lower = wLower @@ faceArgs; upper = wUpper @@ faceArgs];
 register[ns, {1, -4, 0}]; register[us, {0, 1, -1}]; register[js, {1, -3, -1}];
 register[cs, {0, 0, 0}]; register[{gs, grs, gds}, {0, 0, -1}]; register[{Head[lower], Head[upper]}, {0, 1, 0}];
 ns = (Apply[#, args] &) /@ ns; us = Map[Apply[#, args] &, us, {2}]; js = Map[Apply[#, args] &, js, {2}];
 gs = (Apply[#, args] &) /@ gs; grs = (Apply[#, args] &) /@ grs; gds = (Apply[#, args] &) /@ gds;
 cs = {0, cs @@ args, 0};
 density = {ns, cs ns};
 flux = {MapThread[Times, {ns, us}], MapThread[Plus, {MapThread[Times, {cs ns, us}], js}]};
 rhs = {ConstantArray[0, 3], ns gs};
 orderIdentity = MapThread[Inactive[Equal], {gs, grs - gds}];
 factor = Table[Last[us[[k]]] - D[face, t] - Most[us[[k]]].(D[face, #] & /@ Most[spatial]), {k, 3}, {face, {lower, upper}}];
 <|"coords" -> Prepend[spatial, t], "faces" -> {lower, upper}, "density" -> density, "flux" -> flux, "rhs" -> rhs, "gamma_identity" -> orderIdentity, "face_factor" -> factor, "order" -> cs, "field_arguments" -> args|>
];

(* Geometry: the angular measure is computed from the ordinary spherical
   coordinate map, independent of the physical fields. *)
working["spherical specialization / geometric measure"];
spatialMap = {r Sin[theta] Cos[phi], r Sin[theta] Sin[phi], r Cos[theta]};
spatialJacobian = Table[D[spatialMap[[i]], v], {i, 3}, {v, {r, theta, phi}}];
spatialJacobianDet = FullSimplify[Det[spatialJacobian]];
angularDensity = FullSimplify[Abs[spatialJacobianDet], Assumptions -> r > 0 && 0 < theta < Pi];
sphereMeasure = Integrate[angularDensity, {theta, 0, Pi}, {phi, 0, 2 Pi}, Assumptions -> r > 0];
angularConstant = Cancel[sphereMeasure/r^2];
alone["spherical specialization / spatial Jacobian", spatialJacobian];
alone["spherical specialization / sphere measure", sphereMeasure];

(* Compute the one-dimensional integration operator on an arbitrary current.
   Substitution of a computed current is done only after this integration. *)
working["transverse integration operator"];
ftcAntiderivative = Integrate[Derivative[1][primitive][q], q];
ftcTemplate = (ftcAntiderivative /. q -> bb) - (ftcAntiderivative /. q -> aa);
emit["transverse integration / operand", Inactive[Integrate][Derivative[1][primitive][q], {q, aa, bb}]];
emit["transverse integration / object", ftcTemplate];
emit["transverse integration / comparison", "There is nothing to compare this object with."];
ClearAll[boundaryIntegral, boundaryContributions, formalIntegral, formalD, additiveTerms, normalizeArguments];
normalizeArguments[expr_] := expr /. {
 Derivative[orders__][head_][args__] :> Apply[Derivative[orders][head], (Expand[Cancel[Together[#]]] &) /@ {args}],
 (head_[args__] /; KeyExistsQ[fieldDimensions, head]) :> Apply[head, (Expand[Cancel[Together[#]]] &) /@ {args}]
};
boundaryContributions[expr_, var_, a_, b_] := Module[{template, slots = {lowerEndpointSlot, upperEndpointSlot}, ends = {a, b}},
 template = ftcTemplate /. {primitive[aa] -> First[slots], primitive[bb] -> Last[slots]};
 Table[Coefficient[template, slots[[k]]] (expr /. var -> ends[[k]]), {k, 2}]
];
boundaryIntegral[expr_, var_, a_, b_] := Total[boundaryContributions[expr, var, a, b]];
additiveTerms[expr_] := If[Head[expr] === Plus, List @@ expr, {expr}];
formalIntegral[expr_, var_, a_, b_] := Total[Map[Function[term,
 Module[{factors, outside, inside},
  factors = If[Head[term] === Times, List @@ term, {term}];
  outside = Times @@ Select[factors, FreeQ[#, var] &];
  inside = Times @@ Select[factors, !FreeQ[#, var] &];
  outside Inactive[Integrate][inside, {var, a, b}]]], additiveTerms[Expand[expr]]]];
(* Differentiation of formal integrals uses the chain rule and the
   fundamental theorem, retaining all moving endpoint contributions. *)
formalD[expr_, var_] := Module[{entries = {}, replaced, token, ans},
 replaced = expr /. int:Inactive[Integrate][body_, {v_, a_, b_}] :>
   (token = Unique["integralSlot"]; AppendTo[entries, {token, body, v, a, b}]; token);
 ans = D[replaced, var];
 ans = ans + Total[Function[entry, D[replaced, entry[[1]]] (
   formalIntegral[formalD[entry[[2]], var], entry[[3]], entry[[4]], entry[[5]]] +
   (entry[[2]] /. entry[[3]] -> entry[[5]]) D[entry[[5]], var] -
   (entry[[2]] /. entry[[3]] -> entry[[4]]) D[entry[[4]], var]) ] /@ entries];
 ans /. (Function[entry, entry[[1]] -> Inactive[Integrate][entry[[2]], {entry[[3]], entry[[4]], entry[[5]]}]] /@ entries)
];

(* A transverse delta at an endpoint is integrated as its coefficient, not
   with a half weight. This is the contract's closed-region measure. *)
ClearAll[projection];
projection[data_, measure_, prefix_] := Module[
 {coords = data["coords"], faces = data["faces"], transverse, parallel, chart, map, jac, jacDet, inverse, masks, maskGradient, currents, nativeDiv, packed, piola, piolaDiv, density, planeFlux, transverseFace, atoms, concentrated, regular, closedFace, left, required, suppliedSource, bulkSource, projectedDerivative, closedOperand, result = <||>, eq, k, s, j, d, pulledNative, pulledDiv, sourceResidual, integratePullback, nativeIntegral, partial, directIntegrated},
 transverse = Last[coords]; parallel = Most[coords]; d = Length[coords];
 chart = Append[parallel, z];
 map = Append[parallel, First[faces] + z (Last[faces] - First[faces])];
 jac = Table[D[map[[k]], chart[[j]]], {k, d}, {j, d}];
 jacDet = Factor[Det[jac]]; inverse = Simplify[Inverse[jac]];
 masks = {1 - HeavisideTheta[transverse - First[faces]], HeavisideTheta[transverse - First[faces]] - HeavisideTheta[transverse - Last[faces]], HeavisideTheta[transverse - Last[faces]]};
 maskGradient = Table[D[masks[[k]], coords[[j]]], {k, 3}, {j, d}];
 alone[prefix <> "/native field arguments", data["field_arguments"]];
 alone[prefix <> "/U1 order fields", data["order"]];
 alone[prefix <> "/face ordering", Inactive[Less] @@ faces];
 alone[prefix <> "/supplied gamma identity", data["gamma_identity"]];
 alone[prefix <> "/supplied material face factor", data["face_factor"]];
 alone[prefix <> "/region masks", masks];
 alone[prefix <> "/differentiated masks", maskGradient];
 alone[prefix <> "/coordinate map", map];
 alone[prefix <> "/coordinate Jacobian", jac];
 alone[prefix <> "/coordinate Jacobian determinant", jacDet];
 nativeIntegral[expr_] := formalIntegral[expr, transverse, First[faces], Last[faces]];
 integratePullback[expr_] := nativeIntegral[Cancel[normalizeArguments[(expr /. z -> (transverse - First[faces])/(Last[faces] - First[faces]))/jacDet]]];
 Do[
  eq = {"material", "order"}[[s]];
  working[prefix <> "/" <> eq <> "/regional currents and distribution"];
  currents = MapThread[Prepend, {data["flux"][[s]], data["density"][[s]]}];
  nativeDiv = Table[Total[Table[D[measure currents[[k, j]], coords[[j]]], {j, d}]]/measure, {k, 3}];
  packed = Total[MapThread[Times, {masks, currents}]];
  concentrated = Expand[Total[Table[D[measure packed[[j]], coords[[j]]], {j, d}]]/measure];
  atoms = Table[Expand[Coefficient[concentrated, DiracDelta[transverse - faces[[k]]]]] /. transverse -> faces[[k]], {k, 2}];
  regular = nativeDiv[[2]];
  alone[prefix <> "/" <> eq <> "/regional currents", currents];
  alone[prefix <> "/" <> eq <> "/regional balance operands", {nativeDiv, data["rhs"][[s]]}];
  alone[prefix <> "/" <> eq <> "/regional balance residuals", Expand[nativeDiv - data["rhs"][[s]]]];
  alone[prefix <> "/" <> eq <> "/distribution divergence", concentrated];
  alone[prefix <> "/" <> eq <> "/face concentrated coefficients", atoms];
  alone[prefix <> "/" <> eq <> "/face concentrated measures", MapThread[Times, {atoms, DiracDelta[transverse - #] & /@ faces}]];
  (* Piola pullback of the interior current, with its geometry density. *)
  working[prefix <> "/" <> eq <> "/fixed interval transport"];
  piola = normalizeArguments[Simplify[jacDet inverse . ((measure currents[[2]]) /. transverse -> Last[map])]];
  pulledNative = Expand[normalizeArguments[jacDet ((measure regular) /. transverse -> Last[map])]];
  pulledDiv = Expand[normalizeArguments[Total[Table[D[piola[[j]], chart[[j]]], {j, d}]]]];
  operands[prefix <> "/" <> eq <> "/coordinate divergence", pulledNative, pulledDiv];
  alone[prefix <> "/" <> eq <> "/pulled current", piola];
  density = Cancel[integratePullback[piola[[1]]]/measure];
  planeFlux = (Cancel[integratePullback[#]/measure] &) /@ piola[[2 ;; -2]];
  (* The evaluated transverse integral supplies BOTH endpoints. *)
  transverseFace = Cancel[boundaryIntegral[Last[piola], z, 0, 1]/measure];
  partial = (Cancel[#/measure] &) /@ boundaryContributions[Last[piola], z, 0, 1];
  alone[prefix <> "/" <> eq <> "/transverse derivative integral operand", formalIntegral[D[Last[piola], z], z, 0, 1]/measure];
  alone[prefix <> "/" <> eq <> "/transverse derivative integral", transverseFace];
  alone[prefix <> "/" <> eq <> "/relative interior transport at each face", partial];
  (* Face concentrated parts are included in full, with no value-based branch. *)
  closedFace = Expand[partial + atoms];
  projectedDerivative = formalD[density, t] + Total[MapThread[(formalD[measure #1, #2]/measure) &, {planeFlux, Rest[parallel]}]];
  left = projectedDerivative + Total[closedFace];
  closedOperand = nativeIntegral[regular] + Total[atoms];
  directIntegrated = nativeIntegral[Expand[regular - D[Last[currents[[2]]], transverse]]] + boundaryIntegral[Last[currents[[2]]], transverse, First[faces], Last[faces]] + Total[atoms];
  operands[prefix <> "/" <> eq <> "/closed integration", directIntegrated, left];
  bulkSource = nativeIntegral[data["rhs"][[s, 2]]];
  required = nativeIntegral[regular] + Total[atoms];
  suppliedSource = bulkSource + If[s == 2, Total[atoms], 0];
  alone[prefix <> "/" <> eq <> "/projected density", density];
  alone[prefix <> "/" <> eq <> "/projected in-plane flux", planeFlux];
  alone[prefix <> "/" <> eq <> "/each closed-projection face term", closedFace];
  alone[prefix <> "/" <> eq <> "/closed integral operand", closedOperand];
  alone[prefix <> "/" <> eq <> "/projected derivative and face operand", left];
  If[s == 1,
   alone[prefix <> "/" <> eq <> "/face balance operands", atoms];
   alone[prefix <> "/" <> eq <> "/projected balance", Inactive[Equal][left, bulkSource]],
   alone[prefix <> "/" <> eq <> "/projected order divergence bulk", nativeIntegral[regular]];
   alone[prefix <> "/" <> eq <> "/projected order source face parts", atoms];
   alone[prefix <> "/" <> eq <> "/projected order divergence", required];
   alone[prefix <> "/" <> eq <> "/projected order source", suppliedSource];
   alone[prefix <> "/" <> eq <> "/source law residual", nativeIntegral[Expand[regular - data["rhs"][[s, 2]]]]];
   alone[prefix <> "/" <> eq <> "/projected balance", Inactive[Equal][left, suppliedSource]]];
  AssociateTo[result, eq -> <|"density" -> density, "plane_flux" -> planeFlux, "faces" -> closedFace, "interior_faces" -> partial, "atoms" -> atoms, "derivative" -> projectedDerivative, "left" -> left, "source" -> If[s == 2, suppliedSource, required], "supplied_source" -> suppliedSource, "bulk" -> regular, "native_current" -> currents[[2]]|>],
 {s, 2}];
 result
];

working["general projection"];
generalInput = suppliedData[False];
general = projection[generalInput, 1, "general projection"];
working["spherical specialization / projection"];
sphericalInput = suppliedData[True];
radialMeasure = Cancel[sphereMeasure/angularConstant];
spherical = projection[sphericalInput, radialMeasure, "spherical specialization"];

(* The radial divergence is integrated with the very same FTC operator;
   angular weights were obtained above from the geometric Jacobian. *)
ClearAll[shellIntegral, atRadius];
shellIntegral[expr_] := formalIntegral[Expand[sphereMeasure expr], r, ri, ro];
atRadius[expr_, radius_] := expr /. r -> radius;
working["spherical specialization / shell channels"];
sigma = spherical["material"]["density"];
radialFlux = First[spherical["material"]["plane_flux"]];
faceChannel = shellIntegral[Total[spherical["material"]["faces"]]];
faceChannelInterior = shellIntegral[Total[spherical["material"]["interior_faces"]]];
faceChannelConcentrated = shellIntegral[Total[spherical["material"]["atoms"]]];
radialBoundaryIntegral = boundaryIntegral[sphereMeasure radialFlux, r, ri, ro];
radialBoundaryContributions = boundaryContributions[sphereMeasure radialFlux, r, ri, ro];
coreChannel = First[radialBoundaryContributions];
outerFlux = Last[radialBoundaryContributions];
drainRate = faceChannel + coreChannel;
orderLost = -shellIntegral[spherical["order"]["source"]];
orderLostFromSuppliedBulk = -shellIntegral[spherical["order"]["supplied_source"]];
storedMaterial = shellIntegral[sigma];
storageOperand = Inactive[D][storedMaterial, t];
storage = formalD[storedMaterial, t];
shellBalance = storage + radialBoundaryIntegral + faceChannel;
faceMinusOrderLost = faceChannel - orderLost;
outerFromBalance = Expand[outerFlux - shellBalance];
outerRadiusDerivative = formalD[outerFlux, ro];
nonzeroDensity = Unequal[atRadius[sigma, ro], 0];
projectedVelocity = ConditionalExpression[atRadius[radialFlux/sigma, ro], nonzeroDensity];
velocityVolumeFlux = ConditionalExpression[atRadius[sphereMeasure radialFlux/sigma, ro], nonzeroDensity];
Do[alone["spherical specialization / shell/" <> item[[1]], item[[2]]],
 {item, {{"domain", Inactive[Less][0, ri, ro]}, {"face channel", faceChannel}, {"face channel interior contribution", faceChannelInterior}, {"face channel concentrated contribution", faceChannelConcentrated}, {"core channel", coreChannel}, {"J", drainRate}, {"order lost", orderLost}, {"order lost using supplied bulk source", orderLostFromSuppliedBulk}, {"outer material flux", outerFlux}, {"stored material", storedMaterial}, {"storage derivative operand", storageOperand}, {"storage", storage}, {"radial derivative integral operand", formalIntegral[D[sphereMeasure radialFlux, r], r, ri, ro]}, {"radial derivative integral", radialBoundaryIntegral}, {"balance operands", {storage, outerFlux, coreChannel, faceChannel}}, {"balance residual", shellBalance}, {"balance", Inactive[Equal][shellBalance, 0]}, {"outer flux balance operands", {outerFlux, outerFromBalance}}, {"outer flux balance", Inactive[Equal][outerFlux, outerFromBalance]}, {"face channel minus order lost", faceMinusOrderLost}, {"outer flux as a function of radius", outerFlux}, {"outer flux radius derivative", outerRadiusDerivative}, {"projected velocity domain", nonzeroDensity}, {"projected velocity", projectedVelocity}, {"projected velocity volume flux", velocityVolumeFlux}}}];

(* Structural dimension propagation; heterogeneous additive dimensions remain
   visible as dimensionTerms rather than being discarded. Zero literals have
   no intrinsic dimension and do not constrain an additive expression. *)
ClearAll[dimension, dimensionSum];
dimensionSum[ds_] := Module[{u = DeleteDuplicates[DeleteCases[ds, dimensionOfZero]]}, Which[Length[u] == 0, dimensionOfZero, Length[u] == 1, First[u], True, dimensionTerms[u]]];
dimension[0] := dimensionOfZero;
dimension[e_?NumberQ] := {0, 0, 0};
dimension[Pi] := {0, 0, 0};
dimension[e_Symbol] := Lookup[coordinateDimensions, e, dimensionUnknown[e]];
dimension[e_Plus] := dimensionSum[dimension /@ (List @@ e)];
dimension[e_Times] := Total[DeleteCases[dimension /@ (List @@ e), dimensionOfZero]];
dimension[Power[b_, p_?NumberQ]] := p dimension[b];
dimension[Inactive[Integrate][body_, {v_, a_, b_}]] := dimension[body] + dimension[v];
dimension[Inactive[D][body_, v_]] := dimension[body] - dimension[v];
dimension[ConditionalExpression[body_, condition_]] := ConditionalExpression[dimension[body], condition];
dimension[Derivative[orders__][head_][args__]] := Lookup[fieldDimensions, head, dimensionUnknown[head]] - Total[MapThread[Times, {{orders}, dimension /@ {args}}]];
dimension[e_] /; KeyExistsQ[fieldDimensions, Head[e]] := fieldDimensions[Head[e]];
dimension[e_List] := dimension /@ e;
dimension[e_] := dimensionUnknown[e];
working["dimensional homogeneity"];
alone["dimensions / field assignment", fieldDimensions];
alone["dimensions / coordinate assignment", coordinateDimensions];
Do[
 Do[
  obj = package[[2]][kind];
  Do[alone[package[[1]] <> "/dimensions/" <> kind <> "/" <> key, dimension[obj[key]]],
    {key, {"density", "plane_flux", "faces", "atoms", "derivative", "left", "source", "supplied_source"}}];
  alone[package[[1]] <> "/dimensions/" <> kind <> "/balance terms", dimension /@ additiveTerms[Expand[obj["left"]]]];
  alone[package[[1]] <> "/dimensions/" <> kind <> "/source terms", dimension /@ additiveTerms[Expand[obj["supplied_source"]]]],
 {kind, {"material", "order"}}],
 {package, {{"general projection", general}, {"spherical specialization", spherical}}}];
Do[alone["spherical specialization / dimensions/" <> item[[1]], dimension[item[[2]]]],
 {item, {{"face channel", faceChannel}, {"core channel", coreChannel}, {"J", drainRate}, {"order lost", orderLost}, {"outer material flux", outerFlux}, {"stored material", storedMaterial}, {"storage derivative operand", storageOperand}, {"storage", storage}, {"balance residual", shellBalance}, {"face channel minus order lost", faceMinusOrderLost}, {"outer flux radius derivative", outerRadiusDerivative}, {"projected velocity", projectedVelocity}, {"projected velocity volume flux", velocityVolumeFlux}}}];
