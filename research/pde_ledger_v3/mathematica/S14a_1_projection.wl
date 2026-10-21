SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
(* Independent weak/product-rule construction from the supplied equations and U1.
   No Python source, output, exports, or previously derived balance is read.
   Formal integrals retain arbitrary live fields. No radius limit is taken. *)
$HistoryLength = 0;
transcript = OpenWrite[FileNameJoin[{DirectoryName[$InputFileName], "out", "S14a_1_projection.out"}], PageWidth -> Infinity];
emit[tag_, obj_] := Module[{line = tag <> ": " <> ToString[obj, InputForm, PageWidth -> Infinity]},
  Print[line]; WriteString[transcript, line <> "\n"]];
working[tag_] := emit["WORKING_ON", tag];
object[tag_, obj_] := (emit[tag, obj]; emit[tag <> "/comparison", "No second route; nothing to compare it with."]);
balance[tag_, a_, b_] := (emit[tag <> "/operands", {a, b}];
  emit[tag <> "/residual", Expand[a - b]];
  emit[tag <> "/comparison", "Balance residual, not an independent-route comparison."]);
buildStatus = Check[
emit["contract", "S14a_1_decision_list.md revision 3; supplied balances, U1 and face factor are premises"];
emit["kernel", {$Version, SystemOptions["ParallelOptions"], SystemOptions["InequalitySolvingOptions" -> "ThreadOr"]}];

(* Dimensions are input metadata: exponent order {material number, length, time}. *)
unit = <||>;
numberUnit = {1, 0, 0}; lengthUnit = {0, 1, 0}; timeUnit = {0, 0, 1};
zeroUnit = {0, 0, 0};
register[s_Symbol, v_] := (AssociateTo[unit, SymbolName[s] -> v]; s);
field[name_, args_, v_] := (register[Symbol[name], v]; Symbol[name] @@ args);
Scan[register[#, lengthUnit] &, {x, y, z, w, radius, innerRadius, outerRadius}];
register[time, timeUnit]; register[theta, zeroUnit]; register[phi, zeroUnit];

(* Generic calculus, not a library of physical answers. The linear integral and
   its Leibniz endpoint pieces are kept explicit, including both closed endpoints. *)
linearIntegral[0, lim_] := 0;
linearIntegral[e_, lim_] := If[Head[Expand[e]] === Plus,
  Total[Inactive[Integrate][#, lim] & /@ (List @@ Expand[e])], Inactive[Integrate][Expand[e], lim]];
trace[e_, v_, at_] := e /. v -> at;
leibnizPieces[e_, lim_, q_] := {
  linearIntegral[D[e, q], lim],
  -trace[e, lim[[1]], lim[[2]]] D[lim[[2]], q],
  trace[e, lim[[1]], lim[[3]]] D[lim[[3]], q]};
continuity[d_, f_, vars_, clock_, measure_] := D[d, clock] +
  Total[MapThread[D[measure #1, #2] &, {f, vars}]]/measure;
stepDerivative[e_, v_] := D[e, v];

(* The only native physical products are the supplied equations, below. *)
supplied[d_, v_, c_, extra_, g_, gr_, gd_, vars_, clock_, measure_] := Module[{md, mf, od, of, os},
  md = d; mf = d v; od = c d; of = c d v + extra; os = d g;
  <|"materialDensity" -> md, "materialFlux" -> mf,
    "orderDensity" -> od, "orderFlux" -> of, "orderSource" -> os,
    "materialEquation" -> {continuity[md, mf, vars, clock, measure], 0},
    "orderEquation" -> {continuity[od, of, vars, clock, measure], os},
    "gammaEquation" -> {g, gr - gd}|>];

construct[label_, vars_, clock_, args_, faceArgs_, measure_] := Module[
  {lo, hi, limits, masks, data, names, n, vel, extra, c, g, gr, gd,
   spacetime, density, flux, regular, atomPieces, atoms, source, projected,
   result, pieces, lhs, formal, faceTerms, rhs, local, kind, s, i, j, normal,
   suppliedFactor, contractions, outward, faceTimes},
  working[label <> "/native fields and U1"];
  lo = field[label <> "LowerFace", faceArgs, lengthUnit];
  hi = field[label <> "UpperFace", faceArgs, lengthUnit];
  limits = {Last[vars], lo, hi};
  (* U1 support encoding; zero here is the stipulated exterior order, not a freeze. *)
  masks = {1 - UnitStep[Last[vars] - lo],
    UnitStep[Last[vars] - lo] - UnitStep[Last[vars] - hi], UnitStep[Last[vars] - hi]};
  names = {"Below", "Inside", "Above"};
  data = Table[
    n = field[label <> "n" <> names[[i]], args, numberUnit - 4 lengthUnit];
    vel = Table[field[label <> "u" <> ToString[j] <> names[[i]], args, lengthUnit - timeUnit], {j, Length[vars]}];
    extra = Table[field[label <> "Jchi" <> ToString[j] <> names[[i]], args, numberUnit - 3 lengthUnit - timeUnit], {j, Length[vars]}];
    c = If[i == 2, field[label <> "chiInside", args, zeroUnit], 0];
    g = field[label <> "GammaB" <> names[[i]], args, -timeUnit];
    gr = field[label <> "GammaReturn" <> names[[i]], args, -timeUnit];
    gd = field[label <> "GammaDrain" <> names[[i]], args, -timeUnit];
    supplied[n, vel, c, extra, g, gr, gd, vars, clock, measure], {i, 3}];
  object[label <> "/faces/domain", lo < hi];
  object[label <> "/U1/masks", masks];
  Do[KeyValueMap[(object[label <> "/supplied/" <> names[[i]] <> "/" <> #1, #2]) &, data[[i]]], {i, 3}];
  spacetime = Prepend[vars, clock];
  result = <|"faces" -> {lo, hi}, "limits" -> limits, "data" -> data|>;
  Do[
    working[label <> "/" <> kind <> "/distribution and projection"];
    density = Lookup[data, kind <> "Density"];
    flux = Lookup[data, kind <> "Flux"];
    regular = MapThread[continuity[#1, #2, vars, clock, measure] &, {density, flux}];
    (* Weak product rule: differentiate the region indicators, then contract
       their spacetime covectors with each region's conserved current. *)
    atomPieces = Table[
      local = Total[MapThread[Times, {D[masks[[i]], #] & /@ spacetime, Prepend[flux[[i]], density[[i]]]}]];
      Table[trace[Coefficient[Expand[local], DiracDelta[Last[vars] - {lo, hi}[[s]]]], Last[vars], {lo, hi}[[s]]], {s, 2}], {i, 3}];
    atoms = Total[atomPieces];
    source = If[kind === "order", Lookup[data, "orderSource"], ConstantArray[0, 3]];
    object[label <> "/" <> kind <> "/regional differential operands", regular];
    object[label <> "/" <> kind <> "/face contributions by region", atomPieces];
    object[label <> "/" <> kind <> "/face concentrated parts", atoms];
    object[label <> "/" <> kind <> "/distribution", Total[masks regular] +
      Total[MapThread[Times, {DiracDelta[Last[vars] - #] & /@ {lo, hi}, atoms}]]];
    projected = {linearIntegral[density[[2]], limits], linearIntegral[#, limits] & /@ Most[flux[[2]]]};
    object[label <> "/" <> kind <> "/projected density", projected[[1]]];
    object[label <> "/" <> kind <> "/projected in-plane flux", projected[[2]]];
    pieces = leibnizPieces[density[[2]], limits, clock] +
      Total[MapThread[leibnizPieces[measure #1, limits, #2]/measure &, {Most[flux[[2]]], Most[vars]}]];
    lhs = Total[pieces];
    formal = Inactive[D][projected[[1]], clock] +
      Total[MapThread[Inactive[D][measure #1, #2]/measure &, {projected[[2]], Most[vars]}]];
    (* The integrated regular operator has the same first piece; all face
       terms below are differences of its Leibniz endpoint pieces and atoms. *)
    faceTerms = Expand[atoms - Rest[pieces]];
    rhs = If[kind === "order", linearIntegral[source[[2]], limits] + Total[atoms], 0];
    object[label <> "/" <> kind <> "/projected differential notation", formal];
    object[label <> "/" <> kind <> "/projected differential expansion", lhs];
    Do[object[label <> "/" <> kind <> "/face " <> ToString[s] <> "/leaving term", faceTerms[[s]]], {s, 2}];
    object[label <> "/" <> kind <> "/projected source", rhs];
    balance[label <> "/" <> kind <> "/closed projected balance", lhs + Total[faceTerms], rhs];
    AssociateTo[result, kind -> <|"density" -> projected[[1]], "flux" -> projected[[2]],
      "left" -> lhs, "formal" -> formal, "faces" -> faceTerms, "source" -> rhs,
      "atoms" -> atoms, "regular" -> regular, "projectedBalance" -> {lhs + Total[faceTerms], rhs}|>],
    {kind, {"material", "order"}}];
  (* Supplied face factor, printed with its operands; not an independent check. *)
  Do[
    normal = Join[{-D[{lo, hi}[[s]], clock]}, -D[{lo, hi}[[s]], #] & /@ Most[vars], {1}];
    object[label <> "/face " <> ToString[s] <> "/spacetime conormal", normal];
    Do[
      n = data[[i]]["materialDensity"];
      (* Extract the native velocity from the supplied material flux. *)
      vel = Cancel[#/n] & /@ data[[i]]["materialFlux"];
      suppliedFactor = Last[vel] - D[{lo, hi}[[s]], clock] - Total[Most[vel] (D[{lo, hi}[[s]], #] & /@ Most[vars])];
      object[label <> "/face " <> ToString[s] <> "/" <> names[[i]] <> "/supplied material factor", trace[suppliedFactor, Last[vars], {lo, hi}[[s]]]],
      {i, 3}], {s, 2}];
  result];

cartesian = construct["General", {x, y, z, w}, time, {x, y, z, w, time}, {x, y, z, time}, 1];

working["Spherical specialization/geometry"];
embedding = {radius Sin[theta] Cos[phi], radius Sin[theta] Sin[phi], radius Cos[theta]};
jacobian = Transpose[D[embedding, #] & /@ {radius, theta, phi}];
metric = Simplify[Transpose[jacobian].jacobian];
measure = FullSimplify[Sqrt[Det[metric]], Assumptions -> radius > 0 && 0 < theta < Pi];
area = Integrate[measure, {theta, 0, Pi}, {phi, 0, 2 Pi}, Assumptions -> radius > 0];
object["Spherical specialization/coordinate embedding", embedding];
object["Spherical specialization/metric", metric];
object["Spherical specialization/volume Jacobian", measure];
object["Spherical specialization/sphere area", area];
(* Radial and bulk components are native live functions of {radius,w,time}.
   Angular independence and the radial vector representation are exactly the
   labelled spherical specialization, not frozen Cartesian component fields. *)
spherical = construct["SphericalSpecialization", {radius, w}, time, {radius, w, time}, {radius, time}, area];
emit["SphericalSpecialization/representation", {"u1 and Jchi1 denote radial components; u2 and Jchi2 denote bulk components", "every native field has radius,w,time arguments"}];

working["SphericalSpecialization/shell and outer radius objects"];
shell[e_] := linearIntegral[Expand[area e], {radius, innerRadius, outerRadius}];
sphere[e_, at_] := trace[area e, radius, at];
m = spherical["material"]; o = spherical["order"];
faceChannelParts = shell /@ m["faces"];
faceChannel = Total[faceChannelParts];
coreChannel = -sphere[First[m["flux"]], innerRadius];
drainRate = Total[{faceChannel, coreChannel}];
orderLost = -shell[o["source"]];
outerFlux = sphere[First[m["flux"]], outerRadius];
(* Apply the time derivative to the native density with both live endpoints. *)
storage = shell[Total[leibnizPieces[spherical["data"][[2]]["materialDensity"], spherical["limits"], time]]];
storageNotation = Inactive[D][shell[m["density"]], time];
shellBalanceLeft = Total[{storage, outerFlux, drainRate}];
shellBalanceRight = shell[m["projectedBalance"][[2]]];
faceMinusOrder = Expand[faceChannel - orderLost];
radiusDerivative = Total[leibnizPieces[area spherical["data"][[2]]["materialFlux"][[1]], spherical["limits"], radius]];
radiusDerivative = trace[radiusDerivative, radius, outerRadius];
object["SphericalSpecialization/shell domain", 0 < innerRadius < outerRadius];
object["SphericalSpecialization/face channel by face", faceChannelParts];
object["SphericalSpecialization/face channel", faceChannel];
object["SphericalSpecialization/core channel", coreChannel];
object["SphericalSpecialization/J", drainRate];
object["SphericalSpecialization/order lost", orderLost];
object["SphericalSpecialization/outward outer flux", outerFlux];
object["SphericalSpecialization/storage notation", storageNotation];
object["SphericalSpecialization/storage", storage];
balance["SphericalSpecialization/shell balance", shellBalanceLeft, shellBalanceRight];
balance["SphericalSpecialization/face channel minus order lost", faceChannel, orderLost];
object["SphericalSpecialization/outer flux radius derivative", radiusDerivative];
velocityCondition = m["density"] != 0;
projectedVelocity = First[m["flux"]]/m["density"];
volumeFlux = sphere[projectedVelocity, outerRadius];
object["SphericalSpecialization/projected velocity", ConditionalExpression[projectedVelocity, velocityCondition]];
object["SphericalSpecialization/outer projected velocity", ConditionalExpression[trace[projectedVelocity, radius, outerRadius], trace[velocityCondition, radius, outerRadius]]];
object["SphericalSpecialization/outer volume flux", ConditionalExpression[volumeFlux, trace[velocityCondition, radius, outerRadius]]];

(* Expression-tree dimensional analysis: no target dimensions are supplied for
   the projected objects. Every additive comparison emits operands first. *)
working["dimensional objects"];
dimension[0] := Missing["AdditiveIdentityDimension"];
dimension[e_?NumberQ] := zeroUnit;
dimension[Pi] := zeroUnit;
dimension[e_Symbol] := Lookup[unit, SymbolName[e], Missing["Unregistered", SymbolName[e]]];
dimension[e_Plus] := dimension[First[e]];
dimension[e_Times] := Total[dimension /@ (List @@ e)];
dimension[Power[b_, p_]] := p dimension[b];
dimension[Derivative[orders__][f_][args__]] := Lookup[unit, SymbolName[f]] - Total[MapThread[Times, {{orders}, dimension /@ {args}}]];
dimension[Inactive[Integrate][f_, {v_, a_, b_}]] := dimension[f] + dimension[v];
dimension[Inactive[D][f_, v_]] := dimension[f] - dimension[v];
dimension[e_] := Lookup[unit, SymbolName[Head[e]], Missing["UnregisteredHead", Head[e]]];
dimensionAudit[tag_, e_, dom_: True] := Module[{adds, ds},
  emit[tag <> "/dimension operand", ConditionalExpression[e, dom]]; emit[tag <> "/dimension", dimension[e]];
  adds = Cases[{e}, q_Plus :> List @@ q, Infinity];
  Do[ds = dimension /@ adds[[i]];
    emit[tag <> "/sum " <> ToString[i] <> "/operands", Transpose[{adds[[i]], ds}]];
    emit[tag <> "/sum " <> ToString[i] <> "/dimension residuals", (# - First[ds] &) /@ ds], {i, Length[adds]}];
  emit[tag <> "/comparison", "Only dimensions of additive operands are compared; no second construction of this object."]];
Scan[Function[pair, With[{lab = pair[[1]], obj = pair[[2]]}, Do[
  dimensionAudit[lab <> "/dimensions/" <> kind <> "/density", obj[kind]["density"]];
  MapIndexed[dimensionAudit[lab <> "/dimensions/" <> kind <> "/flux " <> ToString[First[#2]], #1] &, obj[kind]["flux"]];
  MapIndexed[dimensionAudit[lab <> "/dimensions/" <> kind <> "/face " <> ToString[First[#2]], #1] &, obj[kind]["faces"]];
  dimensionAudit[lab <> "/dimensions/" <> kind <> "/balance left", obj[kind]["projectedBalance"][[1]]];
  dimensionAudit[lab <> "/dimensions/" <> kind <> "/balance right", obj[kind]["projectedBalance"][[2]]];
  dimensionAudit[lab <> "/dimensions/" <> kind <> "/balance residual", Subtract @@ obj[kind]["projectedBalance"]],
  {kind, {"material", "order"}}]]], {{"General", cartesian}, {"SphericalSpecialization", spherical}}];
Scan[dimensionAudit["SphericalSpecialization/dimensions/" <> #[[1]], #[[2]]] &,
  {{"face channel", faceChannel}, {"core channel", coreChannel}, {"J", drainRate},
   {"order lost", orderLost}, {"outward outer flux", outerFlux}, {"storage", storage},
   {"shell balance", shellBalanceLeft - shellBalanceRight}, {"face channel minus order lost", faceMinusOrder},
   {"outer flux radius derivative", radiusDerivative}}];
emit["SphericalSpecialization/dimensions/projected velocity/domain", velocityCondition];
dimensionAudit["SphericalSpecialization/dimensions/projected velocity", projectedVelocity, velocityCondition];
emit["SphericalSpecialization/dimensions/volume flux/domain", trace[velocityCondition, radius, outerRadius]];
dimensionAudit["SphericalSpecialization/dimensions/volume flux", volumeFlux, trace[velocityCondition, radius, outerRadius]];
True, False];
Close[transcript];
Quit[If[TrueQ[buildStatus], 0, 1]];
