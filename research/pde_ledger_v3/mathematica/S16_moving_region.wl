SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
(* S16: independently constructed first, before the Python engine.
   Distributional weak differentiation of the shell window, followed by
   integration by parts. There are no imports, reference answers or guards
   on physical values. Integral heads below are linear unevaluated operators. *)
$HistoryLength = 0;
ClearAll["Global`*"];
outfile = FileNameJoin[{DirectoryName[$InputFileName], "out", "S16_moving_region.out"}];
If[!DirectoryQ[DirectoryName[outfile]], CreateDirectory[DirectoryName[outfile]]];
stream = OpenWrite[outfile, PageWidth -> Infinity];
write[tag_, value_] := (WriteString[stream, tag, " = ", ToString[value, InputForm, PageWidth -> Infinity], "\n"]);
working[tag_] := (Print["WORKING_ON = ", tag]; write["WORKING_ON", tag]);
objects = <||>;
emit[tag_, value_] := (write[tag, value]; AssociateTo[objects, tag -> value]; write[tag <> ".comparison", Missing["NoSecondRoute", "nothing to compare it with"]]);

(* The supplied physical equations, live arguments, and control surfaces. *)
xx = {x1, x2, x3};
args = Append[xx, t];
cc = Through[{c1, c2, c3}[t]];
radii = {rin[t], rout[t]};
density = Sigma @@ args;
materialFlux = Through[{F1, F2, F3}[Sequence @@ args]];
momentumDensity = Through[{G1, G2, G3}[Sequence @@ args]];
force = Through[{f1, f2, f3}[Sequence @@ args]];
piHeads = Table[Symbol["Pi" <> ToString[i] <> ToString[j]], {i, 3}, {j, 3}];
momentumFlux = Map[Apply[#, args] &, piHeads, {2}];
faceTerms = {-(n1 @@ args) (aminus @@ args), (n3 @@ args) (aplus @@ args)};
div[v_] := Total[MapThread[D, {v, xx}]];
materialBalance = D[density, t] + div[materialFlux] + Total[faceTerms];
momentumBalance = Table[D[momentumDensity[[i]], t] + div[momentumFlux[[i]]] - force[[i]], {i, 3}];
yy = MapThread[Subtract, {xx, cc}];
distanceSquared = yy.yy;
distance = Sqrt[distanceSquared];
level = {distance - radii[[1]], radii[[2]] - distance};
window = Times @@ (HeavisideTheta /@ level);
direction = {Sin[theta] Cos[phi], Sin[theta] Sin[phi], Cos[theta]};
chart = Table[cc + radii[[k]] direction, {k, 2}];
volumeChart = cc + rho direction;
assumptions = 0 < rin[t] < rout[t] && 0 < theta < Pi && Element[{theta, phi, t}, Reals];
at[k_, expr_] := Cancel[expr /. Power[distanceSquared, p_] :> radii[[k]]^(2 p) /. Thread[xx -> chart[[k]]]];
area = Table[FullSimplify[Sqrt[Cross[D[chart[[k]], theta], D[chart[[k]], phi]].Cross[D[chart[[k]], theta], D[chart[[k]], phi]]], assumptions], {k, 2}];
volumeJacobian = FullSimplify[Det[Transpose[{D[volumeChart, rho], D[volumeChart, theta], D[volumeChart, phi]}]], rho > 0 && 0 < theta < Pi];
normal = Table[at[k, -D[level[[k]], #]] & /@ xx, {k, 2}];
velocity = D[chart, t];
normalSpeed = MapThread[Dot, {normal, velocity}];
(* V[q] integrates q(x,t) in the shell; B[k,q] integrates a chart-pulled
   integrand q(theta,phi,t) with area[k]. Inactive Integrate fixes the measure. *)
write["schema", <|"V" -> "volume integral; x remains a bound Cartesian argument", "B" -> "surface integral; integrand already pulled to chart", "X" -> "dipole/N: displacement from Xc", "physicalPremises" -> "supplied and unfalsifiable", "shapeChoice" -> "concentric spheres, not field symmetry", "momentOrderChoice" -> 2, "forceScope" -> "undivided supplied f", "PiScope" -> "undivided supplied total momentum flux"|>];
write["domain.unnormalized", assumptions];
write["domain.normalized", assumptions && count[t] != 0];
write["integral.V", Inactive[Integrate][q @@ args /. Thread[xx -> volumeChart], {rho, rin[t], rout[t]}, {theta, 0, Pi}, {phi, 0, 2 Pi}] /. Inactive[Integrate][z_, lims__] :> Inactive[Integrate][z volumeJacobian, lims]];
Do[write["integral.B." <> ToString[k], Inactive[Integrate][q[theta, phi, t] area[[k]], {theta, 0, Pi}, {phi, 0, 2 Pi}]], {k, 2}];
emit["supplied.material", materialBalance];
emit["supplied.momentum", momentumBalance];
emit["surfaces.chart", chart]; emit["surfaces.normal", normal];
emit["surfaces.velocity", velocity]; emit["surfaces.normalSpeed", normalSpeed];
emit["surfaces.measure", area]; emit["volume.measure", volumeJacobian];

(* Partial differentiation in independent geometric jet coordinates. Every
   geometric jet present in the operand is included, with no degree bound.
   Restoring the jets retains their live time dependence. These coordinates
   distinguish field evolution from derivatives of the moving reference. *)
geometryHeads = Join[{c1, c2, c3}, {rin, rout}];
geometryChain[z_] := Module[{jets, center, radial, geometry, slots, frozen, restore, rate},
 jets = DeleteDuplicates[Cases[z, Derivative[n_][h_][t] /; MemberQ[geometryHeads, h], {0, Infinity}]];
 center = Join[cc, Select[jets, MemberQ[{c1, c2, c3}, #[[0, 1]]] &]];
 radial = Join[radii, Select[jets, MemberQ[{rin, rout}, #[[0, 1]]] &]];
 geometry = Join[center, radial]; slots = Table[Unique["geometry"], Length[geometry]];
 frozen = z /. Thread[geometry -> slots];
 restore[v_] := v /. Thread[slots -> geometry];
 rate[group_] := Total[Map[Function[g, restore[D[frozen, slots[[First[FirstPosition[geometry, g]]]]]] D[g, t]], group]];
 <|"localEvolution" -> restore[D[frozen, t]], "centerMotion" -> rate[center], "radiusMotion" -> rate[radial]|>];
windowRates = geometryChain[window];

(* Weak differentiation has separate operands for the weighted local law
   and the moving window. Only the former undergoes spatial integration by
   parts. Window motion is differentiated before multiplying the integrand;
   its origin therefore does not depend on the integrand's velocity degree. *)
windowRules = Join[Thread[(HeavisideTheta /@ level) -> {h1, h2}], Thread[(DiracDelta /@ level) -> {d1, d2}]];
windowPolynomial[z_] := Expand[z /. windowRules];
reduceWindow[z_] := Module[{p = windowPolynomial[z]},
 {Coefficient[Coefficient[p, h1], h2] /. {d1 -> 0, d2 -> 0},
  Table[at[k, Coefficient[p, {d1, d2}[[k]]] /. {h1 -> 1, h2 -> 1, d1 -> 0, d2 -> 0}], {k, 2}]}];
weak[w_, den_, flux_, balance_] := Module[{localRate, raw, integrated, jet, coeff, rest, reduced},
 localRate = First[Solve[balance == 0, D[den, t]]];
 raw = Expand[D[w den window, t] /. localRate];
 integrated = Expand[window (D[w den, t] /. localRate)];
 Do[jet = D[flux[[j]], xx[[j]]]; coeff = Coefficient[integrated, jet];
    integrated = Expand[integrated - coeff jet - D[coeff, xx[[j]]] flux[[j]]], {j, 3}];
 rest = reduceWindow[integrated]; reduced = reduceWindow[integrated + w den D[window, t]];
 <|"localRate" -> localRate, "windowDerivative" -> raw, "afterParts" -> Expand[integrated + w den D[window, t]],
   "reduced" -> reduced, "rest" -> rest,
   "centerMotion" -> reduceWindow[w den windowRates["centerMotion"]],
   "radiusMotion" -> reduceWindow[w den windowRates["radiusMotion"]]|>];

terms[z_] := If[Head[Expand[z]] === Plus, List @@ Expand[z], {Expand[z]}];
selectTerms[z_, heads_] := Total[Select[terms[z], !FreeQ[#, Alternatives @@ heads] &]];
rateParts[tag_, w_, den_, flux_, balance_, sources_] := Module[{derivation, body, components, picked, remainder, reference},
 working[tag]; derivation = weak[w, den, flux, balance];
 write[tag <> ".derivation", KeyDrop[derivation, {"reduced", "rest", "centerMotion", "radiusMotion"}]];
 body = derivation["reduced"][[1]]; remainder = body; components = <||>;
 Do[picked = selectTerms[body, {sources[[k]]}]; remainder = Expand[remainder - picked];
    AssociateTo[components, "bulkSource" <> ToString[k] -> V[picked]], {k, Length[sources]}];
 reference = geometryChain[w];
 Do[picked = Expand[den reference[origin]]; remainder = Expand[remainder - picked];
    AssociateTo[components, "bulkReference." <> origin -> V[picked]], {origin, {"centerMotion", "radiusMotion"}}];
 AssociateTo[components, "bulkWeight" -> V[remainder]];
 Do[Do[AssociateTo[components, "surface" <> ToString[k] <> "." <> origin -> B[k, derivation[origin][[2, k]]]],
       {origin, {"rest", "centerMotion", "radiusMotion"}}], {k, 2}];
 KeyValueMap[emit[tag <> "." <> #1, #2] &, components];
 emit[tag <> ".rate", Total[Values[components]]];
 emit[tag <> ".integratedForceCompanion", V /@ force]; components];

working["moments"];
weights = Join[{1}, yy, Flatten[Outer[Times, yy, yy]]];
moments = (V[Expand[# density]] &) /@ weights;
Do[emit["moment." <> ToString[k - 1], moments[[k]]], {k, Length[moments]}];
write["moment.indices", Join[{{}}, List /@ Range[3], Tuples[Range[3], 2]]];
emit["count.definition", count[t] -> First[moments]];
Do[emit["dipole.definition." <> ToString[i], dip[i][t] -> moments[[i + 1]]], {i, 3}];
Do[emit["momentum.definition." <> ToString[i], mom[i][t] -> V[momentumDensity[[i]]]], {i, 3}];
countParts = rateParts["count", 1, density, materialFlux, materialBalance, {n1, n3}];
dipParts = Table[rateParts["dipole." <> ToString[i], yy[[i]], density, materialFlux, materialBalance, {n1, n3}], {i, 3}];
pParts = Table[rateParts["momentum." <> ToString[i], 1, momentumDensity[[i]], momentumFlux[[i]], momentumBalance[[i]], {Head[force[[i]]]}], {i, 3}];
working["momentumReceived"];
received = Table[Total[Lookup[pParts[[i]], {"bulkSource1", "surface1.rest", "surface2.rest"}]], {i, 3}];
emit["momentumReceived", received];
emit["momentumRate", Total[Values[#]] & /@ pParts];
emit["momentumSurfaceMotion", KeyTake[#, {"surface1.centerMotion", "surface1.radiusMotion", "surface2.centerMotion", "surface2.radiusMotion"}] & /@ pParts];

(* Differentiate the unevaluated integrals by their measures. Volume
   transport is derived again from the window derivative with no PDE. *)
volTime[z_] := Module[{red = reduceWindow[D[z window, t]]}, V[red[[1]]] + Sum[B[k, red[[2, k]]], {k, 2}]];
boundTime[k_, z_] := B[k, Expand[D[z area[[k]], t]/area[[k]]]];
time[z_Plus] := Total[time /@ (List @@ z)];
time[z_Times] := Module[{f = List @@ z}, Total[MapIndexed[time[#1] Times @@ Delete[f, First[#2]] &, f]]];
time[z_Power] := z[[2]] z[[1]]^(z[[2]] - 1) time[z[[1]]];
time[V[z_]] := volTime[z];
time[B[k_, z_]] := boundTime[k, z];
time[z_] := D[z, t];
(* Substitute the governed Sigma time jet at arbitrary live arguments;
   rates of F, n, a and f have no supplied evolution and remain derivatives. *)
(* The integrand and the integration domain are differentiated separately.
   A motion jet in z is a coefficient of the window derivative, not another
   attribution for it. This preserves squares and mixed velocity products. *)
timeParts[V[z_]] := Module[{local = geometryChain[z], result, red},
 result = <|"localEvolution" -> V[local["localEvolution"]],
            "bulk.centerMotion" -> V[local["centerMotion"]],
            "bulk.radiusMotion" -> V[local["radiusMotion"]]|>;
 Do[red = reduceWindow[z windowRates[origin]];
    Do[AssociateTo[result, "surface" <> ToString[k] <> "." <> origin -> B[k, red[[2, k]]]], {k, 2}],
    {origin, {"centerMotion", "radiusMotion"}}]; result];
timeParts[B[k_, z_]] := Module[{local = geometryChain[z area[[k]]]},
 <|"localEvolution" -> B[k, Expand[local["localEvolution"]/area[[k]]]],
   "surface" <> ToString[k] <> ".centerMotion" -> B[k, Expand[local["centerMotion"]/area[[k]]]],
   "surface" <> ToString[k] <> ".radiusMotion" -> B[k, Expand[local["radiusMotion"]/area[[k]]]]|>];
(* Linear integral comparison: residuals are integrands on the same printed
   measures. The operands are the direct product derivative and the sum of
   the separately differentiated origins. No expected-value guard is used. *)
integralResidual[left_, right_] := Module[{delta = Expand[left - right], vs, bs, all, residual},
 vs = DeleteDuplicates[Cases[delta, _V, {0, Infinity}]];
 bs = Table[DeleteDuplicates[Cases[delta, B[k, _], {0, Infinity}]], {k, 2}];
 all = Join[vs, Flatten[bs]];
 residual = <|"outsideIntegrals" -> Expand[delta - Total[(Coefficient[delta, #] # &) /@ all]],
               "volumeIntegrand" -> Cancel[Total[(Coefficient[delta, #] #[[1]] &) /@ vs]]|>;
 Do[AssociateTo[residual, "surface" <> ToString[k] <> "Integrand" -> Cancel[Total[(Coefficient[delta, #] #[[2]] &) /@ bs[[k]]]]], {k, 2}]; residual];
partitionReport[tag_, direct_, pieces_Association] := Module[{total = Total[Values[pieces]], residual},
 write[tag <> ".direct", direct]; AssociateTo[objects, tag <> ".direct" -> direct];
 write[tag <> ".partsSum", total]; AssociateTo[objects, tag <> ".partsSum" -> total];
 residual = integralResidual[direct, total];
 write[tag <> ".residual", residual];
 write[tag <> ".routes", {"direct product derivative", "chain-rule origins and window-only motion"}]];
materialRHS = D[density, t] /. First[Solve[materialBalance == 0, D[density, t]]];
governed[z_] := z /. Derivative[0, 0, 0, 1][Sigma][u_, v_, w_, s_] :> (materialRHS /. Thread[args -> {u, v, w, s}]);
working["centerKinematics"];
rateAliases = <||>; derivativeAliases = <||>;
bindParts[p_, prefix_] := Association[KeyValueMap[Function[{key, value}, Module[{name = Symbol[prefix <> StringReplace[key, "." -> ""]][t], dv, pieces},
 emit[prefix <> ".term." <> key, name -> value]; dv = governed[time[value]];
 emit[prefix <> ".termRate." <> key, D[name, t] -> dv];
 pieces = Map[governed, timeParts[value]];
 KeyValueMap[emit[prefix <> ".termRatePart." <> key <> "." <> #1, #2] &, pieces];
 partitionReport[prefix <> ".termRatePartition." <> key, dv, pieces];
 AssociateTo[rateAliases, name -> value]; AssociateTo[derivativeAliases, D[name, t] -> dv]; key -> name]], p]];
namesN = bindParts[countParts, "nrate"];
namesD = Table[bindParts[dipParts[[i]], "drate" <> ToString[i]], {i, 3}];
replacementRates = Join[{Derivative[1][count][t] -> Total[Values[namesN]]}, Table[Derivative[1][dip[i]][t] -> Total[Values[namesD[[i]]]], {i, 3}]];
replacementSeconds = Join[{Derivative[2][count][t] -> Total[D[Values[namesN], t]]}, Table[Derivative[2][dip[i]][t] -> Total[D[Values[namesD[[i]]], t]], {i, 3}]];
relativeCenter = Table[Divide[dip[i][t], count[t]], {i, 3}];
absoluteCenter = MapThread[Plus, {cc, relativeCenter}];
emit["center.relative.definition", relativeCenter];
emit["center.absolute.definition", absoluteCenter];
emit["center.relative.rate", Expand[D[relativeCenter, t] /. replacementRates]];
emit["center.relative.acceleration", Expand[D[relativeCenter, {t, 2}] /. replacementSeconds /. replacementRates]];
emit["center.absolute.rate", Expand[D[absoluteCenter, t] /. replacementRates]];
emit["center.absolute.acceleration", Expand[D[absoluteCenter, {t, 2}] /. replacementSeconds /. replacementRates]];
emit["count.rateForNormalization", Total[Values[namesN]]];
emit["count.accelerationForNormalization", Total[D[Values[namesN], t]]];

(* Dimension vectors use independent length, time, count, momentum units.
   Momentum is an independent unit: no mass or count-to-momentum law.
   Field dimensions are solved from the supplied balance terms; face-factor
   dimensions are velocity dimensions. No numerical value is tested. *)
working["dimensions"];
baseDim = <|x1 -> {1, 0, 0, 0}, x2 -> {1, 0, 0, 0}, x3 -> {1, 0, 0, 0}, t -> {0, 1, 0, 0}, rho -> {1, 0, 0, 0}, theta -> {0, 0, 0, 0}, phi -> {0, 0, 0, 0}|>;
fieldDim = <||>;
Do[AssociateTo[fieldDim, h -> {1, 0, 0, 0}], {h, {c1, c2, c3, rin, rout}}];
AssociateTo[fieldDim, Sigma -> {-3, 0, 1, 0}];
Do[AssociateTo[fieldDim, h -> {-3, 0, 0, 1}], {h, {G1, G2, G3}}];
Do[AssociateTo[fieldDim, h -> {1, -1, 0, 0}], {h, {aminus, aplus}}];
Do[AssociateTo[fieldDim, h -> (fieldDim[Sigma] - baseDim[t] + baseDim[x1])], {h, {F1, F2, F3}}];
Do[AssociateTo[fieldDim, h -> (fieldDim[G1] - baseDim[t] + baseDim[x1])], {h, Flatten[piHeads]}];
Do[AssociateTo[fieldDim, h -> (fieldDim[G1] - baseDim[t])], {h, {f1, f2, f3}}];
Do[AssociateTo[fieldDim, h -> (fieldDim[Sigma] - baseDim[t] - fieldDim[aminus])], {h, {n1, n3}}];
AssociateTo[fieldDim, count -> (fieldDim[Sigma] + 3 baseDim[x1])];
Do[AssociateTo[fieldDim, dip[i] -> (fieldDim[count] + baseDim[x1])]; AssociateTo[fieldDim, mom[i] -> (fieldDim[G1] + 3 baseDim[x1])], {i, 3}];
write["dimensions.basis", {"length", "time", "count", "momentum"}]; write["dimensions.fields", fieldDim];
dimension[0] := {};
dimension[z_?NumberQ] := {{0, 0, 0, 0}};
dimension[Pi] := {{0, 0, 0, 0}};
dimension[z_Symbol] := {Lookup[baseDim, z, Missing["dimension", z]]};
dimension[z_Plus] := DeleteDuplicates[Join @@ (dimension /@ (List @@ z))];
dimension[z_Times] := Fold[DeleteDuplicates[Flatten[Outer[Plus, #1, #2, 1], 1]] &, {{0, 0, 0, 0}}, (dimension /@ (List @@ z))];
dimension[z_Power] := (# z[[2]] &) /@ dimension[z[[1]]];
dimension[Sin[_] | Cos[_] | HeavisideTheta[_]] := {{0, 0, 0, 0}};
dimension[V[z_]] := (# + 3 baseDim[x1] &) /@ dimension[z];
dimension[B[_, z_]] := (# + 2 baseDim[x1] &) /@ dimension[z];
dimension[Derivative[orders__][h_][ar___]] := Module[{dd = Lookup[fieldDim, h, Missing["dimension", h]], av = {ar}},
 {dd - Total[MapThread[Times, {{orders}, If[Length[av] == 1, {baseDim[t]}, Append[ConstantArray[baseDim[x1], 3], baseDim[t]]]}]]}];
dimension[z_Rule] := {dimension[First[z]], dimension[Last[z]]};
dimension[z_List] := dimension /@ z;
dimension[z_Association] := Map[dimension, z];
dimension[z_] := {Lookup[fieldDim, Head[z], Missing["dimension", HoldForm[z]]]};
KeyValueMap[(AssociateTo[fieldDim, Head[#1] -> (If[StringStartsQ[SymbolName[Head[#1]], "nrate"], fieldDim[count], fieldDim[dip[1]]] - baseDim[t])]) &, rateAliases];
(* For every printed object emit the inferred dimensional operands. For
   each sum in its tree also emit dimensional differences, unconditionally;
   these are dimensional homogeneity data, not expected-value guards. *)
KeyValueMap[Function[{tag, value}, Module[{sums, ds},
 write[tag <> ".dimensions", dimension[value]];
 sums = Cases[value, z_Plus :> z, {0, Infinity}];
 Do[ds = dimension /@ (List @@ sums[[k]]);
 write[tag <> ".dimensionSum." <> ToString[k] <> ".operands", ds];
 write[tag <> ".dimensionSum." <> ToString[k] <> ".residuals", If[MemberQ[ds, {}, Infinity], Missing["zeroTerm"], Map[# - First[ds] &, ds]]], {k, Length[sums]}]]], objects];
write["coverage.objects", Keys[objects]];
write["coverage.ungovernedRates", {F1, F2, F3, n1, n3, aminus, aplus}];
write["coverage.fluxDecomposition", Missing["NoSuppliedLaw", "Pi has no separately supplied stress or advection constituents"]];
Close[stream];
Print["TRANSCRIPT = ", outfile];
Exit[0];
