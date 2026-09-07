(* S11c-c2 N6. Author: Codex gpt-6-astra, high.
   Self-contained and input-free. Construction authorities:
   directives/S11b_SHARED_PHYSICS.md; directives/S11c_a_SHARED_PHYSICS.md;
   directives/S11c_b_SHARED_PHYSICS.md; directives/S11c_c1_SHARED_PHYSICS.md;
   directives/S11c_c2_SHARED_PHYSICS.md;
   _measurements/S11c_c2_N6_route2_spec_astra.md; the blind N6 build directive.
   This instrument does not implement the full c2 self-energy fold.
   All quantities are computed before their measurements are classified externally.
*)
$HistoryLength = 0;
ClearAll["Global`*"];
$Messages = {OutputStream["stderr", 2]};
$MaxExtraPrecision = 0;
$IterationLimit = 10000000;

(* Independent builder inputs. Build legs change these in a temporary copy. *)
materialNormalKnife = 0;
actualAdvectionCoefficient = 1;
actualThicknessCoefficient = 1;
actualJunkCoefficient = 0;
actualJunkCase = {"MATERIAL_ADVECTED", "RHO4_CONSTANT"};
predictionAdvectionCoefficient = 1;
predictionThicknessCoefficient = 1;

rcNames = {"CARRIER_EULERIAN", "CARRIER_MATERIAL", "CARRIER_BRIDGE_RESIDUAL",
  "SOURCE_EULERIAN", "SOURCE_MATERIAL", "SOURCE_BRIDGE_RESIDUAL",
  "CARRIER_CHANNEL", "SOURCE_CHANNEL", "CROSS_CHANNEL", "EULERIAN_OPERAND",
  "MATERIAL_OPERAND", "R_N6", "SPLIT_SUM", "SPLIT_CHECK", "FROZEN_RELATIONS",
  "PROVENANCE", "ADVECTION_ABSENCE", "DIMENSIONS"};
covNames = {"SOURCE_ACTUAL", "SOURCE_PREDICTED", "R_COV", "SOURCE_BASELINE",
  "R_COV_BASELINE", "SOURCE_CONTROL_DELTA", "R_COV_CONTROL_DELTA",
  "R_COV_INCREMENT", "FROZEN_PHI", "PHI_DOMAIN_CENSUS",
  "ACTUAL_CONTROL_PARAMETERS", "PROVENANCE"};
guardNames = {"SLOT_GUARD_NATIVE", "SLOT_GUARD_CARRIER", "SLOT_GUARD_RESIDUAL",
  "CLOSURE_GUARD_NATIVE", "CLOSURE_GUARD_CARRIER", "CLOSURE_GUARD_RESIDUAL"};
standardEmissionName[{family_String, name_String}] := Switch[family,
  "RC", "WL_S11CC2_N6RC_" <> name, "COV", "WL_S11CC2_N6COV_" <> name,
  "GUARD", "WL_S11CC2_N6_" <> name, _, "WL_S11CC2_N6_LOCAL_" <> name];
emittedNames = {};
stripConditional[x_] := x /. ConditionalExpression[v_, c_] :>
  <|"CONDITIONAL_VALUE" -> v, "CONDITION_OPERAND" -> HoldForm[c]|>;
relationalObject[l_, r_] := Inactive[Equal][l, r];
(* WriteString dispatches the stream writer/flush callback. This explicit final
   stream write implements the named Flush wrapper without a package dependency;
   redirected visibility is also measured by the external run monitor. *)
Global`Flush[stream_OutputStream] := (WriteString[stream, ""]; Null);
flushOutput[] := Flush[OutputStream["stdout", 1]];
emit[key_, payload_] := Module[{name = standardEmissionName[key]},
  If[!StringQ[name], Quit[90]];
  If[MemberQ[emittedNames, name], Quit[91]];
  AppendTo[emittedNames, name];
  WriteString[$Output, name, " = ",
    ToString[stripConditional[payload], InputForm, PageWidth -> Infinity], "\n"];
  flushOutput[]];
beginAssociationEmission[key_] := Module[{name = standardEmissionName[key]},
  If[!StringQ[name], Quit[90]];
  If[MemberQ[emittedNames, name], Quit[91]];
  AppendTo[emittedNames, name]; streamFirst = True;
  WriteString[$Output, name, " = <|"]; flushOutput[]];
appendAssociationEmission[key_, value_] := (
  If[!streamFirst, WriteString[$Output, ", "]]; streamFirst = False;
  WriteString[$Output, ToString[key, InputForm, PageWidth -> Infinity], " -> "];
  writeCAS[stripConditional[value]]; flushOutput[]);
writeCAS[value_Association] := Module[{first = True},
  WriteString[$Output, "<|"];
  KeyValueMap[Function[{key, val}, If[!first, WriteString[$Output, ", "]]; first = False;
    WriteString[$Output, ToString[key, InputForm, PageWidth -> Infinity], " -> "]; writeCAS[val]; flushOutput[]], value];
  WriteString[$Output, "|>"]];
writeCAS[value_] := WriteString[$Output, ToString[value, InputForm, PageWidth -> Infinity]];
(* A spool is an in-memory list of rules; this engine has no file inputs. *)
appendAssociationEmissionFromSpool[rules_List] := Scan[
  appendAssociationEmission[First[#], Last[#]] &, rules];
endAssociationEmission[] := (WriteString[$Output, "|>\n"]; flushOutput[]);

(* The differential algebra registers atoms when they are constructed. A derivative
   appends a sorted multi-index; no finite hand-written prolongation table is used. *)
jetRegistry = <||>; unitRegistry = <||>;
backgroundProfileAtoms = {w1Profile, m1Profile};
declare[s_Symbol, dim_List] := (AssociateTo[unitRegistry, s -> dim]; s);
fieldDimensions = <|"theta" -> {0, 0, 0}, "eW" -> {0, 0, 0},
  "u1" -> {1, 0, 0}, "u2" -> {1, 0, 0}, "u3" -> {1, 0, 0},
  "zetaC" -> {1, 0, 0}, "WBg" -> {1, 0, 0}, "muRBg" -> {-1, -2, 1}|>;
waveFields = {"theta", "eW", "u1", "u2", "u3", "zetaC"};
jet[field_String, indices_List : {}] := Module[{ii = Sort[indices], atom, dim},
  atom = Symbol[field <> If[ii === {}, "", "Jet" <> StringJoin[ToString /@ ii]]];
  If[!KeyExistsQ[jetRegistry, atom],
    AssociateTo[jetRegistry, atom -> {field, ii}];
    dim = Lookup[fieldDimensions, field, Missing["FieldDimension", field]];
    If[ListQ[dim], declare[atom, dim - {Count[ii, 1 | 2 | 3], Count[ii, 4], 0}]]]; atom];
atomsIn[f_] := Select[DeleteDuplicates[Cases[f, _Symbol, {0, Infinity}]],
  KeyExistsQ[jetRegistry, #] &];
td[f_, direction_Integer] := td[f, direction] = Total[Function[a,
  D[f, a] jet[jetRegistry[a][[1]], Append[jetRegistry[a][[2]], direction]]] /@ atomsIn[f]];
tdMany[f_, indices_List] := Fold[td, f, indices];
theta = jet["theta"]; eW = jet["eW"];
displacement = Table[jet["u" <> ToString[i]], {i, 3}];
zetaC = jet["zetaC"]; WBg = jet["WBg"]; muRBg = jet["muRBg"];
Scan[declare[#, {0, 0, 0}] &, {epsilonShape, etaBg, sigmaW, w1Profile, m1Profile}];
Scan[declare[#, {1, 0, 0}] &, {W0, LW}];
declare[rhoBr, {-3, 0, 1}]; declare[rhoM, {-4, 0, 1}];
declare[rhoBrBgRho4Constant, {-3, 0, 1}];
declare[cS0, {1, -1, 0}]; declare[omega, {0, -1, 0}];
Scan[declare[#, {0, 1, 0}] &, {tauA, tauV, tauX}];
declare[LambdaA0, {-5, 1, 1}]; declare[LambdaV0, {-4, 0, 1}];
declare[LambdaX0, {-4, 0, 1}]; declare[muTheta, {-1, -2, 1}];
declare[rhoFace, {-3, 0, 1}]; declare[pFace, {-2, -2, 1}]; declare[vFace, {1, -1, 0}];
declare[jFace, {-3, -1, 1}]; declare[vBulkFace, {1, -1, 0}]; declare[faceNormalOperand, {0, 0, 0}];
declare[junkMu, {-1, -2, 1}]; declare[muR, {-1, -2, 1}];
declare[bRho, {-2, -2, 1}]; declare[cCoupling, {-2, -2, 1}];
declare[kW, {-3, -2, 1}]; declare[kappaW, {-3, -2, 1}];
lambdaA = LambdaA0/(1 - I omega tauA);
lambdaV = LambdaV0/(1 - I omega tauV);
lambdaX = LambdaX0/(1 - I omega tauX);

(* The inherited four pressure identities are instantiated once and passed by
   value to both builders. No name lookup or route-local pressure construction. *)
pressureSlots = {deltaPPlus, dWDeltaPPlus, deltaPMinus, dWDeltaPMinus};
MapThread[declare, {pressureSlots, {{-2, -2, 1}, {-3, -2, 1},
  {-2, -2, 1}, {-3, -2, 1}}}];
pressureAssumptions = Element[#, Complexes] & /@ pressureSlots;
slotZero = Thread[pressureSlots -> 0];
faces = {1, -1};
gradeIndices = Tuples[{Range[0, 1], Range[0, 1]}];
retainTerm[f_] := retainTerm[f] = Expand[Normal[Series[Normal[Series[f, {etaBg, 0, 1}]], {sigmaW, 0, 1}]]];
retain[f_] := Module[{expanded = Expand[f]}, Total[retainTerm /@
  If[Head[expanded] === Plus, List @@ expanded, {expanded}]]];
profileRules[] := Map[Function[a, With[{f = jetRegistry[a][[1]], ii = jetRegistry[a][[2]]},
  a -> If[ii === {}, If[f === "WBg", W0 (1 + etaBg w1Profile), muR (1 + etaBg m1Profile)],
    With[{v = Symbol[If[f === "WBg", "w1Profile", "m1Profile"] <> "Jet" <> StringJoin[ToString /@ ii]]},
      declare[v, {0, 0, 0}]; backgroundProfileAtoms = Union[backgroundProfileAtoms, {v}]; If[MemberQ[ii, 4], 0,
        sigmaW If[f === "WBg", 1, muR/W0] v/LW^(Length[ii] - 1)]]]]],
  Select[Keys[jetRegistry], MemberQ[{"WBg", "muRBg"}, jetRegistry[#][[1]]] &]];
harmonicRules[f_] := Map[Function[a, With[{field = jetRegistry[a][[1]], indices = jetRegistry[a][[2]]},
  a -> (-I omega)^Count[indices, 4] jet[field, DeleteCases[indices, 4]]]],
  Select[atomsIn[f], MemberQ[waveFields, jetRegistry[#][[1]]] && MemberQ[jetRegistry[#][[2]], 4] &]];
finish[f_] := finish[f] = retain[(f /. harmonicRules[f]) /. profileRules[]];
atPoint[f_, location_String] := Module[{aa, bb},
  aa = Select[backgroundProfileAtoms, !FreeQ[f, #] &];
  bb = Map[Function[a, declare[Symbol[SymbolName[a] <> "At" <> location], unitRegistry[a]]], aa];
  f /. Thread[aa -> bb]];
gradePart[f_, g_List] := Coefficient[Coefficient[f, etaBg, g[[1]]], sigmaW, g[[2]]];
waveScale[f_, degree_Integer] := Module[{aa = Select[atomsIn[f], MemberQ[waveFields, jetRegistry[#][[1]]] &]},
  Coefficient[Expand[f /. Thread[aa -> (waveMarker aa)]], waveMarker, degree]];
el[density_] := D[density, theta] - Sum[td[D[density, jet["theta", {i}]], i], {i, 3}];
(* Spatial differentiation never decreases the number of background-jet factors.
   Pruning products beyond sigma^1 is therefore safe before further derivatives;
   zero-jet WBg factors remain live and are NOT expanded/truncated here. *)
pruneBackgroundProducts[f_] := pruneBackgroundProducts[f] = Module[{aa, expanded, terms},
  aa = Select[atomsIn[f], MemberQ[{"WBg", "muRBg"}, jetRegistry[#][[1]]] && jetRegistry[#][[2]] =!= {} &];
  If[aa === {} || !PolynomialQ[f, aa], Return[f]];
  expanded = Expand[f]; terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  Total[Select[terms, Total[Exponent[#, aa]] <= 1 &]]];

(* Units are inferred from expression trees, including incompatible sums. Empty
   support belongs to zero; it is not assigned an expected physical dimension. *)
unitSupport[0] := {};
unitSupport[x_?NumberQ] := {{0, 0, 0}};
unitSupport[x_Symbol] := If[KeyExistsQ[unitRegistry, x], {unitRegistry[x]}, {Missing["Units", x]}];
unitSupport[x_Plus] := DeleteDuplicates[Flatten[unitSupport /@ (List @@ x), 1]];
unitSupport[x_Times] := Module[{us = unitSupport /@ (List @@ x)},
  If[MemberQ[us, {}], {}, DeleteDuplicates[Map[If[AllTrue[#, ListQ], Total[#], Missing["ProductUnits", #]] &, Tuples[us]]]]];
unitSupport[Power[x_, n_?NumberQ]] := (If[ListQ[#], n #, #] & /@ unitSupport[x]);
unitSupport[x_] := {Missing["Units", HoldForm[x]]};
dimensionRecord[x_] := With[{support = unitSupport[x]},
  <|"SUPPORT" -> support, "CONSISTENCY_OPERAND" -> relationalObject[Length[support], 1]|>];
(* Metadata contains physical operands as well as labels. Inspect the operands
   recursively rather than assigning dimensionless units to a whole metadata tag.
   An unrecognized formal operator remains an explicit unknown. *)
dimensionTree[x_Association] := Map[dimensionTree, x];
dimensionTree[x_List] := dimensionTree /@ x;
dimensionTree[Rule[l_, r_]] := <|"DOMAIN" -> dimensionTree[l], "IMAGE" -> dimensionTree[r]|>;
dimensionTree[x_String] := <|"METADATA_SUPPORT" -> {{0, 0, 0}}|>;
dimensionTree[x_] := dimensionRecord[x];

(* Enumerate O(3)-even scalar contractions. Only the theta-bearing summands can
   enter either held-fixed theta derivative: the pullback of eW and u is theta-free.
   The dependency census records the complementary generated summands as well. *)
pairings[{}] := {{}};
pairings[ll_List] := Flatten[Table[
  (Prepend[#, {First[ll], ll[[j]]}] & /@ pairings[Delete[Rest[ll], j - 1]]),
  {j, 2, Length[ll]}], 1];
tensorValue[{name_, rank_}, indices_List] := Switch[name,
  "theta", theta, "eW", eW, "u", displacement[[First[indices]]],
  "gradTheta", jet["theta", indices], "gradEW", jet["eW", indices],
  "gradU", jet["u" <> ToString[indices[[1]]], {indices[[2]]}],
  "gradW", jet["WBg", indices], "gradMuR", jet["muRBg", indices]];
contractTensors[ts_List] := Module[{ranks = ts[[All, 2]], n, pp},
  n = Total[ranks]; If[OddQ[n], Return[{}]]; pp = pairings[Range[n]];
  DeleteDuplicates[Expand /@ Map[Function[pairing,
    Total[Map[Function[values, Module[{ix = ConstantArray[0, n], offset = 0},
      MapThread[(ix[[#1]] = {#2, #2}) &, {pairing, values}];
      Times @@ Map[Function[t, With[{part = Take[ix, {offset + 1, offset + t[[2]]}]},
        offset += t[[2]]; tensorValue[t, part]]], ts]]], Tuples[Range[3], n/2]]]], pp]]];
energyEulerVector[f_] := Table[D[f, jet[field]] - Sum[td[D[f, jet[field, {i}]], i], {i, 3}],
  {field, {"theta", "eW", "u1", "u2", "u3"}}];
divergenceQuotient[terms_List] := Module[{vectors, variables, signatures, monomials, matrix, rr, pivots},
  vectors = energyEulerVector /@ terms;
  variables = DeleteDuplicates[Flatten[atomsIn /@ Flatten[vectors]]];
  signatures = Map[Function[vector, Association[Flatten[MapIndexed[
    Function[{component, index}, Map[Function[rule,
      ToString[{First[index], First[rule]}, InputForm] -> Last[rule]], CoefficientRules[Expand[component], variables]]], vector]]]], vectors];
  monomials = Union @@ (Keys /@ signatures);
  matrix = Transpose[Lookup[#, monomials, 0] & /@ signatures];
  rr = RowReduce[matrix];
  pivots = DeleteCases[Map[Function[row, FirstPosition[row, x_ /; x =!= 0,
    Missing["EmptyRow"], {1}, Heads -> False]], rr], _Missing];
  If[!AllTrue[Flatten[pivots], IntegerQ[#] && 1 <= # <= Length[terms] &], Quit[101]];
  <|"BASIS" -> terms[[Flatten[pivots]]], "PIVOTS" -> Flatten[pivots],
    "EULER_SIGNATURE_RANK" -> Length[pivots], "MATRIX_DIMENSIONS" -> Dimensions[matrix],
    "SIGNATURE_FINGERPRINT" -> Hash[vectors, "SHA256"]|>];
constructEnergy[] := Module[{ts, candidates, retained, density = 0, coeff, uu, known, thetaAtoms, quotient},
  ts = {{"theta", 0}, {"eW", 0}, {"u", 1}, {"gradTheta", 1}, {"gradEW", 1}, {"gradU", 2}};
  candidates = DeleteDuplicates[Flatten[Table[
    Join[If[MemberQ[{ts[[i, 1]], ts[[j, 1]]}, "u"], {}, contractTensors[{ts[[i]], ts[[j]]}]],
      contractTensors[{{"gradW", 1}, ts[[i]], ts[[j]]}],
      contractTensors[{{"gradMuR", 1}, ts[[i]], ts[[j]]}]], {i, Length[ts]}, {j, i, Length[ts]}]]];
  candidates = candidates /. {theta^2 -> WBg theta^2, theta eW -> WBg theta eW};
  candidates = Join[{WBg theta^2, WBg theta eW}, Complement[candidates, {WBg theta^2, WBg theta eW}]];
  quotient = divergenceQuotient[candidates];
  thetaAtoms = {theta, Sequence @@ Table[jet["theta", {i}], {i, 3}]};
  retained = Select[quotient["BASIS"], !FreeQ[#, Alternatives @@ thetaAtoms] &];
  known = {WBg theta^2 -> bRho/2, WBg theta eW -> cCoupling};
  Do[uu = unitSupport[retained[[i]]];
    coeff = If[MemberQ[First /@ known, retained[[i]]], retained[[i]] /. known,
      declare[Symbol["energyCoefficient" <> ToString[i]], {-1, -2, 1} - First[uu]]];
    density += coeff retained[[i]], {i, Length[retained]}];
  <|"DENSITY" -> density, "GENERATED_CONTRACTIONS" -> candidates,
    "THETA_DEPENDENT_CONTRACTIONS" -> retained, "DIVERGENCE_QUOTIENT" -> quotient,
    "THETA_INDEPENDENT_CENSUS" -> Map[Function[f, <|"TERM" -> f,
      "THETA_DERIVATIVE" -> el[f]|>],
      Select[candidates, FreeQ[#, Alternatives @@ thetaAtoms] &]],
    "QUOTIENT_DISCARDED_CENSUS" -> Map[Function[f, <|"TERM" -> f,
      "EULER_SIGNATURE" -> energyEulerVector[f]|>], Complement[candidates, quotient["BASIS"]]],
    "REPRESENTATIVE" -> Inactive[Modulo][retained, Inactive[Div][compactSupportFlux]]|>];

phiMap[expression_, a_, h_] := Module[{domain, rules, uncovered},
  domain = Select[atomsIn[expression], MemberQ[{"theta", "eW"}, jetRegistry[#][[1]]] &];
  rules = Map[Function[atom, atom -> tdMany[
    If[jetRegistry[atom][[1]] === "theta", theta + a, eW + h], jetRegistry[atom][[2]]]], domain];
  uncovered = Complement[domain, First /@ rules];
  <|"MAP" -> rules, "DOMAIN" -> domain,
    "COVERAGE" -> Association[Map[# -> relationalObject[Count[First /@ rules, #], 1] &, domain]],
    "UNCOVERED" -> uncovered, "MAX_RANK" -> Max[Prepend[Length[jetRegistry[#][[2]]] & /@ domain, 0]]|>];
materialAmplitude[density_, a_, h_, ak_, hk_, jk_] := Module[{mapping, pulled},
  mapping = phiMap[density, ak a, hk h];
  pulled = waveScale[(1 + Sum[jet["u" <> ToString[i], {i}], {i, 3}]) (density /. mapping["MAP"]), 2];
  <|"DENSITY" -> pulled, "MU" -> el[pulled] + jk junkMu eW, "MAP" -> mapping|>];

(* Two independent coordinate builders. Eulerian differentiates the graph level
   set; material differentiates the flattening coordinate and maps its covector
   through the four-dimensional Jacobian INSIDE the builder. *)
linearShape[f_] := Coefficient[Expand[f], shapeParameter, 1];
shapeBackground[f_] := f /. shapeParameter -> 0;
firstBackgroundJet[f_] := Module[{aa, rules},
  aa = Select[atomsIn[f], MemberQ[{"WBg", "muRBg"}, jetRegistry[#][[1]]] &&
    Length[jetRegistry[#][[2]]] > 0 &];
  rules = Thread[aa -> backgroundJetMarker aa];
  Normal[Series[f /. rules, {backgroundJetMarker, 0, 1}]] /. backgroundJetMarker -> 1];
geometryExpansion[f_List] := geometryExpansion /@ f;
geometryExpansion[f_] := firstBackgroundJet[Normal[Series[f, {shapeParameter, 0, 1}]]];
graphGeometry[anchor_, s_] := Module[{height, normalRaw, normal, area, vv, virtual, slope},
  height = s (WBg - If[anchor === "MATERIAL_ADVECTED", shapeParameter displacement.Table[td[WBg, i], {i, 3}], 0])/2 +
    shapeParameter (zetaC + s W0 eW/2);
  slope = Table[td[height, i], {i, 3}];
  normalRaw = s Join[-slope, {1}]; area = Sqrt[1 + slope.slope];
  normal = geometryExpansion[normalRaw/area]; area = geometryExpansion[area];
  vv = Join[shapeParameter Table[td[displacement[[i]], 4], {i, 3}],
    {shapeParameter (td[zetaC + s W0 eW/2, 4] +
       If[anchor === "LAB_HELD", s Table[td[WBg, i], {i, 3}].Table[td[displacement[[i]], 4], {i, 3}]/2, 0])}];
  virtual = Join[virtualU, {virtualCenter + s W0 virtualEW/2 +
    If[anchor === "LAB_HELD", s Table[td[WBg, i], {i, 3}].virtualU/2, 0]}];
  <|"HEIGHT" -> height, "CONORMAL" -> normalRaw, "NORMAL" -> normal,
    "AREA" -> area, "FACE_VECTOR" -> vv, "VIRTUAL" -> virtual,
    "VELOCITY" -> linearShape[normal.vv], "ROUTE" -> "EULERIAN_LEVEL_SET"|>];
materialGeometry[anchor_, s_, knife_] := Module[{deform, invT, thick, center, flattenCovector,
    mapped, norm, area, virtual, vv, height, slope, normal, inverseZero, inverseFirst},
  deform = IdentityMatrix[3] + shapeParameter Table[td[displacement[[i]], j], {i, 3}, {j, 3}];
  thick = WBg + shapeParameter (W0 eW + If[anchor === "LAB_HELD", displacement.Table[td[WBg, i], {i, 3}], 0]);
  center = shapeParameter zetaC;
  (* d_X w' at w'=s/2, d_w w'; lower block of inverse transpose is exact. *)
  flattenCovector = Join[Table[-(td[center, i] + s td[thick, i]/2)/thick, {i, 3}], {1/thick}];
  inverseZero = LinearSolve[Transpose[deform /. shapeParameter -> 0], IdentityMatrix[3]];
  inverseFirst = LinearSolve[Transpose[deform /. shapeParameter -> 0],
    -Transpose[D[deform, shapeParameter] /. shapeParameter -> 0].inverseZero];
  invT = ArrayFlatten[{{inverseZero + shapeParameter inverseFirst, ConstantArray[0, {3, 1}]},
    {ConstantArray[0, {1, 3}], {{1}}}}];
  (* FORM knife: contaminate one mapped covector component by a different slope.
     Only this carrier builder receives knife; source velocity has a fresh builder. *)
  invT[[1, 4]] += knife jet["WBg", {2}];
  mapped = s thick invT.flattenCovector;
  (* Return from X to x before forming source laws. *)
  mapped = Normal[Series[mapped, {shapeParameter, 0, 1}]];
  mapped = mapped /. (a_Symbol /; KeyExistsQ[jetRegistry, a] && jetRegistry[a][[1]] === "WBg") :>
    a - shapeParameter displacement.Table[td[a, i], {i, 3}];
  mapped = Normal[Series[mapped, {shapeParameter, 0, 1}]];
  norm = Sqrt[mapped.mapped]; normal = geometryExpansion[mapped/norm];
  area = geometryExpansion[norm]; (* common Eulerian graph conormal has last component s *)
  virtual = Join[virtualU, {virtualCenter + s W0 virtualEW/2 +
    If[anchor === "LAB_HELD", s Table[td[WBg, i], {i, 3}].virtualU/2, 0]}];
  vv = Join[shapeParameter Table[td[displacement[[i]], 4], {i, 3}],
    {linearShape[td[center + s thick/2, 4]] shapeParameter}];
  <|"FLATTEN_COVECTOR" -> flattenCovector, "INVERSE_TRANSPOSE" -> invT,
    "CONORMAL" -> mapped, "NORMAL" -> normal, "AREA" -> area,
    "FACE_VECTOR" -> vv, "VIRTUAL" -> virtual, "VELOCITY" -> linearShape[normal.vv],
    "ROUTE" -> "MATERIAL_FLATTENING"|>];
virtualU = {virtualU1, virtualU2, virtualU3};
Scan[declare[#, {1, 0, 0}] &, Join[virtualU, {virtualCenter}]];
declare[virtualEW, {0, 0, 0}];
faceLaws[geometry_, s_, mu_, density_] := Module[{p, trace, normal, area, velocity, bulk,
    affinity, flux, closure, traction, work, balance},
  p = If[s === 1, Take[pressureSlots, 2], Take[pressureSlots, -2]];
  trace = shapeParameter (p[[1]] + s (WBg - W0) p[[2]]/2);
  normal = geometry["NORMAL"]; area = geometry["AREA"];
  velocity = geometry["FACE_VECTOR"];
  bulk = shapeParameter Table[declare[Symbol["bulkVelocity" <> ToString[s /. -1 -> 2] <> ToString[i]], {1, -1, 0}], {i, 4}];
  affinity = shapeParameter mu/density - trace/rhoM;
  flux = rhoM (bulk - velocity).normal;
  closure = flux - lambdaA affinity - lambdaV normal.velocity;
  traction = -(trace + lambdaX affinity) normal;
  work = area traction.geometry["VIRTUAL"];
  <|"TRACTION" -> linearShape[traction], "VIRTUAL_WORK" -> linearShape[work],
    "RELATIVE_FLUX" -> linearShape[flux], "TRUE_AREA_FLUX" -> linearShape[area flux],
    "CLOSURE" -> linearShape[closure], "KINEMATIC" -> linearShape[normal.bulk - normal.velocity - flux/rhoM],
    "AREA_DERIVATIVE" -> linearShape[area], "NORMAL" -> shapeBackground[normal],
    "CONORMAL" -> shapeBackground[geometry["CONORMAL"]],
    "VELOCITY" -> geometry["VELOCITY"], "PRESSURE_IDENTITIES" -> p|>];
eulerianSlabFace[laws_, massBase_] := Module[{work, rows},
  work = laws["VIRTUAL_WORK"];
  rows = Join[Table[-D[work, virtualU[[i]]], {i, 3}],
    {massBase + laws["TRUE_AREA_FLUX"] - laws["CLOSURE"], -D[work, virtualEW]}];
  rows];
materialFaceFold[laws_, massBase_] := Module[{correctedOrigins, forces},
  correctedOrigins = laws["TRUE_AREA_FLUX"] - laws["CLOSURE"];
  forces = -Table[Coefficient[Expand[laws["VIRTUAL_WORK"]], v], {v, Join[virtualU, {virtualEW}]}];
  Join[Take[forces, 3], {massBase + correctedOrigins, Last[forces]}]];
carrier[rows_] := Table[D[rows[[r]], p] /. slotZero, {r, 5}, {p, pressureSlots}];
denominatorCircuit[f_] := Times @@ Cases[f, Power[base_, exponent_Integer /; exponent < 0] :>
  base^(-exponent), {0, Infinity}];
massSubstrate[anchor_, density4_] := Module[{deform, jacobian, sigmaE, sigmaM, carriedE,
    massE, massM, vConstraintE, vConstraintM, gradW, gradR, virtualRules},
  deform = IdentityMatrix[3] + shapeParameter Table[td[displacement[[i]], j], {i, 3}, {j, 3}];
  jacobian = Det[deform];
  gradW = Table[td[WBg, i], {i, 3}]; gradR = Table[td[density4, i], {i, 3}];
  sigmaE = (density4 - If[anchor === "MATERIAL_ADVECTED", shapeParameter displacement.gradR, 0])
    (1 + shapeParameter theta) (WBg + shapeParameter (W0 eW -
      If[anchor === "MATERIAL_ADVECTED", displacement.gradW, 0]));
  sigmaM = (density4 + If[anchor === "LAB_HELD", shapeParameter displacement.gradR, 0])
    (1 + shapeParameter theta) (WBg + shapeParameter (W0 eW +
      If[anchor === "LAB_HELD", displacement.gradW, 0])) jacobian;
  carriedE = (sigmaE + shapeParameter displacement.Table[td[sigmaE /. shapeParameter -> 0, i], {i, 3}]) jacobian;
  massE = linearShape[td[sigmaE, 4] + Sum[td[sigmaE shapeParameter td[displacement[[i]], 4], i], {i, 3}]];
  massM = linearShape[td[sigmaM, 4]];
  vConstraintE = linearShape[carriedE]/(density4 WBg);
  vConstraintM = linearShape[sigmaM]/(density4 WBg);
  declare[virtualTheta, {0, 0, 0}];
  virtualRules = Join[{theta -> virtualTheta, eW -> virtualEW, zetaC -> virtualCenter},
    Thread[displacement -> virtualU], Flatten[Table[jet["u" <> ToString[i], {j}] ->
      declare[Symbol["virtualU" <> ToString[i] <> "Jet" <> ToString[j]], {0, 0, 0}], {i, 3}, {j, 3}]]];
  <|"EULERIAN_EVOLUTION_BASE" -> massE, "MATERIAL_EVOLUTION_BASE" -> massM,
    "EULERIAN_VIRTUAL_CONSTRAINT" -> (vConstraintE /. virtualRules),
    "MATERIAL_VIRTUAL_CONSTRAINT" -> (vConstraintM /. virtualRules), "VIRTUAL_FIELD_IDENTITY" -> virtualRules,
    "CENTER_WORK_SLOT" -> virtualCenter, "DENSITY_JACOBIAN" -> jacobian|>];

(* The c1 source is extracted from the supplied closure solve, not re-entered
   as a second formula. Operator order is retained before Fourier normalization. *)
sourceConstruction[] := Module[{fluxLaw, bulkLaw, driving, pRule, resolvent, closedPressure, closedFlux},
  fluxLaw = jFace == lambdaA (muTheta/rhoFace - pFace/rhoM) + lambdaV vFace;
  bulkLaw = vBulkFace == vFace + jFace/rhoM;
  driving = vBulkFace /. First[Solve[{fluxLaw, bulkLaw}, {jFace, vBulkFace}]];
  resolvent = Inactive[Inverse][Inactive[Plus][IdentityOperator,
    Inactive[NonCommutativeMultiply][-Coefficient[driving, pFace], impedanceOperator]]];
  closedPressure = Inactive[NonCommutativeMultiply][resolvent, impedanceOperator, driving /. pFace -> 0];
  closedFlux = (jFace /. First[Solve[fluxLaw, jFace]]) /. pFace -> closedPressure;
  <|"SOURCE" -> (driving /. pFace -> 0), "PRESSURE_COEFFICIENT" -> Coefficient[driving, pFace],
    "DRIVING" -> driving,
    "OPERATOR_EQUATION" -> relationalObject[pFace,
      Inactive[NonCommutativeMultiply][impedanceOperator, driving]],
    "RESOLVENT" -> resolvent, "CLOSED_PRESSURE" -> closedPressure, "CLOSED_FLUX" -> closedFlux,
    "CLOSED_TRACTION" -> -(closedPressure + lambdaX (muTheta/rhoFace - closedPressure/rhoM)) faceNormalOperand|>];
sourceBind[solve_, mu_, velocity_, density_] :=
  waveScale[pruneBackgroundProducts[solve["SOURCE"] /. {muTheta -> mu, vFace -> velocity, rhoFace -> density}], 1];

(* Boundary matching is performed on a generic outgoing input wave. Its pressure
   trace and conormal trace give Z = P N^-1, expanded as ordered operators. *)
momentum = Association[Table[leg -> Table[declare[Symbol[leg <> ToString[i]], {-1, 0, 0}], {i, 3}],
  {leg, {"kOut", "kIn", "kMiddle"}}]];
qLeg = Association[Map[# -> declare[Symbol["q" <> #], {-1, 0, 0}] &, {"Out", "In", "Middle"}]];
constructKernel[s_Integer] := Module[{zz, height, heightJets, p0, n0, p1, n1, z1,
    coordinates = {acousticX1, acousticX2, acousticX3}, wave, flatWave, pressureField,
    normalTrace, pressureTrace, dispersionEquation, dispersionRoot, mixedNormal},
  zz = <||>;
  height = declare[heightTransform, {4, 0, 0}];
  heightJets = Table[declare[Symbol["heightJetTransform" <> ToString[i]], {3, 0, 0}], {i, 3}];
  (* acousticW is the laboratory normal coordinate. The supplied orientation
     and graph map enter separately for each disconnected exterior half-space. *)
  wave = Exp[I (momentum["kIn"].coordinates + s qLeg["In"] acousticW - omega acousticTime)];
  flatWave = wave /. acousticW -> 0;
  pressureField = -rhoM D[wave, acousticTime];
  pressureTrace = pressureField /. acousticW -> s shapeParameter height;
  normalTrace = (s D[wave, acousticW] - shapeParameter Sum[heightJets[[i]] D[wave, coordinates[[i]]], {i, 3}]) /.
    acousticW -> s shapeParameter height;
  p0 = Cancel[(pressureTrace /. shapeParameter -> 0)/flatWave];
  n0 = Cancel[(normalTrace /. shapeParameter -> 0)/flatWave];
  p1 = Cancel[(D[pressureTrace, shapeParameter] /. shapeParameter -> 0)/flatWave];
  n1 = Cancel[(D[normalTrace, shapeParameter] /. shapeParameter -> 0)/flatWave];
  mixedNormal = Cancel[(D[(s D[wave, acousticW] - slopeMarker Sum[heightJets[[i]] D[wave, coordinates[[i]]], {i, 3}]) /.
    acousticW -> s heightMarker height, heightMarker, slopeMarker] /. {heightMarker -> 0, slopeMarker -> 0})/flatWave];
  Do[AssociateTo[zz, leg -> (Cancel[p0/n0] /. qLeg["In"] -> qLeg[leg])], {leg, Keys[qLeg]}];
  dispersionEquation = Expand[Cancel[(D[wave, {acousticTime, 2}] - cS0^2 (
    Sum[D[wave, {coordinates[[i]], 2}], {i, 3}] + D[wave, {acousticW, 2}]))/wave]];
  dispersionRoot = qSquared /. First[Solve[(dispersionEquation /. qLeg["In"]^2 -> qSquared) == 0, qSquared]];
  z1 = p1/n0 - zz["Out"] n1/n0;
  <|"FACE" -> s, "FLAT" -> zz, "SHAPE" -> z1,
    "TRACE_PRESSURE" -> {p0, p1}, "TRACE_CONORMAL" -> {n0, n1}, "TRACE_MIXED_CONORMAL" -> mixedNormal,
    "JET_FOURIER_IDENTITY" -> Thread[heightJets -> I (momentum["kOut"] - momentum["kIn"]) height],
    "BULK_EQUATION_OPERAND" -> dispersionEquation,
    "DISPERSION" -> Association[Map[# -> relationalObject[qLeg[#]^2,
      dispersionRoot /. Thread[momentum["kIn"] -> momentum["k" <> #]]] &, Keys[qLeg]]]|>];

(* Exact, componentwise arithmetic nodes preserve subtraction until sampling.
   No cancellation or full symbolic zero-test is used for comparison nodes. *)
add[a_, b_] := arithmeticPlus[a, b];
add[0, b_] := b;
add[a_, 0] := a;
neg[a_] := mul[-1, a];
sub[a_, b_] := arithmeticPlus[a, neg[b]];
mul[a_, b_] := arithmeticTimes[a, b];
mul[0, b_] := 0;
mul[a_, 0] := 0;
mul[1, b_] := b;
mul[a_, 1] := a;
power[a_, n_Integer] := arithmeticPower[a, n];
circuitExpression[x_arithmeticPlus] := Inactive[Plus] @@ (circuitExpression /@ List @@ x);
circuitExpression[x_arithmeticTimes] := Inactive[Times] @@ (circuitExpression /@ List @@ x);
circuitExpression[arithmeticPower[a_, n_]] := Inactive[Power][circuitExpression[a], n];
circuitExpression[x_] := x;
gAdd[a_List, b_List] := MapThread[add, {a, b}];
gSub[a_List, b_List] := MapThread[sub, {a, b}];
gMul[a_List, b_List] := Table[Fold[add, 0,
  Flatten[Table[If[gradeIndices[[i]] + gradeIndices[[j]] === g, mul[a[[i]], b[[j]]], Nothing],
    {i, 4}, {j, 4}]]], {g, gradeIndices}];
gScale[a_List, b_] := mul[#, b] & /@ a;
graded[f_] := gradePart[f, #] & /@ gradeIndices;
gInverse[a_List] := Module[{out = ConstantArray[0, 4], earlier},
  out[[1]] = power[a[[1]], -1];
  Do[earlier = Fold[add, 0, Flatten[Table[
    If[gradeIndices[[i]] + gradeIndices[[j]] === gradeIndices[[k]], mul[a[[i]], out[[j]]], Nothing],
      {i, 2, 4}, {j, 1, k - 1}]]]; out[[k]] = neg[mul[out[[1]], earlier]], {k, 2, 4}]; out];
gPower[a_List, 0] := {1, 0, 0, 0};
gPower[a_List, n_Integer] := If[n < 0, gPower[gInverse[a], -n], Fold[gMul, {1, 0, 0, 0}, ConstantArray[a, n]]];
gradeCircuit[f_, replacements_Association] := Module[{h = Head[f]},
  If[KeyExistsQ[replacements, f], Return[replacements[f]]];
  If[FreeQ[f, Alternatives @@ Join[{etaBg, sigmaW}, Keys[replacements]]], Return[{f, 0, 0, 0}]];
  Which[f === etaBg, {0, 0, 1, 0}, f === sigmaW, {0, 1, 0, 0},
    h === Plus, Fold[gAdd, {0, 0, 0, 0}, gradeCircuit[#, replacements] & /@ List @@ f],
    h === Times, Fold[gMul, {1, 0, 0, 0}, gradeCircuit[#, replacements] & /@ List @@ f],
    h === Power && IntegerQ[f[[2]]], gPower[gradeCircuit[f[[1]], replacements], f[[2]]],
    True, WriteString[$Messages[[1]], "Unsupported graded arithmetic\n"]; Quit[93]]];

(* Local differential trial parametrizations: u_T = Curl A, u_L = Grad psi.
   Transferring derivatives off arbitrary compact-support trial potentials gives
   (-D_Y + i k_in)^I on each live coefficient. This is not a Helmholtz projector. *)
adjointJet[c_, indices_, leg_] := Fold[Function[{f, i},
  pruneBackgroundProducts[If[i === 4, -I omega f, -td[f, i] + I momentum[leg][[i]] f]]], c, indices];
restrict[expression_, trial_] := restrict[expression, trial] = Module[{aa, out, f, ix, c, j},
  aa = Select[atomsIn[expression], MemberQ[Switch[trial,
    "THETA", {"theta"}, "E_W", {"eW"}, _, {"u1", "u2", "u3"}], jetRegistry[#][[1]]] &];
  out = ConstantArray[0, If[trial === "TRANSVERSE", 3, 1]];
  Do[f = jetRegistry[a][[1]]; ix = jetRegistry[a][[2]]; c = D[expression, a];
    Which[trial === "THETA" && f === "theta", out[[1]] += adjointJet[c, ix, "kIn"],
      trial === "E_W" && f === "eW", out[[1]] += adjointJet[c, ix, "kIn"],
      MemberQ[{"u1", "u2", "u3"}, f], j = ToExpression[StringDrop[f, 1]];
        If[trial === "LONGITUDINAL", out[[1]] += adjointJet[c, Append[ix, j], "kIn"]];
        If[trial === "TRANSVERSE", Do[out[[k]] += Signature[{j, d, k}]
          adjointJet[c, Append[ix, d], "kIn"], {d, 3}, {k, 3}]]], {a, aa}];
  finish /@ out];
testCarrier[cs_, test_, slot_] := testCarrier[cs, test, slot] = Module[{v = cs[[1 ;; 3, slot]], result},
  result = Switch[test, "THETA", {cs[[4, slot]]}, "E_W", {cs[[5, slot]]},
    "LONGITUDINAL", {-Sum[td[v[[i]], i] + I momentum["kOut"][[i]] v[[i]], {i, 3}]},
    "TRANSVERSE", Table[Sum[Signature[{k, i, j}]
      (td[v[[j]], i] + I momentum["kOut"][[i]] v[[j]]), {i, 3}, {j, 3}], {k, 3}]];
  finish /@ result];
blocks = {{"THETA", "TRANSVERSE"}, {"E_W", "TRANSVERSE"}, {"LONGITUDINAL", "TRANSVERSE"},
  {"TRANSVERSE", "THETA"}, {"TRANSVERSE", "E_W"}, {"TRANSVERSE", "LONGITUDINAL"}};
kernelFamilies = {"LOCAL_BARE", "FOURIER_KOUT_Y", "FOURIER_KOUT_KIN_Y", "FOURIER_KOUT_KIN_Y_MIDDLE"};
Scan[declare[#, {-3, 0, 0}] &, {d3KOut, d3KIn, d3KMiddle}];
declare[d3Y, {3, 0, 0}]; declare[Pi, {0, 0, 0}];
xCoordinate = Table[declare[Symbol["xCoordinate" <> ToString[i]], {1, 0, 0}], {i, 3}];
yCoordinate = Table[declare[Symbol["yCoordinate" <> ToString[i]], {1, 0, 0}], {i, 3}];
measureCircuit[family_] := Switch[family, "LOCAL_BARE", 1,
  "FOURIER_KOUT_Y", d3KOut d3Y/(2 Pi)^3,
  "FOURIER_KOUT_KIN_Y", d3KOut d3KIn d3Y/(2 Pi)^6,
  "FOURIER_KOUT_KIN_Y_MIDDLE", d3KOut d3KIn d3Y d3KMiddle/(2 Pi)^9];
keyDimensionRecord[value_, key_List] := Module[{family, measure, pairing = {0, 0, 0}, sectors, support},
  family = SelectFirst[key, MemberQ[kernelFamilies, #] &, "LOCAL_BARE"];
  measure = unitSupport[measureCircuit[family]];
  If[StringQ[First[key]] && StringContainsQ[First[key], "_FROM_"],
    sectors = StringSplit[First[key], "_FROM_"];
    If[family === "LOCAL_BARE", sectors = Take[sectors, 1]];
    pairing = Total[If[MemberQ[{"TRANSVERSE", "LONGITUDINAL"}, #], {2, 0, 0}, {0, 0, 0}] & /@ sectors]];
  support = unitSupport[value];
  <|"COEFFICIENT_SUPPORT" -> support, "MEASURE_SUPPORT" -> measure, "TRIAL_TEST_SUPPORT" -> pairing,
    "PAIRED_SUPPORT" -> Map[If[ListQ[#], # + First[measure] + pairing, #] &, support],
    "ADDITION_CONSISTENCY" -> relationalObject[Length[support], 1]|>];

(* Formal response families are derived by the ordered inverse expansion. The
   middle family retains both cross products of eta and sigma corrections. *)
responseFamilies[kernel_] := Module[{p0, n0, d, a, profile, jets, traceP, traceN, pg, ng, dg,
    first, leftRules, rightRules, leftP, leftD, rightD, mixedTrace, mixed, heightLeft, gradientRight},
  a = lambdaA/rhoM^2; p0 = First[kernel["TRACE_PRESSURE"]]; n0 = First[kernel["TRACE_CONORMAL"]];
  d = Association[Map[# -> ((n0 + a p0) /. qLeg["In"] -> qLeg[#]) &, Keys[qLeg]]];
  profile = declare[profileTransform, {3, 0, 0}];
  jets = Table[Symbol["heightJetTransform" <> ToString[i]], {i, 3}];
  traceP = Last[kernel["TRACE_PRESSURE"]] /. heightTransform -> etaBg W0 profile/2;
  traceN = Last[kernel["TRACE_CONORMAL"]] /. heightTransform -> etaBg W0 profile/2;
  traceN = traceN /. Thread[jets -> sigmaW LW I (momentum["kOut"] - momentum["kIn"]) profile/2];
  pg = graded[retain[traceP]]; ng = graded[retain[traceN]]; dg = pg a + ng;
  (* P (N+a P)^-1 is the boundary-matching form of [I+a Z]^-1 Z.
     Each product below retains its own intermediate momentum; there is no scalar
     replacement of the nonlocal inverse. The inverse series is generated by
     multiplying its defining equation grade by grade. *)
  first = gScale[gSub[pg, gScale[dg, p0/d["Out"]]], 1/d["In"]];
  leftRules = Join[Thread[momentum["kIn"] -> momentum["kMiddle"]],
    {qLeg["In"] -> qLeg["Middle"], profileTransform -> declare[profileTransformLeft, {3, 0, 0}]}];
  rightRules = Join[Thread[momentum["kOut"] -> momentum["kMiddle"]],
    {qLeg["Out"] -> qLeg["Middle"], profileTransform -> declare[profileTransformRight, {3, 0, 0}]}];
  leftP = pg /. leftRules; leftD = dg /. leftRules; rightD = dg /. rightRules;
  heightLeft = etaBg W0 profileTransformLeft/2;
  gradientRight = sigmaW LW I (momentum["kMiddle"] - momentum["kIn"]) profileTransformRight/2;
  mixedTrace = kernel["TRACE_MIXED_CONORMAL"] /. Join[{heightTransform -> heightLeft}, Thread[jets -> gradientRight]];
  mixed = gAdd[gScale[graded[retain[mixedTrace]], -p0/(d["Out"] d["In"])],
    gAdd[gScale[gMul[leftP, rightD], -1/(d["Middle"] d["In"])],
      gScale[gMul[leftD, rightD], p0/(d["Out"] d["Middle"] d["In"])]]];
  <|"FOURIER_KOUT_Y" -> {p0/d["Out"], 0, 0, 0},
    "FOURIER_KOUT_KIN_Y" -> first, "FOURIER_KOUT_KIN_Y_MIDDLE" -> mixed,
    "DENOMINATORS" -> d, "TRACE_PRESSURE_GRADES" -> pg, "TRACE_CONORMAL_GRADES" -> ng,
    "TRACE_MIXED_CONORMAL_GRADES" -> graded[retain[mixedTrace]],
    "ORDERED_MATCHING_OPERAND" -> Inactive[NonCommutativeMultiply][pressureTraceOperator,
      Inactive[Inverse][normalTraceOperator + a pressureTraceOperator]]|>];

formalMeasure[family_] := Switch[family,
  "LOCAL_BARE", 1,
  "FOURIER_KOUT_Y", Inactive[Integral][Inactive[Exp][I momentum["kOut"].(xCoordinate - yCoordinate)],
    {kOut, yCoordinate}, Inactive[Times][d3KOut, d3Y, (2 Pi)^-3]],
  "FOURIER_KOUT_KIN_Y", Inactive[Integral][Inactive[Exp][I (momentum["kOut"].xCoordinate - momentum["kIn"].yCoordinate)],
    {kOut, kIn, yCoordinate}, Inactive[Times][d3KOut, d3KIn, d3Y, (2 Pi)^-6]],
  "FOURIER_KOUT_KIN_Y_MIDDLE", Inactive[Integral][Inactive[Exp][I (momentum["kOut"].xCoordinate - momentum["kIn"].yCoordinate)],
    {kOut, kIn, yCoordinate, kMiddle}, Inactive[Times][d3KOut, d3KIn, d3Y, d3KMiddle, (2 Pi)^-9]]];
normalContinuation[s_] := normalContinuation[s] = Cancel[
  D[Exp[I s qLeg["Out"] continuationW], continuationW]/Exp[I s qLeg["Out"] continuationW]];

(* Assembly keeps the pressure-independent base in the diagnostic circuit. The
   weak increment itself is formed only from P coefficients. *)
buildContraction[cs_, source_, response_, s_, affine_] := Module[{result = <||>, tc, sc, sourceGrade,
    pressureIndex, slot, multiplier, cg, combined, entry, family, pair, block, normalFactor, trialComponents},
  pressureIndex = If[s === 1, {1, 2}, {3, 4}];
  Do[block = pair[[1]] <> "_FROM_" <> pair[[2]];
    sc = atPoint[#, "Y"] & /@ restrict[source, pair[[2]]];
    Do[family = familyName;
      trialComponents = If[family === "LOCAL_BARE", 1, Length[sc]];
      combined = ConstantArray[0, {If[pair[[1]] === "TRANSVERSE", 3, 1], trialComponents, 4}];
      Do[slot = pressureIndex[[slotNumber]]; tc = atPoint[#, "X"] & /@ testCarrier[cs, pair[[1]], slot];
        normalFactor = If[slotNumber === 1, 1, normalContinuation[s]];
        Do[cg = graded[tc[[i]]]; sourceGrade = graded[If[family === "FOURIER_KOUT_Y",
          sc[[j]] /. Thread[momentum["kIn"] -> momentum["kOut"]], sc[[j]]]];
          entry = If[family === "LOCAL_BARE", ConstantArray[0, 4],
            gScale[gMul[gMul[cg, response[family]], sourceGrade], normalFactor]];
          (* Restriction of an AFFINE map retains its constant forcing in the test
             space. It has one output vector, not a trial-component matrix. Taking
             a trial derivative here would incorrectly delete the local -C.p. *)
          If[family === "LOCAL_BARE" && affine,
            entry = gScale[cg, -pressureSlots[[slot]]]];
          combined[[i, j]] = gAdd[combined[[i, j]], entry],
          {i, Length[tc]}, {j, trialComponents}], {slotNumber, 2}];
      Do[AssociateTo[result, {block, family, {1, Sequence @@ gradeIndices[[g]]}, s} ->
        Flatten[combined[[All, All, g]]]], {g, 4}], {familyName, If[affine, kernelFamilies, Rest[kernelFamilies]]}];
    constructionEventIndex++;
    appendAssociationEmission[{caseOrdinal, "WEAK_BLOCK", constructionEventIndex, s}, {block, Length[result]}],
    {pair, blocks}]; result];
mapCombine[a_Association, b_Association, fn_] := Association[Map[# -> MapThread[fn, {a[#], b[#]}] &, Keys[a]]];

emit[{"LOCAL", "SETUP"}, <|"CASES" -> Tuples[{{"LAB_HELD", "MATERIAL_ADVECTED"},
  {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}], "GRADES" -> ({1, Sequence @@ #} & /@ gradeIndices),
  "PRESSURE_SLOTS" -> pressureSlots, "PRESSURE_ASSUMPTIONS" -> pressureAssumptions,
  "KERNEL_FAMILIES" -> Association[# -> formalMeasure[#] & /@ kernelFamilies],
  "WEAK_BLOCKS" -> blocks|>];

(* Remaining run driver follows the symbolic construction. *)

(* A rational circuit is evaluated as numerator/denominator pairs. Addition does
   not cancel denominators; the same arithmetic supplies conservative degrees.
   Comparisons retain their two inputs even when their evaluated numerators vanish. *)
unitSupport[x_arithmeticPlus] := DeleteDuplicates[Flatten[unitSupport /@ List @@ x, 1]];
unitSupport[x_arithmeticTimes] := unitSupport[Times @@ List @@ x];
unitSupport[arithmeticPower[a_, n_]] := (If[ListQ[#], n #, #] & /@ unitSupport[a]);
numberPair[x_Integer, prime_] := {Mod[x, prime], 1};
numberPair[x_Rational, prime_] := {Mod[Numerator[x], prime], Mod[Denominator[x], prime]};
numberPair[x_Complex, prime_] := pairAdd[numberPair[Re[x], prime],
  pairMultiply[{imaginaryRoot, 1}, numberPair[Im[x], prime], prime], prime];
pairAdd[a_, b_, p_] := Mod[{a[[1]] b[[2]] + b[[1]] a[[2]], a[[2]] b[[2]]}, p];
pairMultiply[a_, b_, p_] := Mod[a b, p];
pairPower[a_, n_Integer, p_] := If[n >= 0,
  PowerMod[#, n, p] & /@ a, Reverse[PowerMod[#, -n, p] & /@ a]];
numericPair[x_?NumberQ, p_] := numberPair[x, p];
numericPair[x_Symbol, p_] := {Lookup[sampleValues, x, Missing["Sample", x]], 1};
numericPair[x_Plus | x_arithmeticPlus, p_] := Fold[pairAdd[#1, #2, p] &, {0, 1}, numericPair[#, p] & /@ List @@ x];
numericPair[x_Times | x_arithmeticTimes, p_] := Fold[pairMultiply[#1, #2, p] &, {1, 1}, numericPair[#, p] & /@ List @@ x];
numericPair[Power[x_, n_Integer], p_] := pairPower[numericPair[x, p], n, p];
numericPair[arithmeticPower[x_, n_Integer], p_] := pairPower[numericPair[x, p], n, p];
numericPair[x_, p_] := (WriteString[$Messages[[1]], "Unsupported arithmetic head: ", ToString[Head[x]], "\n"]; Quit[94]);
degreePair[x_?NumberQ] := {0, 0};
degreePair[x_Symbol] := Lookup[parameterDegrees, x, {1, 0}];
degreePair[x_Plus | x_arithmeticPlus] := Fold[
  {Max[#1[[1]] + #2[[2]], #2[[1]] + #1[[2]]], #1[[2]] + #2[[2]]} &, {0, 0}, degreePair /@ List @@ x];
degreePair[x_Times | x_arithmeticTimes] := Total[degreePair /@ List @@ x];
degreePair[Power[x_, n_Integer]] := If[n >= 0, n degreePair[x], -n Reverse[degreePair[x]]];
degreePair[arithmeticPower[x_, n_Integer]] := If[n >= 0, n degreePair[x], -n Reverse[degreePair[x]]];
degreePair[x_] := (WriteString[$Messages[[1]], "Unsupported degree head: ", ToString[Head[x]], "\n"]; Quit[95]);

(* Real branch charts are chosen before modular evaluation. Every leg uses fresh
   chart coordinates, with the same frequency and sound speed. Positive chart
   coordinates select a Zariski-dense real patch in each indicated branch cell.
   P: stereographic sphere with t_i=v_i/(1+sum v), hence sum t_i^2<1.
   E: rational hyperboloid, k.k-lambda^2=(omega/cS0)^2, lambda>0.
   Static E: a sphere of independent radius lambda. These are separate real cells.
   The chosen primes split i^2+1, so i is an exact element of F_p. *)
branchCells = Join[Flatten[Table[{sgn, regimes}, {sgn, {-1, 1}},
  {regimes, Tuples[{"PROPAGATING", "EVANESCENT"}, 3]}], 1],
  {{0, ConstantArray["EVANESCENT", 3]}}];
chartNames = Flatten[Table[Table[Symbol["chart" <> leg <> ToString[i]], {i, 4}], {leg, {"Out", "In", "Middle"}}]];
chartRules[cell_] := Module[{h = frequencyScale/cS0, om, rules, v, t, radius, den, kv, qv, regime},
  om = cell[[1]] frequencyScale;
  rules = {omega -> om};
  Do[v = Table[Symbol["chart" <> leg <> ToString[i]], {i, 4}]; regime = cell[[2, j]];
    If[cell[[1]] === 0,
      radius = v[[4]]; den = 1 + v[[1]]^2 + v[[2]]^2;
      kv = radius {2 v[[1]], 2 v[[2]], 1 - v[[1]]^2 - v[[2]]^2}/den;
      qv = I radius,
      If[regime === "PROPAGATING",
        t = Take[v, 3]/(1 + Total[Take[v, 3]]); den = 1 + t.t;
        kv = 2 h t/den; qv = cell[[1]] h (1 - t.t)/den,
        t = 1/(1 + v[[1]]);
        kv = h {(1 + t^2 - v[[2]]^2 - v[[3]]^2)/(2 t), v[[2]], v[[3]]};
        qv = I h (1 - t^2 + v[[2]]^2 + v[[3]]^2)/(2 t)]];
    rules = Join[rules, Thread[momentum["k" <> leg] -> kv], {qLeg[leg] -> qv}],
    {j, 3}, {leg, {Keys[qLeg][[j]]}}]; rules];

(* Hash-consing gives each scalar arithmetic leaf one evaluation per joint draw.
   This is a WL-native circuit representation, independent of any sibling shape. *)
leafExpressions = {}; leafIndex = <||>;
nodeDefinitions = {}; nodeIndices = <||>;
intern[x_] := Module[{key = ToString[Unevaluated[x], InputForm], n},
  If[KeyExistsQ[leafIndex, key], Return[leafIndex[key]]];
  AppendTo[leafExpressions, x]; n = Length[leafExpressions]; AssociateTo[leafIndex, key -> n]; n];
internNode[type_, children_] := Module[{key = ToString[{type, children}, InputForm], index},
  If[KeyExistsQ[nodeIndices, key], Return[branch[nodeIndices[key]]]];
  AppendTo[nodeDefinitions, {type, children}]; index = Length[nodeDefinitions];
  AssociateTo[nodeIndices, key -> index]; branch[index]];
compileCircuit[x_arithmeticPlus] := internNode["PLUS", compileCircuit /@ List @@ x];
compileCircuit[x_arithmeticTimes] := internNode["TIMES", compileCircuit /@ List @@ x];
compileCircuit[arithmeticPower[x_, n_]] := internNode["POWER", {compileCircuit[x], n}];
compileCircuit[x_] := internNode["LEAF", intern[x]];
evaluateNode[branch[n_], p_] := nodeValueCache[[n]];
nodeDegree[branch[n_]] := nodeDegreeCache[[n]];
evaluateCompiledNodes[p_] := Module[{definition, children},
  nodeValueCache = ConstantArray[{0, 1}, Length[nodeDefinitions]];
  Do[definition = nodeDefinitions[[j]]; children = definition[[2]];
    nodeValueCache[[j]] = Switch[First[definition], "LEAF", leafValues[[children]],
      "PLUS", Fold[pairAdd[#1, #2, p] &, {0, 1}, evaluateNode[#, p] & /@ children],
      "TIMES", Fold[pairMultiply[#1, #2, p] &, {1, 1}, evaluateNode[#, p] & /@ children],
      "POWER", pairPower[evaluateNode[children[[1]], p], children[[2]], p]], {j, Length[nodeDefinitions]}]];
degreeCompiledNodes[] := Module[{definition, children},
  nodeDegreeCache = ConstantArray[{0, 0}, Length[nodeDefinitions]];
  Do[definition = nodeDefinitions[[j]]; children = definition[[2]];
    nodeDegreeCache[[j]] = Switch[First[definition], "LEAF", leafDegrees[[children]],
      "PLUS", Fold[{Max[#1[[1]] + #2[[2]], #2[[1]] + #1[[2]]], #1[[2]] + #2[[2]]} &, {0, 0}, nodeDegree /@ children],
      "TIMES", Total[nodeDegree /@ children],
      "POWER", If[children[[2]] >= 0, children[[2]] nodeDegree[children[[1]]],
        -children[[2]] Reverse[nodeDegree[children[[1]]]]]], {j, Length[nodeDefinitions]}]];
evaluateNode[leaf[n_], p_] := leafValues[[n]];
evaluateNode[x_nodePlus, p_] := Fold[pairAdd[#1, #2, p] &, {0, 1}, evaluateNode[#, p] & /@ List @@ x];
evaluateNode[x_nodeTimes, p_] := Fold[pairMultiply[#1, #2, p] &, {1, 1}, evaluateNode[#, p] & /@ List @@ x];
evaluateNode[nodePower[x_, n_], p_] := pairPower[evaluateNode[x, p], n, p];
nodeDegree[leaf[n_]] := leafDegrees[[n]];
nodeDegree[x_nodePlus] := Fold[{Max[#1[[1]] + #2[[2]], #2[[1]] + #1[[2]]], #1[[2]] + #2[[2]]} &,
  {0, 0}, nodeDegree /@ List @@ x];
nodeDegree[x_nodeTimes] := Total[nodeDegree /@ List @@ x];
nodeDegree[nodePower[x_, n_]] := If[n >= 0, n nodeDegree[x], -n Reverse[nodeDegree[x]]];

outputObjects = <||>; outputMetadata = <||>;
objectId[family_, name_] := family <> ":" <> name;
put[family_, name_, case_, value_] := Module[{id = objectId[family, name], cases},
  cases = Lookup[outputObjects, id, <||>]; AssociateTo[cases, case -> value]; AssociateTo[outputObjects, id -> cases]];
putMeta[family_, name_, case_, value_] := Module[{id = objectId[family, name], cases},
  cases = Lookup[outputMetadata, id, <||>]; AssociateTo[cases, case -> value]; AssociateTo[outputMetadata, id -> cases]];
numericObjects = <||>;
withFaceSum[mapping_, axis_] := Module[{out = mapping, partner, sumKey},
  Do[If[key[[axis]] === 1, partner = ReplacePart[key, axis -> -1]; sumKey = ReplacePart[key, axis -> "SUM"];
    AssociateTo[out, sumKey -> MapThread[add, {mapping[key], mapping[partner]}]]], {key, Keys[mapping]}]; out];
addNumeric[family_, name_, mapping_] := Module[{axis},
  axis = Which[family === "GUARD", 3,
    MemberQ[{"CARRIER_EULERIAN", "CARRIER_MATERIAL", "CARRIER_BRIDGE_RESIDUAL"}, name], 2,
    MemberQ[{"SOURCE_EULERIAN", "SOURCE_MATERIAL", "SOURCE_BRIDGE_RESIDUAL", "SOURCE_ACTUAL", "SOURCE_PREDICTED",
      "R_COV", "SOURCE_BASELINE", "R_COV_BASELINE", "SOURCE_CONTROL_DELTA", "R_COV_CONTROL_DELTA"}, name], 1,
    True, 4];
  AssociateTo[numericObjects, objectId[family, name] -> withFaceSum[mapping, axis]]];
sourceMap[ss_Association] := Association[Flatten[Table[
  With[{f = atPoint[finish[ss[s]], "Y"]}, With[{aa = Select[atomsIn[f], MemberQ[waveFields, jetRegistry[#][[1]]] &]},
    Table[{s, wave, {1, Sequence @@ g}} -> {gradePart[
      Total[Map[If[jetRegistry[#][[1]] === wave, D[f, #] #, 0] &, aa]], g]},
      {wave, waveFields}, {g, gradeIndices}]]], {s, faces}], 2]];
carrierMap[cc_Association] := Association[Flatten[Table[
  {If[r <= 3, "U" <> ToString[r], If[r === 4, "THETA", "E_W"]], s, pressureSlots[[p]],
    {1, Sequence @@ g}} -> {gradePart[atPoint[finish[cc[s][[r, p]]], "X"], g]},
  {s, faces}, {r, 5}, {p, 4}, {g, gradeIndices}], 3]];
joinFaceMaps[builder_] := Join @@ (builder /@ faces);
zeroLike[map_] := Association[Map[# -> ConstantArray[0, Length[map[#]]] &, Keys[map]]];
allFamilies[map_] := Module[{out = map},
  Do[If[!KeyExistsQ[out, {pair[[1]] <> "_FROM_" <> pair[[2]], "LOCAL_BARE", {1, Sequence @@ g}, s}],
    AssociateTo[out, {pair[[1]] <> "_FROM_" <> pair[[2]], "LOCAL_BARE", {1, Sequence @@ g}, s} ->
      ConstantArray[0, If[First[pair] === "TRANSVERSE", 3, 1]]]], {pair, blocks}, {g, gradeIndices}, {s, faces}]; out];

probeCase[case_] := Module[{compiled, flatNodes, degreesByCell, primeList, seeds = {}, tables, metadata,
    pointRows = {}, rejectedRows = {}, bounds = {}, variables, sampleAtoms, chart, chartAtoms, atomStreams,
    d, excluded, numeratorDegree, seed, attempts, valid, samples, rejected, row, nodeValues,
    allValues, prime, cell, point, chartPairs, bad, drawCount = 8, bound, unionCount, metadataCases,
    probePositions, positionCursor = 0, positions},
  leafExpressions = {}; leafIndex = <||>; nodeDefinitions = {}; nodeIndices = <||>;
  compiled = Association[KeyValueMap[Function[{name, entries}, name -> Association[
    KeyValueMap[Function[{key, values}, key -> (compileCircuit /@ values)], entries]]], numericObjects]];
  flatNodes = Flatten[Values /@ Values[compiled]];
  unionCount = Length[flatNodes];
  primeList = {1000000009, 998244353, 1004535809};
  If[!And @@ (PrimeQ /@ primeList), Quit[96]];
  (* Store one joint draw once. Cell positions follow the identical ordered
     flattening used above; emission slices these rows without rebuilding a
     growing nested Association separately for every cell of every draw. *)
  probePositions = Association[KeyValueMap[Function[{name, entries}, name -> Association[
    KeyValueMap[Function[{key, values},
      positions = Range[positionCursor + 1, positionCursor + Length[values]];
      positionCursor += Length[values]; key -> positions], entries]]], compiled]];
  tables = {};
  variables = DeleteDuplicates[Cases[leafExpressions, _Symbol, Infinity]];
  variables = Select[variables, Context[#] === "Global`" &];
  Do[cell = branchCells[[cellIndex]]; chart = chartRules[cell];
    parameterDegrees = <||>;
    Do[AssociateTo[parameterDegrees, First[rule] -> degreePair[Last[rule]]], {rule, chart}];
    leafDegrees = degreePair /@ leafExpressions;
    degreeCompiledNodes[];
    degreesByCell = nodeDegree /@ flatNodes;
    (* Joint singular rejection includes all leaf denominators and chart denominators.
       Product-degree bound is deliberately conservative and is derived from this circuit. *)
    excluded = Total[Last /@ leafDegrees] + Total[Last[degreePair[Last[#]]] & /@ chart] +
      First[degreePair[(momentum["kOut"][[1]] - momentum["kIn"][[1]]) /. chart]] +
      Total[First[degreePair[# /. chart]] & /@ Values[qLeg]] +
      Total[Cases[nodeDefinitions, {"POWER", {base_, n_Integer /; n < 0}} :> First[nodeDegree[base]]]];
    numeratorDegree = Max[Prepend[First /@ degreesByCell, 0]];
    bound = If[excluded < Min[primeList] - 1,
      Min[1, numeratorDegree/(Min[primeList] - 1 - excluded)], 1];
    drawCount = If[0 < bound < 1, Max[8, Ceiling[Log[2^-80/Max[1, unionCount]]/Log[bound]]], 8];
    Do[prime = pp; imaginaryRoot = PowerMod[PrimitiveRoot[prime], (prime - 1)/4, prime];
      If[Mod[imaginaryRoot^2 + 1, prime] =!= 0, Quit[97]];
      seed = 731921 + 10000 caseOrdinal + 100 cellIndex + First[FirstPosition[primeList, prime]];
      SeedRandom[seed, Method -> "MersenneTwister"]; AppendTo[seeds, {cell, prime, seed}];
      samples = {}; rejected = 0; valid = 0;
      sampleAtoms = Complement[Union[variables, chartNames, {frequencyScale, cS0}], First /@ chart];
      (* Named independent streams keep an unaffected atom's draws fixed when a
         one-sided control introduces a previously absent dependency. *)
      atomStreams = Association[Map[Function[atom, atom -> BlockRandom[
        SeedRandom[Hash[{seed, SymbolName[atom]}, "SHA256"], Method -> "MersenneTwister"];
        RandomInteger[{1, prime - 1}, 4 drawCount + 100]]], sampleAtoms]];
      attempts = 0;
      (* The loop conditions are sample validity/count only, never residual values. *)
      While[valid < drawCount && attempts < 4 drawCount + 100,
        attempts++;
        point = atomStreams[#][[attempts]] & /@ sampleAtoms;
        sampleValues = Association[Thread[sampleAtoms -> point]];
        chartPairs = numericPair[Last[#], prime] & /@ chart;
        bad = AnyTrue[chartPairs, Last[#] === 0 &];
        If[!bad,
          MapThread[AssociateTo[sampleValues, First[#1] -> Mod[#2[[1]] PowerMod[#2[[2]], -1, prime], prime]] &,
            {chart, chartPairs}];
          bad = Lookup[sampleValues, momentum["kOut"]] === Lookup[sampleValues, momentum["kIn"]] ||
            MemberQ[Lookup[sampleValues, Values[qLeg]], 0]];
        If[!bad, leafValues = numericPair[#, prime] & /@ leafExpressions;
          bad = AnyTrue[leafValues, Last[#] === 0 || !VectorQ[#, IntegerQ] &]];
        If[!bad, evaluateCompiledNodes[prime]; bad = AnyTrue[nodeValueCache, Last[#] === 0 &]];
        If[bad, rejected++; Continue[]];
        valid++; AppendTo[samples, point];
        nodeValues = evaluateNode[#, prime] & /@ flatNodes;
        AppendTo[tables, {cellIndex, prime, valid, First /@ nodeValues, Last /@ nodeValues}]];
      If[valid < drawCount, WriteString[$Messages[[1]], "Joint sampler exhausted at case ", ToString[case], "\n"]; Quit[98]];
      AppendTo[pointRows, <|"CELL" -> cell, "PRIME" -> prime, "ATOMS" -> sampleAtoms,
        "ATOM_SEEDS" -> (Hash[{seed, SymbolName[#]}, "SHA256"] & /@ sampleAtoms), "DRAWS" -> samples|>];
      AppendTo[rejectedRows, {cell, prime, rejected}];
      bound = If[excluded < prime - 1, Min[1, numeratorDegree/(prime - 1 - excluded)], Missing["ExcludedDegree"]];
      AppendTo[bounds, <|"CELL" -> cell, "PRIME" -> prime, "D" -> numeratorDegree,
        "E" -> excluded, "N" -> prime - 1, "PER_VALID_DRAW" -> bound,
        "DRAWS" -> valid, "FAMILY_BOUND" -> If[NumberQ[bound], Min[1, unionCount bound^valid], bound]|>];
      appendAssociationEmission[{case, "PROBE", cellIndex, prime}, {valid, rejected, Length[leafExpressions]}],
      {pp, primeList}], {cellIndex, Length[branchCells]}];
  KeyValueMap[Function[{name, entries}, Module[{parts = StringSplit[name, ":"], payload},
    payload = Association[KeyValueMap[Function[{key, values}, key -> <|
      "ARITHMETIC" -> (circuitExpression /@ values), "DIMENSIONS" -> (keyDimensionRecord[#, key] & /@ values),
      "COMPONENT_AXES" -> If[name === "GUARD:SLOT_GUARD_RESIDUAL",
        Join[{"G_LINEAR"}, Flatten[Table[{"G_CROSS", p, q}, {p, pressureSlots}, {q, pressureSlots}], 1],
          Table[{"DENOMINATOR_PRESSURE_DEPENDENCY", p}, {p, pressureSlots}]], Range[Length[values]]],
      "PROBE_NUMERATORS" -> SparseArray[tables[[All, 4, probePositions[name][key]]]],
      "PROBE_DENOMINATORS" -> SparseArray[tables[[All, 5, probePositions[name][key]]], Automatic, 1],
      "SAMPLE_INDEX" -> {case, HoldForm[jointCaseSampleIndex]}|>], entries]];
    put[parts[[1]], parts[[2]], case, payload]]], numericObjects];
  putMeta["LOCAL", "PROBE", case, <|"PRIMES" -> Map[# -> Inactive[PrimeQ][#] &, primeList],
    "CELLS" -> branchCells, "SAMPLE_POINTS" -> pointRows, "SEEDS" -> seeds, "REJECTIONS" -> rejectedRows,
    "SAMPLE_INDEX" -> tables[[All, 1 ;; 3]],
    "BOUNDS" -> bounds, "FAMILY_CARDINALITY" -> unionCount,
    "CONDITIONAL_UNION_BOUND" -> Min[1, Total[Table[Max[Lookup[Select[bounds, #["CELL"] === bc &], "FAMILY_BOUND"]], {bc, branchCells}]]],
    "BAD_PRIME_CONDITION" -> Inactive[Exists][goodPrime, Inactive[And][
      Inactive[Element][goodPrime, primeList], Inactive[Unequal][reducedNonzeroCoefficient[goodPrime], 0]]],
    "SCOPE" -> HoldForm[formalWeakKernelRationalCoefficients],
    "ALL_ZERO_TABLE_SEMANTICS" -> Inactive[ConditionalExpression][noNonzeroFound, conditionalFalseNegativeBound],
    "EXCLUDED_LOCI" -> {Inactive[Equal][qOut, 0], Inactive[Equal][clearedDenominator, 0],
      Inactive[Equal][kOut, kIn]}, "COMPLEX_FREQUENCY_COVERAGE" -> Missing["RealAxisPITOnly"]|>];
  Clear[sampleValues, leafValues, parameterDegrees, leafDegrees]; ClearSystemCache[]];

buildCase[case_List] := Module[{anchor = case[[1]], densityKind = case[[2]], density4, density3,
    a, h, energy, muE, mat, baseline, muM, predictedMap, muPred, geometriesE = <||>,
    geometriesM = <||>, geometriesSourceM = <||>, lawE = <||>, lawM = <||>, rowsE = <||>, rowsM = <||>,
    ce = <||>, cm = <||>, es = <||>, ms = <||>, predicted = <||>, baselineSource = <||>,
    velocitiesE = <||>, velocitiesM = <||>, sourceSolve, kernel, response, massBase, eMap, mMap,
    sEMap, sMMap, pMap, baselineMap, bridge, sourceDelta, covDelta, baselineCov,
    eOperand, mOperand, rawResidual, carrierChannel, sourceChannel, crossChannel, splitSum, splitCheck,
    covIncrement, native = <||>, reconstructed = <||>, slotResidual = <||>, closedNative = <||>,
    closedCarrier = <||>, closedResidual = <||>, carriers, rows, reconstructedRow, projectedRow,
    closedRules, directDifference, closureDifference, key, entries, census, fingerprint, materialTagged,
    sourceControl, rCovControl, dimObject, phiFieldMap, pressureCensus, stages,
    imageGrades, imageRules, openGrades, closedGrades, carrierGrades, basePressureGrades,
    familyName, sourceForGuard, pressureImage, virtualMass, virtualConstraint, massData, activeJunk},
  numericObjects = <||>;
  density4 = If[densityKind === "RHO4_CONSTANT", rhoBr/W0, rhoBr/WBg];
  density3 = density4 WBg;
  a = displacement.Table[td[density4, i]/density4, {i, 3}];
  h = If[anchor === "LAB_HELD", displacement.Table[td[WBg, i]/WBg, {i, 3}], 0];
  energy = constructEnergy[]; muE = el[energy["DENSITY"]];
  appendAssociationEmission[{case, "ENERGY"}, {Length[energy["GENERATED_CONTRACTIONS"]], LeafCount[muE]}];
  activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0];
  mat = materialAmplitude[energy["DENSITY"], a, h, actualAdvectionCoefficient,
    actualThicknessCoefficient, activeJunk]; muM = mat["MU"];
  baseline = materialAmplitude[energy["DENSITY"], a, h, 1, 1, 0];
  predictedMap = phiMap[muE, predictionAdvectionCoefficient a, predictionThicknessCoefficient h];
  putMeta["COV", "PHI_DOMAIN_CENSUS", case, KeyDrop[predictedMap, "MAP"]];
  If[predictedMap["UNCOVERED"] =!= {},
    emit[{"COV", "PHI_DOMAIN_CENSUS"}, outputMetadata["COV:PHI_DOMAIN_CENSUS"]]; Quit[92]];
  muPred = muE /. predictedMap["MAP"];
  appendAssociationEmission[{case, "AMPLITUDES"}, LeafCount /@ {muE, muM, muPred}];
  sourceSolve = sourceConstruction[];
  kernel = Association[Table[s -> constructKernel[s], {s, faces}]];
  response = Map[responseFamilies, kernel];
  massData = massSubstrate[anchor, density4];
  Do[AssociateTo[geometriesE, s -> graphGeometry[anchor, s]];
    AssociateTo[geometriesM, s -> materialGeometry[anchor, s, materialNormalKnife]];
    AssociateTo[geometriesSourceM, s -> materialGeometry[anchor, s, 0]];
    AssociateTo[velocitiesE, s -> geometriesE[s]["VELOCITY"]];
    AssociateTo[velocitiesM, s -> geometriesSourceM[s]["VELOCITY"]];
    AssociateTo[lawE, s -> faceLaws[geometriesE[s], s, muE, density3]];
    AssociateTo[lawM, s -> faceLaws[geometriesM[s], s, muM, density3]];
    AssociateTo[rowsE, s -> eulerianSlabFace[lawE[s], massData["EULERIAN_EVOLUTION_BASE"]/2]];
    AssociateTo[rowsM, s -> materialFaceFold[lawM[s], massData["MATERIAL_EVOLUTION_BASE"]/2]];
    AssociateTo[ce, s -> carrier[rowsE[s]]]; AssociateTo[cm, s -> carrier[rowsM[s]]];
    AssociateTo[es, s -> sourceBind[sourceSolve, muE, velocitiesE[s], density3]];
    AssociateTo[ms, s -> sourceBind[sourceSolve, muM, velocitiesM[s], density3]];
    AssociateTo[predicted, s -> sourceBind[sourceSolve, muPred, velocitiesE[s], density3]];
    AssociateTo[baselineSource, s -> sourceBind[sourceSolve, baseline["MU"], velocitiesM[s], density3]];
    appendAssociationEmission[{case, "FACES", s}, {LeafCount[ce[s]], LeafCount[cm[s]]}], {s, faces}];
  eMap = carrierMap[ce]; mMap = carrierMap[cm]; bridge = mapCombine[eMap, mMap, sub];
  addNumeric["RC", "CARRIER_EULERIAN", eMap]; addNumeric["RC", "CARRIER_MATERIAL", mMap];
  addNumeric["RC", "CARRIER_BRIDGE_RESIDUAL", bridge];
  sEMap = sourceMap[es]; sMMap = sourceMap[ms]; pMap = sourceMap[predicted]; baselineMap = sourceMap[baselineSource];
  sourceDelta = mapCombine[sEMap, sMMap, sub]; covDelta = mapCombine[sMMap, pMap, sub];
  baselineCov = mapCombine[baselineMap, pMap, sub];
  addNumeric["RC", "SOURCE_EULERIAN", sEMap]; addNumeric["RC", "SOURCE_MATERIAL", sMMap];
  addNumeric["RC", "SOURCE_BRIDGE_RESIDUAL", sourceDelta];
  addNumeric["COV", "SOURCE_ACTUAL", sMMap]; addNumeric["COV", "SOURCE_PREDICTED", pMap];
  addNumeric["COV", "R_COV", covDelta]; addNumeric["COV", "SOURCE_BASELINE", baselineMap];
  addNumeric["COV", "R_COV_BASELINE", baselineCov];
  addNumeric["COV", "SOURCE_CONTROL_DELTA", mapCombine[sMMap, baselineMap, sub]];
  addNumeric["COV", "R_COV_CONTROL_DELTA", mapCombine[covDelta, baselineCov, sub]];
  appendAssociationEmission[{case, "SOURCES"}, {Length[sEMap], Length[sMMap]}];
  eOperand = joinFaceMaps[buildContraction[ce[#], es[#], response[#], #, True] &];
  mOperand = joinFaceMaps[buildContraction[cm[#], ms[#], response[#], #, True] &];
  rawResidual = mapCombine[eOperand, mOperand, sub];
  carrierChannel = joinFaceMaps[buildContraction[ce[#] - cm[#], ms[#], response[#], #, True] &];
  sourceChannel = joinFaceMaps[buildContraction[cm[#], es[#] - ms[#], response[#], #, False] &];
  crossChannel = joinFaceMaps[buildContraction[ce[#] - cm[#], es[#] - ms[#], response[#], #, False] &];
  splitSum = mapCombine[mapCombine[carrierChannel, allFamilies[sourceChannel], add], allFamilies[crossChannel], add];
  splitCheck = mapCombine[splitSum, rawResidual, sub];
  covIncrement = joinFaceMaps[buildContraction[cm[#], ms[#] - predicted[#], response[#], #, False] &];
  MapThread[addNumeric["RC", #1, #2] &, {{"EULERIAN_OPERAND", "MATERIAL_OPERAND", "R_N6",
    "CARRIER_CHANNEL", "SOURCE_CHANNEL", "CROSS_CHANNEL", "SPLIT_SUM", "SPLIT_CHECK"},
    {eOperand, mOperand, rawResidual, carrierChannel, sourceChannel, crossChannel, splitSum, splitCheck}}];
  addNumeric["COV", "R_COV_INCREMENT", covIncrement];
  appendAssociationEmission[{case, "WEAK_CONTRACTIONS"}, {Length[eOperand], Length[mOperand]}];
  (* The full native circuit is evaluated at the same derived closed-response
     images used in the increment. Each formal integral family is kept separate.
     Graded circuit substitution avoids constructing a second symbolic slab. *)
  Do[rows = If[route === "EULERIAN", rowsE[s], rowsM[s]];
    carriers = If[route === "EULERIAN", ce[s], cm[s]];
    Do[projectedRow = atPoint[finish[rows[[r]] - (rows[[r]] /. slotZero)], "X"];
      reconstructedRow = atPoint[finish[carriers[[r]].pressureSlots], "X"];
      Do[key = {route, r, s, {1, Sequence @@ g}};
        AssociateTo[native, key -> {gradePart[projectedRow, g]}];
        AssociateTo[reconstructed, key -> {gradePart[reconstructedRow, g]}];
        AssociateTo[slotResidual, key -> Join[{sub[gradePart[projectedRow, g], gradePart[reconstructedRow, g]]},
          Flatten[Table[gradePart[atPoint[finish[D[rows[[r]], p, q]], "X"], g], {p, pressureSlots}, {q, pressureSlots}]],
          Table[D[denominatorCircuit[atPoint[finish[rows[[r]]], "X"]], p], {p, pressureSlots}]]],
        {g, gradeIndices}];
      Do[imageRules = <||>;
        Do[sourceForGuard = atPoint[finish[If[route === "EULERIAN", es[face], ms[face]]], "Y"];
          pressureImage = gMul[response[face][familyName], graded[sourceForGuard]];
          AssociateTo[imageRules, pressureSlots[[If[face === 1, 1, 3]]] -> pressureImage];
          AssociateTo[imageRules, pressureSlots[[If[face === 1, 2, 4]]] -> gScale[pressureImage, normalContinuation[face]]],
          {face, faces}];
        openGrades = gradeCircuit[atPoint[finish[rows[[r]]], "X"], <||>];
        closedGrades = gradeCircuit[atPoint[finish[rows[[r]]], "X"], imageRules];
        carrierGrades = Fold[gAdd, {0, 0, 0, 0}, Table[gMul[graded[atPoint[finish[carriers[[r, p]]], "X"]],
          gSub[imageRules[pressureSlots[[p]]], {pressureSlots[[p]], 0, 0, 0}]], {p, 4}]];
        Do[key = {route, r, s, familyName, {1, Sequence @@ gradeIndices[[g]]}};
          AssociateTo[closedNative, key -> {sub[closedGrades[[g]], openGrades[[g]]]}];
          AssociateTo[closedCarrier, key -> {carrierGrades[[g]]}];
          AssociateTo[closedResidual, key -> {sub[First[closedNative[key]], First[closedCarrier[key]]]}], {g, 4}],
        {familyName, Rest[kernelFamilies]}], {r, 5}], {route, {"EULERIAN", "MATERIAL"}}, {s, faces}];
  MapThread[addNumeric["GUARD", #1, #2] &, {guardNames,
    {native, reconstructed, slotResidual, closedNative, closedCarrier, closedResidual}}];
  pressureCensus = Association[Table[pressureSlots[[i]] -> <|"INHERITED" -> pressureSlots[[i]],
    "EULERIAN_BINDING" -> Flatten[Lookup[Values[lawE], "PRESSURE_IDENTITIES"]][[i]],
    "MATERIAL_BINDING" -> Flatten[Lookup[Values[lawM], "PRESSURE_IDENTITIES"]][[i]],
    "ASSUMPTION" -> With[{assumption = pressureAssumptions[[i]]}, HoldForm[assumption]],
    "IDENTITY_COMPARISONS" -> {Inactive[SameQ][pressureSlots[[i]], Flatten[Lookup[Values[lawE], "PRESSURE_IDENTITIES"]][[i]]],
      Inactive[SameQ][pressureSlots[[i]], Flatten[Lookup[Values[lawM], "PRESSURE_IDENTITIES"]][[i]]]}|>, {i, 4}]];
  putMeta["RC", "FROZEN_RELATIONS", case, <|"DENSITY4" -> density4, "DENSITY3" -> density3,
    "LIVE_DENSITY_NAME" -> If[densityKind === "RHO4_CONSTANT", rhoBrBgRho4Constant, rhoBr],
    "A_RHO" -> a, "H_ALPHA" -> h, "DENSITY_GRADIENT" -> Table[td[density4, i], {i, 3}],
    "DENSITY_JACOBIAN" -> 1 + Sum[jet["u" <> ToString[i], {i}], {i, 3}],
    "FIELD_MAP" -> {theta -> theta + a, eW -> eW + h}, "PROLONGATION" -> predictedMap["MAP"],
    "MATERIAL_COVECTOR_MAP" -> Lookup[Values[geometriesM], "INVERSE_TRANSPOSE"],
    "JET_VOCABULARY_BRIDGE" -> Table[<|"FACE_ATOM" -> Symbol["thetaD" <> ToString[i]],
      "ENERGY_ATOM" -> Symbol["gradTheta" <> ToString[i]], "PHYSICAL_JET" -> jet["theta", {i}],
      "IDENTITIES" -> {Symbol["thetaD" <> ToString[i]] -> jet["theta", {i}],
        Symbol["gradTheta" <> ToString[i]] -> jet["theta", {i}]}|>, {i, 3}],
    "PRESSURE_CENSUS" -> pressureCensus, "GRADES" -> gradeIndices|>];
  stages = <|"MU_E" -> {"EL", energy["DENSITY"], muE},
    "MU_M" -> {"PULLBACK_EL", mat["DENSITY"], muM},
    "V_E" -> {"EULERIAN_LEVEL_SET", geometriesE, velocitiesE},
    "V_M" -> {"MATERIAL_FLATTENING", geometriesSourceM, velocitiesM},
    "C_E" -> {"EULERIAN_SLAB_ROWS", rowsE, ce}, "C_M" -> {"MATERIAL_FACE_FOLD", rowsM, cm},
    "ES" -> {"SOURCE_SOLVE", {muE, velocitiesE, sourceSolve}, es},
    "MS" -> {"SOURCE_SOLVE", {muM, velocitiesM, sourceSolve}, ms}|>;
  putMeta["RC", "PROVENANCE", case, <|"STAGES" -> Association[KeyValueMap[Function[{name, data},
    name -> <|"BUILDER" -> First[data], "INPUT_FINGERPRINT" -> Hash[data[[2]], "SHA256"],
      "OUTPUT_FINGERPRINT" -> Hash[data[[3]], "SHA256"]|>], stages]],
    "ENERGY" -> energy, "CONSTITUTIVE_OPERANDS" -> {finish[muE], finish[muM]},
    "FACE_SUBSTRATE_E" -> lawE, "FACE_SUBSTRATE_M" -> lawM,
    "MASS_AND_CONSTRAINT_SUBSTRATE" -> massData,
    "CENTER_WORK_OPERANDS" -> {Table[-D[lawE[s]["VIRTUAL_WORK"], virtualCenter], {s, faces}],
      Table[-D[lawM[s]["VIRTUAL_WORK"], virtualCenter], {s, faces}]},
    "LOCAL_BARE_ROWS_BEFORE_WEAK_RESTRICTION" -> {Map[-#.pressureSlots &, ce], Map[-#.pressureSlots &, cm]},
    "VELOCITY_OPERAND_A" -> Map[finish, velocitiesE], "VELOCITY_OPERAND_B" -> Map[finish, velocitiesM],
    "VELOCITY_DIFFERENCE" -> MapThread[sub, {Values[Map[finish, velocitiesE]], Values[Map[finish, velocitiesM]]}],
    "C1_SOURCE_SOLVE" -> sourceSolve, "C1_KERNEL" -> kernel, "C1_RESPONSE_FAMILIES" -> response,
    "BULK_PRESSURE_DEPENDENCY" -> Table[D[energy["DENSITY"], p], {p, pressureSlots}],
    "INDEPENDENT_BUILDERS" -> KeyTake[stages, {"MU_E", "MU_M", "V_E", "V_M", "C_E", "C_M"}]|>];
  If[densityKind === "RHO4_CONSTANT",
    materialTagged = materialAmplitude[energy["DENSITY"], a, h, advectionTag, 1, 0];
    putMeta["RC", "ADVECTION_ABSENCE", case, <|"DENSITY_GRADIENT" -> Table[td[density4, i], {i, 3}],
      "ADVECTION" -> a, "MATERIAL_MU_TAG_DERIVATIVE" -> finish[D[materialTagged["MU"], advectionTag]]|>]];
  putMeta["COV", "FROZEN_PHI", case, <|"MAP" -> predictedMap["MAP"], "A_RHO" -> a, "H_ALPHA" -> h,
    "PREDICTION_PARAMETERS" -> {predictionAdvectionCoefficient, predictionThicknessCoefficient}|>];
  putMeta["COV", "ACTUAL_CONTROL_PARAMETERS", case, <|"ADVECTION" -> actualAdvectionCoefficient,
    "THICKNESS" -> actualThicknessCoefficient, "JUNK" -> activeJunk, "JUNK_CASE" -> actualJunkCase,
    "JUNK_SYMBOL" -> junkMu, "JUNK_DIMENSIONS" -> unitSupport[junkMu], "JUNK_ASSUMPTION" -> Inactive[Unequal][junkMu, 0],
    "MATERIAL_NORMAL" -> materialNormalKnife|>];
  putMeta["COV", "PROVENANCE", case, <|"ACTUAL_INPUT" -> Hash[mat["DENSITY"], "SHA256"],
    "PREDICTION_INPUT" -> Hash[{muE, predictedMap["MAP"], velocitiesE}, "SHA256"],
    "PREDICTION_BUILDER" -> {HoldForm[el], HoldForm[phiMap], "EULERIAN_VELOCITY"},
    "PREDICTED_MU" -> finish[muPred], "ACTUAL_MU" -> finish[muM],
    "THETA_INDEPENDENT_JUNK" -> activeJunk junkMu eW|>];
  dimObject = Association[KeyValueMap[Function[{name, entries}, name -> Association[
    KeyValueMap[Function[{key, values}, key -> (keyDimensionRecord[#, key] & /@ values)], entries]]], numericObjects]];
  putMeta["RC", "DIMENSIONS", case, <|"OBJECTS" -> dimObject,
    "METADATA_OBJECTS" -> Association[Map[Function[id, id -> dimensionTree[outputMetadata[id][case]]],
      Select[Keys[outputMetadata], # =!= "RC:DIMENSIONS" && KeyExistsQ[outputMetadata[#], case] &]]],
    "BASE_OPERAND" -> deltaPPlus, "CONTROL_OPERAND" -> deltaPPlus + W0 deltaPMinus,
    "EXPRESSION_DIFFERENCE" -> sub[deltaPPlus, deltaPPlus + W0 deltaPMinus],
    "DIMENSION_OPERANDS" -> {unitSupport[deltaPPlus], unitSupport[deltaPPlus + W0 deltaPMinus]},
    "CONTROL_ADDITION_OPERANDS" -> {unitSupport[deltaPPlus], unitSupport[W0 deltaPMinus]},
    "CONTROL_DIMENSION_DIFFERENCE" -> First[unitSupport[deltaPPlus]] - First[unitSupport[W0 deltaPMinus]],
    "CONTROL_ADDITION_CONSISTENCY" -> relationalObject[Length[unitSupport[deltaPPlus + W0 deltaPMinus]], 1],
    "DIMENSION_COMPARISON" -> relationalObject[unitSupport[deltaPPlus], unitSupport[deltaPPlus + W0 deltaPMinus]]|>];
  appendAssociationEmission[{case, "PROBE_INPUT"}, Length[numericObjects]];
  probeCase[case]; ClearSystemCache[]];

caseOrdinal = 0; constructionEventIndex = 0;
clearCaseMemo[] := (
  DownValues[td] = Select[DownValues[td], !FreeQ[First[#], _Pattern] &];
  DownValues[finish] = Select[DownValues[finish], !FreeQ[First[#], _Pattern] &];
  DownValues[retainTerm] = Select[DownValues[retainTerm], !FreeQ[First[#], _Pattern] &];
  DownValues[restrict] = Select[DownValues[restrict], !FreeQ[First[#], _Pattern] &];
  DownValues[testCarrier] = Select[DownValues[testCarrier], !FreeQ[First[#], _Pattern] &];
  DownValues[pruneBackgroundProducts] = Select[DownValues[pruneBackgroundProducts], !FreeQ[First[#], _Pattern] &];
  Clear[leafValues, sampleValues]; ClearSystemCache[]);
beginAssociationEmission[{"LOCAL", "CONSTRUCTION"}];
Do[caseOrdinal++; Check[buildCase[case], Quit[99]];
  appendAssociationEmission[{case, "STAGES"},
    <|"CASE" -> case, "OBJECT_KEYS" -> Keys[numericObjects], "CIRCUIT_LEAVES" -> Length[leafExpressions]|>]; clearCaseMemo[],
  {case, Tuples[{{"LAB_HELD", "MATERIAL_ADVECTED"}, {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}]}];
endAssociationEmission[];
(* Object order is explicit: operands precede differences, and guards follow them. *)
Do[With[{id = objectId[family, name]},
  If[KeyExistsQ[outputObjects, id], beginAssociationEmission[{family, name}];
    appendAssociationEmissionFromSpool[Normal[outputObjects[id]]]; endAssociationEmission[],
    If[KeyExistsQ[outputMetadata, id], emit[{family, name}, outputMetadata[id]]]]],
  {family, {"RC", "COV", "GUARD"}},
  {name, Switch[family, "RC", DeleteCases[rcNames, "DIMENSIONS"], "COV", covNames, "GUARD", guardNames]}];
(* Dimensional guard operands follow the numerical comparisons in every family. *)
emit[{"RC", "DIMENSIONS"}, outputMetadata["RC:DIMENSIONS"]];
emit[{"LOCAL", "PROBE"}, outputMetadata["LOCAL:PROBE"]];
emit[{"LOCAL", "INVENTORY"}, <|"TAGS" -> Append[emittedNames, "WL_S11CC2_N6_LOCAL_INVENTORY"],
  "SYMBOLS" -> Names["Global`*"], "DIMENSION_REGISTRY" -> unitRegistry|>];
Quit[0];
