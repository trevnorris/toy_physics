(* Focused native-source audit. It does not run or modify any production driver.
   Supplied inputs: kinetic action, reference thickness coordinate, external
   traction virtual work, and a stored-energy probe. All native definitions and
   assembly operands are extracted from the pinned b source before evaluation. *)
$HistoryLength = 0;
ClearAll["Global`*"];
$Messages = {OutputStream["stderr", 2]};
auditRoot = Environment["S11C_WOLFRAM_AUDIT_ROOT"];
auditPath = FileNameJoin[{auditRoot, "sources", "mathematica", "S11c_b_brane_operator_mathematica_audit.wl"}];
auditHeld = ToExpression[Import[auditPath, "Text"], InputForm, HoldComplete];
auditKeys = {};
auditRecords = <||>;
auditEmit[key_String, value_, dimension_, waveOrder_:0] := Module[{record},
  record = <|"value" -> value, "dimensionLTM" -> dimension,
    "epsilonOrder" -> waveOrder, "backgroundDomain" -> "arbitrary stationary smooth width; W != 0; rho != 0",
    "lambda" -> "eta homotopy; sigma independent"|>;
  WriteString[$Output, key, " = ", ToString[record, InputForm, PageWidth -> Infinity], "\n"];
  If[MemberQ[auditKeys, key], Quit[91]]; AppendTo[auditKeys, key];
  AssociateTo[auditRecords, key -> record]];
auditTake[held_, pattern_, label_, count_:1, levels_:{2}] := Module[{found},
  found = Cases[held, x:pattern :> HoldComplete[x], levels];
  auditEmit["wlRepairSource" <> label, <|"count" -> Length[found],
    "sha256" -> Hash[found, "SHA256", "HexString"], "definitions" -> found|>, {0, 0, 0}];
  If[Length[found] =!= count, Quit[92]]; found];
auditLoad[pattern_, label_, count_:1] := Scan[ReleaseHold, auditTake[auditHeld, pattern, label, count]];
auditRhs[held_, pattern_, label_] := First[auditTake[held, pattern, label, 1, Infinity]] /.
  HoldComplete[Set[_, value_]] :> HoldComplete[value];

(* The held assembly is a coefficient audit, not an execution of the full N15
   energy census. Both native stored-energy variations are separately probed. *)
auditConstructors = Map[Function[entry, Module[{definition},
  definition = First[auditTake[auditHeld, entry[[2]], entry[[1]]]];
  <|"name" -> entry[[1]],
    "kineticU" -> auditRhs[definition, HoldPattern[Set[kineticU | kineticULive, _]], entry[[1]] <> "KineticU"],
    "kineticEw" -> auditRhs[definition, HoldPattern[Set[kineticEw | kineticEwLive, _]], entry[[1]] <> "KineticEw"],
    "operator" -> auditRhs[definition, HoldPattern[Set[operator | operatorLive, _Association]], entry[[1]] <> "Assembly"]|>]],
  {{"Raw", HoldPattern[SetDelayed[rawModel[___], _]]},
   {"Frozen", HoldPattern[SetDelayed[frozenEvaluatedModel[___], _]]},
   {"Live", HoldPattern[SetDelayed[evaluatedModel[___], _]]}}];

auditLoad[HoldPattern[Set[braneDimension | spatialCoordinates | materialCoordinates | directions | faces | branches | densities, _]], "Coordinates", 7];
auditLoad[HoldPattern[Set[uVector | virtualUVector | thetaField | eWField | zetaCenterField | virtualThetaField | virtualEwField | virtualZetaCenterField, _]], "Fields", 8];
auditLoad[HoldPattern[SetDelayed[(gradient | divergence | inactiveDivergence | variationalSource | linearVirtualVariation | jacobianLinear | sigmaBackground | virtualConstraintSource | faceSources | faceGeneralizedRows | constrainedRows | constrainedRowsWithLiveEnergyEL)[___], _]], "GeometryVariation", 12];
auditLoad[HoldPattern[SetDelayed[(rhoFourBackground | rhoBrBackground | pressureField)[___], _]], "DensityPressure", 6];
auditLoad[HoldPattern[SetDelayed[activateSpatialDivergences[___], _]], "Activation", 4];
auditEmit["wlRepairMechanicalDomain", <|"routes" -> {"EULERIAN", "MATERIAL"},
  "anchorings" -> branches, "densities" -> densities, "faces" -> faces,
  "constructors" -> auditConstructors[[All, "name"]],
  "energyProbe" -> "local quadratic theta/eW plus all three displacement gradients; not full N15 regeneration",
  "excluded" -> {"spectral coverage", "frequency poles", "full c2 self energy"}|>, {0, 0, 0}];
auditEmit["wlRepairOperandOrder", <|
  "inertia" -> {"native", "independent action", "residual", "negative kinetic action residual", "background coordinate action residual"},
  "storedProbe" -> {"native", "independent first variation", "residual"},
  "faceRow" -> {"native stored-row load", "physical generalized force", "sum residual"},
  "faceWorkAndPower" -> {"native or generalized operand", "independent geometry operand", "difference residual"},
  "primitiveUnits" -> <|"muW and rhoBr" -> {-3, 0, 1}, "WZero" -> {1, 0, 0},
    "probe B C K G" -> {-1, -2, 1}, "auditChemical" -> {-1, -2, 1},
    "lambdaXResponse" -> {-4, 0, 1}, "rhoM" -> {-4, 0, 1},
    "auditWidthJet[i,j,k]" -> "[-i-j-k,0,0]", "auditWidthContrast" -> {0, 0, 0}|>|>, {0, 0, 0}];

auditFields = Join[uVector, {eWField}];
auditVirtuals = Join[virtualUVector, {virtualEwField, virtualZetaCenterField}];
auditRowUnits = {{-2, -2, 1}, {-2, -2, 1}, {-2, -2, 1}, {-1, -2, 1}, {-2, -2, 1}};
auditWidth = widthBase[Sequence @@ spatialCoordinates];
auditProfileRules = {
  HoldPattern[Derivative[i_, j_, k_][widthBase][___]] :> sigmaW WZero auditWidthJet[i, j, k],
  HoldPattern[widthBase[___]] :> WZero (1 + etaBg auditWidthContrast)};
auditPhysical[key_, expression_, unit_, epsilonOrder_] := Module[{raw, retained, grades},
  raw = Simplify[expression] /. auditProfileRules;
  retained = Expand[Normal[Series[Normal[Series[raw, {etaBg, 0, 1}]], {sigmaW, 0, 1}]]];
  grades = Select[Tuples[{Range[0, 1], Range[0, 1]}],
    AnyTrue[Flatten[{Coefficient[Coefficient[retained, etaBg, #[[1]]], sigmaW, #[[2]]]}], !TrueQ[# === 0] &] &];
  auditEmit[key, <|"raw" -> raw, "retained" -> retained,
    "discarded" -> Together[raw - retained], "etaSigmaSupport" -> grades|>, unit, epsilonOrder]];

Do[
  auditRho = rhoBrBackground[auditDensity, auditWidth];
  (* Reference-normalized physical coordinate is supplied before varying T. *)
  auditDeltaW = WZero eWField;
  auditT = auditRho D[uVector, time].D[uVector, time]/2 + muW D[auditDeltaW, time]^2/2;
  auditInertia = (D[D[auditT, D[#, time]], time] - D[auditT, #]) & /@ auditFields;
  auditNegativeInertia = (D[D[-auditT, D[#, time]], time] - D[-auditT, #]) & /@ auditFields;
  auditMutatedT = auditRho D[uVector, time].D[uVector, time]/2 + muW D[auditWidth eWField, time]^2/2;
  auditMutatedInertia = (D[D[auditMutatedT, D[#, time]], time] - D[auditMutatedT, #]) & /@ auditFields;
  auditPhysical["wlRepairKineticAction" <> ToString[auditDensityIndex], auditT, {-1, -2, 1}, 2];
  Do[
    auditNative = Block[{rhoBrValue = auditRho, rhoBrValueLive = auditRho},
      Join[ReleaseHold[auditConstructor["kineticU"]], {ReleaseHold[auditConstructor["kineticEw"]]}]];
    Do[auditPhysical["wlRepairInertia" <> ToString[auditDensityIndex] <> auditConstructor["name"] <> ToString[i],
      {auditNative[[i]], auditInertia[[i]], auditNative[[i]] - auditInertia[[i]],
        auditNegativeInertia[[i]] - auditInertia[[i]], auditMutatedInertia[[i]] - auditInertia[[i]]},
      auditRowUnits[[i]], 1], {i, 4}], {auditConstructor, auditConstructors}],
  {auditDensityIndex, Length[densities]}, {auditDensity, {densities[[auditDensityIndex]]}}];

(* Pure assembly coefficients are read from each constructor without replacing
   a computed energy or force by a fabricated zero. Formal symbols here denote
   independent coefficient slots, and are labelled as such. *)
Do[
  auditAssembly = ReleaseHold[auditConstructor["operator"]] /.
    {kineticULive -> kineticU, kineticEwLive -> kineticEw,
     rowsLive -> rows, faceRowsLive -> faceRows};
  auditAssemblyOperands = {auditAssembly["U_MOMENTUM_ROWS"], auditAssembly["THICKNESS_ROW"]};
  auditAssemblySlots = {{kineticU, rows["U_INTERNAL"], faceRows["U_FACE"]},
    {kineticEw, rows["EW_INTERNAL"], faceRows["EW_FACE"]}};
  auditEmit["wlRepairAssemblyCoefficients" <> auditConstructor["name"],
    MapThread[Function[{value, slots}, Coefficient[value, #] & /@ slots],
      {auditAssemblyOperands, auditAssemblySlots}], {0, 0, 0}], {auditConstructor, auditConstructors}];

auditProbeU = auditB thetaField^2/2 + auditC thetaField eWField + auditK eWField^2/2 +
  auditG Total[Flatten[Outer[D, uVector, spatialCoordinates]]^2]/2;
auditCaseNumber = 0;
Do[
  auditCaseNumber++;
  auditConstraint = virtualConstraintSource[auditRoute, auditBranch, auditDensity, widthBase];
  auditThetaRule = First[Solve[auditConstraint == 0, virtualThetaField]];
  (* First variation computed from the probe before spatial integration by parts. *)
  auditVariedProbe = Total[MapThread[Function[{f, v}, D[auditProbeU, f] v +
    Sum[D[auditProbeU, D[f, x]] D[v, x], {x, spatialCoordinates}]],
    {Join[uVector, {thetaField, eWField}], Join[virtualUVector, {virtualThetaField, virtualEwField}]}]];
  auditVariedProbe = auditVariedProbe /. auditThetaRule;
  auditIndependentRows = Table[D[auditVariedProbe, v] - Sum[D[D[auditVariedProbe, D[v, x]], x], {x, spatialCoordinates}], {v, Take[auditVirtuals, 4]}];
  Do[
    auditRows = auditVariationFunction[auditProbeU, auditConstraint];
    auditNativeRows = Activate[Join[auditRows["U_INTERNAL"], {auditRows["EW_INTERNAL"]}]];
    Do[auditPhysical["wlRepairStoredProbe" <> ToString[auditCaseNumber] <> SymbolName[auditVariationFunction] <> ToString[i],
      {auditNativeRows[[i]], auditIndependentRows[[i]], auditNativeRows[[i]] - auditIndependentRows[[i]]}, auditRowUnits[[i]], 1], {i, 4}],
    {auditVariationFunction, {constrainedRows, constrainedRowsWithLiveEnergyEL}}];
  auditEmit["wlRepairMechanicalCase" <> ToString[auditCaseNumber], {auditRoute, auditBranch, auditDensity}, {0, 0, 0}];
  Do[
    auditFace = faceSources[auditRoute, auditBranch, auditFaceSign, widthBase, auditChemical, rhoBrBackground[auditDensity, auditWidth]];
    auditNativeFaceRows = faceGeneralizedRows[<|"singleFace" -> auditFace|>];
    auditNativeForceRow = Activate[Join[auditNativeFaceRows["U_FACE"], {auditNativeFaceRows["EW_FACE"], auditNativeFaceRows["CENTER_FACE"]}]];
    (* Independent traction work from the material face embedding and its area
       covector. The prescribed pressure/chemical response is held fixed. *)
    auditMaterialX = spatialCoordinates + auditVariation virtualUVector;
    auditHeight = auditFaceSign widthBase[Sequence @@ If[auditBranch === "LAB_HELD", auditMaterialX, spatialCoordinates]]/2 +
      auditVariation (virtualZetaCenterField + auditFaceSign WZero virtualEwField/2);
    auditVirtualPosition = D[Join[auditMaterialX, {auditHeight}], auditVariation] /. auditVariation -> 0;
    auditGraph = auditFaceSign auditWidth/2;
    auditCofactor = auditFaceSign Join[-(D[auditGraph, #] & /@ spatialCoordinates), {1}];
    auditTractionPressure = pressureField[auditFaceSign] + lambdaXResponse (auditChemical/rhoBrBackground[auditDensity, auditWidth] - pressureField[auditFaceSign]/rhoM);
    auditIndependentWork = -auditTractionPressure auditCofactor.auditVirtualPosition;
    auditForce = D[auditIndependentWork, #] & /@ auditVirtuals;
    auditVelocityRules = Thread[auditVirtuals -> Join[D[uVector, time], {D[eWField, time], D[zetaCenterField, time]}]];
    auditPhysical["wlRepairFaceWork" <> ToString[auditCaseNumber] <> If[auditFaceSign === 1, "Plus", "Minus"],
      {auditFace["VIRTUAL_WORK"], auditIndependentWork, auditFace["VIRTUAL_WORK"] - auditIndependentWork}, {-1, -2, 1}, 1];
    Do[auditPhysical["wlRepairFaceRow" <> ToString[auditCaseNumber] <> If[auditFaceSign === 1, "Plus", "Minus"] <> ToString[i],
      {auditNativeForceRow[[i]], auditForce[[i]], auditNativeForceRow[[i]] + auditForce[[i]]}, auditRowUnits[[i]], 1], {i, 5}];
    auditPower = auditForce.(auditVirtuals /. auditVelocityRules);
    auditNativePower = -auditFace["MEASURE"] auditFace["TRACTION_PRESSURE"] auditFace["NORMAL_VELOCITY"];
    auditPhysical["wlRepairFacePower" <> ToString[auditCaseNumber] <> If[auditFaceSign === 1, "Plus", "Minus"],
      {auditPower, auditNativePower, auditPower - auditNativePower}, {-1, -3, 1}, 2];
    Do[
      auditMutatedFace = Join[auditFace, <|"VIRTUAL_WORK" -> auditMutation auditFace["VIRTUAL_WORK"]|>];
      auditMutatedRows = faceGeneralizedRows[<|"singleFace" -> auditMutatedFace|>];
      auditMutatedRows = Activate[Join[auditMutatedRows["U_FACE"], {auditMutatedRows["EW_FACE"], auditMutatedRows["CENTER_FACE"]}]];
      Do[auditPhysical["wlRepairLoadMutation" <> ToString[auditCaseNumber] <> If[auditFaceSign === 1, "Plus", "Minus"] <> If[auditMutation === -1, "Sign", "Scale"] <> ToString[i],
        auditMutatedRows[[i]] + auditForce[[i]], auditRowUnits[[i]], 1], {i, 5}], {auditMutation, {-1, 2}}], {auditFaceSign, faces}],
  {auditRoute, {"EULERIAN", "MATERIAL"}}, {auditBranch, branches}, {auditDensity, densities}];
auditEmit["wlRepairMechanicalCompletion", <|"cases" -> auditCaseNumber, "emittedRecords" -> Length[auditKeys]|>, {0, 0, 0}];
auditResidualFamilies = {"wlRepairInertia", "wlRepairStoredProbe", "wlRepairFaceWork", "wlRepairFaceRow", "wlRepairFacePower"};
auditSummary = Association[Table[auditPrefix -> Module[{keys, residuals},
  keys = Select[auditKeys, StringStartsQ[#, auditPrefix] &];
  residuals = auditRecords[#]["value"]["retained"][[3]] & /@ keys;
  <|"records" -> Length[keys], "zeroResiduals" -> Count[residuals, 0],
    "nonzeroResiduals" -> Association[MapThread[Rule, {Pick[keys, (# =!= 0 & /@ residuals)],
      ToString[#, InputForm] & /@ Select[residuals, # =!= 0 &]}]]|>], {auditPrefix, auditResidualFamilies}]];
AssociateTo[auditSummary, "mutationRetainedNonzeroComponents" -> Count[
  (auditRecords[#]["value"]["retained"] & /@ Select[auditKeys, StringStartsQ[#, "wlRepairLoadMutation"] &]), Except[0]]];
AssociateTo[auditSummary, "inertiaSignActionNonzeroComponents" -> Count[
  (auditRecords[#]["value"]["retained"][[4]] & /@ Select[auditKeys, StringStartsQ[#, "wlRepairInertia"] &]), Except[0]]];
AssociateTo[auditSummary, "coordinateActionNonzeroComponents" -> Count[
  (auditRecords[#]["value"]["retained"][[5]] & /@ Select[auditKeys, StringStartsQ[#, "wlRepairInertia"] &]), Except[0]]];
AssociateTo[auditSummary, "assemblyCoefficients" -> Table[auditRecords["wlRepairAssemblyCoefficients" <> name]["value"], {name, {"Raw", "Frozen", "Live"}}]];
Export[FileNameJoin[{auditRoot, "summary.json"}], auditSummary, "RawJSON"];
