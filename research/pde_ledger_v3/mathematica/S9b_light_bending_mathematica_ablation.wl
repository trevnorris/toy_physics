(* One guarded kernel. The canonical engine is loaded unchanged for baseline.
   Each copy has exactly one manifest replacement. No subprocess kernels.
   Run from this repository root; all harness artifacts stay in _scratch. *)
Begin["S9bHarness`"];
$HistoryLength = 0;
(* Optional scope changes execution selection only, never a knife or payload. *)
scope = Replace[Environment["S9B_ABLATION_SCOPE"], $Failed -> "ALL"];
If[!MemberQ[{"ALL", "K11"}, scope], Quit[94]];
root = Directory[];
engine = FileNameJoin[{root, "research", "pde_ledger_v3", "mathematica", "S9b_light_bending_mathematica_audit.wl"}];
work = FileNameJoin[{root, "_scratch", "s9b_wl_repair4", "ablation-final"}];
If[!DirectoryQ[work], CreateDirectory[work, CreateIntermediateDirectories -> True]];
source = Import[engine, "Text"];
S9bHarness`progress = OpenWrite[FileNameJoin[{work, "progress.log"}]];
note[stage_, name_] := (WriteString[S9bHarness`progress,
  ToString[{stage, name, MemoryInUse[], MaxMemoryUsed[]}, InputForm] <> "\n"]; Flush[S9bHarness`progress]);
knives = <|
 "K1" -> {"kineticSymbol = (omega - vVector.kCov)^2;", "kineticSymbol = omega^2;"},
 "K2" -> {"metric2 = DiagonalMatrix[{aMetric, r^2}];", "metric2 = DiagonalMatrix[{1, r^2}];"},
 "K3" -> {"dispersionSpeed = cSquared;", "dispersionSpeed = c0^2;"},
 "K4a" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] /. HoldPattern[Derivative[j_][deltaProfile][arg_]] :> 0;"},
 "K4b" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] /. HoldPattern[Derivative[j_][velocityProfile][arg_]] :> 0;"},
 "K4c" -> {"xiSlope = Simplify[D[xiAnsatz, r], $Assumptions];", "xiSlope = Simplify[D[xiConstant, r], $Assumptions];"},
 "K5" -> {"returnDirection = -1;", "returnDirection = 1;"},
 "K6" -> {"radarSource = roundTrip;", "radarSource = oneWayER;"},
 "K7" -> {"quantifierMode = \"EVERY_B\";", "quantifierMode = \"FIXED_B\";"},
 "K8" -> {"massDensity = rhoBr[r];", "massDensity = rhoBrConstant;"},
 "K9a" -> {"localSoundSquared = bulkSoundSquared /. rho -> rho0 (1 + fLocal);", "localSoundSquared = bulkSoundSquared /. rho -> rho0;"},
 "K9b" -> {"localBulkDensity = rho0 (1 + fLocal);", "localBulkDensity = rho0;"},
 "K10" -> {"branchGateLocal = propagatingLocal;", "branchGateLocal = Not[propagatingLocal];"},
 "K11" -> {"vVector = {vLocal, 0};", "vVector = {vLocal, vLocal/r};"}|>;
print[name_, value_] := (WriteString[First[$Output], "WL_LOCAL_S9B_ABLATION_" <> name <> ": " <>
  ToString[value, InputForm, PageWidth -> Infinity] <> "\n"]; Flush[First[$Output]]);
(* BEGIN_DOMAIN_DIFFERENCE_METHOD
   Domain attachment is a property of the engine's payload schema, not of
   a predicate's head. In particular Closedness, ProfileEquation, Constraint,
   ConstraintOnOtherInputs and NoRealVelocityDomain remain computed VALUES.
   ConditionalExpression conditions and the named attachment fields below
   are collected at their complete structural paths, including nested cases.
   A predicate emitted alone is a value, not an attached domain. Coverage is
   never passed to this projection (K11's sample domain is not an attachment).
   Baseline and corrupted payloads are always printed in full separately. *)
attachedDomainKeys = {"Domain", "Domains", "BranchDomain", "AlongRayDomain",
  "ConstantDomain", "GeneralConstantDomain", "FunctionalDomain", "Premise",
  "Component", "ZeroSetJoining"};
splitPayload[payload_] := Module[{domains = <||>, walk, value},
  walk[z_Association, path_] := Association[KeyValueMap[
    Function[{key, item}, If[MemberQ[attachedDomainKeys,key],
      AssociateTo[domains, ToString[Append[path,key],InputForm] -> item]; Nothing,
      key -> walk[item,Append[path,key]]]], z]];
  walk[z_List, path_] := MapIndexed[walk[#1,Join[path,#2]] &,z];
  walk[z_ConditionalExpression, path_] := (
    AssociateTo[domains, ToString[Append[path,"ConditionalExpression"],InputForm] -> z[[2]]];
    walk[z[[1]],path]);
  (* Traverse other heads as well: nested ConditionalExpression attachments
     are separated without changing a relation into a domain. *)
  walk[z_, path_] := If[AtomQ[z] || FreeQ[z,_ConditionalExpression],z,
    MapIndexed[walk[#1,Join[path,#2]] &,z]];
  value = walk[payload,{}]; <|"Value" -> value,"Domain" -> domains|>];
logicalObjectQ[z_] := MemberQ[{True,False},z] ||
  MemberQ[{Equal,Unequal,Less,LessEqual,Greater,GreaterEqual,Inequality,
    And,Or,Not,Xor,Equivalent,Implies,Element,Exists,ForAll},Head[z]] ||
  MemberQ[{Inactive[ForAll],Inactive[Exists]},Head[z]];
(* Canonical exact subtraction, with an equality short circuit before any
   traversal. Logical leaves use signed indicator subtraction, including
   relations/predicates used as values. No assumption, sampling, global
   Simplify/Reduce or exponential Boolean expansion is introduced. *)
partDifference[a_,b_] := Which[
  SameQ[a,b], 0,
  AssociationQ[a] && AssociationQ[b], AssociationMap[
    partDifference[Lookup[a,#,Missing["Unemitted"]],Lookup[b,#,Missing["Unemitted"]]] &,
    Union[Keys[a],Keys[b]]],
  ListQ[a] && ListQ[b] && Length[a] == Length[b], MapThread[partDifference,{a,b}],
  logicalObjectQ[a] && logicalObjectQ[b], Boole[a]-Boole[b],
  Head[a] === Rule && Head[b] === Rule && SameQ[First[a],First[b]],
    First[a] -> partDifference[Last[a],Last[b]],
  ListQ[a] || ListQ[b] || StringQ[a] || StringQ[b] ||
    !FreeQ[{a,b}, _String | _Association | _Rule | _Missing | _Piecewise],
    <|"Changed" -> UnsameQ[a,b],"Corrupted" -> a,"Baseline" -> b|>,
  True,a-b];
(* Attachment predicates may be a bare Boolean parameter. Metadata-valued
   domains retain structural comparison; symbolic predicates use indicators. *)
attachedDifference[a_,b_] := Which[SameQ[a,b],0,
  (logicalObjectQ[a] || Head[a] === Symbol) &&
    (logicalObjectQ[b] || Head[b] === Symbol), Boole[a]-Boole[b],
  True,partDifference[a,b]];
domainDifference[a_Association,b_Association] := If[SameQ[a,b],0,
  AssociationMap[attachedDifference[Lookup[a,#,Missing["Unemitted"]],
    Lookup[b,#,Missing["Unemitted"]]] &,Union[Keys[a],Keys[b]]]];
difference[a_,b_] := Module[{partsA = splitPayload[a],partsB = splitPayload[b]},
  <|"ValueDifference" -> partDifference[partsA["Value"],partsB["Value"]],
    "DomainDifference" -> domainDifference[partsA["Domain"],partsB["Domain"]]|>];
(* END_DOMAIN_DIFFERENCE_METHOD *)
run[path_, name_] := Module[{capture, result},
  ClearSystemCache[]; note["engine begin", name];
  capture = OpenWrite[FileNameJoin[{work, name <> ".stdout"}]];
  Block[{$Context = "Global`", $ContextPath = {"Global`", "System`"}, $Output = {capture}}, Get[path]];
  Close[capture]; result = Global`emitted; note["engine end", name]; result];
mutate[name_] := Module[{pair = knives[name], target},
  print[name <> "_SITE", <|"Construction" -> First[pair], "Mutation" -> Last[pair],
    "Occurrences" -> StringCount[source, First[pair]]|>];
  If[StringCount[source, First[pair]] != 1, Quit[92]];
  target = FileNameJoin[{work, name <> ".wl"}];
  Export[target, StringReplace[source, Rule @@ pair], "Text"]; target];
baseline = run[engine, "BASELINE"];
print["MANIFEST", knives];
Do[copy = mutate[knife]; corrupted = run[copy, knife];
  Do[note["difference begin", knife <> "_" <> tag]; print[knife <> "_" <> tag, <|"Baseline" -> Lookup[baseline, tag, Missing["Unemitted"]],
    "Corrupted" -> Lookup[corrupted, tag, Missing["Unemitted"]],
    "Difference" -> difference[Lookup[corrupted, tag, Missing["Unemitted"]],
      Lookup[baseline, tag, Missing["Unemitted"]]]|>];
    note["difference end", knife <> "_" <> tag], {tag, Union[Keys[baseline], Keys[corrupted]]}],
  {knife, If[scope == "K11", {}, DeleteCases[Keys[knives], "K11"]]}];
(* Compact evaluation is exact symbolic prefix evaluation. No points, seeds,
   finite precision or extra grade truncation. Both copies have the same
   extraction boundary. Later tags are explicitly unreached. *)
copy = mutate["K11"];
(* Fully qualified harness state survives ClearAll["Global`*"]. Prepare
   and validate BOTH extracts before executing either. No fallback to a full
   source is possible when the marker is missing or extraction is malformed. *)
S9bHarness`compactBoundary = "(* COMPACT_K11_BOUNDARY *)";
S9bHarness`prepareCompact[path_, name_] := Module[{text, target, position, prefix},
  text = Import[path, "Text"];
  If[!StringQ[S9bHarness`compactBoundary] ||
    StringCount[text, S9bHarness`compactBoundary] != 1, Quit[93]];
  position = First[First[StringPosition[text, S9bHarness`compactBoundary]]];
  prefix = StringTake[text, position - 1];
  If[!StringQ[prefix] || StringLength[prefix] == 0 ||
    StringLength[prefix] >= StringLength[text] ||
    StringContainsQ[prefix, S9bHarness`compactBoundary], Quit[93]];
  target = FileNameJoin[{work, name <> ".wl"}];
  Export[target, prefix, "Text"];
  <|"Path" -> target, "Text" -> prefix|>];
S9bHarness`compactInputs = {
  S9bHarness`prepareCompact[engine, "K11_COMPACT_BASELINE"],
  S9bHarness`prepareCompact[copy, "K11_COMPACT_CORRUPTED"]};
If[S9bHarness`compactInputs[[2]]["Text"] =!=
  StringReplace[S9bHarness`compactInputs[[1]]["Text"], Rule @@ knives["K11"]], Quit[93]];
compactBaseline = run[S9bHarness`compactInputs[[1]]["Path"], "K11_COMPACT_BASELINE"];
compactCorrupted = run[S9bHarness`compactInputs[[2]]["Path"], "K11_COMPACT_CORRUPTED"];
coverage = <|"Method" -> "exact symbolic mechanical prefix extraction",
  "Boundary" -> S9bHarness`compactBoundary, "RetainedGrades" -> Global`grades,
  "LocalConstructions" -> "untruncated dispersion/Fermat and retained-ring nonreciprocal one-form",
  "AdditionalTruncation" -> None, "SampleDomain" -> "symbolic constructor domains",
  "Arithmetic" -> "exact", "Samples" -> {}, "Seeds" -> {}|>;
Do[note["difference begin", "K11_" <> tag]; print["K11_" <> tag, If[KeyExistsQ[compactBaseline, tag] && KeyExistsQ[compactCorrupted, tag],
  <|"Baseline" -> baseline[tag], "CompactBaseline" -> compactBaseline[tag],
    "CompactCorrupted" -> compactCorrupted[tag],
    "Difference" -> difference[compactCorrupted[tag], compactBaseline[tag]],
    "MethodResidual" -> difference[compactBaseline[tag], baseline[tag]], "Coverage" -> coverage|>,
  <|"Baseline" -> baseline[tag], "CompactBaseline" -> Missing["NOT_EVALUATED"],
    "CompactCorrupted" -> Missing["NOT_EVALUATED"], "Difference" -> Missing["NOT_EVALUATED"],
    "MethodResidual" -> Missing["NOT_EVALUATED"], "Coverage" -> coverage|>]]; note["difference end", "K11_" <> tag], {tag, Keys[baseline]}];
(* Harness dead-path ablation: both inputs are the actual baseline. *)
If[scope == "ALL", Do[note["difference begin", "HARNESS_DEAD_PATH_" <> tag]; print["HARNESS_DEAD_PATH_" <> tag, <|"Baseline" -> baseline[tag],
 "Corrupted" -> baseline[tag], "Difference" -> difference[baseline[tag], baseline[tag]]|>]; note["difference end", "HARNESS_DEAD_PATH_" <> tag], {tag, Keys[baseline]}]];
note["complete", scope]; Close[S9bHarness`progress];
End[];
