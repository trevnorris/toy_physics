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
work = FileNameJoin[{root, "_scratch", "s9b_wl_repair1", "k11_repair", "ablation"}];
If[!DirectoryQ[work], CreateDirectory[work, CreateIntermediateDirectories -> True]];
source = Import[engine, "Text"];
progress = OpenWrite[FileNameJoin[{work, "progress.log"}]];
note[stage_, name_] := (WriteString[progress,
  ToString[{stage, name, MemoryInUse[], MaxMemoryUsed[]}, InputForm] <> "\n"]; Flush[progress]);
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
(* Exact structural subtraction, with no global Simplify/Reduce call.
   Associations/lists are traversed without omitting leaves. Identical
   leaves compute to zero; arithmetic leaves use exact canonical subtraction.
   Nonnumeric leaves (including string-valued Piecewise and logical gates)
   retain the complete ordered operands in Inactive[Subtract], as do domains.
   This is a full difference representation, not an approximation or a
   compact evaluation. No profile, grade, term or condition is discarded. *)
difference[a_Association, b_Association] := AssociationMap[
  difference[Lookup[a, #, Missing["Unemitted"]], Lookup[b, #, Missing["Unemitted"]]] &,
  Union[Keys[a], Keys[b]]];
difference[a_List, b_List] /; Length[a] == Length[b] := MapThread[difference, {a,b}];
difference[a_, b_] := Which[SameQ[a,b], 0,
  StringQ[a] || StringQ[b] || MemberQ[{True, False}, a] || MemberQ[{True, False}, b] ||
    !FreeQ[{a,b}, _String | _Association | _Rule | _Missing | _ConditionalExpression |
      _Piecewise | _Equal | _Unequal | _Less | _LessEqual | _Greater | _GreaterEqual |
      _Inequality | _And | _Or | _Not | _Xor | _Element | _Exists | _ForAll],
    Inactive[Subtract][a,b],
  True, a-b];
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
Do[print["K11_" <> tag, If[KeyExistsQ[compactBaseline, tag] && KeyExistsQ[compactCorrupted, tag],
  <|"Baseline" -> baseline[tag], "CompactBaseline" -> compactBaseline[tag],
    "CompactCorrupted" -> compactCorrupted[tag],
    "Difference" -> difference[compactCorrupted[tag], compactBaseline[tag]],
    "MethodResidual" -> difference[compactBaseline[tag], baseline[tag]], "Coverage" -> coverage|>,
  <|"Baseline" -> baseline[tag], "CompactBaseline" -> Missing["NOT_EVALUATED"],
    "CompactCorrupted" -> Missing["NOT_EVALUATED"], "Difference" -> Missing["NOT_EVALUATED"],
    "MethodResidual" -> Missing["NOT_EVALUATED"], "Coverage" -> coverage|>]], {tag, Keys[baseline]}];
(* Harness dead-path ablation: both inputs are the actual baseline. *)
If[scope == "ALL", Do[print["HARNESS_DEAD_PATH_" <> tag, <|"Baseline" -> baseline[tag],
 "Corrupted" -> baseline[tag], "Difference" -> difference[baseline[tag], baseline[tag]]|>], {tag, Keys[baseline]}]];
note["complete", scope]; Close[progress];
End[];
