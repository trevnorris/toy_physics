(* One guarded kernel. The canonical engine is loaded unchanged for baseline.
   Each copy has exactly one manifest replacement. No subprocess kernels.
   Run from this repository root; all harness artifacts stay in _scratch. *)
Begin["S9bHarness`"];
$HistoryLength = 0;
root = Directory[];
engine = FileNameJoin[{root, "research", "pde_ledger_v3", "mathematica", "S9b_light_bending_mathematica_audit.wl"}];
work = FileNameJoin[{root, "_scratch", "s9b_wl_ablation"}];
If[!DirectoryQ[work], CreateDirectory[work, CreateIntermediateDirectories -> True]];
source = Import[engine, "Text"];
knives = <|
 "K1" -> {"kineticSymbol = (omega - vVector.kCov)^2;", "kineticSymbol = omega^2;"},
 "K2" -> {"metric2 = DiagonalMatrix[{aMetric, r^2}];", "metric2 = DiagonalMatrix[{1, r^2}];"},
 "K3" -> {"dispersionSpeed = cSquared;", "dispersionSpeed = c0^2;"},
 "K4a" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] + powD ampD D[z, ampD]/x;"},
 "K4b" -> {"rayDerivative[z_, x_] := D[z, x];", "rayDerivative[z_, x_] := D[z, x] + powV ampV D[z, ampV]/x;"},
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
(* Structural differences retain nonnumeric leaves; zero is emitted for
   identical leaves, including metadata. There is no expected outcome. *)
difference[a_Association, b_Association] := AssociationMap[
  difference[Lookup[a, #, Missing["Unemitted"]], Lookup[b, #, Missing["Unemitted"]]] &,
  Union[Keys[a], Keys[b]]];
difference[a_List, b_List] /; Length[a] == Length[b] := MapThread[difference, {a,b}];
difference[a_, b_] := Which[SameQ[a,b], 0,
  StringQ[a] || StringQ[b] || MemberQ[{True, False}, a] || MemberQ[{True, False}, b] ||
    !FreeQ[{a,b}, _Association | _Rule | _Missing | _ConditionalExpression], Inactive[Subtract][a,b],
  True, Simplify[a-b]];
run[path_, name_] := Module[{capture, result},
  capture = OpenWrite[FileNameJoin[{work, name <> ".stdout"}]];
  Block[{$Context = "Global`", $ContextPath = {"Global`", "System`"}, $Output = {capture}}, Get[path]];
  Close[capture]; result = Global`emitted; result];
mutate[name_] := Module[{pair = knives[name], target},
  print[name <> "_SITE", <|"Construction" -> First[pair], "Mutation" -> Last[pair],
    "Occurrences" -> StringCount[source, First[pair]]|>];
  If[StringCount[source, First[pair]] != 1, Quit[92]];
  target = FileNameJoin[{work, name <> ".wl"}];
  Export[target, StringReplace[source, Rule @@ pair], "Text"]; target];
baseline = run[engine, "BASELINE"];
print["MANIFEST", knives];
Do[copy = mutate[knife]; corrupted = run[copy, knife];
  Do[print[knife <> "_" <> tag, <|"Baseline" -> Lookup[baseline, tag, Missing["Unemitted"]],
    "Corrupted" -> Lookup[corrupted, tag, Missing["Unemitted"]],
    "Difference" -> difference[Lookup[corrupted, tag, Missing["Unemitted"]],
      Lookup[baseline, tag, Missing["Unemitted"]]]|>], {tag, Union[Keys[baseline], Keys[corrupted]]}],
  {knife, DeleteCases[Keys[knives], "K11"]}];
(* Compact evaluation is exact symbolic prefix evaluation. No points, seeds,
   finite precision or extra grade truncation. Both copies have the same
   extraction boundary. Later tags are explicitly unreached. *)
copy = mutate["K11"];
boundary = "(* COMPACT_K11_BOUNDARY *)";
compact[path_, name_] := Module[{text, target},
  text = Import[path, "Text"];
  target = FileNameJoin[{work, name <> ".wl"}];
  Export[target, First[StringSplit[text, boundary]], "Text"];
  run[target, name]];
compactBaseline = compact[engine, "K11_COMPACT_BASELINE"];
compactCorrupted = compact[copy, "K11_COMPACT_CORRUPTED"];
coverage = <|"Method" -> "exact symbolic mechanical prefix extraction",
  "Boundary" -> boundary, "RetainedGrades" -> Global`grades,
  "LocalConstructions" -> "untruncated dispersion and Fermat one-form",
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
Do[print["HARNESS_DEAD_PATH_" <> tag, <|"Baseline" -> baseline[tag],
 "Corrupted" -> baseline[tag], "Difference" -> difference[baseline[tag], baseline[tag]]|>], {tag, Keys[baseline]}];
End[];
