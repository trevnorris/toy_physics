(* Fixed four-engine harness, WL leg. The Python driver supplies the exact
   production, budget and knife sources in a fresh /tmp tree for each label.
   Held complete-source validation precedes Get of the complete production
   entrypoint. No extractor, emitter, serializer, or engine object is patched.
   BUDGET_RESTRICTION / NON-KNIFE: only the buildCase driver iterator and the
   adaptive probeCase draw assignment. All labels use the same restrictions. *)
Begin["S11CC2Harness`"];
plan = Import[Environment["S11CC2_HARNESS_PLAN"], "RawJSON"];
parseHeld[s_String] := Quiet[Check[ToExpression[s, InputForm, HoldComplete], $Failed]];
production = Import[plan["production"], "Text"];
budget = Import[plan["budget"], "Text"];
executed = Import[plan["executed"], "Text"];
pHeld = parseHeld[production]; bHeld = parseHeld[budget]; eHeld = parseHeld[executed];
If[MemberQ[{pHeld, bHeld, eHeld}, $Failed], Exit[81]];
(* Global context is explicit: this file has a private driver context. *)
End[];
(* Parse the source in Global`, exactly as Get will. HoldComplete prevents
   ClearAll, Do, Quit and every construction from evaluating during matching. *)
S11CC2Harness`pHeld = Quiet[Check[ToExpression[S11CC2Harness`production, InputForm, HoldComplete], $Failed]];
S11CC2Harness`bHeld = Quiet[Check[ToExpression[S11CC2Harness`budget, InputForm, HoldComplete], $Failed]];
S11CC2Harness`eHeld = Quiet[Check[ToExpression[S11CC2Harness`executed, InputForm, HoldComplete], $Failed]];
S11CC2Harness`drivers = Cases[S11CC2Harness`pHeld,
  x : HoldPattern[Do[body_, iterator_]] /; !FreeQ[Unevaluated[body], HoldPattern[buildCase[case]]] :> HoldComplete[x], {2}];
S11CC2Harness`probes = Cases[S11CC2Harness`pHeld,
  x : HoldPattern[SetDelayed[probeCase[___], _]] :> HoldComplete[x], {2}];
S11CC2Harness`builds = Cases[S11CC2Harness`pHeld,
  x : HoldPattern[SetDelayed[buildCase[___], _]] :> HoldComplete[x], {2}];
If[Length[S11CC2Harness`drivers] != 1 || Length[S11CC2Harness`probes] != 1 ||
   Length[S11CC2Harness`builds] != 1, Exit[82]];
S11CC2Harness`drawAssignments = Cases[First[S11CC2Harness`probes],
  HoldPattern[Set[drawCount, If[0 < bound < 1,
    Max[8, Ceiling[Log[2^-80/Max[1, unionCount]]/Log[bound]]], 8]]] :> 1, Infinity];
S11CC2Harness`iteratorCount = Count[First[S11CC2Harness`drivers],
  HoldPattern[{case, Tuples[{{"LAB_HELD", "MATERIAL_ADVECTED"}, {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}]}], Infinity];
If[Length[S11CC2Harness`drawAssignments] != 1 || S11CC2Harness`iteratorCount != 1, Exit[83]];
(* The driver applies literal replacements only in balanced top-level scopes.
   Match knife ownership in the held tree, without evaluating replacement Set
   expressions or rewriting any production definition in the running kernel. *)
S11CC2Harness`knifeScopeCount = Switch[S11CC2Harness`plan["scope"],
  "materialNormalKnife", Count[S11CC2Harness`pHeld, HoldPattern[Set[materialNormalKnife, 0]], {2}],
  "actualJunkCoefficient", Count[S11CC2Harness`pHeld, HoldPattern[Set[actualJunkCoefficient, 0]], {2}],
  "buildCase", Count[First[S11CC2Harness`builds],
    HoldPattern[Set[sourceChannel, joinFaceMaps[buildContraction[cm[#], es[#] - ms[#], response[#], #, False] &]]], Infinity],
  "NONE", 0, _, -1];
If[S11CC2Harness`knifeScopeCount != S11CC2Harness`plan["knife_count"], Exit[85]];
Export[S11CC2Harness`plan["guard"], <|"wolfram_version" -> $Version,
  "held_parseable" -> 1, "driver_count" -> Length[S11CC2Harness`drivers],
  "probe_assignment_count" -> Length[S11CC2Harness`drawAssignments],
  "knife_count" -> S11CC2Harness`plan["knife_count"]|>, "RawJSON"];
(* The engine clears Global` and exits itself. The path is substituted first. *)
With[{path = S11CC2Harness`plan["executed"]}, Get[path]];
