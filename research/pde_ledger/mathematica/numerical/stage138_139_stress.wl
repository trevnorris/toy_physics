(* Retired Stage 138/139 plain harness wrapper.

   The exploratory plain harness was weaker than the fixed-point harness and is
   kept only as a compatibility entry point. Delegate to the authoritative
   fixed-point stress script instead of maintaining a second implementation. *)

ClearAll["Global`*"];

rootDir = DirectoryName[$InputFileName];
targetScript = FileNameJoin[{rootDir, "stage138_139_fixedpoint_stress.wl"}];
scriptArgs = If[Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine], {}];

If[Length[scriptArgs] >= 1,
  stage138139FixedpointConfigOverride = scriptArgs[[1]];
];

Print[
  "stage138_139_stress.wl is retired; delegating to "
  <> "stage138_139_fixedpoint_stress.wl"
];
Get[targetScript];
