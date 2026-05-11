(* Retired Stage 155/156 plain harness wrapper.

   The exploratory plain harness was weaker than the fixed-point harness and is
   kept only as a compatibility entry point. Delegate to the authoritative
   fixed-point stress script instead of maintaining a second implementation. *)

ClearAll["Global`*"];

rootDir = DirectoryName[$InputFileName];
targetScript = FileNameJoin[{rootDir, "stage155_156_fixedpoint_stress.wl"}];
scriptArgs = If[Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine], {}];

If[Length[scriptArgs] >= 1,
  stage155139FixedpointConfigOverride = scriptArgs[[1]];
];

Print[
  "stage155_156_stress.wl is retired; delegating to "
  <> "stage155_156_fixedpoint_stress.wl"
];
Get[targetScript];
