(* Path-A B2b Spike-1 VSH Maxwell operator input ledger.

   Dispatched by mt15_03_spike1_vsh_maxwell_operator.wls --patha-b2b.
   This worker verifies that the Path-A route reads the self-consistent A0
   background and frozen geometry before Spike-2 performs the Green transfer.
*)

ClearAll["Global`pathaB2BSpike1`*"];

args = $ScriptCommandLine;

argValue[flag_String, default_] := Module[{pos},
  pos = FirstPosition[args, flag, Missing["NotFound"]];
  If[MissingQ[pos] || pos[[1]] >= Length[args], default, args[[pos[[1]] + 1]]]
];

toJsonNumber[x_] := Module[{n = Quiet[Check[N[x], $Failed]]},
  If[n === $Failed || ! NumericQ[n],
    Print["ERROR: nonnumeric JSON value: ", InputForm[x]];
    Exit[4]
  ];
  N[Re[n], 17]
];
toJsonVector[x_List] := toJsonNumber /@ x;

clipToGrid[x_, xs_List] := Min[Max[x, First[xs]], Last[xs]];
interp1Clamped[xs_List, ys_List, x_] := Module[{xc = clipToGrid[x, xs], k, t},
  If[xc <= First[xs], Return[First[ys]]];
  If[xc >= Last[xs], Return[Last[ys]]];
  k = First[FirstPosition[UnitStep[xs - xc], 1]] - 1;
  If[k < 1, k = 1];
  If[k >= Length[xs], k = Length[xs] - 1];
  t = (xc - xs[[k]])/(xs[[k + 1]] - xs[[k]]);
  (1 - t) ys[[k]] + t ys[[k + 1]]
];

backgroundPath = argValue["--background", ""];
outPath = argValue["--out", ""];
nr = ToExpression[argValue["--nr", "7"]];
nw = ToExpression[argValue["--nw", "7"]];
radialScale = ToExpression[argValue["--radial-scale", "1.0"]];

If[backgroundPath == "" || outPath == "",
  Print["ERROR: --background and --out are required for --patha-b2b Spike-1"];
  Exit[2]
];

background = Import[backgroundPath, "RawJSON"];
geometry = background["geometry"];
constants = background["constants"];

a = N[geometry["a"], 17];
L = N[geometry["L"], 17];
rOuter = N[background["grid"]["r_max"] radialScale, 17];
rMin = N[0.22 a, 17];
dr = (rOuter - rMin)/(nr + 1);
dw = L/(nw + 1);
rCenters = N[Table[rMin + i dr, {i, 1, nr}], 17];
wCenters = N[Table[i dw, {i, 1, nw}], 17];

wOld = N[background["grid"]["w_centers"], 17];
zOld = N[background["derived"]["Z_w"], 17];
r0Old = N[background["fields"]["R0_w"], 17];
A00 = N[Flatten[background["fields"]["A_00"]], 17];
Ar0 = N[Flatten[background["fields"]["A_r0"]], 17];
Aw0 = N[Flatten[background["fields"]["A_w0"]], 17];
zCenters = N[interp1Clamped[wOld, zOld, #] & /@ wCenters, 17];
r0Centers = N[interp1Clamped[wOld, r0Old, #] & /@ wCenters, 17];

packet = <|
  "schema" -> "stage1_patha_b2b_spike1_vsh_operator_ledger/v1",
  "scope" -> "Path-A B2b Spike-1 input ledger; full Green transfer is Spike-2",
  "background_path" -> backgroundPath,
  "background_content_hash" -> background["content_hash"],
  "tau" -> toJsonNumber[constants["tau"]],
  "geometry" -> <|
    "a" -> toJsonNumber[a],
    "L" -> toJsonNumber[L],
    "radial_scale" -> toJsonNumber[radialScale],
    "r_min_vsh" -> toJsonNumber[rMin],
    "r_outer_vsh" -> toJsonNumber[rOuter]
  |>,
  "grid" -> <|
    "nr" -> nr,
    "nw" -> nw,
    "r_centers" -> toJsonVector[rCenters],
    "w_centers" -> toJsonVector[wCenters],
    "Z_w_on_vsh_grid" -> toJsonVector[zCenters],
    "R0_w_on_vsh_grid" -> toJsonVector[r0Centers]
  |>,
  "A0_linf" -> <|
    "A_00" -> toJsonNumber[Max[Abs[A00]]],
    "A_r0" -> toJsonNumber[Max[Abs[Ar0]]],
    "A_w0" -> toJsonNumber[Max[Abs[Aw0]]]
  |>,
  "operator_form" -> <|
    "ell" -> 2,
    "Lharm" -> 6,
    "lanes" -> {"phi", "a_r", "a_E", "a_B", "a_w"},
    "gauge_weight" -> "H=Z_w from Path-A background",
    "smoke_geometry_used" -> False
  |>,
  "term_map" -> {
    <|"term" -> "VSH lane split", "implementation" -> "phi,a_r,a_E,a_B,a_w with angular weight L=6", "source" -> "decision 05 D1"|>,
    <|"term" -> "Path-A background", "implementation" -> "reads A_00,A_r0,A_w0,Z_w,R0_w from B2a JSON", "source" -> "pathA_11 B2b requirement 1"|>,
    <|"term" -> "frozen geometry", "implementation" -> "a=1,L=37/20", "source" -> "pathA Gate-A freeze / B2a exporter"|>
  }
|>;

outDir = DirectoryName[outPath];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
Export[outPath, packet, "JSON"];

Print["patha_b2b_spike1_operator_ledger: ", outPath];
Print["background_hash: ", background["content_hash"]];
Print["geometry a/L: ", a, " / ", L];
Print["A0 linf: ", packet["A0_linf"]];
