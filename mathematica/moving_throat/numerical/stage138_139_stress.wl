(* Numerical stress harness for the Family-1 co-evolving mouth/core map (Stages 138-139). *)

ClearAll["Global`*"];

rootDir = DirectoryName[$InputFileName];
configPath = FileNameJoin[{
   rootDir, "..", "..", "..", "scripts", "moving_throat", "numerical",
   "stage138_139_config.json"
}];
cfg = Import[ExpandFileName[configPath], "RawJSON"];

tol = 10^-10;
momTol = 10^-10;
nearQ[lhs_, rhs_, t_: tol] := Abs[lhs - rhs] <= t (1 + Abs[rhs]);
fmt[x_] := ToString @ NumberForm[N[x, 14], {Infinity, 12}, ExponentFunction -> (Null &)];

require[label_, condition_, detail_] := Module[{status},
  status = If[TrueQ[condition], "PASS", "FAIL"];
  Print["[", status, "] ", label, ": ", detail];
  If[! TrueQ[condition], Throw[$Failed]]
];

gridPts = cfg["grid"]["points"];
xs = N[Range[0, gridPts - 1]/(gridPts - 1)];
weights = ConstantArray[1./(gridPts - 1), gridPts];
weights[[1]] /= 2; weights[[-1]] /= 2;

rF1 = cfg["constants"]["r_F1"];
gStar = cfg["constants"]["g_star"];
sigma0Star = cfg["frozen_traction"]["Sigma0"];
stage138Ref = cfg["stage138_reference"];
stage139Ref = cfg["stage139_reference"];
fixedPointTol = cfg["tolerances"]["fixed_point_tol"];
maxIter = cfg["tolerances"]["max_iter"];
relax = cfg["tolerances"]["relaxation"];
rootTol = cfg["tolerances"]["root_tol"];

cosVec = Cos[Pi xs/2];
sVec = Cosh[Pi (1 - xs)/2]/Cosh[Pi/2];
kernelTs = Outer[Min, xs, xs];
kernelTq = Outer[
  (Sinh[Pi Min[#1, #2]/2] Cosh[Pi (1 - Max[#1, #2])/2])/((Pi/2) Cosh[Pi/2]) &,
  xs, xs, 1
];

normalize[v_] := v/(Total[v weights]);
profileFromSeed[seed_] := Module[{name = seed["name"], vals},
  Switch[name,
    "canonical_exponential",
      vals = Exp[-seed["Pi"] xs],
    "uniform_broad",
      vals = ConstantArray[1., Length[xs]],
    "derivative_like",
      vals = Cos[Pi xs/2],
    "endpoint_biased",
      vals = Exp[-6. seed["mix"] xs] + 0.15,
    _,
      vals = ConstantArray[1., Length[xs]]
  ];
  normalize[vals]
];

moments[prof_] := <|
  "g" -> Total[prof cosVec weights],
  "S" -> Total[prof sVec weights]
|>;

applyKernel[kernel_, prof_] := kernel . (prof weights);

stepMap[prof_, sigma0_] := Module[{mom, g, S, R, Ts, Tq, phi, vals},
  mom = moments[prof];
  g = mom["g"]; S = mom["S"];
  R = (g - rF1)^2/(1 + rF1^2);
  Ts = applyKernel[kernelTs, prof];
  Tq = applyKernel[kernelTq, prof];
  phi = sigma0 (Ts - R Tq);
  vals = Exp[-phi];
  normalize[vals]
];

fixedPoint[seed_, sigma0_] := Module[{prof = profileFromSeed[seed], next, change, mom, lastMom = <||>, it},
  For[it = 1, it <= maxIter, it++,
    next = (1 - relax) prof + relax stepMap[prof, sigma0];
    next = normalize[next];
    change = Max[Abs[next - prof]];
    mom = moments[next];
    If[change < fixedPointTol && (lastMom === <||> || Max[Abs[{mom["g"], mom["S"]} - {lastMom["g"], lastMom["S"]}]] < momTol),
      Return[<|"profile" -> next, "iterations" -> it, "g" -> mom["g"], "S" -> mom["S"], "R" -> (mom["g"] - rF1)^2/(1 + rF1^2)|>]
    ];
    prof = next;
    lastMom = mom;
  ];
  Throw[$Failed]
];

Print["Loaded config from ", configPath];
Print["Grid points: ", gridPts];

frozenSeeds = cfg["frozen_traction"]["seeds"];
frozenFP = Table[fixedPoint[seed, sigma0Star], {seed, frozenSeeds}];

Print[""];
Print["=== Frozen traction multi-seed solve ==="];
Do[
  Print[
    "[PASS] seed ", seed["name"], ": iterations=", fp["iterations"],
    ", g_fp=", fmt[fp["g"]], ", S_fp=", fmt[fp["S"]], ", R_fp=", fmt[fp["R"]]
  ],
  {seed, frozenSeeds}, {fp, {frozenFP[[First@FirstPosition[frozenSeeds, seed]]}}]
];

gVals = frozenFP[[All, "g"]];
sVals = frozenFP[[All, "S"]];
rVals = frozenFP[[All, "R"]];
gAvg = Mean[gVals];
sAvg = Mean[sVals];
rAvg = Mean[rVals];
Scan[require["frozen-traction g consistency", nearQ[#, gAvg, momTol], "g_fp=" <> fmt[#] <> ", mean=" <> fmt[gAvg]] &, gVals];
Scan[require["frozen-traction S consistency", nearQ[#, sAvg, momTol], "S_fp=" <> fmt[#] <> ", mean=" <> fmt[sAvg]] &, sVals];
Scan[require["frozen-traction R consistency", nearQ[#, rAvg, momTol], "R_fp=" <> fmt[#] <> ", mean=" <> fmt[rAvg]] &, rVals];
require["frozen traction stays below compensated target", gAvg < gStar, "g_fp(mean)=" <> fmt[gAvg] <> ", target=" <> fmt[gStar]];
require["frozen traction remains above compensated R threshold", rAvg > 1/4, "R_fp(mean)=" <> fmt[rAvg]];
require["frozen traction matches Stage 138 reference g", nearQ[gAvg, stage138Ref["g_fp"], 10^-5], "g_fp(mean)=" <> fmt[gAvg] <> ", ref=" <> fmt[stage138Ref["g_fp"]]];
require["frozen traction matches Stage 138 reference S", nearQ[sAvg, stage138Ref["S_fp"], 10^-5], "S_fp(mean)=" <> fmt[sAvg] <> ", ref=" <> fmt[stage138Ref["S_fp"]]];
require["frozen traction matches Stage 138 reference R", nearQ[rAvg, stage138Ref["R_fp"], 10^-5], "R_fp(mean)=" <> fmt[rAvg] <> ", ref=" <> fmt[stage138Ref["R_fp"]]];

rootSeeds = cfg["root_search"]["seeds"];
rootSeedFP[sigma0_] := Mean[(fixedPoint[#, sigma0]["g"] & /@ rootSeeds)];

Print[""];
Print["=== Root search samples ==="];
samples = Table[{sigma0, rootSeedFP[sigma0]}, {sigma0, cfg["root_search"]["samples"]}];
Scan[(Print["sigma0=", fmt[#[[1]]], " -> g_fp(mean)=", fmt[#[[2]]]]) &, samples];
Do[
  require[
    "sample monotonicity on analyzed interval",
    samples[[i + 1, 2]] > samples[[i, 2]],
    "g_left=" <> fmt[samples[[i, 2]]] <> ", g_right=" <> fmt[samples[[i + 1, 2]]]
  ],
  {i, 1, Length[samples] - 1}
];

low = First@First@Select[Partition[samples, 2, 1], #[[1, 2]] < gStar < #[[2, 2]] &];
high = Last@First@Select[Partition[samples, 2, 1], #[[1, 2]] < gStar < #[[2, 2]] &];
require["root bracket order", low < high, "lo=" <> fmt[low] <> ", hi=" <> fmt[high]];

target = cfg["root_search"]["target_g"];
left = low; right = high;
Do[
  mid = (left + right)/2;
  midVal = rootSeedFP[mid];
  If[(rootSeedFP[left] - target) (midVal - target) <= 0,
    right = mid,
    left = mid
  ];
  If[Abs[right - left] < rootTol, Break[]],
  {50}
];

sigma0Can = (left + right)/2;
rootFP = fixedPoint[#, sigma0Can] & /@ rootSeeds;
rootG = Mean[rootFP[[All, "g"]]];
rootS = Mean[rootFP[[All, "S"]]];
rootR = Mean[rootFP[[All, "R"]]];
Scan[require["root-seed consistency", nearQ[#[["g"]], rootG, momTol], "g_fp=" <> fmt[#[["g"]]] <> ", mean=" <> fmt[rootG]] &, rootFP];

Print[""];
Print["=== Renormalized canonical solve ==="];
Print["sigma0_can = ", fmt[sigma0Can]];
Print["g_can(mean) = ", fmt[rootG]];
Print["S_can(mean) = ", fmt[rootS]];
Print["R_can(mean) = ", fmt[rootR]];
require["root moment close to target", nearQ[rootG, target, 10^-7], "g_can=" <> fmt[rootG] <> ", target=" <> fmt[target]];
require["compensated branch recovered", nearQ[rootR, 1/4, 10^-7], "R_can=" <> fmt[rootR]];
require["root traction matches Stage 139 reference", nearQ[sigma0Can, stage139Ref["Sigma0_can"], 10^-4], "sigma0_can=" <> fmt[sigma0Can] <> ", ref=" <> fmt[stage139Ref["Sigma0_can"]]];
require["root mouth traction matches Stage 139 reference", nearQ[tmCan, stage139Ref["Tm_can"], 10^-4], "Tm_can=" <> fmt[tmCan] <> ", ref=" <> fmt[stage139Ref["Tm_can"]]];

tmCan = Sqrt[9 sigma0Can/20];
Print["Tm_can = ", fmt[tmCan]];
Print["Pi_can estimated = ", fmt[sigma0Can (1 - rootR rootS)]];

Print[""];
Print["All Family-1 co-evolution stress checks passed."];
