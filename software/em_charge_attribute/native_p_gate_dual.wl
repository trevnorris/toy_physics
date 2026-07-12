(* Genuine independent Wolfram Language engine for the native-P gate.          *)
(* Every native/control input follows buildH2 -> diracSearch -> searchG.       *)

ClearAll[fail, check, pb, independentIndices, buildH2, diracSearch, searchG,
  executeModel, nativeModel, maxwellModel, gaugedModel, sigmaModel,
  gaugeFixedModel, globalModel, couplingGuard, tunedRejectionGuard, nativeAnalysis];

fail[msg_] := (Print["FAIL: " <> ToString[msg, InputForm]]; Exit[1]);
check[test_, msg_] := If[! TrueQ[test], fail[msg]];
pb[f_, g_, pairs_List] := Factor[Total[(D[f, #[[1]]] D[g, #[[2]]] -
      D[f, #[[2]]] D[g, #[[1]]]) & /@ pairs]];
str[x_] := ToString[Factor[x], InputForm];
strList[x_List] := str /@ x;
strMatrix[x_List] := Map[str, x, {2}];
kvec = {1, 2, 3}; ksq = kvec.kvec;

independentIndices[rows_List] := Module[{chosen = {}, rank = 0, trial},
  Do[trial = Append[chosen, i];
    If[MatrixRank[rows[[trial]]] > rank, chosen = trial; rank++],
    {i, Length[rows]}]; chosen];

buildH2[model_Association] := Module[
  {q, v, lag, n, mom, momMap, hessian, nulls, residual, primaries, rank, safe,
   rowIds, colIds, minor, rhs, solved, rules, hc, pairs},
  q = model["q"]; v = model["v"]; lag = model["L"]; n = Length[q];
  safe = StringReplace[model["key"], "_" -> ""];
  mom = Table[Unique["Pi" <> safe <> "$"] , {n}];
  momMap = D[lag, #] & /@ v;
  hessian = Table[D[momMap[[i]], v[[j]]], {i, n}, {j, n}];
  nulls = NullSpace[Transpose[hessian]];
  residual = mom - (momMap /. Thread[v -> ConstantArray[0, n]]);
  primaries = Factor[#.residual] & /@ nulls;
  rank = MatrixRank[hessian];
  rowIds = independentIndices[hessian];
  colIds = independentIndices[Transpose[hessian]];
  check[Length[rowIds] == rank && Length[colIds] == rank,
    model["key"] <> " Hessian minor selection failed"];
  rules = Thread[v -> ConstantArray[0, n]];
  If[rank > 0,
    minor = hessian[[rowIds, colIds]];
    rhs = residual[[rowIds]];
    solved = Inverse[minor].rhs;
    Do[rules[[colIds[[a]]]] = v[[colIds[[a]]]] -> solved[[a]], {a, rank}]
  ];
  (* ReplaceAll uses the last duplicate rule, so pivot rules override zeros. *)
  hc = Factor[(mom.v - lag) /. rules];
  check[FreeQ[hc, Alternatives @@ v], model["key"] <> " retained a velocity"];
  pairs = Transpose[{q, mom}];
  <|"model" -> model, "momenta" -> mom, "momentum_map" -> momMap,
    "hessian" -> hessian, "hessian_rank" -> rank,
    "hessian_nullspace" -> nulls, "primary_constraints" -> primaries,
    "H2" -> hc, "pairs" -> pairs|>
];

addsIndependent[expr_, constraints_List, z_List] := Module[{row, old},
  row = {D[expr, #] & /@ z};
  If[! AnyTrue[First[row], Simplify[# != 0] &], Return[False]];
  (* Explicit table avoids evaluation ambiguities in nested pure functions. *)
  old = Table[D[constraints[[i]], z[[j]]], {i, Length[constraints]}, {j, Length[z]}];
  MatrixRank[Join[old, row]] > MatrixRank[old]
];

diracSearch[build_Association] := Module[
  {constraints, stages, primaryCount, z, bmat, hflow, compat, new,
   candidate, m, rank, kernel, first, classes, stage},
  constraints = build["primary_constraints"]; stages = ConstantArray[0, Length[constraints]];
  primaryCount = Length[constraints]; z = Join[build["model"]["q"], build["momenta"]];
  For[stage = 1, stage <= 2 Length[z], stage++,
    If[Length[constraints] == 0, Break[]];
    bmat = Table[pb[constraints[[i]], build["primary_constraints"][[j]], build["pairs"]],
      {i, Length[constraints]}, {j, primaryCount}];
    hflow = pb[#, build["H2"], build["pairs"]] & /@ constraints;
    compat = NullSpace[Transpose[bmat]]; new = {};
    Do[candidate = Factor[left.hflow];
      If[addsIndependent[candidate, Join[constraints, new], z], AppendTo[new, candidate]],
      {left, compat}];
    If[Length[new] == 0, Break[]];
    constraints = Join[constraints, new]; stages = Join[stages, ConstantArray[stage, Length[new]]];
  ];
  check[stage <= 2 Length[z], build["model"]["key"] <> " Dirac iteration did not close"];
  m = Table[pb[constraints[[i]], constraints[[j]], build["pairs"]],
    {i, Length[constraints]}, {j, Length[constraints]}];
  rank = If[Length[constraints] == 0, 0, MatrixRank[m]];
  kernel = If[Length[constraints] == 0, {}, NullSpace[m]]; first = Length[kernel];
  check[first == Length[constraints] - rank, build["model"]["key"] <> " rank/nullity mismatch"];
  classes = Table[If[AllTrue[m[[i]], TrueQ[Factor[#] === 0] &], "FIRST_CLASS", "SECOND_CLASS_COMPONENT"],
    {i, Length[constraints]}];
  <|"constraints" -> constraints, "constraint_stages" -> stages,
    "primary_count" -> primaryCount, "bracket_matrix" -> m, "bracket_rank" -> rank,
    "first_class_count" -> first, "second_class_count" -> rank,
    "kernel" -> kernel, "constraint_classes_list" -> classes|>
];

proportionalK[action_List] := Length[action] == 3 && ! AllTrue[action, Simplify[# == 0] &] &&
  And @@ Flatten[Table[Simplify[action[[i]] kvec[[j]] - action[[j]] kvec[[i]]] == 0,
    {i, 3}, {j, i + 1, 3}]];

searchG[build_Association, d_Association] := Module[
  {constraints, m, r, primaryFC, accepted = {}, rejected = {}, coeff, primary,
   descendant, paction, saction, fcdesc, scalar, onlyscalar, spatialBlocks, spatial,
   nontrivial, desczero, nongradient, certified, reason, item},
  constraints = d["constraints"]; m = d["bracket_matrix"]; r = d["primary_count"];
  primaryFC = If[Length[constraints] == 0 || r == 0, {}, NullSpace[m[[All, 1 ;; r]]]];
  Do[
    primary = Factor[Sum[coeff[[i]] constraints[[i]], {i, r}]];
    descendant = Factor[pb[primary, build["H2"], build["pairs"]]];
    paction = pb[#, primary, build["pairs"]] & /@ build["model"]["q"];
    saction = pb[#, descendant, build["pairs"]] & /@ build["model"]["q"];
    fcdesc = AllTrue[constraints, Simplify[pb[descendant, #, build["pairs"]] == 0] &];
    scalar = AnyTrue[build["model"]["scalar_connection"], Simplify[paction[[#]] != 0] &];
    onlyscalar = scalar && And @@ Table[
      If[MemberQ[Join[build["model"]["scalar_connection"], build["model"]["multipliers"]], i],
        True, Simplify[paction[[i]] == 0]], {i, Length[paction]}];
    spatialBlocks = (<|"field_indices" -> #,
        "action" -> strList[saction[[#]]],
        "proportional_to_k" -> proportionalK[saction[[#]]]|> &) /@
      build["model"]["spatial_connections"];
    spatial = AnyTrue[spatialBlocks, TrueQ[#["proportional_to_k"]] &];
    nontrivial = AnyTrue[saction, Simplify[# != 0] &];
    desczero = TrueQ[Simplify[descendant == 0]];
    nongradient = ! TrueQ[spatial]; certified = TrueQ[desczero || nongradient];
    reason = Which[desczero, "DESCENDANT_ZERO", nongradient,
      "SECONDARY_ACTION_NON_GRADIENT", True, "OTHER_GAUSS_CRITERION"];
    item = <|"primary" -> str[primary], "descendant" -> str[descendant],
      "primary_action" -> strList[paction], "secondary_action" -> strList[saction],
      "spatial_action_blocks" -> spatialBlocks,
      "first_class_descendant" -> fcdesc, "scalar_primary_end" -> onlyscalar,
      "spatial_gradient_action" -> spatial, "descendant_zero" -> desczero,
      "secondary_action_non_gradient" -> nongradient,
      "descendant_rejection_certified" -> certified,
      "computed_rejection_reason" -> reason|>;
    If[TrueQ[fcdesc && onlyscalar && spatial && nontrivial], AppendTo[accepted, item], AppendTo[rejected, item]],
    {coeff, primaryFC}];
  <|"computed_kernel_dimension" -> Length[d["kernel"]],
    "computed_first_class_primaries" -> Length[primaryFC],
    "gauss_candidates" -> Length[accepted], "additional_G_exists" -> (Length[accepted] > 0),
    "candidates" -> accepted, "rejected_first_class_primaries" -> rejected|>
];

tunedRejectionGuard[theory_String, condition_String, d_Association, g_Association] := Module[
  {fc, rejected, candidates, primaryFC, records},
  fc = d["first_class_count"];
  If[fc == 0, Return[<|"status" -> "NOT_APPLICABLE_FC0", "computed" -> True,
    "checked_directions" -> 0|>]];
  rejected = g["rejected_first_class_primaries"]; candidates = g["candidates"];
  primaryFC = g["computed_first_class_primaries"];
  check[primaryFC == Length[rejected] + Length[candidates],
    "HARDENING-TUNED-DESCENDANT-REJECTION THEORY-" <> theory <> " " <> condition <>
      ": first-class primary accounting mismatch"];
  Do[check[TrueQ[item["descendant_rejection_certified"]],
      "HARDENING-TUNED-DESCENDANT-REJECTION THEORY-" <> theory <> " " <> condition <>
        " direction " <> ToString[i] <>
        ": rejected primary has a nonzero proportional-to-k descendant"],
    {i, Length[rejected]}, {item, {rejected[[i]]}}];
  records = MapIndexed[<|"direction" -> First[#2], "primary" -> #1["primary"],
      "descendant" -> #1["descendant"], "secondary_action" -> #1["secondary_action"],
      "spatial_action_blocks" -> #1["spatial_action_blocks"],
      "descendant_zero" -> #1["descendant_zero"],
      "secondary_action_non_gradient" -> #1["secondary_action_non_gradient"],
      "computed_rejection_reason" -> #1["computed_rejection_reason"]|> &, rejected];
  <|"status" -> "PASS", "computed" -> True, "condition" -> condition,
    "first_class_count" -> fc, "first_class_primaries" -> primaryFC,
    "accepted_gauss_directions" -> Length[candidates],
    "checked_directions" -> Length[records], "records" -> records|>
];

serialize[build_, d_, g_] := <|
  "model" -> build["model"]["key"], "lagrangian" -> str[build["model"]["L"]],
  "coordinates" -> strList[build["model"]["q"]], "momentum_map" -> strList[build["momentum_map"]],
  "hessian" -> strMatrix[build["hessian"]], "hessian_rank" -> build["hessian_rank"],
  "hessian_nullity" -> Length[build["hessian_nullspace"]],
  "primary_constraints" -> strList[build["primary_constraints"]], "H2" -> str[build["H2"]],
  "constraints" -> strList[d["constraints"]], "constraint_stages" -> d["constraint_stages"],
  "constraint_classes" -> d["constraint_classes_list"],
  "class_signature" -> Counts[d["constraint_classes_list"]],
  "bracket_matrix" -> strMatrix[d["bracket_matrix"]], "bracket_rank" -> d["bracket_rank"],
  "first_class_count" -> d["first_class_count"], "second_class_count" -> d["second_class_count"],
  "G_search" -> g|>;

executeModel[model_] := Module[{b = buildH2[model], d, g},
  d = diracSearch[b]; g = searchG[b, d]; {b, d, g, serialize[b, d, g]}];

nativeModel[theory_String, boundary_: 0, tuning_: "generic", extraRules_: {}] := Module[
  {p, u, s, lam, sigma, b, xi, q, v, vp, vu, vs, vb, gt, gs, gd, gb,
   ru, ap, au, ab, trans, kinetic, potential, lag, couplings},
  p = Array[Symbol["p" <> theory <> ToString[#]] &, 3];
  u = Array[Symbol["u" <> theory <> ToString[#]] &, 3];
  s = Symbol["s" <> theory]; lam = Symbol["lambda" <> theory];
  sigma = Symbol["sigma" <> theory]; b = Symbol["b" <> theory]; xi = Symbol["xi" <> theory];
  q = Join[p, u, {s, lam, sigma, b}, If[theory == "C", {xi}, {}]];
  v = Symbol["d" <> SymbolName[#]] & /@ q; vp = v[[1 ;; 3]]; vu = v[[4 ;; 6]];
  vs = v[[7]]; vb = v[[10]];
  gt = Symbol["gt" <> theory]; gs = Symbol["gs" <> theory];
  gd = Symbol["gd" <> theory]; gb = Symbol["gb" <> theory];
  ru = Symbol["rhou" <> theory]; ap = Symbol["ap" <> theory];
  au = Symbol["au" <> theory]; ab = Symbol["ab" <> theory];
  trans = ksq IdentityMatrix[3] - Outer[Times, kvec, kvec];
  kinetic = (vp.vp + vs^2 + ru vu.vu)/2 + gt vp.vu + If[theory == "A", vb^2/2, 0];
  potential = ap p.p/2 + au u.trans.u/2 + gs ksq p.u + gd (kvec.p)^2/2 + ab b^2/2 +
    If[theory == "A", gb b kvec.p, 0];
  lag = kinetic - potential + lam s + sigma kvec.u + If[theory == "C", xi b, 0];
  couplings = Join[{gt, gs, gd}, If[theory == "A", {gb}, {}]];
  If[boundary != 0,
    lag = Factor[lag /. Switch[tuning,
      "generic", {gt -> boundary, ru -> 1},
      "rankdrop", {gt -> boundary, ru -> 1, gs -> boundary (ap + 14 au)/28},
      "common", {gt -> boundary, ru -> 1, ap -> 14 au, gs -> boundary au}]];
    couplings = If[tuning == "generic", Rest[couplings], Drop[couplings, 2]]
  ];
  If[Length[extraRules] > 0,
    lag = Factor[lag /. extraRules];
    couplings = Select[couplings, ! MemberQ[First /@ extraRules, #] &]
  ];
  <|"key" -> "native_" <> theory, "q" -> q, "v" -> v, "L" -> lag,
    "couplings" -> couplings, "scalar_connection" -> {},
    "spatial_connections" -> {{1, 2, 3}, {4, 5, 6}},
    "multipliers" -> Join[{8, 9}, If[theory == "C", {11}, {}]]|>
];

maxwellModel[key_String : "maxwell", kineticA0_: False, defect_: False] := Module[
  {q, v, av, va, electric, magnetic, lag, safe, currentData = <||>},
  safe = StringReplace[key, "_" -> ""];
  q = Array[Symbol[safe <> "A" <> ToString[# - 1]] &, 4]; v = Symbol["d" <> SymbolName[#]] & /@ q;
  av = q[[2 ;; 4]]; va = v[[2 ;; 4]]; electric = va - kvec q[[1]];
  magnetic = ksq IdentityMatrix[3] - Outer[Times, kvec, kvec];
  lag = electric.electric/2 - av.magnetic.av/2 + If[kineticA0, v[[1]]^2/2, 0];
  If[key == "nonconserved_current",
    lag = lag + ncrho q[[1]] + {ncj1, ncj2, ncj3}.av;
    currentData = <|"rho" -> ncrho, "rho_dot" -> ncrhodot,
      "j" -> {ncj1, ncj2, ncj3}, "impose_conservation" -> (! defect)|>
  ];
  <|"key" -> key, "q" -> q, "v" -> v, "L" -> lag, "couplings" -> {},
    "scalar_connection" -> {1}, "spatial_connections" -> {{2, 3, 4}}, "multipliers" -> {},
    "current_data" -> currentData|>
];

gaugedModel[removeHard_: False] := Module[{base, q, v, av, va, phi, vf, rot, electric, magnetic, lag},
  base = maxwellModel["gaugedhard"]; phi = {ghphi1, ghphi2};
  q = Join[base["q"], phi, If[removeHard, {}, {ghlambda}]]; v = Symbol["d" <> SymbolName[#]] & /@ q;
  av = q[[2 ;; 4]]; va = v[[2 ;; 4]]; electric = va - kvec q[[1]];
  magnetic = ksq IdentityMatrix[3] - Outer[Times, kvec, kvec]; vf = v[[5 ;; 6]]; rot = {-phi[[2]], phi[[1]]};
  lag = electric.electric/2 - av.magnetic.av/2 + (vf - q[[1]] rot).(vf - q[[1]] rot)/2 - phi.phi/2 +
    If[removeHard, 0, ghlambda (phi.phi - 1)/2];
  <|"key" -> "gauged_hard_unit", "q" -> q, "v" -> v, "L" -> lag, "couplings" -> {},
    "scalar_connection" -> {1}, "spatial_connections" -> {{2, 3, 4}},
    "multipliers" -> If[removeHard, {}, {7}]|>
];

sigmaModel[removeHard_: False] := Module[{phi, q, v, lag},
  phi = {bsphi1, bsphi2, bsphi3, bsphi4}; q = Join[phi, If[removeHard, {}, {bslambda}]];
  v = Symbol["d" <> SymbolName[#]] & /@ q; lag = v[[1 ;; 4]].v[[1 ;; 4]]/2 - phi.phi/2 +
    If[removeHard, 0, bslambda (phi.phi - 1)/2];
  <|"key" -> "bare_sigma", "q" -> q, "v" -> v, "L" -> lag, "couplings" -> {},
    "scalar_connection" -> {}, "spatial_connections" -> {}, "multipliers" -> If[removeHard, {}, {5}]|>
];

gaugeFixedModel[] := Module[{q, v, av, va, electric, magnetic, lag},
  q = {gfA0, gfA1, gfA2, gfA3, gfeta, gfzeta}; v = Symbol["d" <> SymbolName[#]] & /@ q;
  av = q[[2 ;; 4]]; va = v[[2 ;; 4]]; electric = va - kvec q[[1]];
  magnetic = ksq IdentityMatrix[3] - Outer[Times, kvec, kvec];
  lag = electric.electric/2 - av.magnetic.av/2 + gfeta gfA0 + gfzeta kvec.av;
  <|"key" -> "gauge_fixed_maxwell", "q" -> q, "v" -> v, "L" -> lag, "couplings" -> {},
    "scalar_connection" -> {1}, "spatial_connections" -> {{2, 3, 4}}, "multipliers" -> {5, 6}|>
];

globalModel[] := <|"key" -> "global_only", "q" -> {globalx, globaly}, "v" -> {dglobalx, dglobaly},
  "L" -> (dglobalx^2 + dglobaly^2)/2 - (globalx^2 + globaly^2)/2, "couplings" -> {},
  "scalar_connection" -> {}, "spatial_connections" -> {}, "multipliers" -> {}|>;

couplingGuard[build_, d_] := Module[{locations = <||>, hits, c},
  Do[hits = {};
    If[! FreeQ[build["momentum_map"], c], AppendTo[hits, "momentum_map"]];
    If[! FreeQ[build["hessian"], c], AppendTo[hits, "hessian"]];
    If[! FreeQ[d["constraints"], c], AppendTo[hits, "constraints"]];
    If[! FreeQ[d["bracket_matrix"], c], AppendTo[hits, "bracket_matrix"]];
    check[Length[hits] > 0, "GUARD-COUPLINGS-ENTER dropped " <> SymbolName[c]];
    AssociateTo[locations, SymbolName[c] -> hits], {c, build["model"]["couplings"]}];
  <|"status" -> "PASS", "computed" -> True, "locations" -> locations|>
];

nativeAnalysis[theory_String] := Module[
  {run, b, d, g, data, guard, boundary, br, anyg, scan, condition, rejection,
   sweep, sweepSeed, sample, auv, apv, extra, tuned, tunedGuard, physicalDet,
   componentDet, gt, ru, ap, au, gd, gb, scanned},
  run = executeModel[nativeModel[theory]]; {b, d, g, data} = run; guard = couplingGuard[b, d];
  gt = Symbol["gt" <> theory]; ru = Symbol["rhou" <> theory];
  ap = Symbol["ap" <> theory]; au = Symbol["au" <> theory];
  gd = Symbol["gd" <> theory]; gb = Symbol["gb" <> theory];
  physicalDet = Factor[Det[b["hessian"][[1 ;; 6, 1 ;; 6]]]];
  componentDet = Factor[Det[b["hessian"][[{1, 4}, {1, 4}]]]];
  check[TrueQ[Simplify[physicalDet == (ru - gt^2)^3]],
    "THEORY-" <> theory <> ": unexpected physical kinetic-Hessian degeneracy"];
  check[TrueQ[Simplify[componentDet == ru - gt^2]],
    "THEORY-" <> theory <> ": kinetic component determinant mismatch"];
  boundary = Flatten[Table[Table[
    br = executeModel[nativeModel[theory, sign, scan]];
    condition = "rhou=1,gt=" <> ToString[sign] <> ";" <> scan;
    rejection = tunedRejectionGuard[theory, condition, br[[2]], br[[3]]];
    <|"condition" -> condition,
      "coordinates" -> strList[br[[1]]["model"]["q"]],
      "first_class_count" -> br[[2]]["first_class_count"],
      "first_class_primaries" -> br[[3]]["computed_first_class_primaries"],
      "gauss_candidates" -> br[[3]]["gauss_candidates"],
      "additional_G_exists" -> br[[3]]["additional_G_exists"],
      "hessian_nullity" -> br[[4]]["hessian_nullity"],
      "rejected_first_class_primaries" -> br[[3]]["rejected_first_class_primaries"],
      "tuned_descendant_rejection_guard" -> rejection|>,
    {scan, {"generic", "rankdrop", "common"}}], {sign, {-1, 1}}], 1];
  sweepSeed = 260712 + If[theory == "A", 1, 3];
  sweep = BlockRandom[SeedRandom[sweepSeed]; Flatten[Table[Table[
    auv = RandomInteger[{1, 5}]; apv = RandomInteger[{1, 40}];
    While[apv == 14 auv, apv = RandomInteger[{1, 40}]];
    extra = {ap -> apv, au -> auv, gd -> 0}; If[theory == "A", AppendTo[extra, gb -> 0]];
    br = executeModel[nativeModel[theory, sign, "rankdrop", extra]];
    condition = "rank-drop randomized sample sign=" <> ToString[sign] <>
      " index=" <> ToString[sample];
    rejection = tunedRejectionGuard[theory, condition, br[[2]], br[[3]]];
    <|"condition" -> condition,
      "sample_values" -> <|"a_p" -> apv, "a_u" -> auv, "g_d" -> 0,
        "g_b" -> If[theory == "A", 0, "not-applicable"]|>,
      "first_class_count" -> br[[2]]["first_class_count"],
      "first_class_primaries" -> br[[3]]["computed_first_class_primaries"],
      "gauss_candidates" -> br[[3]]["gauss_candidates"],
      "additional_G_exists" -> br[[3]]["additional_G_exists"],
      "hessian_nullity" -> br[[4]]["hessian_nullity"],
      "tuned_descendant_rejection_guard" -> rejection|>,
    {sample, 3}], {sign, {-1, 1}}], 1]];
  scanned = Join[boundary, sweep];
  tuned = Select[#["tuned_descendant_rejection_guard"] & /@ scanned,
    #["status"] == "PASS" &];
  check[Length[tuned] > 0, "THEORY-" <> theory <> ": tuned FC hardening audit empty"];
  tunedGuard = <|"status" -> "PASS", "computed" -> True,
    "checked_strata" -> Length[tuned],
    "checked_directions" -> Total[#["checked_directions"] & /@ tuned],
    "records" -> tuned|>;
  anyg = TrueQ[g["additional_G_exists"]] || AnyTrue[scanned, TrueQ[#["additional_G_exists"]] &];
  data = Join[data, <|
    "couplings" -> (SymbolName /@ b["model"]["couplings"]), "coupling_guard" -> guard,
    "tuned_descendant_rejection_guard" -> tunedGuard,
    "coupling_scan" -> <|"physical_kinetic_hessian_determinant" -> str[physicalDet],
      "physical_kinetic_determinant_per_component" -> str[componentDet],
      "only_kinetic_hessian_degeneracy" -> str[componentDet] <> "=0",
      "regular_result" -> <|"first_class_count" -> d["first_class_count"], "gauss_candidates" -> g["gauss_candidates"],
        "additional_G_exists" -> g["additional_G_exists"]|>, "semidefinite_boundary" -> boundary,
      "rankdrop_representative_sweep" -> <|
        "method" -> "fixed-seed randomized representative points on the solved non-common rank-drop locus",
        "seed" -> sweepSeed, "sample_count" -> Length[sweep], "samples" -> sweep|>|>,
    "decision_order" -> "quadratic (field-amplitude order 2)",
    "verdict" -> If[anyg, If[TrueQ[g["additional_G_exists"]], "FIRST_CLASS_GENERIC_EM_CANDIDATE", "FIRST_CLASS_TUNED_INVERSE_DESIGN"], "NATIVE_P_NO_EMERGENT_GAUSS"]|>];
  data
];

controlModels = <|"maxwell" -> maxwellModel[], "gauged_hard_unit" -> gaugedModel[],
  "bare_sigma" -> sigmaModel[], "nonconserved_current" -> maxwellModel["nonconserved_current", False, True],
  "gauge_fixed_maxwell" -> gaugeFixedModel[], "global_only" -> globalModel[]|>;
controlRuns = Map[executeModel, controlModels];
controls = Map[Last, controlRuns];
ncRun = controlRuns["nonconserved_current"];
ncPrimaryDirections = NullSpace[ncRun[[2]]["bracket_matrix"][[All, 1 ;; ncRun[[2]]["primary_count"]]]];
check[Length[ncPrimaryDirections] > 0, "nonconserved control lost its first-class primary"];
ncPrimary = Factor[Sum[ncPrimaryDirections[[1, i]] ncRun[[2]]["constraints"][[i]],
  {i, ncRun[[2]]["primary_count"]}]];
ncGauss = Factor[pb[ncPrimary, ncRun[[1]]["H2"], ncRun[[1]]["pairs"]]];
ncPreservation = Factor[D[ncGauss, ncrho] ncrhodot + pb[ncGauss, ncRun[[1]]["H2"], ncRun[[1]]["pairs"]]];
check[! FreeQ[ncPreservation, ncrhodot] && ! TrueQ[ncPreservation === 0],
  "nonconserved current did not produce a continuity defect"];
controls["nonconserved_current"] = Join[controls["nonconserved_current"],
  <|"current_inconsistent" -> True,
    "current_preservation" -> <|"gauss" -> str[ncGauss], "raw" -> str[ncPreservation]|>|>];
classifications = <|
  "maxwell" -> If[controls["maxwell"]["G_search"]["gauss_candidates"] > 0 && controls["maxwell"]["second_class_count"] == 0, "FIRST_CLASS_GAUSS", "FAIL"],
  "gauged_hard_unit" -> If[controls["gauged_hard_unit"]["G_search"]["gauss_candidates"] > 0 && controls["gauged_hard_unit"]["first_class_count"] > 0 && controls["gauged_hard_unit"]["second_class_count"] > 0, "MIXED", "FAIL"],
  "bare_sigma" -> If[controls["bare_sigma"]["first_class_count"] == 0 && controls["bare_sigma"]["second_class_count"] > 0, "SECOND_CLASS_RADIAL_NO_GAUSS", "FAIL"],
  "nonconserved_current" -> If[controls["nonconserved_current"]["G_search"]["gauss_candidates"] > 0 && TrueQ[controls["nonconserved_current"]["current_inconsistent"]], "INCONSISTENT_PRESERVATION", "FAIL"],
  "gauge_fixed_maxwell" -> If[controls["gauge_fixed_maxwell"]["first_class_count"] == 0 && controls["gauge_fixed_maxwell"]["G_search"]["gauss_candidates"] == 0, "SECOND_CLASS_NO_LOCAL_GAUGE", "FAIL"],
  "global_only" -> If[controls["global_only"]["hessian_nullity"] == 0 && controls["global_only"]["first_class_count"] == 0, "GLOBAL_CHARGE_NO_LOCAL_GAUSS", "FAIL"]|>;
check[FreeQ[Values[classifications], "FAIL"], "a shared-path control failed"];
check[controls["maxwell"]["G_search"]["gauss_candidates"] > 0 &&
  controls["gauged_hard_unit"]["G_search"]["gauss_candidates"] > 0, "GUARD-SEARCH-CAPABLE failed"];

theoryA = nativeAnalysis["A"]; theoryC = nativeAnalysis["C"];
baseDir = DirectoryName[$InputFileName]; outDir = FileNameJoin[{baseDir, "reports", "native_p_gate_artifacts"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
result = <|"schema" -> "native_p_constraint_gate/mathematica/rebuild-v1", "engine" -> "WolframLanguage",
  "engine_version" -> ToString[$Version], "pipeline" -> "buildH2 -> diracSearch -> searchG",
  "theories" -> <|"A" -> theoryA, "C" -> theoryC|>, "controls" -> controls,
  "control_classifications" -> classifications,
  "guards" -> <|"GUARD-COUPLINGS-ENTER" -> <|"A" -> theoryA["coupling_guard"], "C" -> theoryC["coupling_guard"]|>,
    "GUARD-SEARCH-CAPABLE" -> <|"status" -> "PASS", "computed" -> True|>,
    "HARDENING-TUNED-DESCENDANT-REJECTION" -> <|
      "A" -> theoryA["tuned_descendant_rejection_guard"],
      "C" -> theoryC["tuned_descendant_rejection_guard"]|>|>, "algebra_status" -> "PASS"|>;
Export[FileNameJoin[{outDir, "mathematica_results.json"}], result, "RawJSON"];
Print["MATHEMATICA_GENUINE_QUADRATIC_DIRAC: PASS"];
Print["GUARD-COUPLINGS-ENTER: PASS"];
Print["GUARD-SEARCH-CAPABLE: PASS"];
Print["HARDENING-TUNED-DESCENDANT-REJECTION: PASS"];
Print["THEORY-A: FC=" <> ToString[theoryA["first_class_count"]] <> " G=" <> ToString[theoryA["G_search"]["gauss_candidates"]] <> " VERDICT=" <> theoryA["verdict"]];
Print["THEORY-C: FC=" <> ToString[theoryC["first_class_count"]] <> " G=" <> ToString[theoryC["G_search"]["gauss_candidates"]] <> " VERDICT=" <> theoryC["verdict"]];
Print["SIX_CONTROLS_SHARED_PIPELINE: PASS"];
Exit[0];
