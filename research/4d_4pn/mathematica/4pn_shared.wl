If[! ValueQ[$FourPNSharedLoaded],
  $FourPNSharedLoaded = True;

  banner[title_] := Module[{line = StringRepeat["=", 88]},
    Print[""];
    Print[line];
    Print[title];
    Print[line];
  ];

  subbanner[title_] := Module[{line = StringRepeat["-", 88]},
    Print[""];
    Print[line];
    Print[title];
    Print[line];
  ];

  checkZero[name_, expr_] := Module[{res = FullSimplify[Expand[expr]]},
    Print[name, " = ", res];
    If[res === 0,
      passCount++,
      failCount++;
      Print["FAIL: ", name];
      Throw[$Failed]
    ];
  ];

  checkEqual[name_, lhs_, rhs_] := checkZero[name, lhs - rhs];

  checkTrue[name_, cond_] := Module[{ok = TrueQ[cond]},
    Print[name, " = ", ok];
    If[ok,
      passCount++,
      failCount++;
      Print["FAIL: ", name];
      Throw[$Failed]
    ];
  ];

  finalizeAudit[] := (
    Print["Passes: ", passCount];
    Print["Fails: ", failCount];
    If[failCount != 0, Exit[1]];
  );

  aa = Symbol["aa"];
  bb = Symbol["bb"];
  cc = Symbol["cc"];
  dd = Symbol["dd"];
  ee = Symbol["ee"];
  pp = Symbol["pp"];
  qq = Symbol["qq"];
  V2sym = Symbol["V2"];
  rdSym = Symbol["rd"];
  nuSym = Symbol["nu"];
  delSym = Symbol["del4pn"];

  Xa4PN = (1 + delSym)/2;
  Xb4PN = (1 - delSym)/2;

  blockRanges4PN = <|
    "Q" -> Range[7, 11],
    "T" -> Range[12, 15],
    "S" -> Range[16, 18],
    "U" -> Range[19, 20],
    "W" -> {21}
  |>;

  comSubs4PN = {
    pp -> Xa4PN,
    qq -> Xb4PN,
    aa -> Xb4PN^2 V2sym,
    bb -> Xa4PN^2 V2sym,
    cc -> -Xa4PN Xb4PN V2sym,
    dd -> Xb4PN rdSym,
    ee -> -Xa4PN rdSym
  };

  swapFull4PN[expr_] := Expand[
    expr /. {aa -> bb, bb -> aa, cc -> cc, dd -> ee, ee -> dd, pp -> qq, qq -> pp}
  ];

  canonicalSym4PN[expr_] := Module[
    {vars = {pp, qq, aa, bb, cc, dd, ee}, s, rules, coeffs, dens, lcm, ints, gcd, scaled, lead},
    s = Expand[expr + swapFull4PN[expr]];
    rules = CoefficientRules[s, vars];
    coeffs = Last /@ rules;
    dens = Denominator[Together /@ coeffs];
    lcm = Fold[LCM, 1, dens];
    ints = Numerator[Together[lcm #]] & /@ coeffs;
    gcd = Fold[GCD, 0, Abs /@ ints];
    If[gcd === 0, gcd = 1];
    scaled = Expand[lcm s/gcd];
    lead = First @ SortBy[CoefficientRules[scaled, vars], First];
    If[Last[lead] < 0, scaled = -scaled];
    Expand[scaled]
  ];

  generateBasis4PN[massDeg_, velDeg_] := Module[
    {basis = <||>, maxpow = Quotient[velDeg, 2] + 1, expr, sym, tm},
    Do[
      If[2 pa + 2 pb + 2 pc + pd + pe =!= velDeg, Continue[]];
      expr = pp^mp qq^mq aa^pa bb^pb cc^pc dd^pd ee^pe;
      sym = canonicalSym4PN[expr];
      tm = Expand[sym /. {bb -> 0, cc -> 0, ee -> 0, pp -> 0, qq -> 1}];
      If[tm === 0,
        basis[ToString[sym, InputForm]] = sym;
      ],
      {mp, 0, massDeg}, {mq, {massDeg - mp}},
      {pa, 0, maxpow}, {pb, 0, maxpow}, {pc, 0, maxpow},
      {pd, 0, velDeg}, {pe, 0, velDeg}
    ];
    SortBy[Values[basis], ToString[#, InputForm] &]
  ];

  toNu4PN[expr_] := Module[{res},
    res = Expand[expr /. comSubs4PN];
    res = Expand[(res + (res /. delSym -> -delSym))/2];
    res = FixedPoint[Expand[# /. delSym^2 -> 1 - 4 nuSym] &, res];
    Expand[res /. delSym -> 0]
  ];

  toNuEven4PN[expr_] := Module[{res},
    res = Expand[(expr + (expr /. delSym -> -delSym))/2];
    res = FixedPoint[Expand[# /. delSym^2 -> 1 - 4 nuSym] &, res];
    Expand[res /. delSym -> 0]
  ];

  blockSlots4PN[expr_, block_] := Module[{res = Expand[expr]},
    Switch[block,
      "Q",
        {
          Expand[Coefficient[res, V2sym, 4] /. rdSym -> 0],
          Expand[Coefficient[Coefficient[res, V2sym, 3], rdSym, 2]],
          Expand[Coefficient[Coefficient[res, V2sym, 2], rdSym, 4]],
          Expand[Coefficient[Coefficient[res, V2sym, 1], rdSym, 6]],
          Expand[(Coefficient[res, rdSym, 8]) /. V2sym -> 0]
        },
      "T",
        {
          Expand[(Coefficient[res, V2sym, 3]) /. rdSym -> 0],
          Expand[Coefficient[Coefficient[res, V2sym, 2], rdSym, 2]],
          Expand[Coefficient[Coefficient[res, V2sym, 1], rdSym, 4]],
          Expand[(Coefficient[res, rdSym, 6]) /. V2sym -> 0]
        },
      "S",
        {
          Expand[(Coefficient[res, V2sym, 2]) /. rdSym -> 0],
          Expand[Coefficient[Coefficient[res, V2sym, 1], rdSym, 2]],
          Expand[(Coefficient[res, rdSym, 4]) /. V2sym -> 0]
        },
      "U",
        {
          Expand[(Coefficient[res, V2sym, 1]) /. rdSym -> 0],
          Expand[(Coefficient[res, rdSym, 2]) /. V2sym -> 0]
        },
      "W",
        {Expand[res /. {V2sym -> 0, rdSym -> 0}]},
      _,
        Throw[$Failed]
    ]
  ];

  coeffDict4PN[expr_] := Association[CoefficientRules[Expand[expr], {pp, qq, aa, bb, cc, dd, ee}]];

  coordinateMatrix4PN[basis_List] := Module[{mons, pos, mat},
    mons = DeleteDuplicates[Flatten[Keys /@ (coeffDict4PN /@ basis), 1]];
    pos = AssociationThread[mons, Range[Length[mons]]];
    mat = ConstantArray[0, {Length[mons], Length[basis]}];
    Do[
      KeyValueMap[(mat[[pos[#1], j]] = Expand[#2]) &, coeffDict4PN[basis[[j]]]],
      {j, 1, Length[basis]}
    ];
    {mat, pos}
  ];

  coordsInBasis4PN[expr_, basisMat_, monPos_] := Module[{vec, vars, sol},
    vec = ConstantArray[0, Length[monPos]];
    KeyValueMap[(vec[[monPos[#1]]] = Expand[#2]) &, coeffDict4PN[expr]];
    vars = Array[x, Dimensions[basisMat][[2]]];
    sol = Solve[Thread[basisMat.vars == vec], vars];
    If[sol === {}, Throw[$Failed]];
    Expand[vars /. (First[sol] /. C[_] -> 0)]
  ];

  imageMatrixPolynomial4PN[block_, basis_List, maxdeg_] := Module[{rows, slots, poly},
    rows = Table[
      slots = blockSlots4PN[toNu4PN[basis[[j]]], block];
      Flatten@Table[
        poly = Expand[slots[[i]]];
        Table[Expand[Coefficient[poly, nuSym, k]], {k, 1, maxdeg}],
        {i, 1, Length[slots]}
      ],
      {j, 1, Length[basis]}
    ];
    Transpose[rows]
  ];

  targetVector4PN[slots_List, maxdeg_] := Flatten@Table[
    Table[Expand[Coefficient[Expand[slots[[i]]], nuSym, k]], {k, 1, maxdeg}],
    {i, 1, Length[slots]}
  ];

  particularSolution4PN[m_, vec_] := Module[{vars, sol},
    vars = Array[x, Dimensions[m][[2]]];
    sol = Solve[Thread[m.vars == vec], vars];
    If[sol === {}, Throw[$Failed]];
    Expand[vars /. (First[sol] /. C[_] -> 0)]
  ];

  nonzeroTerms4PN[coords_List, basis_List] := Select[
    Table[{i, Expand[coords[[i]]], basis[[i]]}, {i, 1, Min[Length[coords], Length[basis]]}],
    #[[2]] =!= 0 &
  ];

  canonicalExpr4PN[coords_List, basis_List] := Expand[Total[MapThread[Times, {coords, basis}]]];

  unionExprs4PN[lists__List] := SortBy[
    Values @ Merge[
      Association /@ (Map[ToString[#, InputForm] -> # &, #] & /@ {lists}),
      Last
    ],
    ToString[#, InputForm] &
  ];

  evenize4PN[expr_, pr_, pt_, p2_, pr2_] := Expand[
    expr /. {
      pr^(n_Integer?EvenQ) :> pr2^(n/2),
      pt^(n_Integer?EvenQ) :> (p2 - pr2)^(n/2)
    }
  ];

  carriedLowerOrderBlocks4PN[] := Module[
    {nu = nuSym, u, rd, vt, v2, l1, l2, h3, l3c, l3},
    u = Symbol["u"];
    rd = Symbol["rdL"];
    vt = Symbol["vtL"];
    v2 = rd^2 + vt^2;

    l1 =
      ((1 - 3 nu)/8) v2^2 +
      u (((3 + nu)/2) v2 + (nu/2) rd^2) -
      u^2/2;

    l2 =
      ((1 - 5 nu + 5 nu^2)/16) v2^3 +
      u (((7/8) - (7 nu)/4 - nu^2/8) v2^2 + (nu/4 - nu^2/4) rd^2 v2 + (3 nu^2/8) rd^4) +
      u^2 ((2 - (7 nu)/8) v2 + (15 nu/8) rd^2) +
      u^3 (1/4 + (3 nu)/4);

    h3 = <|
      1 -> (-5 + 35 nu - 70 nu^2 + 35 nu^3)/128,
      2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0,
      6 -> (-7 + 42 nu - 53 nu^2 - 5 nu^3)/16,
      7 -> (2 - 3 nu) nu^2/16,
      8 -> 3 (1 - nu) nu^2/16,
      9 -> -5 nu^3/16,
      10 -> (-27 + 136 nu + 109 nu^2)/16,
      11 -> (17 + 30 nu) nu/16,
      12 -> (5 + 43 nu) nu/12,
      13 -> (-600 + (3 Pi^2 - 1340) nu - 552 nu^2)/192,
      14 -> -(340 + 3 Pi^2 + 112 nu) nu/64,
      15 -> (12 + (872 - 63 Pi^2) nu)/96
    |>;

    l3c = <|
      1 -> FullSimplify[(3 nu)/16 - (21 nu^2)/16 + (9 nu^3)/4 - h3[1]],
      2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0,
      6 -> FullSimplify[1/4 + (7 nu)/8 - (35 nu^2)/8 - (21 nu^3)/4 - h3[6]],
      7 -> FullSimplify[(11 nu^2)/8 - (9 nu^3)/2 - h3[7]],
      8 -> FullSimplify[(3 nu^2)/4 - (9 nu^3)/4 - h3[8]],
      9 -> -h3[9],
      10 -> FullSimplify[5/4 + (15 nu)/8 + (123 nu^2)/8 + (13 nu^3)/4 - h3[10]],
      11 -> FullSimplify[(7 nu)/8 + (41 nu^2)/8 + (31 nu^3)/4 - h3[11]],
      12 -> FullSimplify[(9 nu^2)/2 + 4 nu^3 - h3[12]],
      13 -> FullSimplify[-3/2 - (59 nu)/4 - (25 nu^2)/4 - nu^3/2 - h3[13]],
      14 -> FullSimplify[(7 nu)/4 - (31 nu^2)/4 - (7 nu^3)/2 - h3[14]],
      15 -> -h3[15]
    |>;

    l3 =
      l3c[1] v2^4 + l3c[2] v2^3 rd^2 + l3c[3] v2^2 rd^4 + l3c[4] v2 rd^6 + l3c[5] rd^8 +
      u (l3c[6] v2^3 + l3c[7] v2^2 rd^2 + l3c[8] v2 rd^4 + l3c[9] rd^6) +
      u^2 (l3c[10] v2^2 + l3c[11] v2 rd^2 + l3c[12] rd^4) +
      u^3 (l3c[13] v2 + l3c[14] rd^2) +
      u^4 l3c[15];

    <|"l1" -> l1, "l2" -> l2, "l3" -> l3, "u" -> u, "rd" -> rd, "vt" -> vt, "v2" -> v2|>
  ];

  quarticComMap4PN[] := quarticComMap4PN[] = Module[
    {data, l1, l2, l3, u, rd, vt, pr, pt, p2, pr2, coeffs, l4, A0, B0, D0, C0, E0, tcontr,
     h4expr, poly, rules, hmap, indexToMonom},
    data = carriedLowerOrderBlocks4PN[];
    {l1, l2, l3, u, rd, vt} = Lookup[data, {"l1", "l2", "l3", "u", "rd", "vt"}];
    pr = Symbol["pr"];
    pt = Symbol["pt"];
    p2 = Symbol["p2"];
    pr2 = Symbol["pr2"];
    coeffs = Table[Symbol["L" <> ToString[i]], {i, 1, 21}];

    l4 =
      coeffs[[1]] data["v2"]^5 +
      coeffs[[2]] data["v2"]^4 rd^2 +
      coeffs[[3]] data["v2"]^3 rd^4 +
      coeffs[[4]] data["v2"]^2 rd^6 +
      coeffs[[5]] data["v2"] rd^8 +
      coeffs[[6]] rd^10 +
      u (coeffs[[7]] data["v2"]^4 + coeffs[[8]] data["v2"]^3 rd^2 + coeffs[[9]] data["v2"]^2 rd^4 + coeffs[[10]] data["v2"] rd^6 + coeffs[[11]] rd^8) +
      u^2 (coeffs[[12]] data["v2"]^3 + coeffs[[13]] data["v2"]^2 rd^2 + coeffs[[14]] data["v2"] rd^4 + coeffs[[15]] rd^6) +
      u^3 (coeffs[[16]] data["v2"]^2 + coeffs[[17]] data["v2"] rd^2 + coeffs[[18]] rd^4) +
      u^4 (coeffs[[19]] data["v2"] + coeffs[[20]] rd^2) +
      u^5 coeffs[[21]];

    A0 = {D[l1, rd], D[l1, vt]} /. {rd -> pr, vt -> pt};
    B0 = {D[l2, rd], D[l2, vt]} /. {rd -> pr, vt -> pt};
    D0 = {D[l3, rd], D[l3, vt]} /. {rd -> pr, vt -> pt};
    C0 = D[l1, {{rd, vt}, 2}] /. {rd -> pr, vt -> pt};
    E0 = D[l2, {{rd, vt}, 2}] /. {rd -> pr, vt -> pt};

    tcontr = Sum[
      (D[l1, {rd, vt}[[i]], {rd, vt}[[j]], {rd, vt}[[k]]] /. {rd -> pr, vt -> pt}) *
      A0[[i]] * A0[[j]] * A0[[k]],
      {i, 1, 2}, {j, 1, 2}, {k, 1, 2}
    ];

    h4expr = Expand[
      -(l4 /. {rd -> pr, vt -> pt}) +
      A0.D0 +
      (B0.B0)/2 -
      (B0.C0.A0) -
      (A0.E0.A0)/2 +
      (A0.C0.C0.A0)/2 +
      tcontr/6
    ];

    poly = CoefficientRules[evenize4PN[h4expr, pr, pt, p2, pr2], {p2, pr2, u}];
    rules = Association[poly];

    indexToMonom = Association@Join[
      AssociationThread[Range[1, 6], {{5, 0, 0}, {4, 1, 0}, {3, 2, 0}, {2, 3, 0}, {1, 4, 0}, {0, 5, 0}}],
      AssociationThread[Range[7, 11], {{4, 0, 1}, {3, 1, 1}, {2, 2, 1}, {1, 3, 1}, {0, 4, 1}}],
      AssociationThread[Range[12, 15], {{3, 0, 2}, {2, 1, 2}, {1, 2, 2}, {0, 3, 2}}],
      AssociationThread[Range[16, 18], {{2, 0, 3}, {1, 1, 3}, {0, 2, 3}}],
      AssociationThread[Range[19, 20], {{1, 0, 4}, {0, 1, 4}}],
      <|21 -> {0, 0, 5}|>
    ];

    hmap = Association@Table[
      i -> FullSimplify[If[KeyExistsQ[rules, indexToMonom[i]], rules[indexToMonom[i]], 0]],
      {i, 1, 21}
    ];
    <|"map" -> hmap, "coeffs" -> coeffs|>
  ];

  localHamiltonianTarget4PN[] := <|
    1 -> 7/256 - (63 nuSym)/256 + (189 nuSym^2)/256 - (105 nuSym^3)/128 + (63 nuSym^4)/256,
    2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0, 6 -> 0,
    7 -> 45/128 - (45 nuSym)/16 + (423 nuSym^2)/64 - (1013 nuSym^3)/256 - (35 nuSym^4)/128,
    8 -> -(3 nuSym^2)/32 + (23 nuSym^3)/64 - (5 nuSym^4)/32,
    9 -> -(9 nuSym^2)/64 + (69 nuSym^3)/128 - (9 nuSym^4)/64,
    10 -> -(5 nuSym^3)/64 - (5 nuSym^4)/32,
    11 -> (35 nuSym^3)/256 - (35 nuSym^4)/128,
    12 -> 13/8 - (791 nuSym)/64 + (4857 nuSym^2)/256 + (2335 nuSym^3)/256,
    13 -> (49 nuSym)/16 - (545 nuSym^2)/64 + (1135 nuSym^3)/256,
    14 -> -(889 nuSym)/192 + (9475 nuSym^2)/768 - (1649 nuSym^3)/768,
    15 -> (369 nuSym)/160 - (1151 nuSym^2)/128 + (10353 nuSym^3)/1280,
    16 -> 105/32 + (2749 Pi^2 nuSym)/8192 - (589189 nuSym)/19200 + (18491 Pi^2 nuSym^2)/16384 - (1189789 nuSym^2)/28800 - (553 nuSym^3)/128,
    17 -> (63347 nuSym)/1600 - (1059 Pi^2 nuSym)/1024 - (127 nuSym^2)/3 - (4035 Pi^2 nuSym^2)/2048 - (225 nuSym^3)/64,
    18 -> (375 Pi^2 nuSym)/8192 - (23533 nuSym)/1280 + (57563 nuSym^2)/1920 - (38655 Pi^2 nuSym^2)/16384 - (381 nuSym^3)/128,
    19 -> 105/32 + (185761 nuSym)/19200 - (21837 Pi^2 nuSym)/8192 + (672811 nuSym^2)/19200 - (158177 Pi^2 nuSym^2)/49152,
    20 -> (3401779 nuSym)/57600 - (28691 Pi^2 nuSym)/24576 + (110099 Pi^2 nuSym^2)/49152 - (21827 nuSym^2)/3840,
    21 -> -1/16 + (6237 Pi^2 nuSym)/1024 - (169199 nuSym)/2400 + (7403 Pi^2 nuSym^2)/3072 - (1256 nuSym^2)/45
  |>;

  ordinaryTargetFromHamiltonian4PN[] := ordinaryTargetFromHamiltonian4PN[] = Module[
    {mapData, hmap, coeffs, feedback, targetH, targetL},
    mapData = quarticComMap4PN[];
    hmap = mapData["map"];
    coeffs = mapData["coeffs"];
    feedback = Association@Table[i -> FullSimplify[hmap[i] /. Thread[coeffs -> ConstantArray[0, 21]]], {i, 1, 21}];
    targetH = localHamiltonianTarget4PN[];
    targetL = Association@Table[i -> FullSimplify[feedback[i] - targetH[i]], {i, 1, 21}];
    <|"feedback" -> feedback, "targetH" -> targetH, "targetL" -> targetL|>
  ];

  naturalSeedLocalOrdinary4PN[] := Join[
    Association[Table[i -> 0, {i, 1, 21}]],
    <|
      1 -> FullSimplify[(7/256) toNuEven4PN[Xa4PN^9 + Xb4PN^9]],
      7 -> FullSimplify[(75/128) toNuEven4PN[Xa4PN^8 + Xb4PN^8]],
      12 -> FullSimplify[(59/16) toNuEven4PN[Xa4PN^7 + Xb4PN^7]],
      16 -> FullSimplify[(203/32) toNuEven4PN[Xa4PN^6 + Xb4PN^6]],
      19 -> FullSimplify[(31/32) toNuEven4PN[Xa4PN^5 + Xb4PN^5]],
      21 -> FullSimplify[(1/16) toNuEven4PN[Xa4PN^4 + Xb4PN^4]]
    |>
  ];

  hamiltonianSeed4PN[] := <|
    "K" -> toNuEven4PN[(7/256) (Xa4PN^9 + Xb4PN^9)],
    "Q1" -> toNuEven4PN[(45/128) (Xa4PN^8 + Xb4PN^8)],
    "Q2" -> 0, "Q3" -> 0, "Q4" -> 0, "Q5" -> 0,
    "T1" -> toNuEven4PN[(13/8) (Xa4PN^7 + Xb4PN^7)],
    "T2" -> 0, "T3" -> 0, "T4" -> 0,
    "S1" -> toNuEven4PN[(105/32) (Xa4PN^6 + Xb4PN^6)],
    "S2" -> 0, "S3" -> 0,
    "U1" -> toNuEven4PN[(105/32) (Xa4PN^5 + Xb4PN^5)],
    "U2" -> 0,
    "W1" -> toNuEven4PN[-(1/16) (Xa4PN^4 + Xb4PN^4)]
  |>;

  hamiltonianSeedIndexed4PN[] := Module[{named = hamiltonianSeed4PN[]},
    <|
      1 -> named["K"], 2 -> 0, 3 -> 0, 4 -> 0, 5 -> 0, 6 -> 0,
      7 -> named["Q1"], 8 -> named["Q2"], 9 -> named["Q3"], 10 -> named["Q4"], 11 -> named["Q5"],
      12 -> named["T1"], 13 -> named["T2"], 14 -> named["T3"], 15 -> named["T4"],
      16 -> named["S1"], 17 -> named["S2"], 18 -> named["S3"],
      19 -> named["U1"], 20 -> named["U2"],
      21 -> named["W1"]
    |>
  ];

  residualSlots4PN[target_Association, seed_Association] := Module[{res},
    res = AssociationMap[Expand[target[#] - seed[#]] &, Keys[target]];
    <|
      "free" -> {res[1]},
      "Q" -> Table[res[i], {i, 7, 11}],
      "T" -> Table[res[i], {i, 12, 15}],
      "S" -> Table[res[i], {i, 16, 18}],
      "U" -> Table[res[i], {i, 19, 20}],
      "W" -> {res[21]}
    |>
  ];

  canonicalHamiltonianData4PN[] := canonicalHamiltonianData4PN[] = Module[
    {qBasis, tBasis, sBasis, uBasis, wBasis, basisAssoc, targetH, seedIdx, residual, mats, maxdeg,
     coords, blocks, nullspaces, nullPolys},
    qBasis = generateBasis4PN[0, 8];
    tBasis = generateBasis4PN[1, 6];
    sBasis = generateBasis4PN[2, 4];
    uBasis = generateBasis4PN[3, 2];
    wBasis = generateBasis4PN[4, 0];
    basisAssoc = <|"Q" -> qBasis, "T" -> tBasis, "S" -> sBasis, "U" -> uBasis, "W" -> wBasis|>;

    targetH = localHamiltonianTarget4PN[];
    seedIdx = hamiltonianSeedIndexed4PN[];
    residual = residualSlots4PN[targetH, seedIdx];
    maxdeg = <|"Q" -> 4, "T" -> 3, "S" -> 3, "U" -> 2, "W" -> 2|>;
    mats = AssociationMap[imageMatrixPolynomial4PN[#, basisAssoc[#], maxdeg[#]] &, Keys[basisAssoc]];
    coords = AssociationMap[
      particularSolution4PN[mats[#], targetVector4PN[residual[#], maxdeg[#]]] &,
      Keys[basisAssoc]
    ];
    blocks = AssociationMap[canonicalExpr4PN[coords[#], basisAssoc[#]] &, Keys[basisAssoc]];
    nullspaces = AssociationMap[NullSpace[mats[#]] &, Keys[basisAssoc]];
    nullPolys = AssociationMap[
      Expand[canonicalExpr4PN[#, basisAssoc[#2]]] & @@@
        Table[{nullspaces[key][[i]], key}, {key, Keys[nullspaces]}, {i, 1, Length[nullspaces[key]]}] // Flatten[#, 1] &,
      {}
    ];
    <|
      "basis" -> basisAssoc,
      "targetH" -> targetH,
      "seedHIndexed" -> seedIdx,
      "residualSlots" -> residual,
      "matrices" -> mats,
      "coords" -> coords,
      "blocks" -> blocks,
      "nullspaces" -> nullspaces
    |>
  ];

  blockTargetSlots4PN[targetL_Association, block_] := targetL /@ blockRanges4PN[block];

  ordinaryTranslationData4PN[] := ordinaryTranslationData4PN[] = Module[
    {ordData, targetL, feedback, naturalSeed, hData, hBlocks, residualH, lBlocks, seedHIdx,
     alignedSeed, ordResidualAligned, misalign, misBlocks, oldBases, seedBases, coords, deltaExprs,
     qNat, tNat, sNat, uNat, wNat, alignedBlocks, fullBlocks},
    ordData = ordinaryTargetFromHamiltonian4PN[];
    targetL = ordData["targetL"];
    feedback = ordData["feedback"];
    naturalSeed = naturalSeedLocalOrdinary4PN[];
    hData = canonicalHamiltonianData4PN[];
    hBlocks = hData["blocks"];
    residualH = hData["residualSlots"];
    lBlocks = AssociationMap[Expand[-hBlocks[#]] &, Keys[hBlocks]];
    seedHIdx = hamiltonianSeedIndexed4PN[];
    alignedSeed = Association@Table[i -> FullSimplify[feedback[i] - seedHIdx[i]], {i, 1, 21}];
    ordResidualAligned = <|
      "Q" -> Table[FullSimplify[targetL[i] - alignedSeed[i]], {i, 7, 11}],
      "T" -> Table[FullSimplify[targetL[i] - alignedSeed[i]], {i, 12, 15}],
      "S" -> Table[FullSimplify[targetL[i] - alignedSeed[i]], {i, 16, 18}],
      "U" -> Table[FullSimplify[targetL[i] - alignedSeed[i]], {i, 19, 20}],
      "W" -> {FullSimplify[targetL[21] - alignedSeed[21]]}
    |>;
    misalign = Association@Table[i -> FullSimplify[alignedSeed[i] - naturalSeed[i]], {i, 1, 21}];
    misBlocks = AssociationMap[Table[misalign[i], {i, blockRanges4PN[#]}] &, Keys[blockRanges4PN]];

    oldBases = <|
      "Q" -> generateBasis4PN[0, 8],
      "T" -> generateBasis4PN[1, 6],
      "S" -> generateBasis4PN[2, 4],
      "U" -> generateBasis4PN[3, 2],
      "W" -> generateBasis4PN[4, 0]
    |>;
    seedBases = <|
      "Q" -> oldBases["Q"],
      "T" -> unionExprs4PN[oldBases["T"], Expand[pp qq #] & /@ oldBases["T"]],
      "S" -> unionExprs4PN[oldBases["S"], Expand[pp qq #] & /@ oldBases["S"]],
      "U" -> unionExprs4PN[oldBases["U"], Expand[(pp qq)^2 #] & /@ oldBases["U"]],
      "W" -> {}
    |>;

    coords = <|
      "Q" -> particularSolution4PN[imageMatrixPolynomial4PN["Q", seedBases["Q"], 4], targetVector4PN[misBlocks["Q"], 4]],
      "T" -> particularSolution4PN[imageMatrixPolynomial4PN["T", seedBases["T"], 4], targetVector4PN[misBlocks["T"], 4]],
      "S" -> particularSolution4PN[imageMatrixPolynomial4PN["S", seedBases["S"], 4], targetVector4PN[misBlocks["S"], 4]],
      "U" -> particularSolution4PN[imageMatrixPolynomial4PN["U", seedBases["U"], 4], targetVector4PN[misBlocks["U"], 4]]
    |>;

    deltaExprs = <|
      "Q" -> canonicalExpr4PN[coords["Q"], seedBases["Q"]],
      "T" -> canonicalExpr4PN[coords["T"], seedBases["T"]],
      "S" -> canonicalExpr4PN[coords["S"], seedBases["S"]],
      "U" -> canonicalExpr4PN[coords["U"], seedBases["U"]],
      "W" -> 0
    |>;

    qNat = (75/128) (aa^4 + bb^4);
    tNat = (59/16) (qq aa^3 + pp bb^3);
    sNat = (203/32) (qq^2 aa^2 + pp^2 bb^2);
    uNat = (31/32) (qq^3 aa + pp^3 bb);
    wNat = (1/16) (qq^4 + pp^4);

    alignedBlocks = <|
      "Q" -> Expand[qNat + deltaExprs["Q"]],
      "T" -> Expand[tNat + deltaExprs["T"]],
      "S" -> Expand[sNat + deltaExprs["S"]],
      "U" -> Expand[uNat + deltaExprs["U"]],
      "W" -> Expand[wNat]
    |>;

    fullBlocks = <|
      "Q" -> Expand[alignedBlocks["Q"] + lBlocks["Q"]],
      "T" -> Expand[alignedBlocks["T"] + lBlocks["T"]],
      "S" -> Expand[alignedBlocks["S"] + lBlocks["S"]],
      "U" -> Expand[alignedBlocks["U"] + lBlocks["U"]],
      "W" -> Expand[alignedBlocks["W"] + lBlocks["W"]]
    |>;

    <|
      "targetL" -> targetL,
      "feedback" -> feedback,
      "naturalSeed" -> naturalSeed,
      "alignedSeed" -> alignedSeed,
      "hamiltonianBlocks" -> hBlocks,
      "ordinaryResidualBlocks" -> lBlocks,
      "ordinaryResidualSlotsAligned" -> ordResidualAligned,
      "misalignment" -> misalign,
      "misBlocks" -> misBlocks,
      "seedBases" -> seedBases,
      "seedCoords" -> coords,
      "deltaExprs" -> deltaExprs,
      "alignedBlocks" -> alignedBlocks,
      "fullBlocks" -> fullBlocks
    |>
  ];
];
