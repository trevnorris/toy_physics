(* PathA-43 density quadrupole port, Mathematica engine.

   Standalone print-only gate:

     timeout 600 math -script software/stage1_solver/tools/pathA_43_density_quadrupole_port.wl

   This engine derives the density two-port numerator by Green-function/DtN
   interface matching and a different elimination order from the SymPy engine.
   It does not read SymPy output.
*)

ClearAll["Global`*"];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
assertTrue[cond_, msg_] := If[! TrueQ[cond], fail[msg]];
boolText[x_] := If[TrueQ[x], "PASS", "FAIL"];
exprText[x_] := ToString[FullSimplify[x], InputForm];
listText[x_] := StringRiffle[ToString[#, InputForm] & /@ x, ", "];

a = Symbol["a"]; cs = Symbol["cS"]; cc = Symbol["c"]; gravG = Symbol["G"];
d0 = Symbol["D0"]; rho = Symbol["rho"]; i25 = Symbol["I25"];
iw2 = Symbol["Iwrong2"]; xiQ = Symbol["XiQ"]; xiDeferred = Symbol["Xideferred"];
etaq = Symbol["etaq"]; etaphi = Symbol["etaphi"];
varpiq = Symbol["varpiq2"]; varpiphi = Symbol["varpiPhi2"];
lambdac = Symbol["lambdac"]; q2sym = Symbol["q2"]; phi2sym = Symbol["Phi2"];
sigmaHidden = Symbol["sigmahidden"]; freeCarrier = Symbol["freecarrier"];
omegaU = Symbol["OmegaU"]; omegaW = Symbol["OmegaW"]; rmix = Symbol["Rmix"];
gU = Symbol["gU"]; gW = Symbol["gW"];
aw = Symbol["Aw"]; fmuw = Symbol["Fmuw"]; jwvec = Symbol["Jwvec"];
uvec = Symbol["Uvec"]; wvec = Symbol["Wvec"];
omegaWall = Symbol["omegawall"]; omegaRho = Symbol["omegarho"];
rMixRelabel = Symbol["rmixrel"]; gRho = Symbol["grho"]; gQOld = Symbol["gqold"];
z = Symbol["z"]; omega = Symbol["omega"];

pretty[s_String] := StringReplace[s, {
  "cS" -> "c_s", "Iwrong2" -> "I_wrong2", "XiQ" -> "Xi_Q", "Xideferred" -> "Xi_deferred",
  "etaq" -> "eta_q", "etaphi" -> "eta_phi", "varpiq2" -> "varpi_q2",
  "varpiPhi2" -> "varpi_Phi2", "lambdac" -> "lambda_c", "sigmahidden" -> "sigma_hidden",
  "freecarrier" -> "free_carrier", "OmegaU" -> "Omega_U", "OmegaW" -> "Omega_W",
  "Rmix" -> "R_mix", "gU" -> "g_U", "gW" -> "g_W", "omegawall" -> "omega_wall",
  "Aw" -> "A_w", "Fmuw" -> "F_muw", "Jwvec" -> "Jw", "Uvec" -> "U", "Wvec" -> "W",
  "omegarho" -> "omega_rho", "rmixrel" -> "r_mix", "grho" -> "g_rho", "gqold" -> "g_qold"
}];
exprText[x_] := pretty[ToString[FullSimplify[x], InputForm]];

$Assumptions = And @@ Thread[
    {a, cs, cc, gravG, d0, rho, i25, iw2, xiQ, xiDeferred, etaq, etaphi,
     varpiq, varpiphi, lambdac} > 0
  ];

zeroDim = {0, 0, 0};
dimAdd[x_, y_] := x + y;
dimSub[x_, y_] := x - y;
dimScale[x_, q_] := q x;

dimToText[dim_] := Module[{pairs, parts, expText},
  expText[e_] := If[Denominator[Together[e]] === 1, ToString[Numerator[Together[e]]],
    ToString[Numerator[Together[e]]] <> "/" <> ToString[Denominator[Together[e]]]];
  pairs = {{"L", dim[[1]]}, {"T", dim[[3]]}, {"M", dim[[2]]}};
  parts = (If[#2 === 1, #1, #1 <> "^" <> expText[#2]] &) @@@ Select[pairs, ! TrueQ[#[[2]] == 0] &];
  If[Length[parts] == 0, "1", StringRiffle[parts, " "]]
];

dimOf[expr_, dims_Association] := Module[{head, args, termDims, first, base, pow},
  Which[
    TrueQ[expr === 0] || NumberQ[expr], zeroDim,
    AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr],
    AtomQ[expr], fail["missing sourced dimension for " <> ToString[expr, InputForm]],
    Head[expr] === Times, Total[dimOf[#, dims] & /@ (List @@ expr)],
    Head[expr] === Power,
      base = expr[[1]]; pow = expr[[2]];
      If[! NumberQ[pow], fail["non-numeric dimension exponent in " <> ToString[expr, InputForm]]];
      dimScale[dimOf[base, dims], pow],
    Head[expr] === Plus,
      args = Select[List @@ expr, ! TrueQ[# == 0] &];
      termDims = dimOf[#, dims] & /@ args;
      If[Length[termDims] == 0, zeroDim,
        first = First[termDims];
        If[! And @@ (TrueQ[# == first] & /@ termDims),
          fail["dimension mismatch in sum " <> ToString[expr, InputForm]]
        ];
        first
      ],
    True, fail["unsupported dimension expression " <> ToString[expr, InputForm]]
  ]
];

scalePower[expr_, powers_Association] := Module[{args, vals, first, base, pow},
  Which[
    TrueQ[expr === 0] || NumberQ[expr], 0,
    TrueQ[expr === a], 1,
    AtomQ[expr] && KeyExistsQ[powers, expr], powers[expr],
    AtomQ[expr], 0,
    Head[expr] === Times,
      vals = scalePower[#, powers] & /@ (List @@ expr);
      If[MemberQ[vals, Missing["unknown"]], Missing["unknown"], Total[vals]],
    Head[expr] === Power,
      base = expr[[1]]; pow = expr[[2]];
      If[! NumberQ[pow], Missing["unknown"],
        With[{p = scalePower[base, powers]},
          If[p === Missing["unknown"], Missing["unknown"], pow p]
        ]
      ],
    Head[expr] === Plus,
      args = Select[List @@ expr, ! TrueQ[# == 0] &];
      vals = scalePower[#, powers] & /@ args;
      If[Length[vals] == 0, 0,
        If[MemberQ[vals, Missing["unknown"]], Missing["unknown"],
          first = First[vals];
          If[And @@ (TrueQ[# == first] & /@ vals), first, Missing["unknown"]]
        ]
      ],
    True, Missing["unknown"]
  ]
];

baseDims = <|
  a -> {1, 0, 0}, cs -> {1, 0, -1}, cc -> {1, 0, -1},
  gravG -> {3, -1, -2}, d0 -> {-1, 1, -2}, rho -> {-3, 1, 0},
  i25 -> {5/2, 0, 0}, iw2 -> {2, 0, 0}, xiQ -> zeroDim,
  xiDeferred -> zeroDim, etaq -> zeroDim, etaphi -> zeroDim,
  varpiq -> {0, 0, -2}, varpiphi -> {0, 0, -2}, lambdac -> {0, 0, -2},
  sigmaHidden -> zeroDim, freeCarrier -> zeroDim,
  aw -> zeroDim, fmuw -> zeroDim, jwvec -> zeroDim, uvec -> zeroDim, wvec -> zeroDim,
  omegaU -> {0, 0, -1}, omegaW -> {0, 0, -1}, rmix -> {0, 0, -2},
  gU -> {-1/2, 1/2, -2}, gW -> {-1/2, 1/2, -2},
  omegaWall -> {0, 0, -1}, omegaRho -> {0, 0, -1}, rMixRelabel -> {0, 0, -2},
  gRho -> {-1/2, 1/2, -2}, gQOld -> {-1/2, 1/2, -2}
|>;

basePowers = <|
  cs -> 0, cc -> 0, gravG -> 0, d0 -> 0, rho -> 0, i25 -> 0, iw2 -> 0,
  xiQ -> 0, xiDeferred -> 0, etaq -> 0, etaphi -> 0,
  varpiq -> -2, varpiphi -> -2, lambdac -> -2,
  sigmaHidden -> 0, freeCarrier -> 0,
  aw -> 0, fmuw -> 0, jwvec -> 0, uvec -> 0, wvec -> 0,
  omegaU -> -1, omegaW -> -1, rmix -> -2, gU -> -1, gW -> -1,
  omegaWall -> -1, omegaRho -> -1, rMixRelabel -> -2, gRho -> -1, gQOld -> -1
|>;

cfg[name_, rules___Rule] := Join[
  <|
    "name" -> name, "kind" -> "density", "couplingZero" -> False,
    "corruptDimension" -> False, "freeCarrierRider" -> False,
    "freeCarrierTagged" -> False,
    "incomingSign" -> False, "couplingAPower" -> -7/2,
    "deferredUncertified" -> False, "provenDeferred" -> False,
    "fakeContinuity" -> False, "hiddenVectorIntermediate" -> False,
    "corruptContinuityMoment" -> False, "vectorInjectedDensity" -> False
  |>,
  <|rules|>
];

vectorSymbols = {aw, fmuw, jwvec, uvec, wvec, omegaU, omegaW, rmix, gU, gW};
continuityOperatorId = "pathA_29_projected_continuity";
continuityL0 = "M0 = Integral(S_leak d3x)";
continuityL1 = "D1_i = Integral(x_i*S_leak d3x) + Integral(j_i d3x)";
continuityL2 = "Q2_m = Integral(Y2_m_star*S_leak d3x)";
continuityL2Kernel = "Y2_m_star*S_leak";
badContinuityL2 = "GARBAGE_NOT_A_CONTINUITY_MOMENT_AT_ALL";

baseSourceTags = <|
  a -> {"pathA_29_bulk", "pathA_32_wall"},
  cs -> {"pathA_29_bulk"}, cc -> {"calibrated_anchor"}, gravG -> {"calibrated_anchor"},
  d0 -> {"calibrated_anchor"}, rho -> {"pathA_29_bulk"},
  xiQ -> {"continuity_interface"}, xiDeferred -> {"continuity_interface"},
  etaq -> {"continuity_interface"}, etaphi -> {"continuity_interface"},
  varpiq -> {"pathA_32_wall"}, varpiphi -> {"pathA_29_bulk"},
  lambdac -> {"continuity_interface", "pathA_29_bulk", "pathA_32_wall"},
  sigmaHidden -> {"vector_port"}, freeCarrier -> {},
  aw -> {"vector_port"}, fmuw -> {"vector_port"}, jwvec -> {"vector_port"},
  uvec -> {"vector_port"}, wvec -> {"vector_port"},
  omegaU -> {"vector_port"}, omegaW -> {"vector_port"}, rmix -> {"vector_port"},
  gU -> {"vector_port"}, gW -> {"vector_port"},
  omegaWall -> {"vector_port"}, omegaRho -> {"vector_port"},
  rMixRelabel -> {"vector_port"}, gRho -> {"vector_port"}, gQOld -> {"vector_port"}
|>;

containsAll[text_, tokens_] := And @@ (StringContainsQ[ToString[text], #] & /@ tokens);

lineageValidQ[lineage_Association] := Module[
  {moments = Lookup[lineage, "moments", <||>], l0, l1, l2, l2Kernel},
  l0 = Lookup[moments, "l0", ""];
  l1 = Lookup[moments, "l1", ""];
  l2 = Lookup[moments, "l2", ""];
  l2Kernel = Lookup[lineage, "l2Kernel", ""];
  TrueQ[Lookup[lineage, "operatorId", ""] === continuityOperatorId] &&
  containsAll[l0, {"M0", "Integral", "S_leak", "d3x"}] &&
  containsAll[l1, {"D1_i", "Integral", "x_i", "S_leak", "j_i", "d3x"}] &&
  containsAll[l2, {"Q2_m", "Integral", "Y2_m_star", "S_leak", "d3x"}] &&
  containsAll[l2Kernel, {"Y2_m_star", "S_leak"}]
];

continuityMomentSymbol[c_Association, lineage_Association] := Module[{valid = lineageValidQ[lineage]},
  {If[TrueQ[valid] && TrueQ[c["couplingAPower"] == -7/2], i25, iw2], valid}
];

sourceTagMap[momentSym_, momentValid_, freeCarrierTagged_:False] := Module[{out = baseSourceTags},
  If[TrueQ[freeCarrierTagged],
    out = Join[out, <|freeCarrier -> {"pathA_34_dimensionless_free_carrier"}|>]
  ];
  If[TrueQ[momentValid],
    Join[out, <|momentSym -> {"continuity_interface", "pathA_32_wall"}|>],
    Join[out, <|momentSym -> {}|>]
  ]
];

taintForExpr[expr_, tagMap_Association] := Module[{symbols, tags = {}, missing = {}},
  symbols = Variables[{expr}];
  Scan[
    If[KeyExistsQ[tagMap, #],
      tags = Join[tags, tagMap[#]],
      missing = Append[missing, #]
    ] &,
    symbols
  ];
  {Sort[DeleteDuplicates[tags]], missing}
];

sourceGraphForExpr[expr_, tagMap_Association, coordinateHosts_] := Module[
  {symbols, tagsMissing, tags, missing, empty, vectorHosts, symbolTags},
  symbols = SortBy[Variables[{expr}], pretty[ToString[#, InputForm]] &];
  tagsMissing = taintForExpr[expr, tagMap]; tags = tagsMissing[[1]]; missing = tagsMissing[[2]];
  empty = Select[symbols, KeyExistsQ[tagMap, #] && Length[Lookup[tagMap, #, {}]] == 0 &];
  vectorHosts = Intersection[symbols, vectorSymbols];
  symbolTags = Association@Table[
    pretty[ToString[sym, InputForm]] -> Sort[Lookup[tagMap, sym, {}]],
    {sym, symbols}
  ];
  <|"symbolTags" -> symbolTags, "taintSet" -> tags,
    "missingSourceSymbols" -> (pretty[ToString[#, InputForm]] & /@ missing),
    "emptySourceSymbols" -> (pretty[ToString[#, InputForm]] & /@ empty),
    "vectorHostSymbols" -> (pretty[ToString[#, InputForm]] & /@ vectorHosts),
    "coordinateHosts" -> Sort[coordinateHosts]|>
];

vectorAblatedExpr[expr_, tagMap_Association] := Module[
  {vectorTaintedSymbols, ablateSymbols},
  vectorTaintedSymbols = Keys[Select[tagMap, MemberQ[#, "vector_port"] &]];
  ablateSymbols = Intersection[Variables[{expr}], DeleteDuplicates[Join[vectorSymbols, vectorTaintedSymbols]]];
  Quiet[FullSimplify[expr /. Thread[ablateSymbols -> 0], $Assumptions]]
];

lineageFor[c_Association] := Which[
  TrueQ[c["fakeContinuity"]],
    <|"operatorId" -> "mis_tagged_vector_formula",
      "moments" -> <|"l2" -> "Omega_U^2*g_W + R_mix*g_U, relabeled as continuity"|>,
      "failure" -> "no l=0->1->2 projected-continuity lineage; no Integral(Y2star*S_leak*d3x)"|>,
  c["kind"] === "density",
    <|"operatorId" -> continuityOperatorId,
      "moments" -> <|
        "l0" -> continuityL0,
        "l1" -> continuityL1,
        "l2" -> If[TrueQ[c["corruptContinuityMoment"]], badContinuityL2, continuityL2]|>,
      "l2Kernel" -> continuityL2Kernel,
      "lineage" -> "same projected-continuity operator extended to the l=2 Y2* source moment"|>,
  True,
    <|"operatorId" -> "none", "moments" -> <||>,
      "failure" -> "vector port has no continuity moment lineage"|>
];

densityExpression[c_Association, lineage_Association] := Module[
  {xi, momentPair, q2Moment, momentValid, gBase, gqDen, gphiDen, qRule, phiEq, phiResponse, pden, delta, n0, trace},
  xi = If[TrueQ[c["deferredUncertified"]] || TrueQ[c["provenDeferred"]], xiDeferred, xiQ];
  momentPair = continuityMomentSymbol[c, lineage];
  q2Moment = momentPair[[1]]; momentValid = momentPair[[2]];
  gBase = q2Moment xi cs^2 Sqrt[rho]/a^(7/2);
  If[! TrueQ[c["couplingAPower"] == -7/2],
    assertTrue[TrueQ[c["couplingAPower"] == -3], "only isolated scaling control power is implemented"];
    gBase = q2Moment xi cs^2 Sqrt[rho]/a^3
  ];
  If[TrueQ[c["vectorInjectedDensity"]], gBase = gBase omegaU/omegaW];
  If[TrueQ[c["couplingZero"]], gBase = 0];
  gqDen = FullSimplify[gBase etaq, $Assumptions];
  gphiDen = FullSimplify[gBase etaphi, $Assumptions];

  (* DtN/interface matching route: solve wall equation for q2 first, then
     insert it into the bulk Phi2 interface equation. *)
  qRule = q2sym -> FullSimplify[(lambdac phi2sym + gqDen)/varpiq, $Assumptions];
  phiEq = FullSimplify[varpiphi phi2sym - lambdac q2sym == gphiDen /. qRule, $Assumptions];
  phiResponse = FullSimplify[phi2sym /. First[Solve[phiEq, phi2sym]], $Assumptions];
  delta = FullSimplify[varpiq varpiphi - lambdac^2, $Assumptions];
  pden = FullSimplify[varpiq gphiDen + lambdac gqDen, $Assumptions];
  n0 = FullSimplify[phiResponse^2, $Assumptions];
  If[TrueQ[c["hiddenVectorIntermediate"]], n0 = FullSimplify[n0 sigmaHidden, $Assumptions]];
  trace = <|
    "method" -> "Mathematica Green-function/DtN interface matching",
    "qEliminated" -> qRule,
    "phiEquationAfterWallElimination" -> phiEq,
    "Phi2Response" -> phiResponse,
    "DeltaDen" -> delta,
    "PDen" -> pden,
    "gBase" -> gBase,
    "gqDen" -> gqDen,
    "gphiDen" -> gphiDen,
    "continuityMomentSymbol" -> q2Moment,
    "continuityMomentValid" -> momentValid
  |>;
  {FullSimplify[n0, $Assumptions], trace}
];

vectorExpression[c_Association] := Module[{p, delta},
  If[c["kind"] === "relabel",
    p = omegaWall^2 gRho + rMixRelabel gQOld;
    delta = omegaWall^2 omegaRho^2 - rMixRelabel^2,
    p = omegaU^2 gW + rmix gU;
    delta = omegaU^2 omegaW^2 - rmix^2
  ];
  {FullSimplify[p^2/delta^2, $Assumptions], <|"method" -> "old vector port fixture", "Pold" -> p, "DeltaOld" -> delta|>}
];

derive[c_Association] := Module[
  {lineage, momentPair, momentSym, momentValid, tagMap, exprTrace, expr, trace,
   hosts, coordinateHosts, sourceGraph},
  lineage = lineageFor[c];
  momentPair = continuityMomentSymbol[c, lineage];
  momentSym = momentPair[[1]]; momentValid = momentPair[[2]];
  tagMap = sourceTagMap[momentSym, momentValid, c["freeCarrierTagged"]];
  Which[
    MemberQ[{"vector", "relabel"}, c["kind"]],
      (exprTrace = vectorExpression[c]; expr = exprTrace[[1]]; trace = exprTrace[[2]]),
    TrueQ[c["fakeContinuity"]],
      (exprTrace = vectorExpression[cfg[c["name"], "kind" -> "vector"]];
       expr = exprTrace[[1]]; trace = exprTrace[[2]]),
    True,
      (exprTrace = densityExpression[c, lineage]; expr = exprTrace[[1]]; trace = exprTrace[[2]])
  ];
  If[TrueQ[c["freeCarrierRider"]], expr = FullSimplify[expr freeCarrier, $Assumptions]];
  hosts = Sort[pretty /@ (ToString[#, InputForm] & /@ Variables[{expr}])];
  coordinateHosts = If[c["kind"] === "density", {"q2", "Phi2"}, {}];
  If[TrueQ[c["fakeContinuity"]], coordinateHosts = {"q2_fake", "Phi2_fake"}];
  hosts = Sort[DeleteDuplicates[Join[hosts, coordinateHosts]]];
  sourceGraph = sourceGraphForExpr[expr, tagMap, coordinateHosts];
  <|"expr" -> FullSimplify[expr, $Assumptions], "trace" -> trace,
    "tags" -> sourceGraph["taintSet"], "sourceTagMap" -> tagMap, "sourceGraph" -> sourceGraph,
    "lineage" -> lineage, "lineageValid" -> lineageValidQ[lineage],
    "continuityMomentSymbol" -> momentSym, "continuityMomentValid" -> momentValid,
    "hostSymbols" -> hosts|>
];

dtnSign[kind_] := Module[{j2, y2, h, lam, yhat, ser, coeff},
  j2 = (3/z^3 - 1/z) Sin[z] - 3 Cos[z]/z^2;
  y2 = (1/z - 3/z^3) Cos[z] - 3 Sin[z]/z^2;
  h = If[kind === "outgoing", j2 + I y2, j2 - I y2];
  lam = FullSimplify[z D[h, z]/h];
  yhat = FullSimplify[-3/lam];
  ser = Collect[Normal[Series[yhat, {z, 0, 6}]], z, FullSimplify];
  coeff = FullSimplify[Coefficient[ser, z, 5]/I];
  <|"kind" -> kind, "YhatSeries" -> ser, "radiativeCoeff" -> coeff,
    "matchesOutgoing" -> TrueQ[FullSimplify[coeff == 1/27]]|>
];

closureOverlay[n0_] := Module[{p0phys, kbar0, kbar2, kbar4, kbar4res, gamma5, gammares},
  p0phys = FullSimplify[(cs/a)^2 n0/d0, $Assumptions];
  kbar0 = 54 gravG cs^5/(5 a^5 cc^5);
  kbar2 = 6 gravG cs^3/(5 a^3 cc^5);
  kbar4 = 8 gravG cs/(15 a cc^5);
  kbar4res = FullSimplify[kbar4 - 4 kbar2^2/kbar0, $Assumptions];
  gamma5 = FullSimplify[kbar0 a^5/(27 cs^5), $Assumptions];
  gammares = FullSimplify[gamma5 - 2 gravG/(5 cc^5), $Assumptions];
  <|"P0Physical" -> p0phys, "Kbar0" -> kbar0, "Kbar2" -> kbar2, "Kbar4" -> kbar4,
    "Kbar4Residual" -> kbar4res, "Gamma5Residual" -> gammares,
    "ok" -> TrueQ[kbar4res == 0] && TrueQ[gammares == 0]|>
];

evaluate[c_Association] := Module[
  {base, expr, tags, sourceGraph, tagMap, ablatedExpr, dims, powers, dimVal, dimText,
   dimOk, n0Power, p0Power, scalingOk, scalingWrong, scalingUndecidable, sign,
   signOk, vectorHostSymbols, sourceMapComplete, continuityDependencyOk,
   vanishedContinuityCoupling, originOk, vectorOk, nonzeroOk, closure, verdict, bad},
  base = derive[c];
  expr = base["expr"]; tags = base["tags"]; sourceGraph = base["sourceGraph"]; tagMap = base["sourceTagMap"];
  ablatedExpr = vectorAblatedExpr[expr, tagMap];
  dims = baseDims;
  If[TrueQ[c["corruptDimension"]], dims = Join[dims, <|i25 -> baseDims[i25] + {1, 0, 0}|>]];
  powers = basePowers;
  If[TrueQ[c["deferredUncertified"]], powers = Join[powers, <|xiDeferred -> Missing["unknown"]|>]];
  If[TrueQ[c["provenDeferred"]], powers = Join[powers, <|xiDeferred -> 0|>]];
  dimVal = Quiet[Check[dimOf[expr, dims], "inhomogeneous"]];
  dimText = If[ListQ[dimVal], dimToText[dimVal], "inhomogeneous"];
  dimOk = TrueQ[dimVal == {-1, 1, 0}];
  n0Power = scalePower[expr, powers];
  p0Power = If[n0Power === Missing["unknown"], Missing["unknown"], FullSimplify[n0Power - 2 - powers[d0]]];
  scalingOk = TrueQ[p0Power == -5];
  scalingWrong = p0Power =!= Missing["unknown"] && ! TrueQ[p0Power == -5];
  scalingUndecidable = p0Power === Missing["unknown"];
  sign = dtnSign[If[TrueQ[c["incomingSign"]], "incoming", "outgoing"]];
  signOk = TrueQ[sign["matchesOutgoing"]];
  vectorHostSymbols = sourceGraph["vectorHostSymbols"];
  sourceMapComplete = TrueQ[
    Length[sourceGraph["missingSourceSymbols"]] == 0 &&
    Length[sourceGraph["emptySourceSymbols"]] == 0
  ];
  continuityDependencyOk = TrueQ[
    base["lineageValid"] && base["continuityMomentValid"] &&
    (MemberQ[Variables[{expr}], base["continuityMomentSymbol"]] ||
      (TrueQ[FullSimplify[expr, $Assumptions] == 0] && TrueQ[c["couplingZero"]]))
  ];
  vanishedContinuityCoupling = TrueQ[
    FullSimplify[expr, $Assumptions] == 0 && TrueQ[c["couplingZero"]] &&
    Length[vectorHostSymbols] == 0 && sourceMapComplete && continuityDependencyOk
  ];
  originOk = TrueQ[
    (MemberQ[tags, "continuity_interface"] && ! MemberQ[tags, "vector_port"] &&
      Length[vectorHostSymbols] == 0 && sourceMapComplete && continuityDependencyOk) ||
    vanishedContinuityCoupling
  ];
  vectorOk = TrueQ[
    Length[vectorHostSymbols] == 0 && ! MemberQ[tags, "vector_port"] &&
    FullSimplify[expr - ablatedExpr, $Assumptions] == 0
  ];
  nonzeroOk = ! TrueQ[FullSimplify[expr, $Assumptions] == 0];
  closure = closureOverlay[expr];
  verdict = Which[
    ! originOk || ! vectorOk, "FAIL_NOT_DENSITY_DERIVED",
    ! nonzeroOk, "FAIL_PORT_VANISHES",
    ! dimOk || scalingWrong || ! signOk,
      bad = {};
      If[! dimOk, bad = Append[bad, "dimensional"]];
      If[scalingWrong, bad = Append[bad, "scaling"]];
      If[! signOk, bad = Append[bad, "sign"]];
      "FAIL_PORT_MALFORMED(" <> StringRiffle[bad, ","] <> ")",
    originOk && nonzeroOk && dimOk && scalingOk && signOk && vectorOk && TrueQ[closure["ok"]],
      "DENSITY_PORT_HOSTED",
    True, "PORT_INCONCLUSIVE_SIM_DEFERRED"
  ];
  <|
    "expr" -> expr, "trace" -> base["trace"], "hostSymbols" -> base["hostSymbols"], "taintSet" -> tags,
    "sourceGraph" -> sourceGraph, "lineage" -> base["lineage"], "lineageValid" -> base["lineageValid"],
    "continuityMomentSymbol" -> pretty[ToString[base["continuityMomentSymbol"], InputForm]],
    "continuityDependencyOk" -> continuityDependencyOk, "vectorHostSymbols" -> vectorHostSymbols,
    "sourceMapComplete" -> sourceMapComplete, "emptySourceSymbols" -> sourceGraph["emptySourceSymbols"],
    "ablationExpr" -> ablatedExpr,
    "ablationDelta" -> FullSimplify[expr - ablatedExpr, $Assumptions],
    "checks" -> <|"origin_derivation_ok" -> originOk, "nonzero" -> nonzeroOk, "dimension" -> dimOk,
      "a_scaling_provenance_ok" -> scalingOk, "radiative_sign" -> signOk,
      "vector_independence_ok" -> vectorOk|>,
    "dimensions" -> <|"N0_den" -> dimText|>,
    "scaling" -> <|"N0_den_a_power" -> n0Power, "P0_physical_a_power" -> p0Power|>,
    "sign" -> sign, "closure" -> closure, "verdict" -> verdict
  |>
];

controls[] := Module[{fixtures, out},
  fixtures = <|
    "vector_hosted" -> cfg["vector_hosted", "kind" -> "vector"],
    "relabel_rig" -> cfg["relabel_rig", "kind" -> "relabel"],
    "fake_continuity" -> cfg["fake_continuity", "fakeContinuity" -> True],
    "ablation_isolation" -> cfg["ablation_isolation", "hiddenVectorIntermediate" -> True],
    "attack2_continuity_corruption" -> cfg["attack2_continuity_corruption", "corruptContinuityMoment" -> True],
    "attack5_vector_injection" -> cfg["attack5_vector_injection", "vectorInjectedDensity" -> True],
    "zero_coupling" -> cfg["zero_coupling", "couplingZero" -> True],
    "dimensional" -> cfg["dimensional", "corruptDimension" -> True],
    "sign" -> cfg["sign", "incomingSign" -> True],
    "scaling" -> cfg["scaling", "couplingAPower" -> -3],
    "deferred_scalar" -> cfg["deferred_scalar", "deferredUncertified" -> True],
    "deferred_scalar_proven_converse" -> cfg["deferred_scalar_proven_converse", "provenDeferred" -> True],
    "free_carrier_dimension_corruption" -> cfg["free_carrier_dimension_corruption",
      "freeCarrierRider" -> True, "freeCarrierTagged" -> True],
    "provenance_less_rider" -> cfg["provenance_less_rider", "freeCarrierRider" -> True]
  |>;
  Association @ KeyValueMap[(#1 -> evaluate[#2]) &, fixtures]
];

assertGate[baseline_, cps_] := Module[{expected},
  assertTrue[baseline["verdict"] === "DENSITY_PORT_HOSTED", "baseline verdict " <> baseline["verdict"]];
  Scan[assertTrue[#, "baseline check failed"] &, Values[baseline["checks"]]];
  assertTrue[baseline["taintSet"] === {"continuity_interface", "pathA_29_bulk", "pathA_32_wall"}, "bad taint set"];
  assertTrue[baseline["sourceGraph"]["emptySourceSymbols"] === {}, "baseline empty provenance symbol"];
  Scan[assertTrue[Length[#] > 0, "baseline empty symbol tag set"] &, Values[baseline["sourceGraph"]["symbolTags"]]];
  assertTrue[baseline["dimensions"]["N0_den"] === "L^-1 M", "bad N0 dimension"];
  assertTrue[TrueQ[baseline["scaling"]["P0_physical_a_power"] == -5], "bad P0 a-scaling"];
  assertTrue[TrueQ[baseline["closure"]["Kbar4Residual"] == 0], "bad Kbar4 residual"];
  expected = <|
    "vector_hosted" -> "FAIL_NOT_DENSITY_DERIVED",
    "relabel_rig" -> "FAIL_NOT_DENSITY_DERIVED",
    "fake_continuity" -> "FAIL_NOT_DENSITY_DERIVED",
    "ablation_isolation" -> "FAIL_NOT_DENSITY_DERIVED",
    "attack2_continuity_corruption" -> "FAIL_NOT_DENSITY_DERIVED",
    "attack5_vector_injection" -> "FAIL_NOT_DENSITY_DERIVED",
    "zero_coupling" -> "FAIL_PORT_VANISHES",
    "dimensional" -> "FAIL_PORT_MALFORMED(dimensional)",
    "sign" -> "FAIL_PORT_MALFORMED(sign)",
    "scaling" -> "FAIL_PORT_MALFORMED(scaling)",
    "deferred_scalar" -> "PORT_INCONCLUSIVE_SIM_DEFERRED",
    "deferred_scalar_proven_converse" -> "DENSITY_PORT_HOSTED",
    "free_carrier_dimension_corruption" -> "DENSITY_PORT_HOSTED",
    "provenance_less_rider" -> "FAIL_NOT_DENSITY_DERIVED"|>;
  KeyValueMap[
    (assertTrue[cps[#1]["verdict"] === #2,
      "control " <> #1 <> ": expected " <> #2 <> ", got " <> cps[#1]["verdict"]]) &,
    expected
  ];
];

printResult[baseline_, cps_] := Module[{},
  Print["PATHA_43_DENSITY_QUADRUPOLE_PORT_MATHEMATICA"];
  Print["method: ", baseline["trace"]["method"]];
  Print["verdict: ", baseline["verdict"]];
  Print["N0_den: ", exprText[baseline["expr"]]];
  Print["port_picture: ii two-port(q2,Phi2)"];
  Print["host_set: ", StringRiffle[baseline["hostSymbols"], ", "]];
  Print["taint_set: ", StringRiffle[baseline["taintSet"], ", "]];
  Print["vector_host_symbols: ", If[Length[baseline["vectorHostSymbols"]] == 0, "empty", StringRiffle[baseline["vectorHostSymbols"], ", "]]];
  Print["continuity_moment_symbol: ", baseline["continuityMomentSymbol"]];
  Print["continuity_dependency_ok: ", boolText[baseline["continuityDependencyOk"]]];
  Print["source_graph_symbol_tags:"];
  KeyValueMap[(Print["  ", #1, ": ", If[Length[#2] == 0, "empty", StringRiffle[#2, ","]]]) &, baseline["sourceGraph"]["symbolTags"]];
  Print["lineage:"];
  KeyValueMap[(Print["  ", #1, ": ", #2]) &, baseline["lineage"]["moments"]];
  Print["lineage_valid: ", boolText[baseline["lineageValid"]]];
  Print["derivation_trace:"];
  Print["  g_base: ", exprText[baseline["trace"]["gBase"]]];
  Print["  g_q_den: ", exprText[baseline["trace"]["gqDen"]]];
  Print["  g_phi_den: ", exprText[baseline["trace"]["gphiDen"]]];
  Print["  P_den: ", exprText[baseline["trace"]["PDen"]]];
  Print["  Delta_den: ", exprText[baseline["trace"]["DeltaDen"]]];
  Print["  Phi2_response: ", exprText[baseline["trace"]["Phi2Response"]]];
  Print["checks:"];
  KeyValueMap[(Print["  ", #1, ": ", boolText[#2]]) &, baseline["checks"]];
  Print["vector_ablation_delta: ", exprText[baseline["ablationDelta"]]];
  Print["scaling:"];
  Print["  N0_den_a_power: ", baseline["scaling"]["N0_den_a_power"]];
  Print["  P0_physical_a_power: ", baseline["scaling"]["P0_physical_a_power"]];
  Print["closure:"];
  Print["  P0_physical: ", exprText[baseline["closure"]["P0Physical"]]];
  Print["  Kbar4_residual: ", exprText[baseline["closure"]["Kbar4Residual"]]];
  Print["  Gamma5_residual_to_2G_over_5c5: ", exprText[baseline["closure"]["Gamma5Residual"]]];
  Print["scope_tags:"];
  Print["  CALIBRATED: G, 2/5, 54/5"];
  Print["  SIM_DEFERRED: Xi_Q exact branch magnitude, eta_q, eta_phi, lambda_c literal throat value"];
  Print["controls:"];
  KeyValueMap[
    (Print["  ", #1, ": verdict=", #2["verdict"], "; ablation_delta=", exprText[#2["ablationDelta"]],
      "; taint_set=", If[Length[#2["taintSet"]] == 0, "empty", StringRiffle[#2["taintSet"], ","]],
      "; source_map_complete=", boolText[#2["sourceMapComplete"]],
      "; empty_source_symbols=", If[Length[#2["emptySourceSymbols"]] == 0, "empty", StringRiffle[#2["emptySourceSymbols"], ","]]]) &,
    cps
  ];
];

baseline = evaluate[cfg["baseline"]];
cps = controls[];
assertGate[baseline, cps];
printResult[baseline, cps];
