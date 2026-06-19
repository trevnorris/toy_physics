#!/usr/bin/env wolframscript

(* MT15-03 / Spike-1: ell=2 vector-spherical-harmonic localized Maxwell
   operator, gauge deflation, MMS, and assembly bookkeeping checks.

   Scope is deliberately limited to Decision-05 Spike-1:
   - no outgoing DtN/open-impedance solve,
   - no wall->matter->Maxwell transfer,
   - no N0/Z0 extraction,
   - no V2 packet emission.
*)

$HistoryLength = 0;

ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/mt15_03_spike1_vsh_maxwell_operator.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
repoRoot = AbsoluteFileName[FileNameJoin[{scriptDir, "..", "..", ".."}]];

If[MemberQ[$ScriptCommandLine, "--patha-b2b"],
  Get[FileNameJoin[{scriptDir, "mt15_03_patha_b2b_vsh_maxwell_operator.wls"}]];
  Exit[0]
];

runDir = FileNameJoin[{scriptDir, "runs", "mt15_03_spike1_vsh_maxwell_operator"}];
If[! DirectoryQ[runDir], CreateDirectory[runDir, CreateIntermediateDirectories -> True]];

diagnosticsPath = FileNameJoin[{runDir, "mt15_03_spike1_vsh_maxwell_operator_diagnostics.json"}];
reportPath = FileNameJoin[{scriptDir, "mt15_03_spike1_vsh_maxwell_operator_report.md"}];

ClearAll[
  smootherStep, zWeight, r0Geometry, toJsonNumber, toJsonVector,
  diagSparse, firstDerivativeMatrix, secondDerivativeMatrix, gridData, blockRow, weightedGram,
  buildSystem, sourceCommitString, sampleExpr, weightedRelativeLaneError,
  maskedWeightedRelativeLaneError, interiorMask,
  estimateOrders, vectorNormWeighted, maxAbsFieldStrength, fieldNormFromComponents,
  makeGaugeVector, deflatedVector, formatScientific, formatFixed, markdownTable,
  buildMmsSourceExpressions, buildStrongMatrixFromEuler
];

ell = 2;
Lharm = ell (ell + 1);
xiGauge = 1.0;

(* Target-blind geometry/localization placeholders, matching the M1a/M1b
   smooth open-throat geometry style.  R0 is recorded for reuse/provenance;
   the closed MMS harness uses a rectangular interior r-grid away from r=0 so
   the singular VSH coefficients are tested without origin-parity machinery. *)
wLength = 2.0;
aMouth = 1.0;
rExit = 1.65;
rMinClosed = 0.35;
zFloor = 0.08;
zAmp = 0.92;
zCenter = 0.5 wLength;
zSigma = 0.38 wLength;

smootherStep[x_] := 3 x^2 - 2 x^3;
r0Geometry[w_] := aMouth + (rExit - aMouth) smootherStep[w/wLength];
zWeight[w_] := zFloor + zAmp Exp[-((w - zCenter)/zSigma)^2];

toJsonNumber[x_] := Module[{n = Quiet[Check[N[x], $Failed]]},
  If[n === $Failed || ! NumericQ[n],
    Print["ERROR: nonnumeric JSON export value: ", InputForm[x]];
    Exit[4]
  ];
  If[Abs[Im[n]] < 10^-13, N[Re[n], 17], N[n, 17]]
];
toJsonVector[x_List] := toJsonNumber /@ x;

sourceCommitString[] := Module[
  {gitDir, headPath, head, refName, refPath, packedPath, packedLines, match, sha},
  gitDir = FileNameJoin[{repoRoot, ".git"}];
  headPath = FileNameJoin[{gitDir, "HEAD"}];
  If[! FileExistsQ[headPath], Return["unavailable"]];

  head = Quiet[Check[StringTrim[Import[headPath, "Text"]], ""]];
  If[! StringQ[head] || StringLength[head] == 0, Return["unavailable"]];

  sha = If[StringStartsQ[head, "ref:"],
    refName = StringTrim[StringDrop[head, 4]];
    refPath = FileNameJoin[Join[{gitDir}, FileNameSplit[refName]]];
    If[FileExistsQ[refPath],
      Quiet[Check[StringTrim[Import[refPath, "Text"]], ""]],
      packedPath = FileNameJoin[{gitDir, "packed-refs"}];
      If[! FileExistsQ[packedPath],
        "",
        packedLines = StringSplit[Quiet[Check[Import[packedPath, "Text"], ""]], "\n"];
        match = Select[packedLines, StringMatchQ[#, Except["#"] ~~ ___ ~~ " " ~~ refName] &];
        If[Length[match] > 0, First[StringSplit[First[match]]], ""]
      ]
    ],
    head
  ];

  If[StringQ[sha] && StringLength[sha] >= 7, StringTake[sha, 7], "unavailable"]
];

diagSparse[v_List] := SparseArray[Band[{1, 1}] -> N[v, 17], {Length[v], Length[v]}];

firstDerivativeMatrix[n_Integer, h_?NumericQ] := Module[{rules},
  rules = Join[
    Table[If[i < n, {i, i + 1} -> 1.0/(2.0 h), Nothing], {i, 1, n}],
    Table[If[i > 1, {i, i - 1} -> -1.0/(2.0 h), Nothing], {i, 1, n}]
  ];
  SparseArray[rules, {n, n}]
];

secondDerivativeMatrix[n_Integer, h_?NumericQ] := Module[{rules},
  rules = Join[
    Table[{i, i} -> -2.0/h^2, {i, 1, n}],
    Table[If[i < n, {i, i + 1} -> 1.0/h^2, Nothing], {i, 1, n}],
    Table[If[i > 1, {i, i - 1} -> 1.0/h^2, Nothing], {i, 1, n}]
  ];
  SparseArray[rules, {n, n}]
];

gridData[nr_Integer, nw_Integer] := Module[
  {dr, dw, rVals1, wVals1, rVals, wVals, zVals, r0Vals, quadWeights},
  dr = (rExit - rMinClosed)/(nr + 1);
  dw = wLength/(nw + 1);
  rVals1 = N[Table[rMinClosed + i dr, {i, 1, nr}], 17];
  wVals1 = N[Table[i dw, {i, 1, nw}], 17];
  rVals = N[Flatten[Table[rVals1[[i]], {j, 1, nw}, {i, 1, nr}]], 17];
  wVals = N[Flatten[Table[wVals1[[j]], {j, 1, nw}, {i, 1, nr}]], 17];
  zVals = N[zWeight /@ wVals, 17];
  r0Vals = N[r0Geometry /@ wVals, 17];
  quadWeights = N[dr dw rVals^2 zVals, 17];
  <|
    "nr" -> nr, "nw" -> nw, "n" -> nr nw,
    "dr" -> dr, "dw" -> dw,
    "r1" -> rVals1, "w1" -> wVals1,
    "r" -> rVals, "w" -> wVals,
    "z" -> zVals, "R0" -> r0Vals,
    "w0" -> quadWeights,
    "h" -> Max[dr, dw]
  |>
];

blockRow[blocks_List] := ArrayFlatten[{blocks}];
weightedGram[m_, wdiag_, sign_] := sign ConjugateTranspose[m].wdiag.m;

buildSystem[nr_Integer, nw_Integer, omega_?NumericQ] := Module[
  {
    grid, n, total, drMat1, dwMat1, idR, idW, idN, zN, drMat, dwMat,
    rVals, invR, invR2, rDiag, r2Diag, invRDiag, invR2Diag, w0Diag, lw0Diag,
    wLaneVals, wLane, wLaneInv, divR, mG, mEr, mEE, mEB, mEw, mCr, mCE,
    mCB, mBB, mBr, mBE, kAction, kFieldPositive, aAction, aSmoke,
    crossRules, crossBlock, gaugeMap, gaugeAnchorRow, offsets
  },
  grid = gridData[nr, nw];
  n = grid["n"];
  total = 5 n;
  drMat1 = firstDerivativeMatrix[nr, grid["dr"]];
  dwMat1 = firstDerivativeMatrix[nw, grid["dw"]];
  idR = SparseArray[Band[{1, 1}] -> 1.0, {nr, nr}];
  idW = SparseArray[Band[{1, 1}] -> 1.0, {nw, nw}];
  idN = SparseArray[Band[{1, 1}] -> 1.0, {n, n}];
  zN = SparseArray[{}, {n, n}];
  drMat = KroneckerProduct[idW, drMat1];
  dwMat = KroneckerProduct[dwMat1, idR];

  rVals = grid["r"];
  invR = N[1.0/rVals, 17];
  invR2 = N[invR^2, 17];
  rDiag = diagSparse[rVals];
  r2Diag = diagSparse[rVals^2];
  invRDiag = diagSparse[invR];
  invR2Diag = diagSparse[invR2];
  w0Diag = diagSparse[grid["w0"]];
  lw0Diag = Lharm w0Diag;
  wLaneVals = Join[grid["w0"], grid["w0"], Lharm grid["w0"], Lharm grid["w0"], grid["w0"]];
  wLane = diagSparse[wLaneVals];
  wLaneInv = diagSparse[1.0/wLaneVals];

  divR = invR2Diag . drMat . r2Diag;

  mG = blockRow[{I omega idN, divR, -Lharm invRDiag, zN, dwMat}];
  mEr = blockRow[{-drMat, I omega idN, zN, zN, zN}];
  mEE = blockRow[{-invRDiag, zN, I omega idN, zN, zN}];
  mEB = blockRow[{zN, zN, zN, I omega idN, zN}];
  mEw = blockRow[{-dwMat, zN, zN, zN, I omega idN}];
  mCr = blockRow[{zN, -dwMat, zN, zN, drMat}];
  mCE = blockRow[{zN, zN, -dwMat, zN, invRDiag}];
  mCB = blockRow[{zN, zN, zN, -dwMat, zN}];
  mBB = blockRow[{zN, -invRDiag, invRDiag . drMat . rDiag, zN, zN}];
  mBr = blockRow[{zN, zN, zN, -Lharm invRDiag, zN}];
  mBE = blockRow[{zN, zN, zN, -invRDiag . drMat . rDiag, zN}];

  kAction =
    weightedGram[mEr, w0Diag, 1.0] +
    weightedGram[mEE, lw0Diag, 1.0] +
    weightedGram[mEB, lw0Diag, 1.0] +
    weightedGram[mEw, w0Diag, 1.0] -
    weightedGram[mCr, w0Diag, 1.0] -
    weightedGram[mCE, lw0Diag, 1.0] -
    weightedGram[mCB, lw0Diag, 1.0] -
    weightedGram[mBB, lw0Diag, 1.0] -
    weightedGram[mBr, w0Diag, 1.0] -
    weightedGram[mBE, lw0Diag, 1.0] -
    (1.0/xiGauge) weightedGram[mG, w0Diag, 1.0];

  kFieldPositive =
    weightedGram[mEr, w0Diag, 1.0] +
    weightedGram[mEE, lw0Diag, 1.0] +
    weightedGram[mEB, lw0Diag, 1.0] +
    weightedGram[mEw, w0Diag, 1.0] +
    weightedGram[mCr, w0Diag, 1.0] +
    weightedGram[mCE, lw0Diag, 1.0] +
    weightedGram[mCB, lw0Diag, 1.0] +
    weightedGram[mBB, lw0Diag, 1.0] +
    weightedGram[mBr, w0Diag, 1.0] +
    weightedGram[mBE, lw0Diag, 1.0];

  aAction = wLaneInv . kAction;

  (* The scalarized smoke wrapper is represented by deleting exactly the D1
     vector-Laplacian cross block from the ar/aE strong operator. *)
  crossRules = Join[
    Table[{n + k, n + k} -> -2.0 invR2[[k]], {k, 1, n}],
    Table[{n + k, 2 n + k} -> 2.0 Lharm invR2[[k]], {k, 1, n}],
    Table[{2 n + k, n + k} -> 2.0 invR2[[k]], {k, 1, n}]
  ];
  crossBlock = SparseArray[crossRules, {total, total}];
  aSmoke = aAction - crossBlock;

  gaugeMap = ArrayFlatten[{
      {I omega idN},
      {drMat},
      {invRDiag},
      {zN},
      {dwMat}
    }];

  gaugeAnchorRow = SparseArray[{{1, Ceiling[n/2]} -> 1.0}, {1, n}];
  offsets = <|"phi" -> 0, "a_r" -> n, "a_E" -> 2 n, "a_B" -> 3 n, "a_w" -> 4 n|>;

  Join[
    grid,
    <|
      "omega" -> omega,
      "total" -> total,
      "Dr" -> drMat,
      "Dw" -> dwMat,
      "W0" -> w0Diag,
      "WLane" -> wLane,
      "WLaneInv" -> wLaneInv,
      "WLaneValues" -> wLaneVals,
      "KAction" -> kAction,
      "KFieldPositive" -> kFieldPositive,
      "AAction" -> aAction,
      "ASmoke" -> aSmoke,
      "CrossBlock" -> crossBlock,
      "GaugeMap" -> gaugeMap,
      "GaugeAnchorRow" -> gaugeAnchorRow,
      "FieldMatrices" -> <|
        "G" -> mG,
        "E_r" -> mEr, "E_E" -> mEE, "E_B" -> mEB, "E_w" -> mEw,
        "C_r" -> mCr, "C_E" -> mCE, "C_B" -> mCB,
        "B_B" -> mBB, "B_r" -> mBr, "B_E" -> mBE
      |>,
      "Offsets" -> offsets
    |>
  ]
];

buildMmsSourceExpressions[] := Module[
  {
    r, w, ph, ar, aE, aB, aw, zExpr, g, er, ee, eb, ew, cr, ce, cb,
    bb, br, be, density, euler, x, y, s1, s2, t1, t2, t3,
    phiExact, arExact, aEExact, aBExact, awExact, depRules, allRules,
    strongExprs, sourceExprs
  },
  Clear[r, w, ph, ar, aE, aB, aw];
  zExpr = zFloor + zAmp Exp[-((w - zCenter)/zSigma)^2];
  g = (1/r^2) D[r^2 ar[r, w], r] - Lharm aE[r, w]/r + D[aw[r, w], w];
  er = -D[ph[r, w], r];
  ee = -ph[r, w]/r;
  eb = 0 ph[r, w];
  ew = -D[ph[r, w], w];
  cr = D[aw[r, w], r] - D[ar[r, w], w];
  ce = aw[r, w]/r - D[aE[r, w], w];
  cb = -D[aB[r, w], w];
  bb = D[r aE[r, w], r]/r - ar[r, w]/r;
  br = -Lharm aB[r, w]/r;
  be = -D[r aB[r, w], r]/r;

  density = 0.5 r^2 zExpr (
      er^2 + Lharm ee^2 + Lharm eb^2 + ew^2 -
      cr^2 - Lharm ce^2 - Lharm cb^2 -
      Lharm bb^2 - br^2 - Lharm be^2 -
      (1/xiGauge) g^2
    );

  euler[f_] :=
    D[density, f[r, w]] -
    D[D[density, Derivative[1, 0][f][r, w]], r] -
    D[D[density, Derivative[0, 1][f][r, w]], w];

  x = (r - rMinClosed)/(rExit - rMinClosed);
  y = w/wLength;
  s1 = Sin[Pi x];
  s2 = Sin[2 Pi x];
  t1 = Sin[Pi y];
  t2 = Sin[2 Pi y];
  t3 = Sin[3 Pi y];

  phiExact = 0.17 s1 t1 + 0.05 s1^2 t2;
  arExact = 0.31 s1^2 t1 + 0.04 s2 t3;
  aEExact = -0.27 s1 t1^2 + 0.06 s2 t1;
  aBExact = 0.22 s2 t1 + 0.03 s1 t3;
  awExact = -0.19 s1^3 t1 + 0.05 s1 t2;

  depRules[f_, expr_] := Join[
    {f[r, w] -> expr},
    Flatten[Table[Derivative[i, j][f][r, w] -> D[expr, {r, i}, {w, j}], {i, 0, 3}, {j, 0, 3}]]
  ];
  allRules = Join[
    depRules[ph, phiExact],
    depRules[ar, arExact],
    depRules[aE, aEExact],
    depRules[aB, aBExact],
    depRules[aw, awExact]
  ];
  strongExprs = Simplify[
    {
      euler[ph]/(r^2 zExpr),
      euler[ar]/(r^2 zExpr),
      euler[aE]/(Lharm r^2 zExpr),
      euler[aB]/(Lharm r^2 zExpr),
      euler[aw]/(r^2 zExpr)
    }
  ];
  sourceExprs = Simplify[strongExprs /. allRules];
  <|
    "variables" -> {r, w},
    "heads" -> {ph, ar, aE, aB, aw},
    "exact" -> {phiExact, arExact, aEExact, aBExact, awExact},
    "strong" -> strongExprs,
    "source" -> sourceExprs
  |>
];

buildStrongMatrixFromEuler[sys_Association, strongExprs_List, heads_List, rSym_, wSym_] := Module[
  {
    n = sys["n"], nr = sys["nr"], nw = sys["nw"], idR, idW, idN,
    dr1, dw1, drr1, dww1, drMat, dwMat, drrMat, dwwMat, derivMatrices,
    derivativeKeys, zBlock, rowBlocks, coeffExpr, coeffVals
  },
  idR = SparseArray[Band[{1, 1}] -> 1.0, {nr, nr}];
  idW = SparseArray[Band[{1, 1}] -> 1.0, {nw, nw}];
  idN = SparseArray[Band[{1, 1}] -> 1.0, {n, n}];
  dr1 = firstDerivativeMatrix[nr, sys["dr"]];
  dw1 = firstDerivativeMatrix[nw, sys["dw"]];
  drr1 = secondDerivativeMatrix[nr, sys["dr"]];
  dww1 = secondDerivativeMatrix[nw, sys["dw"]];
  drMat = KroneckerProduct[idW, dr1];
  dwMat = KroneckerProduct[dw1, idR];
  drrMat = KroneckerProduct[idW, drr1];
  dwwMat = KroneckerProduct[dww1, idR];
  derivMatrices = <|
    {0, 0} -> idN,
    {1, 0} -> drMat,
    {0, 1} -> dwMat,
    {2, 0} -> drrMat,
    {0, 2} -> dwwMat,
    {1, 1} -> drMat . dwMat
  |>;
  derivativeKeys = Keys[derivMatrices];
  rowBlocks = Table[
    Table[
      Total[
        Table[
          coeffExpr = If[key === {0, 0},
            Coefficient[strongExprs[[outLane]], heads[[inLane]][rSym, wSym]],
            Coefficient[
              strongExprs[[outLane]],
              Derivative[key[[1]], key[[2]]][heads[[inLane]]][rSym, wSym]
            ]
          ];
          If[TrueQ[coeffExpr === 0],
            SparseArray[{}, {n, n}],
            coeffVals = sampleExpr[coeffExpr, rSym, wSym, sys["r"], sys["w"]];
            diagSparse[coeffVals] . derivMatrices[key]
          ],
          {key, derivativeKeys}
        ]
      ],
      {inLane, 1, Length[heads]}
    ],
    {outLane, 1, Length[strongExprs]}
  ];
  ArrayFlatten[rowBlocks]
];

sampleExpr[expr_, rSym_, wSym_, rVals_List, wVals_List] :=
  N[MapThread[(expr /. {rSym -> #1, wSym -> #2}) &, {rVals, wVals}], 17];

vectorNormWeighted[vec_List, weights_List] := Sqrt[Max[0.0, Re[Conjugate[vec].(weights vec)]]];

weightedRelativeLaneError[err_List, ref_List, weights_List] := Module[{den},
  den = vectorNormWeighted[ref, weights];
  vectorNormWeighted[err, weights]/Max[den, 10^-30]
];

maskedWeightedRelativeLaneError[err_List, ref_List, weights_List, mask_List] :=
  weightedRelativeLaneError[Pick[err, mask], Pick[ref, mask], Pick[weights, mask]];

interiorMask[sys_Association, layer_Integer] := Flatten[
  Table[
    i > layer && i <= sys["nr"] - layer && j > layer && j <= sys["nw"] - layer,
    {j, 1, sys["nw"]}, {i, 1, sys["nr"]}
  ]
];

estimateOrders[hs_List, errs_List] := Table[
  N[Log[errs[[i]]/errs[[i + 1]]]/Log[hs[[i]]/hs[[i + 1]]], 10],
  {i, 1, Length[errs] - 1}
];

fieldNormFromComponents[sys_Association, vec_List] := Module[
  {fm = sys["FieldMatrices"], w0, lw0, terms, compVals},
  w0 = sys["W0"];
  lw0 = Lharm w0;
  terms = {
    {fm["E_r"], w0}, {fm["E_E"], lw0}, {fm["E_B"], lw0}, {fm["E_w"], w0},
    {fm["C_r"], w0}, {fm["C_E"], lw0}, {fm["C_B"], lw0},
    {fm["B_B"], lw0}, {fm["B_r"], w0}, {fm["B_E"], lw0}
  };
  compVals = Table[
    With[{component = terms[[i, 1]] . vec, weight = terms[[i, 2]]},
      Conjugate[component] . (weight . component)
    ],
    {i, 1, Length[terms]}
  ];
  Sqrt[Max[0.0, Re[Total[compVals]]]]
];

maxAbsFieldStrength[sys_Association, vec_List] := Module[{fm = sys["FieldMatrices"], vals},
  vals = (Abs[# . vec] & /@ Values[KeyDrop[fm, "G"]]);
  Max[Flatten[vals]]
];

makeGaugeVector[sys_Association, chi_List] := Module[
  {omega = sys["omega"], drMat = sys["Dr"], dwMat = sys["Dw"], rVals = sys["r"], n = sys["n"]},
  Join[
    I omega chi,
    drMat . chi,
    chi/rVals,
    ConstantArray[0.0, n],
    dwMat . chi
  ]
];

deflatedVector[sys_Association, vec_List] := Module[
  {b = sys["GaugeMap"], w = sys["WLane"], gram, rhs, coeff},
  gram = ConjugateTranspose[b].w.b;
  rhs = ConjugateTranspose[b].w.vec;
  coeff = LinearSolve[gram, rhs];
  vec - b.coeff
];

formatScientific[x_] := StringReplace[ToString[FortranForm[N[x, 8]]], " " -> ""];
formatFixed[x_] := ToString[NumberForm[N[x, 8], {10, 6}]];

markdownTable[headers_List, rows_List] := StringRiffle[
  Join[
    {"| " <> StringRiffle[headers, " | "] <> " |"},
    {"| " <> StringRiffle[ConstantArray["---", Length[headers]], " | "] <> " |"},
    ("| " <> StringRiffle[ToString /@ #, " | "] <> " |") & /@ rows
  ],
  "\n"
];

Print["MT15-03 Spike-1 VSH Maxwell operator"];
Print["building MMS source expressions from the D1 weak density"];
mmsExpr = buildMmsSourceExpressions[];
{rSym, wSym} = mmsExpr["variables"];
fieldHeads = mmsExpr["heads"];
exactExprs = mmsExpr["exact"];
strongExprs = mmsExpr["strong"];
sourceExprs = mmsExpr["source"];

mmsGrids = {{8, 8}, {12, 12}, {18, 18}, {26, 26}};
laneOrder = {"phi", "a_r", "a_E", "a_B", "a_w"};
cleanSmokeLanes = {"phi", "a_B", "a_w"};
weakMmsClosureLayer = 2;
weakMmsRows = {};
strongMmsRows = {};
consistencyRows = {};
mmsErrorsWeak = AssociationThread[laneOrder, ConstantArray[{}, Length[laneOrder]]];
mmsErrorsWeakSmoke = AssociationThread[laneOrder, ConstantArray[{}, Length[laneOrder]]];
mmsErrorsStrong = AssociationThread[laneOrder, ConstantArray[{}, Length[laneOrder]]];
mmsErrorsStrongSmoke = AssociationThread[laneOrder, ConstantArray[{}, Length[laneOrder]]];
weakStrongRatios = AssociationThread[laneOrder, ConstantArray[{}, Length[laneOrder]]];
hs = {};

Do[
  {nr, nw} = gridPair;
  sys = buildSystem[nr, nw, 0.0];
  strongOp = buildStrongMatrixFromEuler[sys, strongExprs, fieldHeads, rSym, wSym];
  strongSmokeOp = strongOp - sys["CrossBlock"];
  weakOp = sys["AAction"];
  weakSmokeOp = sys["ASmoke"];
  n = sys["n"];
  weightsByLane = {
    sys["w0"], sys["w0"], Lharm sys["w0"], Lharm sys["w0"], sys["w0"]
  };
  weakMask = interiorMask[sys, weakMmsClosureLayer];
  uExact = Join @@ (sampleExpr[#, rSym, wSym, sys["r"], sys["w"]] & /@ exactExprs);
  fExact = Join @@ (sampleExpr[#, rSym, wSym, sys["r"], sys["w"]] & /@ sourceExprs);
  weakAction = weakOp.uExact;
  weakSmokeAction = weakSmokeOp.uExact;
  strongAction = strongOp.uExact;
  strongSmokeAction = strongSmokeOp.uExact;
  residualWeak = weakAction - fExact;
  residualWeakSmoke = weakSmokeAction - fExact;
  residualStrong = strongAction - fExact;
  residualStrongSmoke = strongSmokeAction - fExact;
  weakStrongDelta = weakAction - strongAction;
  AppendTo[hs, sys["h"]];

  laneErrWeak = Table[
    maskedWeightedRelativeLaneError[
      residualWeak[[(lane - 1) n + 1 ;; lane n]],
      fExact[[(lane - 1) n + 1 ;; lane n]],
      weightsByLane[[lane]],
      weakMask
    ],
    {lane, 1, 5}
  ];
  laneErrWeakSmoke = Table[
    maskedWeightedRelativeLaneError[
      residualWeakSmoke[[(lane - 1) n + 1 ;; lane n]],
      fExact[[(lane - 1) n + 1 ;; lane n]],
      weightsByLane[[lane]],
      weakMask
    ],
    {lane, 1, 5}
  ];
  laneErrStrong = Table[
    weightedRelativeLaneError[
      residualStrong[[(lane - 1) n + 1 ;; lane n]],
      fExact[[(lane - 1) n + 1 ;; lane n]],
      weightsByLane[[lane]]
    ],
    {lane, 1, 5}
  ];
  laneErrStrongSmoke = Table[
    weightedRelativeLaneError[
      residualStrongSmoke[[(lane - 1) n + 1 ;; lane n]],
      fExact[[(lane - 1) n + 1 ;; lane n]],
      weightsByLane[[lane]]
    ],
    {lane, 1, 5}
  ];
  laneWeakStrongRatio = Table[
    maskedWeightedRelativeLaneError[
      weakStrongDelta[[(lane - 1) n + 1 ;; lane n]],
      strongAction[[(lane - 1) n + 1 ;; lane n]],
      weightsByLane[[lane]],
      weakMask
    ],
    {lane, 1, 5}
  ];

  Do[
    mmsErrorsWeak[laneOrder[[lane]]] = Append[mmsErrorsWeak[laneOrder[[lane]]], laneErrWeak[[lane]]];
    mmsErrorsWeakSmoke[laneOrder[[lane]]] = Append[mmsErrorsWeakSmoke[laneOrder[[lane]]], laneErrWeakSmoke[[lane]]];
    mmsErrorsStrong[laneOrder[[lane]]] = Append[mmsErrorsStrong[laneOrder[[lane]]], laneErrStrong[[lane]]];
    mmsErrorsStrongSmoke[laneOrder[[lane]]] = Append[mmsErrorsStrongSmoke[laneOrder[[lane]]], laneErrStrongSmoke[[lane]]];
    weakStrongRatios[laneOrder[[lane]]] = Append[weakStrongRatios[laneOrder[[lane]]], laneWeakStrongRatio[[lane]]];
    ,
    {lane, 1, 5}
  ];
  AppendTo[weakMmsRows, Join[{nr, nw, formatFixed[sys["h"]]}, formatScientific /@ laneErrWeak, formatScientific /@ laneErrWeakSmoke]];
  AppendTo[strongMmsRows, Join[{nr, nw, formatFixed[sys["h"]]}, formatScientific /@ laneErrStrong, formatScientific /@ laneErrStrongSmoke]];
  AppendTo[consistencyRows, Join[{nr, nw, formatFixed[sys["h"]]}, formatScientific /@ laneWeakStrongRatio]];
  Print["MMS grid ", nr, "x", nw,
    " AAction max rel err: ", formatScientific[Max[laneErrWeak]],
    " ASmoke max rel err: ", formatScientific[Max[laneErrWeakSmoke]],
    " strongOp max rel err: ", formatScientific[Max[laneErrStrong]],
    " weak-vs-strong max ratio: ", formatScientific[Max[laneWeakStrongRatio]]];
  ,
  {gridPair, mmsGrids}
];

mmsOrdersWeak = AssociationMap[estimateOrders[hs, mmsErrorsWeak[#]] &, laneOrder];
mmsOrdersWeakSmoke = AssociationMap[estimateOrders[hs, mmsErrorsWeakSmoke[#]] &, laneOrder];
mmsOrdersStrong = AssociationMap[estimateOrders[hs, mmsErrorsStrong[#]] &, laneOrder];
mmsOrdersStrongSmoke = AssociationMap[estimateOrders[hs, mmsErrorsStrongSmoke[#]] &, laneOrder];
weakStrongRatioOrders = AssociationMap[estimateOrders[hs, weakStrongRatios[#]] &, laneOrder];
mmsFinalWeak = AssociationMap[Last[mmsErrorsWeak[#]] &, laneOrder];
mmsFinalWeakSmoke = AssociationMap[Last[mmsErrorsWeakSmoke[#]] &, laneOrder];
mmsFinalStrong = AssociationMap[Last[mmsErrorsStrong[#]] &, laneOrder];
mmsFinalStrongSmoke = AssociationMap[Last[mmsErrorsStrongSmoke[#]] &, laneOrder];
weakStrongFinalRatios = AssociationMap[Last[weakStrongRatios[#]] &, laneOrder];

selfSys = buildSystem[8, 7, 0.37];
selfAsymMatrix = selfSys["WLane"].selfSys["AAction"] - ConjugateTranspose[selfSys["AAction"]].selfSys["WLane"];
selfAsymNorm = N[Norm[Normal[selfAsymMatrix], "Frobenius"]/
    Max[Norm[Normal[selfSys["WLane"].selfSys["AAction"]], "Frobenius"], 10^-30], 17];

gaugeSys = buildSystem[10, 9, 0.37];
chiVals = N[
  MapThread[
    Sin[Pi ((#1 - rMinClosed)/(rExit - rMinClosed))]^2 Sin[Pi #2/wLength]^2 &,
    {gaugeSys["r"], gaugeSys["w"]}
  ],
  17
];
gaugeVec = makeGaugeVector[gaugeSys, chiVals];
fieldMaxGauge = maxAbsFieldStrength[gaugeSys, gaugeVec];
fieldNormGauge = fieldNormFromComponents[gaugeSys, gaugeVec];
gResidualNorm = vectorNormWeighted[gaugeSys["FieldMatrices"]["G"].gaugeVec, gaugeSys["w0"]];
deflatedGaugeVec = deflatedVector[gaugeSys, gaugeVec];
deflatedGaugeNorm = vectorNormWeighted[deflatedGaugeVec, gaugeSys["WLaneValues"]];
deflatedOperatorNorm = vectorNormWeighted[gaugeSys["AAction"].deflatedGaugeVec, gaugeSys["WLaneValues"]];
nondeflatedOperatorNorm = vectorNormWeighted[gaugeSys["AAction"].gaugeVec, gaugeSys["WLaneValues"]];
gaugeGram = ConjugateTranspose[gaugeSys["GaugeMap"]].gaugeSys["WLane"].gaugeSys["GaugeMap"];
gaugeRank = MatrixRank[Normal[gaugeGram], Tolerance -> 10^-10];
gResidualRank = MatrixRank[Normal[gaugeSys["FieldMatrices"]["G"]], Tolerance -> 10^-10];
rawDim = gaugeSys["total"];
physicalQuotientDim = rawDim - gaugeRank - gResidualRank;
anchorCell = Ceiling[gaugeSys["n"]/2];

crossAuditSys = buildSystem[9, 9, 0.0];
crossAuditInput = Join @@ (sampleExpr[#, rSym, wSym, crossAuditSys["r"], crossAuditSys["w"]] & /@ exactExprs);
crossAuditVec = crossAuditSys["CrossBlock"].crossAuditInput;
crossAuditNorms = <|
  "a_r_cross_block_norm" -> vectorNormWeighted[
    crossAuditVec[[crossAuditSys["n"] + 1 ;; 2 crossAuditSys["n"]]],
    crossAuditSys["w0"]
  ],
  "a_E_cross_block_norm" -> vectorNormWeighted[
    crossAuditVec[[2 crossAuditSys["n"] + 1 ;; 3 crossAuditSys["n"]]],
    Lharm crossAuditSys["w0"]
  ],
  "a_B_cross_block_norm" -> vectorNormWeighted[
    crossAuditVec[[3 crossAuditSys["n"] + 1 ;; 4 crossAuditSys["n"]]],
    Lharm crossAuditSys["w0"]
  ]
|>;

aActionMmsPass = And[
  Max[Values[mmsFinalWeak]] < 0.04,
  Min[Last /@ Values[mmsOrdersWeak]] > 1.7
];
aActionTeethPass = And[
  mmsFinalWeakSmoke["a_r"] > 0.05,
  mmsFinalWeakSmoke["a_E"] > 0.05,
  Last[mmsOrdersWeakSmoke["a_r"]] < 0.2,
  Last[mmsOrdersWeakSmoke["a_E"]] < 0.2,
  Max[mmsFinalWeakSmoke /@ cleanSmokeLanes] < 0.04,
  Min[Last /@ (mmsOrdersWeakSmoke /@ cleanSmokeLanes)] > 1.7
];
strongMmsCrossRefPass = And[
  Max[Values[mmsFinalStrong]] < 0.02,
  Min[Flatten[Values[mmsOrdersStrong]]] > 1.7
];
strongTeethCrossRefPass = And[
  mmsFinalStrongSmoke["a_r"] > 0.05,
  mmsFinalStrongSmoke["a_E"] > 0.05,
  Max[mmsOrdersStrongSmoke["a_r"]] < 0.2,
  Max[mmsOrdersStrongSmoke["a_E"]] < 0.2
];
weakStrongConsistencyPass = And @@ Table[
  With[{ratios = weakStrongRatios[lane]},
    And @@ Thread[Rest[ratios] < Most[ratios]] && Last[ratios] < 0.03
  ],
  {lane, laneOrder}
];
mmsGatePass = And[
  aActionMmsPass,
  aActionTeethPass,
  weakStrongConsistencyPass,
  strongMmsCrossRefPass,
  strongTeethCrossRefPass
];

diagnostics = <|
  "schema" -> "mt15_03_spike1_vsh_maxwell_operator/v1",
  "scope" -> "Spike-1 only: operator, gauge deflation, MMS, assembly bookkeeping identity",
  "source_commit" -> sourceCommitString[],
  "target_blind" -> True,
  "geometry" -> <|
    "R0_placeholder" -> "R0(w)=a+(R_exit-a)(3x^2-2x^3), x=w/L",
    "a" -> toJsonNumber[aMouth],
    "R_exit" -> toJsonNumber[rExit],
    "closed_test_r_min" -> toJsonNumber[rMinClosed],
    "L_w" -> toJsonNumber[wLength],
    "Z_placeholder" -> "Z(w)=Z_floor+Z_amp Exp[-((w-w0)/sigma)^2]",
    "Z_floor" -> toJsonNumber[zFloor],
    "Z_amp" -> toJsonNumber[zAmp],
    "Z_center" -> toJsonNumber[zCenter],
    "Z_sigma" -> toJsonNumber[zSigma]
  |>,
  "operator" -> <|
    "ell" -> ell,
    "L" -> Lharm,
    "lanes" -> {"phi", "a_r", "a_E", "a_B", "a_w"},
    "assembly" -> "K=sum sign M^dagger W M from G,E,C,B component matrices; AAction=W_lane^-1 K; primary MMS applies AAction directly to u_exact",
    "closed_boundary_condition" -> "homogeneous Dirichlet/manufactured values on all closed test boundaries",
    "smoke_teeth_operator" -> "ASmoke=AAction-cross_block, deleting D1 vector-Laplacian cross terms from the weak operator"
  |>,
  "mms" -> <|
    "primary_gate" -> "AAction weak Gram assembly against analytic Euler-Lagrange source",
    "secondary_cross_reference" -> "strongOp collocation matrix generated from the Euler variation of the same weak density",
    "grids" -> mmsGrids,
    "h" -> toJsonVector[hs],
    "weak_aaction_residual_norm_region" -> <|
      "type" -> "interior strong-residual rows",
      "dirichlet_closure_layers_excluded_per_side" -> weakMmsClosureLayer,
      "note" -> "The weak Gram operator remains WLaneInv.KAction on the full closed grid. The direct pointwise MMS residual is evaluated away from centered-difference Dirichlet closure rows; the full-grid strongOp cross-reference is reported separately and unchanged."
    |>,
    "aaction" -> <|
      "relative_errors" -> AssociationMap[toJsonVector[mmsErrorsWeak[#]] &, laneOrder],
      "orders" -> AssociationMap[toJsonVector[mmsOrdersWeak[#]] &, laneOrder],
      "final_errors" -> AssociationMap[toJsonNumber[mmsFinalWeak[#]] &, laneOrder],
      "pass" -> aActionMmsPass
    |>,
    "aaction_smoke_cross_terms_dropped" -> <|
      "relative_errors" -> AssociationMap[toJsonVector[mmsErrorsWeakSmoke[#]] &, laneOrder],
      "orders" -> AssociationMap[toJsonVector[mmsOrdersWeakSmoke[#]] &, laneOrder],
      "final_errors" -> AssociationMap[toJsonNumber[mmsFinalWeakSmoke[#]] &, laneOrder],
      "teeth_pass" -> aActionTeethPass
    |>,
    "strongop_secondary_cross_reference" -> <|
      "relative_errors" -> AssociationMap[toJsonVector[mmsErrorsStrong[#]] &, laneOrder],
      "orders" -> AssociationMap[toJsonVector[mmsOrdersStrong[#]] &, laneOrder],
      "final_errors" -> AssociationMap[toJsonNumber[mmsFinalStrong[#]] &, laneOrder],
      "pass" -> strongMmsCrossRefPass
    |>,
    "strongop_smoke_cross_terms_dropped_secondary" -> <|
      "relative_errors" -> AssociationMap[toJsonVector[mmsErrorsStrongSmoke[#]] &, laneOrder],
      "orders" -> AssociationMap[toJsonVector[mmsOrdersStrongSmoke[#]] &, laneOrder],
      "final_errors" -> AssociationMap[toJsonNumber[mmsFinalStrongSmoke[#]] &, laneOrder],
      "teeth_pass" -> strongTeethCrossRefPass
    |>,
    "weak_vs_strong_consistency" -> <|
      "ratio" -> "norm(AAction.u_exact-strongOp.u_exact)/norm(strongOp.u_exact), by lane",
      "ratios" -> AssociationMap[toJsonVector[weakStrongRatios[#]] &, laneOrder],
      "orders" -> AssociationMap[toJsonVector[weakStrongRatioOrders[#]] &, laneOrder],
      "final_ratios" -> AssociationMap[toJsonNumber[weakStrongFinalRatios[#]] &, laneOrder],
      "pass" -> weakStrongConsistencyPass
    |>,
    "gate_pass" -> mmsGatePass
  |>,
  "gauge_deflation" -> <|
    "grid" -> {gaugeSys["nr"], gaugeSys["nw"]},
    "omega" -> toJsonNumber[gaugeSys["omega"]],
    "raw_dim" -> rawDim,
    "pure_gauge_rank" -> gaugeRank,
    "lorenz_residual_rank" -> gResidualRank,
    "scalar_anchor" -> <|"type" -> "chi(anchor_cell)=0", "cell" -> anchorCell|>,
    "physical_quotient_dim_G_equals_0_mod_gauge" -> physicalQuotientDim,
    "expected_physical_quotient_dim_3N" -> 3 gaugeSys["n"],
    "pure_gauge_max_abs_F" -> toJsonNumber[fieldMaxGauge],
    "pure_gauge_physical_field_norm" -> toJsonNumber[fieldNormGauge],
    "pure_gauge_G_residual_norm_before_constraint" -> toJsonNumber[gResidualNorm],
    "nondeflated_gauge_fixed_operator_norm" -> toJsonNumber[nondeflatedOperatorNorm],
    "deflated_vector_weighted_norm" -> toJsonNumber[deflatedGaugeNorm],
    "deflated_operator_weighted_norm" -> toJsonNumber[deflatedOperatorNorm]
  |>,
  "assembly_hermiticity_identity_construction_guaranteed" -> <|
    "note" -> "WLane.AAction = KAction = Sum[M^dagger W M] is Hermitian by construction; this verifies WLaneInv bookkeeping, NOT operator physics. Operator correctness is covered by the primary AAction MMS.",
    "grid" -> {selfSys["nr"], selfSys["nw"]},
    "omega" -> toJsonNumber[selfSys["omega"]],
    "relative_W_asymmetry_frobenius" -> toJsonNumber[selfAsymNorm]
  |>,
  "cross_term_audit" -> Join[
    <|
      "D1_cross_block" -> {
        "(Delta_3 A)_r contains -2 a_r/r^2 + 2 L a_E/r^2",
        "(Delta_3 A)_E contains +2 a_r/r^2",
        "(Delta_3 A)_B has no cross block"
      }
    |>,
    Map[toJsonNumber, crossAuditNorms]
  ],
  "term_map" -> {
    <|"term" -> "parent localized Maxwell action", "source" -> "notes/moving_throat_pde_program_compact.md:592-620", "implementation" -> "weighted Z(w) weak density with H=Z gauge fixing"|>,
    <|"term" -> "exact localized Maxwell equation", "source" -> "notes/moving_throat_pde_program_compact.md:674-689", "implementation" -> "operator is the Hessian of the VSH-reduced weak form, A=W^-1 K"|>,
    <|"term" -> "mixed gauge invariants E_w,C_a", "source" -> "notes/moving_throat_pde_program_compact.md:769-786", "implementation" -> "E_w=i omega a_w-d_w phi; C_r/C_E/C_B include d_w vector lanes"|>,
    <|"term" -> "five-lane VSH decomposition", "source" -> "software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9829-9843 and decisions/05_mixed_maxwell_spike_design.md:12-19", "implementation" -> "lane order (phi,a_r,a_E,a_B,a_w)"|>,
    <|"term" -> "G,E,C,B reduced components", "source" -> "software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9845-9862", "implementation" -> "component matrices mG,mEr,mEE,mEB,mEw,mCr,mCE,mCB,mBB,mBr,mBE"|>,
    <|"term" -> "quadratic VSH weak form and angular weights", "source" -> "software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9864-9872 and decisions/05_mixed_maxwell_spike_design.md:14-16", "implementation" -> "K=sum sign M^dagger W M with L weights on Psi/Phi components"|>,
    <|"term" -> "pure-gauge quotient", "source" -> "software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9874-9879", "implementation" -> "GaugeMap chi=(i omega chi,d_r chi,chi/r,0,d_w chi); W-orthogonal deflation"|>,
    <|"term" -> "vector-Laplacian cross terms", "source" -> "software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9881-9890 and decisions/05_mixed_maxwell_spike_design.md:27-29", "implementation" -> "CrossBlock rows a_r,a_E delete-test block; smoke A=A-crossBlock"|>,
    <|"term" -> "closed Spike-1 BC and physical norm", "source" -> "software/stage1_solver/_scratch/mixed_maxwell_spike_design.log:9912-9925 and software/stage1_solver/directives/mt15_03_spike1_vsh_maxwell_operator.md:29-48", "implementation" -> "closed Dirichlet MMS/assembly identity; positive F norm for quotient diagnostics"|>
  }
|>;

Export[diagnosticsPath, diagnostics, "JSON"];

weakOrdersRows = Table[
  {
    lane,
    StringRiffle[formatFixed /@ mmsOrdersWeak[lane], ", "],
    formatScientific[mmsFinalWeak[lane]],
    StringRiffle[formatFixed /@ mmsOrdersWeakSmoke[lane], ", "],
    formatScientific[mmsFinalWeakSmoke[lane]]
  },
  {lane, laneOrder}
];

strongOrdersRows = Table[
  {
    lane,
    StringRiffle[formatFixed /@ mmsOrdersStrong[lane], ", "],
    formatScientific[mmsFinalStrong[lane]],
    StringRiffle[formatFixed /@ mmsOrdersStrongSmoke[lane], ", "],
    formatScientific[mmsFinalStrongSmoke[lane]]
  },
  {lane, laneOrder}
];

consistencyOrdersRows = Table[
  {
    lane,
    StringRiffle[formatFixed /@ weakStrongRatios[lane], ", "],
    StringRiffle[formatFixed /@ weakStrongRatioOrders[lane], ", "],
    formatScientific[weakStrongFinalRatios[lane]]
  },
  {lane, laneOrder}
];

weakMmsTable = markdownTable[
  Join[{"nr", "nw", "h"}, ("AAction " <> #) & /@ laneOrder, ("ASmoke " <> #) & /@ laneOrder],
  weakMmsRows
];
weakOrdersTable = markdownTable[
  {"lane", "AAction orders", "AAction final rel err", "ASmoke orders", "ASmoke final rel err"},
  weakOrdersRows
];
strongMmsTable = markdownTable[
  Join[{"nr", "nw", "h"}, ("strongOp " <> #) & /@ laneOrder, ("strong smoke " <> #) & /@ laneOrder],
  strongMmsRows
];
strongOrdersTable = markdownTable[
  {"lane", "strongOp orders", "strongOp final rel err", "strong smoke orders", "strong smoke final rel err"},
  strongOrdersRows
];
consistencyTable = markdownTable[
  Join[{"nr", "nw", "h"}, ("ratio " <> #) & /@ laneOrder],
  consistencyRows
];
consistencyOrdersTable = markdownTable[
  {"lane", "weak-vs-strong ratios", "ratio orders", "final ratio"},
  consistencyOrdersRows
];
termMapTable = markdownTable[
  {"operator term", "implementation", "parent/design source"},
  ({#["term"], #["implementation"], #["source"]} & /@ diagnostics["term_map"])
];

report = StringRiffle[
  {
    "# MT15-03 Spike-1 VSH Maxwell Operator Report",
    "",
    "Status: Spike-1 only. No outgoing DtN/open-impedance solve, no wall-to-current transfer, no N0/Z0 extraction, and no V2 packet emission.",
    "",
    "## Five-lane operator",
    "",
    "Lane order: `(phi, a_r, a_E, a_B, a_w)` with `ell=2`, `L=ell(ell+1)=6`.",
    "",
    "Target-blind placeholders: `R0(w)=a+(R_exit-a)(3x^2-2x^3)` with `a=1`, `R_exit=1.65`, and `Z(w)=Z_floor+Z_amp Exp[-((w-w0)/sigma)^2]` with `Z_floor=0.08`, `Z_amp=0.92`, `w0=L_w/2`, `sigma=0.38 L_w`.  These are geometry/localization placeholders only; the closed test grid is rectangular and avoids the origin (`r_min=0.35`) so Spike-1 can isolate the VSH operator without origin-parity or open-exit machinery.",
    "",
    "The assembled action operator is `AAction=W_lane^-1 K`, where `K=sum sign M^dagger W M` over the D1 component matrices `G`, `E_r/E_E/E_B/E_w`, `C_r/C_E/C_B`, and `B_B/B_r/B_E`.  The primary MMS applies this weak Gram assembly directly to `u_exact`; the collocation `strongOp` generated from the Euler variation is retained as a secondary cross-reference.  The physical quotient norm uses the positive field-strength form `|E|^2+|B|^2+|C|^2`.",
    "",
    "Cross-term carrying lanes: `a_r` row contains `-2 a_r/r^2 + 2 L a_E/r^2`; `a_E` row contains `+2 a_r/r^2`; `a_B` has the scalar vector-harmonic block only.  The smoke test deletes exactly this cross block.",
    "",
    "## MMS",
    "",
    "Primary gate: `AAction.u_exact-f_exact` is measured directly against the analytic Euler-Lagrange source on the interior rows where the centered weak Gram assembly represents a pointwise strong residual; two Dirichlet closure layers per side are excluded from this pointwise norm.  The full-grid `strongOp` cross-reference below is unchanged.  Manufactured solution has all five lanes nonzero and has coupled nonzero `a_r` and `a_E`; the source is derived by Euler variation of the same VSH weak density used by the matrix assembly.",
    "",
    weakMmsTable,
    "",
    weakOrdersTable,
    "",
    "Weak-operator teeth result: `aaction_teeth_pass=" <> ToString[aActionTeethPass] <> "`.  The genuine weak operator converges; deleting the D1 cross block from `AAction` leaves an O(1) residual on `a_r`/`a_E`, while `phi`/`a_B`/`a_w` remain second-order.",
    "",
    "Weak-vs-strong consistency ratio `||AAction.u_exact-strongOp.u_exact||/||strongOp.u_exact||` decreases under refinement:",
    "",
    consistencyTable,
    "",
    consistencyOrdersTable,
    "",
    "Secondary `strongOp` cross-reference MMS (kept to preserve the original collocation Euler check):",
    "",
    strongMmsTable,
    "",
    strongOrdersTable,
    "",
    "## Gauge Deflation",
    "",
    "- Gauge test grid: " <> ToString[gaugeSys["nr"]] <> "x" <> ToString[gaugeSys["nw"]] <> ", raw dimension `" <> ToString[rawDim] <> "`.",
    "- Pure-gauge rank: `" <> ToString[gaugeRank] <> "`; Lorenz/Gauss residual rank: `" <> ToString[gResidualRank] <> "`.",
    "- Scalar gauge anchor: `chi(anchor_cell)=0` at cell `" <> ToString[anchorCell] <> "`.",
    "- Physical quotient dimension (`G=0` modulo pure gauge): `" <> ToString[physicalQuotientDim] <> "`; expected `3N=" <> ToString[3 gaugeSys["n"]] <> "`.",
    "- Pure-gradient max `|F|`: `" <> formatScientific[fieldMaxGauge] <> "`; positive physical field norm: `" <> formatScientific[fieldNormGauge] <> "`.",
    "- Non-deflated gauge-fixed operator norm: `" <> formatScientific[nondeflatedOperatorNorm] <> "`; after W-orthogonal gauge deflation, vector norm `" <> formatScientific[deflatedGaugeNorm] <> "`, operator norm `" <> formatScientific[deflatedOperatorNorm] <> "`.",
    "",
    "The nonzero pre-constraint `G` norm (`" <> formatScientific[gResidualNorm] <> "`) is expected: this is a pure gauge field-strength test, not a Lorenz-satisfying representative.  The physical quotient deflates it.",
    "",
    "## Assembly Hermiticity Identity",
    "",
    "`WLane.AAction = KAction = sum M^dagger W M` is Hermitian by construction, so this is an assembly/bookkeeping identity check for `WLaneInv`, not a physics gate.  Operator correctness is covered by the primary `AAction` MMS.  At `omega=0.37`, relative Frobenius norm of `W A - A^dagger W` is `" <> formatScientific[selfAsymNorm] <> "`.",
    "",
    "## Cross-Term Audit",
    "",
    "- `a_r` cross-block norm on the MMS profile: `" <> formatScientific[crossAuditNorms["a_r_cross_block_norm"]] <> "`.",
    "- `a_E` cross-block norm on the MMS profile: `" <> formatScientific[crossAuditNorms["a_E_cross_block_norm"]] <> "`.",
    "- `a_B` cross-block norm on the MMS profile: `" <> formatScientific[crossAuditNorms["a_B_cross_block_norm"]] <> "`.",
    "",
    "## Term Map",
    "",
    termMapTable,
    "",
    "Diagnostics JSON: `" <> diagnosticsPath <> "`."
  },
  "\n"
];

Export[reportPath, report, "Text"];

Print["diagnostics_path: ", diagnosticsPath];
Print["report_path: ", reportPath];
Print["assembly_hermiticity_identity_construction_guaranteed_relative_W_asymmetry: ", formatScientific[selfAsymNorm]];
Print["physical_quotient_dim: ", physicalQuotientDim, " expected_3N: ", 3 gaugeSys["n"]];
Print["pure_gauge_max_abs_F: ", formatScientific[fieldMaxGauge]];
Print["pure_gauge_physical_norm: ", formatScientific[fieldNormGauge]];
Print["deflated_operator_norm: ", formatScientific[deflatedOperatorNorm]];
Print["AAction MMS final max relative error: ", formatScientific[Max[Values[mmsFinalWeak]]]];
Print["AAction MMS orders by lane: ", AssociationMap[formatFixed /@ mmsOrdersWeak[#] &, laneOrder]];
Print["AAction smoke final a_r/a_E relative errors: ",
  formatScientific[mmsFinalWeakSmoke["a_r"]], " / ", formatScientific[mmsFinalWeakSmoke["a_E"]]];
Print["weak_vs_strong_final_ratios: ", AssociationMap[formatScientific[weakStrongFinalRatios[#]] &, laneOrder]];
Print["strongOp secondary MMS final max relative error: ", formatScientific[Max[Values[mmsFinalStrong]]]];
Print["aaction_mms_pass: ", aActionMmsPass];
Print["aaction_teeth_pass: ", aActionTeethPass];
Print["weak_vs_strong_consistency_pass: ", weakStrongConsistencyPass];
Print["strongop_secondary_cross_reference_pass: ", strongMmsCrossRefPass];
Print["mms_gate_pass: ", mmsGatePass];

If[! TrueQ[mmsGatePass],
  Print["ERROR: primary AAction MMS, weak teeth, weak-vs-strong consistency, or secondary strongOp cross-reference failed thresholds."];
  Exit[5]
];
If[selfAsymNorm > 10^-12,
  Print["ERROR: assembly hermiticity identity asymmetry exceeds tolerance."];
  Exit[6]
];
If[fieldMaxGauge > 10^-10 || fieldNormGauge > 10^-10 || deflatedOperatorNorm > 10^-8,
  Print["ERROR: pure-gauge deflation check exceeds tolerance."];
  Exit[7]
];
If[physicalQuotientDim =!= 3 gaugeSys["n"],
  Print["ERROR: physical quotient dimension mismatch."];
  Exit[8]
];

Print["MT15-03 Spike-1 checks passed."];
