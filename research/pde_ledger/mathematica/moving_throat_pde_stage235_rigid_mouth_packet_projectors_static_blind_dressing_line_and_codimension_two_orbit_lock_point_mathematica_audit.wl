ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

fmt[expr_] := StringReplace[ToString[InputForm[expr]], "Global`" -> ""];
pass[name_String] := Print["PASS: ", name];
fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanExpr[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectMatrixZero[name_String, mat_] := Module[{depth, res},
  depth = ArrayDepth[mat];
  res = Map[cleanExpr, mat, {depth}];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === ConstantArray[0, Dimensions[res]]],
    pass[name],
    fail[name, res]
  ];
];

expectTrue[name_String, cond_] := Module[{res},
  res = stripConditional[FullSimplify[cond, Assumptions -> $Assumptions]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

solveRules[eqns_] := stripConditional[
  FullSimplify[Solve[eqns, {R1, E1}, Reals], Assumptions -> $Assumptions]
];

$Assumptions = (
  cEta > 0 && L > 0 &&
  Element[{R1, E1, qNt, qEta}, Reals]
);

banner["Stage 235 Mathematica audit: rigid-mouth packet projectors and orbit lock"];

mouthCompiler = {{-1, -cEta}, {0, 1}};
directState = {R1, E1};
packetState = {qNt, qEta};
xiPacket = -R1 - cEta E1;

subbanner["M1. Involutive rigid-mouth compiler"];

compiledDirect = mouthCompiler . directState;
inversePacketMap = mouthCompiler . packetState;

expectMatrixZero["M1 M_rm.{R1,E1} - {Xi1,E1}", compiledDirect - {xiPacket, E1}];
expectMatrixZero["M1 explicit vector target", compiledDirect - {-R1 - cEta E1, E1}];
expectMatrixZero["M1 MatrixPower[M_rm,2] - I", MatrixPower[mouthCompiler, 2] - IdentityMatrix[2]];
expectMatrixZero["M1 inverse packet map", inversePacketMap - {-qNt - cEta qEta, qEta}];
expectZero["M1 Det[M_rm] + 1", Det[mouthCompiler] + 1];

Print["M_rm = ", fmt[mouthCompiler]];
Print["M_rm . {R1,E1} = ", fmt[compiledDirect]];
Print["M_rm . {q_nt,q_eta} = ", fmt[inversePacketMap]];

subbanner["M2. Direct-space projectors from similarity"];

axisNt = DiagonalMatrix[{1, 0}];
axisEta = DiagonalMatrix[{0, 1}];
projNt = mouthCompiler . axisNt . mouthCompiler;
projEta = mouthCompiler . axisEta . mouthCompiler;
zeroTwo = ConstantArray[0, {2, 2}];

expectMatrixZero["M2 P_nt closed form", projNt - {{1, cEta}, {0, 0}}];
expectMatrixZero["M2 P_eta closed form", projEta - {{0, -cEta}, {0, 1}}];
expectMatrixZero["M2 P_nt idempotent", projNt . projNt - projNt];
expectMatrixZero["M2 P_eta idempotent", projEta . projEta - projEta];
expectMatrixZero["M2 P_nt.P_eta zero", projNt . projEta - zeroTwo];
expectMatrixZero["M2 P_eta.P_nt zero", projEta . projNt - zeroTwo];
expectMatrixZero["M2 projector sum identity", projNt + projEta - IdentityMatrix[2]];

Print["P_nt = ", fmt[projNt]];
Print["P_eta = ", fmt[projEta]];

subbanner["M3. Direct-space packet decomposition"];

ntPiece = projNt . directState;
etaPiece = projEta . directState;

expectMatrixZero["M3 P_nt.x_rm closed form", ntPiece - {R1 + cEta E1, 0}];
expectMatrixZero["M3 P_nt.x_rm equals {-Xi1,0}", ntPiece - {-xiPacket, 0}];
expectMatrixZero["M3 P_eta.x_rm closed form", etaPiece - {-cEta E1, E1}];
expectMatrixZero["M3 x_rm recomposes", directState - ntPiece - etaPiece];

Print["P_nt . x_rm = ", fmt[ntPiece]];
Print["P_eta . x_rm = ", fmt[etaPiece]];

subbanner["M4. Codimension-two orbit lock"];

lockFromQ = solveRules[Thread[compiledDirect == {0, 0}]];
originFromDirect = mouthCompiler . {0, 0};
solXiE = solveRules[{xiPacket == 0, E1 == 0}];
solXiR = solveRules[{xiPacket == 0, R1 == 0}];

expectTrue["M4 Det[M_rm] != 0", Det[mouthCompiler] != 0];
expectTrue["M4 q_nt=q_eta=0 implies R1=E1=0", lockFromQ === {{R1 -> 0, E1 -> 0}}];
expectMatrixZero["M4 R1=E1=0 implies q_nt=q_eta=0", originFromDirect - {0, 0}];
expectTrue["M4 Solve[{Xi1==0,E1==0}] unique origin", solXiE === {{R1 -> 0, E1 -> 0}}];
expectTrue["M4 Solve[{Xi1==0,R1==0}] unique origin", solXiR === {{R1 -> 0, E1 -> 0}}];
expectTrue["M4 static strip Xi1==0 has one free direct coordinate", Length[Solve[{xiPacket == 0}, {R1}, Reals]] == 1];

Print["Solve[q_rm==0] = ", fmt[lockFromQ]];
Print["Solve[Xi1==0,E1==0] = ", fmt[solXiE]];
Print["Solve[Xi1==0,R1==0] = ", fmt[solXiR]];

subbanner["M5. Static-blind dressing line and exact norm"];

blindLine = {-cEta qEta, qEta};
blindXi = -blindLine[[1]] - cEta blindLine[[2]];
blindNormSquared = blindLine . blindLine;
lengthRule = qEta -> L/Sqrt[1 + cEta^2];

expectZero["M5 Xi1 on blind line", blindXi];
expectZero["M5 blind-line norm squared", blindNormSquared - (1 + cEta^2) qEta^2];
expectZero["M5 norm with q_eta = L/Sqrt[1+c_eta^2]", (blindNormSquared /. lengthRule) - L^2];
expectZero["M5 Xi1 remains zero at arbitrary size L", blindXi /. lengthRule];

Print["x_blind = ", fmt[blindLine]];
Print["x_blind.x_blind = ", fmt[cleanExpr[blindNormSquared]]];

subbanner["M6. Correction compilers"];

packetDirect = mouthCompiler . packetState;
staticCorrection = -mouthCompiler . axisNt . packetState;
orbitCorrection = -packetDirect;

expectMatrixZero["M6 x_q direct compiler", packetDirect - {-qNt - cEta qEta, qEta}];
expectMatrixZero["M6 Delta_x_static", staticCorrection - {qNt, 0}];
expectMatrixZero["M6 Delta_x_orbit", orbitCorrection - {qNt + cEta qEta, -qEta}];
expectMatrixZero[
  "M6 additive static-to-orbit relation",
  orbitCorrection - (staticCorrection + {cEta qEta, -qEta})
];
expectMatrixZero["M6 full-lock identity", packetDirect + orbitCorrection];

Print["Delta_x_static = ", fmt[staticCorrection]];
Print["Delta_x_orbit = ", fmt[orbitCorrection]];

Print[""];
Print["All Stage 235 Mathematica checks passed."];
