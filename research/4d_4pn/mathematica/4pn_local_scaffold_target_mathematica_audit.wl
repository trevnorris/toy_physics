#!/usr/bin/env wolframscript

ClearAll["Global`*"];
passCount = 0;
failCount = 0;

Get[FileNameJoin[{DirectoryName[$InputFileName], "4pn_shared.wl"}]];

oneBodyHamiltonianGate4PN[] := <|
  "K" -> 7/256,
  "Q1" -> 45/128,
  "T1" -> 13/8,
  "S1" -> 105/32,
  "U1" -> 105/32,
  "W1" -> -1/16
|>;

banner["PART I — RAW LOCAL 4PN GENERIC-FRAME SCAFFOLD"];

qBasis = generateBasis4PN[0, 8];
tBasis = generateBasis4PN[1, 6];
sBasis = generateBasis4PN[2, 4];
uBasis = generateBasis4PN[3, 2];
wBasis = generateBasis4PN[4, 0];

qM = imageMatrixPolynomial4PN["Q", qBasis, 4];
tM = imageMatrixPolynomial4PN["T", tBasis, 3];
sM = imageMatrixPolynomial4PN["S", sBasis, 3];
uM = imageMatrixPolynomial4PN["U", uBasis, 2];
wM = imageMatrixPolynomial4PN["W", wBasis, 2];

checkZero["Q basis size - 52", Length[qBasis] - 52];
checkZero["T basis size - 46", Length[tBasis] - 46];
checkZero["S basis size - 29", Length[sBasis] - 29];
checkZero["U basis size - 10", Length[uBasis] - 10];
checkZero["W basis size - 2", Length[wBasis] - 2];

checkZero["Q image rank - 5", MatrixRank[qM] - 5];
checkZero["T image rank - 4", MatrixRank[tM] - 4];
checkZero["S image rank - 3", MatrixRank[sM] - 3];
checkZero["U image rank - 2", MatrixRank[uM] - 2];
checkZero["W image rank - 1", MatrixRank[wM] - 1];

checkZero["Q nullity - 47", Length[qBasis] - MatrixRank[qM] - 47];
checkZero["T nullity - 42", Length[tBasis] - MatrixRank[tM] - 42];
checkZero["S nullity - 26", Length[sBasis] - MatrixRank[sM] - 26];
checkZero["U nullity - 8", Length[uBasis] - MatrixRank[uM] - 8];
checkZero["W nullity - 1", Length[wBasis] - MatrixRank[wM] - 1];
checkZero["total raw directions - 139", Total[Length /@ {qBasis, tBasis, sBasis, uBasis, wBasis}] - 139];
checkZero["total COM rank - 15", Total[MatrixRank /@ {qM, tM, sM, uM, wM}] - 15];
checkZero["total COM nullity - 124", (Total[Length /@ {qBasis, tBasis, sBasis, uBasis, wBasis}] - Total[MatrixRank /@ {qM, tM, sM, uM, wM}]) - 124];

banner["PART II — STRICT ONE-BODY 4PN HAMILTONIAN GATE AND IMPORTED TARGET"];

gate = oneBodyHamiltonianGate4PN[];
target = localHamiltonianTarget4PN[];

checkZero["free slot K at nu->0", (target[1] /. nuSym -> 0) - gate["K"]];
checkZero["G/r leading slot Q1 at nu->0", (target[7] /. nuSym -> 0) - gate["Q1"]];
checkZero["G^2/r^2 leading slot T1 at nu->0", (target[12] /. nuSym -> 0) - gate["T1"]];
checkZero["G^3/r^3 leading slot S1 at nu->0", (target[16] /. nuSym -> 0) - gate["S1"]];
checkZero["G^4/r^4 leading slot U1 at nu->0", (target[19] /. nuSym -> 0) - gate["U1"]];
checkZero["G^5/r^5 leading slot W1 at nu->0", (target[21] /. nuSym -> 0) - gate["W1"]];

Do[
  checkZero["subleading one-body slot " <> ToString[i], target[i] /. nuSym -> 0],
  {i, {8, 9, 10, 11, 13, 14, 15, 17, 18, 20}}
];

checkZero["G^4/r^4 max nu-degree - 2", Max[Exponent[Expand[target[#]], nuSym] & /@ {19, 20}] - 2];
checkZero["G^5/r^5 max nu-degree - 2", Exponent[Expand[target[21]], nuSym] - 2];

banner["FINAL LOCAL SCAFFOLD / TARGET-IMPORT LEDGER"];
finalizeAudit[];

(*"
Output:

========================================================================================
PART I — RAW LOCAL 4PN GENERIC-FRAME SCAFFOLD
========================================================================================
Q basis size - 52 = 0
T basis size - 46 = 0
S basis size - 29 = 0
U basis size - 10 = 0
W basis size - 2 = 0
Q image rank - 5 = 0
T image rank - 4 = 0
S image rank - 3 = 0
U image rank - 2 = 0
W image rank - 1 = 0
Q nullity - 47 = 0
T nullity - 42 = 0
S nullity - 26 = 0
U nullity - 8 = 0
W nullity - 1 = 0
total raw directions - 139 = 0
total COM rank - 15 = 0
total COM nullity - 124 = 0

========================================================================================
PART II — STRICT ONE-BODY 4PN HAMILTONIAN GATE AND IMPORTED TARGET
========================================================================================
free slot K at nu->0 = 0
G/r leading slot Q1 at nu->0 = 0
G^2/r^2 leading slot T1 at nu->0 = 0
G^3/r^3 leading slot S1 at nu->0 = 0
G^4/r^4 leading slot U1 at nu->0 = 0
G^5/r^5 leading slot W1 at nu->0 = 0
subleading one-body slot 8 = 0
subleading one-body slot 9 = 0
subleading one-body slot 10 = 0
subleading one-body slot 11 = 0
subleading one-body slot 13 = 0
subleading one-body slot 14 = 0
subleading one-body slot 15 = 0
subleading one-body slot 17 = 0
subleading one-body slot 18 = 0
subleading one-body slot 20 = 0
G^4/r^4 max nu-degree - 2 = 0
G^5/r^5 max nu-degree - 2 = 0

========================================================================================
FINAL LOCAL SCAFFOLD / TARGET-IMPORT LEDGER
========================================================================================
Passes: 36
Fails: 0
"*)
