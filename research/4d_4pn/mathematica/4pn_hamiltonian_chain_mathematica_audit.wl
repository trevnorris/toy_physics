#!/usr/bin/env wolframscript

ClearAll["Global`*"];
passCount = 0;
failCount = 0;

Get[FileNameJoin[{DirectoryName[$InputFileName], "4pn_shared.wl"}]];

banner["PART I — EXACT 21-SLOT QUARTIC COM MAP"];

mapData = quarticComMap4PN[];
hmap = mapData["map"];
coeffs = mapData["coeffs"];
jac = Table[D[hmap[i], coeffs[[j]]], {i, 1, 21}, {j, 1, 21}];
checkTrue["Jacobian + I is identically zero", And @@ Flatten[Map[FullSimplify[# === 0] &, jac + IdentityMatrix[21], {2}]]];

checkZero["h1 map identity", hmap[1] + coeffs[[1]] + (2217 nuSym^4 - 1936 nuSym^3 + 473 nuSym^2 - 4 nuSym - 7)/128];
checkZero["h7 map identity", hmap[7] + coeffs[[7]] + (-625 nuSym^4 - 517 nuSym^3 + 350 nuSym^2 + 13 nuSym - 15)/16];
checkZero["h12 map identity", hmap[12] + coeffs[[12]] + (938 nuSym^4 + 2510 nuSym^3 - 195 nuSym^2 + 438 nuSym - 170)/32];
checkZero["h16 map identity", hmap[16] + coeffs[[16]] - (1584 nuSym^4 + 5880 nuSym^3 + 5556 nuSym^2 + 9 Pi^2 nuSym^2 - 3 Pi^2 nuSym - 9868 nuSym + 1848)/192];
checkZero["h19 map identity", hmap[19] + coeffs[[19]] + (48 nuSym^4 + 108 nuSym^3 - 239 nuSym^2 + 3 Pi^2 nuSym^2 + 9 Pi^2 nuSym + 840 nuSym - 408)/96];
checkZero["h21 map identity", hmap[21] + coeffs[[21]]];

banner["PART II — EXACT ORDINARY LOCAL 4PN TARGET"];

ordData = ordinaryTargetFromHamiltonian4PN[];
feedback = ordData["feedback"];
targetH = ordData["targetH"];
targetL = ordData["targetL"];

checkZero["L1 target formula", targetL[1] - (-(4497 nuSym^4 - 4082 nuSym^3 + 1135 nuSym^2 - 71 nuSym - 7)/256)];
checkZero["L7 target formula", targetL[7] - ((10070 nuSym^4 + 9285 nuSym^3 - 7292 nuSym^2 + 512 nuSym + 150)/256)];
checkZero["L12 target formula", targetL[12] - (-(7504 nuSym^4 + 22415 nuSym^3 + 3297 nuSym^2 + 340 nuSym - 944)/256)];
checkZero["L16 target formula", targetL[16] - ((30412800 nuSym^4 + 128822400 nuSym^3 - 3987675 Pi^2 nuSym^2 + 258968192 nuSym^2 - 76341312 nuSym - 1294650 Pi^2 nuSym + 23385600)/3686400)];
checkZero["L19 target formula", targetL[19] + (614400 nuSym^4 + 1382400 nuSym^3 - 3916025 Pi^2 nuSym^2 + 40000704 nuSym^2 - 3160350 Pi^2 nuSym + 22640704 nuSym - 1190400)/1228800];
checkZero["L21 target formula", targetL[21] + (-6430720 nuSym^2 + 555225 Pi^2 nuSym^2 - 16243104 nuSym + 1403325 Pi^2 nuSym - 14400)/230400];

checkZero["H target G/r degree ceiling - 4", Max[Exponent[Expand[targetH[#]], nuSym] & /@ Range[7, 11]] - 4];
checkZero["H target G^2/r^2 degree ceiling - 3", Max[Exponent[Expand[targetH[#]], nuSym] & /@ Range[12, 15]] - 3];
checkZero["H target G^3/r^3 degree ceiling - 3", Max[Exponent[Expand[targetH[#]], nuSym] & /@ Range[16, 18]] - 3];
checkZero["H target G^4/r^4 degree ceiling - 2", Max[Exponent[Expand[targetH[#]], nuSym] & /@ Range[19, 20]] - 2];
checkZero["H target G^5/r^5 degree ceiling - 2", Exponent[Expand[targetH[21]], nuSym] - 2];

checkZero["L target G/r degree ceiling - 4", Max[Exponent[Expand[targetL[#]], nuSym] & /@ Range[7, 11]] - 4];
checkZero["L target G^2/r^2 degree ceiling - 4", Max[Exponent[Expand[targetL[#]], nuSym] & /@ Range[12, 15]] - 4];
checkZero["L target G^3/r^3 degree ceiling - 4", Max[Exponent[Expand[targetL[#]], nuSym] & /@ Range[16, 18]] - 4];
checkZero["L target G^4/r^4 degree ceiling - 4", Max[Exponent[Expand[targetL[#]], nuSym] & /@ Range[19, 20]] - 4];
checkZero["L target G^5/r^5 degree ceiling - 2", Exponent[Expand[targetL[21]], nuSym] - 2];

banner["PART III — HAMILTONIAN-CHART GENERIC-FRAME LIFT AND CANONICAL SLICE"];

hData = canonicalHamiltonianData4PN[];
basis = hData["basis"];
mats = hData["matrices"];
coords = hData["coords"];
blocks = hData["blocks"];
residualSlots = hData["residualSlots"];
nullspaces = hData["nullspaces"];

checkZero["Q rank - 20", MatrixRank[mats["Q"]] - 20];
checkZero["T rank - 12", MatrixRank[mats["T"]] - 12];
checkZero["S rank - 9", MatrixRank[mats["S"]] - 9];
checkZero["U rank - 4", MatrixRank[mats["U"]] - 4];
checkZero["W rank - 2", MatrixRank[mats["W"]] - 2];

checkZero["Q nullity - 32", Length[basis["Q"]] - MatrixRank[mats["Q"]] - 32];
checkZero["T nullity - 34", Length[basis["T"]] - MatrixRank[mats["T"]] - 34];
checkZero["S nullity - 20", Length[basis["S"]] - MatrixRank[mats["S"]] - 20];
checkZero["U nullity - 6", Length[basis["U"]] - MatrixRank[mats["U"]] - 6];
checkZero["W nullity - 0", Length[basis["W"]] - MatrixRank[mats["W"]]];

Do[
  checkTrue["canonical solve " <> block <> " is exact",
    And @@ Map[
      FullSimplify[# === 0] &,
      Expand[mats[block].coords[block] - targetVector4PN[residualSlots[block], <|"Q" -> 4, "T" -> 3, "S" -> 3, "U" -> 2, "W" -> 2|>[block]]]
    ]
  ],
  {block, {"Q", "T", "S", "U", "W"}}
];

Do[
  checkTrue["canonical COM verification " <> block,
    And @@ MapThread[FullSimplify[Expand[#1 - #2]] === 0 &, {blockSlots4PN[toNu4PN[blocks[block]], block], residualSlots[block]}]
  ],
  {block, {"Q", "T", "S", "U", "W"}}
];

linIdeal = {pp aa + qq cc, pp cc + qq bb, pp dd + qq ee};
fullIdeal = Join[linIdeal, {aa bb - cc^2, aa ee - cc dd, bb dd - cc ee}];

Do[
  checkTrue["Q-null polynomial in full COM ideal " <> ToString[i],
    Last[PolynomialReduce[canonicalExpr4PN[nullspaces["Q"][[i]], basis["Q"]], fullIdeal, {pp, qq, aa, bb, cc, dd, ee}]] === 0
  ],
  {i, 1, Length[nullspaces["Q"]]}
];
Do[
  checkTrue[block <> "-null polynomial in linear COM ideal " <> ToString[i],
    Last[PolynomialReduce[canonicalExpr4PN[nullspaces[block][[i]], basis[block]], linIdeal, {pp, qq, aa, bb, cc, dd, ee}]] === 0
  ],
  {block, {"T", "S", "U"}}, {i, 1, Length[nullspaces[block]]}
];
checkZero["W nullspace dimension", Length[nullspaces["W"]]];

banner["FINAL HAMILTONIAN 4PN CHAIN LEDGER"];
finalizeAudit[];

(*"
Output:

========================================================================================
PART I — EXACT 21-SLOT QUARTIC COM MAP
========================================================================================
Jacobian + I is identically zero = True
h1 map identity = 0
h7 map identity = 0
h12 map identity = 0
h16 map identity = 0
h19 map identity = 0
h21 map identity = 0

========================================================================================
PART II — EXACT ORDINARY LOCAL 4PN TARGET
========================================================================================
L1 target formula = 0
L7 target formula = 0
L12 target formula = 0
L16 target formula = 0
L19 target formula = 0
L21 target formula = 0
H target G/r degree ceiling - 4 = 0
H target G^2/r^2 degree ceiling - 3 = 0
H target G^3/r^3 degree ceiling - 3 = 0
H target G^4/r^4 degree ceiling - 2 = 0
H target G^5/r^5 degree ceiling - 2 = 0
L target G/r degree ceiling - 4 = 0
L target G^2/r^2 degree ceiling - 4 = 0
L target G^3/r^3 degree ceiling - 4 = 0
L target G^4/r^4 degree ceiling - 4 = 0
L target G^5/r^5 degree ceiling - 2 = 0

========================================================================================
PART III — HAMILTONIAN-CHART GENERIC-FRAME LIFT AND CANONICAL SLICE
========================================================================================
"*)
