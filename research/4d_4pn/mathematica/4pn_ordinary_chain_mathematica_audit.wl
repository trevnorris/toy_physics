#!/usr/bin/env wolframscript

ClearAll["Global`*"];
passCount = 0;
failCount = 0;

Get[FileNameJoin[{DirectoryName[$InputFileName], "4pn_shared.wl"}]];

data = ordinaryTranslationData4PN[];

targetL = data["targetL"];
naturalSeed = data["naturalSeed"];
alignedSeed = data["alignedSeed"];
hBlocks = data["hamiltonianBlocks"];
lBlocks = data["ordinaryResidualBlocks"];
ordResidual = data["ordinaryResidualSlotsAligned"];
misBlocks = data["misBlocks"];
seedBases = data["seedBases"];
seedCoords = data["seedCoords"];
deltaExprs = data["deltaExprs"];
alignedBlocks = data["alignedBlocks"];
fullBlocks = data["fullBlocks"];

banner["PART I — EXACT ORDINARY SIGN-FLIP AND ALIGNED SEED"];

Do[
  checkZero["sign flip " <> block, lBlocks[block] + hBlocks[block]],
  {block, {"Q", "T", "S", "U", "W"}}
];

Do[
  checkTrue["ordinary residual COM verification " <> block,
    And @@ MapThread[FullSimplify[Expand[#1 - #2]] === 0 &, {blockSlots4PN[toNu4PN[lBlocks[block]], block], ordResidual[block]}]
  ],
  {block, {"Q", "T", "S", "U", "W"}}
];

Do[
  checkZero["aligned free slot " <> ToString[i], alignedSeed[i] - targetL[i]],
  {i, 1, 6}
];

checkZero["Q misalignment degree ceiling - 4", Max[Exponent[Expand[#], nuSym] & /@ misBlocks["Q"]] - 4];
checkZero["T misalignment degree ceiling - 4", Max[Exponent[Expand[#], nuSym] & /@ misBlocks["T"]] - 4];
checkZero["S misalignment degree ceiling - 4", Max[Exponent[Expand[#], nuSym] & /@ misBlocks["S"]] - 4];
checkZero["U misalignment degree ceiling - 4", Max[Exponent[Expand[#], nuSym] & /@ misBlocks["U"]] - 4];
checkZero["W misalignment is zero", First[misBlocks["W"]]];

banner["PART II — STRUCTURED GENERIC-FRAME LIFT OF THE SEED-ALIGNMENT CORRECTION"];

Do[
  checkTrue["structured solve " <> block,
    And @@ Flatten[Expand[
      imageMatrixPolynomial4PN[block, seedBases[block], 4].seedCoords[block] - targetVector4PN[misBlocks[block], 4]
    ] == 0]
  ],
  {block, {"Q", "T", "S", "U"}}
];

Do[
  checkTrue["delta-seed COM verification " <> block,
    And @@ MapThread[FullSimplify[Expand[#1 - #2]] === 0 &, {blockSlots4PN[toNu4PN[deltaExprs[block]], block], misBlocks[block]}]
  ],
  {block, {"Q", "T", "S", "U"}}
];
checkZero["W delta-seed block", deltaExprs["W"]];

banner["PART III — FULL LOCAL REFEREE RECONSTRUCTION"];

Do[
  checkTrue["aligned block COM verification " <> block,
    And @@ MapThread[FullSimplify[Expand[#1 - #2]] === 0 &, {
      blockSlots4PN[toNu4PN[alignedBlocks[block]], block],
      If[block === "W", {alignedSeed[21]}, alignedSeed /@ blockRanges4PN[block]]
    }]
  ],
  {block, {"Q", "T", "S", "U", "W"}}
];

Do[
  checkTrue["full local candidate COM verification " <> block,
    And @@ MapThread[FullSimplify[Expand[#1 - #2]] === 0 &, {
      blockSlots4PN[toNu4PN[fullBlocks[block]], block],
      blockTargetSlots4PN[targetL, block]
    }]
  ],
  {block, {"Q", "T", "S", "U", "W"}}
];

Do[
  checkTrue["target decomposition " <> block,
    And @@ MapThread[FullSimplify[Expand[#1 - #2]] === 0 &, {
      blockTargetSlots4PN[targetL, block],
      If[block === "W",
        {naturalSeed[21] + First[misBlocks["W"]] + First[ordResidual["W"]]},
        (naturalSeed /@ blockRanges4PN[block]) + misBlocks[block] + ordResidual[block]
      ]
    }]
  ],
  {block, {"Q", "T", "S", "U", "W"}}
];

banner["FINAL ORDINARY 4PN CHAIN LEDGER"];
finalizeAudit[];

(*"
Output:

"*)
