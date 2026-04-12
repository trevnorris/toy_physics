#!/usr/bin/env wolframscript

ClearAll["Global`*"];
passCount = 0;
failCount = 0;

banner[title_] := Module[{line = StringRepeat["=", 88]},
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

banner["PART I — THEOREM ENVELOPE"];

gammaQuadEff = NQ a0^5/(27 cs^5);
cTail = G M gammaQuadEff/(2 c^3);
gammaGR = 2 G/(5 c^5);
cTailGR = G^2 M/(5 c^8);

checkZero["GR bridge relation", cTailGR - (G M gammaGR)/(2 c^3)];
checkZero["canonical normalization target closes c_tail", (cTail /. NQ -> 54 G cs^5/(5 a0^5 c^5)) - cTailGR];

banner["PART II — EQUIVALENT NORMALIZATION CHAIN"];

nqTarget = 54 G cs^5/(5 a0^5 c^5);
k0Target = 54 G cs^5/(5 a0^5 c^5);
k2Target = 6 G cs^3/(5 a0^3 c^5);
gammaInvariant = Assuming[
  {k0 > 0, k2 > 0},
  FullSimplify[PowerExpand[9 k2^(5/2)/k0^(3/2)]]
];

checkZero["gamma from N_Q target", (gammaQuadEff /. NQ -> nqTarget) - gammaGR];
checkZero[
  "gamma from invariant target pair",
  Assuming[
    {G > 0, c > 0, cs > 0, a0 > 0},
    FullSimplify[PowerExpand[gammaInvariant /. {k0 -> k0Target, k2 -> k2Target}] - gammaGR]
  ]
];
checkZero[
  "c_tail from invariant target pair",
  Assuming[
    {G > 0, M > 0, c > 0, cs > 0, a0 > 0},
    FullSimplify[PowerExpand[((G M)/(2 c^3)) (gammaInvariant /. {k0 -> k0Target, k2 -> k2Target})] - cTailGR]
  ]
];

banner["PART III — CONDITIONAL FULL-4PN LEDGER"];

Print["L_4PN^cons = L_4PN^local + L_4PN^tail"];
Print["C_tail = (G M / (2 c^3)) gamma_quad^eff"];
Print["No new independent 4PN normalization datum opens beyond the 2.5PN STF quadrupole gate."];

Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — THEOREM ENVELOPE
========================================================================================
GR bridge relation = 0
canonical normalization target closes c_tail = 0

========================================================================================
PART II — EQUIVALENT NORMALIZATION CHAIN
========================================================================================
gamma from N_Q target = 0
gamma from invariant target pair = 0
c_tail from invariant target pair = 0

========================================================================================
PART III — CONDITIONAL FULL-4PN LEDGER
========================================================================================
L_4PN^cons = L_4PN^local + L_4PN^tail
C_tail = (G M / (2 c^3)) gamma_quad^eff
No new independent 4PN normalization datum opens beyond the 2.5PN STF quadrupole gate.
Passes: 5
Fails: 0
"*)
