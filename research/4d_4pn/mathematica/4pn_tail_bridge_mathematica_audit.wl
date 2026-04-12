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

stfOuter[a_, b_] := Module[{mat},
  mat = (Outer[Times, a, b] + Outer[Times, b, a])/2;
  FullSimplify[mat - Tr[mat] IdentityMatrix[3]/3]
];

banner["PART I — EXACT GR 4PN TAIL COEFFICIENT"];

gammaGR = 2 G/(5 c^5);
alphaTailGR = G^2 M/(5 c^8);

checkZero["alpha_tail_GR - (G M / (2 c^3)) gamma_GR", alphaTailGR - (G M gammaGR)/(2 c^3)];

banner["PART II — NEWTONIAN ORDER-REDUCED STF QUADRUPOLE KERNEL"];

x = {r, 0, 0};
v = {rd, vt, 0};
a = -G M x/r^3;
j = -G M (v/r^3 - 3 rd x/r^4);

i3Direct = FullSimplify[2 mu (3 stfOuter[a, v] + stfOuter[x, j])];
A = stfOuter[x, v];
B = stfOuter[x, x];
i3Formula = FullSimplify[-2 G M mu/r^3 (4 A - 3 rd B/r)];

Do[
  checkZero["I3 direct - formula (" <> ToString[i] <> "," <> ToString[jj] <> ")", i3Direct[[i, jj]] - i3Formula[[i, jj]]],
  {i, 1, 3}, {jj, 1, 3}
];

i3Sq = FullSimplify[Expand[Total[Flatten[i3Formula^2]]]];
i3SqExpected = (8 G^2 M^2 mu^2 (12 v2 - 11 rd^2))/(3 r^4);
checkZero["I3_sq - expected", Expand[i3Sq /. vt^2 -> v2 - rd^2] - i3SqExpected];

banner["PART III — UNIVERSAL LOCAL LOGARITHMIC TAIL SHADOW"];

fTail = FullSimplify[(2 G^2 M i3SqExpected)/(5 c^8)];
fTailPerMu = FullSimplify[((fTail/mu) /. mu -> nu M) /. (G^4 M^4/r^4) -> U^4];
checkZero["tail shadow law", fTailPerMu - (16 nu U^4 (12 v2 - 11 rd^2))/(15 c^8)];

banner["PART IV — TOY-MODEL HEREDITARY BRIDGE"];

gammaQuadEff = NQ a0^5/(27 cs^5);
cTailToy = FullSimplify[(G M gammaQuadEff)/(2 c^3)];
nqTarget = 54 G cs^5/(5 a0^5 c^5);
checkZero["GR recovery from N_Q target", (cTailToy /. NQ -> nqTarget) - alphaTailGR];

k0Target = 54 G cs^5/(5 a0^5 c^5);
k2Target = 6 G cs^3/(5 a0^3 c^5);
gammaBranch = Assuming[
  {k0 > 0, k2 > 0},
  FullSimplify[PowerExpand[9 k2^(5/2)/k0^(3/2)]]
];
checkZero[
  "branch gamma target",
  Assuming[
    {G > 0, c > 0, cs > 0, a0 > 0},
    FullSimplify[PowerExpand[gammaBranch /. {k0 -> k0Target, k2 -> k2Target}] - 2 G/(5 c^5)]
  ]
];
checkZero[
  "tail bridge from invariant pair",
  Assuming[
    {G > 0, M > 0, c > 0, cs > 0, a0 > 0},
    FullSimplify[
      PowerExpand[((G M gammaBranch)/(2 c^3)) /. {k0 -> k0Target, k2 -> k2Target}] - alphaTailGR
    ]
  ]
];

banner["FINAL 4PN TAIL / HEREDITARY BRIDGE LEDGER"];
Print["Passes: ", passCount];
Print["Fails: ", failCount];
If[failCount != 0, Exit[1]];

(*"
Output:

========================================================================================
PART I — EXACT GR 4PN TAIL COEFFICIENT
========================================================================================
alpha_tail_GR - (G M / (2 c^3)) gamma_GR = 0

========================================================================================
PART II — NEWTONIAN ORDER-REDUCED STF QUADRUPOLE KERNEL
========================================================================================
I3 direct - formula (1,1) = 0
I3 direct - formula (1,2) = 0
I3 direct - formula (1,3) = 0
I3 direct - formula (2,1) = 0
I3 direct - formula (2,2) = 0
I3 direct - formula (2,3) = 0
I3 direct - formula (3,1) = 0
I3 direct - formula (3,2) = 0
I3 direct - formula (3,3) = 0
I3_sq - expected = 0

========================================================================================
PART III — UNIVERSAL LOCAL LOGARITHMIC TAIL SHADOW
========================================================================================
tail shadow law = 0

========================================================================================
PART IV — TOY-MODEL HEREDITARY BRIDGE
========================================================================================
GR recovery from N_Q target = 0
branch gamma target = 0
tail bridge from invariant pair = 0

========================================================================================
FINAL 4PN TAIL / HEREDITARY BRIDGE LEDGER
========================================================================================
Passes: 15
Fails: 0
"*)
