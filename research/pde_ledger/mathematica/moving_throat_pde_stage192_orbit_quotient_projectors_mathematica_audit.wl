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

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

stripCE[expr_] := expr /. ConditionalExpression[e_, _] :> e;

cleanScalar[expr_] := FullSimplify[
  Together[Expand[stripCE[expr]]],
  Assumptions -> $Assumptions
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {Length[Dimensions[expr]]}],
  cleanScalar[expr]
];

zeroTensorQ[expr_] := If[
  ListQ[expr],
  And @@ (TrueQ[# === 0] & /@ Flatten[expr]),
  TrueQ[expr === 0]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanTensor[expr];
  Print[name, " = ", fmt[res]];
  If[zeroTensorQ[res], pass[name], fail[name, res]];
];

expectTrue[name_String, claim_] := Module[{res},
  res = FullSimplify[stripCE[claim], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 192 -- MATHEMATICA ORBIT/QUOTIENT PROJECTOR AUDIT"];

Clear[chi, delU, ee, ff, dl, dc, dg, du, dk, dw, dm, dt, qT, qN, qE];

$Assumptions =
  Element[{chi, delU, ee, ff, dl, dc, dg, du, dk, dw, dm, dt, qT, qN, qE}, Reals] &&
  chi > 0 && delU > 0 && ee > 0 && ff > 0;

drift = {dl, dc, dg, du, dk, dw, dm, dt};

(* Ordered drift basis: lambda, c, gamma, U, K_eta, W, mu, T. *)
quotientRows = {
  {0, 1 + delU, 1 + delU, -(2 + chi + delU), 0, 0, 0, 1 + chi},
  {2*(1 + ee), 0, 2*ee, ff - ee, -1, -(2 + ee), 1, -ff},
  {0, 2, 0, -1, -1, 0, 0, 0}
};

depSlots = {8, 5, 7};       (* T, K_eta, mu *)
freeSlots = {1, 2, 3, 4, 6};
freeSelector = UnitVector[8, #] & /@ freeSlots;
constrainedSystem = Join[quotientRows, freeSelector];

subbanner["I. Pivot block and native constrained section"];

pivotBlock = quotientRows[[All, depSlots]];
pivotExpected = {
  {1 + chi, 0, 0},
  {-ff, -1, 1},
  {0, -1, 0}
};
pivotInverse = cleanTensor[Inverse[pivotBlock]];
pivotInverseExpected = {
  {1/(1 + chi), 0, 0},
  {0, 0, -1},
  {ff/(1 + chi), 1, -1}
};

Print["P_(T,K_eta,mu) = ", fmt[pivotBlock]];
expectZero["pivot block - SymPy target block", pivotBlock - pivotExpected];
Print["det(P_(T,K_eta,mu)) = ", fmt[Factor[Det[pivotBlock]]]];
expectZero["det(P) - (1+chi)", Factor[Det[pivotBlock]] - (1 + chi)];
expectZero["pivot inverse - SymPy target inverse", pivotInverse - pivotInverseExpected];
expectZero["P^(-1) P - I_3", pivotInverse . pivotBlock - IdentityMatrix[3]];
expectZero["P P^(-1) - I_3", pivotBlock . pivotInverse - IdentityMatrix[3]];

(* Independent route: each section column is the unique vector with zero free
   coordinates whose quotient packet is the corresponding unit basis vector. *)
sectionFromSolve = Transpose[
  Table[
    LinearSolve[
      constrainedSystem,
      Join[UnitVector[3, col], ConstantArray[0, Length[freeSlots]]]
    ],
    {col, 3}
  ]
];
sectionFromSolve = cleanTensor[sectionFromSolve];

sectionExpected = {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, -1},
  {0, 0, 0},
  {ff/(1 + chi), 1, -1},
  {1/(1 + chi), 0, 0}
};

Print["S from constrained LinearSolve = ", fmt[sectionFromSolve]];
expectZero["section - SymPy target section", sectionFromSolve - sectionExpected];
expectZero["free rows of S vanish", sectionFromSolve[[freeSlots, All]]];
expectZero["M S - I_3", quotientRows . sectionFromSolve - IdentityMatrix[3]];

subbanner["II. Complementary projectors"];

quotientProjector = cleanTensor[sectionFromSolve . quotientRows];
orbitProjector = cleanTensor[IdentityMatrix[8] - quotientProjector];

expectZero["Q^2 - Q", quotientProjector . quotientProjector - quotientProjector];
expectZero["O^2 - O", orbitProjector . orbitProjector - orbitProjector];
expectZero["Q O", quotientProjector . orbitProjector];
expectZero["O Q", orbitProjector . quotientProjector];
expectZero["M O", quotientRows . orbitProjector];
expectZero["M Q - M", quotientRows . quotientProjector - quotientRows];
expectZero["Q free-coordinate rows", quotientProjector[[freeSlots, All]]];

subbanner["III. Quotient failure from a constrained solve"];

packet = cleanTensor[quotientRows . drift];
failureByProjection = cleanTensor[quotientProjector . drift];
failureByConstrainedSolve = cleanTensor[
  LinearSolve[
    constrainedSystem,
    Join[packet, ConstantArray[0, Length[freeSlots]]]
  ]
];
failureExpected = {
  0,
  0,
  0,
  0,
  -packet[[3]],
  0,
  ff*packet[[1]]/(1 + chi) + packet[[2]] - packet[[3]],
  packet[[1]]/(1 + chi)
};

Print["q = M Delta x = ", fmt[packet]];
Print["Delta x_fail = ", fmt[failureByProjection]];
expectZero["Q Delta x - constrained failure solve", failureByProjection - failureByConstrainedSolve];
expectZero["Q Delta x - SymPy dependent-triple support", failureByProjection - failureExpected];
expectZero["M Delta x_fail - q", quotientRows . failureByProjection - packet];

subbanner["IV. Orbit representative from the kernel constraints"];

alpha = (1 + delU)/(1 + chi);
orbitByProjection = cleanTensor[orbitProjector . drift];
orbitByKernelSolve = cleanTensor[
  LinearSolve[
    constrainedSystem,
    Join[ConstantArray[0, 3], drift[[freeSlots]]]
  ]
];
orbitExpected = {
  dl,
  dc,
  dg,
  du,
  2*dc - du,
  dw,
  2*dc - du + 2*dw - 2*dl
    - ee*(2*dg + 2*dl - du - dw)
    - ff*alpha*(dg + dc - du),
  du - alpha*(dg + dc - du)
};

Print["Delta x_orbit = ", fmt[orbitByProjection]];
expectZero["O Delta x - constrained kernel solve", orbitByProjection - orbitByKernelSolve];
expectZero["O Delta x - SymPy orbit law", orbitByProjection - orbitExpected];
expectZero["M Delta x_orbit", quotientRows . orbitByProjection];
expectZero["Delta x - (orbit + fail)", drift - orbitByProjection - failureByProjection];

subbanner["V. Orbit-lock equivalence"];

independentPacket = {qT, qN, qE};
quotientRepresentative = cleanTensor[sectionFromSolve . independentPacket];
representativeExpected = {
  0,
  0,
  0,
  0,
  -qE,
  0,
  ff*qT/(1 + chi) + qN - qE,
  qT/(1 + chi)
};

Print["Canonical representative S q = ", fmt[quotientRepresentative]];
expectZero["S q - SymPy representative", quotientRepresentative - representativeExpected];
expectZero["M (S q) - q", quotientRows . quotientRepresentative - independentPacket];
expectZero["Q (S q) - S q", quotientProjector . quotientRepresentative - quotientRepresentative];
expectZero["O (S q)", orbitProjector . quotientRepresentative];
expectZero["S * 0", sectionFromSolve . {0, 0, 0}];
expectZero["orbit-law branch lies in ker(M)", quotientRows . orbitExpected];
expectZero["Q Delta x - S (M Delta x)", failureByProjection - sectionFromSolve . packet];
expectZero["M (Q Delta x) - M Delta x", quotientRows . failureByProjection - quotientRows . drift];

lockByMap = First[Solve[
  Thread[quotientRows . drift == ConstantArray[0, 3]],
  drift[[depSlots]],
  Reals
]];
lockByProjector = First[Solve[
  Thread[quotientProjector . drift == ConstantArray[0, 8]],
  drift[[depSlots]],
  Reals
]];
expectZero[
  "Solve laws for Q Delta x == 0 and M Delta x == 0 agree",
  (drift[[depSlots]] /. lockByProjector) - (drift[[depSlots]] /. lockByMap)
];
expectZero[
  "Q Delta x vanishes on the M-kernel law",
  quotientProjector . (drift /. lockByMap)
];
expectZero[
  "M Delta x vanishes on the Q-kernel law",
  quotientRows . (drift /. lockByProjector)
];

banner["STAGE 192 MATHEMATICA LEDGER"];
Print["1. The dependent triple (Delta_T, Delta_Keta, Delta_mu) has determinant 1+chi > 0."];
Print["2. A constrained LinearSolve constructs the canonical quotient section S independently."];
Print["3. The native projectors Q = S M and O = I - Q are complementary and exact."];
Print["4. Q Delta x is supported only on (Delta_T, Delta_Keta, Delta_mu)."];
Print["5. O Delta x is the unique kernel representative with unchanged free coordinates."];
Print["6. Solve confirms Q Delta x = 0 iff M Delta x = 0."];

Exit[0];
