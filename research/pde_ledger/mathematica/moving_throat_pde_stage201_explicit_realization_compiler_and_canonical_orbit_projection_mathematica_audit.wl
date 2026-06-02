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
  PowerExpand[Together[Expand[stripCE[expr]]]],
  Assumptions -> $Assumptions
];

cleanTensor[expr_] := If[
  ListQ[expr],
  Map[cleanScalar, expr, {ArrayDepth[expr]}],
  cleanScalar[expr]
];

zeroTensorQ[expr_] := If[
  ListQ[expr],
  And @@ (TrueQ[# == 0] & /@ Flatten[expr]),
  TrueQ[expr == 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = cleanTensor[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[zeroTensorQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[zeroTensorQ[res], pass[name], fail[name, res]]
  ];
];

expectTrue[name_String, claim_] := Module[{res},
  res = FullSimplify[stripCE[claim], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 201 -- EXPLICIT REALIZATION COMPILER AND CANONICAL ORBIT PROJECTION"];

Clear[
  chi0s, deltaUs, Estar, Fstar,
  Rtr, Rnt, Reta, r,
  chiQ, CtrStar, CntStar, epsEtaStar,
  Ctr1, Cnt1, eps1, Ctr2, Cnt2, eps2,
  lambda, c, gamma, KU, Keta, KW, mu, Tstate,
  dT, dKeta, dmu,
  Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT
];

$Assumptions =
  Element[
    {
      chi0s, deltaUs, Estar, Fstar,
      Rtr, Rnt, Reta, r,
      chiQ, CtrStar, CntStar, epsEtaStar,
      Ctr1, Cnt1, eps1, Ctr2, Cnt2, eps2,
      lambda, c, gamma, KU, Keta, KW, mu, Tstate,
      dT, dKeta, dmu,
      Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT
    },
    Reals
  ] &&
  chi0s > 0 && deltaUs > 0 && Estar > 0 && Fstar > 0 &&
  Rtr > 0 && Rnt > 0 && Reta > 0 && r > 0 && r != 1 &&
  CtrStar > 0 && CntStar > 0 && epsEtaStar > 0 &&
  Ctr1 > 0 && Cnt1 > 0 && eps1 > 0 &&
  Ctr2 > 0 && Cnt2 > 0 && eps2 > 0 &&
  lambda > 0 && c > 0 && gamma > 0 && KU > 0 &&
  Keta > 0 && KW > 0 && mu > 0 && Tstate > 0;

driftBasis = {
  "Delta_lambda", "Delta_c", "Delta_gamma", "Delta_U",
  "Delta_Keta", "Delta_W", "Delta_mu", "Delta_T"
};

trRow = AssociationThread[
  driftBasis,
  {0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s}
];
ntRow = AssociationThread[
  driftBasis,
  {2 (1 + Estar), 0, 2 Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar}
];
etaRow = AssociationThread[
  driftBasis,
  {0, 2, 0, -1, -1, 0, 0, 0}
];

Mstar = ({trRow, ntRow, etaRow}[[All, #]] & /@ driftBasis) // Transpose;
q = {Log[Rtr], Log[Rnt], Log[Reta]};
depSlots = {8, 5, 7};
freeSlots = {1, 2, 3, 4, 6};
freeSelector = UnitVector[8, #] & /@ freeSlots;
constrainedSystem = Join[Mstar, freeSelector];

Print["M_* ="];
Print[MatrixForm[Mstar]];

subbanner["M1. Right-section identity and repair cancellation"];

Ssection = cleanTensor[
  Transpose[
    Table[
      LinearSolve[
        constrainedSystem,
        Join[UnitVector[3, col], ConstantArray[0, Length[freeSlots]]]
      ],
      {col, 3}
    ]
  ]
];
dxRep = cleanTensor[
  LinearSolve[
    constrainedSystem,
    Join[-q, ConstantArray[0, Length[freeSlots]]]
  ]
];

Print["S from constrained LinearSolve ="];
Print[MatrixForm[Ssection]];
Print["Delta x_rep from constrained repair solve ="];
Print[prettyArray[dxRep]];
expectZero["M1 residual: M_* . S - I_3", Mstar . Ssection - IdentityMatrix[3]];
expectZero["M1 residual: M_* . Delta x_rep + q", Mstar . dxRep + q];

subbanner["M2. Mismatch chart by exponent arithmetic"];

mT = Exp[q[[1]]/(1 + chi0s)];
mK = Exp[-q[[3]]];
mMu = Exp[q[[2]] - q[[3]] + Fstar q[[1]]/(1 + chi0s)];

expectZero[
  "M2 residual: m_T - R_tr^(1/(1+chi0_star))",
  mT - Rtr^(1/(1 + chi0s))
];
expectZero["M2 residual: m_K - R_eta^(-1)", mK - Reta^(-1)];
expectZero[
  "M2 residual: m_mu - R_nt R_eta^(-1) R_tr^(F_star/(1+chi0_star))",
  mMu - Rnt Reta^(-1) Rtr^(Fstar/(1 + chi0s))
];

subbanner["M3. Repair vector closed form and support"];

dxRepExpected = {
  0,
  0,
  0,
  0,
  q[[3]],
  0,
  -q[[2]] + q[[3]] - Fstar q[[1]]/(1 + chi0s),
  -q[[1]]/(1 + chi0s)
};

expectZero["M3 residual: Delta x_rep - closed form", dxRep - dxRepExpected];
expectZero["M3 residual: free-coordinate support rows vanish", dxRep[[freeSlots]]];
expectZero[
  "M3 residual: dependent entries match (Keta, mu, T)",
  dxRep[[{5, 7, 8}]] - dxRepExpected[[{5, 7, 8}]]
];

subbanner["M4. Repair vanishes exactly on the target ratios"];

targetRatioSubs = {Rtr -> 1, Rnt -> 1, Reta -> 1};
etaOnlyPerturbation = {Rtr -> 1, Rnt -> 1, Reta -> r};

expectZero[
  "M4 residual: Delta x_rep at R_tr=R_nt=R_eta=1",
  dxRep /. targetRatioSubs
];
expectTrue[
  "M4 check: eta-only perturbation has nonzero Keta repair",
  Exp[dxRep[[5]] /. etaOnlyPerturbation] != 1
];

subbanner["M5. Canonical orbit projection"];

xState = {lambda, c, gamma, KU, Keta, KW, mu, Tstate};
xProj = {
  lambda,
  c,
  gamma,
  KU,
  Keta Reta,
  KW,
  mu Rnt^(-1) Reta Rtr^(-Fstar/(1 + chi0s)),
  Tstate Rtr^(-1/(1 + chi0s))
};
projectionLogDrift = Table[Log[xProj[[i]]/xState[[i]]], {i, 8}];

Print["Pi_O^can(x) ="];
Print[prettyArray[xProj]];
expectZero[
  "M5 residual: Log[x_proj[[i]]/x[[i]]] - Delta x_rep[[i]]",
  projectionLogDrift - dxRep
];
expectZero[
  "M5 residual: fixed-point reduction x_proj - x at target ratios",
  (xProj - xState) /. targetRatioSubs
];

subbanner["M6. Same-free-quintuple uniqueness by direct dependent solve"];

dxDep = {0, 0, 0, 0, dKeta, 0, dmu, dT};
dependentBlock = Mstar[[All, depSlots]];
depDet = Det[dependentBlock];
depSolutions = Solve[
  Thread[Mstar . dxDep == -q],
  {dT, dKeta, dmu},
  Reals
];
depSolution = First[depSolutions];

Print["dependent block columns (T, Keta, mu) ="];
Print[MatrixForm[dependentBlock]];
Print["dependent solution = ", fmt[depSolution]];
expectZero["M6 residual: Det[dependent block] - (1+chi0_star)", depDet - (1 + chi0s)];
expectTrue["M6 check: dependent block determinant is nonzero", depDet != 0];
expectTrue["M6 check: Solve returns a unique dependent triple", Length[depSolutions] == 1];
expectZero[
  "M6 residual: dT solution + q_tr/(1+chi0_star)",
  (dT /. depSolution) + q[[1]]/(1 + chi0s)
];
expectZero["M6 residual: dKeta solution - q_eta", (dKeta /. depSolution) - q[[3]]];
expectZero[
  "M6 residual: dmu solution - closed form",
  (dmu /. depSolution) - (-q[[2]] + q[[3]] - Fstar q[[1]]/(1 + chi0s))
];
expectZero[
  "M6 residual: direct dependent solve satisfies M_* . Delta x_dep + q",
  (Mstar . dxDep + q) /. depSolution
];

subbanner["M7. Intrinsic packet equals pairwise-witness packet"];

pairwisePacket = {
  chiQ - 1,
  Log[Ctr2/Ctr1],
  Log[Cnt2/Cnt1],
  Log[eps2/eps1]
};
intrinsicPacket = {
  chiQ - 1,
  Log[Ctr2/CtrStar],
  Log[Cnt2/CntStar],
  Log[eps2/epsEtaStar]
};
witnessSubs = {Ctr1 -> CtrStar, Cnt1 -> CntStar, eps1 -> epsEtaStar};

expectZero[
  "M7 residual: pairwise witness packet - intrinsic packet",
  (pairwisePacket /. witnessSubs) - intrinsicPacket
];

subbanner["M8. First-order linearized compiler"];

generalDrift = {Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT};
qLin = cleanTensor[Mstar . generalDrift];
dxRepLin = cleanTensor[-Ssection . qLin];

Print["q^lin = M_* . Delta x ="];
Print[prettyArray[qLin]];
Print["Delta x_rep^lin = -S . q^lin ="];
Print[prettyArray[dxRepLin]];
expectZero[
  "M8 residual: M_* . (Delta x + Delta x_rep^lin)",
  Mstar . (generalDrift + dxRepLin)
];
expectZero["M8 residual: q^lin - M_* . Delta x", qLin - Mstar . generalDrift];

banner["STAGE 201 MATHEMATICA LEDGER"];
Print["1. The Stage-192 quotient map has a constrained right section on (T, Keta, mu)."];
Print["2. The canonical repair vector cancels the quotient packet and is supported only on the dependent triple."];
Print["3. The mismatch chart and orbit projection follow by exact exponent arithmetic."];
Print["4. Direct Solve gives the unique same-free-quintuple dependent repair."];
Print["5. The linearized compiler projects any first-order drift back into ker(M_*)."];

Exit[0];
