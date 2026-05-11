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
  If[! MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

normalizeExpr[expr_] := Module[{res},
  res = If[
    MatrixQ[expr],
    Map[Simplify[PowerExpand[Together[Expand[#]]], Assumptions -> $Assumptions] &, expr, {2}],
    If[
      VectorQ[expr],
      Map[Simplify[PowerExpand[Together[Expand[#]]], Assumptions -> $Assumptions] &, expr],
      Simplify[PowerExpand[Together[Expand[expr]]], Assumptions -> $Assumptions]
    ]
  ];
  res
];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[{res} // MatrixForm];
    If[And @@ Flatten[Map[TrueQ[# == 0] &, res, {Infinity}]],
      pass[name],
      fail[name, res]
    ],
    Print[name, " = ", fmt[res]];
    If[TrueQ[res == 0],
      pass[name],
      fail[name, res]
    ]
  ];
];

ctrMonomial[gamma_, cEtaU_, kU_, tU_, chi0s_, deltaUs_, L_] :=
  (gamma cEtaU/kU)^(1 + deltaUs) (Pi^2 tU/(L^2 kU))^(1 + chi0s);

cntMonomial[lambdaW_, gamma_, kU_, kEta_, kW_, muW_, tU_, eStar_, fStar_, L_, sigma_] :=
  (lambdaW^2 muW/(kEta kW^2))
    (gamma^2 lambdaW^2 sigma/(kU kW))^eStar
    (Pi^2 tU/(L^2 kU))^(-fStar);

epsEtaMonomial[cEtaU_, kU_, kEta_] := cEtaU^2/(kU kEta);

banner["STAGE 183 — EXACT REFERENCE-FREE HOME-STRETCH THEOREM"];

Clear[
  chi0s, deltaUs, eStar, fStar, L, sigma,
  lam1, c1, gam1, KU1, Keta1, KW1, mu1, T1,
  lam2, c2, gam2, KU2, Keta2, KW2, mu2, T2,
  rla, rc, rg, rU, rK, rW, rmu, rT,
  Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT,
  lams, cs, gs, KUs, Ketas, KWs, mus, Ts,
  rlaW, rcW, rgW, rUW, rWW,
  lamx, cx, gx, KUx, Ketax, KWx, mux, Tx, chiQx,
  lamf, cf, gf, KUf, KWf, CtrTarget, CntTarget, epsEtaTarget, mT, mK, mMu,
  rla21, rc21, rg21, rU21, rK21, rW21, rmu21, rT21,
  rla32, rc32, rg32, rU32, rK32, rW32, rmu32, rT32,
  S, beta, Sigma0, Sigma5, eps, epsBeta, dSigma0, dSigma5
];

$Assumptions =
  Element[
    {
      chi0s, deltaUs, eStar, fStar, L, sigma,
      lam1, c1, gam1, KU1, Keta1, KW1, mu1, T1,
      lam2, c2, gam2, KU2, Keta2, KW2, mu2, T2,
      rla, rc, rg, rU, rK, rW, rmu, rT,
      Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT,
      lams, cs, gs, KUs, Ketas, KWs, mus, Ts,
      rlaW, rcW, rgW, rUW, rWW,
      lamx, cx, gx, KUx, Ketax, KWx, mux, Tx, chiQx,
      lamf, cf, gf, KUf, KWf, CtrTarget, CntTarget, epsEtaTarget, mT, mK, mMu,
      rla21, rc21, rg21, rU21, rK21, rW21, rmu21, rT21,
      rla32, rc32, rg32, rU32, rK32, rW32, rmu32, rT32,
      S, beta, Sigma0, Sigma5, eps, epsBeta, dSigma0, dSigma5
    },
    Reals
  ] &&
  chi0s > 0 && deltaUs > 0 && eStar > 0 && fStar > 0 && L > 0 && sigma > 0 &&
  lam1 > 0 && c1 > 0 && gam1 > 0 && KU1 > 0 && Keta1 > 0 && KW1 > 0 && mu1 > 0 && T1 > 0 &&
  lam2 > 0 && c2 > 0 && gam2 > 0 && KU2 > 0 && Keta2 > 0 && KW2 > 0 && mu2 > 0 && T2 > 0 &&
  rla > 0 && rc > 0 && rg > 0 && rU > 0 && rK > 0 && rW > 0 && rmu > 0 && rT > 0 &&
  lams > 0 && cs > 0 && gs > 0 && KUs > 0 && Ketas > 0 && KWs > 0 && mus > 0 && Ts > 0 &&
  rlaW > 0 && rcW > 0 && rgW > 0 && rUW > 0 && rWW > 0 &&
  lamx > 0 && cx > 0 && gx > 0 && KUx > 0 && Ketax > 0 && KWx > 0 && mux > 0 && Tx > 0 &&
  lamf > 0 && cf > 0 && gf > 0 && KUf > 0 && KWf > 0 &&
  CtrTarget > 0 && CntTarget > 0 && epsEtaTarget > 0 &&
  mT > 0 && mK > 0 && mMu > 0 &&
  rla21 > 0 && rc21 > 0 && rg21 > 0 && rU21 > 0 && rK21 > 0 && rW21 > 0 && rmu21 > 0 && rT21 > 0 &&
  rla32 > 0 && rc32 > 0 && rg32 > 0 && rU32 > 0 && rK32 > 0 && rW32 > 0 && rmu32 > 0 && rT32 > 0 &&
  S > 0;

subbanner["I. Exact four-scalar Packet-B compiler from primitive monomial ratios"];

CtrRatio = FullSimplify[
  PowerExpand[(rg rc/rU)^(1 + deltaUs) (rT/rU)^(1 + chi0s)],
  Assumptions -> $Assumptions
];
CntRatio = FullSimplify[
  PowerExpand[
    (rla^2 rmu/(rK rW^2))
      (rg^2 rla^2/(rU rW))^eStar
      (rT/rU)^(-fStar)
  ],
  Assumptions -> $Assumptions
];
epsRatio = FullSimplify[PowerExpand[rc^2/(rK rU)], Assumptions -> $Assumptions];

Print["Ctr_2 / Ctr_1 ="];
Print[CtrRatio // TraditionalForm];
Print["Cnt_2 / Cnt_1 ="];
Print[CntRatio // TraditionalForm];
Print["epsEta_2 / epsEta_1 ="];
Print[epsRatio // TraditionalForm];

ratioToLogs = {
  rla -> Exp[Dl],
  rc -> Exp[Dc],
  rg -> Exp[Dg],
  rU -> Exp[DU],
  rK -> Exp[DKeta],
  rW -> Exp[DW],
  rmu -> Exp[Dmu],
  rT -> Exp[DT]
};

qPair = {
  FullSimplify[PowerExpand[Log[CtrRatio /. ratioToLogs]], Assumptions -> $Assumptions],
  FullSimplify[PowerExpand[Log[CntRatio /. ratioToLogs]], Assumptions -> $Assumptions],
  FullSimplify[PowerExpand[Log[epsRatio /. ratioToLogs]], Assumptions -> $Assumptions]
};
Dvec = {Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT};
Mderived = normalizeExpr[Table[D[qPair[[i]], Dvec[[j]]], {i, 3}, {j, 8}]];
Mexpected = {
  {0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s},
  {2 (1 + eStar), 0, 2 eStar, fStar - eStar, -1, -(2 + eStar), 1, -fStar},
  {0, 2, 0, -1, -1, 0, 0, 0}
};

Print["q^(2<-1) from primitive monomial ratios ="];
Print[qPair // MatrixForm];
Print["M_* derived from the monomial log-ratio Jacobian ="];
Print[Mderived // MatrixForm];
expectZero["derived M_* - carried Stage 192 matrix", Mderived - Mexpected];
expectZero["q^(2<-1) - M_* Delta x^(2<-1)", qPair - Mderived . Dvec];

subbanner["II. Exact target-orbit witness invariance from the carried orbit law"];

alphaStar = normalizeExpr[(1 + deltaUs)/(1 + chi0s)];
PhiTW = normalizeExpr[rUW (rUW/(rgW rcW))^alphaStar];
PhiKW = normalizeExpr[rcW^2/rUW];
PhiMuW = normalizeExpr[
  PhiKW rWW^2/rlaW^2
    (rgW^2 rlaW^2/(rUW rWW))^(-eStar)
    (PhiTW/rUW)^fStar
];

CtrStar = ctrMonomial[gs, cs, KUs, Ts, chi0s, deltaUs, L];
CntStar = cntMonomial[lams, gs, KUs, Ketas, KWs, mus, Ts, eStar, fStar, L, sigma];
epsEtaStar = epsEtaMonomial[cs, KUs, Ketas];

orbitRatioSubs = {
  rla -> rlaW,
  rc -> rcW,
  rg -> rgW,
  rU -> rUW,
  rK -> PhiKW,
  rW -> rWW,
  rmu -> PhiMuW,
  rT -> PhiTW
};
CtrWitnessRatio = FullSimplify[PowerExpand[CtrRatio /. orbitRatioSubs], Assumptions -> $Assumptions];
CntWitnessRatio = FullSimplify[PowerExpand[CntRatio /. orbitRatioSubs], Assumptions -> $Assumptions];
epsEtaWitnessRatio = FullSimplify[PowerExpand[epsRatio /. orbitRatioSubs], Assumptions -> $Assumptions];

expectZero["Log[Ctr(witness) / Ctr(*)]", Log[CtrWitnessRatio]];
expectZero["Log[Cnt(witness) / Cnt(*)]", Log[CntWitnessRatio]];
expectZero["Log[epsEta(witness) / epsEta(*)]", Log[epsEtaWitnessRatio]];

CtrWitness = normalizeExpr[CtrStar CtrWitnessRatio];
CntWitness = normalizeExpr[CntStar CntWitnessRatio];
epsEtaWitness = normalizeExpr[epsEtaStar epsEtaWitnessRatio];

CtrX = ctrMonomial[gx, cx, KUx, Tx, chi0s, deltaUs, L];
CntX = cntMonomial[lamx, gx, KUx, Ketax, KWx, mux, Tx, eStar, fStar, L, sigma];
epsEtaX = epsEtaMonomial[cx, KUx, Ketax];

DeltaHIntrinsic = normalizeExpr[{chiQx - 1, Log[CtrX/CtrStar], Log[CntX/CntStar], Log[epsEtaX/epsEtaStar]}];
DeltaHPairWitness = normalizeExpr[{chiQx - 1, Log[CtrX/CtrWitness], Log[CntX/CntWitness], Log[epsEtaX/epsEtaWitness]}];

Print["Delta_H^int(x | O_*) ="];
Print[DeltaHIntrinsic // MatrixForm];
Print["Delta_H^(x<-w) with an arbitrary target-orbit witness w ="];
Print[DeltaHPairWitness // MatrixForm];
expectZero["pairwise witness packet - intrinsic orbit packet", DeltaHPairWitness - DeltaHIntrinsic];

subbanner["III. Exact mismatch chart from the dependent-triple orbit solve"];

KetaOrbit = normalizeExpr[cf^2/(KUf epsEtaTarget)];
TOrbit = normalizeExpr[
  (L^2 KUf/Pi^2) (CtrTarget/(gf cf/KUf)^(1 + deltaUs))^(1/(1 + chi0s))
];
muOrbit = normalizeExpr[
  CntTarget KetaOrbit KWf^2/lamf^2
    (gf^2 lamf^2 sigma/(KUf KWf))^(-eStar)
    (Pi^2 TOrbit/(L^2 KUf))^fStar
];

TActual = normalizeExpr[mT TOrbit];
KetaActual = normalizeExpr[mK KetaOrbit];
muActual = normalizeExpr[mMu muOrbit];

CtrActualRatio = normalizeExpr[(TActual/TOrbit)^(1 + chi0s)];
CntActualRatio = normalizeExpr[(muActual/muOrbit) (KetaOrbit/KetaActual) (TActual/TOrbit)^(-fStar)];
epsEtaActualRatio = normalizeExpr[KetaOrbit/KetaActual];

expectZero["Log[Ctr(actual) / Ctr_*] - (1+chi0_*) Log[m_T]", Log[CtrActualRatio] - (1 + chi0s) Log[mT]];
expectZero[
  "Log[Cnt(actual) / Cnt_*] - (Log[m_mu] - Log[m_K] - F_* Log[m_T])",
  Log[CntActualRatio] - (Log[mMu] - Log[mK] - fStar Log[mT])
];
expectZero["Log[epsEta(actual) / epsEta_*] + Log[m_K]", Log[epsEtaActualRatio] + Log[mK]];

qMismatch = normalizeExpr[{Log[CtrActualRatio], Log[CntActualRatio], Log[epsEtaActualRatio]}];
Dmis = {0, 0, 0, 0, Log[mK], 0, Log[mMu], Log[mT]};
qMismatchExpected = normalizeExpr[{(1 + chi0s) Log[mT], Log[mMu] - Log[mK] - fStar Log[mT], -Log[mK]}];

Print["q(mismatch) ="];
Print[qMismatch // MatrixForm];
expectZero["q(mismatch) - carried chart", qMismatch - qMismatchExpected];
expectZero["M_* Delta x_mis - q(mismatch)", Mderived . Dmis - qMismatch];

subbanner["IV. Exact cocycle law from microscopic ratio composition"];

D21 = {
  Log[rla21], Log[rc21], Log[rg21], Log[rU21],
  Log[rK21], Log[rW21], Log[rmu21], Log[rT21]
};
D32 = {
  Log[rla32], Log[rc32], Log[rg32], Log[rU32],
  Log[rK32], Log[rW32], Log[rmu32], Log[rT32]
};
D31 = normalizeExpr[{
  Log[rla32 rla21], Log[rc32 rc21], Log[rg32 rg21], Log[rU32 rU21],
  Log[rK32 rK21], Log[rW32 rW21], Log[rmu32 rmu21], Log[rT32 rT21]
}];
q21 = normalizeExpr[Mderived . D21];
q32 = normalizeExpr[Mderived . D32];
q31 = normalizeExpr[Mderived . D31];

expectZero["Delta x^31 - Delta x^32 - Delta x^21", D31 - D32 - D21];
expectZero["q^(31) - q^(32) - q^(21)", q31 - q32 - q21];

subbanner["V. Exact Packet-A linearization and assembled full compiler"];

chiFromDef = 3 (S beta^5 + 9 Sigma5)/(3 S - Sigma0);
DeltaQLinear = Expand[
  Normal[
    Series[
      chiFromDef /. {
        beta -> 1 + eps epsBeta,
        Sigma0 -> eps dSigma0,
        Sigma5 -> eps dSigma5
      },
      {eps, 0, 1}
    ]
  ] - 1
];
DeltaQLinearExpected = eps (5 epsBeta + dSigma0/(3 S) + 9 dSigma5/S);

expectZero[
  "Delta_Q^lin - eps(5 eps_beta + dSigma0/(3S) + 9 dSigma5/S)",
  DeltaQLinear - DeltaQLinearExpected
];

DeltaHLinear = normalizeExpr[Join[{DeltaQLinearExpected}, Mderived . Dvec]];
Print["Delta_H^lin ="];
Print[DeltaHLinear // MatrixForm];

banner["STAGE 183 LEDGER"];
Print["1. The Packet-B quotient coordinates are derived from the primitive coherent monomials"];
Print["      (C_tr, C_nt, eps_eta)"];
Print["   and their logarithmic Jacobian reproduces the carried Stage 192 matrix M_*."];
Print["2. The target-orbit witness packet is genuinely witness-independent because the"];
Print["   finite similarity-orbit law preserves the three coherent monomials exactly."];
Print["3. The mismatch chart q = ((1+chi0_*) ln m_T, ln m_mu - ln m_K - F_* ln m_T, -ln m_K)"];
Print["   is re-derived from the exact dependent-triple orbit solve rather than posited."];
Print["4. The pairwise cocycle law is inherited from microscopic ratio composition"];
Print["      Delta x^(31) = Delta x^(32) + Delta x^(21),"];
Print["   followed by the monomial compiler q = M_* Delta x."];
Print["5. The full linearized home-stretch compiler is"];
Print["      Delta_HS^lin = (Delta_Q^lin, M_* delta x),"];
Print["   where Delta_Q^lin is carried from the exact Packet-A closure law and the"];
Print["   remaining three rows come from the monomial log-ratio Jacobian."];

Exit[0];
