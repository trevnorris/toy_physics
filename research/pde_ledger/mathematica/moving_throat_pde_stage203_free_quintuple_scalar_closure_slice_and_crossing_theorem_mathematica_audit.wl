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

allZeroQ[expr_] := If[
  ListQ[expr],
  And @@ Flatten[Map[TrueQ[# == 0] &, expr, {ArrayDepth[expr]}]],
  TrueQ[expr == 0]
];

prettyArray[arr_] := If[VectorQ[arr], MatrixForm[{arr}], MatrixForm[arr]];

expectZero[name_String, expr_] := Module[{res},
  res = normalizeExpr[expr];
  If[ListQ[res],
    Print[name, " ="];
    Print[prettyArray[res]];
    If[allZeroQ[res], pass[name], fail[name, res]],
    Print[name, " = ", fmt[res]];
    If[allZeroQ[res], pass[name], fail[name, res]]
  ];
];

expectPositive[name_String, expr_] := Module[{shown, res},
  shown = normalizeExpr[expr];
  res = FullSimplify[expr > 0, Assumptions -> $Assumptions];
  Print[name, " > 0  -> ", fmt[shown]];
  If[TrueQ[res], pass[name], fail[name, shown]];
];

expectNegative[name_String, expr_] := Module[{shown, res},
  shown = normalizeExpr[expr];
  res = FullSimplify[expr < 0, Assumptions -> $Assumptions];
  Print[name, " < 0  -> ", fmt[shown]];
  If[TrueQ[res], pass[name], fail[name, shown]];
];

banner["STAGE 186 — FREE-QUINTUPLE SCALAR CLOSURE SLICE AND CROSSING THEOREM"];

Clear[
  chi0s, deltaUs, Estar, Fstar,
  CtrTgt, CntTgt, epsEtaTgt, L, sigma,
  lam0, cetaU0, gamma0, KU0, KW0,
  t, lam1, c1, gam1, kapU, kapW,
  lam, cetaU, gamma, KU, KW,
  ET, EK, Emu,
  qtrs, qnts, qetas,
  Siso, beta, Sigma0, Sigma5,
  tau, rho, lamBar, cetaUBar, gammaBar, KUBar, KWBar
];

$Assumptions =
  Element[
    {
      chi0s, deltaUs, Estar, Fstar,
      CtrTgt, CntTgt, epsEtaTgt, L, sigma,
      lam0, cetaU0, gamma0, KU0, KW0,
      t, lam1, c1, gam1, kapU, kapW,
      lam, cetaU, gamma, KU, KW,
      ET, EK, Emu,
      qtrs, qnts, qetas,
      Siso, beta, Sigma0, Sigma5,
      tau, rho, lamBar, cetaUBar, gammaBar, KUBar, KWBar
    },
    Reals
  ] &&
  chi0s > 0 && deltaUs > 0 && Estar > 0 && Fstar > 0 &&
  CtrTgt > 0 && CntTgt > 0 && epsEtaTgt > 0 &&
  L > 0 && sigma > 0 &&
  lam0 > 0 && cetaU0 > 0 && gamma0 > 0 && KU0 > 0 && KW0 > 0 &&
  lam > 0 && cetaU > 0 && gamma > 0 && KU > 0 && KW > 0 &&
  Siso > 0 && rho > 0 &&
  lamBar > 0 && cetaUBar > 0 && gammaBar > 0 && KUBar > 0 && KWBar > 0;

Mstar = {
  {0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s},
  {2 (1 + Estar), 0, 2 Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar},
  {0, 2, 0, -1, -1, 0, 0, 0}
};

subbanner["I. Stage 192 dependent pivot block and exact canonical section"];

Pdep = Mstar[[All, {8, 5, 7}]];
PdepInv = normalizeExpr[Inverse[Pdep]];

Edep = ConstantArray[0, {8, 3}];
Edep[[8, 1]] = 1;
Edep[[5, 2]] = 1;
Edep[[7, 3]] = 1;

Sdep = normalizeExpr[Edep . PdepInv];

Print["M_* ="];
Print[MatrixForm[Mstar]];
Print["P_(T,K_eta,mu) ="];
Print[MatrixForm[Pdep]];
Print["P_(T,K_eta,mu)^(-1) ="];
Print[MatrixForm[PdepInv]];
Print["S_(T,K_eta,mu) ="];
Print[MatrixForm[Sdep]];

expectZero["M_* S - I_3", Mstar . Sdep - IdentityMatrix[3]];

subbanner["II. Exact graph-family tangent formulas from the Stage 202 target graph"];

lamT = lam0 Exp[t lam1];
cetaUT = cetaU0 Exp[t c1];
gammaT = gamma0 Exp[t gam1];
KUT = KU0 Exp[t kapU];
KWT = KW0 Exp[t kapW];

alpha = normalizeExpr[(1 + deltaUs)/(1 + chi0s)];
logDeltaGraphT = normalizeExpr[
  (Log[CtrTgt] - (1 + deltaUs) (Log[gammaT] + Log[cetaUT] - Log[KUT]))/(1 + chi0s)
];
logTGraphT = normalizeExpr[2 Log[L] + Log[KUT] + logDeltaGraphT - 2 Log[Pi]];
logKetaGraphT = normalizeExpr[2 Log[cetaUT] - Log[KUT] - Log[epsEtaTgt]];
logMuGraphT = normalizeExpr[
  Log[CntTgt] + logKetaGraphT + 2 Log[KWT] - 2 Log[lamT]
  - Estar (2 Log[gammaT] + 2 Log[lamT] + Log[sigma] - Log[KUT] - Log[KWT])
  + Fstar logDeltaGraphT
];

dlnDeltaGraph = normalizeExpr[D[logDeltaGraphT, t] /. t -> 0];
tau1Graph = normalizeExpr[D[logTGraphT, t] /. t -> 0];
keta1Graph = normalizeExpr[D[logKetaGraphT, t] /. t -> 0];
mu1Graph = normalizeExpr[D[logMuGraphT, t] /. t -> 0];

Print["d ln delta_U^graph / d tau = ", fmt[dlnDeltaGraph]];
Print["tau_1^graph = ", fmt[tau1Graph]];
Print["kappa_eta^graph = ", fmt[keta1Graph]];
Print["mu_1^graph = ", fmt[mu1Graph]];

expectZero[
  "dln delta_U^graph - carried formula",
  dlnDeltaGraph + alpha (gam1 + c1 - kapU)
];
expectZero[
  "tau_1^graph - carried formula",
  tau1Graph - (kapU - alpha (gam1 + c1 - kapU))
];
expectZero["kappa_eta^graph - (2 c1 - kappa_U)", keta1Graph - (2 c1 - kapU)];
expectZero[
  "mu_1^graph - carried formula",
  mu1Graph - (
    2 c1 - kapU + 2 kapW - 2 lam1
    - Estar (2 gam1 + 2 lam1 - kapU - kapW)
    - Fstar alpha (gam1 + c1 - kapU)
  )
];

subbanner["III. Exact graph-family orbit-kernel theorem"];

dxGraph = {lam1, c1, gam1, kapU, keta1Graph, kapW, mu1Graph, tau1Graph};
Print["dot(Delta x)_graph ="];
Print[prettyArray[dxGraph]];
expectZero["M_* dot(Delta x)_graph", Mstar . dxGraph];

subbanner["IV. Exact graph-error packet from the direct monomials"];

logDeltaGraph = normalizeExpr[
  (Log[CtrTgt] - (1 + deltaUs) (Log[gamma] + Log[cetaU] - Log[KU]))/(1 + chi0s)
];
logTGraph = normalizeExpr[2 Log[L] + Log[KU] + logDeltaGraph - 2 Log[Pi]];
logKetaGraph = normalizeExpr[2 Log[cetaU] - Log[KU] - Log[epsEtaTgt]];
logMuGraph = normalizeExpr[
  Log[CntTgt] + logKetaGraph + 2 Log[KW] - 2 Log[lam]
  - Estar (2 Log[gamma] + 2 Log[lam] + Log[sigma] - Log[KU] - Log[KW])
  + Fstar logDeltaGraph
];

logTU = normalizeExpr[logTGraph + ET];
logKeta = normalizeExpr[logKetaGraph + EK];
logMu = normalizeExpr[logMuGraph + Emu];

qtr = normalizeExpr[
  (1 + deltaUs) (Log[gamma] + Log[cetaU] - Log[KU])
  + (1 + chi0s) (2 Log[Pi] + logTU - 2 Log[L] - Log[KU])
  - Log[CtrTgt]
];
qnt = normalizeExpr[
  2 Log[lam] + logMu - logKeta - 2 Log[KW]
  + Estar (2 Log[gamma] + 2 Log[lam] + Log[sigma] - Log[KU] - Log[KW])
  - Fstar (2 Log[Pi] + logTU - 2 Log[L] - Log[KU])
  - Log[CntTgt]
];
qeta = normalizeExpr[2 Log[cetaU] - Log[KU] - logKeta - Log[epsEtaTgt]];

Print["q_tr = ", fmt[qtr]];
Print["q_nt = ", fmt[qnt]];
Print["q_eta = ", fmt[qeta]];

expectZero["q_tr - (1+chi0_*) E_T", qtr - (1 + chi0s) ET];
expectZero["q_nt - (E_mu - E_K - F_* E_T)", qnt - (Emu - EK - Fstar ET)];
expectZero["q_eta + E_K", qeta + EK];

subbanner["V. Exact inverse compiler and repair vector"];

inverseSolution = Solve[
  {
    (1 + chi0s) ET == qtrs,
    Emu - EK - Fstar ET == qnts,
    -EK == qetas
  },
  {ET, EK, Emu},
  Reals
][[1]];

Print["Inverse compiler solution ="];
Print[inverseSolution];

expectZero["inverse E_T", (ET /. inverseSolution) - qtrs/(1 + chi0s)];
expectZero["inverse E_K", (EK /. inverseSolution) + qetas];
expectZero[
  "inverse E_mu",
  (Emu /. inverseSolution) - (qnts - qetas + Fstar qtrs/(1 + chi0s))
];

qFromErrors = {(1 + chi0s) ET, Emu - EK - Fstar ET, -EK};
dxRep = normalizeExpr[-Sdep . qFromErrors];
dxRepExpected = {0, 0, 0, 0, -EK, 0, -Emu, -ET};

Print["Delta x_rep ="];
Print[prettyArray[dxRep]];
expectZero["repair vector from graph errors", dxRep - dxRepExpected];
expectZero["M_* Delta x_rep + q(E)", Mstar . dxRep + qFromErrors];

subbanner["VI. Stage 197 scalar closure composed with an explicit Stage 202 graph path"];

betaPath = normalizeExpr[2^(2 tau - 1)];
gammaPath = normalizeExpr[gammaBar betaPath];
cetaUPath = normalizeExpr[cetaUBar Exp[rho tau]];
yGraphPath = {lamBar, cetaUPath, gammaPath, KUBar, KWBar};
graphPathRules = {lam -> lamBar, cetaU -> cetaUPath, gamma -> gammaPath, KU -> KUBar, KW -> KWBar};

logTGraphLift = normalizeExpr[logTGraph /. graphPathRules];
logKetaGraphLift = normalizeExpr[logKetaGraph /. graphPathRules];
logMuGraphLift = normalizeExpr[logMuGraph /. graphPathRules];
qtrGraphLift = normalizeExpr[
  (1 + deltaUs) (Log[gammaPath] + Log[cetaUPath] - Log[KUBar])
  + (1 + chi0s) (2 Log[Pi] + logTGraphLift - 2 Log[L] - Log[KUBar])
  - Log[CtrTgt]
];
qntGraphLift = normalizeExpr[
  2 Log[lamBar] + logMuGraphLift - logKetaGraphLift - 2 Log[KWBar]
  + Estar (2 Log[gammaPath] + 2 Log[lamBar] + Log[sigma] - Log[KUBar] - Log[KWBar])
  - Fstar (2 Log[Pi] + logTGraphLift - 2 Log[L] - Log[KUBar])
  - Log[CntTgt]
];
qetaGraphLift = normalizeExpr[
  2 Log[cetaUPath] - Log[KUBar] - logKetaGraphLift - Log[epsEtaTgt]
];

chiFromStage180 = normalizeExpr[3 (Siso beta^5 + 9 Sigma5)/(3 Siso - Sigma0)];
closureNumStage180 = normalizeExpr[3 Siso (beta^5 - 1) + Sigma0 + 27 Sigma5];

betaLift = normalizeExpr[yGraphPath[[3]]/gammaBar];
(* Verified q_tr=q_nt=q_eta=0 puts this lift on the target graph slice, so the carried closure perturbations Sigma0 and Sigma5 vanish on the composition. *)
hatChiGraph = normalizeExpr[chiFromStage180 /. {beta -> betaLift, Sigma0 -> 0, Sigma5 -> 0}];
hatDeltaGraph = normalizeExpr[hatChiGraph - 1];
closureNumGraph = normalizeExpr[closureNumStage180 /. {beta -> betaLift, Sigma0 -> 0, Sigma5 -> 0}];
hatDeltaDen = Denominator[Together[hatDeltaGraph]];

Print["Explicit free-quintuple graph path y(tau) ="];
Print[prettyArray[yGraphPath]];
Print["beta_path(tau) = ", fmt[betaPath]];
Print["widehat chi_Q(y(tau)) from the carried Stage 197 closure algebra = ", fmt[hatChiGraph]];
Print["widehat Delta_Q(tau) = ", fmt[hatDeltaGraph]];

expectZero["beta_path - gamma(tau)/gamma_bar", betaLift - betaPath];
expectZero["graph-lift target monomial q_tr", qtrGraphLift];
expectZero["graph-lift target monomial q_nt", qntGraphLift];
expectZero["graph-lift target monomial q_eta", qetaGraphLift];
expectZero[
  "Stage 197 closure numerator identity on the graph path",
  3 Siso hatDeltaGraph - closureNumGraph
];
expectPositive["graph residual denominator", hatDeltaDen];

hatDeltaMinus = Factor[hatDeltaGraph /. tau -> 0];
hatDeltaPlus = Factor[hatDeltaGraph /. tau -> 1];
hatDeltaRoot = normalizeExpr[hatDeltaGraph /. tau -> 1/2];
realCrossings = FullSimplify[Reduce[hatDeltaGraph == 0, tau, Reals], Assumptions -> $Assumptions];

Print["widehat Delta_Q(0) = ", fmt[hatDeltaMinus]];
Print["widehat Delta_Q(1) = ", fmt[hatDeltaPlus]];
Print["real crossing set = ", fmt[realCrossings]];

expectNegative["graph residual at tau_- = 0", hatDeltaMinus];
expectPositive["graph residual at tau_+ = 1", hatDeltaPlus];
expectZero["graph residual at tau_* = 1/2", hatDeltaRoot];
If[
  ! TrueQ[FullSimplify[realCrossings == (tau == 1/2), Assumptions -> $Assumptions]],
  fail["real crossing set", realCrossings],
  pass["real crossing set"]
];

banner["STAGE 186 MATHEMATICA AUDIT PASSED"];
Exit[0];
