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
    Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr, {2}],
    If[
      VectorQ[expr],
      Map[FullSimplify[#, Assumptions -> $Assumptions] &, expr],
      FullSimplify[expr, Assumptions -> $Assumptions]
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

banner["STAGE 222 — RIGID-MOUTH PHYSICAL NORMAL FORM"];

Clear[
  U, V, T2, T2ref, epsEta, epsEtaRef, Rtarget, RtargetRef, Lambda0,
  chi0, deltaU, ZW, OmegaW2, eps, Bcoeff, Ccoeff, Rtr, RtrRef,
  zeta, Mmix, h, dlnT2, dlnepsEta
];

$Assumptions =
  Element[{U, V, zeta, Mmix, h, dlnT2, dlnepsEta}, Reals] &&
  Element[
    {
      T2, T2ref, epsEta, epsEtaRef, Rtarget, RtargetRef, Lambda0,
      chi0, deltaU, ZW, OmegaW2, eps, Rtr, RtrRef
    },
    Reals
  ] &&
  Element[{Bcoeff, Ccoeff}, Reals] &&
  T2 > 0 && T2ref > 0 && 0 < epsEta < 1 && 0 < epsEtaRef < 1 &&
  Rtarget > 0 && RtargetRef > 0 && Lambda0 > 0 &&
  chi0 > 0 && deltaU > 0 && ZW > 0 && OmegaW2 > 0 && eps > 0 &&
  Rtr > 0 && RtrRef > 0 && Bcoeff != 0 && Ccoeff != 0;

subbanner["I. Exact rigid-mouth physical logarithmic chart"];

Udef = Log[T2/T2ref];
Vdef = Log[epsEta/epsEtaRef];

expectZero["physical log coordinate U", Exp[Udef] - T2/T2ref];
expectZero["physical log coordinate V", Exp[Vdef] - epsEta/epsEtaRef];

Mphys = IdentityMatrix[2];
xPhys = {U, V};
qPhys = Mphys . xPhys;
expectZero["diagonal rigid-mouth packet compiler", qPhys - {U, V}];

RtargetFromIdentity = FullSimplify[Lambda0 (1 - epsEta)/T2, Assumptions -> $Assumptions];
RtargetRefFromIdentity = FullSimplify[Lambda0 (1 - epsEtaRef)/T2ref, Assumptions -> $Assumptions];
ratioFromIdentity = FullSimplify[RtargetFromIdentity/RtargetRefFromIdentity, Assumptions -> $Assumptions];
ratioFromUV = FullSimplify[
  ((1 - epsEtaRef Exp[V])/(1 - epsEtaRef)) Exp[-U],
  Assumptions -> $Assumptions
];

expectZero[
  "target-ratio reconstruction in physical chart",
  FullSimplify[
    ratioFromUV /. {U -> Udef, V -> Vdef},
    Assumptions -> $Assumptions
  ] - ratioFromIdentity
];

subbanner["II. Exact physical projectors and commuting finite legs"];

PT = {{1, 0}, {0, 0}};
Peta = {{0, 0}, {0, 1}};

expectZero["P_T idempotence", PT . PT - PT];
expectZero["P_eta idempotence", Peta . Peta - Peta];
expectZero["P_T P_eta = 0", PT . Peta];
expectZero["P_eta P_T = 0", Peta . PT];
expectZero["P_T + P_eta = I", PT + Peta - IdentityMatrix[2]];

xT = PT . xPhys;
xEta = Peta . xPhys;
expectZero["physical packet decomposition", xPhys - xT - xEta];

ratioTransfer = FullSimplify[ratioFromUV /. V -> 0, Assumptions -> $Assumptions];
ratioDressing = FullSimplify[ratioFromUV /. U -> 0, Assumptions -> $Assumptions];

expectZero["pure transfer-shape leg", ratioTransfer - Exp[-U]];
expectZero[
  "pure dressing leg",
  ratioDressing - (1 - epsEtaRef Exp[V])/(1 - epsEtaRef)
];
expectZero[
  "exact commutativity / factorization of transfer and dressing legs",
  ratioFromUV - ratioTransfer ratioDressing
];

subbanner["III. Stage-219/221 dependent-plane compiler in the physical chart"];

T2Stage221 = ZW (1 + chi0)^2/(OmegaW2 (1 - eps)^2);
RtargetStage221 =
  Lambda0 OmegaW2 (1 - epsEta) (1 - eps)^2/(ZW (1 + chi0)^2);
expectZero[
  "Stage-221 transfer-shape identity",
  RtargetStage221 T2Stage221 - Lambda0 (1 - epsEta)
];

qNtGeneral =
  Bcoeff Log[Rtr/RtrRef] +
  Log[(1 - epsEta)/(1 - epsEtaRef)] -
  Log[Rtarget/RtargetRef];
qNtRigidStage221 = FullSimplify[
  ExpandAll @ PowerExpand[
    qNtGeneral /. {
      Rtr -> RtrRef,
      Rtarget -> RtargetStage221,
      RtargetRef -> Lambda0 (1 - epsEtaRef)/T2ref
    }
  ],
  Assumptions -> $Assumptions
];
UStage221 = Log[T2Stage221/T2ref];
VStage221 = Log[epsEta/epsEtaRef];

expectZero[
  "Stage-221 rigid-mouth packet factorization for q_nt",
  Exp[qNtRigidStage221] - T2Stage221/T2ref
];
expectZero[
  "Stage-221 rigid-mouth identification q_nt = U",
  FullSimplify[
    ExpandAll @ PowerExpand[qNtRigidStage221 - UStage221],
    Assumptions -> $Assumptions
  ]
];
expectZero[
  "Stage-221 rigid-mouth identification q_eta = V",
  FullSimplify[
    ExpandAll @ PowerExpand[VStage221 - Vdef],
    Assumptions -> $Assumptions
  ]
];

SrmDep = {{0, 0}, {0, -1}, {1, -1}};
CphysDep = SrmDep . Mphys;
yDep = CphysDep . xPhys;

expectZero[
  "physical-to-microscopic compiler inherited from Stage 219",
  CphysDep - {{0, 0}, {0, -1}, {1, -1}}
];
expectZero[
  "physical-to-microscopic dependent compiler",
  yDep - {0, -V, U - V}
];

xPhysStage221 = {UStage221, VStage221};
yDepStage221 = SrmDep . {qNtRigidStage221, VStage221};
expectZero[
  "Stage-219/221 dependent compiler propagated into the physical chart",
  FullSimplify[
    ExpandAll @ PowerExpand[yDepStage221 - CphysDep . xPhysStage221],
    Assumptions -> $Assumptions
  ]
];

inverseSolution = Solve[
  {
    yK == -qEtaVar,
    yMu == qNtVar - qEtaVar
  },
  {qNtVar, qEtaVar},
  Reals
];
expectedInverse = {{qNtVar -> yMu - yK, qEtaVar -> -yK}};
If[inverseSolution =!= expectedInverse,
  Print["FAIL: unexpected dependent-plane inverse"];
  Print["  actual -> ", fmt[inverseSolution]];
  Exit[1];
];

LphysDep = {{0, -1, 1}, {0, -1, 0}};
expectZero[
  "left inverse reconstructed from Stage-219 dependent-plane equations",
  LphysDep . {0, yK, yMu} - {qNtVar, qEtaVar} /. expectedInverse[[1]]
];
expectZero["left inverse of physical compiler", LphysDep . CphysDep - IdentityMatrix[2]];
expectZero["recovery of U,V from dependent correction", LphysDep . yDep - xPhys];

subbanner["IV. Exact microscopic images of the two physical axes"];

yT = CphysDep . xT;
yEta = CphysDep . xEta;

expectZero["pure transfer-shape microscopic image", yT - {0, 0, U}];
expectZero["pure dressing microscopic image", yEta + V {0, 1, 1}];
expectZero["dependent correction splits into the two physical axes", yDep - yT - yEta];

subbanner["V. Exact correction compilers"];

deltaYStatic = -yT;
deltaYEtaRest = V {0, 1, 1};
deltaYOrbit = -yDep;

expectZero["static-only correction", deltaYStatic - {0, 0, -U}];
expectZero["post-static dressing correction", deltaYEtaRest - {0, V, V}];
expectZero["full orbit-lock correction", deltaYOrbit - {0, V, V - U}];
expectZero["full correction splits into static + dressing parts", deltaYOrbit - deltaYStatic - deltaYEtaRest];
expectZero["full orbit-lock correction cancels the dependent defect", yDep + deltaYOrbit];

subbanner["VI. Propagation of Stage-221 support-blindness"];

T2sb = T2sbFn[zeta, Mmix];
epsEtaSb = epsEtaSbFn[zeta, Mmix];
Ublind = Log[T2sb/T2ref];
Vblind = Log[epsEtaSb/epsEtaRef];
supportBlindRules = {
  Derivative[1, 0][T2sbFn][zeta, Mmix] -> 0,
  Derivative[0, 1][T2sbFn][zeta, Mmix] -> 0,
  Derivative[1, 0][epsEtaSbFn][zeta, Mmix] -> 0,
  Derivative[0, 1][epsEtaSbFn][zeta, Mmix] -> 0
};

expectZero[
  "Stage-221 branch formula for T^2 is support-blind w.r.t. zeta",
  D[T2Stage221, zeta]
];
expectZero[
  "Stage-221 branch formula for T^2 is support-blind w.r.t. M_mix",
  D[T2Stage221, Mmix]
];
expectZero[
  "Stage-221 support-blind T^2 propagates to U w.r.t. zeta",
  (D[Ublind, zeta] /. supportBlindRules)
];
expectZero[
  "Stage-221 support-blind eps_eta propagates to V w.r.t. zeta",
  (D[Vblind, zeta] /. supportBlindRules)
];
expectZero[
  "Stage-221 support-blind T^2 propagates to U w.r.t. M_mix",
  (D[Ublind, Mmix] /. supportBlindRules)
];
expectZero[
  "Stage-221 support-blind eps_eta propagates to V w.r.t. M_mix",
  (D[Vblind, Mmix] /. supportBlindRules)
];

yDepBlind = CphysDep . {Ublind, Vblind};
deltaYStaticBlind = -(CphysDep . {Ublind, 0});
deltaYOrbitBlind = -yDepBlind;

expectZero[
  "Stage-221 support-blindness propagates to the dependent correction w.r.t. zeta",
  (D[#, zeta] & /@ yDepBlind) /. supportBlindRules
];
expectZero[
  "Stage-221 support-blindness propagates to the dependent correction w.r.t. M_mix",
  (D[#, Mmix] & /@ yDepBlind) /. supportBlindRules
];
expectZero[
  "Stage-221 support-blindness propagates to the static correction w.r.t. zeta",
  (D[#, zeta] & /@ deltaYStaticBlind) /. supportBlindRules
];
expectZero[
  "Stage-221 support-blindness propagates to the orbit correction w.r.t. zeta",
  (D[#, zeta] & /@ deltaYOrbitBlind) /. supportBlindRules
];
expectZero[
  "Stage-221 support-blindness propagates to the static correction w.r.t. M_mix",
  (D[#, Mmix] & /@ deltaYStaticBlind) /. supportBlindRules
];
expectZero[
  "Stage-221 support-blindness propagates to the orbit correction w.r.t. M_mix",
  (D[#, Mmix] & /@ deltaYOrbitBlind) /. supportBlindRules
];

subbanner["VII. Cartesian orbit-lock equivalence and first-order form"];

orbitLockSolution = Solve[
  {
    yDep[[2]] == 0,
    yDep[[3]] == 0
  },
  {U, V},
  Reals
];
If[orbitLockSolution =!= {{U -> 0, V -> 0}},
  Print["FAIL: unexpected orbit-lock solution set"];
  Print["  actual -> ", fmt[orbitLockSolution]];
  Exit[1];
];

expectZero["U = V = 0 cancels the dependent defect", yDep /. {U -> 0, V -> 0}];
expectZero["T^2 = T_ref^2 implies U = 0", Udef /. T2 -> T2ref];
expectZero["eps_eta = eps_eta_ref implies V = 0", Vdef /. epsEta -> epsEtaRef];

transferEquilibriumSolution = Solve[T2/T2ref == 1, T2, Reals];
dressingEquilibriumSolution = Solve[epsEta/epsEtaRef == 1, epsEta, Reals];
If[transferEquilibriumSolution =!= {{T2 -> T2ref}},
  Print["FAIL: unexpected transfer-equilibrium solve"];
  Print["  actual -> ", fmt[transferEquilibriumSolution]];
  Exit[1];
];
If[dressingEquilibriumSolution =!= {{epsEta -> epsEtaRef}},
  Print["FAIL: unexpected dressing-equilibrium solve"];
  Print["  actual -> ", fmt[dressingEquilibriumSolution]];
  Exit[1];
];

expectZero["U = V = 0 restores the target ratio", (ratioFromUV /. {U -> 0, V -> 0}) - 1];

T2pert = T2 Exp[h dlnT2];
epsPert = epsEta Exp[h dlnepsEta];

Ufirst = FullSimplify[D[Log[T2pert/T2], h] /. h -> 0, Assumptions -> $Assumptions];
Vfirst = FullSimplify[D[Log[epsPert/epsEta], h] /. h -> 0, Assumptions -> $Assumptions];

expectZero["first-order U compiler", Ufirst - dlnT2];
expectZero["first-order V compiler", Vfirst - dlnepsEta];

yFirst = CphysDep . {Ufirst, Vfirst};
expectZero[
  "first-order dependent correction",
  yFirst - {0, -dlnepsEta, dlnT2 - dlnepsEta}
];

yStaticBlind = yDep /. U -> 0;
expectZero["static-blind line maps to equal-drift dressing ray", yStaticBlind + V {0, 1, 1}];

Print[""];
Print["All Stage-222 symbolic checks passed."];
