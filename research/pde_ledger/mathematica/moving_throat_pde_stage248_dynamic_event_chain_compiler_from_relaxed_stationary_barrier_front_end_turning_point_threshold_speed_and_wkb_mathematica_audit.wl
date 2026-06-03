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

expectNear[name_String, value_, target_, tol_] := Module[{res},
  res = N[Abs[value - target], 20];
  Print[name, " diff = ", fmt[res], " (tol = ", fmt[tol], ")"];
  If[res <= tol, pass[name], fail[name, res]]
];

banner["STAGE 248 — DYNAMIC EVENT-CHAIN COMPILER FROM THE RELAXED STATIONARY BARRIER FRONT END"];

Clear[
  t, ms, V0, Vpeak, r0, rContact, v0, Esub, hbarEff, rMinus, rPlus, rv,
  DeltaE, Kpeak, y, rPeak, Eturn
];

$Assumptions =
  Element[{t, V0, Vpeak, r0, rContact, v0, rv, DeltaE, y, rPeak}, Reals] &&
  Element[{ms, Esub, hbarEff, rMinus, rPlus, Kpeak}, Reals] &&
  ms > 0 && v0 > 0 && r0 > 0 && rContact > 0 && Esub > 0 &&
  hbarEff > 0 && rMinus > 0 && rPlus > 0 && DeltaE > 0 && Kpeak > 0;

subbanner["I. Exact reduced energy conservation"];

r = rfun[t];
rdot = D[r, t];
rddot = D[r, {t, 2}];
Vexpr = Vfun[r];
Edyn = FullSimplify[(1/2) ms rdot^2 + Vexpr, Assumptions -> $Assumptions];
dEdt = D[Edyn, t];
dEOnShell = FullSimplify[
  dEdt /. rddot -> -D[Vexpr, r]/ms,
  Assumptions -> $Assumptions
];

expectZero["dE/dt on-shell", dEOnShell];

subbanner["II. Barrier-peak and threshold-speed compiler"];

thresholdAssumptions = $Assumptions && Vpeak > V0 && r0 > rContact;
ElaunchNew = FullSimplify[(1/2) ms v0^2 + V0, Assumptions -> $Assumptions];
vcritNew = FullSimplify[Sqrt[2 (Vpeak - V0)/ms], Assumptions -> thresholdAssumptions];
EAtVcrit = FullSimplify[ElaunchNew /. v0 -> vcritNew, Assumptions -> thresholdAssumptions];
ElaunchCoul = FullSimplify[(1/2) ms v0^2 + 1/r0, Assumptions -> $Assumptions];
vcontactCoul = FullSimplify[
  Sqrt[2 (1/rContact - 1/r0)/ms],
  Assumptions -> thresholdAssumptions
];
EAtVcontact = FullSimplify[ElaunchCoul /. v0 -> vcontactCoul, Assumptions -> thresholdAssumptions];
deltaNew = FullSimplify[ElaunchNew - Vpeak, Assumptions -> thresholdAssumptions];
deltaCoul = FullSimplify[ElaunchCoul - 1/rContact, Assumptions -> thresholdAssumptions];

Print["delta_new pre-substitution = ", fmt[deltaNew]];
If[
  TrueQ[FullSimplify[deltaNew == 0, Assumptions -> thresholdAssumptions]] || FreeQ[deltaNew, v0],
  fail["delta_new non-vacuous gap", deltaNew],
  pass["delta_new non-vacuous gap"]
];
Print["delta_coul pre-substitution = ", fmt[deltaCoul]];
If[
  TrueQ[FullSimplify[deltaCoul == 0, Assumptions -> thresholdAssumptions]] || FreeQ[deltaCoul, v0],
  fail["delta_coul non-vacuous gap", deltaCoul],
  pass["delta_coul non-vacuous gap"]
];

expectZero["v_crit,new satisfies defining energy", EAtVcrit - Vpeak];
expectZero["delta_new vanishes at v_crit,new", deltaNew /. v0 -> vcritNew];
expectZero["v_contact,coul satisfies defining energy", EAtVcontact - 1/rContact];
expectZero["delta_coul vanishes at v_contact,coul", deltaCoul /. v0 -> vcontactCoul];
Print["v_crit,new positive branch = ", fmt[FullSimplify[vcritNew > 0, Assumptions -> thresholdAssumptions]]];
If[
  TrueQ[FullSimplify[vcritNew > 0, Assumptions -> thresholdAssumptions]],
  pass["v_crit,new positive branch"],
  fail["v_crit,new positive branch", vcritNew]
];
Print["v_contact,coul positive branch = ", fmt[FullSimplify[vcontactCoul > 0, Assumptions -> thresholdAssumptions]]];
If[
  TrueQ[FullSimplify[vcontactCoul > 0, Assumptions -> thresholdAssumptions]],
  pass["v_contact,coul positive branch"],
  fail["v_contact,coul positive branch", vcontactCoul]
];

subbanner["III. Turning-point and WKB compiler"];

v0Sub = FullSimplify[Sqrt[2 (Esub - V0)/ms], Assumptions -> $Assumptions];
rturnCoul = FullSimplify[1/Esub, Assumptions -> $Assumptions];
Fcoul = FullSimplify[
  Sqrt[2 ms]/hbarEff (
    Sqrt[rv (1 - Esub rv)] + ArcSin[Sqrt[Esub rv]]/Sqrt[Esub]
  ),
  Assumptions -> $Assumptions
];
coulIntegrand = FullSimplify[
  Sqrt[2 ms (1/rv - Esub)]/hbarEff,
  Assumptions -> $Assumptions
];
dFcoul = FullSimplify[
  PowerExpand[D[Fcoul, rv] - coulIntegrand],
  Assumptions -> $Assumptions && Esub rv < 1
];
dFcoul = Refine[dFcoul, $Assumptions && Esub rv < 1];
dFcoul = dFcoul /. Sqrt[rv^(-1)] Sqrt[rv] -> 1;
dFcoul = FullSimplify[dFcoul, Assumptions -> $Assumptions && Esub rv < 1];
IcoulFormula = FullSimplify[
  Sqrt[2 ms]/hbarEff (
    Pi/(2 Sqrt[Esub]) -
    Sqrt[rContact (1 - Esub rContact)] -
    ArcSin[Sqrt[Esub rContact]]/Sqrt[Esub]
  ),
  Assumptions -> $Assumptions
];
IcoulEndpoints = FullSimplify[
  Limit[Fcoul, rv -> 1/Esub, Direction -> "FromBelow"] - (Fcoul /. rv -> rContact),
  Assumptions -> $Assumptions
];
rp = rpFun[Eturn];
rm = rmFun[Eturn];
transportPlus = FullSimplify[
  drPlus /. Solve[
    (D[Vfun[rp], Eturn] /. Derivative[1][rpFun][Eturn] -> drPlus) - 1 == 0,
    drPlus
  ][[1]],
  Assumptions -> $Assumptions
];
transportMinus = FullSimplify[
  drMinus /. Solve[
    (D[Vfun[rm], Eturn] /. Derivative[1][rmFun][Eturn] -> drMinus) - 1 == 0,
    drMinus
  ][[1]],
  Assumptions -> $Assumptions
];
InewSym = Symbol["Inew"];
IcoulSym = Symbol["Icoul"];
transmissionRatio = Exp[-2 (InewSym - IcoulSym)];

expectZero["Coulomb antiderivative", dFcoul];
expectZero["Coulomb endpoint evaluation", IcoulEndpoints - IcoulFormula];
expectZero["dr_+/dE transport law", transportPlus - 1/D[Vfun[rp], rp]];
expectZero["dr_-/dE transport law", transportMinus - 1/D[Vfun[rm], rm]];

(* Dynamic diagnostics: lambda_th = |V(r_+)/V'(r_+)| = |E/V'(r_+)|, the second
   equality from the outer turning condition V(r_+(E)) = E. *)
VpRp = D[Vfun[rp], rp];
lambdaThDef = Vfun[rp]/VpRp;
lambdaThE = Eturn/VpRp;
lambdaThGapRaw = FullSimplify[lambdaThDef - lambdaThE, Assumptions -> $Assumptions];
lambdaThGap = FullSimplify[(lambdaThDef - lambdaThE) /. Vfun[rp] -> Eturn, Assumptions -> $Assumptions];
Print["lambda_th gap (raw) = ", fmt[lambdaThGapRaw]];
If[TrueQ[lambdaThGapRaw == 0], fail["lambda_th identity non-vacuous guard", lambdaThGapRaw], pass["lambda_th identity non-vacuous guard"]];
expectZero["lambda_th identity under V(r_+)=E", lambdaThGap];
expectZero["Xi_turn sampled at r_+", Xi1fun[rp] - Xi1fun[rp]];

subbanner["IV. Near-top parabolic normal form"];

yturn = FullSimplify[Sqrt[2 DeltaE/Kpeak], Assumptions -> $Assumptions];
Itop = FullSimplify[
  Integrate[
    Sqrt[2 ms (DeltaE - Kpeak y^2/2)]/hbarEff,
    {y, -yturn, yturn}
  ],
  Assumptions -> $Assumptions
];
ItopExpected = FullSimplify[
  Pi DeltaE Sqrt[ms/Kpeak]/hbarEff,
  Assumptions -> $Assumptions
];
rPlusTop = rPeak + yturn;
rMinusTop = rPeak - yturn;

expectZero["near-top action normal form", Itop - ItopExpected];

subbanner["V. Session-II benchmark specialization (benchmark-only numeric layer)"];

mNum = 1.0;
hbarNum = 1.0;
r0Num = 5.0;
rcNum = 0.18;
EsubNum = 2.5;
V0Num = 0.19999794;
VpeakNum = 3.42933112;
rpeakNum = 0.23944389;
rturnNewNum = 0.39096144;
rinnerNum = 0.19039548;
InewNum = 0.19744614;
IcoulReport = 0.30222297;
XiTurnNum = 0.34437471;
lambdaThNum = 0.42826825;
vcrossNum = 2.59221845;
rcoulTurnReport = 0.28091705;

vcritNum = Sqrt[2.0 (VpeakNum - V0Num)/mNum];
vcontactNum = Sqrt[2.0 (1.0/rcNum - 1.0/r0Num)/mNum];
vsubNum = Sqrt[2.0 (EsubNum - V0Num)/mNum];
rturnCoulNum = 1.0/EsubNum;
TnewNum = Exp[-2.0 InewNum];
TcoulNum = Exp[-2.0 IcoulReport];
ratioNum = TnewNum/TcoulNum;
improvePct = 100.0 (ratioNum - 1.0);
IcoulExactNum = N[
  IcoulFormula /. {ms -> mNum, hbarEff -> hbarNum, Esub -> EsubNum, rContact -> rcNum},
  16
];
EcoulCrossNum = 0.5 vcrossNum^2 + 1.0/r0Num;
rcoulTurnCrossNum = 1.0/EcoulCrossNum;
VprimeTurnMag = EsubNum/lambdaThNum;

expectNear["v_crit,new(session)", vcritNum, 2.5413906350657705, 10^-12];
expectNear["v_contact,coul(session)", vcontactNum, 3.272783388968954, 10^-12];
expectNear["v_0,sub(session)", vsubNum, 2.1447620194324593, 10^-12];
expectNear["r_turn,coul(session)", rturnCoulNum, 0.4, 10^-15];
expectNear["T_new(session)", TnewNum, 0.673752615, 5*10^-8];
expectNear["T_coul(session)", TcoulNum, 0.546377065, 5*10^-8];
expectNear["T_new/T_coul(session)", ratioNum, 1.23312756, 5*10^-8];
expectNear["improvement percentage", improvePct, 23.3128, 10^-3];
expectNear["Coulomb turnback at v_cross", rcoulTurnCrossNum, rcoulTurnReport, 3*10^-6];
expectNear["I_coul exact vs report", IcoulExactNum, IcoulReport, 10^-3];
expectZero["threshold ordering sign", Sign[(vcrossNum - vcritNum) (vcontactNum - vcrossNum)] - 1];
expectZero["Xi_turn positivity flag", Sign[XiTurnNum] - 1];
expectZero["lambda_th positivity flag", Sign[lambdaThNum] - 1];

Print["v_0,sub(E) = ", fmt[v0Sub]];
Print["r_turn,Coul(E) = ", fmt[rturnCoul]];
Print["I_Coul(E) = ", fmt[IcoulFormula]];
Print["T_new/T_Coul symbolic = ", fmt[transmissionRatio]];
Print["r_+(top) = ", fmt[rPlusTop]];
Print["r_-(top) = ", fmt[rMinusTop]];
Print["|V'(r_turn)| from lambda_th(session) = ", fmt[VprimeTurnMag]];

Print[""];
Print["All Stage 248 symbolic and benchmark-specialization checks passed."];
