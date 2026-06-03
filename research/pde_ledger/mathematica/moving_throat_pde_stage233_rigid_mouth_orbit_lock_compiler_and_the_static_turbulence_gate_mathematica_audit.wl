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
fmt[expr_] := StringReplace[ToString[InputForm[expr]], "Global`" -> ""];
stripConditional[expr_] := expr /. ConditionalExpression[e_, _] :> e;

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  detail -> ", fmt[detail]]];
  Exit[1];
);

cleanExpr[expr_] := Module[{res},
  res = stripConditional[expr];
  res = FullSimplify[Cancel[Together[res]], Assumptions -> $Assumptions];
  res = stripConditional[res];
  FullSimplify[res, Assumptions -> $Assumptions]
];

expectZero[name_String, expr_] := Module[{res},
  res = cleanExpr[expr];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectSmall[name_String, expr_, tol_] := Module[{res, ok},
  res = cleanExpr[expr];
  ok = FullSimplify[Abs[res] <= tol, Assumptions -> $Assumptions];
  Print[name, " residual = ", fmt[res]];
  If[TrueQ[ok], pass[name], fail[name, res]];
];

exactDecimal[s_String] := Module[
  {txt, sign, parts, whole, frac, scale, value},
  txt = s;
  sign = 1;
  If[StringStartsQ[txt, "-"],
    sign = -1;
    txt = StringDrop[txt, 1]
  ];
  parts = StringSplit[txt, "."];
  whole = ToExpression[parts[[1]]];
  If[Length[parts] == 1,
    value = whole,
    frac = parts[[2]];
    scale = 10^StringLength[frac];
    value = (whole scale + ToExpression[frac])/scale
  ];
  sign Rationalize[value, 0]
];

Clear[
  dlnRtr, dlnNstar, dlnEpsEta, bstar, epsEtaStar,
  n0, d0, n01, d01, deltaNorm, tQuad, mhat0, pcrit, pbar
];

$Assumptions = (
  Element[
    {
      dlnRtr, dlnNstar, dlnEpsEta, bstar, epsEtaStar,
      n0, d0, n01, d01, deltaNorm, tQuad, mhat0, pcrit, pbar
    },
    Reals
  ]
  && epsEtaStar != 1
  && n0 != 0 && d0 != 0
  && deltaNorm >= 0 && tQuad > 0 && mhat0 > 0 && pcrit > 0 && pbar > 0
);

banner["STAGE 233 -- RIGID-MOUTH ORBIT-LOCK COMPILER MATHEMATICA AUDIT"];

subbanner["M1-M2. Coefficient compiler and track lock"];

theta1 = dlnRtr;
xi1 = dlnNstar - bstar dlnRtr;
ceta = epsEtaStar/(1 - epsEtaStar);
r1 = -ceta dlnEpsEta - xi1;

thetaCoeff = Coefficient[theta1, dlnRtr];
thetaConstant = Coefficient[theta1, dlnRtr, 0];
xiCoeff = Coefficient[xi1, dlnRtr];
xiConstant = Coefficient[xi1, dlnRtr, 0];
r1TrackFromCoeff = -ceta dlnEpsEta - xiConstant;

expectZero["M1 Theta1 dlnRtr coefficient", thetaCoeff - 1];
expectZero["M1 Xi1 dlnRtr coefficient", xiCoeff + bstar];
expectZero["M2 Theta1 constant term", thetaConstant];
expectZero["M2 Theta1 track lock from coefficient data", thetaConstant + thetaCoeff 0];
expectZero["M2 Xi1 track lock from coefficient data", xiConstant - dlnNstar];
expectZero["M2 R1 track lock from coefficient data", r1TrackFromCoeff + ceta dlnEpsEta + dlnNstar];

subbanner["M3-M4. Prefactor quotient and operator rigidity"];

p0 = n0/d0;
p1 = (n01 d0 - n0 d01)/d0^2;
xiLoad = n01/n0 - d01/d0;
quotientResidual = Together[p1/p0 - xiLoad];

expectZero["M3 prefactor identity", quotientResidual];
expectZero["M4 operator-rigid load scalar", (xiLoad /. d01 -> 0) - n01/n0];

subbanner["M5-M7. Transported ceiling and Pbar forms"];

pbarExpr = (deltaNorm + tQuad)/mhat0^2;
gateRhs = pcrit mhat0^2/(deltaNorm + tQuad) - 1;
pbarGate = pcrit/pbarExpr - 1;
abstractPbarGate = pcrit/pbar - 1;

expectZero["M5 transported ceiling equivalence", gateRhs - pbarGate];

calibratedFromPbar = FullSimplify[pbarGate /. deltaNorm -> 0, Assumptions -> $Assumptions];
expectZero["M6 calibrated branch through Pbar form", calibratedFromPbar - (pcrit mhat0^2/tQuad - 1)];

expectZero["M7 abstract Pbar form specializes to gate RHS", (abstractPbarGate /. pbar -> pbarExpr) - gateRhs];

subbanner["M8. Carried budget recovery from independent critical pressures"];

(* Critical-pressure literals are the same independent inputs used in the Stage 233 SymPy fix:
   scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py ceilings dict. *)
pbarCompat = exactDecimal["0.002069792318062885"];
pcritRobust = exactDecimal["0.0028313316855593175"];
pcritNonempty = exactDecimal["0.0035965105896846573"];
budgetRobust = exactDecimal["0.367930328492646"];
budgetNonempty = exactDecimal["0.737619063660757"];
budgetTol = exactDecimal["0.000000000000001"];

expectZero["M8 symbolic budget form", (pcrit/pbar - 1) - abstractPbarGate];
expectSmall["M8 robust numeric budget", pcritRobust/pbarCompat - 1 - budgetRobust, budgetTol];
expectSmall["M8 nonempty numeric budget", pcritNonempty/pbarCompat - 1 - budgetNonempty, budgetTol];

Print[""];
Print["All Stage 233 Mathematica checks passed."];
