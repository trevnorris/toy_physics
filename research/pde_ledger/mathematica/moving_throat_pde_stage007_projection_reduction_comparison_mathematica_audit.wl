ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Print["STAGE 007 PROJECTION REDUCTION COMPARISON MATHEMATICA AUDIT"];

Clear[w, x, k, lambda, rho, xi, mu0, sigma, tau, epsilon, etaSym];
$Assumptions =
  lambda > 0 && rho > 0 && xi > 0 && mu0 > 0 && sigma > 0 && tau > 0 && epsilon > 0 &&
    Element[{x, k, etaSym}, Reals] && etaSym != 0;

Z[w_] := Exp[-w^2/lambda^2];
Hprofile[w_] := Exp[-w^2/rho^2];

Zint = Integrate[Z[w], {w, -Infinity, Infinity}];
Z2int = Integrate[Z[w]^2, {w, -Infinity, Infinity}];
Hint = Integrate[Hprofile[w], {w, -Infinity, Infinity}];

(* M1: Gaussian profile area. *)
m1Residual = FullSimplify[Zint - Sqrt[Pi]*lambda];
Print["M1 residual = ", fmt[m1Residual]];
If[!(FullSimplify[Zint - Sqrt[Pi]*lambda] === 0),
  Print["FAIL: M1 Gaussian profile area residual = ", fmt[m1Residual]]; Exit[1]
];
Print["PASS: M1 Gaussian profile area"];

(* M2: Gaussian squared profile area. *)
m2Residual = FullSimplify[Z2int - Sqrt[Pi/2]*lambda];
Print["M2 residual = ", fmt[m2Residual]];
If[!(FullSimplify[Z2int - Sqrt[Pi/2]*lambda] === 0),
  Print["FAIL: M2 Gaussian squared profile area residual = ", fmt[m2Residual]]; Exit[1]
];
Print["PASS: M2 Gaussian squared profile area"];

(* M2b: Gaussian gauge-weight profile area. *)
m2bResidual = FullSimplify[Hint - Sqrt[Pi]*rho];
Print["M2b residual = ", fmt[m2bResidual]];
If[!(FullSimplify[Hint - Sqrt[Pi]*rho] === 0),
  Print["FAIL: M2b Gaussian gauge-weight area residual = ", fmt[m2bResidual]]; Exit[1]
];
Print["PASS: M2b Gaussian gauge-weight area"];

Wsmooth[w_] := Exp[-w^2/sigma^2]/(Sqrt[Pi]*sigma);
Ssmooth[w_] := Exp[-w^2/tau^2]/(Sqrt[Pi]*tau);

IWZsmooth = Integrate[Wsmooth[w] Z[w], {w, -Infinity, Infinity}];
IWHsmooth = Integrate[Wsmooth[w] Hprofile[w], {w, -Infinity, Infinity}];
IWSsmooth = Integrate[Wsmooth[w] Ssmooth[w], {w, -Infinity, Infinity}];

(* M3: smooth observer/profile overlap. *)
m3Residual = FullSimplify[IWZsmooth - lambda/Sqrt[lambda^2 + sigma^2]];
Print["M3 residual = ", fmt[m3Residual]];
If[!(FullSimplify[IWZsmooth - lambda/Sqrt[lambda^2 + sigma^2]] === 0),
  Print["FAIL: M3 smooth observer/profile overlap residual = ", fmt[m3Residual]]; Exit[1]
];
Print["PASS: M3 smooth observer/profile overlap"];

(* M4: smooth observer/source overlap. *)
m4Residual = FullSimplify[IWSsmooth - 1/(Sqrt[Pi]*Sqrt[sigma^2 + tau^2])];
Print["M4 residual = ", fmt[m4Residual]];
If[!(FullSimplify[IWSsmooth - 1/(Sqrt[Pi]*Sqrt[sigma^2 + tau^2])] === 0),
  Print["FAIL: M4 smooth observer/source overlap residual = ", fmt[m4Residual]]; Exit[1]
];
Print["PASS: M4 smooth observer/source overlap"];

(* M4b: smooth observer/gauge-weight overlap. *)
m4bResidual = FullSimplify[IWHsmooth - rho/Sqrt[rho^2 + sigma^2]];
Print["M4b residual = ", fmt[m4bResidual]];
If[!(FullSimplify[IWHsmooth - rho/Sqrt[rho^2 + sigma^2]] === 0),
  Print["FAIL: M4b smooth observer/gauge-weight overlap residual = ", fmt[m4bResidual]]; Exit[1]
];
Print["PASS: M4b smooth observer/gauge-weight overlap"];

(* M4c: smooth projected gauge parameter. *)
m4cTarget = xi*lambda*Sqrt[rho^2 + sigma^2]/(rho*Sqrt[lambda^2 + sigma^2]);
m4cResidual = FullSimplify[xi*IWZsmooth/IWHsmooth - m4cTarget];
Print["M4c residual = ", fmt[m4cResidual]];
If[!(FullSimplify[xi*IWZsmooth/IWHsmooth - m4cTarget] === 0),
  Print["FAIL: M4c smooth projected gauge parameter residual = ", fmt[m4cResidual]]; Exit[1]
];
Print["PASS: M4c smooth projected gauge parameter"];

Fmut[x_, w_] := Sin[k x] + x^2 + etaSym*x*w^2;
fieldMutationLHS =
  Integrate[Wsmooth[w] Z[w] D[Fmut[x, w], x], {w, -Infinity, Infinity}];
fieldMutationDelta =
  FullSimplify[fieldMutationLHS - IWZsmooth*D[Sin[k x] + x^2, x]];

(* M5: field-mutation Gaussian moment. *)
m5Target = etaSym*lambda^3*sigma^2/(2*(lambda^2 + sigma^2)^(3/2));
m5Residual = FullSimplify[fieldMutationDelta - m5Target];
Print["M5 residual = ", fmt[m5Residual]];
If[!(FullSimplify[
    fieldMutationDelta
      - etaSym*lambda^3*sigma^2/(2*(lambda^2 + sigma^2)^(3/2))] === 0),
  Print["FAIL: M5 field-mutation Gaussian moment residual = ", fmt[m5Residual]]; Exit[1]
];
Print["PASS: M5 field-mutation Gaussian moment"];

Jmut[x_, w_] := Ssmooth[w]*(Cos[k x] + x) + etaSym*x*w^2*Ssmooth[w];
sourceMutationRHS =
  Integrate[Wsmooth[w] mu0 Jmut[x, w], {w, -Infinity, Infinity}];
sourceMutationDelta =
  FullSimplify[sourceMutationRHS - mu0*IWSsmooth*(Cos[k x] + x)];

(* M6: source-mutation Gaussian moment. *)
m6Target =
  etaSym*mu0*x*sigma^2*tau^2/(2*Sqrt[Pi]*(sigma^2 + tau^2)^(3/2));
m6Residual = FullSimplify[sourceMutationDelta - m6Target];
Print["M6 residual = ", fmt[m6Residual]];
If[!(FullSimplify[
    sourceMutationDelta
      - etaSym*mu0*x*sigma^2*tau^2/(2*Sqrt[Pi]*(sigma^2 + tau^2)^(3/2))] === 0),
  Print["FAIL: M6 source-mutation Gaussian moment residual = ", fmt[m6Residual]]; Exit[1]
];
Print["PASS: M6 source-mutation Gaussian moment"];

Wmatch[w_] := Z[w]/Zint;
IWZmatch = Integrate[Wmatch[w] Z[w], {w, -Infinity, Infinity}];
IWHmatch = Integrate[Wmatch[w] Hprofile[w], {w, -Infinity, Infinity}];

(* M7: matched-observer overlap. *)
m7Residual = FullSimplify[IWZmatch - 1/Sqrt[2]];
Print["M7 residual = ", fmt[m7Residual]];
If[!(FullSimplify[IWZmatch - 1/Sqrt[2]] === 0),
  Print["FAIL: M7 matched-observer overlap residual = ", fmt[m7Residual]]; Exit[1]
];
Print["PASS: M7 matched-observer overlap"];

(* M7b: matched-observer gauge-weight overlap. *)
m7bResidual = FullSimplify[IWHmatch - rho/Sqrt[lambda^2 + rho^2]];
Print["M7b residual = ", fmt[m7bResidual]];
If[!(FullSimplify[IWHmatch - rho/Sqrt[lambda^2 + rho^2]] === 0),
  Print["FAIL: M7b matched-observer gauge overlap residual = ", fmt[m7bResidual]]; Exit[1]
];
Print["PASS: M7b matched-observer gauge overlap"];

IWSmatch = Wmatch[0];
mu0ProjMatch = mu0*IWSmatch/IWZmatch;
mu0Red = mu0/Zint;
xiProjMatch = xi*IWZmatch/IWHmatch;
xiRed = xi*Zint/Hint;

(* M8: matched delta-source projection/reduction ratio. *)
m8Residual = FullSimplify[mu0ProjMatch/mu0Red - Sqrt[2]];
Print["M8 residual = ", fmt[m8Residual]];
If[!(FullSimplify[mu0ProjMatch/mu0Red - Sqrt[2]] === 0),
  Print["FAIL: M8 matched projection/reduction ratio residual = ", fmt[m8Residual]]; Exit[1]
];
Print["PASS: M8 matched projection/reduction ratio"];

(* M8b: matched projected gauge parameter against reduction-first reference. *)
m8bTarget = Sqrt[lambda^2 + rho^2]/(Sqrt[2]*lambda);
m8bResidual = FullSimplify[xiProjMatch/xiRed - m8bTarget];
Print["M8b residual = ", fmt[m8bResidual]];
If[!(FullSimplify[xiProjMatch/xiRed - m8bTarget] === 0),
  Print["FAIL: M8b matched gauge projection/reduction ratio residual = ", fmt[m8bResidual]]; Exit[1]
];
Print["PASS: M8b matched gauge projection/reduction ratio"];

(* M8c: H=Z specialization aligns the gauge parameter. *)
m8cResidual = FullSimplify[(xiProjMatch/xiRed /. rho -> lambda) - 1];
Print["M8c residual = ", fmt[m8cResidual]];
If[!(FullSimplify[(xiProjMatch/xiRed /. rho -> lambda) - 1] === 0),
  Print["FAIL: M8c H=Z matched gauge alignment residual = ", fmt[m8cResidual]]; Exit[1]
];
Print["PASS: M8c H=Z matched gauge alignment"];

Weps[w_] := Exp[-w^2/epsilon^2]/(Sqrt[Pi]*epsilon);
IWZeps = Integrate[Weps[w] Z[w], {w, -Infinity, Infinity}];
IWHeps = Integrate[Weps[w] Hprofile[w], {w, -Infinity, Infinity}];
IWSeps = Integrate[Weps[w]^2, {w, -Infinity, Infinity}];

(* M9: regulated observer/profile overlap. *)
m9Residual = FullSimplify[IWZeps - lambda/Sqrt[epsilon^2 + lambda^2]];
Print["M9 residual = ", fmt[m9Residual]];
If[!(FullSimplify[IWZeps - lambda/Sqrt[epsilon^2 + lambda^2]] === 0),
  Print["FAIL: M9 regulated observer/profile overlap residual = ", fmt[m9Residual]]; Exit[1]
];
Print["PASS: M9 regulated observer/profile overlap"];

(* M10: regulated observer self-overlap. *)
m10Residual = FullSimplify[IWSeps - Sqrt[2]/(2*Sqrt[Pi]*epsilon)];
Print["M10 residual = ", fmt[m10Residual]];
If[!(FullSimplify[IWSeps - Sqrt[2]/(2*Sqrt[Pi]*epsilon)] === 0),
  Print["FAIL: M10 regulated observer self-overlap residual = ", fmt[m10Residual]]; Exit[1]
];
Print["PASS: M10 regulated observer self-overlap"];

(* M10b: regulated observer/gauge-weight overlap. *)
m10bResidual = FullSimplify[IWHeps - rho/Sqrt[epsilon^2 + rho^2]];
Print["M10b residual = ", fmt[m10bResidual]];
If[!(FullSimplify[IWHeps - rho/Sqrt[epsilon^2 + rho^2]] === 0),
  Print["FAIL: M10b regulated gauge overlap residual = ", fmt[m10bResidual]]; Exit[1]
];
Print["PASS: M10b regulated gauge overlap"];

xiProjEps = FullSimplify[xi*IWZeps/IWHeps];
m10cTarget = xi*lambda*Sqrt[epsilon^2 + rho^2]/(rho*Sqrt[epsilon^2 + lambda^2]);
m10cResidual = FullSimplify[xiProjEps - m10cTarget];
Print["M10c residual = ", fmt[m10cResidual]];
If[!(FullSimplify[xiProjEps - m10cTarget] === 0),
  Print["FAIL: M10c regulated projected gauge parameter residual = ", fmt[m10cResidual]]; Exit[1]
];
Print["PASS: M10c regulated projected gauge parameter"];

(* M11: sharp-sampling limit. *)
m11Residual =
  FullSimplify[
    Limit[IWZeps, epsilon -> 0, Direction -> "FromAbove",
      Assumptions -> lambda > 0] - 1,
    Assumptions -> lambda > 0];
Print["M11 residual = ", fmt[m11Residual]];
If[!(Limit[IWZeps, epsilon -> 0, Direction -> "FromAbove",
      Assumptions -> lambda > 0] === 1),
  Print["FAIL: M11 sharp-sampling limit residual = ", fmt[m11Residual]]; Exit[1]
];
Print["PASS: M11 sharp-sampling limit"];

(* M11b: sharp-sampling gauge limit. *)
m11bResidual =
  FullSimplify[
    Limit[IWHeps, epsilon -> 0, Direction -> "FromAbove",
      Assumptions -> rho > 0] - 1,
    Assumptions -> rho > 0];
Print["M11b residual = ", fmt[m11bResidual]];
If[!(Limit[IWHeps, epsilon -> 0, Direction -> "FromAbove",
      Assumptions -> rho > 0] === 1),
  Print["FAIL: M11b sharp gauge-sampling limit residual = ", fmt[m11bResidual]]; Exit[1]
];
Print["PASS: M11b sharp gauge-sampling limit"];

m11cResidual =
  FullSimplify[
    Limit[xiProjEps, epsilon -> 0, Direction -> "FromAbove",
      Assumptions -> lambda > 0 && rho > 0 && xi > 0] - xi,
    Assumptions -> lambda > 0 && rho > 0 && xi > 0];
Print["M11c residual = ", fmt[m11cResidual]];
If[!(FullSimplify[m11cResidual] === 0),
  Print["FAIL: M11c sharp projected gauge-parameter limit residual = ", fmt[m11cResidual]]; Exit[1]
];
Print["PASS: M11c sharp projected gauge-parameter limit"];

Print["STATUS: PASS"];
Exit[0];
