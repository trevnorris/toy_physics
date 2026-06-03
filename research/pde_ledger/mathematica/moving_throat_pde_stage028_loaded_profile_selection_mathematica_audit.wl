ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 028 — LOADED PROFILE SELECTION"];

Clear[alpha, theta, Keta, TOmega, Tw, L, xvar, eps];
$Assumptions =
  Element[{alpha, theta, Keta, TOmega, Tw, L, xvar, eps}, Reals] &&
  alpha > 0 && Tw > 0 && L > 0 && eps > 0;

kappa0 = FullSimplify[2*Sqrt[2]/Pi, Assumptions -> $Assumptions];
kappa1 = FullSimplify[-4/(3*Pi), Assumptions -> $Assumptions];
deltaK = FullSimplify[Tw*Pi^2/L^2, Assumptions -> $Assumptions];
K0 = Keta + 6*TOmega;
K1 = K0 + deltaK;
v = {kappa0, kappa1};

kBare = {{K0, 0}, {0, K1}};
kEff = FullSimplify[kBare - alpha*Outer[Times, v, v], Assumptions -> $Assumptions];

Print["v = ", fmt[v]];
Print["K_bare = ", fmt[kBare]];
Print["K_eff = ", fmt[kEff]];

Print["kappa0^2 - kappa1^2 = ", fmt[FullSimplify[kappa0^2 - kappa1^2, Assumptions -> $Assumptions]]];
Print["2 kappa0 kappa1 = ", fmt[FullSimplify[2*kappa0*kappa1, Assumptions -> $Assumptions]]];

trEff = FullSimplify[Tr[kEff], Assumptions -> $Assumptions];
detEff = FullSimplify[Det[kEff], Assumptions -> $Assumptions];
trExpected = FullSimplify[K0 + K1 - alpha*(kappa0^2 + kappa1^2), Assumptions -> $Assumptions];
detExpected = FullSimplify[K0*K1 - alpha*(K1*kappa0^2 + K0*kappa1^2), Assumptions -> $Assumptions];

Print["trace(K_eff) = ", fmt[trEff]];
Print["det(K_eff) = ", fmt[detEff]];
expectZero["trace - expected", trEff - trExpected];
expectZero["det - expected", detEff - detExpected];

disc = FullSimplify[
  (deltaK + alpha*(kappa0^2 - kappa1^2))^2 + 4*alpha^2*kappa0^2*kappa1^2,
  Assumptions -> $Assumptions
];
lambdaMinus = FullSimplify[(trExpected - Sqrt[disc])/2, Assumptions -> $Assumptions];
lambdaPlus = FullSimplify[(trExpected + Sqrt[disc])/2, Assumptions -> $Assumptions];
eigvalsDirect = Eigenvalues[kEff];
expectZero[
  "Eigenvalues[kEff] sum vs trace",
  FullSimplify[Total[eigvalsDirect] - (lambdaMinus + lambdaPlus), Assumptions -> $Assumptions]
];
expectZero[
  "Eigenvalues[kEff] product vs determinant",
  FullSimplify[Times @@ eigvalsDirect - lambdaMinus*lambdaPlus, Assumptions -> $Assumptions]
];
charResidual = Expand[(xvar - lambdaMinus)*(xvar - lambdaPlus) - (xvar^2 - trEff*xvar + detEff)];
expectZero["characteristic factorization", charResidual];

Print["lambda_- = ", fmt[lambdaMinus]];
Print["lambda_+ = ", fmt[lambdaPlus]];

q = {Cos[theta], Sin[theta]};
energy = FullSimplify[(q.kEff.q)/2, Assumptions -> $Assumptions];
dEnergy = FullSimplify[TrigExpand[D[energy, theta]], Assumptions -> $Assumptions];
rhs = FullSimplify[
  2*alpha*kappa0*kappa1/(deltaK + alpha*(kappa0^2 - kappa1^2)),
  Assumptions -> $Assumptions
];
stationarityExpected = FullSimplify[
  (deltaK + alpha*(kappa0^2 - kappa1^2))*Sin[2*theta] - 2*alpha*kappa0*kappa1*Cos[2*theta],
  Assumptions -> $Assumptions
];

expectZero["dE/dtheta - stationarity_expected/2", dEnergy - stationarityExpected/2];
Print["tan(2 theta) = ", fmt[rhs]];
Print["-tan(2 theta) = ", fmt[FullSimplify[-rhs, Assumptions -> $Assumptions]]];
(* Manifest positivity of -kappa0*kappa1: kappa0 > 0, kappa1 < 0. *)
expectZero["kappa0 sign check (kappa0 > 0)", kappa0 - Abs[kappa0]];
expectZero["kappa1 sign check (kappa1 < 0)", kappa1 + Abs[kappa1]];

weakCoefficient = FullSimplify[SeriesCoefficient[rhs/2, {alpha, 0, 1}], Assumptions -> Tw > 0 && L > 0];
Print["weak-loading theta coefficient = ", fmt[weakCoefficient]];
expectZero["weak-loading coefficient - kappa0 kappa1/deltaK", weakCoefficient - kappa0*kappa1/deltaK];

strongLimit = Block[
  {$Assumptions = Element[{Tw, L}, Reals] && Tw > 0 && L > 0},
  FullSimplify[Limit[rhs, alpha -> Infinity], Assumptions -> $Assumptions]
];
tMax = FullSimplify[kappa1/kappa0, Assumptions -> $Assumptions];
tan2TMax = FullSimplify[2*tMax/(1 - tMax^2), Assumptions -> $Assumptions];
Print["lim_{alpha->Infinity} tan(2 theta) = ", fmt[strongLimit]];
Print["tan(theta_max) = ", fmt[tMax]];
expectZero["strong-loading limit - tan(2 theta_max)", strongLimit - tan2TMax];
expectZero["tan(theta_max) + Sqrt[2]/3", tMax + Sqrt[2]/3];

alphaCritSolved = alpha /. First[Solve[detEff == 0, alpha]] /. ConditionalExpression[value_, _] :> value;
expectZero[
  "alphaCrit solved vs ratio closed form",
  FullSimplify[alphaCritSolved - K0*K1/(K1*kappa0^2 + K0*kappa1^2), Assumptions -> $Assumptions]
];
alphaCrit = FullSimplify[K0*K1/(K1*kappa0^2 + K0*kappa1^2), Assumptions -> $Assumptions];
alphaCritClosed = FullSimplify[9*Pi^2*K0*K1/(8*(11*K0 + 9*deltaK)), Assumptions -> $Assumptions];
Print["alpha_crit = ", fmt[alphaCrit]];
expectZero["alpha_crit - finite-throat closed form", alphaCrit - alphaCritClosed];
expectZero["det(alpha_crit)", detEff /. alpha -> alphaCrit];

detBelow = FullSimplify[Factor[detEff /. alpha -> alphaCrit*(1 - eps)], Assumptions -> $Assumptions];
Print["det(alpha_crit*(1-eps)) = ", fmt[detBelow]];

Print[""];
Print["Stage 028 Mathematica audit passed."];

Exit[0];
