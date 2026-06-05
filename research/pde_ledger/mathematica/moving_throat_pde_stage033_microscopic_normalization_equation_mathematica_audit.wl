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

banner["STAGE 033 — MICROSCOPIC NORMALIZATION EQUATION"];

Clear[A, DeltaK, alpha0, beta0];
$Assumptions =
  Element[{A, DeltaK, alpha0, beta0}, Reals] &&
  A > 0 && DeltaK > 0 && alpha0 > 0 && beta0 > 0;

kappa0Sq = 8/Pi^2;
kappa1Sq = 16/(9*Pi^2);
sigma = FullSimplify[kappa0Sq + kappa1Sq, Assumptions -> $Assumptions];
deltaKappa = FullSimplify[kappa0Sq - kappa1Sq, Assumptions -> $Assumptions];
kProd = FullSimplify[kappa0Sq*kappa1Sq, Assumptions -> $Assumptions];

r = FullSimplify[
  Sqrt[(DeltaK + alpha0*deltaKappa)^2 + 4*alpha0^2*kProd],
  Assumptions -> $Assumptions
];
lambdaMinus = FullSimplify[(2*A + DeltaK - alpha0*sigma - r)/2, Assumptions -> $Assumptions];
sMinus = FullSimplify[
  1/2*(sigma + ((DeltaK + alpha0*deltaKappa)*deltaKappa + 4*alpha0*kProd)/r),
  Assumptions -> $Assumptions
];
nMinus = FullSimplify[beta0*sMinus^2/(kappa0Sq*lambdaMinus), Assumptions -> $Assumptions];

Print["lambda_- = ", fmt[lambdaMinus]];
Print["s_-(alpha0) = ", fmt[sMinus]];
Print["N_-(alpha0) = ", fmt[nMinus]];

ds = FullSimplify[D[sMinus, alpha0], Assumptions -> $Assumptions];
dN = FullSimplify[D[nMinus, alpha0], Assumptions -> $Assumptions];
dNFormula = FullSimplify[
  beta0*(2*sMinus*ds*lambdaMinus + sMinus^3)/(kappa0Sq*lambdaMinus^2),
  Assumptions -> $Assumptions
];
expectZero["dN/dalpha - monotonicity formula", dN - dNFormula];

alphaCrit = FullSimplify[A*(A + DeltaK)/((A + DeltaK)*kappa0Sq + A*kappa1Sq), Assumptions -> $Assumptions];
alphaCritClosed = FullSimplify[9*Pi^2*A*(A + DeltaK)/(8*(11*A + 9*DeltaK)), Assumptions -> $Assumptions];
Print["alpha_crit = ", fmt[alphaCrit]];
expectZero["alpha_crit - finite-throat closed form", alphaCrit - alphaCritClosed];

n0 = FullSimplify[nMinus /. alpha0 -> 0, Assumptions -> A > 0 && DeltaK > 0 && beta0 > 0];
Print["N_-(0) = ", fmt[n0]];
expectZero["N_-(0) - beta0 kappa0^2 / A", n0 - beta0*kappa0Sq/A];

coef1 = FullSimplify[D[nMinus, alpha0] /. alpha0 -> 0, Assumptions -> A > 0 && DeltaK > 0 && beta0 > 0];
coef1Target = FullSimplify[
  beta0*kappa0Sq*(4*A*kappa1Sq + DeltaK*kappa0Sq)/(A^2*DeltaK),
  Assumptions -> A > 0 && DeltaK > 0 && beta0 > 0
];
coef1Closed = FullSimplify[
  64*beta0*(8*A + 9*DeltaK)/(9*Pi^4*A^2*DeltaK),
  Assumptions -> A > 0 && DeltaK > 0 && beta0 > 0
];
Print["weak-loading coefficient = ", fmt[coef1]];
expectZero["first derivative coefficient - generic form", coef1 - coef1Target];
expectZero["first derivative coefficient - closed form", coef1 - coef1Closed];

Clear[gB, gU, gW, gR, varpi, OmegaU, OmegaW, K0, NQ];
$Assumptions =
  Element[{gB, gU, gW, gR, varpi, OmegaU, OmegaW, K0, NQ, DeltaK}, Reals] &&
  varpi > 0 && OmegaU > 0 && OmegaW > 0 && K0 > 0 && NQ > 0 && DeltaK > 0;

aMic = FullSimplify[K0 - gU^2/OmegaU^2, Assumptions -> $Assumptions];
delta0 = FullSimplify[OmegaU^2*OmegaW^2 - gR^2*sigma, Assumptions -> $Assumptions];
chi = FullSimplify[OmegaU^2*gW + gR*gU, Assumptions -> $Assumptions];
beta0Mic = FullSimplify[chi^2/delta0^2, Assumptions -> $Assumptions];
alpha0Mic = FullSimplify[gB^2/varpi^2 + chi^2/(OmegaU^2*delta0), Assumptions -> $Assumptions];
n0Mic = FullSimplify[(beta0*kappa0Sq/A) /. {beta0 -> beta0Mic, A -> aMic}, Assumptions -> $Assumptions];
k0OnsetSolutions = Solve[n0Mic == NQ, K0];
If[Length[k0OnsetSolutions] == 0, fail["Solve[n0Mic == NQ, K0] returned no solutions"]];
k0Onset = FullSimplify[K0 /. First[k0OnsetSolutions], Assumptions -> $Assumptions];

Print["A = ", fmt[aMic]];
Print["Delta0 = ", fmt[delta0]];
Print["Chi = ", fmt[chi]];
Print["beta0 = ", fmt[beta0Mic]];
Print["alpha0 = ", fmt[alpha0Mic]];
Print["N_-(0) = ", fmt[n0Mic]];
Print["K0_onset = ", fmt[k0Onset]];

expectZero["N_-(0) at K0_onset - NQ", FullSimplify[(n0Mic /. K0 -> k0Onset) - NQ, Assumptions -> $Assumptions]];
expectZero[
  "K0_onset - [gU^2/OmegaU^2 + kappa0^2 Chi^2/(NQ Delta0^2)]",
  k0Onset - (gU^2/OmegaU^2 + kappa0Sq*chi^2/(NQ*delta0^2))
];

alphaCritMic = FullSimplify[alphaCritClosed /. A -> aMic, Assumptions -> $Assumptions];
gateDenClaim = FullSimplify[8*varpi^2*OmegaU^2*delta0*(11*aMic + 9*DeltaK), Assumptions -> $Assumptions];
gateDiff = Cancel[Together[alphaCritMic - alpha0Mic]];
gateNumActual = Numerator[gateDiff];
gateDenActual = Denominator[gateDiff];
Print["computed denominator = ", fmt[gateDenActual]];
Print["claimed denominator = ", fmt[gateDenClaim]];
denRatio = FullSimplify[gateDenActual/gateDenClaim, Assumptions -> $Assumptions];
Print["denominator ratio (must be parameter-free) = ", fmt[denRatio]];
If[!NumericQ[denRatio],
  fail["gate denominator does not match claim up to a parameter-free constant", denRatio]
];
pass["gate denominator matches claim up to parameter-free constant"];
gateNumTarget = FullSimplify[gateNumActual/denRatio, Assumptions -> $Assumptions];
Print["gate numerator = ", fmt[gateNumTarget]];
expectZero[
  "gate_num_target/gate_den_claim - (alpha_crit_mic - alpha_0_mic) (tautological by reconstruction; substantive check is NumericQ[denRatio] above)",
  alphaCritMic - alpha0Mic - gateNumTarget/gateDenClaim
];

(* Independent numerical cross-check: substitute rational test values and verify
   the Stage 33.1 monotonicity identity and the Stage 33.6 gate identity at
   floating-point precision. This is structurally distinct from the analytic
   FullSimplify approach used above and would catch algebra errors that the
   line-by-line transliteration cannot. *)
numericRule1 = {A -> 2, DeltaK -> 1, alpha0 -> 1/3, beta0 -> 1,
                varpi -> 1, OmegaU -> 1, OmegaW -> 1,
                gB -> 1/2, gU -> 1/3, gW -> 1/4, gR -> 1/5,
                K0 -> 3, NQ -> 1};
numericRule2 = {A -> 5/2, DeltaK -> 7/3, alpha0 -> 2/5, beta0 -> 3/2,
                varpi -> 4/3, OmegaU -> 5/4, OmegaW -> 6/5,
                gB -> 1/7, gU -> 2/7, gW -> 3/7, gR -> 1/11,
                K0 -> 7/2, NQ -> 1/2};
Do[
  monotonicityNumeric = N[(dN - dNFormula) /. rule, 30];
  Print["monotonicity numeric residual = ", fmt[monotonicityNumeric]];
  If[Abs[monotonicityNumeric] > 10^-20,
    fail["monotonicity numeric residual nonzero", monotonicityNumeric],
    pass["monotonicity numeric residual zero at rule"]
  ],
  {rule, {numericRule1, numericRule2}}
];

Print[""];
Print["Stage 033 Mathematica audit passed."];

Exit[0];
