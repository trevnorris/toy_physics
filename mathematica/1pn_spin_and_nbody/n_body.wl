(* ===================================================================== *)
(* SUPERFLUID DEFECT TOY MODEL: PAPER III N-BODY / VECTOR SECTOR (UPDATED) *)
(* ===================================================================== *)
(* v2: Fixes Solve::ivar by introducing an explicit symbol a2 = alpha^2.  *)
(* ===================================================================== *)

ClearAll["Global`*"];

Print["\n======================================================="];
Print[" PAPER III: N-BODY / VECTOR SECTOR (UPDATED TRANSLATIONAL WAKE)"];
Print["======================================================="];

(* --------------------------------------------------------------------- *)
(* PART 0: Targets from the EIH 1PN Lagrangian (cross terms only)         *)
(* --------------------------------------------------------------------- *)

TargetPara = -7/2;  (* coefficient of vA·vB *)
TargetLong = -1/2;  (* coefficient of (vA·n)(vB·n) *)

Print["\nEIH Targets (cross terms only):"];
Print["  TargetPara (vA·vB) = ", TargetPara];
Print["  TargetLong ((vA·n)(vB·n)) = ", TargetLong];
Print["  Target ratio (Long/Para) = ", N[TargetLong/TargetPara]];

(* --------------------------------------------------------------------- *)
(* PART 1: Translational wake basis in Fourier space                      *)
(* --------------------------------------------------------------------- *)

PL[k_] := Outer[Times, k, k]/(k.k);
PT[k_] := IdentityMatrix[3] - PL[k];

uWake[k_, v_, aH_, alpha_] := Module[{kk = Sqrt[k.k]},
  I*( (PT[k].v)/kk + aH*Cross[k, v]/(kk^2) + alpha*(k*(k.v))/(kk^3) )
];

Print["\nWake basis:"];
Print["  uT  ~ PT[v]/k  (translation/shear transverse)"];
Print["  uH  ~ (k×v)/k^2 (helical transverse; optional)"];
Print["  uL  ~ k(k·v)/k^3 (longitudinal)"];

(* --------------------------------------------------------------------- *)
(* PART 2: Overlap integral -> cross-term coefficients                    *)
(* --------------------------------------------------------------------- *)

Print["\n... Computing overlap integral (symbolic) ..."];

veck = {kx, ky, kz};
vecvA = {vAx, vAy, vAz};
vecvB = {vBx, vBy, vBz};

aHsym = aH;
alphasym = alpha;

uAk = uWake[veck, vecvA, aHsym, alphasym];
uBnegk = uWake[-veck, vecvB, aHsym, alphasym];
dotProduct = Simplify[uBnegk . uAk];

integrand = dotProduct /. {
  kx -> k*Sin[t]*Cos[p],
  ky -> k*Sin[t]*Sin[p],
  kz -> k*Cos[t]
};

fullIntegrand = integrand * k^2 * Sin[t] * Exp[I*k*d*Cos[t]];

angInt = Integrate[fullIntegrand, {p, 0, 2*Pi}, {t, 0, Pi},
  Assumptions -> {d > 0, k > 0, Element[{aHsym, alphasym}, Reals]}
];

radInt = Integrate[angInt, {k, 0, Infinity},
  Assumptions -> {d > 0, Element[{aHsym, alphasym}, Reals]}
];

CoeffVV = Coefficient[radInt, vAx*vBx];
CoeffZZ = Coefficient[radInt, vAz*vBz];

VecParaShape = Simplify[CoeffVV * d];
VecLongShape = Simplify[(CoeffZZ - CoeffVV) * d];

Print["\nDerived geometric shapes (up to overall coupling K):"];
Print["  VecParaShape  (vA·vB)      = ", VecParaShape];
Print["  VecLongShape  ((vA·n)(vB·n)) = ", VecLongShape];

RatioModel = Simplify[VecLongShape/VecParaShape];
RatioTarget = Simplify[TargetLong/TargetPara];

Print["\nModel ratio (Long/Para) = ", RatioModel];
Print["Target ratio (Long/Para) = ", RatioTarget];

(* --------------------------------------------------------------------- *)
(* PART 3: Solve ratio constraint and choose a minimal real solution      *)
(* --------------------------------------------------------------------- *)

Print["\n... Solving ratio constraint and choosing minimal wake ..."];

(* IMPORTANT: Mathematica cannot Solve[...] for alpha^2 directly because alpha^2
   is not a valid Solve variable. Introduce a2 := alpha^2 as an independent symbol. *)
a2 = Unique["a2"];

RatioModelA2 = FullSimplify[RatioModel /. alphasym^2 -> a2, Assumptions -> Element[aHsym, Reals]];
solA2 = FullSimplify[Solve[RatioModelA2 == RatioTarget, a2], Assumptions -> Element[aHsym, Reals]];

Print["\nConstraint solutions for alpha^2 (via a2 = alpha^2):"];
Print[solA2 /. a2 -> alphasym^2];

(* Choose minimal wake: minimize aH^2 + alpha^2 subject to ratio constraint *)
min = Minimize[
  {aHsym^2 + alphasym^2, RatioModel == RatioTarget && 0 <= aHsym^2 <= 1 && alphasym^2 >= 0},
  {aHsym, alphasym},
  Reals
];

Print["\nMinimize result {minValue, {aH->..., alpha->...}}:"];
Print[min];

rules = min[[2]];
aHval = aHsym /. rules;
alphaval = alphasym /. rules;

Print["\nSelected real parameters:"];
Print["  aH = ", aHval];
Print["  alpha = ", alphaval];
Print["  alpha^2 = ", Simplify[alphaval^2]];

(* --------------------------------------------------------------------- *)
(* PART 4: Fix overall coupling K to match TargetPara                     *)
(* --------------------------------------------------------------------- *)

Kval = FullSimplify[TargetPara/VecParaShape /. rules];
Print["\nOverall coupling K needed to match TargetPara:"];
Print["  K = ", Kval];

PredPara = FullSimplify[Kval * VecParaShape /. rules];
PredLong = FullSimplify[Kval * VecLongShape /. rules];

Print["\nCheck (should reproduce EIH cross coefficients):"];
Print["  Predicted Para = ", PredPara, "   (target ", TargetPara, ")"];
Print["  Predicted Long = ", PredLong, "   (target ", TargetLong, ")"];

If[Element[alphaval, Reals] && Element[aHval, Reals],
  Print["\nSUCCESS: Found a real-parameter wake that matches EIH cross terms."],
  Print["\nFAILURE: Parameters are not real (unexpected)."]
];

(* --------------------------------------------------------------------- *)
(* PART 5: Notes / placeholders for static nonlinearity (cavitation)      *)
(* --------------------------------------------------------------------- *)

Print["\n--- NOTE: Static (G^2) / cavitation sector ---"];
Print["This script focuses on the VECTOR cross-term tensor match."];
Print["Static nonlinearity (density-dependent masses, cavitation, etc.)"];
Print["should be treated separately and does not alter the EIH cross tensor."];
Print["\nDone."];

(*"
Output:

EIH Targets (cross terms only):
  TargetPara (vA·vB) = -7/2
  TargetLong ((vA·n)(vB·n)) = -1/2
  Target ratio (Long/Para) = 0.14285714285714285

Wake basis:
  uT  ~ PT[v]/k  (translation/shear transverse)
  uH  ~ (k×v)/k^2 (helical transverse; optional)
  uL  ~ k(k·v)/k^3 (longitudinal)

... Computing overlap integral (symbolic) ...

Derived geometric shapes (up to overall coupling K):
  VecParaShape  (vA·vB)      = (-1 + aH^2 - alpha^2)*Pi^2
  VecLongShape  ((vA·n)(vB·n)) = (-1 + aH^2 + alpha^2)*Pi^2

Model ratio (Long/Para) = (-1 + aH^2 + alpha^2)/(-1 + aH^2 - alpha^2)
Target ratio (Long/Para) = 1/7

... Solving ratio constraint and choosing minimal wake ...

Constraint solutions for alpha^2 (via a2 = alpha^2):
{{alpha^2 -> (-3*(-1 + aH^2))/4}}

Minimize result {minValue, {aH->..., alpha->...}}:
{3/4, {aH -> 0, alpha -> -1/2*Sqrt[3]}}

Selected real parameters:
  aH = 0
  alpha = -1/2*Sqrt[3]
  alpha^2 = 3/4

Overall coupling K needed to match TargetPara:
  K = 2/Pi^2

Check (should reproduce EIH cross coefficients):
  Predicted Para = -7/2   (target -7/2)
  Predicted Long = -1/2   (target -1/2)

SUCCESS: Found a real-parameter wake that matches EIH cross terms.
"*)
