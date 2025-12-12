(* ===================================================================== *)
(* PAPER III REBUILD: TRANSLATIONAL WAKE -> EIH CROSS-TERM MATCHING       *)
(* ===================================================================== *)
(* What this script does                                                  *)
(*   1) Computes the pairwise interaction kernel from an *isotropic*       *)
(*      linear-response wake decomposition in Fourier space:              *)
(*         u(k;v) = u_T(k;v) + u_H(k;v) + u_L(k;v)                         *)
(*      where:                                                            *)
(*         u_T  : transverse projector wake   (translation / shear)       *)
(*         u_H  : helical transverse wake     (vorticity-like)            *)
(*         u_L  : longitudinal (compressible) wake                         *)
(*   2) Extracts the two EIH *cross* coefficients (vA·vB and (vA·n)(vB·n)).*)
(*   3) Solves for real wake-mixing parameters that reproduce the EIH      *)
(*      cross-tensor structure {-7/2, -1/2} up to an overall coupling K.  *)
(*                                                                        *)
(* IMPORTANT (physics bookkeeping):                                        *)
(*   - This script matches ONLY the cross-velocity tensor structure        *)
(*     (vA·vB) and (vA·n)(vB·n).                                           *)
(*   - The self terms (vA^2 + vB^2) in Eq.(EIH_target) do NOT arise from   *)
(*     an overlap uA·uB integral (which is bilinear in vA,vB). Those must  *)
(*     come from the scalar/metric dressing of each body's kinetic term.  *)
(*   - Therefore: do NOT add “optical offsets” to the longitudinal cross   *)
(*     coefficient. That was the category error in older scripts.         *)
(* ===================================================================== *)

ClearAll["Global`*"];

Print["\n======================================================="];
Print[" PAPER III: TRANSLATIONAL WAKE REBUILD (CROSS-TERM MATCH)"];
Print["======================================================="];

(* --------------------------------------------------------------------- *)
(* 0. Targets from Paper III Eq.(EIH_target)                               *)
(* --------------------------------------------------------------------- *)
TargetPara = -7/2;   (* coefficient of vA·vB *)
TargetLong = -1/2;   (* coefficient of (vA·n)(vB·n) *)

Print["\nEIH Targets (cross terms only):"];
Print["  TargetPara (vA·vB) = ", TargetPara];
Print["  TargetLong ((vA·n)(vB·n)) = ", TargetLong];
Print["  Target ratio (Long/Para) = ", N[TargetLong/TargetPara]];

(* --------------------------------------------------------------------- *)
(* 1. Fourier-space wake basis (isotropic linear response)                 *)
(* --------------------------------------------------------------------- *)
(* We align separation along z in real space; the Fourier integral will    *)
(* still produce an isotropic tensor that can be read off from components. *)

veck = {kx, ky, kz};
vecvA = {vAx, vAy, vAz};
vecvB = {vBx, vBy, vBz};

(* Projectors (written so PT[veck,v] expands directly into kx,ky,kz) *)
PT[kVec_, v_] := Module[{k2loc = kVec.kVec}, v - kVec*(kVec.v)/k2loc];
PL[kVec_, v_] := Module[{k2loc = kVec.kVec}, kVec*(kVec.v)/k2loc];

(* Mode amplitudes (dimensionless). We will solve for these. *)
aT = 1;              (* set overall transverse-projector amplitude as normalization *)
aH = aHsym;          (* helical transverse amplitude (free) *)
aL = aLsym;          (* longitudinal amplitude (free) *)

(* Each mode is chosen to scale like v/k so that uA·uB ~ 1/k^2 and the k-integral gives ~1/d. *)

uT[kVec_, vVec_] := Module[{k2loc = kVec.kVec, knloc},
  knloc = Sqrt[k2loc];
  I * aT * (PT[kVec, vVec]/knloc)
];

uH[kVec_, vVec_] := Module[{k2loc = kVec.kVec},
  I * aH * (Cross[kVec, vVec]/k2loc)
];

uL[kVec_, vVec_] := Module[{k2loc = kVec.kVec, knloc},
  knloc = Sqrt[k2loc];
  I * aL * (kVec*(kVec.vVec)/(knloc^3))
];

uTotal[kVec_, vVec_] := uT[kVec, vVec] + uH[kVec, vVec] + uL[kVec, vVec];

Print["\nWake basis:"];
Print["  uT  ~ PT[v]/k  (translation/shear transverse)"];
Print["  uH  ~ (k×v)/k^2 (helical transverse; optional)"];
Print["  uL  ~ k(k·v)/k^3 (longitudinal)"];

(* --------------------------------------------------------------------- *)
(* 2. Overlap integral I(d) ~ ∫ d^3k uA(k)·uB(-k) e^{ik·d}                 *)
(* --------------------------------------------------------------------- *)
Print["\n... Computing overlap integral (symbolic) ..."];

uAk = uTotal[veck, vecvA];
uBnegk = uTotal[-veck, vecvB];

(* Dot product (bilinear in vA and vB) *)
dotProduct = Simplify[uAk . uBnegk];

(* Convert to spherical coordinates for k, with separation vector along z. *)
(* kx = k sin t cos p, ky = k sin t sin p, kz = k cos t                   *)
(* exp(i k·d) = exp(i k d cos t)                                          *)

integrand = dotProduct /. {
    kx -> k*Sin[t]*Cos[p],
    ky -> k*Sin[t]*Sin[p],
    kz -> k*Cos[t],
    k2 -> k^2,
    kn -> k
};

fullIntegrand = integrand * k^2 * Sin[t] * Exp[I*k*d*Cos[t]];

angInt = Integrate[
  fullIntegrand,
  {p, 0, 2*Pi}, {t, 0, Pi},
  Assumptions -> {d > 0, k > 0, Element[aHsym, Reals], Element[aLsym, Reals]},
  GenerateConditions -> False
];

radInt = Integrate[
  angInt,
  {k, 0, Infinity},
  Assumptions -> {d > 0, Element[aHsym, Reals], Element[aLsym, Reals]},
  GenerateConditions -> False
];

(* --------------------------------------------------------------------- *)
(* 3. Extract the isotropic cross-term tensor coefficients                *)
(* --------------------------------------------------------------------- *)
(* With separation along z, n = z-hat. Then (vA·n)(vB·n) = vAz vBz.        *)
(* The isotropic decomposition is:                                        *)
(*   I(d) = (const/d) * [ Cpara (vA·vB) + Clong (vA·n)(vB·n) ]             *)

CoeffVV = Coefficient[radInt, vAx*vBx];
CoeffZZ = Coefficient[radInt, vAz*vBz];

(* Normalize out the expected 1/d scaling. *)
VecParaShape = Simplify[CoeffVV * d];
VecLongShape = Simplify[(CoeffZZ - CoeffVV) * d];

Print["\nDerived geometric shapes (up to overall coupling K):"];
Print["  VecParaShape  (vA·vB)      = ", VecParaShape];
Print["  VecLongShape  ((vA·n)(vB·n)) = ", VecLongShape];

RatioModel = Simplify[VecLongShape/VecParaShape];
Print["\nModel ratio (Long/Para) = ", RatioModel];

(* --------------------------------------------------------------------- *)
(* 4. Solve for real wake parameters that match the EIH ratio             *)
(* --------------------------------------------------------------------- *)
RatioTarget = TargetLong/TargetPara;
Print["Target ratio (Long/Para) = ", RatioTarget];

(* We have two free real parameters (aHsym, aLsym) and one ratio equation. *)
(* Choose the *simplest* real solution by minimizing aH^2 + aL^2 subject  *)
(* to the ratio constraint. This avoids arbitrary hand-picking.           *)

Print["\n... Solving ratio constraint and minimizing wake complexity ..."];

minSol = Quiet @ Check[
  Minimize[
    {
      aHsym^2 + aLsym^2,
      RatioModel == RatioTarget,
      Element[aHsym, Reals], Element[aLsym, Reals]
    },
    {aHsym, aLsym}
  ],
  $Failed
];

If[minSol === $Failed || Head[minSol] =!= List,
  Print["Minimize did not return a symbolic solution. Falling back to FindInstance..."];
  inst = Quiet @ FindInstance[
    {
      RatioModel == RatioTarget,
      Element[aHsym, Reals], Element[aLsym, Reals]
    },
    {aHsym, aLsym},
    Reals,
    1
  ];
  If[inst === {},
    Print["\nFAILURE: No real solution found within this wake basis."];
    Print["  -> Add additional response channels (e.g., k-dependent weights A(k),B(k)) or revise the wake ansatz."];
    Abort[];
  ];
  params = inst[[1]];
  Print["\nOne real solution instance:"];
  Print[params];
,
  Print["\nMinimize result {minValue, {aHsym->..., aLsym->...}}:"];
  Print[minSol];
  params = Last[minSol];
];
aHval = aHsym /. params;
aLval = aLsym /. params;

Print["\nSelected real parameters:"];
Print["  aT (fixed) = ", aT];
Print["  aH = ", aHval];
Print["  aL = ", aLval];

(* Determine the overall coupling constant K that matches the parallel coefficient. *)
Kval = Simplify[TargetPara / (VecParaShape /. params)];

Print["\nOverall coupling K needed to match TargetPara:"];
Print["  K = ", Kval];

(* Predict the longitudinal coefficient with this K (should match TargetLong). *)
LongPred = Simplify[Kval * (VecLongShape /. params)];
ParaPred = Simplify[Kval * (VecParaShape /. params)];

Print["\nCheck (should reproduce EIH cross coefficients):"];
Print["  Predicted Para = ", ParaPred, "   (target ", TargetPara, ")"];
Print["  Predicted Long = ", LongPred, "   (target ", TargetLong, ")"];

(* --------------------------------------------------------------------- *)
(* 5. Stability sanity check                                              *)
(* --------------------------------------------------------------------- *)
(* The underlying quadratic functional is positive if the mode amplitudes *)
(* are real and the medium response weights for PT and PL are positive.   *)
(* Here we at least confirm real-valued mixing parameters.               *)

If[Element[aHval, Reals] && Element[aLval, Reals],
  Print["\nSUCCESS: Found a real-parameter wake that matches EIH cross terms."],
  Print["\nFAILURE: No real-parameter solution found by this criterion."]
];

Print["\nNOTE:"];
Print["  - If Minimize returns no real solution, the restricted basis is still too small,"];
Print["    and you must add additional response channels (e.g., k-dependent weights A(k))."];
Print["  - If a real solution exists, Paper III should treat these as *response parameters*"];
Print["    of the stiff vacuum (determined by EOS / microphysics), not as an imaginary alpha." ];


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
  VecParaShape  (vA·vB)      = (-1 + aHsym^2 - aLsym^2)*Pi^2
  VecLongShape  ((vA·n)(vB·n)) = (-1 + aHsym^2 + aLsym^2)*Pi^2

Model ratio (Long/Para) = (-1 + aHsym^2 + aLsym^2)/(-1 + aHsym^2 - aLsym^2)
Target ratio (Long/Para) = 1/7

... Solving ratio constraint and minimizing wake complexity ...

Minimize result {minValue, {aHsym->..., aLsym->...}}:
{3/4, {aHsym -> 0, aLsym -> -1/2*Sqrt[3]}}

Selected real parameters:
  aT (fixed) = 1
  aH = 0
  aL = -1/2*Sqrt[3]

Overall coupling K needed to match TargetPara:
  K = 2/Pi^2

Check (should reproduce EIH cross coefficients):
  Predicted Para = -7/2   (target -7/2)
  Predicted Long = -1/2   (target -1/2)

SUCCESS: Found a real-parameter wake that matches EIH cross terms.
"*)
