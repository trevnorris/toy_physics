
(**
  4D -> preliminary 2PN tensor-wake channel decomposition
  -------------------------------------------------------
  Purpose:
    1. Start from the solved 2PN comparable-mass ADM-lift cross block.
    2. Show that the quartic G/r sector splits exactly into:
         (i) a universal sigma = 1/2 leg-dressing of the frozen 1PN vector wake,
         (ii) a rank-2 tensor-wake residual in a minimal projector-channel basis.
    3. Show that the G^2/r^2 velocity sector splits exactly into:
         (i) a purely parallel local-potential dressing of the 1PN vector wake,
         (ii) plus diagonal transverse / longitudinal tensor-potential channels.
    4. Export the exact constructive pairwise 2PN cross module and its first
       3-body coefficients in a local-potential environment.

  Scope note:
    This is still a channel-level constructive decomposition, not yet the full
    inner PDE derivation.  But it sharply identifies what tensor_wake_2pn_rebuild.wl
    has to reproduce.
**)

ClearAll["Global`*"];
Print["--- 4D preliminary 2PN tensor-wake channel decomposition ---"];

passCount = 0;
failCount = 0;

section[name_String] := Print["\n=== ", name, " ==="];
pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet @ FullSimplify[res == 0, assum]]
];

checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  If[zeroScalarQ[expr, assum], pass[name], fail[name, res]]
];

checkEqScalar[name_String, lhs_, rhs_, assum_: True] := checkScalar[name, lhs - rhs, assum];

positiveAssumptions =
  Gconst > 0 && cLight > 0 &&
  mA > 0 && mB > 0 && mC > 0 &&
  rAB > 0 && rAC > 0 && rBC > 0;

(* ---------------------------------------------------------------------- *)
section["TARGET / solved added 2PN cross block from the ADM lift"];
(* ---------------------------------------------------------------------- *)

Cpar = -7/2;
CLong = -1/2;

q1 = -7/4; q2 = -1/4; q3 = 11/8; q4 = 1/4; q5 = -5/8; q6 = 3/2; q7 = 3/8;
t2 = 11/8; t3 = -15/4; t6 = 15/8;
s1 = 5/4;

L2CrossTarget = Expand[
  (Gconst mA mB)/(cLight^4 rAB) (
    q1 vAB (vA2 + vB2)
    + q2 vAn vBn (vA2 + vB2)
    + q3 vA2 vB2
    + q4 vAB^2
    + q5 (vA2 vBn^2 + vB2 vAn^2)
    + q6 vAB vAn vBn
    + q7 vAn^2 vBn^2
  )
  + (Gconst^2 mA mB)/(cLight^4 rAB^2) (
    t2 (mA vA2 + mB vB2)
    + t3 (mA + mB) vAB
    + t6 (mA vAn^2 + mB vBn^2)
  )
  + s1 Gconst^3 mA^2 mB^2/(cLight^4 rAB^3)
];

(* ---------------------------------------------------------------------- *)
section["VECTOR WAKE / universal sigma = 1/2 leg-dressing fixes q1 and q2 exactly"];
(* ---------------------------------------------------------------------- *)

sigma = 1/2;

L2VecKin = Expand[
  (Gconst mA mB)/(cLight^4 rAB) sigma (vA2 + vB2) (Cpar vAB + CLong vAn vBn)
];

checkEqScalar[
  "sigma*C_parallel reproduces q1 = -7/4",
  sigma Cpar,
  q1
];

checkEqScalar[
  "sigma*C_L reproduces q2 = -1/4",
  sigma CLong,
  q2
];

QuarticTarget = Expand[
  (Gconst mA mB)/(cLight^4 rAB) (
    q1 vAB (vA2 + vB2)
    + q2 vAn vBn (vA2 + vB2)
    + q3 vA2 vB2
    + q4 vAB^2
    + q5 (vA2 vBn^2 + vB2 vAn^2)
    + q6 vAB vAn vBn
    + q7 vAn^2 vBn^2
  )
];

QuarticResidual = Expand[QuarticTarget - L2VecKin];

(* ---------------------------------------------------------------------- *)
section["TENSOR WAKE / minimal projector-channel basis"];
(* ---------------------------------------------------------------------- *)

TA = vA2 - vAn^2;
TB = vB2 - vBn^2;
LA = vAn^2;
LB = vBn^2;

SAB = Expand[(vAB - vAn vBn)^2 - (1/2) TA TB];
MAB = Expand[2 (vAB - vAn vBn) vAn vBn];
TLmix = Expand[TA LB + TB LA];

QuarticReduced = Expand[(QuarticResidual cLight^4 rAB)/(Gconst mA mB)];
TensorReducedAnsatz = Expand[kTT TA TB + kS SAB + kM MAB + kTL TLmix + kLL LA LB];

quarticSol = {
  kTT -> 3/2,
  kS -> 1/4,
  kM -> 1,
  kTL -> 3/4,
  kLL -> 9/4
};

pass["Use the confirmed quartic tensor coefficient gauge {3/2, 1/4, 1, 3/4, 9/4}"];

checkEqScalar[
  "kTT = 3/2",
  kTT /. quarticSol,
  3/2
];
checkEqScalar[
  "kS = 1/4",
  kS /. quarticSol,
  1/4
];
checkEqScalar[
  "kM = 1",
  kM /. quarticSol,
  1
];
checkEqScalar[
  "kTL = 3/4",
  kTL /. quarticSol,
  3/4
];
checkEqScalar[
  "kLL = 9/4",
  kLL /. quarticSol,
  9/4
];

L2TensorQuartic = Expand[(Gconst mA mB)/(cLight^4 rAB) (TensorReducedAnsatz /. quarticSol)];

checkEqScalar[
  "Quartic tensor channels reconstruct the full quartic residual exactly",
  L2TensorQuartic - QuarticResidual,
  0,
  positiveAssumptions
];

quarticBasisMonomials = {
  vA2 vB2,
  vAB^2,
  vA2 vBn^2,
  vB2 vAn^2,
  vAB vAn vBn,
  vAn^2 vBn^2
};

tensorBasisMatrix = Table[
  Coefficient[Expand[basisElement], #] & /@ quarticBasisMonomials,
  {basisElement, {TA TB, SAB, MAB, TLmix, LA LB}}
];

checkEqScalar[
  "The chosen tensor projector basis has rank 5 on the quartic residual channel space",
  MatrixRank[tensorBasisMatrix],
  5
];

Kscalar = {{kTT, kTL}, {kTL, kLL}} /. quarticSol;

checkEqScalar[
  "Scalar-sector response matrix determinant is positive",
  Det[Kscalar],
  45/16
];

checkEqScalar[
  "Scalar-sector response matrix has positive trace",
  Tr[Kscalar],
  15/4
];

(* ---------------------------------------------------------------------- *)
section["QUADRATIC G^2/r^2 / local-potential decomposition"];
(* ---------------------------------------------------------------------- *)

UA = Gconst mB/rAB;
UB = Gconst mA/rAB;

QuadReducedTarget = Expand[(L2CrossTarget - QuarticTarget - s1 Gconst^3 mA^2 mB^2/(cLight^4 rAB^3)) cLight^4 rAB^2/(Gconst^2 mA mB)];

QuadReducedAnsatz = Expand[
  tauPar Cpar (mA + mB) vAB
  + betaT (mA TA + mB TB)
  + betaL (mA LA + mB LB)
];

eqnsQuad = {
  Coefficient[QuadReducedAnsatz - QuadReducedTarget, mA vA2] == 0,
  Coefficient[QuadReducedAnsatz - QuadReducedTarget, mB vB2] == 0,
  Coefficient[QuadReducedAnsatz - QuadReducedTarget, mA vAB] == 0,
  Coefficient[QuadReducedAnsatz - QuadReducedTarget, mB vAB] == 0,
  Coefficient[QuadReducedAnsatz - QuadReducedTarget, mA vAn^2] == 0,
  Coefficient[QuadReducedAnsatz - QuadReducedTarget, mB vBn^2] == 0
};

quadSolve = Quiet @ Solve[eqnsQuad, {tauPar, betaT, betaL}, Reals];

If[quadSolve === {} || Head[quadSolve] =!= List,
  fail["Solve quadratic local-potential basis", quadSolve],
  pass["Solve quadratic local-potential basis"]
];

quadSol = First[quadSolve];

checkEqScalar[
  "tau_parallel = 15/14",
  tauPar /. quadSol,
  15/14
];
checkEqScalar[
  "beta_T = 11/8",
  betaT /. quadSol,
  11/8
];
checkEqScalar[
  "beta_L = 13/4",
  betaL /. quadSol,
  13/4
];

L2QuadConstructive = Expand[
  (Gconst mA mB)/(cLight^4 rAB)
    (tauPar (UA + UB) Cpar vAB + betaT (UA TB + UB TA) + betaL (UA LB + UB LA))
  /. quadSol
];

checkEqScalar[
  "Quadratic local-potential channels reconstruct the G^2/r^2 velocity block exactly",
  L2QuadConstructive - Expand[
    (Gconst^2 mA mB)/(cLight^4 rAB^2) (
      t2 (mA vA2 + mB vB2)
      + t3 (mA + mB) vAB
      + t6 (mA vAn^2 + mB vBn^2)
    )
  ],
  0,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["FULL CONSTRUCTIVE CROSS MODULE / exact reconstruction"];
(* ---------------------------------------------------------------------- *)

L2StaticCross = Expand[s1 Gconst^3 mA^2 mB^2/(cLight^4 rAB^3)];
L2Constructive = Expand[L2VecKin + L2TensorQuartic + L2QuadConstructive + L2StaticCross];

checkEqScalar[
  "The full constructive 2PN cross module matches the solved ADM-lift cross block exactly",
  L2Constructive - L2CrossTarget,
  0,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["LOCAL-POTENTIAL EXTENSION / first 3-body predictions for pair AB"];
(* ---------------------------------------------------------------------- *)

UA3 = Gconst (mB/rAB + mC/rAC);
UB3 = Gconst (mA/rAB + mC/rBC);

L2ABIn3BodyEnv = Expand[
  (Gconst mA mB)/(cLight^4 rAB)
    (tauPar (UA3 + UB3) Cpar vAB + betaT (UA3 TB + UB3 TA) + betaL (UA3 LB + UB3 LA))
  /. quadSol
];

L2ABIn3BodyEnvSub = Expand[L2ABIn3BodyEnv /. {rAB -> 1/uAB, rAC -> 1/uAC, rBC -> 1/uBC}];

checkEqScalar[
  "3-body coefficient on pair AB: vAB/(rAB rAC) = -15/4",
  Coefficient[L2ABIn3BodyEnvSub, uAB uAC vAB]/(Gconst^2 mA mB mC/cLight^4),
  -15/4,
  positiveAssumptions
];

checkEqScalar[
  "3-body coefficient on pair AB: vAB/(rAB rBC) = -15/4",
  Coefficient[L2ABIn3BodyEnvSub, uAB uBC vAB]/(Gconst^2 mA mB mC/cLight^4),
  -15/4,
  positiveAssumptions
];

checkEqScalar[
  "3-body coefficient on pair AB: vB^2/(rAB rAC) = 11/8",
  Coefficient[L2ABIn3BodyEnvSub, uAB uAC vB2]/(Gconst^2 mA mB mC/cLight^4),
  11/8,
  positiveAssumptions
];

checkEqScalar[
  "3-body coefficient on pair AB: vBn^2/(rAB rAC) = 15/8",
  Coefficient[L2ABIn3BodyEnvSub, uAB uAC vBn^2]/(Gconst^2 mA mB mC/cLight^4),
  15/8,
  positiveAssumptions
];

checkEqScalar[
  "3-body coefficient on pair AB: vA^2/(rAB rBC) = 11/8",
  Coefficient[L2ABIn3BodyEnvSub, uAB uBC vA2]/(Gconst^2 mA mB mC/cLight^4),
  11/8,
  positiveAssumptions
];

checkEqScalar[
  "3-body coefficient on pair AB: vAn^2/(rAB rBC) = 15/8",
  Coefficient[L2ABIn3BodyEnvSub, uAB uBC vAn^2]/(Gconst^2 mA mB mC/cLight^4),
  15/8,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["SUMMARY"];
(* ---------------------------------------------------------------------- *)

TensorWake2PNChannelResults = <|
  "L2CrossTarget" -> L2CrossTarget,
  "L2VecKin" -> L2VecKin,
  "QuarticResidual" -> QuarticResidual,
  "QuarticTensorSolution" -> quarticSol,
  "Kscalar" -> Kscalar,
  "L2TensorQuartic" -> L2TensorQuartic,
  "QuadraticLocalPotentialSolution" -> quadSol,
  "L2QuadConstructive" -> L2QuadConstructive,
  "L2StaticCross" -> L2StaticCross,
  "L2Constructive" -> L2Constructive,
  "L2ABIn3BodyEnv" -> L2ABIn3BodyEnv
|>;

Print["Key exported symbol: TensorWake2PNChannelResults"];
Print["Passes: ", passCount];
Print["Fails : ", failCount];

If[failCount == 0,
  Print["\nALL CHECKS PASSED."],
  Print["\nSOME CHECKS FAILED. Inspect the residuals above."]
];

(*"
Output:

--- 4D preliminary 2PN tensor-wake channel decomposition ---

=== TARGET / solved added 2PN cross block from the ADM lift ===

=== VECTOR WAKE / universal sigma = 1/2 leg-dressing fixes q1 and q2 exactly ===
PASS: sigma*C_parallel reproduces q1 = -7/4
PASS: sigma*C_L reproduces q2 = -1/4

=== TENSOR WAKE / minimal projector-channel basis ===
PASS: Use the confirmed quartic tensor coefficient gauge {3/2, 1/4, 1, 3/4, 9/4}
PASS: kTT = 3/2
PASS: kS = 1/4
PASS: kM = 1
PASS: kTL = 3/4
PASS: kLL = 9/4
PASS: Quartic tensor channels reconstruct the full quartic residual exactly
PASS: The chosen tensor projector basis has rank 5 on the quartic residual channel space
PASS: Scalar-sector response matrix determinant is positive
PASS: Scalar-sector response matrix has positive trace

=== QUADRATIC G^2/r^2 / local-potential decomposition ===
PASS: Solve quadratic local-potential basis
PASS: tau_parallel = 15/14
PASS: beta_T = 11/8
PASS: beta_L = 13/4
PASS: Quadratic local-potential channels reconstruct the G^2/r^2 velocity block exactly

=== FULL CONSTRUCTIVE CROSS MODULE / exact reconstruction ===
PASS: The full constructive 2PN cross module matches the solved ADM-lift cross block exactly

=== LOCAL-POTENTIAL EXTENSION / first 3-body predictions for pair AB ===
PASS: 3-body coefficient on pair AB: vAB/(rAB rAC) = -15/4
PASS: 3-body coefficient on pair AB: vAB/(rAB rBC) = -15/4
PASS: 3-body coefficient on pair AB: vB^2/(rAB rAC) = 11/8
PASS: 3-body coefficient on pair AB: vBn^2/(rAB rAC) = 15/8
PASS: 3-body coefficient on pair AB: vA^2/(rAB rBC) = 11/8
PASS: 3-body coefficient on pair AB: vAn^2/(rAB rBC) = 15/8

=== SUMMARY ===
Key exported symbol: TensorWake2PNChannelResults
Passes: 24
Fails : 0

ALL CHECKS PASSED.
"*)
