
(**
  4D -> preliminary 2PN comparable-mass ADM lift
  ----------------------------------------------
  Purpose:
    Continue the 2PN program past the one-body DtN closure and determine the
    minimal additional comparable-mass cross sector needed to match the
    standard generic-frame ADM 2PN Hamiltonian.

  Strategy:
    1) Freeze the already-passed sectors:
         - Newtonian + full frozen 1PN Lagrangian,
         - DtN-corrected 2PN one-body/self sector,
         - quadratic local mass scaling with lambda_rho = 1/2.
    2) Use the generic-frame ADM Hamiltonians H_1PN and H_2PN from Appendix C
       (Eqs. C.7-C.8) of the Schäfer-Jaranowski 2018 review as the target.
    3) Use the perturbative Legendre transform for
          L = L0 + eps L1 + eps^2 L2
       with quadratic L0:
          H1 = -L1(v0),
          H2 = -L2(v0) + (1/2) A0^T M^-1 A0,
       where A0 = (∂L1/∂v)|_{v=v0}, v0 = p/m.
    4) Compute the H_2PN residual of the DtN-corrected self/static candidate.
    5) Solve that residual against a compact invariant 2PN cross-sector basis.
    6) Rebuild the full candidate and verify exact H_2PN ADM matching.

  Scope note:
    This is a target-matching / bookkeeping result.  It does not yet supply the
    constructive 2PN wake derivation that would explain the solved coefficients
    from the inner-throat geometry.
**)

ClearAll["Global`*"];
Print["--- 4D preliminary 2PN comparable-mass ADM lift ---"];

passCount = 0;
failCount = 0;
section[name_String] := Print["\n=== ", name, " ==="];
pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);
info[msg_String] := Print["INFO: ", msg];

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet @ FullSimplify[res == 0, assum]]
];
checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  If[zeroScalarQ[expr, assum], pass[name], fail[name, res]]
];
checkEqScalar[name_String, lhs_, rhs_, assum_: True] := checkScalar[name, lhs - rhs, assum];
checkCondition[name_String, cond_] := If[TrueQ[cond], pass[name], fail[name, cond]];

positiveAssumptions = Gconst > 0 && rAB > 0 && mA > 0 && mB > 0;

(* ---------------------------------------------------------------------- *)
section["SETUP / kinematic variables and frozen Lagrangian sectors"];
(* ---------------------------------------------------------------------- *)

vA = {vAx, vAy, vAz};
vB = {vBx, vBy, vBz};
pA = {pAx, pAy, pAz};
pB = {pBx, pBy, pBz};

vA2 = Expand[vA.vA];
vB2 = Expand[vB.vB];
vAB = Expand[vA.vB];
vAn = vAz;
vBn = vBz;

pA2 = Expand[pA.pA];
pB2 = Expand[pB.pB];
pAB = Expand[pA.pB];
pAn = pAz;
pBn = pBz;

v0Rules = {
  vAx -> pAx/mA, vAy -> pAy/mA, vAz -> pAz/mA,
  vBx -> pBx/mB, vBy -> pBy/mB, vBz -> pBz/mB
};

pToMvRules = {
  pAx -> mA vAx, pAy -> mA vAy, pAz -> mA vAz,
  pBx -> mB vBx, pBy -> mB vBy, pBz -> mB vBz
};

L0 = (1/2) mA vA2 + (1/2) mB vB2 + Gconst mA mB/rAB;

L1 = (
  (mA vA2^2 + mB vB2^2)/8
  + (Gconst mA mB/rAB) ((3/2) (vA2 + vB2) - (7/2) vAB - (1/2) vAn vBn)
  - (Gconst^2 mA mB (mA + mB))/(2 rAB^2)
);

(* DtN-corrected self/static 2PN block:
   free v^6 coefficient = 1/16
   self U v^4 coefficient = 7/8
   self U^2 v^2 coefficient = 2
   static U^3 coefficient = 1/4  (lambda_rho = 1/2)
*)
L2SelfStaticDtN = (
  (mA vA2^3 + mB vB2^3)/16
  + (7 Gconst mA mB/(8 rAB)) (vA2^2 + vB2^2)
  + (2 Gconst^2 mA mB/rAB^2) (mB vA2 + mA vB2)
  + (Gconst^3 mA mB (mA^2 + mB^2))/(4 rAB^3)
);

(* ---------------------------------------------------------------------- *)
section["TARGET / generic-frame ADM Hamiltonians through 2PN"];
(* ---------------------------------------------------------------------- *)

H0Model = Expand[pA2/(2 mA) + pB2/(2 mB) - Gconst mA mB/rAB];

A0 = Expand /@ ((D[L1, #] /. v0Rules) & /@ {vAx, vAy, vAz, vBx, vBy, vBz});
quadLegendre = Expand[
  (A0[[1]]^2 + A0[[2]]^2 + A0[[3]]^2)/(2 mA)
  + (A0[[4]]^2 + A0[[5]]^2 + A0[[6]]^2)/(2 mB)
];

H1Model = Expand[-(L1 /. v0Rules)];
H2ModelNoCross = Expand[-(L2SelfStaticDtN /. v0Rules) + quadLegendre];

H1Target = Expand[
  -(pA2^2)/(8 mA^3) - (pB2^2)/(8 mB^3)
  + (Gconst mA mB/(4 rAB)) (
      -6 pA2/mA^2
      - 6 pB2/mB^2
      + 14 pAB/(mA mB)
      + 2 pAn pBn/(mA mB)
    )
  + (Gconst^2 mA mB (mA + mB))/(2 rAB^2)
];

H2TargetBase = Expand[
  pA2^3/(16 mA^5)
  + (Gconst mA mB/(8 rAB)) (
      5 pA2^2/mA^4
      - (11/2) pA2 pB2/(mA^2 mB^2)
      - pAB^2/(mA^2 mB^2)
      + 5 pA2 pBn^2/(mA^2 mB^2)
      - 6 pAB pAn pBn/(mA^2 mB^2)
      - (3/2) pAn^2 pBn^2/(mA^2 mB^2)
    )
  + (Gconst^2 mA mB/(4 rAB^2)) (
      mB (10 pA2/mA^2 + 19 pB2/mB^2)
      - (1/2) (mA + mB) (27 pAB + 6 pAn pBn)/(mA mB)
    )
  - (Gconst^3 mA mB (mA^2 + 5 mA mB + mB^2))/(8 rAB^3)
];

swapRules = {
  mA -> mB, mB -> mA,
  pAx -> pBx, pAy -> pBy, pAz -> pBz,
  pBx -> pAx, pBy -> pAy, pBz -> pAz
};

H2Target = Expand[H2TargetBase + (H2TargetBase /. swapRules)];

checkEqScalar[
  "Frozen 1PN Lagrangian still Legendre-transforms to the generic-frame ADM H1PN target",
  H1Model,
  H1Target,
  positiveAssumptions
];

(* ---------------------------------------------------------------------- *)
section["CONTROLLED / residual of the DtN-corrected self/static 2PN candidate"];
(* ---------------------------------------------------------------------- *)

H2ResidualNoCross = Expand[H2ModelNoCross - H2Target];
info["H2ResidualNoCross is the exact 2PN Hamiltonian residual of the DtN-corrected self/static candidate against the generic-frame ADM target."];

checkCondition[
  "The no-cross comparable-mass 2PN residual is nonzero (so a genuine 2PN cross sector is still needed)",
  Not[zeroScalarQ[H2ResidualNoCross, positiveAssumptions]]
];

L2CrossRequiredMapped = Expand[H2ResidualNoCross /. pToMvRules];
info["Because H2 = -L2(v0) + (1/2) A0^T M^-1 A0, any added 2PN cross block enters the Hamiltonian only as -L2Cross(v0); therefore the required added Lagrangian block is H2ResidualNoCross with p -> m v."];

(* ---------------------------------------------------------------------- *)
section["SOLVE / compact invariant 2PN cross-sector basis"];
(* ---------------------------------------------------------------------- *)

quarticBasis = {
  mA mB vAB (vA2 + vB2),
  mA mB vAn vBn (vA2 + vB2),
  mA mB vA2 vB2,
  mA mB vAB^2,
  mA mB (vA2 vBn^2 + vB2 vAn^2),
  mA mB vAB vAn vBn,
  mA mB vAn^2 vBn^2
};

quadraticBasis = {
  mA mB (mB vA2 + mA vB2),
  mA mB (mA vA2 + mB vB2),
  mA mB (mA + mB) vAB,
  mA mB (mA + mB) vAn vBn,
  mA mB (mB vAn^2 + mA vBn^2),
  mA mB (mA vAn^2 + mB vBn^2)
};

staticBasis = mA^2 mB^2;

crossCoeffVars = {q1, q2, q3, q4, q5, q6, q7, t1, t2, t3, t4, t5, t6, s1};

L2CrossAnsatz = Expand[
  (Gconst/rAB) (
    q1 quarticBasis[[1]] + q2 quarticBasis[[2]] + q3 quarticBasis[[3]]
    + q4 quarticBasis[[4]] + q5 quarticBasis[[5]] + q6 quarticBasis[[6]]
    + q7 quarticBasis[[7]]
  )
  + (Gconst^2/rAB^2) (
    t1 quadraticBasis[[1]] + t2 quadraticBasis[[2]] + t3 quadraticBasis[[3]]
    + t4 quadraticBasis[[4]] + t5 quadraticBasis[[5]] + t6 quadraticBasis[[6]]
  )
  + (Gconst^3/rAB^3) s1 staticBasis
];

basisVars = {vAx, vAy, vAz, vBx, vBy, vBz};
crossSolveRules = Module[{coeffs, sol},
  coeffs = Values @ CoefficientRules[Expand[L2CrossRequiredMapped - L2CrossAnsatz], basisVars];
  sol = Quiet @ Solve[Thread[coeffs == 0], crossCoeffVars];
  If[sol === {} || Head[sol] =!= List, $Failed, First[sol]]
];

checkCondition[
  "The chosen 7+6+1 invariant cross basis admits a coefficient solution",
  crossSolveRules =!= $Failed
];

If[crossSolveRules =!= $Failed,
  checkEqScalar["q1 = -7/4", q1 /. crossSolveRules, -7/4, positiveAssumptions];
  checkEqScalar["q2 = -1/4", q2 /. crossSolveRules, -1/4, positiveAssumptions];
  checkEqScalar["q3 = 11/8", q3 /. crossSolveRules, 11/8, positiveAssumptions];
  checkEqScalar["q4 = 1/4", q4 /. crossSolveRules, 1/4, positiveAssumptions];
  checkEqScalar["q5 = -5/8", q5 /. crossSolveRules, -5/8, positiveAssumptions];
  checkEqScalar["q6 = 3/2", q6 /. crossSolveRules, 3/2, positiveAssumptions];
  checkEqScalar["q7 = 3/8", q7 /. crossSolveRules, 3/8, positiveAssumptions];

  checkEqScalar["t1 = 0", t1 /. crossSolveRules, 0, positiveAssumptions];
  checkEqScalar["t2 = 11/8", t2 /. crossSolveRules, 11/8, positiveAssumptions];
  checkEqScalar["t3 = -15/4", t3 /. crossSolveRules, -15/4, positiveAssumptions];
  checkEqScalar["t4 = 0", t4 /. crossSolveRules, 0, positiveAssumptions];
  checkEqScalar["t5 = 0", t5 /. crossSolveRules, 0, positiveAssumptions];
  checkEqScalar["t6 = 15/8", t6 /. crossSolveRules, 15/8, positiveAssumptions];
  checkEqScalar["s1 = 5/4", s1 /. crossSolveRules, 5/4, positiveAssumptions];

  L2CrossAdded2PN = Expand[L2CrossAnsatz /. crossSolveRules];

  checkEqScalar[
    "The solved added cross block equals the required mapped residual exactly",
    L2CrossAdded2PN,
    L2CrossRequiredMapped,
    positiveAssumptions
  ];
];

(* ---------------------------------------------------------------------- *)
section["FULL / exact generic-frame ADM H2PN match after adding the solved cross block"];
(* ---------------------------------------------------------------------- *)

If[crossSolveRules =!= $Failed,
  L2Full2PN = Expand[L2SelfStaticDtN + L2CrossAdded2PN];
  H2FullCandidate = Expand[-(L2Full2PN /. v0Rules) + quadLegendre];

  checkEqScalar[
    "Full DtN-corrected 2PN candidate now Legendre-transforms to the generic-frame ADM H2PN target exactly",
    H2FullCandidate,
    H2Target,
    positiveAssumptions
  ];

  L2CrossAdded2PNClean = Expand[
    (Gconst mA mB/rAB) (
      -(7/4) vAB (vA2 + vB2)
      -(1/4) vAn vBn (vA2 + vB2)
      +(11/8) vA2 vB2
      +(1/4) vAB^2
      -(5/8) (vA2 vBn^2 + vB2 vAn^2)
      +(3/2) vAB vAn vBn
      +(3/8) vAn^2 vBn^2
    )
    + (Gconst^2 mA mB/rAB^2) (
      +(11/8) (mA vA2 + mB vB2)
      -(15/4) (mA + mB) vAB
      +(15/8) (mA vAn^2 + mB vBn^2)
    )
    + (5 Gconst^3 mA^2 mB^2)/(4 rAB^3)
  ];

  checkEqScalar[
    "Clean invariant formula reproduces the solved cross block",
    L2CrossAdded2PNClean,
    L2CrossAdded2PN,
    positiveAssumptions
  ];

  L2CrossPerUnitTestMass = FullSimplify[
    Limit[(L2CrossAdded2PNClean /. {vBx -> 0, vBy -> 0, vBz -> 0})/mA, mA -> 0],
    Gconst > 0 && rAB > 0 && mB > 0
  ];

  checkEqScalar[
    "Added comparable-mass cross block vanishes in the strict test-mass limit with the heavy body at rest",
    L2CrossPerUnitTestMass,
    0,
    Gconst > 0 && rAB > 0 && mB > 0
  ];

  L2CombinedStatic = FullSimplify[
    Coefficient[Expand[L2Full2PN /. {vAx -> 0, vAy -> 0, vAz -> 0, vBx -> 0, vBy -> 0, vBz -> 0}], Gconst^3/rAB^3],
    positiveAssumptions
  ];

  checkEqScalar[
    "The combined 2PN static mass polynomial becomes +(mA mB/4)(mA^2 + 5 mA mB + mB^2)",
    L2CombinedStatic,
    (mA mB (mA^2 + 5 mA mB + mB^2))/4,
    positiveAssumptions
  ];
];

(* ---------------------------------------------------------------------- *)
section["SUMMARY"];
(* ---------------------------------------------------------------------- *)

If[crossSolveRules =!= $Failed,
  TwoPNComparableMassADMLiftResults = <|
    "H2ResidualNoCross" -> H2ResidualNoCross,
    "CrossBasisSolution2PN" -> crossSolveRules,
    "L2CrossAdded2PN" -> L2CrossAdded2PN,
    "L2CrossAdded2PNClean" -> L2CrossAdded2PNClean,
    "L2Full2PN" -> L2Full2PN,
    "H2TargetADM" -> H2Target,
    "H2FullCandidate" -> H2FullCandidate
  |>;
];

Print["Key exported symbols: H2ResidualNoCross, CrossBasisSolution2PN, L2CrossAdded2PN, L2CrossAdded2PNClean, L2Full2PN, H2TargetADM, H2FullCandidate, TwoPNComparableMassADMLiftResults."];
Print["Passes: ", passCount];
Print["Fails : ", failCount];

If[failCount == 0,
  Print["\nALL CHECKS PASSED."],
  Print["\nSOME CHECKS FAILED."]
];

(*"
Output:

--- 4D preliminary 2PN comparable-mass ADM lift ---

=== SETUP / kinematic variables and frozen Lagrangian sectors ===

=== TARGET / generic-frame ADM Hamiltonians through 2PN ===
PASS: Frozen 1PN Lagrangian still Legendre-transforms to the generic-frame ADM H1PN target

=== CONTROLLED / residual of the DtN-corrected self/static 2PN candidate ===
INFO: H2ResidualNoCross is the exact 2PN Hamiltonian residual of the DtN-corrected self/static candidate against the generic-frame ADM target.
PASS: The no-cross comparable-mass 2PN residual is nonzero (so a genuine 2PN cross sector is still needed)
INFO: Because H2 = -L2(v0) + (1/2) A0^T M^-1 A0, any added 2PN cross block enters the Hamiltonian only as -L2Cross(v0); therefore the required added Lagrangian block is H2ResidualNoCross with p -> m v.

=== SOLVE / compact invariant 2PN cross-sector basis ===
PASS: The chosen 7+6+1 invariant cross basis admits a coefficient solution
PASS: q1 = -7/4
PASS: q2 = -1/4
PASS: q3 = 11/8
PASS: q4 = 1/4
PASS: q5 = -5/8
PASS: q6 = 3/2
PASS: q7 = 3/8
PASS: t1 = 0
PASS: t2 = 11/8
PASS: t3 = -15/4
PASS: t4 = 0
PASS: t5 = 0
PASS: t6 = 15/8
PASS: s1 = 5/4
PASS: The solved added cross block equals the required mapped residual exactly

=== FULL / exact generic-frame ADM H2PN match after adding the solved cross block ===
PASS: Full DtN-corrected 2PN candidate now Legendre-transforms to the generic-frame ADM H2PN target exactly
PASS: Clean invariant formula reproduces the solved cross block
PASS: Added comparable-mass cross block vanishes in the strict test-mass limit with the heavy body at rest
PASS: The combined 2PN static mass polynomial becomes +(mA mB/4)(mA^2 + 5 mA mB + mB^2)

=== SUMMARY ===
Key exported symbols: H2ResidualNoCross, CrossBasisSolution2PN, L2CrossAdded2PN, L2CrossAdded2PNClean, L2Full2PN, H2TargetADM, H2FullCandidate, TwoPNComparableMassADMLiftResults.
Passes: 22
Fails : 0

ALL CHECKS PASSED.
"*)
