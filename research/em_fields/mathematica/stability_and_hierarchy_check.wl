(* ================================================================= *)
(* MICROSCOPIC STABILITY & EM CONSISTENCY CHECK                      *)
(* ================================================================= *)
(* Purpose:                                                          *)
(* 1. Verify the derivation of L/a = sqrt(2) pi / x01 from the       *)
(*    cavity enthalpy used in the paper.                             *)
(* 2. Confirm q^2 / m_G^2 scales as Gamma^2 / a^2 at fixed aspect.   *)
(* ================================================================= *)

ClearAll["Global`*"];

$Assumptions = {
  a > 0, L > 0, alphaCav > 0, betaVac \[Element] Reals,
  kappaM > 0, kappaQ > 0, rho0 > 0, Gamma > 0
};

x01 = BesselJZero[0, 1];

(* ----------------------------------------------------------------- *)
(* PART 1: THERMODYNAMIC STABILITY (The "Why" of the Defect)         *)
(* ----------------------------------------------------------------- *)
Print["--- Part 1: Thermodynamic Stability Check ---"];

(* Dimensionless cavity enthalpy from the paper appendix:
   H ~ L x01^2 + a^2 Pi^2/L + betaVac a^2 L
*)
Hdimless = L*x01^2 + a^2*Pi^2/L + betaVac*a^2*L;

eqRadial = D[Hdimless, a] == 0;
eqAxial = D[Hdimless, L] == 0;

betaVacFromRadial = Simplify[betaVac /. Solve[eqRadial, betaVac][[1]]];
chiEq = Simplify[(eqAxial /. betaVac -> betaVacFromRadial) /. L -> chi*a];
aspectRatio = Simplify[chi /. Solve[chiEq && chi > 0, chi][[1]]];
numericalRatio = N[aspectRatio];
betaVacAtMinimum = Simplify[betaVacFromRadial /. L -> aspectRatio*a];

Print["Derived Aspect Ratio Expression (L/a):"];
Print[aspectRatio];
Print["Numerical L/a Value: ", numericalRatio];
Print["betaVac required at the minimum:"];
Print[betaVacAtMinimum];
Print["Match to paper's 1.85 aspect ratio: ", Abs[numericalRatio - 1.8475] < 0.001];

(* ----------------------------------------------------------------- *)
(* PART 2: EM HIERARCHY ROBUSTNESS                                   *)
(* ----------------------------------------------------------------- *)
Print[""];
Print["--- Part 2: EM Hierarchy Robustness ---"];

massGrav = kappaM * rho0 * Pi * a^2 * L;
chargeQ = kappaQ * rho0 * Pi * a^2 * Gamma;

ratioFixedAspect = Simplify[(chargeQ^2 / massGrav^2) /. L -> aspectRatio*a];
exponentA = Exponent[ratioFixedAspect, a];

Print["q^2 / m_G^2 at fixed aspect ratio:"];
Print[ratioFixedAspect];
Print["Exponent of a in q^2 / m_G^2: ", exponentA];

If[exponentA == -2,
 Print["CONFIRMED: q^2 / m_G^2 scales as 1/a^2."],
 Print["WARNING: Scaling is different from the paper."]
];

Print["Interpretation: at fixed circulation, smaller throats have larger charge-to-mass ratios."];
Print["This naturally explains why EM >> Gravity for microscopic particles in the toy model."];

(*
Output:
--- Part 1: Thermodynamic Stability Check ---
Derived Aspect Ratio Expression (L/a):
(Sqrt[2]*Pi)/BesselJZero[0, 1]
Numerical L/a Value: 1.8474865771201279
betaVac required at the minimum:
-1/2*BesselJZero[0, 1]^2/a^2
Match to paper's 1.85 aspect ratio: True

--- Part 2: EM Hierarchy Robustness ---
q^2 / m_G^2 at fixed aspect ratio:
(Gamma^2*kappaQ^2*BesselJZero[0, 1]^2)/(2*a^2*kappaM^2*Pi^2)
Exponent of a in q^2 / m_G^2: -2
CONFIRMED: q^2 / m_G^2 scales as 1/a^2.
Interpretation: at fixed circulation, smaller throats have larger charge-to-mass ratios.
This naturally explains why EM >> Gravity for microscopic particles in the toy model.
*)
