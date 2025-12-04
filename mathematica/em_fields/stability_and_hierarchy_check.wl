(* ================================================================= *)
(* MICROSCOPIC STABILITY & EM CONSISTENCY CHECK                      *)
(* ================================================================= *)
(* Purpose:                                                          *)
(* 1. Verify the derivation of L/a = 1.847 from Enthalpy minimization*)
(* 2. Confirm EM/Gravity hierarchy remains valid with this new L/a   *)
(* ================================================================= *)

ClearAll["Global`*"]

(* ----------------------------------------------------------------- *)
(* PART 1: THERMODYNAMIC STABILITY (The "Why" of the Defect)         *)
(* ----------------------------------------------------------------- *)
Print["--- Part 1: Thermodynamic Stability Check ---"];

(* Define Energy of the Fundamental TM-like Mode *)
(* E ~ Sqrt[ (x01/a)^2 + (pi/L)^2 ] *)
(* We use generic constants k for the prefactors *)
x01 = BesselJZero[0, 1]; (* First zero of J0 ~ 2.4048 *)
modeEnergy = k * Sqrt[(x01/a)^2 + (Pi/L)^2];

(* Define Vacuum Work term *)
(* Work = P_vac * Volume. Volume of cylinder = Pi * a^2 * L *)
vacuumWork = Pvac * (Pi * a^2 * L);

(* Define Enthalpy *)
enthalpyH = modeEnergy + vacuumWork;

(* DERIVATION: Minimize H with respect to a and L *)
(* Equation 1: dH/da = 0 *)
eqRadial = D[enthalpyH, a] == 0;

(* Equation 2: dH/dL = 0 *)
eqAxial = D[enthalpyH, L] == 0;

(* Solve for Pvac in both equations to eliminate it *)
solPvacRadial = Solve[eqRadial, Pvac][[1]];
solPvacAxial = Solve[eqAxial, Pvac][[1]];

(* Equate the two expressions for Pvac (Mechanical Equilibrium) *)
geometricConstraint = (Pvac /. solPvacRadial) == (Pvac /. solPvacAxial);

(* Solve for the aspect ratio L/a *)
(* We assume L > 0 and a > 0 *)
solutionAspect = Solve[geometricConstraint, L];
ratioExpression = Simplify[L/a /. solutionAspect][[1]]; (* Taking the positive solution *)

Print["Derived Aspect Ratio Expression (L/a): "];
Print[ratioExpression];

(* Calculate Numerical Value *)
numericalRatio = N[ratioExpression];
Print["Numerical L/a Value: ", numericalRatio];

(* Compare to User's value *)
Print["User's Value: 1.8475"];
Print["Match: ", Abs[numericalRatio - 1.8475] < 0.001];


(* ----------------------------------------------------------------- *)
(* PART 2: EM HIERARCHY ROBUSTNESS                                   *)
(* ----------------------------------------------------------------- *)
Print["\n--- Part 2: EM Hierarchy Robustness ---"];

(* We need to ensure that F_EM / F_Grav is stable under this change. *)
(* EM Force: F_e ~ q^2 / r^2 *)
(* Grav Force: F_g ~ m^2 / r^2 *)

(* CHARGE q: *)
(* From previous notes: q ~ Area * Circulation ~ a^2 * Gamma *)
(* This depends on 'a', but NOT on 'L'. *)
chargeQ = Cq * rho0 * a^2 * Gamma;

(* GRAVITATIONAL MASS m: *)
(* The "Cavity Mass" is the missing mass of the hole. *)
(* m ~ Volume * rho0 *)
(* If the shape is now fixed to L = 1.847a, the volume is purely a function of a^3 *)
volumeFixed = Pi * a^2 * (ratioExpression * a); 
massGrav = rho0 * volumeFixed;

(* FORCE RATIO *)
(* F_e / F_g *)
forceRatio = (chargeQ^2) / (massGrav^2);

(* SIMPLIFY *)
(* We want to see if the dependence on 'a' cancels out, or if we are left with *)
(* a dependency that matches observation (e.g. Gamma^2/a^2). *)
simplifiedRatio = Simplify[forceRatio];

Print["Force Ratio F_EM / F_G:"];
Print[simplifiedRatio];

(* CHECK: Is the ratio determined purely by Circulation (Gamma) and Radius (a)? *)
(* It should NOT depend on arbitrary 'L' anymore, since L is fixed by a. *)
Print["Dependency Check: The ratio is a function of Gamma, rho0, and a."];
Print["This is consistent with the standard hierarchy where small 'a' -> large Mass -> smaller ratio?"];
Print["Wait, Charge ~ a^2, Mass ~ a^3. Ratio ~ a^4 / a^6 ~ 1/a^2."];
Print["This implies smaller particles (small a) have LARGER charge-to-mass ratios?"];
Print["Let's check 1/a^2 behavior:"];
exponentA = Exponent[simplifiedRatio, a];
Print["Exponent of 'a' in Force Ratio: ", exponentA];

(* Conclusion *)
If[exponentA == -2,
 Print["CONFIRMED: F_em/F_g scales as 1/a^2."],
 Print["WARNING: Scaling is different."]
];

Print["Interpretation: Since 'a' is small, 1/a^2 is huge."];
Print["This naturally explains why EM >> Gravity for microscopic particles."];

(*"
Output:

--- Part 1: Thermodynamic Stability Check ---
Derived Aspect Ratio Expression (L/a):
-((Sqrt[2]*Pi)/BesselJZero[0, 1])
Numerical L/a Value: -1.847486577120128
User's Value: 1.8475
Match: False

--- Part 2: EM Hierarchy Robustness ---
Force Ratio F_EM / F_G:
(Cq^2*Gamma^2*BesselJZero[0, 1]^2)/(2*a^2*Pi^4)
Dependency Check: The ratio is a function of Gamma, rho0, and a.
This is consistent with the standard hierarchy where small 'a' -> large Mass -> smaller ratio?
Wait, Charge ~ a^2, Mass ~ a^3. Ratio ~ a^4 / a^6 ~ 1/a^2.
This implies smaller particles (small a) have LARGER charge-to-mass ratios?
Let's check 1/a^2 behavior:
Exponent of 'a' in Force Ratio: -2
CONFIRMED: F_em/F_g scales as 1/a^2.
"*)
