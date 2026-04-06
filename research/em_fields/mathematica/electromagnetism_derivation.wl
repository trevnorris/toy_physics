(* ================================================================= *)
(* SUPERFLUID VACUUM: EM & GRAVITY UNIFICATION DERIVATION            *)
(* ================================================================= *)
(* Purpose:                                                          *)
(* 1. Collect the main symbolic checks used in the paper.            *)
(* 2. Keep the omnibus script aligned with the split derivations.    *)
(* ================================================================= *)

ClearAll["Global`*"];

$Assumptions = {
  a > 0, L > 0, c > 0, rho0 > 0, Gamma > 0,
  kappaM > 0, kappaQ > 0, betaVac \[Element] Reals
};

x01 = BesselJZero[0, 1];

(* ================================================================= *)
(* PART 1: DEFECT STABILITY & GEOMETRY                               *)
(* ================================================================= *)
Print["================================================="];
Print["PART 1: DEFECT STABILITY & GEOMETRY"];
Print["================================================="];

Hdimless = L*x01^2 + a^2*Pi^2/L + betaVac*a^2*L;
eqRadial = D[Hdimless, a] == 0;
eqAxial = D[Hdimless, L] == 0;

betaVacFromRadial = Simplify[betaVac /. Solve[eqRadial, betaVac][[1]]];
chiEq = Simplify[(eqAxial /. betaVac -> betaVacFromRadial) /. L -> chi*a];
aspectRatio = Simplify[chi /. Solve[chiEq && chi > 0, chi][[1]]];
numericalRatio = N[aspectRatio];

Print["Derived Stability Condition:"];
Print["Geometric Aspect Ratio (L/a) = ", aspectRatio];
Print["Numerical Value: ", numericalRatio];
Print["betaVac at the minimum = ", Simplify[betaVacFromRadial /. L -> aspectRatio*a]];

(* ================================================================= *)
(* PART 2: HIERARCHY PROBLEM RESOLUTION                              *)
(* ================================================================= *)
Print[""];
Print["================================================="];
Print["PART 2: HIERARCHY PROBLEM RESOLUTION"];
Print["================================================="];

massG = kappaM * rho0 * Pi * a^2 * L /. L -> aspectRatio*a;
chargeQ = kappaQ * rho0 * Pi * a^2 * Gamma;
ratioFixedAspect = Simplify[chargeQ^2 / massG^2];
exponentA = Exponent[ratioFixedAspect, a];

Print["Force Ratio (q^2 / m_G^2) at fixed aspect ratio:"];
Print[ratioFixedAspect];
Print["Exponent of a in q^2 / m_G^2: ", exponentA];

If[exponentA == -2,
   Print["SUCCESS: Ratio scales as 1/a^2."],
   Print["FAILURE: Ratio does not scale as 1/a^2."]
];

Print["Interpretation: As a -> 0, the electromagnetic-to-gravitational hierarchy grows."];

(* ================================================================= *)
(* PART 3: HYDRODYNAMIC MAXWELL DICTIONARY                           *)
(* ================================================================= *)
Print[""];
Print["================================================="];
Print["PART 3: HYDRODYNAMIC MAXWELL DICTIONARY"];
Print["================================================="];

vField = {vx[x, y, z], vy[x, y, z], vz[x, y, z]};
vecB = {
  D[vField[[3]], y] - D[vField[[2]], z],
  D[vField[[1]], z] - D[vField[[3]], x],
  D[vField[[2]], x] - D[vField[[1]], y]
};
divB = Simplify[D[vecB[[1]], x] + D[vecB[[2]], y] + D[vecB[[3]], z]];

Print["Divergence of B (Hydrodynamic Identity): ", divB];
If[divB == 0, Print["CONFIRMED: No magnetic monopoles (div B = 0)."]];

Print[""];
Print["--- Testing Coulomb Potential Origin ---"];

rExpr = Sqrt[x^2 + y^2 + z^2];
potentialFlow = (Q/rExpr) * Cos[omega*t];
vFlow = Grad[potentialFlow, {x, y, z}];
vMagSquared = Simplify[vFlow.vFlow];
vMagSquaredR = Simplify[vMagSquared /. {x^2 + y^2 + z^2 -> r^2}];
termUnsteady = Simplify[(-D[potentialFlow, t]) /. {rExpr -> r}];
termKinetic = Simplify[-(1/2) * vMagSquaredR];

Print["Unsteady Term (Acoustic Pressure) Scaling:"];
Print[Normal[Series[termUnsteady, {r, Infinity, 2}]]];
Print["Kinetic Term (Dynamic Pressure) Scaling:"];
Print[Normal[Series[termKinetic, {r, Infinity, 5}]]];
Print["CONCLUSION:"];
Print["The kinetic term falls as 1/r^4 and is negligible in the far field."];
Print["The unsteady term falls as 1/r and provides the Coulomb potential."];

(* ================================================================= *)
(* PART 4: FORCE LAW MATCHING                                        *)
(* ================================================================= *)
Print[""];
Print["================================================="];
Print["PART 4: FORCE LAW MATCHING"];
Print["================================================="];

uVec = {ux, uy, uz};
vInf = {v0, 0, 0};
kVec = {0, 0, 1};

fMagnus = rho0 * Gamma * Cross[kVec, (uVec - vInf)];
fMagnusU = rho0 * Gamma * Cross[kVec, uVec];
qDef = kappaQ * rho0 * Pi * a^2 * Gamma;
fLorentzMag = (qDef/L) * Cross[uVec, B0*kVec];
solutionB0 = Solve[Simplify[fMagnusU] == Simplify[fLorentzMag], B0];

Print["Magnus force per unit length: ", fMagnus];
Print["Magnetic Lorentz term per unit length: ", Simplify[fLorentzMag]];
Print["B0 required for matching: ", solutionB0];
Print["Result: the u-dependent Magnus term matches the magnetic Lorentz term in the paper's straight-vortex geometry."];

(* ================================================================= *)
(* FINAL OUTPUT SUMMARY                                              *)
(* ================================================================= *)
Print[""];
Print["--- SUMMARY OF DERIVATION ---"];
Print["1. Stability: L/a = ", numericalRatio, " from the cavity enthalpy."];
Print["2. Hierarchy: q^2/m_G^2 ~ 1/a^2 at fixed aspect ratio."];
Print["3. Coulomb sector: the unsteady breathing mode dominates as 1/r."];
Print["4. Magnetic sector: the Magnus and Lorentz terms match after fixing B0."];

(*
Output:
=================================================
PART 1: DEFECT STABILITY & GEOMETRY
=================================================
Derived Stability Condition:
Geometric Aspect Ratio (L/a) = (Sqrt[2]*Pi)/BesselJZero[0, 1]
Numerical Value: 1.8474865771201279
betaVac at the minimum = -1/2*BesselJZero[0, 1]^2/a^2

=================================================
PART 2: HIERARCHY PROBLEM RESOLUTION
=================================================
Force Ratio (q^2 / m_G^2) at fixed aspect ratio:
(Gamma^2*kappaQ^2*BesselJZero[0, 1]^2)/(2*a^2*kappaM^2*Pi^2)
Exponent of a in q^2 / m_G^2: -2
SUCCESS: Ratio scales as 1/a^2.
Interpretation: As a -> 0, the electromagnetic-to-gravitational hierarchy grows.

=================================================
PART 3: HYDRODYNAMIC MAXWELL DICTIONARY
=================================================
Divergence of B (Hydrodynamic Identity): 0
CONFIRMED: No magnetic monopoles (div B = 0).

--- Testing Coulomb Potential Origin ---
Unsteady Term (Acoustic Pressure) Scaling:
(omega*Q*Sin[omega*t])/Sqrt[x^2 + y^2 + z^2]
Kinetic Term (Dynamic Pressure) Scaling:
-1/2*(Q^2*Cos[omega*t]^2)/r^4
CONCLUSION:
The kinetic term falls as 1/r^4 and is negligible in the far field.
The unsteady term falls as 1/r and provides the Coulomb potential.

=================================================
PART 4: FORCE LAW MATCHING
=================================================
Magnus force per unit length: {-(Gamma*rho0*uy), Gamma*rho0*(ux - v0), 0}
Magnetic Lorentz term per unit length: {(a^2*B0*Gamma*kappaQ*Pi*rho0*uy)/L, -((a^2*B0*Gamma*kappaQ*Pi*rho0*ux)/L), 0}
B0 required for matching: {{B0 -> -(L/(a^2*kappaQ*Pi))}}
Result: the u-dependent Magnus term matches the magnetic Lorentz term in the paper's straight-vortex geometry.

--- SUMMARY OF DERIVATION ---
1. Stability: L/a = 1.8474865771201279 from the cavity enthalpy.
2. Hierarchy: q^2/m_G^2 ~ 1/a^2 at fixed aspect ratio.
3. Coulomb sector: the unsteady breathing mode dominates as 1/r.
4. Magnetic sector: the Magnus and Lorentz terms match after fixing B0.
*)
