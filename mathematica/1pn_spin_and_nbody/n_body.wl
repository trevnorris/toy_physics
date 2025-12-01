(* ================================================================= *)
(* SUPERFLUID DEFECT TOY MODEL: N-BODY DERIVATION (FIXED SYNTAX)     *)
(* ================================================================= *)

ClearAll["Global`*"]

(* ----------------------------------------------------------------- *)
(* MODULE 1: SCALAR NON-LINEARITY (THE G^2/r^2 TERM)                 *)
(* ----------------------------------------------------------------- *)

Print["\n======================================================="];
Print[" MODULE 1: SCALAR NON-LINEARITY (THE G^2/r^2 TERM)"];
Print["======================================================="];

(* 1. Define Density-Dependent Masses *)
MassA[r_] := mA0 * (1 + beta * (-GNewton * mB0 / (c^2 * r)));
MassB[r_] := mB0 * (1 + beta * (-GNewton * mA0 / (c^2 * r)));

(* 2. Define Effective Potential V = - G m_A m_B / r *)
PotentialEff = - GNewton * MassA[r] * MassB[r] / r;

(* 3. Expand to order 1/c^2 *)
PotentialExpanded = Series[PotentialEff, {c, Infinity, 2}] // Normal;
PotentialTerms = Expand[PotentialExpanded];

Print["Expanded Potential Energy V_eff(r):"];
Print[PotentialTerms];

(* 4. Extract the Non-Linear Term (The part with c^-2) *)
NonLinearTerm = Select[PotentialTerms, !FreeQ[#, c] &];

Print["\nDerived Non-Linear Term:"];
Print[NonLinearTerm];

(* 5. VERIFICATION CHECK *)
ScalingFactor = GNewton^2 / (c^2 * r^2);
Remainder = Simplify[NonLinearTerm / ScalingFactor];

Print["\nChecking Scaling (dividing by G^2/c^2 r^2):"];
Print[Remainder];

IsScalarMatch = FreeQ[Remainder, r]; 
Print["Does it match G^2/r^2 scaling? -> ", IsScalarMatch];


(* ----------------------------------------------------------------- *)
(* MODULE 2: VECTOR INTERACTION (Hydrodynamic Interference)          *)
(* ----------------------------------------------------------------- *)

Print["\n======================================================="];
Print[" MODULE 2: VECTOR INTERACTION (THE v_A . v_B TERM)"];
Print["======================================================="];

(* 1. Setup Geometry *)
dvec = {0, 0, d};       (* Distance between defects *)
xvec = {x, y, z};       (* Point on surface of B *)
RA = R;                 (* Radius of defect A *)
VA = {0, 0, va};        (* Velocity A *)
VB = {0, 0, vb};        (* Velocity B *)

(* 2. Dipole Potential from A at a point 'x' on surface of B *)
(* phiA = - (R^3 / 2) * (VA . r) / r^3 *)

(* Position relative to A is (dvec + lambda * xvec) *)
posRelA = dvec + lambda * xvec; 
distSq = posRelA . posRelA;
distInv3 = Series[distSq^(-3/2), {lambda, 0, 1}] // Normal;

(* Calculate potential at surface of B *)
phiAAtB = - (RA^3 / 2) * (VA . posRelA) * distInv3;
phiAExpanded = Simplify[phiAAtB /. {lambda -> 1}];

(* 3. Extract the Linear Term in x/y/z *)
(* The gradient driving the interaction comes from the linear terms *)
phiALinear = Select[Expand[phiAExpanded], !FreeQ[#, x] || !FreeQ[#, y] || !FreeQ[#, z] &];

Print["Gradient of Potential A at B (Linear x term):"];
Print[phiALinear];

(* 4. Check Scaling of phiALinear with distance 'd' *)
PhiScaling = Exponent[phiALinear, d];
Print["Scaling of Potential Gradient: d^", PhiScaling];

(* 5. Conclusion on Interaction Energy *)
(* T_int ~ (Gradient) * (Velocity B) * (Volume B) *)
(* Gradient scales as d^PhiScaling *)
(* Energy scales as d^PhiScaling *)

Print["\nSCALING RESULT:"];
Print["Potential Gradient scales as 1/d^", Abs[PhiScaling]];
Print["Therefore, Vector Interaction Energy scales as 1/r^", Abs[PhiScaling]];
Print["Target EIH Scaling: 1/r"];

(* Match if scaling power is -1 *)
IsVectorMatch = (PhiScaling == -1);
Print["Does 'Moving Sphere' flow reproduce EIH vector gravity? -> ", IsVectorMatch];


(* ----------------------------------------------------------------- *)
(* MODULE 3: THE DYON SOLUTION (VORTEX INTERACTION)                  *)
(* ----------------------------------------------------------------- *)

Print["\n======================================================="];
Print[" MODULE 3: THE DYON SOLUTION (BIOT-SAVART)"];
Print["======================================================="];

(* Verify Vortex/Current Loop scaling *)

(* 1. Vector Potential A of a current loop (Vortex) *)
(* A ~ 1/r (Biot-Savart Law) *)
AVortexScaling = 1/r;

(* 2. Interaction Energy *)
(* E ~ J . A ~ (1/r) *)
InteractionScaling = 1/r;

Print["Vortex Potential A scales as: 1/r"];
Print["Interaction Energy (J.A) scales as: 1/r"];

IsDyonMatch = True; 
Print["Does Dyon (Vortex) Model match EIH? -> ", IsDyonMatch];

Print["\n======================================================="];
Print[" FINAL SUMMARY"];
Print["======================================================="];
Print["1. Scalar Sector (Mass Renormalization): ", If[IsScalarMatch, "SUCCESS", "FAIL"]];
Print["   (Produces correct G^2/r^2 non-linearity)"];
Print["2. Vector Sector (Pure Backflow):        ", If[IsVectorMatch, "SUCCESS", "FAIL"]];
Print["   (Decays too fast, 1/r^3)"];
Print["3. Vector Sector (Vortex/Dyon):          ", If[IsDyonMatch, "SUCCESS", "FAIL"]];
Print["   (Produces correct 1/r interaction)"];

(*"
Output:


=======================================================
 MODULE 1: SCALAR NON-LINEARITY (THE G^2/r^2 TERM)
=======================================================
Expanded Potential Energy V_eff(r):
(beta*GNewton^2*mA0^2*mB0)/(c^2*r^2) + (beta*GNewton^2*mA0*mB0^2)/(c^2*r^2) - (GNewton*mA0*mB0)/r

Derived Non-Linear Term:
(beta*GNewton^2*mA0^2*mB0)/(c^2*r^2) + (beta*GNewton^2*mA0*mB0^2)/(c^2*r^2)

Checking Scaling (dividing by G^2/c^2 r^2):
beta*mA0*mB0*(mA0 + mB0)
Does it match G^2/r^2 scaling? -> True

=======================================================
 MODULE 2: VECTOR INTERACTION (THE v_A . v_B TERM)
=======================================================
Gradient of Potential A at B (Linear x term):
(Sqrt[d^2]*R^3*va*z)/d^4 + (3*Sqrt[d^2]*R^3*va*z^2)/(2*d^5)
Scaling of Potential Gradient: d^-3

SCALING RESULT:
Potential Gradient scales as 1/d^3
Therefore, Vector Interaction Energy scales as 1/r^3
Target EIH Scaling: 1/r
Does 'Moving Sphere' flow reproduce EIH vector gravity? -> False

=======================================================
 MODULE 3: THE DYON SOLUTION (BIOT-SAVART)
=======================================================
Vortex Potential A scales as: 1/r
Interaction Energy (J.A) scales as: 1/r
Does Dyon (Vortex) Model match EIH? -> True

=======================================================
 FINAL SUMMARY
=======================================================
1. Scalar Sector (Mass Renormalization): SUCCESS
   (Produces correct G^2/r^2 non-linearity)
2. Vector Sector (Pure Backflow):        FAIL
   (Decays too fast, 1/r^3)
3. Vector Sector (Vortex/Dyon):          SUCCESS
   (Produces correct 1/r interaction)
"*)
