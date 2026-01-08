(* ---------------------------------------------------------------------- *)
(* Script: Density_Cross_Talk_Derivation.wl *)
(* Purpose: Derive the effective Lagrangian terms arising from "Density Starvation" *)
(* and compare them to the standard EIH 1PN/2PN requirements. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Density Cross-Talk & 2PN Verification"];
Print["======================================================="];

(* 1. Define Physics Constants & Parameters *)
(* G: Gravitational Constant, c: Speed of Light *)
(* m1, m2: Rest masses of the bodies *)
(* v1: Velocity of body 1 *)
(* r: Separation distance *)
(* k: 'Coupling Strength' of the density effect. 
      Hypothesis: rho_local ~ rho0 * (1 - Phi/c^2) implies k = 1 *)

(* 2. Define The Density Starvation Hypothesis *)
(* The effective mass of Body 1 depends on the gravitational potential of Body 2 *)
(* Potential Phi = -G * m2 / r *)
(* Variable Mass M1(r) = m1 * (1 - k * G * m2 / (r * c^2)) *)
(* Note: We use a minus sign because 'density drops' near the other body. *)

MassEff[mSelf_, mOther_, r_, k_] := mSelf * (1 - k * (G * mOther / (r * c^2)));

Print["\n1. Defined Variable Mass Model:"];
Print["   M_eff(r) = m * (1 - k * G * M_other / (r * c^2))"];

(* 3. Construct the Relativistic Lagrangian *)
(* L = -M_eff(r) * c^2 * Sqrt[1 - v^2/c^2] *)
(* This couples the position-dependent mass to the relativistic kinetic term. *)

LagrangianToy = -MassEff[m1, m2, r, k] * c^2 * Sqrt[1 - v1^2/c^2];

(* 4. Perform PN Expansion (Series expansion in 1/c) *)
(* We expand to order (1/c)^2 relative to the Newtonian term. *)
(* This effectively captures terms up to c^-2 (1PN) and c^-4 (2PN precursor). *)

LExpanded = Series[LagrangianToy, {c, Infinity, 2}];
LExpandedPoly = Normal[LExpanded]; (* Convert series to polynomial *)

Print["\n2. Expanded Lagrangian (Terms up to 1/c^2):"];
(* Filter out the rest mass constant for clarity *)
LInteractionOnly = Simplify[LExpandedPoly + m1*c^2]; 
Print[LInteractionOnly];

(* 5. Extract Specific Interaction Terms *)
(* We are looking for the v^2/r term (The 'Kinetic Interaction') *)
(* Format: Coeff * (G m1 m2 / r) * v1^2 *)

(* Isolate the coefficient of (G m1 m2 v1^2 / r) *)
(* We divide the expression by the target variables to find the dimensionless scalar factor *)
TermOfInterest = G * m1 * m2 * v1^2 / (r * c^2);
ComputedCoeff = Simplify[ Coefficient[LExpandedPoly, TermOfInterest] * TermOfInterest / TermOfInterest ];

(* Note: Coefficient[] might be tricky with division, so we use a robust extraction: *)
(* Differentiate twice with respect to v1, then look at the potential part. *)
ExtractedCoeff = Simplify[ D[D[LExpandedPoly, v1], v1] / 2 / (G * m2 / (r * c^2)) / m1 ] * -1; 
(* The * -1 adjusts for the potential sign conventions if needed. 
   Let's just look at the raw expanded output manually in the print loop below for safety. *)

(* Let's collect terms by power of c *)
TermsOrder0 = Coefficient[LExpandedPoly, c, 0]; (* Newtonian Energy *)
TermsOrderMinus2 = Coefficient[LExpandedPoly, c, -2]; (* 1PN Correction *)

Print["\n3. Coefficient Analysis:"];
Print["   Newtonian Term (Order c^0):"];
Print[TermsOrder0]; 
(* Should be 1/2 m v^2 + G m m / r *)

Print["   1PN Term (Order c^-2):"];
Print[TermsOrderMinus2];

(* 6. Define The EIH Target (Harmonic Coordinates) *)
(* Standard GR 1PN Lagrangian Interaction part: *)
(* L_EIH_Interact = (G m1 m2 / r) * [ 3/2 * (v1^2 + v2^2) / c^2 ... ] *)
(* We strictly look at the v1^2 term coefficient. Target is +3/2. *)

TargetCoeff = 3/2;

(* Extract the actual v1^2 coeff from our TermsOrderMinus2 *)
(* TermsOrderMinus2 contains: 1/8 m1 v1^4 (Relativistic Kinetic Correction) *)
(* + Coeff * G m1 m2 v1^2 / r (The interaction) *)

ActualInteractionTerm = Coefficient[TermsOrderMinus2, v1, 2];
(* This effectively isolates the term linear in G and m2 *)
ToyCoeff = Simplify[ ActualInteractionTerm / (G * m1 * m2 / r) ];

Print["\n4. Results vs EIH Target:"];
Print["   Target Coefficient (GR): ", TargetCoeff];
Print["   Toy Model Coefficient:   ", ToyCoeff];

(* 7. Calculate the Gap *)
Gap = TargetCoeff - ToyCoeff;

Print["\n5. The Gap (Missing Physics):"];
Print["   Delta = Target - Toy = ", Gap];

Print["\n6. Conclusion for Paper VI:"];
If[Gap == 0,
    Print["SUCCESS: Density Starvation explains 100% of the dynamics."],
    Print["PARTIAL: Density Starvation explains part. The Vector Kernel must provide the remaining coefficient: ", Gap]
];

(* Output results for easy reading *)
{TargetCoeff, ToyCoeff, Gap}

(*"
Output:

=======================================================
   Paper VI: Density Cross-Talk & 2PN Verification
=======================================================

1. Defined Variable Mass Model:
   M_eff(r) = m * (1 - k * G * M_other / (r * c^2))

2. Expanded Lagrangian (Terms up to 1/c^2):
(m1*(-4*G*k*m2*v1^2 + r*v1^4 + 4*c^2*(2*G*k*m2 + r*v1^2)))/(8*c^2*r)

3. Coefficient Analysis:
   Newtonian Term (Order c^0):
m1*((G*k*m2)/r + v1^2/2)
   1PN Term (Order c^-2):
(m1*v1^2*(-4*G*k*m2 + r*v1^2))/(8*r)

4. Results vs EIH Target:
   Target Coefficient (GR): 3/2
   Toy Model Coefficient:   -1/2*k

5. The Gap (Missing Physics):
   Delta = Target - Toy = 3/2 + k/2
"*)
