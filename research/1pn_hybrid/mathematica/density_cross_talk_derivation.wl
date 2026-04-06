(* ---------------------------------------------------------------------- *)
(* Script: density_cross_talk_derivation.wl                                *)
(* Purpose: Verify the scalar density-starvation contribution to the       *)
(* mixed Phi v^2/c^2 coefficient and quantify the gap to the EIH target.  *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Density Cross-Talk Verification"];
Print["======================================================="];

MassEff[mSelf_, mOther_, r_, k_] := mSelf * (1 - k * (G * mOther/(r * c^2)));

Print[""];
Print["1. Variable-mass ansatz:"];
Print["   M_eff(r) = m * (1 - k * G M_other/(r c^2))"];

LagrangianToy = -MassEff[m1, m2, r, k] * c^2 * Sqrt[1 - v1^2/c^2];
LExpandedPoly = Normal[Series[LagrangianToy, {c, Infinity, 2}]];
LInteractionOnly = Simplify[LExpandedPoly + m1 * c^2];

Print[""];
Print["2. Expanded Lagrangian through 1/c^2:"];
Print[LInteractionOnly];

TermsOrder0 = Coefficient[LExpandedPoly, c, 0];
TermsOrderMinus2 = Coefficient[LExpandedPoly, c, -2];
TargetCoeff = 3/2;
ActualInteractionTerm = Coefficient[TermsOrderMinus2, v1, 2];
ToyCoeff = Simplify[ActualInteractionTerm/(G * m1 * m2/r)];
Gap = Simplify[TargetCoeff - ToyCoeff];
GapAtK1 = Simplify[Gap /. k -> 1];

Print[""];
Print["3. Newtonian term:"];
Print[TermsOrder0];
Print["1PN interaction term:"];
Print[TermsOrderMinus2];

Print[""];
Print["4. Mixed Phi v^2/c^2 coefficient:"];
Print["   Target coefficient (GR): ", TargetCoeff];
Print["   Scalar-sector coefficient: ", ToyCoeff];
Print["   Remaining gap: ", Gap];
Print["   For k = 1, the scalar sector leaves a gap of ", GapAtK1];

Print[""];
Print["5. Conclusion:"];
Print["   Density starvation alone does not saturate the EIH coefficient."];
Print["   The missing positive contribution must come from the vector/flow sector."];

(*"
Output:

=======================================================
   Paper VI: Density Cross-Talk Verification
=======================================================

1. Variable-mass ansatz:
   M_eff(r) = m * (1 - k * G M_other/(r c^2))

2. Expanded Lagrangian through 1/c^2:
(m1*(-4*G*k*m2*v1^2 + r*v1^4 + 4*c^2*(2*G*k*m2 + r*v1^2)))/(8*c^2*r)

3. Newtonian term:
m1*((G*k*m2)/r + v1^2/2)
1PN interaction term:
(m1*v1^2*(-4*G*k*m2 + r*v1^2))/(8*r)

4. Mixed Phi v^2/c^2 coefficient:
   Target coefficient (GR): 3/2
   Scalar-sector coefficient: -1/2*k
   Remaining gap: (3 + k)/2
   For k = 1, the scalar sector leaves a gap of 2

5. Conclusion:
   Density starvation alone does not saturate the EIH coefficient.
   The missing positive contribution must come from the vector/flow sector.
"*)
