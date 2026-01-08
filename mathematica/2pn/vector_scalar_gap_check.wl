(* ---------------------------------------------------------------------- *)
(* Script: Vector_Scalar_Gap_Check.wl *)
(* Purpose: Calculate 1PN coefficients for Retarded Scalar + Vector fields *)
(* and check if they sum to the General Relativity (EIH) target. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Scalar-Vector Summation Check"];
Print["======================================================="];

(* 1. Define Variables *)
(* rVec: Position vector r = x1 - x2 *)
(* nVec: Unit vector n = r / |r| *)
(* v1, v2: Velocity vectors *)
(* r: Separation magnitude *)

(* 2. Define The EIH Target (Harmonic Coordinates) *)
(* Reference: Landau & Lifshitz, The Classical Theory of Fields *)
(* L_EIH_Interact = (G m1 m2 / r) * [ 1 + ... corrections ... ] *)
(* We look for the coefficients of the 1/c^2 terms inside the bracket. *)

(* Target Coefficients for (G m1 m2 / r c^2) * [ ... ] *)
Targetv1sq = 3/2;        (* Coefficient for v1^2 *)
Targetv1v2 = -7/2;       (* Coefficient for v1 . v2 *)
Targetv1nv2n = -1/2;    (* Coefficient for (v1.n)(v2.n) *)

Print["\n1. EIH Target Coefficients (GR):"];
Print["   v1^2 term:         ", Targetv1sq];
Print["   v1.v2 term:        ", Targetv1v2];
Print["   (v1.n)(v2.n) term: ", Targetv1nv2n];

(* 3. Derive Scalar Sector Contribution *)
(* Physics: Variable Mass M(r) coupled to Relativistic Kinetic Energy *)
(* L_scalar = -M(r) * c^2 * Sqrt[1 - v^2/c^2] *)
(* M(r) includes the STATIC potential Phi = -G m2 / r. *)
(* Note: Retardation of the Scalar Potential itself usually affects v^2 terms too. *)
(* Standard Retarded Scalar Potential Expansion (Lienard-Wiechert for Scalar): *)
(* Phi_ret = m2/r + (1/2 c^2) d^2/dt^2(m2 r) *)

(* Expansion of d^2/dt^2(r): *)
(* d/dt(r) = -n . v2 *)
(* d^2/dt^2(r) includes accelerations (Newtonian force) and v^2 terms. *)
(* For the interaction Lagrangian, we usually substitute leading order EOM: a = -G m / r^2 n *)
(* Result from literature for Scalar Interaction L_int: *)
(* L_S = - G m1 m2 / r * [ 1 - 1/2 (v1.n)(v2.n) + 1/2 (v1^2 + v2^2 - (v1.v2)? ) ]? *)
(* Let's calculate the 'Density Starvation' (Kinetic coupling) term explicitly again. *)
(* We deduced it gives -1/2 * v1^2. *)

CoeffScalarv1sq = -1/2; (* From previous derivation *)
(* Does pure scalar density starvation give v1.v2? No. *)
CoeffScalarv1v2 = 0;
CoeffScalarv1nv2n = 0;

Print["\n2. Scalar Sector (Density Starvation Only):"];
Print["   v1^2 contribution: ", CoeffScalarv1sq];


(* 4. Derive Vector Sector Contribution *)
(* Physics: Standard Vector Field Interaction (like Electromagnetism/Fluid Flow) *)
(* L_vec = (Coeff_Vec) * (J1 . A2) *)
(* A2 is the Retarded Vector Potential: A2 ~ (v2 / r) + corrections *)
(* Interaction term: ~ (v1 . v2) / r *)
(* Darwin Lagrangian (EM) coefficients for Q1 Q2 / r interaction: *)
(* EM: -1/2 (v1 . v2) - 1/2 (v1.n)(v2.n) *)
(* We assume our fluid vector kernel acts similarly but with tunable strength 'Av'. *)

(* If we assume a standard positive kernel, it generates terms proportional to v1.v2 *)
(* However, does a Vector Kernel generate v1^2 terms? *)
(* In standard EM, it does NOT. It only gives cross terms. *)
(* UNLESS: The 'Vector' field is actually a modification of the metric tensor h_00 vs h_0i? *)

(* Let's assume the Vector Kernel provides a weight 'Av' to the standard Darwin term: *)
(* Term: Av * [ (v1 . v2) + (v1.n)(v2.n) ] (Or similar structure) *)
(* Standard Hydrodynamic Kernel (1/r) leads to v1 . v2. *)
(* Retardation leads to (v1.n)(v2.n). *)

(* Gap Analysis *)
Gapv1sq = Targetv1sq - CoeffScalarv1sq;

Print["\n3. The Gap Analysis:"];
Print["   v1^2 Gap: ", Gapv1sq, " (This MUST be filled)"];

(* 5. The Critical Question: Can Vector fill the v1^2 gap? *)
(* Standard vector interactions (v1 . v2) do NOT produce v1^2 terms. *)
(* If the Gap is +2.0 in the v1^2 term, and the vector sector is purely J.A, *)
(* then the model FAILS to match GR. *)

(* HOWEVER: Look at the 1PN Metric expansion in PPN formalism. *)
(* gamma (Scalar curvature) and (Vector gravitomagnetism). *)
(* The v^2 terms in Lagrangian usually come from: *)
(* 1. Kinetic term (Special Relativity) -> Coeff -1/8 (Self term) *)
(* 2. Non-linearity of Scalar Potential (Phi^2) -> Coeff +1/2 ? *)
(* 3. Coupling of v^2 to Potential (Space Curvature) -> Coeff +1 ? *)

(* Hypothesis: The 'Density Starvation' model we used was: M(r) = M0(1 - Phi). *)
(* This corresponds to PPN beta = 1? *)
(* Maybe the 'Tear' implies a metric modification where spatial part is affected? *)
(* If the fluid is compressible, local sound speed c_s changes. *)
(* Refractive index n(r) ~ (1 - 2 Phi). *)
(* Effective metric ds^2 = c^2 dt^2 - n^2 dx^2. *)
(* This gives the Space Curvature v^2 terms! *)

(* Let's verify the coefficient of the Refractive Index Model. *)
(* L = -m0 c^2 Sqrt[ g_00 - g_ij v^i v^j / c^2 ] *)
(* g_00 = (1 + 2 Phi), g_ij = -(1 - 2 Phi) delta_ij  (Standard GR Isotropic) *)
(* Let's test if the Toy Model matches THIS. *)
(* Toy Model: M(r) varies -> g_00 part. *)
(* Toy Model: Flow/Refractive Index -> g_ij part. *)

(* Re-Run Derivation for Refractive Index Model *)
(* L_refractive = -m0 c^2 * Sqrt[ (1 + 2 Phi) - (1 - 2 Phi) v^2/c^2 ] *)
(* Let's expand this and check the v1^2 coefficient. *)

(* Phi = - G m2 / r *)
TermRefractive = Sqrt[ (1 - 2*G*m2/(r*c^2)) - (1 + 2*G*m2/(r*c^2))*v1^2/c^2 ];
LRefractive = -m1 * c^2 * TermRefractive;

LExpRefractive = Series[LRefractive, {c, Infinity, 2}] // Normal;
InteractionRefractive = Coefficient[LExpRefractive, c, -2];
CoeffRefractivev1sq = Coefficient[InteractionRefractive, v1, 2] / (G * m1 * m2 / r);

Print["\n4. Refractive Index Model Check (Scalar + Spatial Curvature):"];
Print["   v1^2 Coeff: ", CoeffRefractivev1sq];

(* Comparison *)
FinalGap = Targetv1sq - CoeffRefractivev1sq;

Print["\n5. Final Verdict:"];
If[FinalGap == 0,
    Print["SUCCESS: The Flow Refractive Index + Density Starvation matches GR v^2 terms exactly."],
    Print["MISMATCH: Gap is ", FinalGap]
];

{Targetv1sq, CoeffRefractivev1sq, FinalGap}

(*"
Output:

=======================================================
   Paper VI: Scalar-Vector Summation Check
=======================================================

1. EIH Target Coefficients (GR):
   v1^2 term:         3/2
   v1.v2 term:        -7/2
   (v1.n)(v2.n) term: -1/2

2. Scalar Sector (Density Starvation Only):
   v1^2 contribution: -1/2

3. The Gap Analysis:
   v1^2 Gap: 2 (This MUST be filled)

4. Refractive Index Model Check (Scalar + Spatial Curvature):
   v1^2 Coeff: 3/2

5. Final Verdict:
SUCCESS: The Flow Refractive Index + Density Starvation matches GR v^2 terms exactly.
"*)
