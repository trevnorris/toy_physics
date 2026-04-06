(* ---------------------------------------------------------------------- *)
(* Script: Fundamental_Fluid_Check.wl *)
(* Purpose: Derive 2PN coefficients from the Fundamental Equation of State *)
(* P ~ rho^n without assuming the metric form. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Zero-Assumption Fluid Derivation"];
Print["======================================================="];

(* 1. Define Fundamentals *)
(* Variables: n (Polytropic Index), q (Mass Scaling Power) *)
(* phi: Gravitational Potential (Phi/c^2), small parameter *)
(* v: Velocity (v/c), small parameter *)

(* Equation of State: P = K * rho^n *)
(* Speed of Sound c_s^2 = dP/drho ~ rho^(n-1) *)
(* We linearize around background density rho0 *)
(* rho = rho0 * (1 + delta) *)

(* Relation between Density and Potential *)
(* Enthalpy h = Integral(dP/rho) ~ Phi *)
(* h ~ rho^(n-1) *)
(* Therefore: rho^(n-1) ~ Phi -> rho ~ Phi^(1/(n-1)) ?? *)
(* Careful: Linearized perturbation. *)
(* dh = c_s^2 * drho/rho *)
(* h = c_s0^2 * delta / (n-1)? No. *)
(* h = n/(n-1) * P/rho = c_s^2 / (n-1). *)
(* Perturbation: dh = d(c_s^2)/(n-1) = d(rho^(n-1))/(n-1) ~ rho^(n-2) drho *)
(* Let's stick to: dP = -rho dPhi (Hydrostatic/Bernoulli balance) *)
(* c_s^2 drho = -rho dPhi *)
(* drho/rho = - dPhi / c_s^2 *)
(* Integrating: log(rho) = -Phi/c_s^2 (for isothermal) *)
(* For polytrope: density perturbation delta = - Phi / c_s0^2 *)
(* Wait, potential Phi is NEGATIVE (-GM/r). *)
(* So delta = - (-|Phi|) = +|Phi|. Density INCREASES in a potential well? *)
(* NO. The defect is a SINK. The flow generates the potential. *)
(* Paper I/II: Phi_Newton ~ - Enthalpy. *)
(* High density = High Enthalpy. *)
(* Low density (near sink) = Low Enthalpy = Low Potential. *)
(* So Density drops near the sink. *)
(* Relation: delta_rho / rho0 = Phi / c_s0^2 *)
(* (Phi is negative, so delta_rho is negative. Correct.) *)

(* 2. Define Perturbed Quantities *)
(* Let small parameter eps = Phi_Newton / c^2 *)
(* We assume c_s0 = c (Speed of light/sound background) *)

(* Density Perturbation: rho(r) = rho0 * (1 + eps) *)
(* Note: eps is negative (-GM/rc^2) *)

(* Speed of Sound Perturbation: c_s(r) *)
(* c_s^2 ~ rho^(n-1) *)
(* c_s(r) = c * (1 + eps)^((n-1)/2) *)

(* Mass Perturbation: M(r) *)
(* Hypothesis: M ~ rho^q *)
(* M(r) = M0 * (1 + eps)^q *)

(* 3. Construct the Lagrangian *)
(* L = - M(r) * c^2 * Sqrt[ 1 - v^2 / c_s(r)^2 ] *)
(* Note: The 'c' in the Gamma factor is the LOCAL speed of sound c_s(r). *)
(* The 'c^2' prefactor is the rest energy scale. *)

c = 1; (* Normalize background speed to 1 for algebra *)
cs[eps_, n_] := c * (1 + eps)^((n-1)/2);
Mass[eps_, q_] := M0 * (1 + eps)^q;

Lagrangian = - Mass[eps, q] * c^2 * Sqrt[ 1 - v^2 / cs[eps, n]^2 ];

(* 4. Expand to 2PN Order *)
(* We need terms of order: eps (Newtonian), v^2 (Kinetic), eps*v^2 (Interaction) *)
SeriesExp = Series[Lagrangian, {eps, 0, 1}, {v, 0, 2}];
PolyExp = Normal[SeriesExp];

Print["\n1. Expanded Lagrangian (terms of eps and v^2):"];
Print[Simplify[PolyExp]];

(* 5. Extract Coefficients *)
(* Term 1: Kinetic Energy (coeff of v^2) *)
(* Should be 1/2 M0 v^2 *)
KineticTerm = Coefficient[PolyExp, v, 2];
KineticTerm0 = Coefficient[KineticTerm, eps, 0]; 
Print["\n2. Newtonian Kinetic Term (should be M0/2):"];
Print[KineticTerm0];

(* Term 2: Potential Energy (coeff of eps) *)
(* Should be M0 * eps (Since eps = Phi = -GM/r, and V = M Phi) *)
(* Lag = T - V = ... - M Phi. So coeff should be -M0. *)
PotentialTerm = Coefficient[PolyExp, eps, 1];
PotentialTerm0 = Coefficient[PotentialTerm, v, 0];
Print["\n3. Newtonian Potential Term (should be -M0):"];
Print[PotentialTerm0];

(* Term 3: Interaction Energy (coeff of eps * v^2) *)
(* This is the 2PN term G M M / r * v^2 *)
(* GR Requirement: + 3/2 * (G M M / r) * v^2 *)
(* In our notation: + 3/2 * (-eps) * M0 * v^2 ? No. *)
(* Interaction in L_EIH is + 3/2 G M M / r v^2. *)
(* Our eps is - G M / r. *)
(* So we look for term: - 3/2 * eps * M0 * v^2. *)
(* Wait, check signs carefully. *)
(* L_EIH ~ G M M / r ( 1 + 3/2 v^2/c^2 ) *)
(* L ~ - V ( 1 + ... ) = + G M M / r + 3/2 G M M / r v^2. *)
(* eps = - G M / r. *)
(* So L ~ - M0 * eps + 3/2 * M0 * (-eps) * v^2 ? NO. *)
(* + 3/2 * (G M M / r) * v^2 = + 3/2 * (-eps * M0) * v^2 = -3/2 M0 eps v^2. *)
(* So we need the coefficient of (M0 * eps * v^2) to be -3/2. *)

InteractionTerm = Coefficient[KineticTerm, eps, 1];
Print["\n4. Interaction Term Coeff (target -3/2):"];
Print[Simplify[InteractionTerm / M0]];

(* 6. Solve for n and q *)
(* Condition 1: Potential Term = -M0 (Recovers Newton) *)
Eq1 = PotentialTerm0 == -M0;

(* Condition 2: Interaction Term = -3/2 M0 (Recovers GR) *)
Eq2 = InteractionTerm == -3/2 * M0;

Sol = Solve[{Eq1, Eq2}, {n, q}];

Print["\n5. SOLUTION FOR PARAMETERS:"];
Print[Sol];

(* 7. Check n=5 case *)
Print["\n6. Checking Stiff Superfluid (n=5):"];
If[Length[Sol] > 0,
   qSol = q /. Sol[[1]];
   nSol = n /. Sol[[1]];
   Print["   Required q: ", qSol];
   Print["   Required n: ", nSol];
   
   If[nSol == 5, 
      Print["   RESULT: n=5 MATCHES EXACTLY."],
      Print["   RESULT: n=5 does not match. Requires n=", nSol]
   ];
   ,
   Print["   No solution found."]
];

(*"
Output:

=======================================================
   Paper VI: Zero-Assumption Fluid Derivation
=======================================================

1. Expanded Lagrangian (terms of eps and v^2):
(M0*(-2 + v^2 - eps*(-1 + n)*v^2 + eps*q*(-2 + v^2)))/2

2. Newtonian Kinetic Term (should be M0/2):
M0/2

3. Newtonian Potential Term (should be -M0):
-(M0*q)

4. Interaction Term Coeff (target -3/2):
(1 - n + q)/2

5. SOLUTION FOR PARAMETERS:
{{n -> 5, q -> 1}}

6. Checking Stiff Superfluid (n=5):
   Required q: 1
   Required n: 5
   RESULT: n=5 MATCHES EXACTLY.
"*)
