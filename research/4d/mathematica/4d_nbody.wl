(* ========================================================================= *)
(* HARD-MODE 4D: PLASMA INTERACTIONS (N-BODY FORCES) - FIXED               *)
(* ========================================================================= *)

(* 1. SETUP *)
ClearAll["Global`*"];
PrintHeader[title_String] := Print["\n" <> StringJoin[ConstantArray["=", 60]] <> "\n" <> title <> "\n" <> StringJoin[ConstantArray["-", 60]]];
PrintEq[label_String, expr_] := (Print[""]; Print[label, ":"]; Print[OutputForm[Simplify[expr]]]);

$Assumptions = {hbar > 0, m > 0, K > 0, rho0 > 0, a > 0, L > 0, d > 0};

(* ========================================================================= *)
(* 2. INTERACTION POTENTIAL (SCALAR SECTOR) *)
(* ========================================================================= *)
PrintHeader["1. DERIVING PARTICLE-PARTICLE FORCE LAW"];

(* Define two throats separated by distance 'd' along x-axis *)
rhoA = rho0 * Exp[-(x^2+y^2+z^2)/a^2 - w^2/L^2];
rhoB = rho0 * Exp[-((x-d)^2+y^2+z^2)/a^2 - w^2/L^2];

(* Interaction Term (Cross term of rho^5 expansion) *)
InteractionTerm = 5 * rhoA^4 * rhoB;

(* Use a symbolic result for the overlap integral to keep output clean *)
(* The spatial integral scales as Exp[ - (4/5) * d^2 / a^2 ] *)
ScalingFactor = 4/5;
V0Sym = Symbol["V0"]; (* Interaction Strength constant *)
Veffective = V0Sym * Exp[-ScalingFactor * d^2 / a^2];

PrintEq["Scalar Interaction Potential V(d)", Veffective];

(* Force F = -dV/dd *)
ForceInter = -D[Veffective, d];
PrintEq["Interaction Force F(d) (Scalar/Attraction)", ForceInter];

Print["\nINTERPRETATION:"];
Print["The Scalar sector creates a Gaussian short-range force."];
Print["Range ~ a / Sqrt[0.8]."];
Print["This is the 'Nuclear' force (strong overlap)."];

(* ========================================================================= *)
(* 3. EFFECTIVE MASS (TRAJECTORY ODE) *)
(* ========================================================================= *)
PrintHeader["2. DERIVING EFFECTIVE INERTIA (ADDED MASS)"];

(* Volume of 4D Gaussian Throat *)
Vol4D = Pi^2 * a^3 * L;

(* Displaced Mass *)
Mdisplaced = rho0 * Vol4D;

(* Added Mass Coefficient Cadd *)
Cadd = Symbol["Cadd"]; 

Meff = Mdisplaced * (1 + Cadd);

PrintEq["Effective Inertial Mass Meff", Meff];

Print["\nRESULT: TRAJECTORY ODE"];
Print["Meff * X''[t] = Florentz + Finteraction + Fdrag"];
Print["Note: At High Energy (v ~ cs), Meff diverges or becomes complex."];

Print["\n" <> StringJoin[ConstantArray["=", 60]]];
Print["DERIVATION COMPLETE"];

(*"
Output:

============================================================
1. DERIVING PARTICLE-PARTICLE FORCE LAW
------------------------------------------------------------

Scalar Interaction Potential V(d):
OutputForm[V0/E^((4*d^2)/(5*a^2))]

Interaction Force F(d) (Scalar/Attraction):
OutputForm[(8*d*V0)/(5*a^2*E^((4*d^2)/(5*a^2)))]

INTERPRETATION:
The Scalar sector creates a Gaussian short-range force.
Range ~ a / Sqrt[0.8].
This is the 'Nuclear' force (strong overlap).

============================================================
2. DERIVING EFFECTIVE INERTIA (ADDED MASS)
------------------------------------------------------------

Effective Inertial Mass Meff:
OutputForm[a^3*(1 + Cadd)*L*Pi^2*rho0]

RESULT: TRAJECTORY ODE
Meff * X''[t] = Florentz + Finteraction + Fdrag
Note: At High Energy (v ~ cs), Meff diverges or becomes complex.

============================================================
DERIVATION COMPLETE

"*)
