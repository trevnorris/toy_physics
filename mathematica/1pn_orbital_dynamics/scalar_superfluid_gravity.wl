(* ================================================================= *)
(* SCALAR SUPERFLUID GRAVITY: LENSING (CORRECTED)                  *)
(* ================================================================= *)
ClearAll["Global`*"]

(* 1. PARAMETERS *)
(* mu = GM *)
(* c0 = Speed of light/sound in vacuum *)
(* rho0 = Vacuum density *)
(* polyN = Polytropic Index (P ~ rho^polyN) *)

(* 2. THE TENSION FIELD (4D BULK) *)
(* We postulate the mass acts as a drain/sink, creating a Pressure Deficit. *)
(* Scaling is 1/r (Flux Tube geometry). *)
(* Sign is NEGATIVE (Suction). *)
deltaP = - (mu * rho0) / r;

(* 3. EQUATION OF STATE RESPONSE *)
(* P ~ rho^n  =>  dP/drho = cs^2 *)
(* Linearized relationships: *)
(* deltaP = c0^2 * deltaRho *)
deltaRho = deltaP / c0^2;

(* Sound speed change: d(cs)/cs = ((n-1)/2) * (d(rho)/rho) *)
fractionalCsChange = ((polyN - 1) / 2) * (deltaRho / rho0);

(* 4. REFRACTIVE INDEX *)
(* Local sound speed *)
localCs = c0 * (1 + fractionalCsChange);

(* Refractive Index n = c0 / c_local *)
(* Note: If localCs drops (slower), n rises (bending IN) *)
nRef = c0 / localCs;

(* Linearize to 1st order in mu *)
nRefLinear = Normal[Series[nRef, {mu, 0, 1}]];

Print["Refractive Index n(r):"]
Print[nRefLinear]

(* Extract coefficient alpha: n(r) = 1 + alpha * mu/(r*c0^2) *)
(* We divide by mu/(r*c0^2) to isolate alpha *)
alpha = Simplify[Coefficient[nRefLinear, mu / (c0^2 * r)]];

Print["Refraction Coefficient (alpha):"]
Print[alpha]

(* 5. CALCULATE BENDING *)
(* Component A: Inflow Drag (Time Dilation) = 2 GM / b c^2 *)
thetaTime = 2 * mu / (b * c0^2);

(* Component B: Refraction (Space Curvature) = 2 alpha GM / b c^2 *)
thetaSpace = 2 * alpha * mu / (b * c0^2);

thetaTotal = Simplify[thetaTime + thetaSpace];
thetaGR    = 4 * mu / (b * c0^2);

(* 6. SOLVE FOR REQUIRED FLUID *)
Print["\n-------------------------------------------------"]
Print["RESULTS"]
Print["-------------------------------------------------"]
Print["Total Scalar Bending: ", thetaTotal]
Print["Target GR Bending:    ", thetaGR]

solution = Solve[thetaTotal == thetaGR, polyN];

Print["\nREQUIRED EQUATION OF STATE (polyN):"]
Print[solution]

(*"
Output:

Refractive Index n(r):
1 + (mu*(-1 + polyN))/(2*c0^2*r)
Refraction Coefficient (alpha):
(-1 + polyN)/2

-------------------------------------------------
RESULTS
-------------------------------------------------
Total Scalar Bending: (mu*(1 + polyN))/(b*c0^2)
Target GR Bending:    (4*mu)/(b*c0^2)

REQUIRED EQUATION OF STATE (polyN):
{{polyN -> 3}}
"*)
