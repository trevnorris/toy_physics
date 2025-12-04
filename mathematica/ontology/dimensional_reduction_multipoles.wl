(* ---------------------------------------------------------------------- *)
(* 3. FINAL ROBUST SCRIPT: Explicit Step-by-Step Integration *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

(* Define the Parameters explicitly as assumptions for simplifications *)
assumePositive = {a > 0, L > 0, r > 0, rho0 > 0};

(* 1. Define Legendre Polynomial P2 *)
P2 = (3 * Cos[theta]^2 - 1) / 2;

(* ---------------------------------------------------------------------- *)
(* 2. Calculate Effective 3D Density (rho3) *)
(* ---------------------------------------------------------------------- *)
Print["--- Step 1: Calculating rho3 ---"];

(* 4D density *)
rho4 = rho0 * Exp[-r^2/a^2] * Exp[-w^2/L^2] * (1 + eps * P2);

(* Integrate over w (-inf to inf) *)
(* We use Assuming[] to pass assumptions cleanly *)
rho3 = Assuming[assumePositive, 
    Integrate[rho4, {w, -Infinity, Infinity}]
];
rho3 = Simplify[rho3];

Print["rho3(r, theta) = "];
Print[rho3 // TraditionalForm];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 3. Total Mass M *)
(* ---------------------------------------------------------------------- *)
Print["--- Step 2: Calculating Mass M ---"];

(* STRATEGY: Integrate Angular parts first, then Radial part separately. *)
(* Integrand = rho3 * r^2 * Sin[theta] *)

(* A. Isolate the angular part of the density *)
(* rho3 looks like: Factor(r) * (Factor(theta)) *)
(* Let's just integrate the whole expression over angles *)

integrandMassAng = rho3 * Sin[theta];

(* Integrate over phi (0 to 2pi) -> just multiplies by 2pi since no phi dependence *)
massPhi = 2 * Pi * integrandMassAng; 

(* Integrate over theta (0 to pi) *)
(* We use 'TrigExpand' to help Mathematica see it's simple polynomials of Sin/Cos *)
massTheta = Integrate[TrigExpand[massPhi], {theta, 0, Pi}];

(* B. Integrate Radial part (0 to inf) *)
(* Don't forget the r^2 from the volume element! *)
(* massTheta currently contains the r dependence from rho3 *)
integrandMassRad = massTheta * r^2; 

M = Assuming[assumePositive, 
    Integrate[integrandMassRad, {r, 0, Infinity}]
];

Print["Total Mass M = "];
Print[M // TraditionalForm];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 4. Quadrupole Moment Q20 *)
(* ---------------------------------------------------------------------- *)
Print["--- Step 3: Calculating Quadrupole Q20 ---"];

(* Integrand = rho3 * r^2 * P2 * (r^2 Sin[theta]) *)
(* = rho3 * r^4 * P2 * Sin[theta]       *)

(* A. Angular Integral *)
integrandQuadAng = rho3 * P2 * Sin[theta];
quadPhi = 2 * Pi * integrandQuadAng; (* Phi integral again just 2pi *)

(* Theta integral *)
quadTheta = Integrate[TrigExpand[quadPhi], {theta, 0, Pi}];

(* B. Radial Integral *)
(* Multiply by r^4 (r^2 from definition, r^2 from volume element) *)
integrandQuadRad = quadTheta * r^2; (* Note: r^2 from volume. The other r^2 is inside quadTheta? *)
(* Wait, let's be careful. Q20 definition is Integral[ rho * r^2 * P2 * dV ] *)
(* dV = r^2 dr sin(theta) dtheta dphi *)
(* So total powers of r is r (from rho) * r^2 (from Q def) * r^2 (from dV) = r^? *)
(* Actually, let's just use the result of Angular integral *)
(* quadTheta has the r-dependence from rho3. We multiply by r^4 (r^2 from Q, r^2 from dV) *)
(* Let's re-assemble carefully: *)
(* Q20 = Int[ (rho3) * (r^2 P2) * (r^2 Sin Theta) ] *)
(* We integrated (rho3 * P2 * Sin Theta) over angles. *)
(* This leaves us with the radial part of rho3. *)
(* We must multiply by r^4 before integrating dr. *)

Q20 = Assuming[assumePositive,
    Integrate[quadTheta * r^2 * r^2, {r, 0, Infinity}]
];

Q20 = Simplify[Q20];

Print["Quadrupole Q20 = "];
Print[Q20 // TraditionalForm];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 5. Ratio *)
(* ---------------------------------------------------------------------- *)
ratio = Simplify[Q20 / M];

Print["Ratio Q20 / M = "];
Print[ratio // TraditionalForm];

(*"
Output:

--- Step 1: Calculating rho3 ---
rho3(r, theta) =
TraditionalForm[(L*Sqrt[Pi]*rho0*(4 + eps + 3*eps*Cos[2*theta]))/(4*E^(r^2/a^2))]

--- Step 2: Calculating Mass M ---
Total Mass M =
TraditionalForm[a^3*L*Pi^2*rho0]

--- Step 3: Calculating Quadrupole Q20 ---
Quadrupole Q20 =
TraditionalForm[(3*a^5*eps*L*Pi^2*rho0)/10]

Ratio Q20 / M =
TraditionalForm[(3*a^2*eps)/10]
"*)
