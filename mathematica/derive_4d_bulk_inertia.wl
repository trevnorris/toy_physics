(* 4D Bulk Inertia Derivation - Corrected *)
(* Goal: Calculate the Added Mass coefficient for a 3-sphere throat moving in 4D Bulk *)

ClearAll["Global`*"]

(* --------------------------------------------------------- *)
(* 1. Geometric Setup & Definitions *)
(* --------------------------------------------------------- *)
Print["--- 1. Geometric Setup ---"]

(* Coordinates: x, y, z are brane dimensions. w is the bulk dimension. *)
(* The throat is a sphere S^3 of radius 'a' centered at origin. *)
(* Motion is along the x-axis with velocity V. *)
(* r4 is the radial coordinate in 4D *)

r4 = Sqrt[x^2 + y^2 + z^2 + w^2];

(* The Potential Ansatz for 4D Dipole *)
(* phi = - (A * V * x) / r4^4 *)
(* We use 'Astrength' instead of A to avoid conflicts *)

phi[x_, y_, z_, w_] := - (Astrength * V * x) / r4^4

(* Verify it satisfies Laplace Equation in 4D *)
laplacian = D[phi[x,y,z,w], {x,2}] + D[phi[x,y,z,w], {y,2}] + 
            D[phi[x,y,z,w], {z,2}] + D[phi[x,y,z,w], {w,2}];
Print["Laplacian Check (should be 0): ", Simplify[laplacian]];


(* --------------------------------------------------------- *)
(* 2. Boundary Conditions to Fix 'Astrength' *)
(* --------------------------------------------------------- *)
Print["--- 2. Solving for Dipole Strength Astrength ---"]

(* Boundary Condition: At r=a, radial velocity of fluid = radial velocity of surface *)
(* v_rad = d(phi)/dr *)
(* For phi ~ x/r^4 = cos(theta)/r^3, d/dr is -3/r^4 *)

(* Using the Ansatz phi = - A * V * x / r^4 *)
(* d/dr (1/r^4) = -4/r^5 *)
(* Wait, d/dr (x/r^4) = x * (-4/r^5) ?? No. *)
(* x = r * cos(theta). *)
(* phi = - A * V * cos(theta) / r^3 *)
(* dphi/dr = - A * V * cos(theta) * (-3/r^4) = 3 * A * V * cos(theta) / r^4 *)

(* Surface velocity v_surf = V * cos(theta) *)
(* Match at r=a: *)
(* 3 * A * V * cos / a^4 = V * cos *)
(* 3 * A / a^4 = 1  =>  A = a^4 / 3 *)

Asol = a^4 / 3;
Print["Dipole Strength Asol: ", Asol];

(* Update Potential *)
phiFinal = phi[x,y,z,w] /. Astrength -> Asol;


(* --------------------------------------------------------- *)
(* 3. Calculate Kinetic Energy in Bulk *)
(* --------------------------------------------------------- *)
Print["--- 3. Integrating Bulk Kinetic Energy ---"]

(* Energy Density T_dens = 1/2 * rho * (grad phi)^2 *)
(* Green's Identity: T = - 1/2 * rho * Integral_Surface (phi * dphi/dr) dA *)

(* On surface r=a: *)
(* phiSurf = - (Asol * V * x) / a^4 *)
(* dphiDrSurf = V * (x/a)   <-- From BC v_rad = v_surf *)

(* Integrand = phi * dphi/dr *)
(* integrand = - (a^4/3 * V * x / a^4) * (V * x / a) *)
(* integrand = - (1/3) * V^2 * (x^2 / a) *)

(* Integral over Surface of S^3 *)
(* Integral(x^2) = (1/4) * a^2 * SurfaceArea(S^3) by symmetry *)
(* Surface Area of unit 3-sphere S^3 is 2*pi^2 * a^3 *)

AreaS3 = 2 * Pi^2 * a^3;
IntX2 = (1/4) * a^2 * AreaS3;

(* Surface Integral Value *)
SurfaceInt = - (1/3) * V^2 * (1/a) * IntX2;

(* Total Kinetic Energy T *)
(* T = - 1/2 * rho * SurfaceInt *)
Tbulk = Simplify[- 1/2 * rho * SurfaceInt];

Print["Total Bulk Kinetic Energy Tbulk: ", Tbulk];


(* --------------------------------------------------------- *)
(* 4. Extract Kappa *)
(* --------------------------------------------------------- *)
Print["--- 4. Extracting Kappa ---"]

(* T = 1/2 * m_eff * V^2 *)
mEff = 2 * Coefficient[Tbulk, V^2];

Print["Effective Bulk Mass mEff: ", mEff];

(* Reference Cavitation Mass (3D Sphere Volume) *)
(* mCav = rho * (4/3) * pi * a^3 *)

mCav = rho * (4/3) * Pi * a^3;
Print["Reference Cavitation Mass (3D): ", mCav];

(* Resulting Kappa *)
kappaPV = Simplify[mEff / mCav];

Print["===================================="];
Print["RESULTING KAPPA_PV: ", kappaPV];
Print["===================================="];

(* Check relative to 4D Displaced Mass *)
mDisp4D = rho * (1/2) * Pi^2 * a^4;
kappa4DRef = Simplify[mEff / mDisp4D];
Print["(Check: Kappa relative to 4D displaced mass): ", kappa4DRef];
