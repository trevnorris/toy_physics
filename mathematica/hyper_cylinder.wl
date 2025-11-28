(* 4D Hyper-Cylinder Added Mass Derivation - Corrected *)
(* Goal: Calculate Added Mass for a Sphere extruded into a Cylinder in 4D *)

ClearAll["Global`*"]

Print["================================================="];
Print["   4D HYPER-CYLINDER (TUBE) INERTIA TEST       "];
Print["================================================="];

(* --------------------------------------------------------- *)
(* 1. Geometric & Potential Setup *)
(* --------------------------------------------------------- *)
Print["\n--- 1. Setup ---"]

(* Coordinates: *)
(* 3D radial coordinate: r3 = Sqrt[x^2 + y^2 + z^2] *)
(* w is the axial coordinate of the cylinder (-L/2 to L/2) *)

r3 = Sqrt[x^2 + y^2 + z^2];

(* Potential Ansatz: *)
(* 3D Dipole Potential: phi = - A * x / r3^3 *)
(* Using 'Astrength' to avoid underscore or single-letter issues *)

phi[x_, y_, z_] := - (Astrength * V * x) / r3^3;

(* Verify Laplacian in 4D *)
(* d^2/dw^2 is 0 since phi is independent of w *)
lap3D = D[phi[x,y,z], {x,2}] + D[phi[x,y,z], {y,2}] + D[phi[x,y,z], {z,2}];
Print["Laplacian Check (should be 0): ", Simplify[lap3D]];

(* --------------------------------------------------------- *)
(* 2. Fix Dipole Strength *)
(* --------------------------------------------------------- *)
Print["\n--- 2. Boundary Condition ---"]

(* Radial velocity in 3D slice: v_r = dphi/dr3 *)
(* Surface normal velocity: v_surf = V * cos(theta) *)
(* Asol = a^3 / 2 derived from standard sphere BC *)

Asol = a^3 / 2;
Print["Dipole Strength Asol (Standard 3D Sphere): ", Asol];

phiFinal = phi[x,y,z] /. Astrength -> Asol;

(* --------------------------------------------------------- *)
(* 3. Integrate Kinetic Energy *)
(* --------------------------------------------------------- *)
Print["\n--- 3. Integrating Bulk Kinetic Energy ---"]

(* Energy Density Tdens = 1/2 * rho * (grad phi)^2 *)
gradPhi = {D[phiFinal, x], D[phiFinal, y], D[phiFinal, z]};
vSquared = Simplify[gradPhi . gradPhi];

(* Integration Strategy: *)
(* Tbulk = L * Tslice *)
(* Tslice = - 1/2 * rho * Integral_Surface (phi * dphi/dr) dA *)

(* On Surface r=a: *)
(* phiSurf = - (a^3/2 * V * x) / a^3 = - V * x / 2 *)
(* dphiDrSurf = V * (x/a) *)

(* Integrand = - ( - V x / 2 ) * ( V x / a ) = V^2 * x^2 / (2a) *)

(* Surface Integral over Sphere S^2: *)
(* Int(x^2) dA = (4/3) * Pi * a^4 *)

IntX2 = (4/3) * Pi * a^4;
SurfaceInt = (V^2 / (2*a)) * IntX2;
Tslice = 1/2 * rho * SurfaceInt;

(* Total Tbulk = L * Tslice *)
Tbulk = Simplify[L * Tslice];

Print["Total Kinetic Energy (Infinite Cylinder Approx): ", Tbulk];

(* --------------------------------------------------------- *)
(* 4. Extract Kappa *)
(* --------------------------------------------------------- *)
Print["\n--- 4. Extracting Kappa ---"]

(* T = 1/2 * mEff * V^2 *)
mEff = 2 * Coefficient[Tbulk, V^2];
Print["Effective Bulk Mass mEff: ", mEff];

(* Displaced Mass of Hyper-Cylinder *)
(* V_disp = Volume(Sphere) * Length = (4/3 pi a^3) * L *)
mDisp = rho * (4/3) * Pi * a^3 * L;

Print["Reference Displaced Mass (Tube): ", mDisp];

kappa = Simplify[mEff / mDisp];

Print["===================================="];
Print["RESULTING KAPPA (Tube): ", kappa];
Print["===================================="];

(* --------------------------------------------------------- *)
(* 5. Interpretation *)
(* --------------------------------------------------------- *)
Print["\n--- Interpretation ---"];
If[kappa === 1/2,
    Print["RESULT: 1/2. The tube behaves like a stack of spheres."],
    Print["RESULT: " <> ToString[kappa]]
];
