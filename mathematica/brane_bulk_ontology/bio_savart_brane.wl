(* ================================================================= *)
(* DERIVATION OF BIOT-SAVART LAW FROM BRANE VORTICITY                *)
(* ================================================================= *)
(* Purpose:                                                          *)
(* 1. Demonstrate that Bulk Potential Flow yields B = 0.             *)
(* 2. Demonstrate that a "Brane Wake" (Stokeslet) yields B ~ u x r.  *)
(* 3. Confirm this matches the Biot-Savart Law for a point charge.   *)
(* ================================================================= *)

ClearAll["Global`*"]

(* ----------------------------------------------------------------- *)
(* 1. SETUP & PARAMETERS                                             *)
(* ----------------------------------------------------------------- *)
Print["--- 1. Setup & Definitions ---"];

(* Cartesian Coordinates *)
rVec = {x, y, z};
rMag = Sqrt[x^2 + y^2 + z^2];

(* Defect Velocity (Current Source) *)
(* Assume defect is at origin, moving with velocity uVec *)
uVec = {ux, uy, uz};

(* Standard Biot-Savart Law for a Point Charge q moving with v *)
(* B_BS ~ (v x r) / r^3 *)
biotSavartStandard = Cross[uVec, rVec] / rMag^3;

Print["Target Form (Standard Biot-Savart ~ u x r / r^3):"];
Print[biotSavartStandard];


(* ----------------------------------------------------------------- *)
(* 2. THE BULK: IRROTATIONAL POTENTIAL FLOW                          *)
(* ----------------------------------------------------------------- *)
Print["\n--- 2. The Bulk: Irrotational Potential Flow ---"];

(* A moving body in an inviscid fluid creates a Dipole Potential *)
(* Phi ~ u . r / r^3 *)
phiBulk = (uVec . rVec) / rMag^3;

(* Velocity v = Grad(Phi) *)
vBulk = Grad[phiBulk, {x, y, z}];

(* Magnetic Field B = Curl(v) *)
bBulk = Curl[vBulk, {x, y, z}];

Print["Bulk Vorticity (B-field):"];
Print[Simplify[bBulk]];

If[AllTrue[Simplify[bBulk], # == 0 &],
   Print["RESULT: B = 0 in the Bulk (As expected for potential flow)."],
   Print["RESULT: Non-zero B in Bulk."]
];


(* ----------------------------------------------------------------- *)
(* 3. THE BRANE: INDUCED VORTICITY (STOKESLET)                       *)
(* ----------------------------------------------------------------- *)
Print["\n--- 3. The Brane: Induced Vorticity (Stokeslet Model) ---"];
Print["Hypothesis: The Brane has an effective viscosity/structure."];
Print["A moving defect creates a drag wake scaling like 1/r (Stokeslet)."];

(* The Stokeslet Velocity Field (Flow due to a point force in viscous fluid) *)
(* v_brane ~ (u/r) + ((u.r)r / r^3) *)
(* We calculate the Curl of the dominant 1/r term (Green's function) *)
(* Effective Vector Potential A ~ v_brane *)

(* Simplified Ansatz: A ~ u / r  (The Liénard-Wiechert potential) *)
(* Note: This matches the 'Stokeslet' leading order term *)
vBraneSimple = uVec / rMag;

(* Calculate Derived B-field *)
bBrane = Curl[vBraneSimple, {x, y, z}];

Print["Calculated Brane Vorticity (B-field):"];
Print[Simplify[bBrane]];

(* ----------------------------------------------------------------- *)
(* 4. COMPARISON & VERIFICATION                                      *)
(* ----------------------------------------------------------------- *)
Print["\n--- 4. Comparison with Biot-Savart ---"];

(* Check if bBrane is proportional to u x r / r^3 *)
(* We compute the ratio or cross product check *)

(* Standard Identity: Curl(u/r) = Grad(1/r) x u = -(r/r^3) x u = (u x r)/r^3 *)
(* Let's see if Mathematica agrees *)

matchCheck = Simplify[bBrane - (-Cross[rVec, uVec] / rMag^3)];
(* Note: Cross[u, r] = -Cross[r, u]. The signs should align. *)

Print["Difference between Derived B and Target Biot-Savart (scaling):"];
Print[matchCheck];

If[AllTrue[Simplify[matchCheck], # == 0 &],
   Print["SUCCESS: Brane Vorticity matches Biot-Savart structure exactly!"],
   Print["FAILURE: Structure mismatch."]
];

(* ----------------------------------------------------------------- *)
(* 5. PHYSICAL INTERPRETATION                                        *)
(* ----------------------------------------------------------------- *)
Print["\n--- Conclusion ---"];
Print["In the Bulk, v ~ Grad(1/r^2), so Curl(v) = 0."];
Print["On the Brane, drag induces v ~ 1/r, so Curl(v) ~ 1/r^2 (Biot-Savart)."];
Print["The 'Magnetic Field' is the vorticity of the Brane Wake."];
