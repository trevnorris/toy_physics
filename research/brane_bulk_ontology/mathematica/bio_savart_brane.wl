(* ================================================================= *)
(* BRANE TRANSVERSE WAKE AND THE BIOT-SAVART STRUCTURE               *)
(* ================================================================= *)
(* Purpose:                                                          *)
(* 1. Verify that the irrotational bulk potential gives Curl[v] = 0. *)
(* 2. Verify that a brane-localized transverse wake A ~ u/r gives    *)
(*    B = Curl[A] ~ (u x r)/r^3.                                     *)
(* 3. Match the result to the paper's vector-Poisson construction.   *)
(* ================================================================= *)

ClearAll["Global`*"]

Print["--- 1. Setup & Definitions ---"];

rVec = {x, y, z};
rMag = Sqrt[x^2 + y^2 + z^2];
uVec = {ux, uy, uz};
kappaA = Symbol["kappaA"];
qEff = Symbol["Q"];

biotSavartTarget =
  Simplify[(kappaA*qEff/(4*Pi)) * Cross[uVec, rVec] / rMag^3];

Print["Target Biot-Savart form B = (kappaA Q / 4 Pi) (u x r)/r^3:"];
Print[biotSavartTarget];


Print["\n--- 2. The Bulk: Irrotational Potential Flow ---"];

phiBulk = (uVec . rVec) / rMag^3;
vBulk = Grad[phiBulk, {x, y, z}];
bBulk = Simplify[Curl[vBulk, {x, y, z}]];

Print["Bulk vorticity Curl[v_bulk]:"];
Print[bBulk];

If[AllTrue[bBulk, # == 0 &],
   Print["RESULT: Curl[v_bulk] = 0 outside the core, as required for the gravitational sector."],
   Print["RESULT: Non-zero bulk vorticity found."]
];


Print["\n--- 3. The Brane: Transverse Wake Potential ---"];

aBrane = (kappaA*qEff/(4*Pi)) * uVec / rMag;
bBrane = Simplify[Curl[aBrane, {x, y, z}]];

Print["Brane vector potential A(r) = (kappaA Q / 4 Pi) u/r:"];
Print[Simplify[aBrane]];
Print["Derived magnetic field B = Curl[A]:"];
Print[bBrane];


Print["\n--- 4. Comparison with the Target Structure ---"];

matchCheck = Simplify[bBrane - biotSavartTarget];

Print["Difference between derived B and target Biot-Savart form:"];
Print[matchCheck];

If[AllTrue[matchCheck, # == 0 &],
   Print["SUCCESS: The brane transverse wake reproduces the Biot-Savart structure exactly."],
   Print["FAILURE: Structure mismatch."]
];


Print["\n--- Conclusion ---"];
Print["The bulk potential flow remains irrotational, so it cannot supply far-zone magnetostatics."];
Print["The brane-localized transverse wake A ~ u/r yields B ~ (u x r)/r^3 with the expected normalization."];
Print["This matches the paper's vector-Poisson interpretation and does not rely on a Stokeslet/viscous-drag picture."];

(* 
Output:
--- 1. Setup & Definitions ---
Target Biot-Savart form B = (kappaA Q / 4 Pi) (u x r)/r^3:
{(kappaA*Q*(-(uz*y) + uy*z))/(4*Pi*(x^2 + y^2 + z^2)^(3/2)), (kappaA*Q*(uz*x - ux*z))/(4*Pi*(x^2 + y^2 + z^2)^(3/2)), (kappaA*Q*(-(uy*x) + ux*y))/(4*Pi*(x^2 + y^2 + z^2)^(3/2))}

--- 2. The Bulk: Irrotational Potential Flow ---
Bulk vorticity Curl[v_bulk]:
{0, 0, 0}
RESULT: Curl[v_bulk] = 0 outside the core, as required for the gravitational sector.

--- 3. The Brane: Transverse Wake Potential ---
Brane vector potential A(r) = (kappaA Q / 4 Pi) u/r:
{(kappaA*Q*ux)/(4*Pi*Sqrt[x^2 + y^2 + z^2]), (kappaA*Q*uy)/(4*Pi*Sqrt[x^2 + y^2 + z^2]), (kappaA*Q*uz)/(4*Pi*Sqrt[x^2 + y^2 + z^2])}
Derived magnetic field B = Curl[A]:
{(kappaA*Q*(-(uz*y) + uy*z))/(4*Pi*(x^2 + y^2 + z^2)^(3/2)), (kappaA*Q*(uz*x - ux*z))/(4*Pi*(x^2 + y^2 + z^2)^(3/2)), (kappaA*Q*(-(uy*x) + ux*y))/(4*Pi*(x^2 + y^2 + z^2)^(3/2))}

--- 4. Comparison with the Target Structure ---
Difference between derived B and target Biot-Savart form:
{0, 0, 0}
SUCCESS: The brane transverse wake reproduces the Biot-Savart structure exactly.

--- Conclusion ---
The bulk potential flow remains irrotational, so it cannot supply far-zone magnetostatics.
The brane-localized transverse wake A ~ u/r yields B ~ (u x r)/r^3 with the expected normalization.
This matches the paper's vector-Poisson interpretation and does not rely on a Stokeslet/viscous-drag picture.
*)
