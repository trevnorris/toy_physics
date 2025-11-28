(* ================================================================= *)
(* RIGOROUS 2-BODY LAGRANGIAN SOLVER (RECURSION FIXED) *)
(* ================================================================= *)
Print["--- EXACT SYMMETRIC LAGRANGIAN DERIVATION ---"];

ClearAll["Global`*"];

(* 1. Construct Lagrangian directly with functions of t *)
(* r and v_rel definitions *)
rVec = {x1[t] - x2[t], y1[t] - y2[t], z1[t] - z2[t]};
vVec = {x1'[t] - x2'[t], y1'[t] - y2'[t], z1'[t] - z2'[t]};
rMag = Sqrt[rVec . rVec];
vRelSq = vVec . vVec;

(* Kinetic Energy *)
Tkin = 1/2 * m1 * (x1'[t]^2 + y1'[t]^2 + z1'[t]^2) +
       1/2 * m2 * (x2'[t]^2 + y2'[t]^2 + z2'[t]^2);

(* Potential Energy *)
BasePot = -(G * m1 * m2 / rMag) - (chi / rMag^2);
Coupling = (1 + alpha * vRelSq / c^2);
Vpot = BasePot * Coupling;

L = Tkin - Vpot;

(* 2. Euler-Lagrange Equations (X-component) *)
(* Eq: D[ D[L, x'[t]], t ] - D[L, x[t]] == 0 *)

(* Body 1 *)
dLdvx1 = D[L, x1'[t]];
dLdx1  = D[L, x1[t]];
dtDLdvx1 = D[dLdvx1, t]; (* This generates x1''[t] terms *)

(* Body 2 *)
dLdvx2 = D[L, x2'[t]];
dLdx2  = D[L, x2[t]];
dtDLdvx2 = D[dLdvx2, t];

(* 3. Substitute Algebraic Symbols for Accelerations *)
(* We replace the derivatives x1''[t] with symbols ax1, ax2 to solve linearly *)
accelRules = {
    x1''[t] -> ax1, y1''[t] -> ay1, z1''[t] -> az1,
    x2''[t] -> ax2, y2''[t] -> ay2, z2''[t] -> az2
};

Eq1 = ((dtDLdvx1 - dLdx1) /. accelRules) == 0;
Eq2 = ((dtDLdvx2 - dLdx2) /. accelRules) == 0;

(* 4. SOLVE THE COUPLED SYSTEM *)
Print["Solving coupled system for ax1 and ax2..."];
(* We solve for ax1 and ax2 simultaneously. *)
(* This breaks the circular dependency because we treat them as linear unknowns. *)
solutions = Solve[{Eq1, Eq2}, {ax1, ax2}];

If[solutions === {},
   Print["Error: No solution found."],

   (* 5. Extract Body 1 Acceleration *)
   rawA1 = ax1 /. solutions[[1]];

   (* 6. Series Expansion *)
   Print["\n--- DERIVED ACCELERATION a1 (Series Expanded) ---"];

   (* Replace t-functions with readable symbols for the output *)
   readableA1 = rawA1 /. {
       x1[t] -> x1, x2[t] -> x2,
       x1'[t] -> vx1, x2'[t] -> vx2,
       Sqrt[(x1[t]-x2[t])^2 + (y1[t]-y2[t])^2 + (z1[t]-z2[t])^2] -> r
   };

   (* Expand to O(1/c^2) *)
   seriesA1 = Series[readableA1, {c, Infinity, 2}];
   simplifiedA1 = Simplify[Normal[seriesA1]];

   Print[simplifiedA1];

   Print["\n--- INTERPRETATION ---"];
   Print["Check the last term. If it contains (vx1 - vx2)*(x1 - x2),"];
   Print["that is the radial drag/anti-damping term."];
];
