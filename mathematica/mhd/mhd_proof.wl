(* ================================================================= *)
(* SUPERFLUID TOY MODEL -> MHD PROOF SCRIPT (CORRECTED)              *)
(* ================================================================= *)

ClearAll["Global`*"];

(* 1. SETUP *)
coords = {x, y, z};
v = {vx[x, y, z, t], vy[x, y, z, t], vz[x, y, z, t]};
hField = h[x, y, z, t];

(* Operators *)
grad[f_] := Grad[f, coords];
curl[F_] := Curl[F, coords];
laplacian[f_] := Laplacian[f, coords];
vectorLaplacian[F_] := laplacian /@ F; 

(* Advection explicitly *)
advectionTerm = Table[v . grad[v[[i]]], {i, 1, 3}];

(* 2. DICTIONARY *)
\[Omega] = curl[v];
A = \[Lambda] * v;
\[Phi] = \[Lambda] * (hField + (v . v)/2);

B = curl[A]; 
Ee = -grad[\[Phi]] - D[A, t];

(* 3. EULER SUBSTITUTION (THE FIX) *)
(* ----------------------------------------------------------------- *)
eulerRHS = -grad[hField] + \[Nu] * vectorLaplacian[v];

(* We create a list of rules for each component individually *)
dtRules = Thread[ D[v, t] -> (eulerRHS - advectionTerm) ];

(* 4. VERIFY IDEAL OHM'S LAW *)
(* ----------------------------------------------------------------- *)
Print["--- VERIFYING IDEAL OHM'S LAW ---"];

(* E + v x B *)
ohmCheck = Ee + Cross[v, B];

(* Apply the threaded rules. Viscosity \[Nu] -> 0 for ideal check *)
ohmCheckIdeal = ohmCheck /. dtRules /. \[Nu] -> 0;

ohmResult = FullSimplify[ohmCheckIdeal];

If[ohmResult === {0, 0, 0},
  Print["SUCCESS: E + v x B = 0 identically."],
  Print["FAILURE: ", ohmResult]
];

(* 5. VERIFY INDUCTION EQUATION (CORRECTED STRATEGY) *)
(* ----------------------------------------------------------------- *)
Print["\n--- VERIFYING INDUCTION EQUATION ---"];

(* STRATEGY FIX: *)
(* Instead of D[B,t] which traps the time derivative inside a spatial derivative, *)
(* we use the commutation of operators: dt(Curl v) = Curl(dt v). *)
(* We take the Curl of the Euler RHS directly. *)

(* 1. Calculate the time evolution of v from Euler *)
dtVExplicit = eulerRHS - advectionTerm;

(* 2. Calculate the time evolution of B by taking the Curl of dtV *)
(* dtB = dt (curl A) = curl (dt A) = lambda * curl (dt v) *)
dtBHydro = \[Lambda] * curl[dtVExplicit];

(* 3. Calculate the target MHD terms (same as before) *)
mhdTransport = curl[Cross[v, B]];
resistiveTerm = \[Nu] * vectorLaplacian[B];

(* 4. Compare *)
inductionCheck = dtBHydro - (mhdTransport + resistiveTerm);

inductionResult = FullSimplify[inductionCheck];

If[inductionResult === {0, 0, 0},
  Print["SUCCESS: Induction Equation holds! Viscosity maps to Resistivity."],
  Print["FAILURE: ", inductionResult]
];

(*"
Output:

--- VERIFYING IDEAL OHM'S LAW ---
SUCCESS: E + v x B = 0 identically.

--- VERIFYING INDUCTION EQUATION ---
SUCCESS: Induction Equation holds! Viscosity maps to Resistivity.
"*)
