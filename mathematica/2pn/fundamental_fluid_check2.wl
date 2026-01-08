(* ---------------------------------------------------------------------- *)
(* Script: Fundamental_Fluid_Check_2PN.wl *)
(* Purpose: Verify 1PN parameters (n=5, q=1) and perform the *)
(* CRITICAL 2PN TEST: The Static Non-Linearity Check (Phi^2 term) *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: 2PN Static Potential Check"];
Print["======================================================="];

(* 1. Define Variables *)
(* eps: Dimensionless Potential (Phi/c^2). Small parameter. *)
(* v: Dimensionless Velocity (v/c). Small parameter. *)
(* M0: Rest Mass *)
(* n: Polytropic Index *)
(* q: Mass-Density Scaling Exponent *)

(* -------------------------------------------------------- *)
(* CORRECTED MASS FUNCTION (The 2PN Physics) *)
(* -------------------------------------------------------- *)

(* 1. Exact Density-Potential Relation *)
(* From Bernoulli: h + Phi + v^2/2 = const. *)
(* For static fluid (v=0): h(rho) = h0 - Phi *)
(* For polytrope n: h(rho) ~ rho^(n-1) *)
(* Thus: rho/rho0 = (1 + (n-1)*Phi/c^2)^(1/(n-1)) *)

(* 2. Mass-Density Relation *)
(* M ~ rho^q *)
(* Combining them: M(Phi) = M0 * [ (1 + (n-1)*eps)^(1/(n-1)) ]^q *)

Mass[eps_] := M0 * (1 + (n - 1) * eps)^(q / (n - 1));

(* -------------------------------------------------------- *)

(* B. Variable Sound Speed (Vector Sector): Stiff Polytrope *)
(* c_s^2/c^2 = 1 + (n-1) * Phi/c^2 *)
(* Derived from P ~ rho^n and Bernoulli enthalpy h ~ Phi *)
CsSq[eps_] := 1 + (n - 1) * eps;

(* 3. Construct Hybrid Lagrangian *)
(* L = -M(Phi) * c^2 * Sqrt[1 - v^2 / c_s^2(Phi)] *)
(* We drop c^2 prefactor for expansions, keeping it dimensionless relative to mc^2 *)
Lagrangian = -Mass[eps] * Sqrt[1 - v^2 / CsSq[eps]];

(* 4. Perform Series Expansion (Now to 2nd Order!) *)
(* We assume v ~ sqrt(eps) scaling for ordering, but expanding strictly in eps works for static check *)
Print["\n--- Expanding Lagrangian to Order epsilon^2 ---"];
LExpansion = Series[Lagrangian, {eps, 0, 2}, {v, 0, 4}] // Normal;
LExpansion = Simplify[LExpansion];

(* Display the expansion (formatted) *)
Print["Expanded L (approx):"];
Print[LExpansion];


(* 5. Extract Coefficients *)

(* A. Newtonian Potential (Order eps^1, v^0) *)
(* Target: -M0 (Standard Newtonian Gravity) *)
CoeffNewton = Coefficient[LExpansion, eps, 1];
CoeffNewton = Coefficient[CoeffNewton, v, 0];

(* B. 1PN Interaction (Order eps^1, v^2) *)
(* Target: -3/2 * M0 (Matches GR EIH) *)
(* Note: The expansion might have a 1/2 factor from v^2 definition, checking 'v^2 eps' term *)
CoeffInteraction = Coefficient[LExpansion, eps, 1];
CoeffInteraction = Coefficient[CoeffInteraction, v, 2];

(* C. 2PN Static Non-Linearity (Order eps^2, v^0) *)
(* THIS IS THE NEW TEST *)
(* GR Target calculation: *)
(* g_00 = 1 - 2*phi + 2*beta*phi^2 (beta=1 for GR) *)
(* L = -sqrt(g_00) = -sqrt(1 - 2eps + 2eps^2) *)
(* Expansion: - (1 - eps + 1/2 eps^2) = -1 + eps - 1/2 eps^2 *)
(* Wait, sign conventions: *)
(* Your previous script defined potential term as -M0*q*eps. *)
(* If Potential V = -1/r, then eps < 0. *)
(* Let's stick to the algebraic target from the previous logic: *)
(* If Newton is -M0, then 2PN static should be +M0/2 for Schwarzschild. *)
CoeffPhi2 = Coefficient[LExpansion, eps, 2];
CoeffPhi2 = Coefficient[CoeffPhi2, v, 0];


Print["\n--- Extracted Coefficients ---"];
Print["1. Newtonian Potential Coeff (Target: -M0): ", CoeffNewton];
Print["2. 1PN Mixed Term Coeff (Target: -3/2 M0): ", CoeffInteraction];
Print["3. 2PN Phi^2 Coeff (Target: +M0/2): ", CoeffPhi2];


(* 6. Solve/Verify *)

Print["\n--- Verification for (n=5, q=1) ---"];
(* Plug in the known solution *)
TestNewton = Simplify[CoeffNewton /. {n -> 5, q -> 1}];
Test1PN = Simplify[CoeffInteraction /. {n -> 5, q -> 1}];
Test2PN = Simplify[CoeffPhi2 /. {n -> 5, q -> 1}];

Print["Newtonian Result (n=5,q=1): ", TestNewton];
Print["1PN Result (n=5,q=1):       ", Test1PN];
Print["2PN Result (n=5,q=1):       ", Test2PN];

(* 7. Final Verdict *)
Print["\n--- FINAL VERDICT ---"];

If[TestNewton == -M0,
    Print["[PASS] Newtonian Limit holds."],
    Print["[FAIL] Newtonian Limit broken."]
];

If[Test1PN == -3/2 * M0,
    Print["[PASS] 1PN Vector Matching holds."],
    Print["[FAIL] 1PN Vector Matching broken."]
];

Target2PN = M0 / 2;
Diff2PN = Simplify[Test2PN - Target2PN];

If[Diff2PN == 0,
    Print["[PASS] 2PN NON-LINEARITY CONFIRMED! The n=5 EOS reproduces the Schwarzschild beta parameter."],
    Print["[FAIL] 2PN Mismatch. The fluid does not reproduce Schwarzschild non-linearity. Diff: ", Diff2PN]
];
