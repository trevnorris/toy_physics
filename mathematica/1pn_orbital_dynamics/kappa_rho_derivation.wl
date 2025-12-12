(* ================================================================== *)
(* RIGOROUS DERIVATION OF KAPPA_RHO (Density Coupling)                *)
(* ================================================================== *)
(* *)
(* Goal: Derive the density response deltaRho(r) of the superfluid    *)
(* to the static potential Phi(r) = -mu/r and extract kappaRho.       *)
(* *)
(* Physics: Hydrostatic Equilibrium in the weak-field (1PN) limit.    *)
(* Justification: We rigorously show dynamic pressure v^2 << Phi.     *)
(* ================================================================== *)

ClearAll["Global`*"];

(* Global assumptions: Define parameters as real and positive *)
(* Note: We purposely do NOT include 'r' here to keep it generic *)
$Assumptions = {mu > 0, cs > 0, rho0 > 0, Q > 0, 
  Element[{mu, cs, rho0, Q}, Reals]};

Print["============================================================="];
Print["STEP 1: Define Potentials and Equation of State"];
Print["============================================================="];

(* 1. The Gravitational/Effective Potential (Static Limit) *)
Phi[r_] := -mu/r;
Print["Potential Phi(r): ", Phi[r]];

(* 2. Equation of State (Barotropic Fluid) *)
(* Linearized around background P0, rho0 *)
PLinear[drho_] := P0 + cs^2 * drho;
Print["Linearized EOS: P(rho) = P0 + cs^2 * (rho - rho0)"];

(* ================================================================== *)
(* STEP 2: Solve Hydrostatic Equilibrium                              *)
(* ================================================================== *)
Print["\n============================================================="];
Print["STEP 2: Hydrostatic Balance Condition"];
Print["============================================================="];

(* Euler Equation: cs^2 * Grad(rho) = - rho * Grad(Phi) *)
(* Linearized: cs^2 * rho'[r] == - rho0 * Phi'[r] *)
eqn = cs^2 * D[rho[r], r] == -rho0 * D[Phi[r], r];

Print["Linearized Hydrostatic Equation:"];
Print["   ", eqn];

(* Solve for rho(r) with boundary condition rho(Infinity) = rho0 *)
sol = DSolve[{eqn, rho[Infinity] == rho0}, rho[r], r];
rhoSol[r_] := Evaluate[rho[r] /. sol[[1]]];

Print["\nSolution for rho(r):"];
Print["   ", rhoSol[r]];

(* ================================================================== *)
(* STEP 3: Analyze Density Perturbation                               *)
(* ================================================================== *)
Print["\n============================================================="];
Print["STEP 3: Extract Fractional Density Perturbation"];
Print["============================================================="];

(* deltaRho = rho(r) - rho0 *)
deltaRho[r_] := Simplify[rhoSol[r] - rho0];
fractionalDelta[r_] := Simplify[deltaRho[r] / rho0];

Print["Fractional Density Change (deltaRho / rho0):"];
Print["   ", fractionalDelta[r]];

(* Check sign: We momentarily assume r > 0 just for this boolean check *)
isPositive = Simplify[fractionalDelta[r] > 0, Assumptions -> r > 0];
Print["\nSign Check:"];
Print["   Is deltaRho positive? ", isPositive];
Print["   (Density INCREASES in the potential well)"];

(* ================================================================== *)
(* STEP 4: Derive Kappa_Rho                                           *)
(* ================================================================== *)
Print["\n============================================================="];
Print["STEP 4: Calculate Effective Mass and KappaRho"];
Print["============================================================="];

(* Defect Mass Model: mEff = rho(r) * VCav *)
(* Reference mass m corresponds to background density rho0 *)
mEff[r_] := m * (1 + fractionalDelta[r]);

(* Inertia correction sigmaRho(r) *)
sigmaRho[r_] := Simplify[(mEff[r] - m)/m];

Print["Inertia correction sigmaRho(r):"];
Print["   ", sigmaRho[r]];

(* Extract coefficient kappaRho from sigma(r) = kappa * mu / (cs^2 * r) *)
kappaRho = Simplify[sigmaRho[r] / (mu/(cs^2 * r))];

Print["\nDERIVED KAPPA_RHO:"];
Print["   ", kappaRho];

(* ================================================================== *)
(* STEP 5: Rigorous Validity Check (v^2 vs Phi)                       *)
(* ================================================================== *)
Print["\n============================================================="];
Print["STEP 5: Rigorous Validity Check (v^2 vs Phi)"];
Print["============================================================="];

(* 1. Velocity Field for a Conserved Sink of strength Q *)
vSink[r_] := Q / (4 Pi r^2);

(* 2. Terms in Bernoulli Equation *)
(* Dynamic Pressure ~ v^2 *)
TermDynamic[r_] := 1/2 * vSink[r]^2;

(* Potential Term ~ 1/r *)
(* Note: using Abs for magnitude comparison *)
TermPotential[r_] := Abs[Phi[r]];

Print["Defining terms:"];
Print["   Dynamic Term ~ v^2 = ", TermDynamic[r]];
Print["   Potential Term ~ 1/r = ", TermPotential[r]];

(* 3. Compute Ratio and Limit *)
Ratio[r_] := TermDynamic[r] / TermPotential[r];

Print["\nRatio (Dynamic / Potential):"];
Print["   ", Simplify[Ratio[r], Assumptions -> r > 0]];

(* FIX: Calculate limit using a dummy variable 'rDummy' to avoid 'alimv' warning *)
LimitRatio = Limit[Ratio[r] /. r -> rDummy, rDummy -> Infinity];

Print["\nLimit of Ratio as r -> Infinity:"];
Print["   ", LimitRatio];

checkApproximation = (LimitRatio == 0);

Print["\nVERDICT:"];
If[checkApproximation && (kappaRho == 1),
    Print["SUCCESS: Dynamic pressure vanishes relative to potential at large r."],
    Print["FAILURE: Approximation invalid or kappaRho != 1."]
];

Print["\nFinal kappaRho value: ", kappaRho];

(*"
Output:

=============================================================
STEP 1: Define Potentials and Equation of State
=============================================================
Potential Phi(r): -(mu/r)
Linearized EOS: P(rho) = P0 + cs^2 * (rho - rho0)

=============================================================
STEP 2: Hydrostatic Balance Condition
=============================================================
Linearized Hydrostatic Equation:
   cs^2*Derivative[1][rho][r] == -((mu*rho0)/r^2)

Solution for rho(r):
   ((mu + cs^2*r)*rho0)/(cs^2*r)

=============================================================
STEP 3: Extract Fractional Density Perturbation
=============================================================
Fractional Density Change (deltaRho / rho0):
   mu/(cs^2*r)

Sign Check:
   Is deltaRho positive? cs^2*mu > 0
   (Density INCREASES in the potential well)

=============================================================
STEP 4: Calculate Effective Mass and KappaRho
=============================================================
Inertia correction sigmaRho(r):
   mu/(cs^2*r)

DERIVED KAPPA_RHO:
   1

=============================================================
STEP 5: Rigorous Validity Check (v^2 vs Phi)
=============================================================
Defining terms:
   Dynamic Term ~ v^2 = Q^2/(32*Pi^2*r^4)
   Potential Term ~ 1/r = Abs[mu/r]

Ratio (Dynamic / Potential):
   Q^2/(32*Pi^2*r^3*Abs[mu])

Limit of Ratio as r -> Infinity:
   0

VERDICT:
SUCCESS: Dynamic pressure vanishes relative to potential at large r.

Final kappaRho value: 1
"*)
