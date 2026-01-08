(* ---------------------------------------------------------------------- *)
(* Script: Throat_Breathing_2PN_Check.wl *)
(* Purpose: Calculate if the "Breathing Mode" (Volume Change) of a Flux Throat *)
(* naturally corrects the 2PN Mass Coefficient to match GR. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Throat Breathing & 2PN Correction"];
Print["======================================================="];

(* 1. Setup Thermodynamics (Fixed from previous check) *)
n = 5;          (* Stiff Superfluid *)
q = 1;          (* Linear Mass-Density Scaling *)

(* 2. Define The Equilibrium Physics *)
(* The radius 'a' is determined by Force Balance: F_open = F_close *)

(* Opening Force: Dynamic/Static Pressure on the 3D->4D entrance *)
(* P ~ rho^n. Force ~ P * Area_Effective. *)
(* From black_hole_scaling.wl: F_open ~ rho^n * a^3 *)
FOpen[rho_, a_] := rho^n * a^3;

(* Closing Force: Surface Tension / Geometry *)
(* We model this as F_close ~ a^lambda *)
(* lambda = 1: Line tension (Tear rim) *)
(* lambda = 2: Surface tension (Bubble) *)
FClose[a_, lambda_] := gamma * a^lambda;

(* 3. Solve for Radius 'a' as function of Density 'rho' *)
(* rho^n * a^3 = gamma * a^lambda  ==>  a^(3-lambda) ~ 1/rho^n *)
(* a(rho) ~ rho^( -n / (3 - lambda) ) *)

(* Define the Breathing Exponent 'k_breath' such that a ~ rho^k_breath *)
kBreath[lambda_] := -n / (3 - lambda);

Print["\n--- 1. Throat Breathing Physics ---"];
Print["Defined Radius Scaling: a ~ rho^k"];
Print["Breathing Exponent k_breath (lambda=1, Line Tension): ", kBreath[1]];
Print["Breathing Exponent k_breath (lambda=2, Surface Tension): ", kBreath[2]];


(* 4. Perform 2PN Expansion with Breathing *)

(* A. Density vs Potential (Bernoulli Exact) *)
(* rho(eps) = rho0 * (1 + (n-1)eps)^(1/(n-1)) *)
Rho[eps_] := (1 + (n - 1) * eps)^(1 / (n - 1));

(* B. Radius vs Potential (The Breathing) *)
(* a(eps) = a0 * (Rho[eps])^k_breath *)
Radius[eps_, lambda_] := (Rho[eps])^kBreath[lambda];

(* C. Effective Mass (Density * Volume) *)
(* M(eps) = M0 * (Rho/rho0)^1 * (a/a0)^3 *)
(* Note: q=1 is implicit in the density factor, volume adds the rest *)
MassEff[eps_, lambda_] := M0 * Rho[eps] * (Radius[eps, lambda])^3;

(* 5. Expand and Extract Coefficients *)

CheckModel[lambda_, Label_] := Module[{LExp, CoeffPhi2, Target, Gap},
    Print["\n--- Testing Model: ", Label, " (lambda=", lambda, ") ---"];
    
    (* Expand Effective Mass to 2nd order in eps *)
    (* We only care about the Mass term expansion because L ~ -M(Phi) *)
    (* The kinetic term corrections are already fixed by n=5 *)
    
    MassExpansion = Series[MassEff[eps, lambda], {eps, 0, 2}] // Normal;
    MassExpansion = Simplify[MassExpansion];
    
    Print["Mass Expansion M(Phi):"];
    Print[MassExpansion];
    
    (* Extract 2PN Coefficient (Coeff of eps^2) *)
    (* Lagrangian L = -M(Phi), so we check the coefficient in M *)
    (* Wait, L ~ +M0/2 * eps^2 in GR? *)
    (* Previous script found Fluid-Only L_2PN = +3/2 M0 * eps^2 *)
    (* We want L_2PN = +1/2 M0 * eps^2 *)
    (* Since L = -M, we want M_2PN = -1/2 M0 *)
    
    CoeffPhi2 = Coefficient[MassExpansion, eps, 2];
    
    Print["2PN Coeff (Computed): ", CoeffPhi2];
    Print["2PN Target (GR):      ", -M0/2]; 
    
    Gap = Simplify[CoeffPhi2 - (-M0/2)];
    
    If[Gap == 0,
        Print["[SUCCESS] The Throat Breathing PERFECTLY cancels the excess mass!"],
        Print["[FAIL] Mismatch. Gap: ", Gap]
    ];
];

(* 6. Run Tests *)

(* Case A: Fixed Radius (The Particle) *)
(* Effectively lambda -> -infinity or just set radius constant *)
Print["\n--- Reference: Rigid Particle (No Breathing) ---"];
MassRigid = Series[M0 * Rho[eps], {eps, 0, 2}] // Normal;
CoeffRigid = Coefficient[MassRigid, eps, 2];
Print["Rigid 2PN Coeff: ", CoeffRigid, " (Matches previous result 3/2 * -1 = -1.5?)"];
(* Note: Previous script output +3/2 for L. L = -M. So M should be -3/2. *)
(* Let's check: Rho ~ 1 + eps + (2-n)/2 eps^2 ? *)
(* For n=5: (1+4eps)^(1/4) ~ 1 + eps - 3/2 eps^2. Correct. *)


(* Case B: Flux Throat with Line Tension (lambda=1) *)
(* Matches the logic in black_hole_scaling.wl *)
CheckModel[1, "Flux Throat (Line Tension)"];

(* Case C: Flux Throat with Area Tension (lambda=2) *)
CheckModel[2, "Flux Throat (Area Tension)"];

(* Case D: Is there a Magic Lambda? *)
(* Solve for lambda that gives Success *)
Print["\n--- Deriving Required Geometry ---"];
(* We need CoeffPhi2 == -1/2 *)
(* M ~ rho * rho^(3k) = rho^(1+3k) *)
(* Expansion of rho^(1+3k): (1+3k) * (-3/2) ? No, careful. *)
(* Let's let Mathematica solve it. *)

GeneralPower = 1 + 3 * k;
(* M ~ (1 + (n-1)eps)^(GeneralPower / (n-1)) *)
(* Expansion 2nd order coeff: *)
(* 1/2 * (GeneralPower/(n-1)) * (GeneralPower/(n-1) - 1) * (n-1)^2 *)
(* = 1/2 * GeneralPower * (GeneralPower - (n-1)) *)
(* Target: -1/2 *)

EqMagic = (1/2 * GeneralPower * (GeneralPower - (n - 1)) == -1/2);
SolMagic = Solve[EqMagic, k];

Print["Required Breathing Exponent k: ", k /. SolMagic];
Print["Corresponding Geometry lambda (from k = -n/(3-lambda)): "];
RequiredLambda = Solve[kBreath[lambda] == (k /. SolMagic)[[1]], lambda];
Print[RequiredLambda];
