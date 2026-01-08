(* ---------------------------------------------------------------------- *)
(* Script: Dark_Star_2PN_Check.wl *)
(* Purpose: Test if a "Dark Star" (Radiation-Supported Throat) naturally *)
(* reproduces the 2PN nonlinearity of General Relativity. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: The Dark Star Hypothesis"];
Print["   (Radiation-Supported Throat vs 2PN Gravity)"];
Print["======================================================="];

(* 1. Setup Thermodynamics *)
n = 5;          (* Stiff Superfluid *)
q = 1;          (* Linear Mass coupling M ~ rho * Vol *)

(* 2. Define The Physics Models *)

(* A. Fluid Density (Bernoulli Exact) *)
(* rho(eps) = rho0 * (1 + (n-1)eps)^(1/(n-1)) *)
Rho[eps_] := (1 + (n - 1) * eps)^(1 / (n - 1));

(* B. Opening Force (Trapped Radiation) *)
(* Energy E is constant (trapped). Volume V ~ a^3. *)
(* Radiation Pressure P_rad = 1/3 * E / V ~ a^-3 *)
(* Force_open ~ P_rad * Area ~ a^-3 * a^2 ~ a^-1 *)
FOpen[a_] := A_rad / a;

(* C. Closing Force Models *)
(* Model 1: Vacuum Pressure (The fluid crushes the bubble) *)
(* F ~ P_vac * Area ~ rho^n * a^2 *)
FClosePressure[rho_, a_] := rho^n * a^2;

(* Model 2: Surface Tension (The rim pulls shut) *)
(* F ~ gamma * Length ~ rho^k * a^1 *)
(* We assume tension scales with density as rho^k_tension *)
FCloseTension[rho_, a_, kTension_] := rho^kTension * a;


(* 3. Solve for Radius Scaling a(rho) *)

(* Case 1: Radiation vs Vacuum Pressure *)
(* A/a = rho^n * a^2  =>  a^3 ~ 1/rho^n  =>  a ~ rho^(-n/3) *)
kRadPressure = -n / 3;

(* Case 2: Radiation vs Surface Tension *)
(* A/a = rho^k * a    =>  a^2 ~ 1/rho^k  =>  a ~ rho^(-k/2) *)
(* We will leave k_tension as a parameter *)

(* 4. Define Mass Calculation Helper *)
(* M(eps) = M0 * (Rho/rho0) * (a/a0)^3 *)
Calculate2PN[k_, ModelName_] := Module[{LExp, MassExpansion, CoeffPhi2, Gap},
    Print["\n--- Testing: ", ModelName, " ---"];
    Print["Radius Scaling: a ~ rho^(", k, ")"];
    
    (* Radius(eps) *)
    RadiusFunc = (Rho[eps])^k;
    
    (* Mass(eps) *)
    MassFunc = M0 * Rho[eps] * RadiusFunc^3;
    
    (* Expansion *)
    MassExpansion = Series[MassFunc, {eps, 0, 2}] // Normal // Simplify;
    Print["Mass Expansion: ", MassExpansion];
    
    CoeffPhi2 = Coefficient[MassExpansion, eps, 2];
    Print["2PN Coeff (Computed): ", CoeffPhi2];
    Print["2PN Target (GR):      ", -M0/2];
    
    Gap = Simplify[CoeffPhi2 - (-M0/2)];
    If[Gap == 0, 
        Print["[SUCCESS] MATCHES GR EXACTLY!"],
        Print["[FAIL] Mismatch. Gap: ", Gap]
    ];
    Return[Gap];
];


(* 5. Run The Tests *)

(* Reference: Rigid Particle *)
Calculate2PN[0, "Rigid Particle (Fixed Radius)"];

(* Test 1: Radiation Pressure vs Vacuum Pressure *)
Calculate2PN[kRadPressure, "Dark Star (Vacuum Crush)"];

(* Test 2: Find the Required Scaling *)
Print["\n--- Solving for Required Geometry ---"];
(* We want CoeffPhi2 == -1/2 M0 *)
(* Mass ~ rho * rho^(3k) = rho^(1+3k) *)
(* Expansion of rho^p is 1 + p*eps + 1/2*p*(p-(n-1))*eps^2 ... wait. *)
(* Let's use the code to solve it exactly. *)

(* Define p = 1+3k *)
(* Coeff of eps^2 for (1+(n-1)eps)^(p/(n-1)) is: *)
(* 1/2 * (p/(n-1)) * (p/(n-1) - 1) * (n-1)^2  =  1/2 * p * (p - (n-1)) *)

(* We want 1/2 * p * (p - (n-1)) == -1/2 *)
(* p^2 - (n-1)p + 1 = 0 *)
(* For n=5: p^2 - 4p + 1 = 0 *)

SolP = Solve[p^2 - 4*p + 1 == 0, p];
Print["Required Power p (where M ~ rho^p): ", p /. SolP];

(* Convert p back to k (radius scaling) *)
(* p = 1 + 3k  =>  k = (p-1)/3 *)
RequiredK = (p - 1)/3 /. SolP;
Print["Required Radius Scaling k (a ~ rho^k): ", RequiredK];

Print["\n--- Interpretation ---"];
Print["k_required_1 approx: ", N[RequiredK[[1]]]];
Print["k_required_2 approx: ", N[RequiredK[[2]]]];

Print["\nChecking Model 1 (Vacuum Crush): k = ", kRadPressure];
Print["Checking Model 2 (Surface Tension): Can we match k?"];
(* If a ~ rho^(-k_tension/2) and we need k ~ 0.26 or 0.73... *)
(* -k_tension/2 = 0.26 => k_tension = -0.52 *)
(* This implies Surface Tension DECREASES as density increases? Inverse tension? *)
