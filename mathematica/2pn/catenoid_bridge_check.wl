(* ---------------------------------------------------------------------- *)
(* Script: Catenoid_Bridge_Check.wl *)
(* Purpose: Test if a "Dual Puncture" Catenoid Bridge supported by *)
(* Trapped Radiation Pressure naturally reproduces the GR 2PN coefficient. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: The Catenoid Bridge Check"];
Print["   (Dual Puncture Geometry vs 2PN Gravity)"];
Print["======================================================="];

(* 1. Setup Thermodynamics *)
n = 5;           (* Stiff Superfluid *)
q = 1;           (* Linear Mass coupling *)

(* 2. Define Physics Models *)

(* A. Fluid Density (Bernoulli Exact) *)
(* rho(eps) = rho0 * (1 + (n-1)eps)^(1/(n-1)) *)
Rho[eps_] := (1 + (n - 1) * eps)^(1 / (n - 1));

(* B. Opening Force: Trapped Radiation in a Catenoid *)
(* A Catenoid throat has effective volume V ~ a^3 (for L ~ a) *)
(* Radiation Pressure P_rad ~ E / V ~ a^-3 *)
(* Force_open ~ P_rad * Area_Throat ~ a^-3 * a^2 ~ a^-1 *)
(* This 1/a scaling is characteristic of "pressure in a tube" *)
FOpen[a_] := A_rad / a;

(* C. Closing Force: Bulk Surface Tension *)
(* For a Catenoid (Minimal Surface), the dominant restoring force is Tension *)
(* Force ~ Gamma * Circumference ~ Gamma * a *)
(* How does Gamma scale with density? *)
(* In a superfluid, Surface Tension gamma ~ Energy/Area ~ rho * c_s^2 * xi ? *)
(* Dimensional analysis for n=5: gamma ~ rho^((n+1)/2)? *)
(* Let's leave the scaling exponent k_gamma as a free parameter to find the "Magic" value *)
FClose[rho_, a_, kGamma_] := rho^kGamma * a;


(* 3. Solve for Radius Scaling a(rho) *)

(* Balance: A/a = rho^k * a  =>  a^2 ~ rho^-k  =>  a ~ rho^(-k/2) *)
(* We define the 'Breathing Exponent' k_breath = -kGamma / 2 *)
kBreath[kGamma_] := -kGamma / 2;


(* 4. Define Mass Calculation *)

(* The critical question: How does Mass scale for a Bridge? *)
(* Hypothesis 1: M ~ rho * Volume ~ rho * a^3 *)
(* Hypothesis 2: M ~ Flux ~ rho * a^2 * c_s *)
(* Hypothesis 3: M ~ Trapped Energy (Constant? No, E redshifts) *)

(* We will test Hypothesis 1 (Standard Defect) first, as used in previous scripts *)
(* M(eps) = M0 * (Rho/rho0) * (a/a0)^3 *)

CheckBridge[kGamma_, Label_] := Module[{k, MassExpansion, CoeffPhi2, Gap},
    Print["\n--- Testing Model: ", Label, " ---"];
    
    (* 1. Determine Radius Scaling *)
    k = kBreath[kGamma];
    Print["Surface Tension Scaling: gamma ~ rho^", kGamma];
    Print["Resulting Radius Scaling: a ~ rho^", k];
    
    (* 2. Expand Mass *)
    (* M ~ rho^(1 + 3k) *)
    
    RadiusFunc = (Rho[eps])^k;
    MassFunc = M0 * Rho[eps] * RadiusFunc^3;
    
    MassExpansion = Series[MassFunc, {eps, 0, 2}] // Normal // Simplify;
    Print["Mass Expansion M(Phi):"];
    Print[MassExpansion];
    
    (* 3. Check 2PN Coefficient *)
    (* We want M_2PN = -1/2 M0 (which corresponds to L_2PN = +1/2 M0) *)
    
    CoeffPhi2 = Coefficient[MassExpansion, eps, 2];
    Print["2PN Coeff (Computed): ", CoeffPhi2];
    Print["2PN Target (GR):      ", -M0/2];
    
    Gap = Simplify[CoeffPhi2 - (-M0/2)];
    
    If[Gap == 0, 
        Print["[SUCCESS] MATCHES GR EXACTLY!"],
        Print["[FAIL] Mismatch. Factor: ", Simplify[CoeffPhi2 / (-M0/2)]]
    ];
];


(* 5. Run The Tests *)

(* Case A: Linear Surface Tension (gamma ~ rho) *)
(* This is the simplest assumption: Tension is proportional to density *)
CheckBridge[1, "Catenoid with Linear Tension (gamma ~ rho)"];

(* Case B: Energy Density Surface Tension (gamma ~ rho * c_s^2) *)
(* For n=5, c_s^2 ~ rho^4. So gamma ~ rho^5. *)
CheckBridge[5, "Catenoid with Stiff Tension (gamma ~ P)"];

(* Case C: The 'Factor of 3' Check *)
(* We assume the 'Rigid' result was 1.5 M0. We want 0.5 M0. *)
(* This requires M ~ rho^(1/3). *)
(* M ~ rho^(1+3k). So 1+3k = 1/3 => 3k = -2/3 => k = -2/9. *)
(* If k = -kGamma/2, then -kGamma/2 = -2/9 => kGamma = 4/9. *)
CheckBridge[4/9, "The 'Factor of 3' Hypothesis (k = -2/9)"];


(* 6. Reverse Engineer the Required Physics *)
Print["\n--- Solving for Required Surface Tension Physics ---"];
(* We need CoeffPhi2 == -1/2 M0 *)
(* Mass ~ rho^p where p = 1 + 3*(-kGamma/2) *)
(* Coeff of eps^2 for rho^p is 1/2 * p * (p - (n-1)) *)
(* We want 1/2 * p * (p - 4) = -1/2  (for n=5) *)
(* p^2 - 4p + 1 = 0 *)

SolP = Solve[p^2 - 4*p + 1 == 0, p];
pVals = p /. SolP;
Print["Required Mass Scaling Power p (M ~ rho^p): ", pVals];
Print["(Note: p=0.268 is close to 1/3, but not exactly 1/3)"];

(* Convert p to kGamma *)
(* p = 1 - 3/2 kGamma  =>  kGamma = 2/3 (1 - p) *)
ReqGamma = 2/3 * (1 - pVals);

Print["Required Surface Tension Scaling gamma ~ rho^k_gamma:"];
Print[ReqGamma];
Print["Numerical values: ", N[ReqGamma]];
