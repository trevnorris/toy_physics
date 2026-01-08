(* ---------------------------------------------------------------------- *)
(* Script: Final_Vortex_Bridge_Check_v2.wl *)
(* Purpose: Test if a Vortex-Supported Bridge with Linear Surface Tension *)
(* solves the "Factor of 3" problem. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: The Vortex Bridge Solution"];
Print["======================================================="];

(* 1. Physics Parameters *)
n = 5;                  (* Stiff Superfluid *)
TargetRatio = 1.8475;   (* The EM Geometric Eigenstate *)

(* 2. Define Exact Geometry (Catenoid) *)
(* V = Pi * a^3 * [ (L/a) + 0.5 * Sinh[2*L/a] ] *)
VolumeExact[a_, L_] := Pi * a^3 * ( (L/a) + 0.5 * Sinh[2*L/a] );

(* 3. Define Thermodynamics *)
Rho[eps_] := (1 + (n - 1) * eps)^(1 / (n - 1));

(* 4. The Solver Module *)
(* We input the Surface Tension Exponent 'kGamma' *)
(* gamma ~ rho^kGamma *)

CheckVortexBridge[kGamma_] := Module[
    {L, a0, TensionVal, FClose, FOpen, k1Val, k2Val, MassSeries, 
     c2, dev, xSeries, BalanceEq, Coeff1, Sol1, Coeff2, Sol2, 
     RadiusExpr, MassExpr},
    
    (* Define Tension Scaling *)
    TensionVal[eps_] := (Rho[eps])^kGamma;
    
    (* Define Geometry *)
    a0 = 1.0; 
    L = TargetRatio * a0;
    
    (* FORCE BALANCE *)
    (* 1. Opening: Vortex Ring (Stiff) -> F ~ 1/a^3 *)
    (* We normalize constants so at eps=0 (a=1, rho=1), F_open = F_close *)
    (* F_close(0) = 1 * 2 * Pi * 1 = 2 Pi *)
    (* So F_open(0) must be 2 Pi. Constant is 2 Pi. *)
    FOpen[a_] := (2 * Pi) / a^3; 
    
    (* 2. Closing: Surface Tension -> F ~ Gamma * Circumference *)
    FClose[a_, eps_] := TensionVal[eps] * 2 * Pi * a;
    
    (* Solve Expansion Coefficients for a(eps) *)
    (* Ansatz: a = a0(1 + k1 eps + k2 eps^2) *)
    xSeries = k1*eps + k2*eps^2;
    
    BalanceEq = Series[
        FOpen[a0*(1+xSeries)] - FClose[a0*(1+xSeries), eps], 
        {eps, 0, 2}
    ];
    
    (* Solve Order 1 *)
    Coeff1 = Coefficient[BalanceEq, eps, 1];
    Sol1 = Solve[Coeff1 == 0, k1];
    k1Val = k1 /. Sol1[[1]];
    
    (* Solve Order 2 *)
    Coeff2 = Coefficient[BalanceEq /. k1 -> k1Val, eps, 2];
    Sol2 = Solve[Coeff2 == 0, k2];
    k2Val = k2 /. Sol2[[1]];
    
    (* CALCULATE EFFECTIVE MASS *)
    (* M = Rho * Volume *)
    RadiusExpr = a0 * (1 + k1Val*eps + k2Val*eps^2);
    MassExpr = Rho[eps] * (VolumeExact[RadiusExpr, L] / VolumeExact[a0, L]);
    MassSeries = Series[MassExpr, {eps, 0, 2}] // Normal;
    
    (* Extract 2PN Coefficient *)
    c2 = Coefficient[MassSeries, eps, 2];
    
    Return[c2];
];


(* 5. Run the Sweep over Tension Physics *)
Print["\nTesting Surface Tension Models (gamma ~ rho^k):"];
Print["k_gamma \t 2PN Coeff \t Target (-0.5) \t Error"];
Print["----------------------------------------------------------"];

(* We test k=0.5 (Wavefunction) to k=1.0 (Classical Energy) *)
Do[
    val = CheckVortexBridge[k];
    err = Abs[val - (-0.5)];
    match = If[err < 0.05, "MATCH! <<", ""];
    
    Print[
        NumberForm[k, {2, 1}], "\t\t", 
        NumberForm[val, {3, 4}], "\t", 
        "-0.5000", "\t\t",
        match
    ];
, {k, 0.5, 1.2, 0.1}];

Print["----------------------------------------------------------"];

(* 6. Interpret the Geometry *)
Print["\n--- Analysis ---"];
Print["If k=1.0 (Linear Tension) provides the match, it implies"];
Print["the Black Hole interface behaves like a classical fluid surface"];
Print["(Energy ~ Area * Density), stabilizing the Vortex Throat."];
