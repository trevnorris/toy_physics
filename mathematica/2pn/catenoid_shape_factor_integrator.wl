(* ---------------------------------------------------------------------- *)
(* Script: Catenoid_Shape_Factor_Integrator_Fixed.wl *)
(* Purpose: Calculate the EXACT 2PN coefficient for a Catenoid Bridge *)
(* by integrating the true shape factor V(a, L) and solving the *)
(* nonlinear force balance equation. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Catenoid Shape Factor Integration"];
Print["   (Exact Geometry Check - FIXED)"];
Print["======================================================="];

(* 1. Physics Parameters *)
n = 5;                  (* Stiff Superfluid *)
q = 1;                  (* Linear Mass *)
gammaExponent = 0.5;    (* "Wavefunction Tension": gamma ~ sqrt(rho) *)

(* 2. Define Exact Geometry *)

(* Profile: r(z) = a * Cosh[z/a] *)
(* Volume of throat from z = -L to z = +L *)
(* V = Integral[ Pi * r(z)^2, {z, -L, L} ] *)
(* Integral[ Cosh[u]^2 ] = u/2 + Sinh[2u]/4 *)
(* Result: V = Pi * a^3 * [ (L/a) + 0.5 * Sinh[2*L/a] ] *)

VolumeExact[a_, L_] := Pi * a^3 * ( (L/a) + 0.5 * Sinh[2*L/a] );

(* 3. Define Thermodynamics *)

(* Density vs Potential (Bernoulli) *)
Rho[eps_] := (1 + (n - 1) * eps)^(1 / (n - 1));

(* Surface Tension (Quantum Stiffness) *)
(* Renamed to TensionVal to avoid conflict with built-in Gamma *)
TensionVal[eps_] := (Rho[eps])^gammaExponent;

(* 4. Define Force Balance at Waist *)
(* The radius 'a' adjusts to balance Radiation Pressure vs Tension *)

(* Opening Force (Radiation Pressure) *)
(* P_rad ~ 1/Volume. Force ~ P * a^2 (Waist Area Scaling) *)
FOpen[a_, L_] := (1.0 / VolumeExact[a, L]) * a^2; 

(* Closing Force (Surface Tension) *)
(* F_close = Tension * Circumference = Tension * 2 * Pi * a *)
FClose[a_, eps_] := TensionVal[eps] * 2 * Pi * a;


(* 5. Solver Module *)

CalculateCoefficient[aspectRatio_] := Module[
    {L, a0, sol, aSol, massSeries, coeff, gap, k1Val, k2Val, Vol0, VolPrime0, S, kLinear},
    
    (* Define Geometry based on Aspect Ratio R = L/a0 *)
    a0 = 1.0; 
    L = aspectRatio * a0;
    
    (* 1. Solve for Linear Radius Scaling k (Breathing Exponent) *)
    
    (* Log-Derivative Matching at eps=0 *)
    (* Log(F_close) = Log(Tension) + Log(a) *)
    (* dLog/deps = (0.5 * dRho/deps) + k *)
    (* dRho/deps = 1 *)
    (* LHS = 0.5 + k *)
    
    (* Log(F_open) = 2 Log(a) - Log(V(a)) *)
    (* dLog/deps = 2k - (a*V'/V)*k = k(2-S) *)
    
    (* Calculate Shape Sensitivity S numerically at a0 *)
    Vol0 = VolumeExact[a0, L];
    VolPrime0 = D[VolumeExact[a, L], a] /. a -> a0;
    S = (a0 * VolPrime0) / Vol0;
    
    (* Equate: 0.5 + k = k(2-S)  =>  0.5 = k(1-S)  => k = 0.5 / (1-S) *)
    kLinear = 0.5 / (1 - S);
    
    (* 2. Calculate 2nd Order Expansion (The 2PN Term) *)
    
    (* Ansatz: a[eps] = a0 * (1 + k1*eps + k2*eps^2) *)
    xSeries = k1*eps + k2*eps^2;
    
    (* Expand Force Imbalance: FOpen - FClose = 0 *)
    BalanceEq = Series[
        FOpen[a0*(1+xSeries), L] - FClose[a0*(1+xSeries), eps], 
        {eps, 0, 2}
    ];
    EqNormal = Normal[BalanceEq];
    
    (* Solve Order 1 for k1 *)
    Coeff1 = Coefficient[EqNormal, eps, 1];
    Sol1 = Solve[Coeff1 == 0, k1];
    k1Val = k1 /. Sol1[[1]];
    
    (* Solve Order 2 for k2 *)
    Coeff2 = Coefficient[EqNormal /. k1 -> k1Val, eps, 2];
    Sol2 = Solve[Coeff2 == 0, k2];
    k2Val = k2 /. Sol2[[1]];
    
    (* 3. Calculate Mass M(eps) *)
    (* M = Rho[eps] * VolumeExact[a[eps], L] *)
    (* We normalize by M0/V0 *)
    
    RadiusExpr = a0 * (1 + k1Val*eps + k2Val*eps^2);
    MassExpr = Rho[eps] * (VolumeExact[RadiusExpr, L] / VolumeExact[a0, L]);
    
    MassSeries = Series[MassExpr, {eps, 0, 2}] // Normal;
    
    (* Extract 2PN Coefficient *)
    c2 = Coefficient[MassSeries, eps, 2];
    
    Return[{c2, k1Val, S}];
];


(* 6. Run Parameter Sweep *)

Print["\n--- Sweeping Aspect Ratio L/a ---"];
Print["Target 2PN Coefficient: -0.500"];
Print["(Corresponds to GR Beta = 1)"];
Print[""];
Print["L/a \t Shape(S) \t k (Breathing) \t 2PN Coeff \t Deviation"];
Print["------------------------------------------------------------------"];

Do[
    result = CalculateCoefficient[ratio];
    c2 = result[[1]];
    k = result[[2]];
    S = result[[3]];
    dev = c2 - (-0.5);
    
    Print[
        NumberForm[ratio, {3, 1}], "\t",
        NumberForm[S, {3, 2}], "\t",
        NumberForm[k, {3, 3}], "\t",
        NumberForm[c2, {3, 4}], "\t",
        If[Abs[dev]<0.005, "MATCH! <<<<<<", NumberForm[dev, {3, 4}]]
    ];
, {ratio, 0.1, 3.0, 0.1}];

Print["------------------------------------------------------------------"];

(* 7. Check Limits *)
Print["\nChecking 'Long Wormhole' (L/a = 10):"];
resultInf = CalculateCoefficient[10.0];
Print["L/a = 10.0 -> 2PN Coeff = ", resultInf[[1]]];
