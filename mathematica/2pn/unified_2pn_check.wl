(* ---------------------------------------------------------------------- *)
(* Script: Unified_2PN_Check.wl *)
(* Purpose: The Final Test. Check if a Vortex-Supported Catenoid Bridge *)
(* with Aspect Ratio L/a ~ 1.85 reproduces GR 2PN Gravity. *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Unified Ontology Check"];
Print["   (Vortex Stiffness + Catenoid Geometry + n=5 Vacuum)"];
Print["======================================================="];

(* 1. Physics Parameters *)
n = 5;                  (* Stiff Superfluid *)
gammaExponent = 0.5;    (* Wavefunction Tension: gamma ~ sqrt(rho) *)

(* 2. Define Exact Geometry (Catenoid) *)
(* V = Pi * a^3 * [ (L/a) + 0.5 * Sinh[2*L/a] ] *)
VolumeExact[a_, L_] := Pi * a^3 * ( (L/a) + 0.5 * Sinh[2*L/a] );

(* 3. Define Thermodynamics *)
Rho[eps_] := (1 + (n - 1) * eps)^(1 / (n - 1));
TensionVal[eps_] := (Rho[eps])^gammaExponent;

(* 4. Define Forces *)

(* Closing Force: Quantum Surface Tension *)
(* F ~ Gamma * Circumference ~ a^1 *)
FClose[a_, eps_] := TensionVal[eps] * 2 * Pi * a;

(* Opening Force: The Critical Update *)
(* Model A: Photon Gas (Previous). Force ~ a^-1 *)
FOpenGas[a_] := C_gas / a;

(* Model B: Vortex Ring (New). Force ~ a^-3 *)
(* Derived from E ~ J^2 / a^2  => F = -dE/da ~ 1/a^3 *)
FOpenVortex[a_] := C_vortex / a^3;

(* 5. Solver Module *)

CheckModel[aspectRatio_, ForceModel_, ModelName_] := Module[
    {L, a0, kLinear, k1Val, k2Val, xSeries, BalanceEq, EqNormal, 
     Coeff1, Sol1, Coeff2, Sol2, RadiusExpr, MassExpr, MassSeries, c2, dev},
    
    (* Define Geometry *)
    a0 = 1.0; 
    L = aspectRatio * a0;
    
    (* Define Force Function based on selection *)
    (* We normalize constants so F_open(a0) = F_close(a0, 0) = 2 Pi *)
    
    FOpenFunc[a_] := If[ForceModel == "Gas",
        (2 * Pi * a0^2) / a,       (* Scales as 1/a *)
        (2 * Pi * a0^4) / a^3      (* Scales as 1/a^3 *)
    ];
    
    (* 1. Calculate 2nd Order Expansion *)
    xSeries = k1*eps + k2*eps^2;
    
    (* Force Balance: FOpen(a) - FClose(a) = 0 *)
    BalanceEq = Series[
        FOpenFunc[a0*(1+xSeries)] - FClose[a0*(1+xSeries), eps], 
        {eps, 0, 2}
    ];
    EqNormal = Normal[BalanceEq];
    
    (* Solve Order 1 *)
    Coeff1 = Coefficient[EqNormal, eps, 1];
    Sol1 = Solve[Coeff1 == 0, k1];
    k1Val = k1 /. Sol1[[1]];
    
    (* Solve Order 2 *)
    Coeff2 = Coefficient[EqNormal /. k1 -> k1Val, eps, 2];
    Sol2 = Solve[Coeff2 == 0, k2];
    k2Val = k2 /. Sol2[[1]];
    
    (* 2. Calculate Mass M(eps) *)
    (* M = Rho[eps] * VolumeExact[a[eps], L] *)
    RadiusExpr = a0 * (1 + k1Val*eps + k2Val*eps^2);
    MassExpr = Rho[eps] * (VolumeExact[RadiusExpr, L] / VolumeExact[a0, L]);
    MassSeries = Series[MassExpr, {eps, 0, 2}] // Normal;
    
    (* Extract 2PN Coefficient *)
    c2 = Coefficient[MassSeries, eps, 2];
    dev = c2 - (-0.5);
    
    Return[{c2, dev}];
];


(* 6. The Grand Finale Test *)

TargetRatio = 1.8475; (* The 'Magic' EM Ratio *)

Print["\n--- TEST 1: The Old Photon Gas Model (Re-check) ---"];
resGas = CheckModel[TargetRatio, "Gas", "Photon Gas"];
Print["Geometry L/a: ", TargetRatio];
Print["2PN Coefficient: ", NumberForm[resGas[[1]], {4, 4}]];
Print["Target (GR):      -0.5000"];
Print["Verdict: ", If[Abs[resGas[[2]]]<0.05, "SUCCESS", "FAIL (Too soft)"]];

Print["\n--- TEST 2: The New Vortex Ring Model ---"];
resVortex = CheckModel[TargetRatio, "Vortex", "Vortex Ring"];
Print["Geometry L/a: ", TargetRatio];
Print["2PN Coefficient: ", NumberForm[resVortex[[1]], {4, 4}]];
Print["Target (GR):      -0.5000"];
Print["Verdict: ", If[Abs[resVortex[[2]]]<0.1, "SUCCESS", "FAIL"]];


(* 7. Scan to see if Vortex matches ANY ratio *)
Print["\n--- Scanning Vortex Model over Aspect Ratios ---"];
Print["L/a \t 2PN Coeff \t Deviation"];
Do[
    res = CheckModel[r, "Vortex", "Vortex"];
    dev = res[[2]];
    marker = If[Abs[dev] < 0.05, "MATCH! <<", ""];
    Print[NumberForm[r, {2, 1}], "\t", NumberForm[res[[1]], {4, 3}], "\t", marker];
, {r, 1.0, 3.0, 0.2}];
