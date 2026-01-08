(* ---------------------------------------------------------------------- *)
(* Script: hybrid_refractive_index_plot_v4.wl *)
(* Purpose: Robustly visualize Scalar vs Vector scaling laws *)
(* Fix: Removed underscores from variable names to prevent syntax errors *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

Print["======================================================="];
Print["   Paper VI: Hybrid Scaling Verification (Robust)"];
Print["======================================================="];

(* 1. Physics Constants *)
n = 5;
GM = 0.5;

(* Horizon Conditions (Sonic Point) *)
(* We set units such that rhoH=1, csH=1 at the horizon *)
rhoH = 1.0;
csH = 1.0; 
vH = 1.0; (* Sonic *)
rH = 1.0;

(* Derived Constants *)
(* K from c_s^2 = n K rho^(n-1) -> 1 = 5 K (1)^4 -> K = 1/5 *)
K = 1/n;
(* Continuity Flux: rho * v * r^2 = C1 *)
C1 = rhoH * vH * rH^2;

(* Bernoulli Constants (Energy Conservation) *)
(* Scalar Potential: Phi = -GM/r *)
Phi[r_] := -GM / r;
(* Enthalpy: h = c_s^2 / (n-1) = n K rho^(n-1) / (n-1) *)
Enthalpy[rho_] := (1/(n-1)) * rho^(n-1);

(* Total Energy (Horizon) *)
(* Renamed variables to remove underscores *)
EnergyVec = Enthalpy[rhoH] + vH^2/2;             (* Vector Only *)
EnergyHyb = Enthalpy[rhoH] + vH^2/2 + Phi[rH];   (* Hybrid *)

(* 2. Analytic Infinity Reference (Crucial for Normalization) *)
(* Vector Only: at inf, v->0. h -> EnergyVec *)
rhoInfVec = ((n-1) * EnergyVec)^(1/(n-1));

(* Hybrid: at inf, v->0, Phi->0. h -> EnergyHyb *)
rhoInfHyb = ((n-1) * EnergyHyb)^(1/(n-1));

(* Scalar Only: Hydrostatic. h + Phi = Const. *)
(* Matches Hybrid at infinity *)
rhoInfScalar = rhoInfHyb;

(* 3. Solvers (Solve for rho directly) *)

(* Scalar Only (Hydrostatic) *)
(* Enthalpy[rho] + Phi[r] = Enthalpy[rhoInf] *)
rhoScalar[r_] := ((n-1) * (Enthalpy[rhoInfScalar] - Phi[r]))^(1/(n-1));
NScalar[r_] := (rhoScalar[r]^2) / (rhoInfScalar^2); (* N ~ rho^2 for n=5 *)

(* Vector Only (Kinematic) *)
(* Enthalpy[rho] + 0.5 (C1/(rho r^2))^2 = EnergyVec *)
rhoVector[r_] := rho /. Quiet[FindRoot[
    Enthalpy[rho] + 0.5*(C1/(rho*r^2))^2 - EnergyVec == 0,
    {rho, 1.0}
]];
NVector[r_] := (rhoVector[r]^2) / (rhoInfVec^2);

(* Hybrid (Full) *)
(* Enthalpy[rho] + 0.5 (C1/(rho r^2))^2 + Phi[r] = EnergyHyb *)
rhoHybrid[r_] := rho /. Quiet[FindRoot[
    Enthalpy[rho] + 0.5*(C1/(rho*r^2))^2 + Phi[r] - EnergyHyb == 0,
    {rho, 1.0}
]];
NHybrid[r_] := (rhoHybrid[r]^2) / (rhoInfHyb^2);

(* 4. Generate Data & Ratio Check *)
r1 = 10.0; r2 = 20.0;
(* Calculate deviations N-1 *)
devV1 = NVector[r1] - 1; devV2 = NVector[r2] - 1;
devS1 = NScalar[r1] - 1; devS2 = NScalar[r2] - 1;

Print["\n--- Scaling Analysis (N-1) ---"];
Print["Scalar (r=10): ", devS1];
Print["Vector (r=10): ", devV1];
Print[""];
Print["Scalar Ratio (10/20): ", devS1/devS2, " (Target: ~2.0)"];
Print["Vector Ratio (10/20): ", devV1/devV2, " (Target: ~16.0)"];

(* 5. Plotting *)
(* Use LogLog to see the power laws *)
data = Table[{r, 
    Abs[NScalar[r]-1], 
    Abs[NVector[r]-1], 
    Abs[NHybrid[r]-1]
}, {r, 1.1, 50, 0.5}];

p1 = ListLogLogPlot[{
    data[[All, {1, 2}]], (* Scalar *)
    data[[All, {1, 3}]], (* Vector *)
    data[[All, {1, 4}]]  (* Hybrid *)
   },
   PlotStyle -> {
     {Blue, Dashed, Thickness[0.007]}, 
     {Red, Dashed, Thickness[0.007]}, 
     {Purple, Thickness[0.004]}
   },
   Joined -> True,
   PlotLegends -> {"Scalar (1/r)", "Vector (1/r^4)", "Hybrid"},
   AxesLabel -> {"Radius r/r_H", "Index Deviation N-1"},
   PlotLabel -> "The Hybrid Resolution",
   GridLines -> Automatic,
   PlotRange -> All
];
Print[p1];

(*"
Output:

--- Scaling Analysis (N-1) ---
Scalar (r=10): 0.09544511501033237
Vector (r=10): -0.00001924556455668025

Scalar Ratio (10/20): 1.9554879614778726 (Target: ~2.0)
Vector Ratio (10/20): 16.00043303798572 (Target: ~16.0)
"*)
