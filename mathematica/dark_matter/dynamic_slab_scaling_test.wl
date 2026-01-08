(* ================================================================= *)
(* DYNAMIC SLAB TEST: BERNOULLI BACK-REACTION & TULLY-FISHER SCALING *)
(* ================================================================= *)
(* Hypothesis: The slab thickness H is not constant. It is defined   *)
(* by a critical pressure isobar in the superfluid bulk.             *)
(* Massive galaxies -> Faster Flow -> Lower Pressure -> Larger H_eff *)
(* ================================================================= *)

ClearAll["Global`*"]

(* 1. SETUP: The Superfluid Sink Model *)
(* 0th order: The flow into a sink in 3D is spherically symmetric *)
(* Velocity u(R) ~ flux / R^2 *)

(* u = Velocity Magnitude, M = Mass parameter (Source Strength) *)
u[R_, M_] := M / R^2

(* 2. PHYSICS: Bernoulli Pressure Field *)
(* P_static = P0 *)
(* P_dynamic(R) = P0 - 1/2 * rho * u(R)^2 *)
(* We look for the "Vacuum Boundary" where P drops to a critical P_c *)

(* The condition for the boundary R_bound: *)
(* 1/2 * rho * u(R_bound)^2 = DeltaP_critical *)

(* Rearranging for R_bound (which acts as our H_eff): *)
(* (M / R^2)^2 ~ Constant *)
(* M^2 / R^4 ~ Constant *)
(* R ~ Sqrt[M] *)

(* Let's verify this scaling numerically and plot the Rotation Curve impact *)

(* 3. NUMERICAL SCALING STUDY *)
Mlist = Table[10^x, {x, 9, 12, 0.5}]; (* Masses from 10^9 to 10^12 *)
Pcrit = 10^-5; (* Arbitrary critical pressure drop threshold *)

(* Calculate H_eff for each Mass *)
(* Solve: (M / H^2)^2 == Pcrit *)
getHeff[mass_] := Module[{h},
  h /. FindRoot[(mass / h^2)^2 == Pcrit, {h, 1.0}]
];

results = Table[
   Module[{h, vslab, vtf},
    h = getHeff[m];
    
    (* The Static Slab Prediction: v^2 = GM / H *)
    (* We use our dynamic H(m) here *)
    vslab = Sqrt[m / h]; 
    
    {m, h, vslab}
   ],
   {m, Mlist}
];

(* 4. ANALYSIS & PLOTTING *)

(* Log-Log Plot of Mass vs Velocity *)
dataPlot = ListLogLogPlot[
   results[[All, {3, 1}]], (* Plot x=Velocity, y=Mass *)
   PlotStyle -> {Red, PointSize[0.02]},
   Frame -> True,
   FrameLabel -> {"Asymptotic Velocity (v)", "Galaxy Mass (M)"},
   PlotLabel -> "Dynamic Slab Scaling (Tully-Fisher Test)",
   GridLines -> Automatic
];

(* Fit the slope: Log[M] = slope * Log[v] + b *)
logData = Log10 /@ results[[All, {3, 1}]];
lm = LinearModelFit[logData, x, x];
slope = lm["ParameterTableEntries"][[2, 1]];

Print["\n================================================="];
Print["SCALING RESULTS (Bernoulli Deformable Slab)"];
Print["================================================="];
Print["If H is constant, slope M vs v should be 2.0 (Newtonian)."];
Print["Observed Tully-Fisher slope is ~4.0."];
Print["\nCALCULATED SLOPE from Dynamic Model: ", NumberForm[slope, 3]];

If[Abs[slope - 4.0] < 0.5, 
 Print["SUCCESS: The hydrodynamic back-reaction recovers the Tully-Fisher relation!"],
 Print["FAIL: The scaling does not match 4.0."]
];

Print["\nInterpretation:"];
Print["Massive galaxies create strong suction."];
Print["This extends the 'low pressure' zone (H_eff) further out."];
Print["Since H scales as Sqrt[M], the velocity v^2 ~ M/H becomes v^2 ~ M/Sqrt[M] ~ Sqrt[M]."];
Print["Therefore v^4 ~ M."];

Show[dataPlot, 
 Plot[Exp[lm[Log[x]]], {x, 10, 1000}, PlotStyle -> Blue]
]

(*"
Output:


=================================================
SCALING RESULTS (Bernoulli Deformable Slab)
=================================================
If H is constant, slope M vs v should be 2.0 (Newtonian).
Observed Tully-Fisher slope is ~4.0.

CALCULATED SLOPE from Dynamic Model: NumberForm[4., 3]
SUCCESS: The hydrodynamic back-reaction recovers the Tully-Fisher relation!

Interpretation:
Massive galaxies create strong suction.
This extends the 'low pressure' zone (H_eff) further out.
Since H scales as Sqrt[M], the velocity v^2 ~ M/H becomes v^2 ~ M/Sqrt[M] ~ Sqrt[M].
Therefore v^4 ~ M.
"*)
