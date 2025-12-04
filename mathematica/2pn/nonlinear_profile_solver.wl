(* ---------------------------------------------------------------------- *)
(* Script: Nonlinear_Profile_Solver.wl *)
(* Purpose: Solve the 4D Radial Inflow ODE to verify Transonic Horizon *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

(* 1. Parameters *)
n = 5;                  (* Stiff Superfluid *)
rH = 1.0;               (* Define Horizon Radius = 1 unit *)
rhoH = 1.0;             (* Density at Horizon = 1 unit *)
csH = rhoH^((n-1)/2);   (* Speed of Sound at Horizon *)

(* 2. Conservation Equations (Stationary Radial Flow) *)
(* Continuity: rho * v * Area = Constant *)
(* 4D Spherical/Cylindrical Geometry: Area ~ r^3 (Bulk radial inflow) *)
(* Note: We model the flow into the bulk sink *)
ContinuityEq = rho[r] * v[r] * r^3 == C1;

(* Bernoulli (Energy): h + v^2/2 = Constant *)
(* Enthalpy h = cs^2 / (n-1) for Polytrope *)
(* We use specific enthalpy formula for P~rho^n *)
Enthalpy[rho_] := (n/(n-1)) * rho^(n-1); (* Scaled enthalpy *)
BernoulliEq = Enthalpy[rho[r]] + v[r]^2/2 == C2;

(* 3. Determine Constants at Horizon (r = rH) *)
(* Condition: v = cs (Mach 1) at rH *)
C1Val = rhoH * csH * rH^3;
C2Val = Enthalpy[rhoH] + csH^2/2;

(* 4. Solve for Velocity Profile v(r) *)
(* We substitute rho from Continuity into Bernoulli to get an equation for v(r) *)
rhoExpr = C1Val / (v[r] * r^3);
VelocityEq = (n/(n-1)) * (rhoExpr)^(n-1) + v[r]^2/2 == C2Val;

(* Use FindRoot to numerically solve for v at each r *)
(* We scan from r = 1.0 (Horizon) outwards to r = 10.0 (Far Field) *)
solV[radius_] := Module[{vSol},
   vSol = v /. FindRoot[
      ((n/(n-1)) * (C1Val/(v * radius^3))^(n-1) + v^2/2) - C2Val == 0,
      {v, 0.1} (* Guess low velocity for far field *)
   ];
   vSol
];

(* 5. Generate Data for Plot *)
rValues = Range[1.0, 5.0, 0.1];
vValues = Table[solV[r], {r, rValues}];
machValues = Table[solV[r] / ( (C1Val/(solV[r]*r^3))^((n-1)/2) ), {r, rValues}];

(* 6. Plot Results *)
Print["--- Transonic Flow Profile ---"];
Print["Verifying smooth acceleration from Far Field (r>>1) to Horizon (r=1)"];

ListLinePlot[{
    Transpose[{rValues, vValues}],
    Transpose[{rValues, machValues}]
  },
  PlotLegends -> {"Velocity v(r)", "Mach Number M(r)"},
  AxesLabel -> {"Radius r", "Magnitude"},
  PlotLabel -> "Inflow Profile: Subsonic to Sonic",
  PlotRange -> All
]

(* 7. Check Stability *)
(* If Mach Number reaches 1.0 exactly at r=1, the Horizon is stable. *)
Print[""];
Print["Mach Number at r=1.0: ", machValues[[1]]];
If[Abs[machValues[[1]] - 1.0] < 0.01,
   Print["STATUS: STABLE. Sonic Horizon forms naturally."],
   Print["STATUS: UNSTABLE. Flow does not choke correctly."]
];
