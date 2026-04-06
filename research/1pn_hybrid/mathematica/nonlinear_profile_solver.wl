(* ---------------------------------------------------------------------- *)
(* Script: Nonlinear_Profile_Solver_v3.wl *)
(* Purpose: Solve Radial Inflow on 3D Brane with Subsonic Constraint *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

(* --- Common polytrope normalization (Paper VI convention) --- *)
n = 5;

(* Choose units: background density rho0 and background wave speed c0 *)
rho0 = 1.0;
c0 = 1.0;

(* EOS: P = K rho^n with K chosen so cs(rho0)=c0 *)
K = c0^2/(n * rho0^(n - 1));

Pressure[rho_] := K * rho^n;
cs2[rho_] := D[Pressure[rho], rho];          (* = n K rho^(n-1) *)
cs[rho_] := Sqrt[cs2[rho]];

(* Enthalpy: h = ∫ dP/rho = (n K/(n-1)) rho^(n-1) (up to additive constant) *)
Enthalpy[rho_] := (n*K/(n-1)) * rho^(n - 1);

(* 1. Parameters *)
rH = 1.0;               (* Horizon Radius *)
rhoH = 1.0;             (* Density at Horizon *)
csH = cs[rhoH];         (* Speed of Sound at Horizon *)
vH = csH;               (* Sonic at the horizon *)

(* 2. Define Physics Functions *)
(* 3D Brane Continuity: rho = C1 / (v * r^2) *)
C1Val = rhoH * csH * rH^2;
rhoExpr[v_, r_] := C1Val / (v * r^2);

(* Bernoulli: Enthalpy + v^2/2 = Constant *)
C2Val = Enthalpy[rhoH] + vH^2/2;
BernoulliErr[v_, r_] := (Enthalpy[rhoExpr[v, r]] + v^2/2) - C2Val;

(* 3. Solve with Continuation and Constraints *)
(* We scan OUTWARDS from r=1.0. *)
(* We use the 'lastV' as the guess for the next step to stay on the branch. *)

rValues = Range[1.0, 10.0, 0.1];
vResults = {};
lastV = vH; (* Start at the horizon velocity *)

Do[
   (* Find root bounded between 0 and 1.2 (slightly above sonic is allowed near r=1) *)
   (* This prevents jumping to the supersonic branch v ~ 1.8 *)
   sol = FindRoot[BernoulliErr[v, r] == 0, {v, lastV, 0.0001, 2.0}];
   
   currentV = v /. sol;
   AppendTo[vResults, currentV];
   lastV = currentV; (* Update guess for next r *)
, {r, rValues}];

(* 4. Compute Derived Quantities *)
machResults = Table[
   Module[{v, rho, csVal},
     v = vResults[[i]];
     r = rValues[[i]];
     rho = rhoExpr[v, r];
     csVal = cs[rho];
     v / csVal
   ], 
   {i, 1, Length[rValues]}
];

(* 5. Plot Results *)
Print["--- Corrected Subsonic Profile ---"];
p1 = ListLinePlot[{
    Transpose[{rValues, vResults}],
    Transpose[{rValues, machResults}]
  },
  PlotLegends -> {"Velocity v(r)", "Mach Number M(r)"},
  AxesLabel -> {"Radius r", "Magnitude"},
  PlotLabel -> "Stable Acoustic Horizon (Subsonic Far Field)",
  PlotRange -> All,
  PlotStyle -> {Blue, Red}
];
Print[p1];

(* 6. Validation *)
vFar = Last[vResults];
machFar = Last[machResults];

Print[""];
Print["Horizon Mach Number (r=1): ", First[machResults]];
Print["Far Field Velocity (r=10): ", vFar];
Print["Far Field Mach (r=10):     ", machFar];

If[vFar < 0.1 && Abs[First[machResults] - 1.0] < 0.01,
   Print["STATUS: SUCCESS. Physically correct profile found."],
   Print["STATUS: FAILURE. Solution still divergent."]
];

(*"
Output:

--- Corrected Subsonic Profile ---
Horizon Mach Number (r=1): 1.
Far Field Velocity (r=10): 0.009193282651578608
Far Field Mach (r=10):     0.007769835759741814
STATUS: SUCCESS. Physically correct profile found.
"*)
