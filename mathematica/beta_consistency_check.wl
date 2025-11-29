(* β = 1.5 Consistency Check Against Superfluid Ontology *)

(* ================================================================== *)
(* C0: Setup                                                          *)
(* ================================================================== *)

ClearAll["Global`*"];

$Assumptions = {
  μ > 0, cs > 0, m > 0, ρ0 > 0, Q > 0, κg > 0,
  r > 0, r0 > 0, rThroat > 0, ξ > 0,
  Element[{μ, cs, m, ρ0, Q, κg, r, r0, rThroat, ξ}, Reals]
};

βVal = 3/2; (* UPDATED from 5/2 *)

Print["C0: Testing consistency of β = ", βVal];

(* ================================================================== *)
(* C1: Regime of validity - where σ(r) becomes O(1)                   *)
(* ================================================================== *)

σ[r_, β_] := β μ/(cs^2 r);

rBreakdownEq = Solve[σ[r, βVal] == 1, r][[1]];
rBreakdown = r /. rBreakdownEq;

Print["C1: σ(r) = 1 at r = ", rBreakdown];

rSchwarzschild = 2 μ/cs^2;

ratioToRs = Simplify[rBreakdown/rSchwarzschild];

Print["C1: r_breakdown = ", ratioToRs, " × r_Schwarzschild"];

(* ================================================================== *)
(* C10: Hydrodynamic origin of β - decomposition                      *)
(* ================================================================== *)

Print["\n============================================================="];
Print["C10: DECOMPOSITION OF β = 1.5"];
Print["============================================================="];

(* β = 1.5 must come from: 
   1. Direct density coupling: κ_ρ
   2. Added mass from flow: κ_add  
   3. Pressure-volume work: κ_PV *)

kappaRho = 1;   (* Standard density coupling *)
kappaAdd = 1/2; (* Standard sphere added mass *)

Print["  κ_ρ (Density)    = ", kappaRho];
Print["  κ_add (Added Mass) = ", kappaAdd];

kappaPV = βVal - kappaRho - kappaAdd;

Print["  -----------------------------"];
Print["  Remains for κ_PV = ", kappaPV];
Print["\nCONCLUSION: No exotic 4D geometry required. PV work is zero."];
Print["============================================================="];

(*"
Output:

C0: Testing consistency of β = 3/2
C1: σ(r) = 1 at r = ConditionalExpression[(3*μ)/(2*cs^2), Element[m, Reals] && Element[Q, Reals] && Element[r0, Reals] && Element[rThroat, Reals] && Element[κg, Reals] && Element[ξ, Reals] && Element[ρ0, Reals]]
C1: r_breakdown = 3/4 × r_Schwarzschild

=============================================================
C10: DECOMPOSITION OF β = 1.5
=============================================================
  κ_ρ (Density)    = 1
  κ_add (Added Mass) = 1/2
  -----------------------------
  Remains for κ_PV = 0

CONCLUSION: No exotic 4D geometry required. PV work is zero.
=============================================================
"*)
