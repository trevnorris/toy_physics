(* ::Package:: *)

(*
Madelung / Hydrodynamic Reduction of the Gross–Pitaevskii Equation (GPE)
======================================================================

What this script gives you (Paper 1 needs):
  1) From the GPE, via ψ = √ρ e^{iθ}, derive:
       (a) Continuity:   ∂t ρ + ∇·(ρ v) = 0
       (b) Phase/Bernoulli:  ħ ∂t θ + (m/2)|v|^2 + V + g ρ + Q = 0
           where Q = -(ħ^2/2m) (∇^2√ρ)/√ρ  (quantum-pressure / quantum potential)
  2) A compact “potential-flow” momentum form from ∇(phase equation):
       ∂t v + ∇(|v|^2/2) = -(1/m)∇(V + gρ + Q)
     and the standard identity connecting it to Euler convective form.
  3) A generic compressible Euler vorticity-transport identity:
       ∂t ω = ∇×(v×ω) + ∇(1/ρ)×∇p
     and barotropic reduction (p=p(ρ) ⇒ baroclinic term = 0).

No file outputs by default:
  • The script prints canonical equations to the CLI.
  • If you ever want TeX strings (still no files), set PrintTeX = True.

Usage:
  wolframscript -file madelung_derivation_noexport.wl

Notes for Paper 1 writeup:
  • These manipulations are exact wherever ρ>0 (i.e., away from vortex cores).
  • In GP, the interaction term gives a barotropic pressure (p ∝ ρ^2), so the
    baroclinic term vanishes away from singular cores.
  • Any “non-ideal” effects you add in the solver (e.g., gated diffusion) can be
    treated as additional terms on the RHS of the induction-like / vorticity equation.

*)

(* ----------------------------- *)
(* 0. Controls & Setup           *)
(* ----------------------------- *)

ClearAll["Global`*"];

PrintTeX = False;   (* if True, prints TeXForm strings to CLI (no file writes) *)
Verbosee = True;     (* if False, prints only the final canonical equations *)

coords = {x, y, z};
tvar = t;

(* Parameters (real) *)
$Assumptions = Element[{x, y, z, t, \[HBar], m, g}, Reals] && m > 0 && \[HBar] > 0;

(* Fields (real-valued functions) *)
rhoF   = \[Rho][x, y, z, t];
thetaF = \[Theta][x, y, z, t];
Vext   = V[x, y, z, t];

(* Helpful assumptions for simplification *)
assumpFields = $Assumptions && rhoF > 0 &&
  Element[{\[Rho][x, y, z, t], \[Theta][x, y, z, t], V[x, y, z, t]}, Reals];

psi = Sqrt[rhoF] Exp[I thetaF];

If[Verbosee,
  Print["--- Madelung ansatz ---"]; 
  Print["psi(x,t) = Sqrt[rho] * Exp[I theta]"];
];

(* ----------------------------- *)
(* 1. Start from the GPE          *)
(* ----------------------------- *)

(* Standard dimensional GPE:
     i ħ ∂t ψ = -(ħ^2/2m) ∇^2 ψ + V ψ + g |ψ|^2 ψ
   Under Madelung, |ψ|^2 = ρ.
*)

eqGPE = I \[HBar] D[psi, tvar] == - (\[HBar]^2/(2 m)) Laplacian[psi, coords] + Vext psi + g rhoF psi;

If[Verbosee,
  Print["\n--- Starting equation (GPE) ---"]; 
  Print[eqGPE];
];

(* Residual of (GPE)/psi = 0, then split into Re/Im *)
res = (I \[HBar] D[psi, tvar] + (\[HBar]^2/(2 m)) Laplacian[psi, coords] - Vext psi - g rhoF psi)/psi;

resRe = FullSimplify[ComplexExpand[Re[res], coords ~Join~ {t, \[HBar], m, g}], assumpFields];
resIm = FullSimplify[ComplexExpand[Im[res], coords ~Join~ {t, \[HBar], m, g}], assumpFields];

If[Verbosee,
  Print["\n--- Residual split (Re, Im) of (GPE)/psi = 0 ---"]; 
  Print["Re(res) = ", resRe];
  Print["Im(res) = ", resIm];
];

(* ----------------------------- *)
(* 2. Continuity equation         *)
(* ----------------------------- *)

(* Imag part → continuity up to an overall factor/sign *)
contCandidate = FullSimplify[(2 rhoF/\[HBar]) resIm, assumpFields];

(* Canonical form: ∂t ρ + (ħ/m) ∇·(ρ ∇θ) = 0 *)
contCanonical = D[rhoF, tvar] + (\[HBar]/m) Div[rhoF Grad[thetaF, coords], coords] == 0;

checkContPlus  = FullSimplify[contCandidate - LeftHandSide[contCanonical], assumpFields];
checkContMinus = FullSimplify[contCandidate + LeftHandSide[contCanonical], assumpFields];

If[Verbosee,
  Print["\n--- Continuity equation ---"]; 
  Print["Candidate (from Im part): ", contCandidate, " == 0"]; 
  Print["Canonical: ", contCanonical];
  Print["Check (candidate - canonical LHS) → ", checkContPlus];
  Print["Check (candidate + canonical LHS) → ", checkContMinus];
];

(* Define velocity field v = (ħ/m) ∇θ *)
v = (\[HBar]/m) Grad[thetaF, coords];

contInV = FullSimplify[
  D[rhoF, tvar] + Div[rhoF v, coords] == 0,
  assumpFields
];

If[Verbosee,
  Print["\nContinuity in v-form: "]; 
  Print[contInV];
];

(* ----------------------------- *)
(* 3. Phase / Bernoulli equation  *)
(* ----------------------------- *)

(* Quantum potential Q = -(ħ^2/2m) (∇^2√ρ)/√ρ *)
Q = - (\[HBar]^2/(2 m)) (Laplacian[Sqrt[rhoF], coords]/Sqrt[rhoF]);

(* Canonical phase equation:
     ħ ∂t θ + (ħ^2/2m)|∇θ|^2 + V + gρ + Q = 0
*)
phaseCanonical = \[HBar] D[thetaF, tvar] + (\[HBar]^2/(2 m)) (Grad[thetaF, coords].Grad[thetaF, coords]) + Vext + g rhoF + Q == 0;

phaseCandidate = FullSimplify[resRe, assumpFields];
checkPhasePlus  = FullSimplify[phaseCandidate - LeftHandSide[phaseCanonical], assumpFields];
checkPhaseMinus = FullSimplify[phaseCandidate + LeftHandSide[phaseCanonical], assumpFields];

If[Verbosee,
  Print["\n--- Phase/Bernoulli equation ---"]; 
  Print["Candidate (from Re part): ", phaseCandidate, " == 0"]; 
  Print["Canonical: ", phaseCanonical];
  Print["Check (candidate - canonical LHS) → ", checkPhasePlus];
  Print["Check (candidate + canonical LHS) → ", checkPhaseMinus];
];

(* Rewrite phase equation in terms of v:
   |∇θ|^2 = (m^2/ħ^2) |v|^2 ⇒ (ħ^2/2m)|∇θ|^2 = (m/2)|v|^2
*)
phaseInV = FullSimplify[
  \[HBar] D[thetaF, tvar] + (m/2) (v.v) + Vext + g rhoF + Q == 0,
  assumpFields
];

If[Verbosee,
  Print["\nPhase equation in v-form: "]; 
  Print[phaseInV];
];

(* ----------------------------- *)
(* 4. Momentum equation (compact) *)
(* ----------------------------- *)

(* Take ∇ of the v-form Bernoulli equation (formal):
     ∂t v + ∇(|v|^2/2) = -(1/m)∇(V + gρ + Q)
   Note: converting ∇(|v|^2/2) to (v·∇)v uses
     ∇(|v|^2/2) = (v·∇)v + v×(∇×v)
   Away from singular cores (or in irrotational regions), v×ω may be dropped.
*)

momCompact = D[v, tvar] + Grad[(v.v)/2, coords] == -(1/m) Grad[Vext + g rhoF + Q, coords];

If[Verbosee,
  Print["\n--- Momentum equation (compact potential-flow form) ---"]; 
  Print[momCompact];
];

(* ----------------------------- *)
(* 5. Generic vorticity transport  *)
(* ----------------------------- *)

(* Generic compressible Euler (for identity):
     ∂t v + (v·∇)v = -(1/ρ)∇p - ∇Φ
   Then ω = ∇×v obeys:
     ∂t ω = ∇×(v×ω) + ∇(1/ρ)×∇p
   If p=p(ρ) (barotropic), the baroclinic term vanishes.
*)

rhoE = \[Rho]E[x, y, z, t];
PhiE = \[CapitalPhi][x, y, z, t];

vE = {u[x, y, z, t], vcomp[x, y, z, t], w[x, y, z, t]};
omegaE = Curl[vE, coords];

pE = p[rhoE];
advE = vE[[1]] D[vE, x] + vE[[2]] D[vE, y] + vE[[3]] D[vE, z];

eulerEq = D[vE, tvar] + advE == -(1/rhoE) Grad[pE, coords] - Grad[PhiE, coords];
baroclinic = Cross[Grad[1/rhoE, coords], Grad[pE, coords]];
vortEqGeneral = D[omegaE, tvar] == Curl[Cross[vE, omegaE], coords] + baroclinic;

baroclinicReduced = FullSimplify[baroclinic, $Assumptions && rhoE > 0];
vortEqBarotropic = D[omegaE, tvar] == Curl[Cross[vE, omegaE], coords];

If[Verbosee,
  Print["\n--- Generic Euler equation (for vorticity identity) ---"]; 
  Print[eulerEq];
  Print["\n--- Vorticity transport (general form) ---"]; 
  Print[vortEqGeneral];
  Print["\nBaroclinic term ∇(1/ρ)×∇p with p=p(ρ) simplifies to: "]; 
  Print[baroclinicReduced];
  Print["(Expected: {0,0,0} since ∇p ∥ ∇ρ for p=p(ρ))"]; 
  Print["\n--- Vorticity transport (barotropic / induction-like) ---"]; 
  Print[vortEqBarotropic];
];

(* ----------------------------- *)
(* 6. Final canonical outputs     *)
(* ----------------------------- *)

Print["\n============================"]; 
Print["FINAL CANONICAL EQUATIONS"]; 
Print["============================"]; 

Print["\n(1) Continuity (v = (ħ/m)∇θ):"]; 
Print[contInV];

Print["\n(2) Phase/Bernoulli:"]; 
Print[phaseCanonical];

Print["\n(3) Quantum potential term Q:"]; 
Print[HoldForm[Q] == Q];

Print["\n(4) Phase in v-form:"]; 
Print[phaseInV];

Print["\n(5) Momentum (compact potential-flow form):"]; 
Print[momCompact];

Print["\n(6) Vorticity transport (general):"]; 
Print[vortEqGeneral];

Print["\n(7) Vorticity transport (barotropic / induction-like):"]; 
Print[vortEqBarotropic];

If[PrintTeX,
  Print["\n--- TeX strings (no files written) ---"]; 
  Print["Continuity TeX: ", ToString[TeXForm[contInV]]];
  Print["Phase TeX: ", ToString[TeXForm[phaseCanonical]]];
  Print["Momentum TeX: ", ToString[TeXForm[momCompact]]];
  Print["Vorticity TeX: ", ToString[TeXForm[vortEqBarotropic]]];
];

Print["\nDone."];

(*
Output:

--- Madelung ansatz ---
psi(x,t) = Sqrt[rho] * Exp[I theta]

--- Starting equation (GPE) ---
I*ℏ*(I*E^(I*θ[x, y, z, t])*Sqrt[ρ[x, y, z, t]]*Derivative[0, 0, 0, 1][θ][x, y, z, t] + (E^(I*θ[x, y, z, t])*Derivative[0, 0, 0, 1][ρ][x, y, z, t])/(2*Sqrt[ρ[x, y, z, t]])) == E^(I*θ[x, y, z, t])*V[x, y, z, t]*Sqrt[ρ[x, y, z, t]] + E^(I*θ[x, y, z, t])*g*ρ[x, y, z, t]^(3/2) - (ℏ^2*(-(E^(I*θ[x, y, z, t])*Sqrt[ρ[x, y, z, t]]*Derivative[0, 0, 1, 0][θ][x, y, z, t]^2) + (I*E^(I*θ[x, y, z, t])*Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t])/Sqrt[ρ[x, y, z, t]] - (E^(I*θ[x, y, z, t])*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2)/(4*ρ[x, y, z, t]^(3/2)) + I*E^(I*θ[x, y, z, t])*Sqrt[ρ[x, y, z, t]]*Derivative[0, 0, 2, 0][θ][x, y, z, t] + (E^(I*θ[x, y, z, t])*Derivative[0, 0, 2, 0][ρ][x, y, z, t])/(2*Sqrt[ρ[x, y, z, t]]) - E^(I*θ[x, y, z, t])*Sqrt[ρ[x, y, z, t]]*Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + (I*E^(I*θ[x, y, z, t])*Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t])/Sqrt[ρ[x, y, z, t]] - (E^(I*θ[x, y, z, t])*Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2)/(4*ρ[x, y, z, t]^(3/2)) + I*E^(I*θ[x, y, z, t])*Sqrt[ρ[x, y, z, t]]*Derivative[0, 2, 0, 0][θ][x, y, z, t] + (E^(I*θ[x, y, z, t])*Derivative[0, 2, 0, 0][ρ][x, y, z, t])/(2*Sqrt[ρ[x, y, z, t]]) - E^(I*θ[x, y, z, t])*Sqrt[ρ[x, y, z, t]]*Derivative[1, 0, 0, 0][θ][x, y, z, t]^2 + (I*E^(I*θ[x, y, z, t])*Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/Sqrt[ρ[x, y, z, t]] - (E^(I*θ[x, y, z, t])*Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2)/(4*ρ[x, y, z, t]^(3/2)) + I*E^(I*θ[x, y, z, t])*Sqrt[ρ[x, y, z, t]]*Derivative[2, 0, 0, 0][θ][x, y, z, t] + (E^(I*θ[x, y, z, t])*Derivative[2, 0, 0, 0][ρ][x, y, z, t])/(2*Sqrt[ρ[x, y, z, t]])))/(2*m)

--- Residual split (Re, Im) of (GPE)/psi = 0 ---
Re(res) = -V[x, y, z, t] - (8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) - 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t]))/(8*m*ρ[x, y, z, t]^2)
Im(res) = (ℏ*(m*Derivative[0, 0, 0, 1][ρ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t]))))/(2*m*ρ[x, y, z, t])

--- Continuity equation ---
Candidate (from Im part): Derivative[0, 0, 0, 1][ρ][x, y, z, t] + (ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t])))/m == 0
Canonical: Derivative[0, 0, 0, 1][ρ][x, y, z, t] + (ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*Derivative[2, 0, 0, 0][θ][x, y, z, t]))/m == 0
Check (candidate - canonical LHS) → -LeftHandSide[m*Derivative[0, 0, 0, 1][ρ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t])) == 0] + Derivative[0, 0, 0, 1][ρ][x, y, z, t] + (ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t])))/m
Check (candidate + canonical LHS) → LeftHandSide[m*Derivative[0, 0, 0, 1][ρ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t])) == 0] + Derivative[0, 0, 0, 1][ρ][x, y, z, t] + (ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t])))/m

Continuity in v-form:
m*Derivative[0, 0, 0, 1][ρ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t])) == 0

--- Phase/Bernoulli equation ---
Candidate (from Re part): -V[x, y, z, t] - (8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) - 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t]))/(8*m*ρ[x, y, z, t]^2) == 0
Canonical: V[x, y, z, t] + g*ρ[x, y, z, t] + ℏ*Derivative[0, 0, 0, 1][θ][x, y, z, t] + (ℏ^2*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2))/(2*m) - (ℏ^2*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]) == 0
Check (candidate - canonical LHS) → -LeftHandSide[8*m*V[x, y, z, t]*ρ[x, y, z, t]^2 + 8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) == 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t])] - V[x, y, z, t] - (8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) - 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t]))/(8*m*ρ[x, y, z, t]^2)
Check (candidate + canonical LHS) → LeftHandSide[8*m*V[x, y, z, t]*ρ[x, y, z, t]^2 + 8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) == 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t])] - V[x, y, z, t] - (8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) - 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t]))/(8*m*ρ[x, y, z, t]^2)

Phase equation in v-form:
8*m*V[x, y, z, t]*ρ[x, y, z, t]^2 + 8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) == 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t])

--- Momentum equation (compact potential-flow form) ---
{(ℏ*Derivative[1, 0, 0, 1][θ][x, y, z, t])/m + ((2*ℏ^2*Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[1, 0, 1, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[1, 1, 0, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[2, 0, 0, 0][θ][x, y, z, t])/m^2)/2, (ℏ*Derivative[0, 1, 0, 1][θ][x, y, z, t])/m + ((2*ℏ^2*Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 1, 1, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 2, 0, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 1, 0, 0][θ][x, y, z, t])/m^2)/2, (ℏ*Derivative[0, 0, 1, 1][θ][x, y, z, t])/m + ((2*ℏ^2*Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 2, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 1, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 1, 0][θ][x, y, z, t])/m^2)/2} == {-((Derivative[1, 0, 0, 0][V][x, y, z, t] + g*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + (ℏ^2*Derivative[1, 0, 0, 0][ρ][x, y, z, t]*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(4*m*ρ[x, y, z, t]^(3/2)) - (ℏ^2*((3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 2, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 2, 0, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + (3*Derivative[1, 0, 0, 0][ρ][x, y, z, t]^3)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[1, 0, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) + Derivative[1, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - (Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[1, 1, 0, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) + Derivative[1, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - (3*Derivative[1, 0, 0, 0][ρ][x, y, z, t]*Derivative[2, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[3, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]))/m), -((Derivative[0, 1, 0, 0][V][x, y, z, t] + g*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + (ℏ^2*Derivative[0, 1, 0, 0][ρ][x, y, z, t]*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(4*m*ρ[x, y, z, t]^(3/2)) - (ℏ^2*((3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2*Derivative[0, 1, 0, 0][ρ][x, y, z, t])/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 2, 0][ρ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]^3)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 1, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) + Derivative[0, 1, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[0, 2, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 3, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) + (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[1, 0, 0, 0][ρ][x, y, z, t]*Derivative[1, 1, 0, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) - (Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[2, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 1, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]))/m), -((Derivative[0, 0, 1, 0][V][x, y, z, t] + g*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + (ℏ^2*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(4*m*ρ[x, y, z, t]^(3/2)) - (ℏ^2*((3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^3)/(8*ρ[x, y, z, t]^(5/2)) - (3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 0, 2, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 0, 3, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) + (3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[0, 1, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 2, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 1, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) + (3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[1, 0, 0, 0][ρ][x, y, z, t]*Derivative[1, 0, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[2, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 1, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]))/m)}

--- Generic Euler equation (for vorticity identity) ---
{Derivative[0, 0, 0, 1][u][x, y, z, t] + w[x, y, z, t]*Derivative[0, 0, 1, 0][u][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 1, 0, 0][u][x, y, z, t] + u[x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t], Derivative[0, 0, 0, 1][vcomp][x, y, z, t] + w[x, y, z, t]*Derivative[0, 0, 1, 0][vcomp][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] + u[x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t], Derivative[0, 0, 0, 1][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t]} == {-((Derivative[1][p][ρE[x, y, z, t]]*Derivative[1, 0, 0, 0][ρE][x, y, z, t])/ρE[x, y, z, t]) - Derivative[1, 0, 0, 0][Φ][x, y, z, t], -((Derivative[1][p][ρE[x, y, z, t]]*Derivative[0, 1, 0, 0][ρE][x, y, z, t])/ρE[x, y, z, t]) - Derivative[0, 1, 0, 0][Φ][x, y, z, t], -((Derivative[1][p][ρE[x, y, z, t]]*Derivative[0, 0, 1, 0][ρE][x, y, z, t])/ρE[x, y, z, t]) - Derivative[0, 0, 1, 0][Φ][x, y, z, t]}

--- Vorticity transport (general form) ---
{-Derivative[0, 0, 1, 1][vcomp][x, y, z, t] + Derivative[0, 1, 0, 1][w][x, y, z, t], Derivative[0, 0, 1, 1][u][x, y, z, t] - Derivative[1, 0, 0, 1][w][x, y, z, t], -Derivative[0, 1, 0, 1][u][x, y, z, t] + Derivative[1, 0, 0, 1][vcomp][x, y, z, t]} == {Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 0, 2, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][vcomp][x, y, z, t] - w[x, y, z, t]*Derivative[0, 1, 1, 0][w][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][w][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t], -(Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t]) - w[x, y, z, t]*Derivative[0, 0, 2, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - u[x, y, z, t]*Derivative[1, 0, 1, 0][u][x, y, z, t] + w[x, y, z, t]*Derivative[1, 0, 1, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[2, 0, 0, 0][w][x, y, z, t], Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - w[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] + u[x, y, z, t]*Derivative[1, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[2, 0, 0, 0][vcomp][x, y, z, t]}

Baroclinic term ∇(1/ρ)×∇p with p=p(ρ) simplifies to:
{0, 0, 0}
(Expected: {0,0,0} since ∇p ∥ ∇ρ for p=p(ρ))

--- Vorticity transport (barotropic / induction-like) ---
{-Derivative[0, 0, 1, 1][vcomp][x, y, z, t] + Derivative[0, 1, 0, 1][w][x, y, z, t], Derivative[0, 0, 1, 1][u][x, y, z, t] - Derivative[1, 0, 0, 1][w][x, y, z, t], -Derivative[0, 1, 0, 1][u][x, y, z, t] + Derivative[1, 0, 0, 1][vcomp][x, y, z, t]} == {Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 0, 2, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][vcomp][x, y, z, t] - w[x, y, z, t]*Derivative[0, 1, 1, 0][w][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][w][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t], -(Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t]) - w[x, y, z, t]*Derivative[0, 0, 2, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - u[x, y, z, t]*Derivative[1, 0, 1, 0][u][x, y, z, t] + w[x, y, z, t]*Derivative[1, 0, 1, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[2, 0, 0, 0][w][x, y, z, t], Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - w[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] + u[x, y, z, t]*Derivative[1, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[2, 0, 0, 0][vcomp][x, y, z, t]}

============================
FINAL CANONICAL EQUATIONS
============================

(1) Continuity (v = (ħ/m)∇θ):
m*Derivative[0, 0, 0, 1][ρ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][θ][x, y, z, t] + Derivative[0, 2, 0, 0][θ][x, y, z, t] + Derivative[2, 0, 0, 0][θ][x, y, z, t])) == 0

(2) Phase/Bernoulli:
V[x, y, z, t] + g*ρ[x, y, z, t] + ℏ*Derivative[0, 0, 0, 1][θ][x, y, z, t] + (ℏ^2*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2))/(2*m) - (ℏ^2*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]) == 0

(3) Quantum potential term Q:
HoldForm[Q] == -1/2*(ℏ^2*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(m*Sqrt[ρ[x, y, z, t]])

(4) Phase in v-form:
8*m*V[x, y, z, t]*ρ[x, y, z, t]^2 + 8*g*m*ρ[x, y, z, t]^3 + 4*ℏ*ρ[x, y, z, t]^2*(2*m*Derivative[0, 0, 0, 1][θ][x, y, z, t] + ℏ*(Derivative[0, 0, 1, 0][θ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][θ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][θ][x, y, z, t]^2)) + ℏ^2*(Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2 + Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2 + Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2) == 2*ℏ^2*ρ[x, y, z, t]*(Derivative[0, 0, 2, 0][ρ][x, y, z, t] + Derivative[0, 2, 0, 0][ρ][x, y, z, t] + Derivative[2, 0, 0, 0][ρ][x, y, z, t])

(5) Momentum (compact potential-flow form):
{(ℏ*Derivative[1, 0, 0, 1][θ][x, y, z, t])/m + ((2*ℏ^2*Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[1, 0, 1, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[1, 1, 0, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[2, 0, 0, 0][θ][x, y, z, t])/m^2)/2, (ℏ*Derivative[0, 1, 0, 1][θ][x, y, z, t])/m + ((2*ℏ^2*Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 1, 1, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 2, 0, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 1, 0, 0][θ][x, y, z, t])/m^2)/2, (ℏ*Derivative[0, 0, 1, 1][θ][x, y, z, t])/m + ((2*ℏ^2*Derivative[0, 0, 1, 0][θ][x, y, z, t]*Derivative[0, 0, 2, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[0, 1, 0, 0][θ][x, y, z, t]*Derivative[0, 1, 1, 0][θ][x, y, z, t])/m^2 + (2*ℏ^2*Derivative[1, 0, 0, 0][θ][x, y, z, t]*Derivative[1, 0, 1, 0][θ][x, y, z, t])/m^2)/2} == {-((Derivative[1, 0, 0, 0][V][x, y, z, t] + g*Derivative[1, 0, 0, 0][ρ][x, y, z, t] + (ℏ^2*Derivative[1, 0, 0, 0][ρ][x, y, z, t]*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(4*m*ρ[x, y, z, t]^(3/2)) - (ℏ^2*((3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 2, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 2, 0, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + (3*Derivative[1, 0, 0, 0][ρ][x, y, z, t]^3)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[1, 0, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) + Derivative[1, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - (Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[1, 1, 0, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) + Derivative[1, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - (3*Derivative[1, 0, 0, 0][ρ][x, y, z, t]*Derivative[2, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[3, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]))/m), -((Derivative[0, 1, 0, 0][V][x, y, z, t] + g*Derivative[0, 1, 0, 0][ρ][x, y, z, t] + (ℏ^2*Derivative[0, 1, 0, 0][ρ][x, y, z, t]*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(4*m*ρ[x, y, z, t]^(3/2)) - (ℏ^2*((3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2*Derivative[0, 1, 0, 0][ρ][x, y, z, t])/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 2, 0][ρ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]^3)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 1, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) + Derivative[0, 1, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[0, 2, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 3, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) + (3*Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[1, 0, 0, 0][ρ][x, y, z, t]*Derivative[1, 1, 0, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) - (Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[2, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 1, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]))/m), -((Derivative[0, 0, 1, 0][V][x, y, z, t] + g*Derivative[0, 0, 1, 0][ρ][x, y, z, t] + (ℏ^2*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*(-1/4*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^2/ρ[x, y, z, t]^(3/2) + Derivative[0, 0, 2, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) - Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 0, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(4*m*ρ[x, y, z, t]^(3/2)) - (ℏ^2*((3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]^3)/(8*ρ[x, y, z, t]^(5/2)) - (3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 0, 2, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 0, 3, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) + (3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 1, 0, 0][ρ][x, y, z, t]^2)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[0, 1, 0, 0][ρ][x, y, z, t]*Derivative[0, 1, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[0, 2, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[0, 2, 1, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]]) + (3*Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[1, 0, 0, 0][ρ][x, y, z, t]^2)/(8*ρ[x, y, z, t]^(5/2)) - (Derivative[1, 0, 0, 0][ρ][x, y, z, t]*Derivative[1, 0, 1, 0][ρ][x, y, z, t])/(2*ρ[x, y, z, t]^(3/2)) - (Derivative[0, 0, 1, 0][ρ][x, y, z, t]*Derivative[2, 0, 0, 0][ρ][x, y, z, t])/(4*ρ[x, y, z, t]^(3/2)) + Derivative[2, 0, 1, 0][ρ][x, y, z, t]/(2*Sqrt[ρ[x, y, z, t]])))/(2*m*Sqrt[ρ[x, y, z, t]]))/m)}

(6) Vorticity transport (general):
{-Derivative[0, 0, 1, 1][vcomp][x, y, z, t] + Derivative[0, 1, 0, 1][w][x, y, z, t], Derivative[0, 0, 1, 1][u][x, y, z, t] - Derivative[1, 0, 0, 1][w][x, y, z, t], -Derivative[0, 1, 0, 1][u][x, y, z, t] + Derivative[1, 0, 0, 1][vcomp][x, y, z, t]} == {Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 0, 2, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][vcomp][x, y, z, t] - w[x, y, z, t]*Derivative[0, 1, 1, 0][w][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][w][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t], -(Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t]) - w[x, y, z, t]*Derivative[0, 0, 2, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - u[x, y, z, t]*Derivative[1, 0, 1, 0][u][x, y, z, t] + w[x, y, z, t]*Derivative[1, 0, 1, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[2, 0, 0, 0][w][x, y, z, t], Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - w[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] + u[x, y, z, t]*Derivative[1, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[2, 0, 0, 0][vcomp][x, y, z, t]}

(7) Vorticity transport (barotropic / induction-like):
{-Derivative[0, 0, 1, 1][vcomp][x, y, z, t] + Derivative[0, 1, 0, 1][w][x, y, z, t], Derivative[0, 0, 1, 1][u][x, y, z, t] - Derivative[1, 0, 0, 1][w][x, y, z, t], -Derivative[0, 1, 0, 1][u][x, y, z, t] + Derivative[1, 0, 0, 1][vcomp][x, y, z, t]} == {Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 0, 2, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][vcomp][x, y, z, t] - w[x, y, z, t]*Derivative[0, 1, 1, 0][w][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][w][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t], -(Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 0, 1, 0][w][x, y, z, t]) - w[x, y, z, t]*Derivative[0, 0, 2, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[0, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] - Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][w][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] + Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - u[x, y, z, t]*Derivative[1, 0, 1, 0][u][x, y, z, t] + w[x, y, z, t]*Derivative[1, 0, 1, 0][w][x, y, z, t] + vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][w][x, y, z, t] + u[x, y, z, t]*Derivative[2, 0, 0, 0][w][x, y, z, t], Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][vcomp][x, y, z, t] + Derivative[0, 0, 1, 0][u][x, y, z, t]*Derivative[0, 1, 0, 0][w][x, y, z, t] + w[x, y, z, t]*Derivative[0, 1, 1, 0][u][x, y, z, t] + vcomp[x, y, z, t]*Derivative[0, 2, 0, 0][u][x, y, z, t] + Derivative[0, 1, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][u][x, y, z, t] - Derivative[0, 1, 0, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[1, 0, 0, 0][u][x, y, z, t]*Derivative[1, 0, 0, 0][vcomp][x, y, z, t] - Derivative[0, 0, 1, 0][vcomp][x, y, z, t]*Derivative[1, 0, 0, 0][w][x, y, z, t] - w[x, y, z, t]*Derivative[1, 0, 1, 0][vcomp][x, y, z, t] + u[x, y, z, t]*Derivative[1, 1, 0, 0][u][x, y, z, t] - vcomp[x, y, z, t]*Derivative[1, 1, 0, 0][vcomp][x, y, z, t] - u[x, y, z, t]*Derivative[2, 0, 0, 0][vcomp][x, y, z, t]}
*)
