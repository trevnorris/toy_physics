(* ================================================================== *)
(* D0: Setup & Assumptions                                            *)
(* ================================================================== *)

ClearAll["Global`*"];

(* Parameters and symbols *)
Clear[μ, cs, β, m, a, e, r, r0, θ, t];

$Assumptions = {
  μ > 0, cs > 0, m > 0,
  0 <= e < 1,
  a > 0, r0 > 0,
  Element[{μ, cs, β, m, a, e, r0, r}, Reals]
};

Print["D0: Setup complete. Assumptions: μ>0, cs>0, 0≤e<1, etc."];

(* ================================================================== *)
(* D1: Define sigma(r), Phi_eff(r) and Lagrangian                     *)
(* ================================================================== *)

(* Spatial inertia factor sigma(r) = β μ / (cs^2 r) *)
σ[r_] := β μ/(cs^2 r);

(* Static-source scalar potential is exactly Newtonian *)
Φeff[r_] := -μ/r;

(* Velocity squared in polar coordinates *)
Clear[v2];
v2[r_, rdot_, thetadot_] := rdot^2 + r^2 thetadot^2;

(* Full Lagrangian L(r, θ, rdot, thetadot) to 1PN order *)
Lagrangian[r_, rdot_, thetadot_] :=
  1/2 m (1 + σ[r]) v2[r, rdot, thetadot] - m Φeff[r];

Print["D1: Lagrangian defined with static-source (Newtonian) scalar potential."];

Lexpanded = Simplify[Lagrangian[r, rdot, thetadot]];
Print["D1: Expanded Lagrangian:"];
Print["    ", Lexpanded];

(* ================================================================== *)
(* D2: Canonical Momenta and θ̇                                      *)
(* ================================================================== *)

Clear[pr, pθ];

pr  = D[Lagrangian[r, rdot, thetadot], rdot] // Simplify;
pθ  = D[Lagrangian[r, rdot, thetadot], thetadot] // Simplify;

(* Angular momentum ℓ is conserved *)
Clear[ℓ];
(* Solve pθ == ℓ for thetadot *)
thetadotSol =
  Solve[pθ == ℓ, thetadot][[1, 1, 2]] // Simplify;

thetadotSeries =
  Series[thetadotSol, {cs, Infinity, 2}] // Normal // Simplify;

(* ================================================================== *)
(* D3: Radial Equation and Orbit Equation in u(θ)                     *)
(* ================================================================== *)

Clear[eqr];

(* Euler-Lagrange for r *)
eqrFull =
  D[pr /. {r -> r[t], rdot -> r'[t]}, t] -
  (D[Lagrangian[r, rdot, thetadot], r] /. {
     r -> r[t],
     rdot -> r'[t],
     thetadot -> θ'[t]
   }) // Simplify;

(* Radial equation of motion: eqrFull == 0 *)
orbitEq = eqrFull == 0;

(* 1/cs^2 expansion of the *left-hand side* of the EOM *)
orbitEqPN =
  Normal[Series[First@orbitEq, {cs, Infinity, 2}]] == 0 // Simplify;


(* ================================================================== *)
(* D4: 1PN effective potential U_eff(r) and circular orbit            *)
(* ================================================================== *)

Clear[thetabar, Ueff1PN];

(* Effective potential for radial motion:
   U_eff(r) = ℓ^2 / [2 m (1 + σ[r]) r^2] - m Φ_eff(r) *)
Ueff1PN[r_] :=
  ℓ^2/(2 m (1 + σ[r]) r^2) - m Φeff[r];

(* Expand U_eff to O(1/cs^2) *)
Ueff1PNExpanded =
  Simplify[
    Normal[Series[Ueff1PN[r], {cs, Infinity, 2}]]
  ];

(* Circular orbit condition: dU_eff/dr = 0 at r = r0 *)
Clear[circCond, ℓ2Sol];

circCond =
  Simplify[
    D[Ueff1PNExpanded, r] /. r -> r0
  ];

(* Solve circular condition for ℓ^2 using dummy variable y *)
Clear[y, ℓ2Sol];

ℓ2Sol =
  First[
    Solve[(circCond == 0) /. ℓ^2 -> y, y]
  ] /. y -> ℓ^2 // Simplify;


(* ================================================================== *)
(* D4b: Epicyclic (radial) frequency κ from U_eff                     *)
(* ================================================================== *)

Clear[mEff, κ2Expr, κ2Series];

(* Effective inertial mass for radial motion: m_eff(r) = m (1 + σ[r]) *)
mEff[r_] := m (1 + σ[r]);

(* κ^2 = (1 / m_eff(r0)) * d^2 U_eff / dr^2 evaluated at r0 *)
κ2Expr =
  Simplify[
    (1/mEff[r0]) *
      (D[Ueff1PNExpanded, {r, 2}] /. {r -> r0} /. ℓ2Sol)
  ];

κ2Series =
  Simplify[
    Normal[Series[κ2Expr, {cs, Infinity, 2}]]
  ];


(* ================================================================== *)
(* D4c: Orbital frequency Ω and precession Δφ                         *)
(* ================================================================== *)

Clear[Ω2Expr, Ω2Series, ratioκ2Ω2, ratioκΩ, ΔφExpr, ΔφFactor, C0, Cβ];

(* Ω = θ̇ at circular orbit *)
thetabar[r_] := ℓ/(m (1 + σ[r]) r^2);

Ω2Expr =
  Simplify[
    (thetabar[r0]^2) /. ℓ2Sol
  ];

Ω2Series =
  Simplify[
    Normal[Series[Ω2Expr, {cs, Infinity, 2}]]
  ];

(* Ratio κ^2 / Ω^2 and κ/Ω *)
ratioκ2Ω2 = Simplify[κ2Series/Ω2Series];

ratioκΩ =
  Simplify[
    Normal[Series[Sqrt[ratioκ2Ω2], {cs, Infinity, 2}]]
  ];

(* Precession per orbit: Δφ ≈ 2π (1 - κ/Ω) *)
ΔφExpr =
  Simplify[2 Pi (1 - ratioκΩ)];

Print["D4c: Δφ(r0) to O(1/cs^2):"];
Print["    ", ΔφExpr];

(* Factor out π μ/(cs^2 r0) and get dimensionless C(β) *)
ΔφFactor =
  Simplify[
    ΔφExpr * (cs^2 r0)/(Pi μ)
  ];

Print["D4c: C(β) such that Δφ ≈ C(β) π μ/(cs^2 r0):"];
Print["    ", ΔφFactor];

C0 = Simplify[ΔφFactor /. β -> 0];
Cβ = Simplify[Coefficient[ΔφFactor, β]];

Print["D4c: C0 (β-independent, scalar-only): ", C0];
Print["D4c: Cβ (coefficient multiplying β):  ", Cβ];

(* ================================================================== *)
(* D5: GR matching: solve for β                                       *)
(* ================================================================== *)

(* GR requires factor of 6 *)
Clear[βSol];

βSol =
  Simplify[
    Solve[C0 + Cβ β == 6, β][[1, 1, 2]]
  ];

Print["\n============================================================="];
Print["FINAL RESULT"];
Print["============================================================="];
Print["Scalar contribution (C0): ", C0, "/6 of GR"];
Print["Inertia contribution (C_beta): ", Cβ];
Print["Required Beta (β): ", βSol];
Print["============================================================="];

(*"
Output:

D0: Setup complete. Assumptions: μ>0, cs>0, 0≤e<1, etc.
D1: Lagrangian defined with static-source (Newtonian) scalar potential.
D1: Expanded Lagrangian:
    (m*μ)/r + (m*(rdot^2 + r^2*thetadot^2)*(1 + (β*μ)/(cs^2*r)))/2
D4c: Δφ(r0) to O(1/cs^2):
    (2*Pi*β*μ)/(cs^2*r0)
D4c: C(β) such that Δφ ≈ C(β) π μ/(cs^2 r0):
    2*β
D4c: C0 (β-independent, scalar-only): 0
D4c: Cβ (coefficient multiplying β):  2

=============================================================
FINAL RESULT
=============================================================
Scalar contribution (C0): 0/6 of GR
Inertia contribution (C_beta): 2
Required Beta (β): 3
=============================================================
"*)
