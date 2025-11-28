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

(* Effective potential from scalar lag: Φ_eff(r) = -μ/r - μ^2/(2 cs^2 r^2) *)
Φeff[r_] := -μ/r - μ^2/(2 cs^2 r^2);

(* Velocity squared in polar coordinates *)
Clear[v2];
v2[r_, rdot_, thetadot_] := rdot^2 + r^2 thetadot^2;

(* Full Lagrangian L(r, θ, rdot, thetadot) to 1PN order *)
Lagrangian[r_, rdot_, thetadot_] :=
  1/2 m (1 + σ[r]) v2[r, rdot, thetadot] - m Φeff[r];

Print["D1: Lagrangian defined with sigma(r) and Φ_eff(r)."];

Lexpanded = Simplify[Lagrangian[r, rdot, thetadot]];
Print["D1: Expanded Lagrangian:"];
Print["    ", Lexpanded];

(* ================================================================== *)
(* D2: Canonical Momenta and θ̇                                      *)
(* ================================================================== *)

Clear[pr, pθ];

pr  = D[Lagrangian[r, rdot, thetadot], rdot] // Simplify;
pθ  = D[Lagrangian[r, rdot, thetadot], thetadot] // Simplify;

Print["D2: Radial momentum pr = ∂L/∂rdot:"];
Print["    ", pr];

Print["D2: Angular momentum pθ = ∂L/∂thetadot:"];
Print["    ", pθ];

(* Angular momentum ℓ is conserved *)
Clear[ℓ];
(* Solve pθ == ℓ for thetadot *)
thetadotSol =
  Solve[pθ == ℓ, thetadot][[1, 1, 2]] // Simplify;

Print["D2: thetadot(ℓ,r,rdot) to 1PN order (exact expression before expansion):"];
Print["    ", thetadotSol];

thetadotSeries =
  Series[thetadotSol, {cs, Infinity, 2}] // Normal // Simplify;

Print["D2: thetadot expanded in 1/cs^2:"];
Print["    ", thetadotSeries];

(* D2: Newtonian limit check *)

thetadotNewtonian = Simplify[thetadotSeries /. β -> 0];

checkD2a = Simplify[thetadotNewtonian == ℓ/(m r^2)];

Print["D2: Newtonian limit thetadot = ℓ/(m r^2)? ", checkD2a];

(* ================================================================== *)
(* D3: Radial Equation and Orbit Equation in u(θ)                     *)
(* ================================================================== *)

Clear[eqr];

(* Euler-Lagrange for r: d/dt (∂L/∂rdot) - ∂L/∂r == 0 *)
(* D3: Correct radial Euler-Lagrange equation *)
eqrFull =
  D[pr /. {r -> r[t], rdot -> r'[t]}, t] -
  (D[Lagrangian[r, rdot, thetadot], r] /. {
     r -> r[t],
     rdot -> r'[t],
     thetadot -> θ'[t]
   }) // Simplify;

Print["D3: Full radial E-L equation (eqrFull == 0):"];
Print["    ", eqrFull];

Clear[u, rOfθ];

(* r(θ) = 1/u(θ) *)
rOfθ[θ_] := 1/u[θ];

(* Convert derivatives: dr/dt = (dr/dθ)(dθ/dt) *)
(* We'll treat r, θ as functions of t but rewrite via chain rule *)

eqrθ =
  eqrFull /. {
    r[t] -> rOfθ[θ[t]],
    r'[t] -> D[rOfθ[θ], θ] θ'[t],
    r''[t] -> D[rOfθ[θ], {θ, 2}] θ'[t]^2 + D[rOfθ[θ], θ] θ''[t]
  } // Simplify;

eqrθSub =
  eqrθ /. {
    θ'[t] -> (thetadotSeries /. r -> rOfθ[θ[t]] /. rdot -> 0),
    θ''[t] -> D[(thetadotSeries /. r -> rOfθ[θ[t]] /. rdot -> 0), t]
  } // Simplify;

(* After cleaning up eqrθSub by hand / stepwise, imagine we reach something like: *)

Clear[orbitEq];

(* ================================================================== *)
(* D3: For now, treat eqrFull == 0 as the radial EOM in time.        *)
(*      We'll later transform this into an orbit equation u(θ).       *)
(* ================================================================== *)

(* orbitEq[θ] == 0 encodes u'' + u - RHS(u) = 0 *)
(* Radial equation of motion: eqrFull == 0 *)
orbitEq = eqrFull == 0;

Print["D3: Radial EOM (time domain) = 0:"];
Print["    ", orbitEq];

(* 1/cs^2 expansion of the *left-hand side* of the EOM *)
orbitEqPN =
  Normal[Series[First@orbitEq, {cs, Infinity, 2}]] == 0 // Simplify;

Print["D3: Radial EOM truncated at O(1/cs^2):"];
Print["    ", orbitEqPN];

(* ================================================================== *)
(* D3b: Radial EOM expression (lhs == 0)                              *)
(* ================================================================== *)

Clear[eqPNExpr];

eqPNExpr =
  Simplify[Subtract @@ (List @@ orbitEqPN)];

Print["D3b: Radial EOM expression eqPNExpr (eqPNExpr == 0):"];
Print["    ", eqPNExpr];

(* ================================================================== *)
(* D4: 1PN effective potential U_eff(r) and circular orbit            *)
(* ================================================================== *)

(* 1. Generalized angular momentum: pθ = m (1 + σ[r]) r^2 thetadot = ℓ *)
(*    -> thetadot = ℓ / (m (1 + σ[r]) r^2) *)
Clear[thetabar, Ueff1PN];

thetabar[r_] := ℓ/(m (1 + σ[r]) r^2);

(* 2. Effective potential for radial motion:
   U_eff(r) = ℓ^2 / [2 m (1 + σ[r]) r^2] - m Φ_eff(r) *)
Ueff1PN[r_] :=
  ℓ^2/(2 m (1 + σ[r]) r^2) - m Φeff[r];

(* Expand U_eff to O(1/cs^2) *)
Ueff1PNExpanded =
  Simplify[
    Normal[Series[Ueff1PN[r], {cs, Infinity, 2}]]
  ];

Print["D4: U_eff,1PN (expanded to O(1/cs^2)):"];
Print["    ", Ueff1PNExpanded];

(* 3. Circular orbit condition: dU_eff/dr = 0 at r = r0 *)
Clear[circCond, ℓ2Sol];

circCond =
  Simplify[
    D[Ueff1PNExpanded, r] /. r -> r0
  ];

Print["D4: Circular condition dU_eff/dr|_(r0)=0:"];
Print["    ", circCond];

(* Solve circular condition for ℓ^2 using dummy variable y *)
Clear[y, ℓ2Sol];

ℓ2Sol =
  First[
    Solve[(circCond == 0) /. ℓ^2 -> y, y]
  ] /. y -> ℓ^2 // Simplify;

Print["D4: ℓ^2 at circular orbit (to O(1/cs^2)):"];
Print["    ", ℓ2Sol];

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

Print["D4b: κ^2 (expanded to O(1/cs^2)):"];
Print["    ", κ2Series];

(* Newtonian sanity: κ^2 and Ω^2 should both → μ/r0^3 *)
κ2Newt = Simplify[κ2Series /. {β -> 0, cs -> Infinity}];
Ω2Newt = Simplify[Ω2Series /. {β -> 0, cs -> Infinity}];

Print["D4b: κ^2 Newtonian limit: ", κ2Newt];
Print["D4c: Ω^2 Newtonian limit: ", Ω2Newt];
Print["D4: κ^2 == Ω^2 in Newtonian limit? ",
      Simplify[κ2Newt == Ω2Newt]];

(* ================================================================== *)
(* D4c: Orbital frequency Ω and precession Δφ                         *)
(* ================================================================== *)

Clear[Ω2Expr, Ω2Series, ratioκ2Ω2, ratioκΩ, ΔφExpr, ΔφFactor, C0, Cβ];

(* Ω = θ̇ at circular orbit: thetabar[r0] with ℓ^2 from ℓ2Sol *)
Ω2Expr =
  Simplify[
    (thetabar[r0]^2) /. ℓ2Sol
  ];

Ω2Series =
  Simplify[
    Normal[Series[Ω2Expr, {cs, Infinity, 2}]]
  ];

Print["D4c: Ω^2 (expanded to O(1/cs^2)):"];
Print["    ", Ω2Series];

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

checkScalarOneSixth = Simplify[C0 == 1];
Print["D4c: Does scalar-only piece give the 1/6 factor (C0=1)? ",
      checkScalarOneSixth];

ΔφSimplified = Simplify[ΔφExpr /. r0 -> a];

(* ================================================================== *)
(* D5: GR matching: solve for β (small-e limit, r0 ≈ a)               *)
(* ================================================================== *)

Clear[βSol];

βSol =
  Simplify[
    Solve[C0 + Cβ β == 6, β][[1, 1, 2]]
  ];

Print["D5: β needed to match GR factor 6 (small-e limit): β = ", βSol];

(* ================================================================== *)
(* D6: Numerical sanity check                                        *)
(* ================================================================== *)

paramsNumeric = {
  μ -> 1.,
  cs -> 20.,   (* large but finite *)
  m -> 1.,
  a -> 10.,
  e -> 0.2,
  β -> βSol /. μ -> 1 /. cs -> 20
};

(* Initial conditions for near-Keplerian orbit: set initial r, rdot, θ, thetadot *)

(* Build numeric EOMs from Lagrangian using Euler-Lagrange or from Hamiltonian *)

(* Use NDSolve to integrate over, say, 20 orbital periods *)

(* Extract times where r[t] is at perihelion (dr/dt = 0, d^2r/dt^2 > 0) *)

(* Compute measured precession per orbit and compare to analytic ΔφSimplified *)

ΔφAnalyticNumeric =
  (ΔφSimplified /. paramsNumeric) // N;

Print["D6: Analytic Δφ with numeric params: ", ΔφAnalyticNumeric, " radians"];

(* After implementing the NDSolve bits, set checkD6 = True/False based on agreement *)

(* ================================================================== *)
(* D7: CHECK SUMMARY                                                 *)
(* ================================================================== *)

allPass =
  checkD2a && checkScalarOneSixth (* && checkD6 once implemented *);

Print["============================================================="];
Print["CHECK SUMMARY"];
Print["============================================================="];
Print["  D2: Newtonian limit of thetadot:                     ", checkD2a];
Print["  D4: Scalar-only limit reproduces 1/6 piece:          ", checkScalarOneSixth];
(* Print["  D6: Numeric orbit vs analytic Δφ agrees:           ", checkD6]; *)
Print["============================================================="];
Print["All passes? ", allPass];
Print["============================================================="];
