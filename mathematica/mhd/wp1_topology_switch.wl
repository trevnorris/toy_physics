(* ::Package:: *)

(* ===================================================================== *)
(* WP1_TOPOLOGY_SWITCH_v2.wl                                             *)
(*                                                                      *)
(* Fast / robust version: avoids heavy FullSimplify explosions.          *)
(*                                                                      *)
(* What it does                                                         *)
(*  1) Verifies the vector identity relating (v·∇)v to v×ω and ∇(v^2/2).  *)
(*  2) Confirms the Rosetta-stone mapping A=λ v, B=∇×A=λ ω,              *)
(*     φ=λ(h+v^2/2) gives ideal Ohm: E + v×B = 0 under inviscid Euler +  *)
(*     incompressibility (Div[v]=0).                                     *)
(*  3) Sets up the "switch slot" term -∇×(η ∇×B) without expanding it.    *)
(*                                                                      *)
(* Notes                                                                *)
(*  - Mathematica's symbol E is protected (Euler's number).              *)
(*    This script uses eField instead of E.                              *)
(*  - The variable-η residual expansion is intentionally NOT performed   *)
(*    by default because it can be extremely slow symbolically.          *)
(* ===================================================================== *)

ClearAll["Global`*"];

(* ---------- USER TOGGLES ---------- *)
ComputeVariableEtaResidual = False;  (* set True only if you really want the full expanded residual *)

(* ---------- ASSUMPTIONS ---------- *)
$Assumptions =
  Element[{x, y, z, t, lam, nu, etaBg, dEta, t0, tau, sigma, x0, y0, z0}, Reals] &&
  lam != 0 && nu >= 0 && etaBg >= 0 && dEta >= 0 && tau > 0 && sigma > 0;

coords = {x, y, z};

grad[f_] := Grad[f, coords];
div[F_]  := Div[F, coords];
curl[F_] := Curl[F, coords];
lap[f_]  := Laplacian[f, coords];
vLap[F_] := lap /@ F;  (* component-wise Laplacian *)

AdvectVector[U_, F_] := Table[
  Sum[U[[i]]*D[F[[j]], coords[[i]]], {i, 1, 3}],
  {j, 1, 3}
];

(* Fast simplifier: try FullSimplify but fall back quickly *)
simp[expr_, seconds_: 1.0] := TimeConstrained[
  FullSimplify[expr, Assumptions -> $Assumptions],
  seconds,
  Simplify[expr, Assumptions -> $Assumptions]
];

(* ---------- SYMBOLIC FIELDS ---------- *)
v = {vx[x, y, z, t], vy[x, y, z, t], vz[x, y, z, t]};
h = hField[x, y, z, t];

omega = curl[v];

(* ---------- 1) VECTOR IDENTITY CHECK ---------- *)
(* Identity:
   (v·∇)v = ∇(v^2/2) - v×ω - v(∇·v)
   So the residual should be exactly v(∇·v). *)
identityResidual =
  simp[
    AdvectVector[v, v] - (grad[1/2 (v.v)] - Cross[v, omega] - v*div[v]),
    0.5
  ];

Print["\n--- CHECK 1: vector identity residual ---"];
Print["Expected: {vx,vy,vz}*(Div[v]). This is NOT an error (it vanishes if Div[v]=0)."];
Print[identityResidual];

identityResidualIncomp =
  simp[identityResidual /. div[v] -> 0, 0.5];

Print["\n--- CHECK 1b: identity residual with Div[v]=0 (should be {0,0,0}) ---"];
Print[identityResidualIncomp];

(* ---------- 2) ROSETTA STONE DICTIONARY ---------- *)
A    = lam v;
B    = curl[A];                         (* = lam omega *)
phiEM = lam (h + 1/2 (v.v));
eField = -grad[phiEM] - D[A, t];

Print["\n--- Sanity: B - lam*omega (should be {0,0,0}) ---"];
Print[simp[B - lam omega, 0.5]];

(* ---------- 3) IDEAL OHM'S LAW CHECK ---------- *)
(* Inviscid Euler: ∂t v + (v·∇)v = -∇h *)
dvdt = D[v, t];
eulerRule = Thread[dvdt -> -AdvectVector[v, v] - grad[h]];

advectToIdentityRule = Thread[
  AdvectVector[v, v] -> (grad[1/2 (v.v)] - Cross[v, omega] - v*div[v])
];

ohmResidual =
  simp[
    (eField + Cross[v, B]) /. eulerRule /. advectToIdentityRule /. div[v] -> 0,
    1.0
  ];

Print["\n--- CHECK 2: Ideal Ohm residual (eField + v×B) under Euler + Div[v]=0 ---"];
Print[ohmResidual];
If[ohmResidual === {0, 0, 0},
  Print["SUCCESS: Ideal Ohm's law holds (frozen-in condition)."],
  Print["NOTE: Nonzero residual. This usually means a definition/assumption mismatch."]
];

(* ---------- 4) DIFFUSION COMMUTATOR CHECK ---------- *)
diffCommResidual =
  simp[curl[vLap[v]] - vLap[curl[v]], 0.5];

Print["\n--- CHECK 3: curl(Laplacian(v)) - Laplacian(curl(v)) (should be {0,0,0}) ---"];
Print[diffCommResidual];

(* ---------- 5) "SWITCH SLOT": resistive/viscous term setup ---------- *)
(* IMPORTANT:
   The induction/vorticity *equations of motion* are not imposed here, so
   "induction residual" is generally nonzero. What we check is consistency:
   B = lam*omega implies (induction residual) = lam*(vorticity residual). *)

vorticityResidual =
  D[omega, t] - curl[Cross[v, omega]] - nu vLap[omega];

inductionResidualConstant =
  D[B, t] - curl[Cross[v, B]] - nu vLap[B];

Print["\n--- CHECK 4: Consistency (inductionResidualConstant - lam*vorticityResidual) ---"];
Print[simp[inductionResidualConstant - lam*vorticityResidual, 1.0]];

Print["\nInterpretation:"];
Print["- The huge expression you saw for 'induction residual' is expected because we did NOT assume"];
Print["  any specific PDE for v beyond the identity checks. It's 'what would need to vanish' if v obeyed"];
Print["  the viscous equations of motion."];
Print["- The important check is that the induction residual matches λ times the vorticity residual (it does)."];

(* ---------- 6) PULSED, LOCALIZED ETA(x,t) (without expanding) ---------- *)
Wt[tt_] := UnitStep[tt - t0] - UnitStep[tt - (t0 + tau)];
Wx[xx_, yy_, zz_] := Exp[-((xx - x0)^2 + (yy - y0)^2 + (zz - z0)^2)/(2 sigma^2)];
eta[xx_, yy_, zz_, tt_] := etaBg + dEta*Wt[tt]*Wx[xx, yy, zz];

ResistiveTerm[Bvec_] := -curl[ eta[x, y, z, t]*curl[Bvec] ];

Print["\n--- Switch slot term (structural form, not expanded) ---"];
Print[HoldForm[D[B, t] == curl[Cross[v, B]] + ResistiveTerm[B]]];

If[ComputeVariableEtaResidual === True,
  Print["\n--- WARNING: Expanding variable-eta induction residual can be very slow ---"];
  inductionResidualVariable =
    simp[D[B, t] - curl[Cross[v, B]] - ResistiveTerm[B], 2.0];
  Print["Variable-eta induction residual (may be large):"];
  Print[Short[inductionResidualVariable, 4]];
];

(* ---------- 7) WP1 SCALARS (formulas you print per run in Python) ---------- *)
Rm[U0_, ell_, eta0_] := (U0*ell)/eta0;
ld[eta0_, tauPulse_] := Sqrt[eta0*tauPulse];
chi[eta0_, tauPulse_, delta0_] := ld[eta0, tauPulse]/delta0;

Print["\n--- WP1 scalars (symbolic formulas) ---"];
Print[HoldForm[Rm[U, δ, η] -> (U*δ)/η]];
Print[HoldForm[ld[η, τ] -> Sqrt[η*τ]]];
Print[HoldForm[χ[η, τ, δ] -> Sqrt[η*τ]/δ]];

(* End *)

(*"
Output:


--- CHECK 1: vector identity residual ---
Expected: {vx,vy,vz}*(Div[v]). This is NOT an error (it vanishes if Div[v]=0).
{vx[x, y, z, t]*(Derivative[0, 0, 1, 0][vz][x, y, z, t] + Derivative[0, 1, 0, 0][vy][x, y, z, t] + Derivative[1, 0, 0, 0][vx][x, y, z, t]), vy[x, y, z, t]*(Derivative[0, 0, 1, 0][vz][x, y, z, t] + Derivative[0, 1, 0, 0][vy][x, y, z, t] + Derivative[1, 0, 0, 0][vx][x, y, z, t]), vz[x, y, z, t]*(Derivative[0, 0, 1, 0][vz][x, y, z, t] + Derivative[0, 1, 0, 0][vy][x, y, z, t] + Derivative[1, 0, 0, 0][vx][x, y, z, t])}

--- CHECK 1b: identity residual with Div[v]=0 (should be {0,0,0}) ---
{0, 0, 0}

--- Sanity: B - lam*omega (should be {0,0,0}) ---
{0, 0, 0}

--- CHECK 2: Ideal Ohm residual (eField + v×B) under Euler + Div[v]=0 ---
{0, 0, 0}
SUCCESS: Ideal Ohm's law holds (frozen-in condition).

--- CHECK 3: curl(Laplacian(v)) - Laplacian(curl(v)) (should be {0,0,0}) ---
{0, 0, 0}

--- CHECK 4: Consistency (inductionResidualConstant - lam*vorticityResidual) ---
{0, 0, 0}

Interpretation:
- The huge expression you saw for 'induction residual' is expected because we did NOT assume
  any specific PDE for v beyond the identity checks. It's 'what would need to vanish' if v obeyed
  the viscous equations of motion.
- The important check is that the induction residual matches λ times the vorticity residual (it does).

--- Switch slot term (structural form, not expanded) ---
HoldForm[D[B, t] == curl[Cross[v, B]] + ResistiveTerm[B]]

--- WP1 scalars (symbolic formulas) ---
HoldForm[Rm[U, δ, η] -> (U*δ)/η]
HoldForm[ld[η, τ] -> Sqrt[η*τ]]
HoldForm[χ[η, τ, δ] -> Sqrt[η*τ]/δ]
"*)
