(* ===================================================================== *)
(* 1PN PERIHELION PRECESSION FROM POSITION-DEPENDENT INERTIA (GENERAL e)  *)
(* --------------------------------------------------------------------- *)
(* Fixes vs v1:                                                          *)
(*  - Avoids u = u[θ] self-recursion.                                     *)
(*  - Avoids FullSimplify recursion.                                      *)
(*  - Series expansion uses a true symbol (η), not μ/(cS^2 p).             *)
(* ===================================================================== *)
(* Model:                                                                *)
(*   L = (1/2) m (1 + σ(r)) ( rdot^2 + r^2 thetadot^2 ) + m μ/r           *)
(*   σ(r) = β μ/(cS^2 r)  (small, 1PN-order)                              *)
(*                                                                        *)
(* Output:                                                                *)
(*   Δϖ = (2β) π μ/(cS^2 a (1 - e^2))   (to first order)                  *)
(* ===================================================================== *)

ClearAll["Global`*"];

(* Small bookkeeping parameter for "keep only O(σ)" *)
ε = Symbol["\[Epsilon]"];

(* Parameters *)
(* IMPORTANT:
   Do NOT set p == a(1-e^2) inside $Assumptions, because FullSimplify may
   rewrite a(1-e^2) back into p. We'll substitute p->a(1-e^2) explicitly
   only when printing the (a,e) form. *)
$Assumptions = $Assumptions && a > 0 && 0 <= e < 1 && p > 0;
assumeAE = (a > 0 && 0 <= e < 1 && μ > 0 && cS > 0 && β ∈ Reals);
pRule = p -> a (1 - e^2);

(* Use u(θ) = 1/r(θ) *)
(* σ expressed as a function of u *)
sigmaU[u_] := ε * β * μ * u / cS^2;
sigmaPrimeU[u_] := D[sigmaU[s], s] /. s -> u;  (* dσ/du *)

(* Constants in the first integral (exact identities) *)
k  = m^2 μ/J^2;   (* Newtonian: k = 1/p *)
C0 = 2 m E/J^2;   (* Newtonian: C0 = -1/(a p) *)

(* --------------------------------------------------------------------- *)
(* Step 1: Binet-type equation from the first integral                    *)
(* --------------------------------------------------------------------- *)
(* Starting point (true for this L):
     u'^2 + u^2 = (1 + σ(u)) ( C0 + 2 k u ).
   Differentiate and divide by 2 u' to get:
     u'' + u = k (1 + σ(u)) + 1/2 (C0 + 2 k u) σ'(u),
   where σ'(u) = dσ/du.
*)

uθ = u[θ];

binetEqExact =
  D[uθ, {θ, 2}] + uθ ==
    k (1 + sigmaU[uθ]) + (C0 + 2 k uθ)/2 * sigmaPrimeU[uθ];

binetEq1PN =
  Normal@Series[binetEqExact, {ε, 0, 1}] /. ε -> 1;

Print["--- Binet equation to O(σ) (general form) ---"];
Print[binetEq1PN];

(* --------------------------------------------------------------------- *)
(* Step 2: Express with Newtonian orbital elements (a,e,p)                *)
(* --------------------------------------------------------------------- *)
ClearAll[a, e, p];
$Assumptions = $Assumptions && a > 0 && 0 <= e < 1 && p == a (1 - e^2);

rulesKepler = {k -> 1/p, C0 -> -1/(a p)};

binetEqAE =
  FullSimplify[binetEq1PN /. rulesKepler, Assumptions -> $Assumptions];

Print["\n--- Binet equation in terms of (a,e,p) ---"];
Print[binetEqAE];

(* --------------------------------------------------------------------- *)
(* Step 3: Extract the frequency shift (the only part that sets Δϖ)       *)
(* --------------------------------------------------------------------- *)
(* The equation has the structure:
     u'' + u == 1/p + const + (2 β μ/(cS^2 p)) u
   Move the u-term to the left:
     u'' + (1 - 2 β μ/(cS^2 p)) u == 1/p + const
   So ω^2 = 1 - 2 β μ/(cS^2 p).
*)

ω2 = 1 - 2 β μ/(cS^2 p);
ω  = Sqrt[ω2];

Print["\n--- Frequency (squared) ω^2 ---"];
Print[ω2];

(* Perihelion advance per orbit: Δϖ = 2π(1/ω - 1) *)
deltaVarpiExact = 2 Pi (1/ω - 1);

(* Make a clean 1st-order series in a SYMBOL η, then substitute η = β μ/(cS^2 p) *)
η = Symbol["\[Eta]"];
deltaVarpiSeries =
  Normal@Series[2 Pi (1/Sqrt[1 - 2 η] - 1), {η, 0, 1}] /. η -> (β μ/(cS^2 p));

deltaVarpiAE = Simplify[deltaVarpiSeries /. (p -> a (1 - e^2)),
  Assumptions -> (a > 0 && 0 <= e < 1 && μ > 0 && cS > 0 && β ∈ Reals)];

Print["\n--- Δϖ per orbit ---"];
Print["Exact (closed form) in p: ", deltaVarpiExact];
Print["Series to first order in μ/(cS^2 p): ", deltaVarpiSeries];
Print["Series to first order in a,e: ", deltaVarpiAE];

Print["\n--- DONE ---"];

(*"
Output:

--- Binet equation to O(σ) (general form) ---
u[θ] + Derivative[2][u][θ] == (m^2*μ)/J^2 + (m*β*μ*(E + 2*m*μ*u[θ]))/(cS^2*J^2)

--- Binet equation in terms of (a,e,p) ---
u[θ] + Derivative[2][u][θ] == p^(-1) + (m*β*μ*(E + 2*m*μ*u[θ]))/(cS^2*J^2)

--- Frequency (squared) ω^2 ---
1 - (2*β*μ)/(cS^2*p)

--- Δϖ per orbit ---
Exact (closed form) in p: 2*Pi*(-1 + 1/Sqrt[1 - (2*β*μ)/(cS^2*p)])
Series to first order in μ/(cS^2 p): (2*Pi*β*μ)/(cS^2*p)
Series to first order in a,e: (-2*Pi*β*μ)/(a*cS^2*(-1 + e^2))
"*)
