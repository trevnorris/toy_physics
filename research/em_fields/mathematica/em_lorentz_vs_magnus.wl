(* em_lorentz_vs_magnus.wl *)

(* ---------------------------------------------------------------------- *)
(* 0. Parameters and Setup *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"];

$Assumptions = {
  rho0 > 0, Gamma > 0, a > 0, L > 0, kappaQ > 0,
  B0 \[Element] Reals, v0 \[Element] Reals,
  ux \[Element] Reals, uy \[Element] Reals, uz \[Element] Reals
};

uVec = {ux, uy, uz};
kVec = {0, 0, 1};
vInf = {v0, 0, 0};

(* ---------------------------------------------------------------------- *)
(* 1. Magnus force on a straight vortex line (per unit length) *)
(* ---------------------------------------------------------------------- *)

FM = rho0 * Gamma * Cross[kVec, (uVec - vInf)];

Print["Magnus force per unit length, f_M = rho0 * Gamma * k x (u - v_inf):"];
Print[ToString[Simplify[FM], InputForm]];
Print[""];

FMuPart = rho0 * Gamma * Cross[kVec, uVec];
FMvinfPart = -rho0 * Gamma * Cross[kVec, vInf];

Print["u-dependent (magnetic-like) part of Magnus force:"];
Print[ToString[Simplify[FMuPart], InputForm]];
Print[""];

Print["Background-flow contribution to Magnus force:"];
Print[ToString[Simplify[FMvinfPart], InputForm]];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 2. Magnetic Lorentz force per unit length *)
(* ---------------------------------------------------------------------- *)

qDef = kappaQ * rho0 * Pi * a^2 * Gamma;
BVec = B0 * kVec;
FLmagPerLength = (qDef/L) * Cross[uVec, BVec];

Print["Toy-model effective charge: q = kappaQ * rho0 * pi * a^2 * Gamma"];
Print["q = ", ToString[qDef, InputForm]];
Print[""];

Print["Magnetic Lorentz force per unit length, f_L,mag = (q/L) * u x B:"];
Print[ToString[Simplify[FLmagPerLength], InputForm]];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 3. Compare with the u-dependent Magnus piece *)
(* ---------------------------------------------------------------------- *)

FMuSimplified = Simplify[FMuPart];
FLuSimplified = Simplify[FLmagPerLength];
solutionB0 = Solve[FMuSimplified == FLuSimplified, B0];

Print["u-dependent Magnus term:"];
Print[ToString[FMuSimplified, InputForm]];
Print[""];

Print["u-dependent Lorentz term (with geometric q):"];
Print[ToString[FLuSimplified, InputForm]];
Print[""];

Print["
We can read off the proportionality by inspection.
Mathematica solution for B0 equating the two forces:
"];
Print[ToString[solutionB0, InputForm]];

Print["
Thus, with the identification
    q = kappaQ * rho0 * pi * a^2 * Gamma,
    B = -(L / (kappaQ * pi * a^2)) k,
the u-dependent Magnus force matches the magnetic part of the Lorentz force
per unit length in the straight-vortex geometry used in the paper.

The remaining background-flow term can be grouped with
pressure-gradient forces and interpreted as part of the effective
electric field contribution qE.
"];

Print["em_lorentz_vs_magnus.wl: algebraic comparison complete."];

(*
Output:
Magnus force per unit length, f_M = rho0 * Gamma * k x (u - v_inf):
{-(Gamma*rho0*uy), Gamma*rho0*(ux - v0), 0}

u-dependent (magnetic-like) part of Magnus force:
{-(Gamma*rho0*uy), Gamma*rho0*ux, 0}

Background-flow contribution to Magnus force:
{0, -(Gamma*rho0*v0), 0}

Toy-model effective charge: q = kappaQ * rho0 * pi * a^2 * Gamma
q = a^2*Gamma*kappaQ*Pi*rho0

Magnetic Lorentz force per unit length, f_L,mag = (q/L) * u x B:
{(a^2*B0*Gamma*kappaQ*Pi*rho0*uy)/L, -((a^2*B0*Gamma*kappaQ*Pi*rho0*ux)/L), 0}

u-dependent Magnus term:
{-(Gamma*rho0*uy), Gamma*rho0*ux, 0}

u-dependent Lorentz term (with geometric q):
{(a^2*B0*Gamma*kappaQ*Pi*rho0*uy)/L, -((a^2*B0*Gamma*kappaQ*Pi*rho0*ux)/L), 0}


We can read off the proportionality by inspection.
Mathematica solution for B0 equating the two forces:

{{B0 -> -(L/(a^2*kappaQ*Pi))}}

Thus, with the identification
    q = kappaQ * rho0 * pi * a^2 * Gamma,
    B = -(L / (kappaQ * pi * a^2)) k,
the u-dependent Magnus force matches the magnetic part of the Lorentz force
per unit length in the straight-vortex geometry used in the paper.

The remaining background-flow term can be grouped with
pressure-gradient forces and interpreted as part of the effective
electric field contribution qE.

em_lorentz_vs_magnus.wl: algebraic comparison complete.
*)
