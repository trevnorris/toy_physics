(* em_lorentz_vs_magnus.wl *)

(* ---------------------------------------------------------------------- *)
(* 0. Parameters and Setup *)
(* ---------------------------------------------------------------------- *)

$Assumptions = {
  rho0 > 0, a > 0, Gamma > 0,
  q \[Element] Reals, B0 \[Element] Reals, v0 \[Element] Reals,
  ux \[Element] Reals, uy \[Element] Reals, uz \[Element] Reals
};

(* Vectors defined as Lists {x, y, z} *)
iVec = {1, 0, 0};
jVec = {0, 1, 0};
kVec = {0, 0, 1};

uVec = {ux, uy, uz};
vInf = v0 * iVec; (* Background flow in +x *)

(* ---------------------------------------------------------------------- *)
(* 1. Magnus force on a straight vortex line (per unit length) *)
(* ---------------------------------------------------------------------- *)

(* Classical Magnus force: F_M = rho0 * Gamma * k x (u - v_inf) *)
FM = rho0 * Gamma * Cross[kVec, uVec - vInf];

Print["Magnus force per unit length, F_M = rho0 * Gamma * k x (u - v_inf):"];
Print[ToString[Simplify[FM], InputForm]];
Print[""];

(* Split into u-dependent (magnetic-like) and v_inf-dependent parts *)
FMuPart = rho0 * Gamma * Cross[kVec, uVec];
FMvinfPart = -rho0 * Gamma * Cross[kVec, vInf];

Print["u-dependent (magnetic-like) part of Magnus force:"];
Print[ToString[Simplify[FMuPart], InputForm]];
Print[""];

Print["Background-flow (v_inf) contribution to Magnus force:"];
Print[ToString[Simplify[FMvinfPart], InputForm]];
Print[""];


(* ---------------------------------------------------------------------- *)
(* 2. Lorentz force with q = rho0 * pi * a^2 * Gamma *)
(* ---------------------------------------------------------------------- *)

qDef = rho0 * Pi * a^2 * Gamma;
Print["Toy-model effective charge: q_defect = rho0 * pi * a^2 * Gamma"];
Print["q = ", ToString[qDef, InputForm]];
Print[""];

(* Effective magnetic field B along z-axis *)
BVec = B0 * kVec;

(* Magnetic part of Lorentz force: F_L,mag = q * u x B *)
FLmag = q * Cross[uVec, BVec];

Print["Magnetic part of Lorentz force, F_L,mag = q * u x B:"];
Print[ToString[Simplify[FLmag], InputForm]];
Print[""];

(* Substitute q_def into Lorentz force *)
FLmagQdef = FLmag /. q -> qDef;

Print["F_L,mag with q = rho0 * pi * a^2 * Gamma:"];
Print[ToString[Simplify[FLmagQdef], InputForm]];
Print[""];


(* ---------------------------------------------------------------------- *)
(* 3. Compare with Magnus u-part *)
(* ---------------------------------------------------------------------- *)

FMuSimplified = Simplify[FMuPart];
FLuSimplified = Simplify[FLmagQdef];

Print["u-dependent Magnus term:"];
Print[ToString[FMuSimplified, InputForm]];
Print[""];

Print["u-dependent Lorentz term (with q_def):"];
Print[ToString[FLuSimplified, InputForm]];
Print[""];

(* Explicit Check: Solve for B0 *)
(* We equate the coefficients of the Cross product or the vectors themselves *)
solutionB0 = Solve[FMuSimplified == FLuSimplified, B0];

Print["
We can read off the proportionality by inspection.
Mathematica solution for B0 equating the two forces:
"];
Print[ToString[solutionB0, InputForm]];

Print["
Thus, with the identification
    q = rho0 * pi * a^2 * Gamma,
    B = -(1 / (pi * a^2)) k,
the u-dependent Magnus force matches the magnetic part of the Lorentz force.

The remaining background-flow term F_M_vinf can be grouped with
pressure-gradient forces and interpreted as part of an effective
electric field E in the Lorentz force qE.
"];

Print["em_lorentz_vs_magnus.wl: algebraic comparison complete."];
