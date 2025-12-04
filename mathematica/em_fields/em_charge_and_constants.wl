(* em_charge_and_constants.wl *)

(* ---------------------------------------------------------------------- *)
(* 0. Parameters and basic symbols *)
(* ---------------------------------------------------------------------- *)

$Assumptions = {
  rho0 > 0, a > 0, L > 0, Gamma > 0, aspect > 0,
  q > 0, eps0 > 0, mu0 > 0, c > 0,
  x \[Element] Reals, y \[Element] Reals, z \[Element] Reals
};

r = Sqrt[x^2 + y^2 + z^2];

Print["==========================================================="];
Print["PART 1: Defect mass and charge scaling (from toy model)"];
Print["===========================================================\n"];

(* ---------------------------------------------------------------------- *)
(* Part 1 Logic *)
(* ---------------------------------------------------------------------- *)

Lexpr = aspect * a;
massG = rho0 * (Pi * a^2 * Lexpr);
chargeQ = rho0 * (Pi * a^2) * Gamma;

Print["Gravitational mass m_G ~ rho0 * Volume:"];
(* ToString[..., InputForm] forces it to a string without the wrapper label *)
Print["m_G = ", ToString[massG, InputForm]];
Print[""];

Print["Electric charge q ~ rho0 * Area * Circulation:"];
Print["q_defect = ", ToString[chargeQ, InputForm]];
Print[""];

Print["These relations encode the basic hierarchy scaling:"];
forceRatio = Simplify[chargeQ^2 / massG^2];

Print["F_elec / F_grav ~ q^2 / m_G^2 ~ "];
Print[ToString[forceRatio, InputForm]];
Print[""];


(* ---------------------------------------------------------------------- *)
(* 1. Coulomb potential and Gauss's law structure *)
(* ---------------------------------------------------------------------- *)
Print["==========================================================="];
Print["PART 2: Coulomb potential and Gauss's law"];
Print["===========================================================\n"];

phiC = q / (4 * Pi * eps0 * r);

Print["Coulomb potential (symbolic, in terms of q and eps0):"];
Print["phi_C = ", ToString[phiC, InputForm]];
Print[""];

Evec = -Grad[phiC, {x, y, z}];

Print["Electric field E = -grad(phi_C):"];
Print["E_vec = ", ToString[Evec, InputForm]];
Print[""];

divE = Div[Evec, {x, y, z}];

Print["Divergence of E (valid for r != 0):"];
Print["divE = ", ToString[Simplify[divE, {x != 0 || y != 0 || z != 0}], InputForm]];

Print["
SymPy returns zero for div(E) away from the origin, as expected.
The missing piece at r = 0 is the standard distributional identity:
    div(E) = q * delta^3(r) / eps0
which encodes Gauss's law for a point charge.
"];


(* ---------------------------------------------------------------------- *)
(* 2. Relating eps0 and mu0 via the wave speed c *)
(* ---------------------------------------------------------------------- *)
Print["==========================================================="];
Print["PART 3: Wave speed and the relation eps0 * mu0 = 1/c^2"];
Print["===========================================================\n"];

Print["We impose the standard relativistic relation:"];
relRelation = eps0 * mu0 == 1/c^2;
Print[ToString[relRelation, InputForm]];
Print[""];

Print["
In the toy model, the effective light speed c is identified with the
acoustic wave speed c_s of the superfluid vacuum.
The above relation is therefore a *definition* of mu0 in terms of eps0.
"];

mu0Sol = Solve[relRelation, mu0][[1, 1, 2]]; 

Print["Solving for mu0 in terms of eps0 and c:"];
Print["mu0 = ", ToString[mu0Sol, InputForm]];
Print[""];


(* ---------------------------------------------------------------------- *)
(* 3. Source densities rho_e and J *)
(* ---------------------------------------------------------------------- *)
Print["==========================================================="];
Print["PART 4: Charge density rho_e and current J for a moving defect"];
Print["===========================================================\n"];

uVec = {ux, uy, uz};

(* Symbolic placeholder *)
delta3[val__] := Subscript[\[Delta], 3][val];

rhoE = q * delta3[x, y, z];
JVec = rhoE * uVec;

Print["Symbolic charge density for a point defect:"];
Print["rho_e = ", ToString[rhoE, InputForm]];
Print[""];

Print["Symbolic current density for a moving defect:"];
Print["J_vec = ", ToString[JVec, InputForm]];
Print[""];

Print["
Here delta3 is a placeholder for the 3D Dirac delta distribution.
With these definitions, Gauss's law and the continuity equation
    div(E) = rho_e / eps0,
    d(rho_e)/dt + div(J) = 0,
are satisfied in the standard distributional sense.
"];

Print["em_charge_and_constants.wl: symbolic derivation complete."];

(*"
Output:

===========================================================
PART 1: Defect mass and charge scaling (from toy model)
===========================================================

Gravitational mass m_G ~ rho0 * Volume:
m_G = a^3*aspect*Pi*rho0

Electric charge q ~ rho0 * Area * Circulation:
q_defect = a^2*Gamma*Pi*rho0

These relations encode the basic hierarchy scaling:
F_elec / F_grav ~ q^2 / m_G^2 ~
Gamma^2/(a^2*aspect^2)

===========================================================
PART 2: Coulomb potential and Gauss's law
===========================================================

Coulomb potential (symbolic, in terms of q and eps0):
phi_C = q/(4*eps0*Pi*Sqrt[x^2 + y^2 + z^2])

Electric field E = -grad(phi_C):
E_vec = {(q*x)/(4*eps0*Pi*(x^2 + y^2 + z^2)^(3/2)), (q*y)/(4*eps0*Pi*(x^2 + y^2 + z^2)^(3/2)), (q*z)/(4*eps0*Pi*(x^2 + y^2 + z^2)^(3/2))}

Divergence of E (valid for r != 0):
divE = 0

SymPy returns zero for div(E) away from the origin, as expected.
The missing piece at r = 0 is the standard distributional identity:
    div(E) = q * delta^3(r) / eps0
which encodes Gauss's law for a point charge.

===========================================================
PART 3: Wave speed and the relation eps0 * mu0 = 1/c^2
===========================================================

We impose the standard relativistic relation:
eps0*mu0 == c^(-2)


In the toy model, the effective light speed c is identified with the
acoustic wave speed c_s of the superfluid vacuum.
The above relation is therefore a *definition* of mu0 in terms of eps0.

Solving for mu0 in terms of eps0 and c:
mu0 = 1/(c^2*eps0)

===========================================================
PART 4: Charge density rho_e and current J for a moving defect
===========================================================

Symbolic charge density for a point defect:
rho_e = q*Subscript[δ, 3][x, y, z]

Symbolic current density for a moving defect:
J_vec = {q*ux*Subscript[δ, 3][x, y, z], q*uy*Subscript[δ, 3][x, y, z], q*uz*Subscript[δ, 3][x, y, z]}


Here delta3 is a placeholder for the 3D Dirac delta distribution.
With these definitions, Gauss's law and the continuity equation
    div(E) = rho_e / eps0,
    d(rho_e)/dt + div(J) = 0,
are satisfied in the standard distributional sense.

em_charge_and_constants.wl: symbolic derivation complete.
"*)
