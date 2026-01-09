(* ================================================================= *)
(* MASTER DERIVATION SCRIPT: 4D THROAT MODEL (N2 + L1 + F1) *)
(* Status: FROZEN for Milestone 2 - FIXED *)
(* ================================================================= *)

(* ZONE 1: GEOMETRY & PRIMITIVES (Family-1) *)
ClearAll[t, x, y, z, w, R3, a, L, dr, dw, epsW, V0, p, m, OmOut, OmIn, Zbulk, q, hbar, mu0, xi, mChi];
R3 = Sqrt[x^2 + y^2 + z^2];

(* Smooth Primitives *)
SmoothStep[u_] := (1 + Tanh[u])/2;
SmoothAbs[val_, eps_] := Sqrt[val^2 + eps^2];

(* Gate G(X) *)
uR = (R3 - a[t])/dr;
uW = (SmoothAbs[w, epsW] - L[t]/2)/dw;
Gate = (1 - SmoothStep[uR]) * (1 - SmoothStep[uW]);

(* Confinement Potential Vconf *)
OmegaSq = OmOut^2 - (OmOut^2 - OmIn^2) * Gate;
Vbrane = (1/2) * m * OmegaSq * w^2;
Vwall = V0 * SmoothStep[uR]^p;
Vcap  = V0 * SmoothStep[uW]^p;
Vconf = Vbrane + Vwall + Vcap;

(* ================================================================= *)
(* ZONE 2: THE UNIFIED ACTION (Two-Field N2 + Dielectric L1) *)

(* Fields *)
psi = \!\(\*OverscriptBox[\(\[Psi]\), \(~\)]\)[t, x, y, z, w];
psiC = \!\(\*OverscriptBox[\(\[Psi]C\), \(~\)]\)[t, x, y, z, w];
chi = \!\(\*OverscriptBox[\(\[Chi]\), \(~\)]\)[t, x, y, z, w];
chiC = \!\(\*OverscriptBox[\(\[Chi]C\), \(~\)]\)[t, x, y, z, w];
A = {A0[t,x,y,z,w], Ax[t,x,y,z,w], Ay[t,x,y,z,w], Az[t,x,y,z,w], Aw[t,x,y,z,w]};

(* Covariant Derivatives *)
Deriv[field_, idx_] := D[field, {t,x,y,z,w}[[idx]]];
CovD[field_, idx_] := D[field, {t,x,y,z,w}[[idx]]] - (I*q/hbar)*A[[idx]]*field;
CovDC[field_, idx_] := D[field, {t,x,y,z,w}[[idx]]] + (I*q/hbar)*A[[idx]]*field;

(* Lagrangians *)
(* Note: Psi is Neutral *)
LPsi = (I*hbar/2)*(psiC*Deriv[psi, 1] - psi*Deriv[psiC, 1]) - (hbar^2/(2*m))*Sum[Deriv[psiC, k]*Deriv[psi, k], {k, 2, 5}] - Vconf*psiC*psi;

(* Note: Chi is Charged (Minimal Coupling) *)
LChi = (I*hbar/2)*(chiC*CovD[chi, 1] - chi*CovDC[chiC, 1]) - (hbar^2/(2*mChi))*Sum[CovDC[chiC, k]*CovD[chi, k], {k, 2, 5}];

(* Note: EM with Dielectric Localization Zloc *)
FieldStrength[m_, n_] := D[A[[n]], {t,x,y,z,w}[[m]]] - D[A[[m]], {t,x,y,z,w}[[n]]];
Zloc = 1 + (Zbulk - 1)*(1 - Gate);
LEM = - (1/(4*mu0)) * Zloc * Sum[FieldStrength[m, n]^2, {m, 1, 5}, {n, 1, 5}];
GaugeFix = - (1/(2*xi*mu0)) * (Sum[D[A[[m]], {t,x,y,z,w}[[m]]], {m, 1, 5}])^2;

LTotal = LPsi + LChi + LEM + GaugeFix;

(* ================================================================= *)
(* ZONE 3: DERIVATION ENGINE *)
Needs["VariationalMethods`"];
vars = {t, x, y, z, w};

(* Equations of Motion & Forces *)
(* EomPsi = EulerEquations[LTotal, psiC, vars]; *) (* Calculated but omitted for speed *)
FaDensity = D[LTotal, a[t]];
FLDensity = D[LTotal, L[t]];

(* Currents *)
(* RENAMED VARIABLES TO REMOVE UNDERSCORES *)
J0Vacuum = D[LPsi, A[[1]]];
J0Source = D[LChi, A[[1]]];

(* ================================================================= *)
(* ZONE 4: PHYSICS CHECKS *)

Print["\n--- PHYSICS CHECKS ---"];
Print["Vacuum Neutrality Verified (Should be 0): ", J0Vacuum === 0];
Print["Radial Force Terms: ", Length[FaDensity]];

(* ================================================================= *)
(* ZONE 5: ROBUST PYTHON CODE GENERATION *)

Print["\n--- CODE GENERATION ---"];

PyReplaceRules = {
   "**" -> "**",
   "a(t)" -> "a_t", "L(t)" -> "L_t",
   "Sqrt" -> "np.sqrt", "Tanh" -> "np.tanh", "Sech" -> "sech", "Complex(0,1)" -> "1j",

   (* Fields - Strip Arguments *)
   "OverTilde(ψ)(t,x,y,z,w)" -> "psi",
   "OverTilde(ψC)(t,x,y,z,w)" -> "np.conj(psi)",
   (* Note: Using special handling for Greek Chi *)
   "OverTilde(χ)(t,x,y,z,w)" -> "chi",
   "OverTilde(χC)(t,x,y,z,w)" -> "np.conj(chi)",
   (* Fallback if Mathematica outputs Greek char directly *)
   "OverTilde(chi)(t,x,y,z,w)" -> "chi",
   "OverTilde(chiC)(t,x,y,z,w)" -> "np.conj(chi)",

   "Az(t,x,y,z,w)" -> "Az", "Aw(t,x,y,z,w)" -> "Aw",
   "Ax(t,x,y,z,w)" -> "Ax", "Ay(t,x,y,z,w)" -> "Ay",
   "A0(t,x,y,z,w)" -> "A0"
};

FormatDerivatives[codeStr_] := StringReplace[codeStr, {
   RegularExpression["Derivative\\(1,0,0,0,0\\)\\((A[0xyzw])\\)\\([^)]+\\)"] :> "grad_t($1)",
   RegularExpression["Derivative\\(0,1,0,0,0\\)\\((A[0xyzw])\\)\\([^)]+\\)"] :> "grad_x($1)",
   RegularExpression["Derivative\\(0,0,1,0,0\\)\\((A[0xyzw])\\)\\([^)]+\\)"] :> "grad_y($1)",
   RegularExpression["Derivative\\(0,0,0,1,0\\)\\((A[0xyzw])\\)\\([^)]+\\)"] :> "grad_z($1)",
   RegularExpression["Derivative\\(0,0,0,0,1\\)\\((A[0xyzw])\\)\\([^)]+\\)"] :> "grad_w($1)"
}];

ExportPy[expr_] := Module[{code},
   code = ToString[FortranForm[expr]];
   code = FormatDerivatives[code];
   code = StringReplace[code, PyReplaceRules];
   code
];

(* ================================================================= *)
(* ADDITION TO ZONE 5: SPATIAL CURRENTS *)
(* ================================================================= *)

(* Calculate Spatial Currents: J^k = dL_Chi / dA_k *)
(* We need Jx, Jy, Jz, Jw to drive the spatial wave equations *)
SpatialCurrents = Table[D[LChi, A[[k]]], {k, 2, 5}]; (* Indices 2-5 are x,y,z,w *)

Print["\n# PYTHON: Spatial Currents (Jx, Jy, Jz, Jw)"];
Print["def compute_spatial_currents(chi, A0, Ax, Ay, Az, Aw):"];
Print["    # Returns tuple (Jx, Jy, Jz, Jw)"];
Print["    # Note: These drive the wave equations for Ax, Ay, Az, Aw"];
Print["    Jx = " <> ExportPy[SpatialCurrents[[1]]]];
Print["    Jy = " <> ExportPy[SpatialCurrents[[2]]]];
Print["    Jz = " <> ExportPy[SpatialCurrents[[3]]]];
Print["    Jw = " <> ExportPy[SpatialCurrents[[4]]]];
Print["    return Jx, Jy, Jz, Jw"];

Print["\n# PYTHON: Radial Force Density"];
Print["def compute_force_radial(psi, chi, A0, Ax, Ay, Az, Aw, a_t, L_t):"];
Print["    # Note: Define sech(x) = 1/np.cosh(x)"];
Print["    return " <> ExportPy[FaDensity]];

Print["\n# PYTHON: Maxwell Source (J^0)"];
Print["    # This is the Right-Hand-Side for the A0 equation"];
Print["    return " <> ExportPy[J0Source]];

(*"
Output:

--- PHYSICS CHECKS ---
Vacuum Neutrality Verified (Should be 0): True
Radial Force Terms: 2

--- CODE GENERATION ---

# PYTHON: Radial Force Density
def compute_force_radial(psi, chi, A0, Ax, Ay, Az, Aw, a_t, L_t):
    # Note: Define sech(x) = 1/np.cosh(x)
    return -((-((p*V0*sech((np.sqrt(x**2 + y**2 + z**2) - a_t)/dr)**2*(1 + np.tanh((np.sqrt(x**2 + y**2 + z**2) - a_t)/dr))**(-1 + p))/(2**p*dr)) - (m*(-OmIn**2 + OmOut**2)*w**2*sech((np.sqrt(x**2 + y**2 + z**2) - a_t)/dr)**2*(1 + (-1 - np.tanh((np.sqrt(epsW**2 + w**2) - L_t/2.)/dw))/2.))/(4.*dr))*psi*np.conj(psi)) + ((-1 + Zbulk)*sech((np.sqrt(x**2 + y**2 + z**2) - a_t)/dr)**2*(1 + (-1 - np.tanh((np.sqrt(epsW**2 + w**2) - L_t/2.)/dw))/2.)*((grad_w(Az) - grad_z(Aw))**2 + (-grad_w(Az) + grad_z(Aw))**2 + (grad_w(Ay) - grad_y(Aw))**2 + (-grad_w(Ay) + grad_y(Aw))**2 + (grad_z(Ay) - grad_y(Az))**2 + (-grad_z(Ay) + grad_y(Az))**2 + (grad_w(Ax) - grad_x(Aw))**2 + (-grad_w(Ax) + grad_x(Aw))**2 + (grad_y(Ax) - grad_x(Ay))**2 + (-grad_y(Ax) + grad_x(Ay))**2 + (grad_z(Ax) - grad_x(Az))**2 + (-grad_z(Ax) + grad_x(Az))**2 + (grad_w(A0) - grad_t(Aw))**2 + (-grad_w(A0) + grad_t(Aw))**2 + (grad_x(A0) - grad_t(Ax))**2 + (-grad_x(A0) + grad_t(Ax))**2 + (grad_y(A0) - grad_t(Ay))**2 + (-grad_y(A0) + grad_t(Ay))**2 + (grad_z(A0) - grad_t(Az))**2 + (-grad_z(A0) + grad_t(Az))**2))/(8.*dr*mu0)

# PYTHON: Maxwell Source (J^0)
    # This is the Right-Hand-Side for the A0 equation
    return q*chi*np.conj(chi)
"*)
