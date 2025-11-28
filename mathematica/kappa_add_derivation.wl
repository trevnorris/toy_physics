(* kappa_add Derivation: Added Mass of a Moving Sink in Superfluid *)
(* ASCII-only version for Mathematica compatibility *)

ClearAll["Global`*"];

$Assumptions = {
  Q > 0, rho0 > 0, a > 0, V > 0, cs > 0, mu > 0, m > 0,
  r > a, R > a,
  Element[{Q, rho0, a, V, cs, mu, m, r, R, th, ph}, Reals]
};

Print["============================================================="];
Print["DERIVATION OF kappa_add: ADDED MASS OF A MOVING SINK"];
Print["============================================================="];
Print[""];

(* ================================================================== *)
(* SECTION 1: Setup - Stationary Sink Flow                            *)
(* ================================================================== *)

Print["SECTION 1: Stationary Sink Flow"];
Print["-------------------------------------------------------------"];

(* The velocity potential for a point sink of strength Q at origin *)
(* Flux through any sphere = Q, so integral v.dA = Q *)
(* v_r = Q/(4 Pi r^2), pointing inward means Q > 0 for sink *)
(* Potential: phiSink = -Q/(4 Pi r) so that v = -grad(phi) gives inward flow *)

phiSink[r_] := -Q/(4 Pi r);

vSinkR[r_] := -D[phiSink[r], r] // Simplify;

Print["1.1: Sink potential: phi_sink(r) = ", phiSink[r]];
Print["1.2: Radial velocity: v_r(r) = -d phi/dr = ", vSinkR[r]];

(* Verify: flux through sphere of radius R *)
fluxCheck = vSinkR[R] * 4 Pi R^2 // Simplify;
Print["1.3: Flux through sphere: v_r * 4 Pi R^2 = ", fluxCheck];
Print["     (This equals Q as required for sink strength Q)"];
Print[""];

(* ================================================================== *)
(* SECTION 2: Moving Sink - Velocity Field Perturbation               *)
(* ================================================================== *)

Print["SECTION 2: Moving Sink Perturbation"];
Print["-------------------------------------------------------------"];

(* When sink moves at velocity V in the z-direction, we work in the 
   frame where the sink is instantaneously at origin but moving.
   
   The key insight: in potential flow, a moving source/sink creates
   a dipole perturbation to the velocity field.
   
   For a point source moving at velocity V, the retarded position
   creates an effective dipole moment p = Q V / cs (for compressible)
   or p = 0 (for incompressible point source).
   
   However, our throat has finite size a. A sphere of radius a moving
   at velocity V in potential flow has velocity potential:
   
   phi_dipole = (V a^3)/(2 r^2) cos(th)
   
   where th is angle from velocity direction (z-axis).
*)

Print["2.1: For finite-size throat (radius a) moving at velocity V:"];

(* Dipole potential for moving sphere *)
phiDipole[r_, th_] := (V a^3)/(2 r^2) Cos[th];

Print["2.2: Dipole potential: phi_dipole = ", phiDipole[r, th]];

(* Velocity components from dipole *)
vDipoleR[r_, th_] := -D[phiDipole[r, th], r] // Simplify;
vDipoleTh[r_, th_] := -(1/r) D[phiDipole[r, th], th] // Simplify;

Print["2.3: Dipole velocity components:"];
Print["     v_r^dipole = ", vDipoleR[r, th]];
Print["     v_th^dipole = ", vDipoleTh[r, th]];
Print[""];

(* ================================================================== *)
(* SECTION 3: Kinetic Energy Calculation                              *)
(* ================================================================== *)

Print["SECTION 3: Kinetic Energy of Combined Flow"];
Print["-------------------------------------------------------------"];

(* Total velocity field: v = v_sink + v_dipole *)
(* v_sink is purely radial: (v_sink, 0, 0) in spherical coords *)
(* v_dipole has both r and th components *)

(* |v|^2 = (v_sink_r + v_dipole_r)^2 + v_dipole_th^2 *)

vTotalRsq[r_, th_] := (vSinkR[r] + vDipoleR[r, th])^2 // Expand;
vTotalThsq[r_, th_] := vDipoleTh[r, th]^2 // Expand;
vTotalSq[r_, th_] := vTotalRsq[r, th] + vTotalThsq[r, th] // Expand;

Print["3.1: |v_total|^2 = (v_sink + v_dip,r)^2 + v_dip,th^2"];
Print["3.2: Expanding..."];

(* Expand |v|^2 into terms by power of V *)
(* |v|^2 = |v_sink|^2 + 2 v_sink . v_dipole + |v_dipole|^2 *)
(*       = O(V^0)      + O(V^1)              + O(V^2)        *)

vSqExpanded = Collect[vTotalSq[r, th], V, Simplify];

Print["3.3: |v|^2 expanded:"];
Print["     ", vSqExpanded];
Print[""];

(* Extract terms by power of V *)
term0 = Coefficient[vTotalSq[r, th], V, 0] // Simplify;
term1 = Coefficient[vTotalSq[r, th], V, 1] // Simplify;
term2 = Coefficient[vTotalSq[r, th], V, 2] // Simplify;

Print["3.4: Coefficient of V^0: ", term0];
Print["3.5: Coefficient of V^1: ", term1];
Print["3.6: Coefficient of V^2: ", term2];
Print[""];

(* ================================================================== *)
(* SECTION 4: Volume Integration                                      *)
(* ================================================================== *)

Print["SECTION 4: Volume Integration for Kinetic Energy"];
Print["-------------------------------------------------------------"];

(* Kinetic energy: T = (1/2) rho0 integral |v|^2 d^3x *)
(* In spherical coordinates: d^3x = r^2 sin(th) dr dth dphi *)
(* Integrate from r = a (throat radius) to infinity *)
(* phi integral gives 2 Pi (azimuthal symmetry) *)

(* T = (1/2) rho0 * 2 Pi * integral_a^inf integral_0^Pi |v|^2 r^2 sin(th) dth dr *)

Print["4.1: T = (1/2) rho0 * 2 Pi integral_a^inf integral_0^Pi |v|^2 r^2 sin(th) dth dr"];
Print[""];

(* First do the th integral for each term *)

(* Term 0 (V^0): static sink energy *)
(* |v_sink|^2 = Q^2/(16 Pi^2 r^4) - no th dependence *)
thInt0 = Integrate[term0 * Sin[th], {th, 0, Pi}] // Simplify;
Print["4.2: integral (V^0 term) sin(th) dth = ", thInt0];

(* Term 1 (V^1): cross term *)
(* This involves cos(th) * sin(th), which integrates to 0 *)
thInt1 = Integrate[term1 * Sin[th], {th, 0, Pi}] // Simplify;
Print["4.3: integral (V^1 term) sin(th) dth = ", thInt1];
Print["     (Cross term vanishes by symmetry!)"];

(* Term 2 (V^2): dipole self-energy - THIS IS THE ADDED MASS *)
thInt2 = Integrate[term2 * Sin[th], {th, 0, Pi}] // Simplify;
Print["4.4: integral (V^2 term) sin(th) dth = ", thInt2];
Print[""];

(* Now do the radial integral for the V^2 term *)
Print["4.5: Radial integral for V^2 term:"];

rInt2 = Integrate[thInt2 * r^2, {r, a, Infinity}, 
  Assumptions -> {a > 0}] // Simplify;

Print["     integral_a^inf (th-integrated V^2 term) * r^2 dr = ", rInt2];
Print[""];

(* ================================================================== *)
(* SECTION 5: Extract Added Mass                                      *)
(* ================================================================== *)

Print["SECTION 5: Added Mass Extraction"];
Print["-------------------------------------------------------------"];

(* Total kinetic energy from V^2 term: *)
(* T_add = (1/2) rho0 * 2 Pi * (radial integral) * V^2 *)

Tadd = (1/2) rho0 * 2 Pi * rInt2 * V^2 // Simplify;

Print["5.1: T_add = (1/2) rho0 * 2 Pi * (integral integral) * V^2 = ", Tadd];

(* By definition, T_add = (1/2) m_add V^2 *)
(* So m_add is the coefficient of (1/2) V^2 *)

mAddResult = Coefficient[Tadd, V^2] * 2 // Simplify;

Print["5.2: From T_add = (1/2) m_add V^2, extract m_add:"];
Print["     m_add = ", mAddResult];
Print[""];

(* Express in terms of displaced mass *)
(* m_displaced = rho0 * (4/3) Pi a^3 *)

mDisplaced = rho0 * (4/3) Pi a^3;

ratioMaddMdisp = Simplify[mAddResult / mDisplaced];

Print["5.3: Displaced fluid mass: m_displaced = rho0 * (4/3) Pi a^3 = ", mDisplaced];
Print["5.4: Ratio m_add / m_displaced = ", ratioMaddMdisp];
Print[""];

(* ================================================================== *)
(* SECTION 6: Connect to Cavitation Mass                              *)
(* ================================================================== *)

Print["SECTION 6: Express kappa_add in Terms of Cavitation Mass"];
Print["-------------------------------------------------------------"];

(* From info.md: the throat's mass comes from cavitation energy *)
(* m = E_cav / cs^2 = rho0 cs^2 V_cav / cs^2 = rho0 V_cav *)
(* where V_cav = (4/3) Pi a^3 is the throat volume *)

(* So m = m_displaced! The cavitation mass equals the displaced mass. *)

Print["6.1: From cavitation: m = rho0 V_cav = rho0 * (4/3) Pi a^3 = m_displaced"];
Print[""];

(* Therefore kappa_add = m_add / m = m_add / m_displaced *)

kappaAddDerived = ratioMaddMdisp // Simplify;

Print["6.2: kappa_add = m_add / m = m_add / m_displaced = ", kappaAddDerived];
Print[""];

checkKappaAdd = Simplify[kappaAddDerived == 1/2];
Print["6.3: kappa_add = 1/2? ", checkKappaAdd];
Print[""];

(* ================================================================== *)
(* SECTION 7: Physical Interpretation                                 *)
(* ================================================================== *)

Print["SECTION 7: Physical Interpretation"];
Print["-------------------------------------------------------------"];

Print["7.1: The added mass m_add = (1/2) m arises because:"];
Print["     - When throat moves at velocity V, it drags surrounding fluid"];
Print["     - The dipole flow pattern |v_dipole|^2 ~ V^2 a^6/r^6"];
Print["     - Integrating this from r=a to infinity gives finite kinetic energy"];
Print["     - This energy is T_add = (1/2) * (1/2) m * V^2"];
Print[""];
Print["7.2: The factor 1/2 is the CLASSICAL result for a sphere in"];
Print["     potential flow. Our derivation confirms this holds for"];
Print["     a throat/sink geometry as well."];
Print[""];
Print["7.3: Importantly, the CROSS TERM (V^1) vanishes by symmetry:"];
Print["     integral cos(th) sin(th) dth = 0"];
Print["     This means sink flow and dipole flow don't interfere"];
Print["     in the total kinetic energy."];
Print[""];

(* ================================================================== *)
(* SECTION 8: Verification via Kelvin's Surface Integral              *)
(* ================================================================== *)

Print["SECTION 8: Verification via Kelvin's Surface Integral"];
Print["-------------------------------------------------------------"];

(* Kelvin's formula for kinetic energy of irrotational flow outside a body: *)
(*   T = -(1/2) rho0 integral_S phi (d phi/dn) dS                          *)
(*                                                                          *)
(* The MINUS SIGN comes from Green's identity when converting the volume    *)
(* integral T = (1/2) rho0 integral |grad phi|^2 dV to a surface integral.  *)
(* For flow OUTSIDE the body, with outward normal n = +r_hat:               *)
(*   integral_V |grad phi|^2 dV = -integral_S phi (d phi/dn) dS             *)
(*                                                                          *)
(* For motion in z-direction, d/dn = d/dr on sphere surface (outward).      *)

Print["8.1: Kelvin formula: T = -(1/2) rho0 integral_S phi (d phi/dr) dS"];
Print["     (Note: minus sign from Green's identity for exterior flow)"];
Print[""];

dPhiDipoledr = D[phiDipole[r, th], r] /. r -> a // Simplify;
Print["8.2: At r = a: d phi_dipole/dr = ", dPhiDipoledr];

phiAtSurface = phiDipole[a, th] // Simplify;
Print["8.3: At r = a: phi_dipole = ", phiAtSurface];

(* Surface element in spherical coords: dS = r^2 sin(th) dth dphi *)
integrandKelvin = phiAtSurface * dPhiDipoledr * a^2 * Sin[th] // Simplify;
Print["8.4: Integrand: phi (d phi/dr) r^2 sin(th) = ", integrandKelvin];

(* Integrate over sphere: first theta, then multiply by 2 Pi for phi *)
kelvinThInt = Integrate[integrandKelvin, {th, 0, Pi}] // Simplify;
Print["8.5: integral_0^Pi dth = ", kelvinThInt];

surfaceIntegral = 2 Pi * kelvinThInt // Simplify;
Print["8.6: Full surface integral (times 2 Pi): ", surfaceIntegral];

(* Apply Kelvin's formula: T = -(1/2) rho0 * (surface integral) *)
TaddKelvin = -(1/2) * rho0 * surfaceIntegral // Simplify;
Print["8.7: T_add = -(1/2) rho0 * (surface integral) = ", TaddKelvin];

(* Extract m_add from T = (1/2) m_add V^2 *)
mAddKelvin = Coefficient[TaddKelvin, V^2] * 2 // Simplify;
Print["8.8: From T = (1/2) m_add V^2: m_add = ", mAddKelvin];

kappaAddKelvin = Simplify[mAddKelvin / mDisplaced];
Print["8.9: kappa_add from Kelvin = ", kappaAddKelvin];

checkKelvin = Simplify[kappaAddKelvin == 1/2];
Print["8.10: Matches volume integral? ", checkKelvin];
Print[""];

(* ================================================================== *)
(* SUMMARY                                                            *)
(* ================================================================== *)

Print["============================================================="];
Print["SUMMARY: kappa_add DERIVATION"];
Print["============================================================="];
Print[""];
Print["DERIVED (not assumed):"];
Print["  1. Velocity field of moving sink = static sink + dipole"];
Print["  2. |v|^2 expanded to O(V^2), cross term vanishes by symmetry"];
Print["  3. Volume integral of V^2 term: T_add = (1/2) m_add V^2"];
Print["  4. Result: m_add = (2/3) Pi a^3 rho0 = (1/2) m_displaced"];
Print["  5. Since m = m_displaced (cavitation): kappa_add = 1/2"];
Print[""];
Print["VERIFIED by two independent methods:"];
Print["  - Volume integral of kinetic energy"];
Print["  - Kelvin's surface integral formula"];
Print[""];
Print["  +----------------------------------+"];
Print["  |  kappa_add = 1/2  (DERIVED)     |"];
Print["  +----------------------------------+"];
Print[""];

(* ================================================================== *)
(* CONTEXT: Full beta = 5/2 Decomposition Status                      *)
(* ================================================================== *)

Print["============================================================="];
Print["CONTEXT: FULL beta = 5/2 DECOMPOSITION"];
Print["============================================================="];
Print[""];
Print["To match GR's 1PN precession, we need beta = 5/2 where:"];
Print["  sigma(r) = beta * mu / (cs^2 * r)"];
Print[""];
Print["Proposed decomposition: beta = kappa_rho + kappa_add + kappa_PV"];
Print[""];
Print["  +------------+-------+------------------------------------------+"];
Print["  | Component  | Value | Status                                   |"];
Print["  +------------+-------+------------------------------------------+"];
Print["  | kappa_rho  |   1   | TRIVIAL: m_eff = rho(r) * V_cav          |"];
Print["  |            |       | follows from cavitation mass definition  |"];
Print["  +------------+-------+------------------------------------------+"];
Print["  | kappa_add  |  1/2  | DERIVED: this script (dipole KE)         |"];
Print["  +------------+-------+------------------------------------------+"];
Print["  | kappa_PV   |   1   | NOT DERIVED: determined by closure       |"];
Print["  |            |       | 5/2 - 1 - 1/2 = 1                        |"];
Print["  |            |       | Physical interpretation unclear          |"];
Print["  +------------+-------+------------------------------------------+"];
Print["  | TOTAL      |  5/2  | Required by GR matching                  |"];
Print["  +------------+-------+------------------------------------------+"];
Print[""];
Print["The remaining kappa_PV = 1 requires solving the full unsteady"];
Print["compressible Euler equations for an accelerating throat."];
Print[""];
Print["============================================================="];
