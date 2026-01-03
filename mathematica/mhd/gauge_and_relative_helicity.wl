(* ::Package:: *)

(*
STEP 5: Gauge invariance + (Relative) Magnetic Helicity notes
Purpose:
  - Show B = curl(A) is gauge-invariant under A -> A + grad(chi)
  - Show helicity density A·B shifts by a pure divergence (so the volume integral is invariant under suitable boundaries)
  - Show E from potentials is gauge-invariant if phi -> phi - d_t chi
  - Show the local helicity balance identity is gauge-consistent
  - (Optional) Show relative helicity density is gauge-invariant up to a divergence
*)

Print["============================"];
Print["STEP 5: Gauge / (Relative) Helicity identities"];
Print["============================\n"];

(* --- Safe vector calculus helpers (do NOT touch protected built-ins) --- *)
ClearAll[vGrad, vDiv, vCurl, vDot, vCross];

vGrad[f_] := {D[f, x], D[f, y], D[f, z]};
vDiv[v_List] := D[v[[1]], x] + D[v[[2]], y] + D[v[[3]], z];
vCurl[v_List] := {
  D[v[[3]], y] - D[v[[2]], z],
  D[v[[1]], z] - D[v[[3]], x],
  D[v[[2]], x] - D[v[[1]], y]
};
vDot[a_List, b_List] := a.b;
vCross[a_List, b_List] := Cross[a, b];

(* --- Fields / potentials --- *)
A    = {Ax[x, y, z, t],  Ay[x, y, z, t],  Az[x, y, z, t]};
phi  =  phiFun[x, y, z, t];
chi  =  chiFun[x, y, z, t];

B  = vCurl[A];

Print["B = vCurl(A): ", B];
Print["Check vDiv(B) (should be 0): ", Simplify[vDiv[B]], "\n"];

(* --- Gauge transform: A -> A + grad chi --- *)
Aprime = A + vGrad[chi];
Bprime = vCurl[Aprime];

Print["============================"];
Print["STEP 5A: B is gauge invariant"];
Print["============================\n"];
Print["Check B' - B (should be {0,0,0}): ", Simplify[Bprime - B], "\n"];

(* --- Helicity density shift under gauge --- *)
h      = vDot[A, B];
hprime = vDot[Aprime, Bprime];

(* Identity target:
   (A + ∇chi)·B = A·B + ∇·(chi B) - chi ∇·B
   Since ∇·B = 0, the difference is a pure divergence.
*)
gaugeShiftCheck = Simplify[hprime - h - vDiv[chi B]];

Print["============================"];
Print["STEP 5B: Helicity density shifts by a divergence"];
Print["============================\n"];
Print["Check: h' - h - vDiv(chi B)  (should be -chi vDiv(B)): ", gaugeShiftCheck];
Print["With vDiv(B)=0, this becomes: ", Simplify[gaugeShiftCheck /. vDiv[B] -> 0], "\n"];

Print["Interpretation:"];
Print["  h' = h + ∇·(chi B)  (since ∇·B=0)."];
Print["  So ∫ h dV is gauge-invariant if the surface flux ∮ chi (B·n) dS vanishes"];
Print["  (e.g., periodic domain, decaying fields, perfectly conducting walls with suitable gauge),\n"];

(* --- Electric field from potentials and gauge invariance --- *)
eFromPot[Avec_List, ph_] := -D[Avec, t] - vGrad[ph];

phiPrime = phi - D[chi, t];

Ee      = eFromPot[A, phi];
Eprime = eFromPot[Aprime, phiPrime];

Print["============================"];
Print["STEP 5C: E from (A,phi) is gauge-invariant if phi -> phi - chi_t"];
Print["============================\n"];
Print["Check E' - E (should be {0,0,0}): ", Simplify[Eprime - Ee], "\n"];

(* --- Local helicity identity consistency under gauge --- *)
(* Standard identity (purely kinematic, no constitutive assumptions):
   If B = curl(A) and E = -A_t - grad(phi), then
     ∂t(A·B) + ∇·(E×A + phi B) = -2 E·B
*)
helicityResidual[Avec_List, ph_] := Module[{Bv, Ev, hv, flux},
  Bv   = vCurl[Avec];
  Ev   = eFromPot[Avec, ph];
  hv   = vDot[Avec, Bv];
  flux = vCross[Ev, Avec] + ph Bv;
  Simplify[D[hv, t] + vDiv[flux] + 2 vDot[Ev, Bv]]
];

res  = helicityResidual[A, phi];
resP = helicityResidual[Aprime, phiPrime];

Print["============================"];
Print["STEP 5D: Local helicity balance is gauge-consistent"];
Print["============================\n"];
Print["Residual for (A,phi)    (should be 0): ", res];
Print["Residual for (A',phi')  (should be 0): ", resP];
Print["Check (resP - res)      (should be 0): ", Simplify[resP - res], "\n"];

Print["Result (paper-ready):"];
Print["  Under A -> A + ∇chi, phi -> phi - ∂t chi, we have B invariant, E invariant,"];
Print["  and A·B shifts only by a divergence term. So the integrated helicity is well-defined"];
Print["  under appropriate boundary conditions.\n"];

(* --- OPTIONAL: Relative helicity density (one common form) --- *)
(* Define a reference potential Aref with Bref = curl(Aref) *)
Aref   = {Ax0[x, y, z, t], Ay0[x, y, z, t], Az0[x, y, z, t]};
Bref   = vCurl[Aref];
chi0   = chi0Fun[x, y, z, t];

(* One common relative helicity density choice:
   hR = (A + Aref) · (B - Bref)
   Gauge shifts A -> A + ∇chi, Aref -> Aref + ∇chi0 change hR by a divergence if div(B)=div(Bref)=0.
*)
hR      = vDot[A + Aref, B - Bref];
hRprime = vDot[(A + vGrad[chi]) + (Aref + vGrad[chi0]), B - Bref];

relShiftCheck = Simplify[hRprime - hR - vDiv[(chi + chi0) (B - Bref)]];

Print["============================"];
Print["STEP 5E (Optional): Relative helicity density shifts by a divergence"];
Print["============================\n"];
Print["Check: hR' - hR - vDiv((chi+chi0)(B-Bref)) (should be -(chi+chi0) vDiv(B-Bref)):"];
Print[relShiftCheck];
Print["Since B and Bref are curls, vDiv(B-Bref) = ", Simplify[vDiv[B - Bref]], " so the check reduces to: ",
  Simplify[relShiftCheck /. vDiv[B - Bref] -> 0], "\n"];

Print["Done."];

(*"
Output:

============================
STEP 5: Gauge / (Relative) Helicity identities
============================

B = vCurl(A): {-Derivative[0, 0, 1, 0][Ay][x, y, z, t] + Derivative[0, 1, 0, 0][Az][x, y, z, t], Derivative[0, 0, 1, 0][Ax][x, y, z, t] - Derivative[1, 0, 0, 0][Az][x, y, z, t], -Derivative[0, 1, 0, 0][Ax][x, y, z, t] + Derivative[1, 0, 0, 0][Ay][x, y, z, t]}
Check vDiv(B) (should be 0): 0

============================
STEP 5A: B is gauge invariant
============================

Check B' - B (should be {0,0,0}): {0, 0, 0}

============================
STEP 5B: Helicity density shifts by a divergence
============================

Check: h' - h - vDiv(chi B)  (should be -chi vDiv(B)): 0
With vDiv(B)=0, this becomes: 0

Interpretation:
  h' = h + ∇·(chi B)  (since ∇·B=0).
  So ∫ h dV is gauge-invariant if the surface flux ∮ chi (B·n) dS vanishes
  (e.g., periodic domain, decaying fields, perfectly conducting walls with suitable gauge),

============================
STEP 5C: E from (A,phi) is gauge-invariant if phi -> phi - chi_t
============================

Check E' - E (should be {0,0,0}): {0, 0, 0}

============================
STEP 5D: Local helicity balance is gauge-consistent
============================

Residual for (A,phi)    (should be 0): 0
Residual for (A',phi')  (should be 0): 0
Check (resP - res)      (should be 0): 0

Result (paper-ready):
  Under A -> A + ∇chi, phi -> phi - ∂t chi, we have B invariant, E invariant,
  and A·B shifts only by a divergence term. So the integrated helicity is well-defined
  under appropriate boundary conditions.

============================
STEP 5E (Optional): Relative helicity density shifts by a divergence
============================

Check: hR' - hR - vDiv((chi+chi0)(B-Bref)) (should be -(chi+chi0) vDiv(B-Bref)):
0
Since B and Bref are curls, vDiv(B-Bref) = 0 so the check reduces to: 0
"*)
