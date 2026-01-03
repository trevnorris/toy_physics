(* ============================================
   STEP 6A: Vector potential formulation consistency
   - B = curl(A) is divergence-free identically
   - A_t = v×B - η J - ∇φ  implies:
       B_t = curl(v×B) - curl(η curl(B))
   - Resistive part dissipates magnetic energy
   - Helicity balance consistency check
   ============================================ *)

ClearAll[
  vGrad, vDiv, vCurl, vCross, vDot, vNorm2,
  x, y, z, t,
  Ax, Ay, Az, u, vcomp, w, phi, eta
];

(* --- Vector calculus helpers (do NOT override built-ins) --- *)
vGrad[f_] := {D[f, x], D[f, y], D[f, z]};

vDiv[v_List] := D[v[[1]], x] + D[v[[2]], y] + D[v[[3]], z];

vCurl[v_List] := {
  D[v[[3]], y] - D[v[[2]], z],
  D[v[[1]], z] - D[v[[3]], x],
  D[v[[2]], x] - D[v[[1]], y]
};

vCross[a_List, b_List] := Cross[a, b];
vDot[a_List, b_List] := a.b;
vNorm2[a_List] := a.a;

section[title_] := Print[
  "\n============================\n",
  title,
  "\n============================\n"
];

(* --- Fields --- *)
A      = {Ax[x, y, z, t], Ay[x, y, z, t], Az[x, y, z, t]};
vField = {u[x, y, z, t],  vcomp[x, y, z, t], w[x, y, z, t]};
ph     = phi[x, y, z, t];
η      = eta[x, y, z, t];

(* --- Derived fields --- *)
B = vCurl[A];
J = vCurl[B];

(* --- Define E and A_t (gauge form) --- *)
eField = -vCross[vField, B] + η J;                (* E = -v×B + ηJ *)
At     =  vCross[vField, B] - η J - vGrad[ph];    (* A_t = -E - ∇φ *)

(* --- Implied B_t from A_t --- *)
Bt = vCurl[At];

section["STEP 6A: Vector potential formulation -> induction equation"];

Print["B = vCurl(A): ", B];
Print["Check vDiv(B) (should be 0): ", Simplify[vDiv[B]]];

expectedBt = vCurl[vCross[vField, B]] - vCurl[η vCurl[B]];

Print["\nCheck Bt - expectedBt (should be {0,0,0}): ",
  Simplify[Bt - expectedBt]
];

Print["Check vDiv(Bt) (should be 0): ", Simplify[vDiv[Bt]]];

section["STEP 6A.1: Gauge consistency checks"];

Print["Check At + eField + vGrad(phi) (should be {0,0,0}): ",
  Simplify[At + eField + vGrad[ph]]
];

Print["Check Bt + vCurl(eField) (should be {0,0,0}): ",
  Simplify[Bt + vCurl[eField]]
];

section["STEP 6A.2: Resistive operator is exactly the recommended one"];

BtRes = vCurl[-η J];

Print["BtRes (resistive part) = -vCurl(η J)"];
Print["Check BtRes + vCurl(η vCurl(B)) (should be {0,0,0}): ",
  Simplify[BtRes + vCurl[η vCurl[B]]]
];

(* FIXED LINE (this is what was broken before) *)
Print["Check vDiv(BtRes) (should be 0): ", Simplify[vDiv[BtRes]]];

section["STEP 6A.3: Magnetic energy dissipation (resistive part)"];

fluxTerm = -vDiv[vCross[η J, B]];
dissTerm = -η vNorm2[J];

Print["Check: B·BtRes - (fluxTerm + dissTerm) (should be 0): ",
  Simplify[vDot[B, BtRes] - (fluxTerm + dissTerm)]
];

Print["Result (local, resistive):  ∂t(|B|^2/2) = -∇·((η J)×B) - η |J|^2"];
Print["Integrated (periodic/decaying): d/dt ∫|B|^2 dV = -2 ∫ η |J|^2 dV <= 0"];

section["STEP 6A.4: Helicity balance consistency (via A_t,B_t)"];

hDot   = vDot[At, B] + vDot[A, Bt];               (* ∂t(A·B) via product rule *)
hFlux  = vCross[eField, A] + ph B;                (* flux = E×A + phi B *)

resHel = FullSimplify[Expand[hDot + vDiv[hFlux] + 2 vDot[eField, B]]];

Print["Check: ht + vDiv(flux) + 2(E·B) (should be 0): ", resHel];
Print["Check: (E·B) - η (J·B) (should be 0): ",
  FullSimplify[Expand[vDot[eField, B] - η vDot[J, B]]]
];

Print["Paper-ready local balance:"];
Print["  ∂t(A·B) + ∇·( E×A + phi B ) = -2 E·B = -2 η (J·B)"];

Print["\nDone."];

(*"
Output:


============================
STEP 6A: Vector potential formulation -> induction equation
============================

B = vCurl(A): {-Derivative[0, 0, 1, 0][Ay][x, y, z, t] + Derivative[0, 1, 0, 0][Az][x, y, z, t], Derivative[0, 0, 1, 0][Ax][x, y, z, t] - Derivative[1, 0, 0, 0][Az][x, y, z, t], -Derivative[0, 1, 0, 0][Ax][x, y, z, t] + Derivative[1, 0, 0, 0][Ay][x, y, z, t]}
Check vDiv(B) (should be 0): 0

Check Bt - expectedBt (should be {0,0,0}): {0, 0, 0}
Check vDiv(Bt) (should be 0): 0

============================
STEP 6A.1: Gauge consistency checks
============================

Check At + eField + vGrad(phi) (should be {0,0,0}): {0, 0, 0}
Check Bt + vCurl(eField) (should be {0,0,0}): {0, 0, 0}

============================
STEP 6A.2: Resistive operator is exactly the recommended one
============================

BtRes (resistive part) = -vCurl(η J)
Check BtRes + vCurl(η vCurl(B)) (should be {0,0,0}): {0, 0, 0}
Check vDiv(BtRes) (should be 0): 0

============================
STEP 6A.3: Magnetic energy dissipation (resistive part)
============================

Check: B·BtRes - (fluxTerm + dissTerm) (should be 0): 0
Result (local, resistive):  ∂t(|B|^2/2) = -∇·((η J)×B) - η |J|^2
Integrated (periodic/decaying): d/dt ∫|B|^2 dV = -2 ∫ η |J|^2 dV <= 0

============================
STEP 6A.4: Helicity balance consistency (via A_t,B_t)
============================

Check: ht + vDiv(flux) + 2(E·B) (should be 0): 0
Check: (E·B) - η (J·B) (should be 0): 0
Paper-ready local balance:
  ∂t(A·B) + ∇·( E×A + phi B ) = -2 E·B = -2 η (J·B)
"*)
