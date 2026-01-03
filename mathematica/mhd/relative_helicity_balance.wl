(*
  STEP 6B (FIXED2): Relative helicity balance (finite-domain / boundary-aware)
  Save as: relative_helicity_balance_step6B_fixed2.wl
*)

Print["============================"];
Print["STEP 6B (FIXED2): Relative helicity balance (finite domain)"];
Print["============================\n"];

ClearAll[vGrad, vDiv, vCurl, vDot, vCross, section];

vGrad[f_] := {D[f, x], D[f, y], D[f, z]};
vDiv[v_List] := D[v[[1]], x] + D[v[[2]], y] + D[v[[3]], z];
vCurl[v_List] := {
  D[v[[3]], y] - D[v[[2]], z],
  D[v[[1]], z] - D[v[[3]], x],
  D[v[[2]], x] - D[v[[1]], y]
};
vDot[a_List, b_List] := Expand[a.b];
vCross[a_List, b_List] := Expand[Cross[a, b]];

section[s_] := (Print["\n============================"]; Print[s]; Print["============================\n"];);

(* ---------------------------
   Define fields
---------------------------- *)
A   = {Ax[x, y, z, t],  Ay[x, y, z, t],  Az[x, y, z, t]};
A0  = {Ax0[x, y, z, t], Ay0[x, y, z, t], Az0[x, y, z, t]};

ph  = phi[x, y, z, t];
ph0 = phi0[x, y, z, t];

eField  = {Ex[x, y, z, t],  Ey[x, y, z, t],  Ez[x, y, z, t]};
eField0 = {Ex0[x, y, z, t], Ey0[x, y, z, t], Ez0[x, y, z, t]};

B   = vCurl[A];
B0  = vCurl[A0];

Print["B  = vCurl(A):  ", B];
Print["B0 = vCurl(A0): ", B0];
Print["Check vDiv(B)  (should be 0):  ", Simplify[vDiv[B]]];
Print["Check vDiv(B0) (should be 0):  ", Simplify[vDiv[B0]]];

(* Gauge definitions: E = -A_t - ∇phi, same for reference *)
At  = -eField  - vGrad[ph];
A0t = -eField0 - vGrad[ph0];

Bt  = vCurl[At];
B0t = vCurl[A0t];

section["STEP 6B.1: Core relative helicity identity (with gauge substitution)"];

hR  = vDot[A + A0, B - B0];
htR = vDot[At + A0t, B - B0] + vDot[A + A0, Bt - B0t];

fluxR = vCross[eField - eField0, A + A0] + (ph + ph0) (B - B0);

resR = Simplify[Expand[ htR + vDiv[fluxR] + 2 (vDot[eField, B] - vDot[eField0, B0]) ]];

Print["Check: ∂t hR + Div(fluxR) + 2(E·B - E0·B0) (should be 0): ", resR];

Print["\nPaper-ready local balance:"];
Print["  hR = (A + A0)·(B - B0)"];
Print["  fluxR = (E - E0)×(A + A0) + (phi + phi0)(B - B0)"];
Print["  ∂t hR + ∇·fluxR = -2 (E·B - E0·B0)"];

section["STEP 6B.2 (Optional): Reduce E·B under Ohm's law (generic, no recursion)"];

(* Generic proof: avoids expanding B=curl(A) which triggers recursion explosions *)
ClearAll[uu, vv, ww, BBx, BBy, BBz, JJx, JJy, JJz, ηsym];

vSym = {uu, vv, ww};
BSym = {BBx, BBy, BBz};
JSym = {JJx, JJy, JJz};

eOhmSym = -Cross[vSym, BSym] + ηsym JSym;

checkEBSym = FullSimplify[
  Expand[eOhmSym.BSym - ηsym (JSym.BSym)]
];

Print["Check: ( (-v×B + ηJ)·B ) - η(J·B) (generic)  (should be 0): ", checkEBSym];

Print["\nInterpretation:"];
Print["  Since (v×B)·B = 0 identically, we get E·B = η (J·B)."];
Print["  If the reference field is chosen current-free (J0=0), then E0·B0 = 0."];
Print["  RHS becomes -2 η (J·B), with boundary transport handled by fluxR."];

Print["\nDone."];

(*"
Output:

============================
STEP 6B (FIXED2): Relative helicity balance (finite domain)
============================

B  = vCurl(A):  {-Derivative[0, 0, 1, 0][Ay][x, y, z, t] + Derivative[0, 1, 0, 0][Az][x, y, z, t], Derivative[0, 0, 1, 0][Ax][x, y, z, t] - Derivative[1, 0, 0, 0][Az][x, y, z, t], -Derivative[0, 1, 0, 0][Ax][x, y, z, t] + Derivative[1, 0, 0, 0][Ay][x, y, z, t]}
B0 = vCurl(A0): {-Derivative[0, 0, 1, 0][Ay0][x, y, z, t] + Derivative[0, 1, 0, 0][Az0][x, y, z, t], Derivative[0, 0, 1, 0][Ax0][x, y, z, t] - Derivative[1, 0, 0, 0][Az0][x, y, z, t], -Derivative[0, 1, 0, 0][Ax0][x, y, z, t] + Derivative[1, 0, 0, 0][Ay0][x, y, z, t]}
Check vDiv(B)  (should be 0):  0
Check vDiv(B0) (should be 0):  0

============================
STEP 6B.1: Core relative helicity identity (with gauge substitution)
============================

Check: ∂t hR + Div(fluxR) + 2(E·B - E0·B0) (should be 0): 0

Paper-ready local balance:
  hR = (A + A0)·(B - B0)
  fluxR = (E - E0)×(A + A0) + (phi + phi0)(B - B0)
  ∂t hR + ∇·fluxR = -2 (E·B - E0·B0)

============================
STEP 6B.2 (Optional): Reduce E·B under Ohm's law (generic, no recursion)
============================

Check: ( (-v×B + ηJ)·B ) - η(J·B) (generic)  (should be 0): 0
"*)
